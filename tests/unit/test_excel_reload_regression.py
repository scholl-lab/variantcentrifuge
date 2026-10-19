"""Regression tests for eliminating openpyxl reloads during Excel generation (#159)."""

import json
from pathlib import Path
from unittest.mock import patch

import openpyxl
import pandas as pd
import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.output_stages import ExcelReportStage


@pytest.fixture
def excel_context(tmp_path: Path) -> PipelineContext:
    """Create a pipeline context with all optional sheets and IGV report mapping."""
    output_dir = tmp_path / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir = output_dir / "report" / "igv"
    report_dir.mkdir(parents=True, exist_ok=True)

    # 1. Setup IGV reports map
    igv_map = {
        "variants": [
            {
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "T",
                "sample_reports": {"Sample1": "igv/Sample1_var1.html"},
            },
            {
                "chrom": "chr2",
                "pos": 200,
                "ref": "C",
                "alt": "G",
                "sample_reports": {
                    "Sample1": "igv/Sample1_var2.html",
                    "Sample2": "igv/Sample2_var2.html",
                },
            },
        ]
    }
    with open(report_dir / "igv_reports_map.json", "w", encoding="utf-8") as f:
        json.dump(igv_map, f)

    # 2. Results DataFrame
    results_df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [100, 200],
            "REF": ["A", "C"],
            "ALT": ["T", "G"],
            "GT": ["Sample1(0/1)", "Sample1(0/1);Sample2(0/1)"],
            "ClinVar": [
                "https://www.ncbi.nlm.nih.gov/clinvar/?term=chr1-100-A-T",
                "https://www.ncbi.nlm.nih.gov/clinvar/?term=chr2-200-C-G",
            ],
        }
    )

    # 3. Create TSV files for optional sheets
    metadata_tsv = output_dir / "metadata.tsv"
    stats_tsv = output_dir / "statistics.tsv"
    gene_burden_tsv = output_dir / "gene_burden.tsv"
    assoc_tsv = output_dir / "association.tsv"
    tsv_output = output_dir / "variants.tsv"

    results_df.to_csv(tsv_output, sep="\t", index=False)
    pd.DataFrame({"Param": ["version", "build"], "Value": ["1.0", "GRCh38"]}).to_csv(
        metadata_tsv, sep="\t", index=False
    )
    pd.DataFrame({"Metric": ["total", "filtered"], "Count": [10, 2]}).to_csv(
        stats_tsv, sep="\t", index=False
    )
    pd.DataFrame({"Gene": ["GENE1"], "P_value": [0.001]}).to_csv(
        gene_burden_tsv, sep="\t", index=False
    )
    pd.DataFrame({"Gene": ["GENE1"], "Assoc_P": [0.005]}).to_csv(assoc_tsv, sep="\t", index=False)

    workspace = Workspace(output_dir, "variants")
    context = PipelineContext(
        args=None,
        config={
            "xlsx": True,
            "output_file": str(tsv_output),
            "no_stats": False,
            "stats_output_file": str(stats_tsv),
            "perform_gene_burden": True,
            "gene_burden_output": str(gene_burden_tsv),
            "perform_association": True,
            "association_output": str(assoc_tsv),
            "igv_enabled": True,
        },
        workspace=workspace,
    )
    context.final_output_path = tsv_output
    context.current_dataframe = results_df
    context.report_paths["metadata"] = metadata_tsv
    context.mark_complete("dataframe_loading")
    context.mark_complete("tsv_output")
    context.mark_complete("metadata_generation")
    context.mark_complete("statistics_generation")

    return context


@pytest.mark.unit
def test_excel_report_stage_performs_no_openpyxl_reload(excel_context: PipelineContext):
    """ExcelReportStage writes Results and optional sheets without loading via openpyxl."""
    stage = ExcelReportStage()

    with patch("openpyxl.load_workbook", wraps=openpyxl.load_workbook) as spy_load:
        stage(excel_context)
        assert spy_load.call_count == 0, (
            f"Expected 0 load_workbook calls, got {spy_load.call_count}"
        )

    excel_file = excel_context.workspace.output_dir / "variants.xlsx"
    assert excel_file.exists()

    # Now verify the produced file using openpyxl in test verification
    wb = openpyxl.load_workbook(excel_file)
    assert wb.sheetnames == ["Results", "Metadata", "Statistics", "Gene Burden", "Association"]

    for name in wb.sheetnames:
        ws = wb[name]
        assert ws.freeze_panes == "A2", f"Sheet {name} missing freeze_panes A2"
        assert ws.auto_filter.ref is not None, f"Sheet {name} missing autofilter"

    # Verify Results sheet contents, hyperlinks, and IGV Report Links
    results_ws = wb["Results"]
    headers = [cell.value for cell in results_ws[1]]
    assert "ClinVar" in headers
    assert "IGV Report Links" in headers

    clinvar_col = headers.index("ClinVar") + 1
    igv_col = headers.index("IGV Report Links") + 1

    # Row 2 (first variant): ClinVar hyperlink and single IGV link
    cell_clinvar_1 = results_ws.cell(row=2, column=clinvar_col)
    assert cell_clinvar_1.value == "ClinVar"
    expected_clinvar_url = "https://www.ncbi.nlm.nih.gov/clinvar/?term=chr1-100-A-T"
    assert expected_clinvar_url in cell_clinvar_1.hyperlink.target

    cell_igv_1 = results_ws.cell(row=2, column=igv_col)
    assert cell_igv_1.value == "Sample1"
    assert cell_igv_1.hyperlink.target.replace("\\", "/") == "report/igv/Sample1_var1.html"

    # Row 3 (second variant): multiple IGV reports with comment
    cell_igv_2 = results_ws.cell(row=3, column=igv_col)
    assert cell_igv_2.value == "Sample1 (+1 others)"
    assert cell_igv_2.comment is not None
    assert "Sample1 (report/igv/Sample1_var2.html)" in cell_igv_2.comment.text
    assert "Sample2 (report/igv/Sample2_var2.html)" in cell_igv_2.comment.text


@pytest.mark.unit
def test_excel_empty_results_and_association_fidelity(tmp_path: Path):
    """Verify empty Results sheet (0 variants) with Association sheet and formatting."""
    output_dir = tmp_path / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_output = output_dir / "empty.tsv"
    assoc_tsv = output_dir / "association.tsv"

    empty_df = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "GT"])
    empty_df.to_csv(tsv_output, sep="\t", index=False)
    pd.DataFrame({"Gene": ["TP53"], "Assoc_P": [0.01]}).to_csv(assoc_tsv, sep="\t", index=False)

    workspace = Workspace(output_dir, "variants")
    context = PipelineContext(
        args=None,
        config={
            "xlsx": True,
            "output_file": str(tsv_output),
            "no_stats": True,
            "perform_association": True,
            "association_output": str(assoc_tsv),
        },
        workspace=workspace,
    )
    context.final_output_path = tsv_output
    context.current_dataframe = empty_df
    context.mark_complete("dataframe_loading")
    context.mark_complete("tsv_output")
    context.mark_complete("metadata_generation")
    context.mark_complete("statistics_generation")

    stage = ExcelReportStage()
    with patch("openpyxl.load_workbook", wraps=openpyxl.load_workbook) as spy_load:
        stage(context)
        assert spy_load.call_count == 0

    excel_file = context.workspace.get_output_path("", ".xlsx")
    assert excel_file.exists()

    wb = openpyxl.load_workbook(excel_file)
    assert wb.sheetnames == ["Results", "Association"]
    assert wb["Results"].max_row == 1  # Header only
    assert wb["Association"].max_row == 2  # Header + 1 row
