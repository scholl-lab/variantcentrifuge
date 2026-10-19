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
            "CHROM": ["chr1", "chr2", "chr3"],
            "POS": [100, 200, 300],
            "REF": ["A", "C", "G"],
            "ALT": ["T", "G", "A"],
            "GT": ["Sample1(0/1)", "Sample1(0/1);Sample2(0/1)", "Sample3(0/1)"],
            "ClinVar": [
                "https://www.ncbi.nlm.nih.gov/clinvar/?term=chr1-100-A-T",
                "https://www.ncbi.nlm.nih.gov/clinvar/?term=chr2-200-C-G",
                "https://www.ncbi.nlm.nih.gov/clinvar/?term=chr3-300-G-A",
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

    # Row 4 (third variant): no IGV report mapped -> N/A
    cell_igv_3 = results_ws.cell(row=4, column=igv_col)
    assert cell_igv_3.value == "N/A"

    # Verify optional sheets cell values
    meta_ws = wb["Metadata"]
    assert [c.value for c in meta_ws[1]] == ["Param", "Value"]
    assert [c.value for c in meta_ws[2]] == ["version", "1.0"]
    assert [c.value for c in meta_ws[3]] == ["build", "GRCh38"]

    stats_ws = wb["Statistics"]
    assert [c.value for c in stats_ws[1]] == ["Metric", "Count"]
    assert [c.value for c in stats_ws[2]] == ["total", 10]
    assert [c.value for c in stats_ws[3]] == ["filtered", 2]

    burden_ws = wb["Gene Burden"]
    assert [c.value for c in burden_ws[1]] == ["Gene", "P_value"]
    assert [c.value for c in burden_ws[2]] == ["GENE1", 0.001]

    assoc_ws = wb["Association"]
    assert [c.value for c in assoc_ws[1]] == ["Gene", "Assoc_P"]
    assert [c.value for c in assoc_ws[2]] == ["GENE1", 0.005]


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


@pytest.mark.unit
def test_excel_legacy_list_format_igv_map(tmp_path: Path):
    """Verify convert_to_excel supports top-level list IGV maps (#159)."""
    from variantcentrifuge.converter import convert_to_excel

    output_dir = tmp_path / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir = output_dir / "report" / "igv"
    report_dir.mkdir(parents=True, exist_ok=True)

    # Legacy top-level list format with flat report_path
    legacy_map = [
        {
            "chrom": "chr1",
            "pos": 100,
            "ref": "A",
            "alt": "T",
            "sample_id": "Sample1",
            "report_path": "igv/Sample1_legacy.html",
        }
    ]
    with open(report_dir / "igv_reports_map.json", "w", encoding="utf-8") as f:
        json.dump(legacy_map, f)

    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["T"],
            "GT": ["Sample1(0/1)"],
        }
    )
    xlsx_file = output_dir / "variants.xlsx"
    convert_to_excel(
        str(xlsx_file),
        cfg={"igv_enabled": True},
        df=df,
        output_file=str(xlsx_file),
    )

    wb = openpyxl.load_workbook(xlsx_file)
    ws = wb["Results"]
    headers = [cell.value for cell in ws[1]]
    assert "IGV Report Links" in headers
    igv_cell = ws.cell(row=2, column=headers.index("IGV Report Links") + 1)
    assert igv_cell.value == "Sample1"
    assert "igv/Sample1_legacy.html" in igv_cell.hyperlink.target.replace("\\", "/")


@pytest.mark.unit
def test_excel_6000_rows_6_urls_no_hyperlink_limit_exhaustion(tmp_path: Path):
    """Verify 6,000 rows x 6 URL columns (36,000 URLs) do not hit XlsxWriter 65,530 limit (#159)."""
    from variantcentrifuge.converter import convert_to_excel

    output_dir = tmp_path / "stress_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    n_rows = 6000
    data = {
        "CHROM": ["chr1"] * n_rows,
        "POS": list(range(1, n_rows + 1)),
        "REF": ["A"] * n_rows,
        "ALT": ["T"] * n_rows,
        "GT": ["Sample1(0/1)"] * n_rows,
    }
    for col_idx in range(1, 7):
        data[f"URL_COL_{col_idx}"] = [
            f"https://example.org/path/{col_idx}/{row}" for row in range(n_rows)
        ]

    df = pd.DataFrame(data)
    xlsx_file = output_dir / "large_urls.xlsx"

    # With strings_to_urls=False, exactly 36,000 hyperlinks are written (under 65,530)
    convert_to_excel(
        str(xlsx_file),
        cfg={"igv_enabled": False},
        df=df,
        output_file=str(xlsx_file),
    )

    assert xlsx_file.exists()
    wb = openpyxl.load_workbook(xlsx_file)
    ws = wb["Results"]
    assert ws.max_row == n_rows + 1

    # Verify first and last rows have proper hyperlinks
    first_url_cell = ws.cell(row=2, column=6)
    assert first_url_cell.value == "URL_COL_1"
    assert first_url_cell.hyperlink.target == "https://example.org/path/1/0"

    last_url_cell = ws.cell(row=n_rows + 1, column=11)
    assert last_url_cell.value == "URL_COL_6"
    assert last_url_cell.hyperlink.target == f"https://example.org/path/6/{n_rows - 1}"


@pytest.mark.unit
def test_excel_generation_benchmark_performance_and_memory(tmp_path: Path):
    """Deterministic structural/memory comparison: multi-pass vs single-pass writer (#159)."""
    import tracemalloc

    from variantcentrifuge.converter import (
        append_tsv_as_sheet,
        convert_to_excel,
        finalize_excel_file,
    )

    n_rows = 1000
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"] * n_rows,
            "POS": list(range(n_rows)),
            "REF": ["A"] * n_rows,
            "ALT": ["T"] * n_rows,
            "GT": ["Sample1(0/1)"] * n_rows,
            "ClinVar": [f"https://example.org/cv/{i}" for i in range(n_rows)],
        }
    )
    meta_tsv = tmp_path / "metadata.tsv"
    stats_tsv = tmp_path / "statistics.tsv"
    pd.DataFrame({"K": ["version"], "V": ["1.0"]}).to_csv(meta_tsv, sep="\t", index=False)
    pd.DataFrame({"Metric": ["count"], "Value": [n_rows]}).to_csv(stats_tsv, sep="\t", index=False)

    extra_sheets = {
        "Metadata": pd.DataFrame({"K": ["version"], "V": ["1.0"]}),
        "Statistics": pd.DataFrame({"Metric": ["count"], "Value": [n_rows]}),
    }

    # 1. Multi-pass reload baseline (historical behavior prior to #159)
    f_baseline = tmp_path / "baseline.xlsx"
    tracemalloc.start()
    df.to_excel(f_baseline, sheet_name="Results", index=False)
    append_tsv_as_sheet(str(f_baseline), str(meta_tsv), sheet_name="Metadata")
    append_tsv_as_sheet(str(f_baseline), str(stats_tsv), sheet_name="Statistics")
    finalize_excel_file(str(f_baseline), cfg={"igv_enabled": False})
    _current, peak_base = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # 2. Single-pass xlsxwriter generator (#159)
    f_single = tmp_path / "single.xlsx"
    tracemalloc.start()
    convert_to_excel(
        str(f_single),
        cfg={"igv_enabled": False},
        df=df,
        extra_sheets=extra_sheets,
        output_file=str(f_single),
    )
    _current, peak_single = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    peak_base_mb = peak_base / (1024 * 1024)
    peak_single_mb = peak_single / (1024 * 1024)

    # Assert deterministic structural and memory improvement without flaky wall-clock caps
    assert peak_single < peak_base, (
        f"Single-pass peak ({peak_single_mb:.2f}MB) should be lower than baseline "
        f"({peak_base_mb:.2f}MB)"
    )
    assert peak_single_mb < 50.0, f"Peak memory {peak_single_mb:.1f}MB, expected < 50MB"

    wb_single = openpyxl.load_workbook(f_single)
    wb_base = openpyxl.load_workbook(f_baseline)
    assert wb_single.sheetnames == ["Results", "Metadata", "Statistics"]
    assert wb_base.sheetnames == ["Results", "Metadata", "Statistics"]


@pytest.mark.unit
def test_finalize_excel_file_missing_coordinates_guard(tmp_path: Path):
    """Verify finalize_excel_file handles missing coordinates without column 0 ValueError (#159)."""
    from variantcentrifuge.converter import finalize_excel_file

    output_dir = tmp_path / "missing_coords"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir = output_dir / "report" / "igv"
    report_dir.mkdir(parents=True, exist_ok=True)

    igv_map = {
        "variants": [
            {
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "T",
                "sample_reports": {"Sample1": "igv/Sample1.html"},
            }
        ]
    }
    with open(report_dir / "igv_reports_map.json", "w", encoding="utf-8") as f:
        json.dump(igv_map, f)

    # DataFrame with only CHROM and POS (missing REF and ALT)
    df = pd.DataFrame({"CHROM": ["chr1"], "POS": [100]})
    xlsx_file = output_dir / "partial.xlsx"
    df.to_excel(xlsx_file, sheet_name="Results", index=False)

    # Should safely skip IGV links rather than calling ws.cell(..., column=0)
    finalize_excel_file(str(xlsx_file), cfg={"igv_enabled": True})

    wb = openpyxl.load_workbook(xlsx_file)
    ws = wb["Results"]
    headers = [cell.value for cell in ws[1]]
    assert "IGV Report Links" not in headers
    assert headers == ["CHROM", "POS"]


@pytest.mark.unit
def test_excel_igv_comment_overflow_creates_overflow_sheet(tmp_path: Path):
    """Verify IGV comment >32,000 chars is truncated and full links saved to IGV_Overflow (#159)."""
    from variantcentrifuge.converter import convert_to_excel

    output_dir = tmp_path / "overflow_test"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir = output_dir / "report" / "igv"
    report_dir.mkdir(parents=True, exist_ok=True)

    # 600 samples with long paths -> comment will exceed 32,000 characters
    n_samples = 600
    subfolder = "igv/long_batch_path_subfolder"
    sample_reports = {
        f"Patient_Sample_ID_{i:04d}": f"{subfolder}/Patient_Sample_ID_{i:04d}_variant_report.html"
        for i in range(n_samples)
    }
    igv_map = {
        "variants": [
            {
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "T",
                "sample_reports": sample_reports,
            }
        ]
    }
    with open(report_dir / "igv_reports_map.json", "w", encoding="utf-8") as f:
        json.dump(igv_map, f)

    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["T"],
            "GT": [";".join(f"Patient_Sample_ID_{i:04d}(0/1)" for i in range(n_samples))],
        }
    )
    xlsx_file = output_dir / "overflow.xlsx"
    convert_to_excel(
        str(xlsx_file),
        cfg={"igv_enabled": True},
        df=df,
        output_file=str(xlsx_file),
    )

    wb = openpyxl.load_workbook(xlsx_file)
    assert "Results" in wb.sheetnames
    assert "IGV_Overflow" in wb.sheetnames

    ws_results = wb["Results"]
    headers = [cell.value for cell in ws_results[1]]
    assert "IGV Report Links" in headers
    igv_col = headers.index("IGV Report Links") + 1

    cell = ws_results.cell(row=2, column=igv_col)
    assert cell.value == f"Patient_Sample_ID_0000 (+{n_samples - 1} others)"
    assert cell.comment is not None

    comment_text = cell.comment.text
    # Cell comment must be at most 32,767 chars (Excel limit) and end with overflow suffix
    assert len(comment_text) <= 32767
    assert len(comment_text) <= 32000
    expected_suffix = "\n... (see IGV_Overflow sheet for all sample links)"
    assert comment_text.endswith(expected_suffix)

    # Verify IGV_Overflow sheet
    ws_overflow = wb["IGV_Overflow"]
    overflow_headers = [c.value for c in ws_overflow[1]]
    expected_headers = ["Row", "CHROM", "POS", "REF", "ALT", "Sample_ID", "IGV_Report_Path"]
    assert overflow_headers == expected_headers

    # 1 header row + 600 sample rows
    assert ws_overflow.max_row == n_samples + 1
    first_row = [
        ws_overflow.cell(row=2, column=i).value for i in range(1, len(expected_headers) + 1)
    ]
    assert first_row == [
        2,
        "chr1",
        100,
        "A",
        "T",
        "Patient_Sample_ID_0000",
        "report/igv/long_batch_path_subfolder/Patient_Sample_ID_0000_variant_report.html",
    ]
