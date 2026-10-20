"""Tests for pseudonymized metadata and Excel report exports (#128)."""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pandas as pd

from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.metadata_helpers import collect_metadata_rows
from variantcentrifuge.stages.output_stages import (
    ExcelReportStage,
    MetadataGenerationStage,
    PseudonymizationStage,
)


def dummy_version(tool: str) -> str:
    return f"{tool}-1.0.0"


def test_pseudonymized_metadata_tsv_and_excel_omits_reidentification(tmp_path: Path):
    """Verify original sample IDs and mapping entries are absent from metadata."""
    workspace = Workspace(tmp_path / "work", "test")
    orig_samples = ["PATIENT_001", "CONTROL_999", "MOTHER_123"]
    mapping = {
        "PATIENT_001": "SAMPLE_A",
        "CONTROL_999": "SAMPLE_B",
        "MOTHER_123": "SAMPLE_C",
    }

    config = {
        "pseudonymize": True,
        "pseudonymize_schema": "sequential",
        "pseudonymize_prefix": "SAMPLE",
        "pseudonymize_table": str(tmp_path / "private_mapping.tsv"),
        "xlsx": True,
        "sample_list": ",".join(orig_samples),
        "case_samples": ["PATIENT_001"],
        "control_samples": ["CONTROL_999"],
        "proband_samples": ["PATIENT_001"],
        "vcf_file": "/data/cohorts/PATIENT_001_trio.vcf.gz",
        "ped_file": "/data/pedigrees/PATIENT_001_family.ped",
        "phenotype_file": "/data/pheno/PATIENT_001_pheno.tsv",
        "filter": "QUAL > 30 && (GEN[PATIENT_001].GT == '0/1')",
        "final_filter": "DP > 10 && PATIENT_001 != '0/0'",
        "tumor_sample": "PATIENT_001",
        "normal_sample": "CONTROL_999",
        "nested_sample_info": {"sample": "PATIENT_001", "status": "affected"},
        "reference": "GRCh38",
        "threads": 4,
    }

    context = PipelineContext(Namespace(), config, workspace)
    context.vcf_samples = list(orig_samples)
    context.current_dataframe = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["BRCA1"],
            "PATIENT_001": ["0/1"],
            "CONTROL_999": ["0/0"],
            "MOTHER_123": ["0/0"],
        }
    )

    # Run PseudonymizationStage to verify private mapping is created and preserved
    pseudo_stage = PseudonymizationStage()
    pseudo_stage._process(context)
    mapping_file = Path(context.config["pseudonymize_table"])
    assert mapping_file.exists(), "Private mapping file must be saved"
    mapping_content = mapping_file.read_text()
    for orig_id in orig_samples:
        assert orig_id in mapping_content, "Private mapping must contain original identifiers"

    mapping = context.config["pseudonymization_mapping"]
    assert mapping is not None

    rows = collect_metadata_rows(context, dummy_version)
    metadata_df = pd.DataFrame(rows, columns=["Parameter", "Value"])

    # Check that NO original sample identifier appears anywhere in the metadata
    for orig_id in orig_samples:
        for val in metadata_df["Value"].astype(str):
            assert orig_id not in val, f"Leaked original sample ID {orig_id} in value: {val}"
        for param in metadata_df["Parameter"].astype(str):
            assert orig_id not in param, f"Leaked ID {orig_id} in parameter: {param}"

    # Verify that the complete lookup mapping is not stored
    assert not any("pseudonymization_mapping" in p for p in metadata_df["Parameter"])
    assert not any("nested_sample_info" in p for p in metadata_df["Parameter"])
    assert not any("sample_list" in p for p in metadata_df["Parameter"])
    assert not any("case_samples" in p for p in metadata_df["Parameter"])
    assert not any("control_samples" in p for p in metadata_df["Parameter"])

    # Verify that tumor/normal has been sanitized with pseudonymized identifiers
    tn_row = metadata_df[metadata_df["Parameter"] == "Tumor/Normal Samples"]
    assert not tn_row.empty
    tn_val = tn_row.iloc[0]["Value"]
    assert f"tumor={mapping['PATIENT_001']}" in tn_val
    assert f"normal={mapping['CONTROL_999']}" in tn_val

    # Verify that filters have been sanitized
    filt_row = metadata_df[metadata_df["Parameter"] == "Filter"]
    assert f"GEN[{mapping['PATIENT_001']}].GT" in filt_row.iloc[0]["Value"]

    # Verify that metadata generation stage creates TSV and Excel reports without leaks
    meta_stage = MetadataGenerationStage()
    meta_stage._process(context)
    meta_tsv_path = context.report_paths["metadata"]
    assert meta_tsv_path.exists()
    tsv_content = meta_tsv_path.read_text()
    for orig_id in orig_samples:
        assert orig_id not in tsv_content

    # Now verify ExcelReportStage with Metadata sheet
    out_tsv = tmp_path / "report.tsv"
    out_tsv.write_text(
        f"CHROM\tPOS\tREF\tALT\tGENE\t{mapping['PATIENT_001']}\t{mapping['CONTROL_999']}\t{mapping['MOTHER_123']}\n"
        "chr1\t100\tA\tT\tBRCA1\t0/1\t0/0\t0/0\n"
    )
    context.final_output_path = out_tsv
    context.config["output_file"] = str(out_tsv)
    excel_stage = ExcelReportStage()
    excel_stage._process(context)
    excel_path = context.report_paths.get("excel")
    assert excel_path is not None and excel_path.exists(), "Excel report must be generated"
    excel_sheets = pd.read_excel(excel_path, sheet_name=None)
    assert "Metadata" in excel_sheets, "Metadata sheet must be present in Excel report"
    excel_meta = excel_sheets["Metadata"]
    excel_text = excel_meta.to_csv(index=False)
    for orig_id in orig_samples:
        assert orig_id not in excel_text

    # Verify internal config was not mutated
    assert context.config["sample_list"] == ",".join(orig_samples)
    assert context.config["case_samples"] == ["PATIENT_001"]
    assert context.config["pseudonymization_mapping"] == mapping


def test_non_pseudonymized_metadata_remains_compatible(tmp_path: Path):
    """Verify that metadata behavior without --pseudonymize remains unchanged and compatible."""
    workspace = Workspace(tmp_path / "work", "test")
    config = {
        "pseudonymize": False,
        "vcf_file": "/data/test.vcf",
        "sample_list": "S1,S2",
        "vcf_samples": ["S1", "S2"],
        "reference": "GRCh38",
        "threads": 2,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.vcf_samples = ["S1", "S2"]
    rows = collect_metadata_rows(context, dummy_version)
    metadata_df = pd.DataFrame(rows, columns=["Parameter", "Value"])

    # In un-pseudonymized mode, config dumps are included and sample_list/vcf_samples preserved
    params = metadata_df["Parameter"].tolist()
    assert "config.threads" in params
    assert "config.reference" in params
    assert "config.sample_list" in params
    assert "config.vcf_samples" in params


def test_pseudonymized_metadata_with_empty_mapping(tmp_path: Path):
    """When pseudonymize=True but mapping is empty/unavailable, identifying strings are redacted."""
    workspace = Workspace(tmp_path / "work", "test")
    config = {
        "pseudonymize": True,
        "vcf_file": "/data/PRIVATE_PATIENT_001.vcf.gz",
        "reference": "GRCh38",
        "filter": "SAMPLE == 'PRIVATE_001'",
        "final_filter": "DP > 10",
        "scoring_config_path": "/configs/scoring.json",
        "phenotype_file": "/data/cases.tsv",
        "stats_config": "/data/PATIENT_001/stats.json",
    }
    context = PipelineContext(Namespace(), config, workspace)
    rows = collect_metadata_rows(context, dummy_version)
    metadata_df = pd.DataFrame(rows, columns=["Parameter", "Value"]).set_index("Parameter")

    assert metadata_df.loc["VCF_file", "Value"] == "<redacted>"
    assert metadata_df.loc["Reference", "Value"] == "<redacted>"
    assert metadata_df.loc["Filter", "Value"] == "<redacted>"
    assert metadata_df.loc["Final Filter", "Value"] == "<redacted>"
    assert metadata_df.loc["Scoring_config", "Value"] == "<redacted>"
    assert metadata_df.loc["Phenotype_file", "Value"] == "<redacted>"
    assert metadata_df.loc["config.reference", "Value"] == "<redacted>"
    assert metadata_df.loc["config.stats_config", "Value"] == "<redacted>"


def test_short_numeric_sample_ids_and_directory_stripping(tmp_path: Path):
    """Short numeric sample IDs ('1', '2') do not corrupt GRCh38.p13 or numeric fields."""
    from variantcentrifuge.stages.metadata_helpers import collect_metadata_rows

    workspace = Workspace(output_dir=tmp_path, base_name="test")
    config = {
        "pseudonymize": True,
        "pseudonymization_mapping": {"1": "STUDY_001", "2": "STUDY_002"},
        "vcf_file": "/secret/patient_dir/cohort_1.vcf.gz",
        "ped_file": "/secret/patient_dir/family_2.ped",
        "reference": "GRCh38.p13",
        "filter": "GEN[1].GT != '0/0'",
        "threads": 1,
        "pipeline_version": "0.15.0",
    }
    context = PipelineContext(Namespace(), config, workspace)
    rows = dict(collect_metadata_rows(context, lambda t: "1.2.3"))

    assert rows["Reference"] == "GRCh38.p13"
    assert rows["Version"] == "0.15.0"
    assert rows["Threads"] == 1
    assert rows["VCF_file"] == "cohort_STUDY_001.vcf.gz"
    assert "/secret/patient_dir" not in rows["VCF_file"]
    assert rows["Pedigree_file"] == "family_STUDY_002.ped"
    assert rows["Filter"] == "GEN[STUDY_001].GT != '0/0'"


def test_generate_html_report_sanitizes_metadata_when_pseudonymized(tmp_path: Path):
    """generate_html_report sanitizes vcf_source, filter_expression, and summary.json metadata."""
    import json

    from variantcentrifuge.generate_html_report import generate_html_report

    variants_path = tmp_path / "variants.json"
    summary_path = tmp_path / "summary.json"
    variants_path.write_text(
        json.dumps([{"CHROM": "chr1", "POS": 100, "REF": "A", "ALT": "T", "GENE": "BRCA1"}]),
        encoding="utf-8",
    )
    summary_path.write_text(
        json.dumps(
            {
                "num_variants": 1,
                "num_genes": 1,
                "impact_distribution": {},
                "metadata": {
                    "vcf_file": "cohort_PATIENT_001.vcf.gz",
                    "filter": "GEN[PATIENT_001].GT != '0/0'",
                },
            }
        ),
        encoding="utf-8",
    )

    cfg = {
        "pseudonymize": True,
        "pseudonymization_mapping": {"PATIENT_001": "STUDY_001"},
        "vcf_file": "/secret/PATIENT_001/cohort_PATIENT_001.vcf.gz",
        "filter": "GEN[PATIENT_001].GT != '0/0'",
        "reference": "GRCh38.p13",
    }
    generate_html_report(str(variants_path), str(summary_path), str(tmp_path), cfg)

    html = (tmp_path / "index.html").read_text(encoding="utf-8")
    assert "PATIENT_001" not in html
    assert "cohort_STUDY_001.vcf.gz" in html
    assert "GEN[STUDY_001].GT" in html

    updated_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert "PATIENT_001" not in json.dumps(updated_summary)
    assert updated_summary["metadata"]["vcf_file"] == "cohort_STUDY_001.vcf.gz"


def test_pseudonymizer_replaces_sample_id_in_inheritance_details():
    """SamplePseudonymizer.pseudonymize_dataframe replaces sample_id inside Inheritance_Details JSON."""
    from variantcentrifuge.pseudonymizer import SamplePseudonymizer, SequentialSchema

    pseudonymizer = SamplePseudonymizer(SequentialSchema(prefix="STUDY"))
    pseudonymizer.create_mapping(["PATIENT_001", "PATIENT_002"])

    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "Inheritance_Details": [
                '{"primary_pattern": "de_novo", "samples_with_pattern": '
                '[{"sample_id": "PATIENT_001", "genotype": "0/1"}, '
                '{"sample_id": "PATIENT_002", "genotype": "0/0"}]}'
            ],
        }
    )
    out = pseudonymizer.pseudonymize_dataframe(df)
    details_str = out.loc[0, "Inheritance_Details"]
    assert "PATIENT_001" not in details_str
    assert "PATIENT_002" not in details_str
    assert "STUDY_001" in details_str
    assert "STUDY_002" in details_str
