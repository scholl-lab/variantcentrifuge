"""Regression tests for Issue #147:
- Wire --fields-file in CLI and load fields from file
- Fix igv_flanking propagation (especially value 0) into config["igv_flanking"]
- Conditionally add PhenotypeLoadingStage in build_pipeline_stages() when phenotype_file is set
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from variantcentrifuge.cli import create_parser
from variantcentrifuge.cli_config import (
    load_fields_from_file,
    populate_cli_config,
    resolve_presets_and_tumor_normal,
)
from variantcentrifuge.pipeline import build_pipeline_stages, create_stages_from_config
from variantcentrifuge.stages.setup_stages import PhenotypeLoadingStage


@pytest.mark.unit
def test_fields_file_cli_loading(tmp_path: Path):
    """--fields-file reads whitespace/newline-separated fields, ignoring comments."""
    fields_file = tmp_path / "fields.txt"
    fields_file.write_text("# Standard VCF fields\nCHROM POS\n\n# Alleles\nREF, ALT\nQUAL\n")

    parser = create_parser()
    args = parser.parse_args(
        [
            "-v",
            "dummy.vcf.gz",
            "-r",
            "GRCh37.p13",
            "-f",
            "PASS",
            "--fields-file",
            str(fields_file),
        ]
    )

    assert hasattr(args, "fields_file")
    assert args.fields_file == str(fields_file)
    loaded = load_fields_from_file(args.fields_file)
    assert loaded == "CHROM POS REF ALT QUAL"


@pytest.mark.unit
def test_fields_precedence_over_fields_file(tmp_path: Path):
    """Explicit --fields should take precedence over --fields-file."""
    fields_file = tmp_path / "fields.txt"
    fields_file.write_text("CHROM POS\n")

    parser = create_parser()
    args = parser.parse_args(
        [
            "-v",
            "dummy.vcf.gz",
            "-r",
            "GRCh37.p13",
            "-f",
            "PASS",
            "-e",
            "CHROM POS REF ALT",
            "--fields-file",
            str(fields_file),
        ]
    )

    # Both args are set, but args.fields takes precedence when extracting
    assert args.fields == "CHROM POS REF ALT"
    assert args.fields_file == str(fields_file)


@pytest.mark.unit
def test_fields_file_not_found_raises():
    """Missing fields_file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_fields_from_file("nonexistent_fields_file.txt")


@pytest.mark.unit
def test_igv_flanking_zero_preserved(tmp_path: Path):
    """--igv-flanking 0 should propagate as 0, not fall back to default 50."""
    bam_map = tmp_path / "bam_map.tsv"
    bam_map.write_text("sample1\tsample1.bam\n")

    parser = create_parser()
    args = parser.parse_args(
        [
            "-v",
            "dummy.vcf.gz",
            "-r",
            "GRCh37.p13",
            "-f",
            "PASS",
            "-e",
            "CHROM POS",
            "--igv",
            "--bam-mapping-file",
            str(bam_map),
            "--igv-reference",
            "hg19",
            "--igv-flanking",
            "0",
        ]
    )

    cfg: dict = {"reference": "GRCh37.p13", "filters": "PASS", "fields_to_extract": "CHROM POS"}
    filters, tn_vars = resolve_presets_and_tumor_normal(args, cfg, "PASS")
    populate_cli_config(args, cfg, parser, "GRCh37.p13", filters, "CHROM POS", tn_vars)

    assert cfg.get("igv_flanking") == 0


@pytest.mark.unit
def test_phenotype_loading_stage_unconditionally_scheduled():
    """PhenotypeLoadingStage is unconditionally scheduled in all stage plans."""
    # When phenotype_file is set in config
    cfg_with_pheno = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS",
        "phenotype_file": "phenotypes.tsv",
    }
    stages_with = create_stages_from_config(cfg_with_pheno)
    stage_types_with = [type(s) for s in stages_with]
    assert PhenotypeLoadingStage in stage_types_with

    # When phenotype_file is set in args
    args_with = argparse.Namespace(
        reference="GRCh37.p13",
        filters="PASS",
        fields_to_extract="CHROM POS",
        phenotype_file="phenotypes.tsv",
        config={},
    )
    stages_args = build_pipeline_stages(args_with)
    assert PhenotypeLoadingStage in [type(s) for s in stages_args]

    # When phenotype_file is NOT set
    cfg_without_pheno = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS",
    }
    stages_without = create_stages_from_config(cfg_without_pheno)
    stage_types_without = [type(s) for s in stages_without]
    assert PhenotypeLoadingStage in stage_types_without


@pytest.mark.unit
def test_configured_analysis_and_report_stages_activated():
    """Configured analyses and reports must activate corresponding stages in
    create_stages_from_config.
    """
    from variantcentrifuge.stages.analysis_stages import (
        AssociationAnalysisStage,
        GeneBurdenAnalysisStage,
    )
    from variantcentrifuge.stages.output_stages import (
        ArchiveCreationStage,
        ExcelReportStage,
        HTMLReportStage,
    )

    cfg = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS",
        "gene_name": "BRCA1",
        "xlsx": True,
        "html_report": True,
        "archive_results": True,
        "perform_gene_burden": True,
        "perform_association": True,
        "association_tests": ["skat"],
    }
    stages = create_stages_from_config(cfg)
    stage_types = [type(s) for s in stages]

    assert ExcelReportStage in stage_types
    assert HTMLReportStage in stage_types
    assert ArchiveCreationStage in stage_types
    assert GeneBurdenAnalysisStage in stage_types
    assert AssociationAnalysisStage in stage_types


@pytest.mark.unit
def test_threads_propagation_to_stage_builders():
    """threads=8 selects parallel_complete_processing; threads=1 or omitted
    retains sequential stages.
    """
    # threads=8 in create_stages_from_config
    cfg_parallel = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS",
        "threads": 8,
    }
    stages_p1 = create_stages_from_config(cfg_parallel)
    names_p1 = [s.name for s in stages_p1]
    assert "parallel_complete_processing" in names_p1
    assert "variant_extraction" not in names_p1

    # threads=8 in build_pipeline_stages
    args_parallel = argparse.Namespace(
        reference="GRCh37.p13",
        filters="PASS",
        fields_to_extract="CHROM POS",
        threads=8,
        config={"threads": 8},
    )
    stages_p2 = build_pipeline_stages(args_parallel)
    names_p2 = [s.name for s in stages_p2]
    assert "parallel_complete_processing" in names_p2
    assert "variant_extraction" not in names_p2

    # threads=1 retains sequential stages
    cfg_seq = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS",
        "threads": 1,
    }
    stages_s1 = create_stages_from_config(cfg_seq)
    names_s1 = [s.name for s in stages_s1]
    assert "parallel_complete_processing" not in names_s1
    assert "variant_extraction" in names_s1


@pytest.mark.unit
def test_cli_main_mocked_execution_preserves_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Real argument parsing with mocked execution preserves configured
    analysis/report flags and gene selectors.
    """
    import json
    from unittest.mock import MagicMock, patch

    from variantcentrifuge.cli import main

    cfg_data = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS",
        "gene_name": "BRCA1",
        "xlsx": True,
        "html_report": True,
        "archive_results": True,
        "perform_gene_burden": True,
        "perform_association": True,
        "association_tests": ["skat"],
    }
    cfg_file = tmp_path / "config.json"
    cfg_file.write_text(json.dumps(cfg_data))

    mock_run = MagicMock()
    monkeypatch.setattr("variantcentrifuge.cli.run_pipeline", mock_run)
    monkeypatch.setattr("variantcentrifuge.cli.validate_vcf_file", lambda *a, **k: None)
    monkeypatch.setattr("variantcentrifuge.cli.validate_required_tools", lambda *a, **k: None)

    with patch("sys.argv", ["variantcentrifuge", "-v", "test.vcf.gz", "-c", str(cfg_file)]):
        ret = main()
        assert ret == 0

    assert mock_run.called
    called_args = mock_run.call_args[0][0]
    assert called_args.gene_name == "BRCA1"
    assert called_args.xlsx is True
    assert called_args.html_report is True
    assert called_args.archive_results is True
    assert called_args.perform_gene_burden is True
    assert called_args.perform_association is True
    assert called_args.association_tests == ["skat"]


@pytest.mark.unit
def test_cli_explicit_flags_override_config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Explicit CLI flags override configuration file settings."""
    import json
    from unittest.mock import MagicMock, patch

    from variantcentrifuge.cli import main

    cfg_data = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS",
        "gene_name": "BRCA1",
    }
    cfg_file = tmp_path / "config.json"
    cfg_file.write_text(json.dumps(cfg_data))

    mock_run = MagicMock()
    monkeypatch.setattr("variantcentrifuge.cli.run_pipeline", mock_run)
    monkeypatch.setattr("variantcentrifuge.cli.validate_vcf_file", lambda *a, **k: None)
    monkeypatch.setattr("variantcentrifuge.cli.validate_required_tools", lambda *a, **k: None)

    with patch(
        "sys.argv",
        ["variantcentrifuge", "-v", "test.vcf.gz", "-c", str(cfg_file), "--gene-name", "TP53"],
    ):
        ret = main()
        assert ret == 0

    assert mock_run.called
    called_args = mock_run.call_args[0][0]
    assert called_args.gene_name == "TP53"


@pytest.mark.unit
def test_config_only_annotate_json_and_igv_pass_validation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Config-only annotate_json_genes and igv pass validation cleanly without CLI flags."""
    import json
    from unittest.mock import MagicMock, patch

    from variantcentrifuge.cli import main

    json_mapping_file = tmp_path / "mapping.json"
    json_mapping_file.write_text('{"GENE": "annotation"}')
    bam_mapping_file = tmp_path / "bams.tsv"
    bam_mapping_file.write_text("sample\t/path/to/sample.bam\n")

    cfg_data = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS",
        "gene_name": "BRCA1",
        "annotate_json_genes": [str(json_mapping_file)],
        "json_gene_mapping": str(json_mapping_file),
        "igv": True,
        "bam_mapping_file": str(bam_mapping_file),
        "igv_reference": "hg19",
    }
    cfg_file = tmp_path / "config.json"
    cfg_file.write_text(json.dumps(cfg_data))

    mock_run = MagicMock()
    monkeypatch.setattr("variantcentrifuge.cli.run_pipeline", mock_run)
    monkeypatch.setattr("variantcentrifuge.cli.validate_vcf_file", lambda *a, **k: None)
    monkeypatch.setattr("variantcentrifuge.cli.validate_required_tools", lambda *a, **k: None)

    with patch(
        "sys.argv",
        ["variantcentrifuge", "-v", "test.vcf.gz", "-c", str(cfg_file)],
    ):
        ret = main()
        assert ret == 0

    assert mock_run.called
    args = mock_run.call_args[0][0]
    assert args.igv is True
    assert args.bam_mapping_file == str(bam_mapping_file)
    assert args.igv_reference == "hg19"
    assert args.annotate_json_genes == [str(json_mapping_file)]
    assert args.json_gene_mapping == str(json_mapping_file)
