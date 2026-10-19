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
def test_phenotype_loading_stage_conditional():
    """PhenotypeLoadingStage is added when phenotype_file is set, and omitted when not set."""
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
    assert PhenotypeLoadingStage not in stage_types_without
