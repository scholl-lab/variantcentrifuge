"""Regression tests for Issue #146: Load supplied PED before CLI inheritance analysis
and preserve JSON config values when CLI flags are omitted.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from variantcentrifuge.cli import create_parser
from variantcentrifuge.cli_config import populate_cli_config, resolve_presets_and_tumor_normal
from variantcentrifuge.pipeline import build_pipeline_stages, create_stages_from_config
from variantcentrifuge.stages.setup_stages import PedigreeLoadingStage


@pytest.mark.unit
def test_cli_ped_schedules_pedigree_loading_stage(tmp_path: Path):
    """Passing --ped on CLI must schedule PedigreeLoadingStage."""
    ped_file = tmp_path / "family.ped"
    ped_file.write_text("fam1\tchild\tfather\tmother\t1\t2\n")

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
            "--ped",
            str(ped_file),
        ]
    )

    cfg: dict = {"reference": "GRCh37.p13", "filters": "PASS", "fields_to_extract": "CHROM POS"}
    filters, tn_vars = resolve_presets_and_tumor_normal(args, cfg, "PASS")
    populate_cli_config(args, cfg, parser, "GRCh37.p13", filters, "CHROM POS", tn_vars)

    # In cli.py, refactored_args is built from args and cfg
    refactored_args = argparse.Namespace(**vars(args))
    refactored_args.config = cfg
    refactored_args.ped_file = cfg.get("ped_file")

    stages = build_pipeline_stages(refactored_args)
    stage_types = [type(s) for s in stages]
    assert PedigreeLoadingStage in stage_types, (
        "PedigreeLoadingStage should be scheduled with --ped"
    )


@pytest.mark.unit
def test_config_ped_file_schedules_pedigree_loading_stage(tmp_path: Path):
    """Providing ped_file in config dict must schedule PedigreeLoadingStage."""
    ped_file = tmp_path / "family.ped"
    ped_file.write_text("fam1\tchild\tfather\tmother\t1\t2\n")

    cfg = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS",
        "ped_file": str(ped_file),
    }
    stages = create_stages_from_config(cfg)
    stage_types = [type(s) for s in stages]
    assert PedigreeLoadingStage in stage_types, (
        "PedigreeLoadingStage should be scheduled when ped_file is in config"
    )


@pytest.mark.unit
def test_no_ped_omits_pedigree_loading_stage():
    """When no PED is provided, PedigreeLoadingStage is omitted."""
    cfg = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS",
    }
    stages = create_stages_from_config(cfg)
    stage_types = [type(s) for s in stages]
    assert PedigreeLoadingStage not in stage_types, (
        "PedigreeLoadingStage should NOT be scheduled without PED"
    )


@pytest.mark.unit
def test_json_config_values_preserved_when_cli_flags_omitted(tmp_path: Path):
    """JSON config values (xlsx, keep_intermediates, archive_results, sort_memory_limit,
    presets, fields_to_extract, Late_filtering, final_filter) must be preserved
    when CLI flags are omitted.
    """
    custom_cfg = {
        "reference": "GRCh37.p13",
        "filters": "PASS",
        "fields_to_extract": "CHROM POS REF ALT",
        "xlsx": True,
        "keep_intermediates": True,
        "archive_results": True,
        "sort_memory_limit": "8G",
        "presets": ["pass"],
        "Late_filtering": True,
        "final_filter": "QUAL > 30",
    }

    parser = create_parser()
    # CLI passes only required args, omitting flags for xlsx, keep_intermediates, etc.
    args = parser.parse_args(["-v", "dummy.vcf.gz", "-r", "GRCh37.p13"])

    cfg = custom_cfg.copy()
    reference = args.reference or cfg.get("reference")
    filters = args.filters or cfg.get("filters")
    fields = args.fields or cfg.get("fields_to_extract")
    filters, tn_vars = resolve_presets_and_tumor_normal(args, cfg, filters)
    populate_cli_config(args, cfg, parser, reference, filters, fields, tn_vars)

    # Mimic cli.py post-parse config updates
    if args.xlsx:
        cfg["xlsx"] = args.xlsx
    elif "xlsx" in custom_cfg:
        cfg["xlsx"] = custom_cfg["xlsx"]

    if args.keep_intermediates:
        cfg["keep_intermediates"] = args.keep_intermediates
    elif "keep_intermediates" in custom_cfg:
        cfg["keep_intermediates"] = custom_cfg["keep_intermediates"]

    if args.archive_results:
        cfg["archive_results"] = args.archive_results
    elif "archive_results" in custom_cfg:
        cfg["archive_results"] = custom_cfg["archive_results"]

    # Verify preserved values
    assert cfg.get("xlsx") is True
    assert cfg.get("keep_intermediates") is True
    assert cfg.get("archive_results") is True
    assert cfg.get("sort_memory_limit") == "8G"
    assert "pass" in cfg.get("preset", [])
    assert cfg.get("late_filtering") is True
    assert "QUAL > 30" in str(cfg.get("final_filter"))
    assert "CHROM POS REF ALT" in str(cfg.get("fields_to_extract"))
