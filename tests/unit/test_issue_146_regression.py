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


@pytest.mark.unit
def test_check_scoring_requires_inheritance(tmp_path: Path):
    """Test check_scoring_requires_inheritance with various scoring configs."""
    from variantcentrifuge.scoring_requirements import check_scoring_requires_inheritance

    # 1. No scoring config -> False
    args = argparse.Namespace()
    assert check_scoring_requires_inheritance(args, {}) is False

    # 2. Scoring config with formula requiring Inheritance_Pattern -> True
    score_dir = tmp_path / "model1"
    score_dir.mkdir()
    (score_dir / "variable_assignment_config.json").write_text("{}")
    (score_dir / "formula_config.json").write_text(
        '{"formulas": [{"test_score": "1.0 if Inheritance_Pattern == \'de_novo\' else 0.0"}]}'
    )
    assert (
        check_scoring_requires_inheritance(
            argparse.Namespace(), {"scoring_config_path": str(score_dir)}
        )
        is True
    )

    # 3. Scoring config with variables containing Inheritance_Details -> True
    score_dir2 = tmp_path / "model2"
    score_dir2.mkdir()
    (score_dir2 / "variable_assignment_config.json").write_text(
        '{"variables": ["Inheritance_Details"]}'
    )
    (score_dir2 / "formula_config.json").write_text('{"formulas": []}')
    assert (
        check_scoring_requires_inheritance(
            argparse.Namespace(), {"scoring_config_path": str(score_dir2)}
        )
        is True
    )

    # 4. Invalid path -> False
    assert (
        check_scoring_requires_inheritance(
            argparse.Namespace(), {"scoring_config_path": str(tmp_path / "nonexistent")}
        )
        is False
    )


@pytest.mark.unit
def test_synthetic_trio_inheritance_classification_with_cli_ped(tmp_path: Path):
    """Synthetic trio passed through CLI configuration and pedigree loading must
    retain parental relationships and produce de novo and compound-heterozygous
    classifications without unrelated-sample fallback.
    """
    from unittest.mock import Mock

    import pandas as pd

    from variantcentrifuge.inheritance.analyzer import analyze_inheritance
    from variantcentrifuge.pipeline_core import PipelineContext

    # Create synthetic trio PED
    ped_file = tmp_path / "trio.ped"
    ped_file.write_text(
        "fam1\tchild\tfather\tmother\t1\t2\nfam1\tfather\t0\t0\t1\t1\nfam1\tmother\t0\t0\t2\t1\n"
    )

    # CLI parse with --ped
    parser = create_parser()
    args = parser.parse_args(["-v", "dummy.vcf.gz", "--ped", str(ped_file)])

    cfg: dict = {}
    filters, tn_vars = resolve_presets_and_tumor_normal(args, cfg, None)
    populate_cli_config(args, cfg, parser, "GRCh37", filters, "CHROM POS", tn_vars)

    # Execute PedigreeLoadingStage
    context = PipelineContext(args=args, config=cfg, workspace=Mock())
    ped_stage = PedigreeLoadingStage()
    context = ped_stage._process(context)

    assert context.pedigree_data is not None
    assert "child" in context.pedigree_data
    assert context.pedigree_data["child"]["father_id"] == "father"
    assert context.pedigree_data["child"]["mother_id"] == "mother"

    # Synthetic variants dataframe:
    # Var 1: De novo (child 0/1, father 0/0, mother 0/0)
    # Var 2: Comp het 1 in GENE2 (child 0/1, father 0/1, mother 0/0)
    # Var 3: Comp het 2 in GENE2 (child 0/1, father 0/0, mother 0/1)
    df = pd.DataFrame(
        [
            {
                "CHROM": "1",
                "POS": 1000,
                "REF": "A",
                "ALT": "G",
                "GENE": "GENE1",
                "father": "0/0",
                "mother": "0/0",
                "child": "0/1",
                "GT_father": "0/0",
                "GT_mother": "0/0",
                "GT_child": "0/1",
            },
            {
                "CHROM": "2",
                "POS": 2000,
                "REF": "C",
                "ALT": "T",
                "GENE": "GENE2",
                "father": "0/1",
                "mother": "0/0",
                "child": "0/1",
                "GT_father": "0/1",
                "GT_mother": "0/0",
                "GT_child": "0/1",
            },
            {
                "CHROM": "2",
                "POS": 3000,
                "REF": "G",
                "ALT": "A",
                "GENE": "GENE2",
                "father": "0/0",
                "mother": "0/1",
                "child": "0/1",
                "GT_father": "0/0",
                "GT_mother": "0/1",
                "GT_child": "0/1",
            },
        ]
    )

    sample_list = ["father", "mother", "child"]

    # 1. With pedigree data loaded via --ped:
    result_with_ped = analyze_inheritance(df.copy(), context.pedigree_data, sample_list)
    assert result_with_ped.iloc[0]["Inheritance_Pattern"] == "de_novo"
    assert result_with_ped.iloc[1]["Inheritance_Pattern"] == "compound_heterozygous"
    assert result_with_ped.iloc[2]["Inheritance_Pattern"] == "compound_heterozygous"

    # 2. Without pedigree (empty pedigree): falls back to unknown / singleton
    result_no_ped = analyze_inheritance(df.copy(), {}, sample_list)
    assert result_no_ped.iloc[0]["Inheritance_Pattern"] == "unknown"
    assert (
        result_no_ped.iloc[1]["Inheritance_Pattern"] == "compound_heterozygous_possible_no_pedigree"
    )
    assert (
        result_no_ped.iloc[2]["Inheritance_Pattern"] == "compound_heterozygous_possible_no_pedigree"
    )
