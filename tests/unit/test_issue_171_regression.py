"""Regression tests for Issue #171: analysis profiles and validated sample-sheet adapter."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from variantcentrifuge import cli
from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.phenotype_case_control_stages import (
    PhenotypeCaseControlAssignmentStage,
)
from variantcentrifuge.stages.sample_and_pedigree_stages import PedigreeLoadingStage
from variantcentrifuge.stages.setup_stages import PhenotypeLoadingStage


def _write_scoring_dir(base_dir: Path) -> Path:
    scoring_dir = base_dir / "scoring_model"
    scoring_dir.mkdir(parents=True, exist_ok=True)
    variable_cfg = {
        "variables": {
            "IMPACT": "impact_val|default:0",
            "dbNSFP_CADD_phred": "cadd_val|default:0",
        },
        "categories": {
            "IMPACT": {"HIGH": "3", "MODERATE": "2", "LOW": "1", "MODIFIER": "0"},
        },
    }
    formula_cfg = {
        "output_scores": ["pathogenicity_score"],
        "formulas": [
            {"pathogenicity_score": "impact_val * 10 + cadd_val"},
        ],
    }
    (scoring_dir / "variable_assignment_config.json").write_text(
        json.dumps(variable_cfg), encoding="utf-8"
    )
    (scoring_dir / "formula_config.json").write_text(json.dumps(formula_cfg), encoding="utf-8")
    return scoring_dir


@pytest.mark.unit
def test_list_profiles_outputs_builtin_analysis_profiles(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """--list-profiles should print built-in analysis profiles and exit cleanly."""
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--list-profiles"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "rare-disease-trio" in captured.out
    assert "rare-disease-singleton" in captured.out
    assert "cohort-burden" in captured.out
    assert "somatic-tumor-normal" in captured.out


@pytest.mark.unit
def test_profile_precedence_cli_over_profile_over_base_config(tmp_path: Path) -> None:
    """Explicit CLI flags must override --profile, which overrides base config.json."""
    from variantcentrifuge.analysis_profiles import apply_analysis_profile
    from variantcentrifuge.cli_sync import sync_cli_and_config
    from variantcentrifuge.config import load_config

    base_config_path = tmp_path / "base_config.json"
    base_cfg_data = load_config()
    base_cfg_data.update(
        {
            "reference": "GRCh38.p14",
            "inheritance_mode": "simple",
            "final_filter": "QUAL > 10",
            "html_report": True,
            "xlsx": False,
            "trait_type": "quantitative",
        }
    )
    base_config_path.write_text(json.dumps(base_cfg_data), encoding="utf-8")

    custom_profile_path = tmp_path / "custom_profile.json"
    custom_profile_path.write_text(
        json.dumps(
            {
                "description": "Custom test profile",
                "inheritance_mode": "columns",
                "final_filter": "QUAL > 50",
                "html_report": False,
                "xlsx": True,
                "preset": ["rare", "coding"],
            }
        ),
        encoding="utf-8",
    )

    argv = [
        "-v",
        "input.vcf.gz",
        "-c",
        str(base_config_path),
        "--profile",
        str(custom_profile_path),
        "--inheritance-mode",
        "full",
        "--preset",
        "super_rare",
    ]
    parser = cli.create_parser()
    args = parser.parse_args(argv)
    cfg = load_config(args.config)

    apply_analysis_profile(args, cfg, parser=parser, argv=argv)
    sync_cli_and_config(args, cfg)

    # From base config (not in profile or CLI)
    assert cfg["reference"] == "GRCh38.p14"
    assert cfg["trait_type"] == "quantitative"
    assert args.trait_type == "quantitative"

    # From profile (overriding base config, not in CLI)
    assert cfg["final_filter"] == "QUAL > 50"
    assert args.final_filter == "QUAL > 50"
    assert cfg["html_report"] is False
    assert args.html_report is False
    assert cfg["xlsx"] is True
    assert args.xlsx is True

    # From explicit CLI (overriding both profile and base config)
    assert cfg["inheritance_mode"] == "full"
    assert args.inheritance_mode == "full"
    assert cfg["preset"] == ["super_rare"]
    assert args.preset == ["super_rare"]


@pytest.mark.unit
def test_unknown_profile_fails_clearly(tmp_path: Path) -> None:
    """An unknown profile name or invalid profile file must fail with ValueError."""
    from variantcentrifuge.analysis_profiles import resolve_analysis_profile
    from variantcentrifuge.config import load_config

    cfg = load_config()
    with pytest.raises(ValueError, match="Analysis profile 'nonexistent-profile' not found"):
        resolve_analysis_profile("nonexistent-profile", cfg)

    bad_profile = tmp_path / "bad_profile.json"
    bad_profile.write_text(json.dumps(["not", "a", "dict"]), encoding="utf-8")
    with pytest.raises(ValueError, match="must be a JSON object"):
        resolve_analysis_profile(str(bad_profile), cfg)


@pytest.mark.unit
def test_save_profile_roundtrip_excludes_sample_specific_paths(tmp_path: Path) -> None:
    """--save-profile must exclude sample-specific/run-specific keys and round-trip cleanly."""
    saved_profile = tmp_path / "saved_profile.json"
    scoring_dir = _write_scoring_dir(tmp_path)

    rc = cli.main(
        [
            "--profile",
            "rare-disease-trio",
            "--inheritance-mode",
            "full",
            "--final-filter",
            "pathogenicity_score >= 20",
            "--scoring-config-path",
            str(scoring_dir),
            "--case-phenotypes",
            "HP:0001250",
            "--output-dir",
            str(tmp_path / "run_out"),
            "--ped",
            str(tmp_path / "family.ped"),
            "--phenotype-file",
            str(tmp_path / "pheno.tsv"),
            "--bam-mapping-file",
            str(tmp_path / "bams.tsv"),
            "--case-samples",
            "PROBAND",
            "--control-samples",
            "FATHER,MOTHER",
            "--save-profile",
            str(saved_profile),
        ]
    )
    assert rc == 0
    assert saved_profile.is_file()

    saved_data = json.loads(saved_profile.read_text(encoding="utf-8"))
    for forbidden_key in (
        "vcf_file",
        "from_tsv",
        "output_dir",
        "output_file",
        "ped",
        "ped_file",
        "phenotype_file",
        "bam_mapping_file",
        "sample_sheet",
        "case_samples",
        "control_samples",
        "case_samples_file",
        "control_samples_file",
        "profiles",
        "field_profiles",
        "presets",
    ):
        assert forbidden_key not in saved_data

    # Reusable settings are preserved
    assert saved_data["inheritance_mode"] == "full"
    assert saved_data["final_filter"] == "pathogenicity_score >= 20"
    assert saved_data["scoring_config_path"] == str(scoring_dir)
    assert saved_data["preset"] == ["rare", "coding"]
    assert saved_data["html_report"] is True
    assert saved_data["xlsx"] is True
    assert saved_data["case_phenotypes"] == "HP:0001250"

    # Validate the saved profile via validate-config subcommand
    rc_val = cli.main(["validate-config", "--profile", str(saved_profile)])
    assert rc_val == 0


@pytest.mark.unit
def test_sample_sheet_validation_errors(tmp_path: Path) -> None:
    """--sample-sheet rejects duplicate sample_id, bad sex/affected_status, bad HPO, conflicts."""
    from variantcentrifuge.sample_sheet import load_sample_sheet

    # 1. Duplicate sample_id
    dup_sheet = tmp_path / "dup.tsv"
    dup_sheet.write_text(
        "family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\taffected_status\n"
        "FAM1\tS1\t0\t0\t1\t2\n"
        "FAM1\tS1\t0\t0\t2\t1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Duplicate sample_id 'S1'"):
        load_sample_sheet(str(dup_sheet))

    # 2. Invalid sex
    bad_sex = tmp_path / "bad_sex.tsv"
    bad_sex.write_text(
        "sample_id\tsex\taffected_status\nS1\t3\t2\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Invalid sex"):
        load_sample_sheet(str(bad_sex))

    # 3. Invalid affected_status
    bad_status = tmp_path / "bad_status.tsv"
    bad_status.write_text(
        "sample_id\tsex\taffected_status\nS1\t1\t9\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Invalid affected_status"):
        load_sample_sheet(str(bad_status))

    # 4. Malformed HPO term
    bad_hpo = tmp_path / "bad_hpo.tsv"
    bad_hpo.write_text(
        "sample_id\thpo_terms\nS1\tHP:0001250;HP:abc\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Malformed HPO term 'HP:abc'"):
        load_sample_sheet(str(bad_hpo))

    # 5. Conflict with --ped, --phenotype-file, or --bam-mapping-file
    valid_sheet = tmp_path / "valid.tsv"
    valid_sheet.write_text(
        "family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\taffected_status\thpo_terms\tbam_path\n"
        "FAM1\tPROBAND\tFATHER\tMOTHER\t1\t2\tHP:0001250\t/path/proband.bam\n",
        encoding="utf-8",
    )
    input_tsv = tmp_path / "in.tsv"
    input_tsv.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\tGT\nchr1\t100\tA\tG\tBRCA1\tPROBAND(0/1)\n",
        encoding="utf-8",
    )
    for conflict_flag, conflict_val in (
        ("--ped", str(tmp_path / "f.ped")),
        ("--phenotype-file", str(tmp_path / "p.tsv")),
        ("--bam-mapping-file", str(tmp_path / "b.tsv")),
    ):
        rc = cli.main(
            [
                "--from-tsv",
                str(input_tsv),
                "--sample-sheet",
                str(valid_sheet),
                conflict_flag,
                conflict_val,
                "-o",
                str(tmp_path / "out.tsv"),
            ]
        )
        assert rc == 1


@pytest.mark.unit
def test_sample_sheet_equivalence_and_checkpoint_restore(tmp_path: Path) -> None:
    """--sample-sheet produces identical outputs and checkpoint restoration to legacy flags."""
    # Create legacy files
    ped_path = tmp_path / "family.ped"
    ped_path.write_text(
        "FAM1\tPROBAND\tFATHER\tMOTHER\t1\t2\nFAM1\tFATHER\t0\t0\t1\t1\nFAM1\tMOTHER\t0\t0\t2\t1\n",
        encoding="utf-8",
    )

    pheno_path = tmp_path / "phenotypes.tsv"
    pheno_path.write_text(
        "sample_id\thpo_terms\nPROBAND\tHP:0001250\nPROBAND\tHP:0001263\nFATHER\t\nMOTHER\t\n",
        encoding="utf-8",
    )

    bam_path = tmp_path / "bams.tsv"
    bam_path.write_text(
        "PROBAND\t/data/proband.bam\nFATHER\t/data/father.bam\nMOTHER\t/data/mother.bam\n",
        encoding="utf-8",
    )

    # Create equivalent sample sheet
    sheet_path = tmp_path / "sample_sheet.tsv"
    sheet_path.write_text(
        "family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\taffected_status\thpo_terms\tbam_path\n"
        "FAM1\tPROBAND\tFATHER\tMOTHER\t1\t2\tHP:0001250;HP:0001263\t/data/proband.bam\n"
        "FAM1\tFATHER\t0\t0\t1\t1\t\t/data/father.bam\n"
        "FAM1\tMOTHER\t0\t0\t2\t1\t\t/data/mother.bam\n",
        encoding="utf-8",
    )

    # Verify checkpoint skip/restore on Pedigree, Phenotype, and Case/Control stages
    ctx = PipelineContext(
        args=argparse.Namespace(sample_sheet=str(sheet_path)),
        config={
            "sample_sheet": str(sheet_path),
            "case_phenotypes": ["HP:0001250"],
            "control_phenotypes": [],
        },
        workspace=Workspace(tmp_path / "ckpt_ws", "test"),
    )
    ctx.vcf_samples = ["PROBAND", "FATHER", "MOTHER"]

    PedigreeLoadingStage()._handle_checkpoint_skip(ctx)
    PhenotypeLoadingStage()._handle_checkpoint_skip(ctx)
    PhenotypeCaseControlAssignmentStage()._handle_checkpoint_skip(ctx)

    assert ctx.pedigree_data is not None
    assert ctx.pedigree_data["PROBAND"]["father_id"] == "FATHER"
    assert ctx.pedigree_data["PROBAND"]["mother_id"] == "MOTHER"
    assert ctx.pedigree_data["PROBAND"]["affected_status"] == "2"
    assert ctx.phenotype_data == {
        "PROBAND": {"HP:0001250", "HP:0001263"},
        "FATHER": set(),
        "MOTHER": set(),
    }
    assert ctx.config["case_samples"] == ["PROBAND"]
    assert set(ctx.config["control_samples"]) == {"FATHER", "MOTHER"}

    # End-to-end pipeline comparison: legacy flags vs --sample-sheet
    input_tsv = tmp_path / "variants.tsv"
    input_tsv.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\tIMPACT\tdbNSFP_CADD_phred\tGT\n"
        "chr1\t1000\tA\tG\tBRCA1\tHIGH\t28.0\tPROBAND(0/1);FATHER(0/0);MOTHER(0/0)\n"
        "chr2\t2000\tC\tT\tTP53\tMODERATE\t22.0\tPROBAND(0/1);FATHER(0/1);MOTHER(0/0)\n",
        encoding="utf-8",
    )
    scoring_dir = _write_scoring_dir(tmp_path)

    legacy_dir = tmp_path / "legacy_out"
    legacy_tsv = legacy_dir / "results.tsv"
    rc_legacy = cli.main(
        [
            "--from-tsv",
            str(input_tsv),
            "--ped",
            str(ped_path),
            "--phenotype-file",
            str(pheno_path),
            "--phenotype-sample-column",
            "sample_id",
            "--phenotype-value-column",
            "hpo_terms",
            "--bam-mapping-file",
            str(bam_path),
            "--case-phenotypes",
            "HP:0001250",
            "--inheritance-mode",
            "columns",
            "--scoring-config-path",
            str(scoring_dir),
            "--output-dir",
            str(legacy_dir),
            "-o",
            str(legacy_tsv),
            "--no-stats",
        ]
    )
    assert rc_legacy == 0

    sheet_dir = tmp_path / "sheet_out"
    sheet_tsv = sheet_dir / "results.tsv"
    rc_sheet = cli.main(
        [
            "--from-tsv",
            str(input_tsv),
            "--sample-sheet",
            str(sheet_path),
            "--case-phenotypes",
            "HP:0001250",
            "--inheritance-mode",
            "columns",
            "--scoring-config-path",
            str(scoring_dir),
            "--output-dir",
            str(sheet_dir),
            "-o",
            str(sheet_tsv),
            "--no-stats",
        ]
    )
    assert rc_sheet == 0

    df_legacy = pd.read_csv(legacy_tsv, sep="\t")
    df_sheet = pd.read_csv(sheet_tsv, sep="\t")

    pd.testing.assert_frame_equal(df_legacy, df_sheet)


@pytest.mark.unit
def test_profile_and_sample_sheet_edge_cases(tmp_path: Path) -> None:
    """Cover edge cases in analysis_profiles.py and sample_sheet.py."""
    from variantcentrifuge.analysis_profiles import (
        apply_analysis_profile,
        get_explicit_cli_dests,
        list_analysis_profiles,
        resolve_analysis_profile,
    )
    from variantcentrifuge.sample_sheet import (
        SampleSheetData,
        classify_samples_from_phenotype_data,
        load_sample_sheet,
        materialize_sample_sheet_bam_mapping,
    )

    assert get_explicit_cli_dests(None, None) == set()
    profs = list_analysis_profiles({})
    assert any(p["name"] == "rare-disease-trio" for p in profs)

    with pytest.raises(ValueError, match="must not be empty"):
        resolve_analysis_profile("   ", {})

    with pytest.raises(ValueError, match="not a file"):
        resolve_analysis_profile(str(tmp_path), {})

    bad_json = tmp_path / "syntax_err.json"
    bad_json.write_text("{not valid json", encoding="utf-8")
    with pytest.raises(ValueError, match="Failed to parse analysis profile JSON"):
        resolve_analysis_profile(str(bad_json), {})

    wrapped_json = tmp_path / "wrapped.json"
    wrapped_json.write_text(
        json.dumps(
            {
                "description": "Wrapped profile",
                "profile": {"inheritance_mode": "columns", "ped_file": "fam.ped"},
            }
        ),
        encoding="utf-8",
    )
    resolved = resolve_analysis_profile(str(wrapped_json), {})
    assert resolved == {"inheritance_mode": "columns", "ped_file": "fam.ped"}

    args = argparse.Namespace(profile=str(wrapped_json), ped=None, ped_file=None)
    cfg: dict[str, object] = {}
    apply_analysis_profile(args, cfg)
    assert args.ped == "fam.ped"

    with pytest.raises(FileNotFoundError, match="Sample sheet file not found"):
        load_sample_sheet(tmp_path / "missing_sheet.tsv")

    empty_sheet = tmp_path / "empty.tsv"
    empty_sheet.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="is empty"):
        load_sample_sheet(empty_sheet)

    comments_only = tmp_path / "comments.tsv"
    comments_only.write_text("# comment only\n", encoding="utf-8")
    with pytest.raises(ValueError, match="contains no header or data"):
        load_sample_sheet(comments_only)

    no_sample_col = tmp_path / "no_col.tsv"
    no_sample_col.write_text("family_id\tsex\nFAM1\t1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="must contain a 'sample_id' column"):
        load_sample_sheet(no_sample_col)

    empty_sid = tmp_path / "empty_sid.tsv"
    empty_sid.write_text("sample_id\tsex\n\t1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Empty sample_id"):
        load_sample_sheet(empty_sid)

    header_only = tmp_path / "header_only.tsv"
    header_only.write_text("sample_id\tsex\n", encoding="utf-8")
    with pytest.raises(ValueError, match="contains header only"):
        load_sample_sheet(header_only)

    assert materialize_sample_sheet_bam_mapping(SampleSheetData(), tmp_path) is None

    pheno_map = {
        "S1_sub": {"HP:0001250"},
        "S2_sub": {"HP:0002000"},
        "S3_sub": set(),
    }
    cases, ctrls = classify_samples_from_phenotype_data(
        ["S1_sub", "S2_sub", "S3_sub"],
        pheno_map,
        case_phenotypes=[],
        control_phenotypes=["HP:0002000"],
        remove_substring="_sub",
    )
    assert ctrls == ["S2"]
    assert set(cases) == {"S1", "S3"}

    cases_both, ctrls_both = classify_samples_from_phenotype_data(
        ["S1_sub", "S2_sub", "S3_sub"],
        pheno_map,
        case_phenotypes=["HP:0001250"],
        control_phenotypes=["HP:0002000"],
        remove_substring="_sub",
    )
    assert cases_both == ["S1"]
    assert ctrls_both == ["S2"]

    # Attached short options (-oattached.tsv, -rGRCh38) and abbreviated long options
    parser = cli.create_parser()
    abbr_argv = [
        "-v",
        "in.vcf.gz",
        "--profile",
        "rare-disease-trio",
        "-oattached.tsv",
        "-rGRCh38",
        "--inheritance-mo",
        "full",
    ]
    abbr_args = parser.parse_args(abbr_argv)
    abbr_cfg: dict[str, object] = {}
    apply_analysis_profile(abbr_args, abbr_cfg, parser=parser, argv=abbr_argv)
    assert abbr_cfg["output_file"] == "attached.tsv"
    assert abbr_cfg["reference"] == "GRCh38"
    assert abbr_cfg["inheritance_mode"] == "full"

    # Chained boolean short flags in get_explicit_cli_dests
    short_parser = argparse.ArgumentParser()
    short_parser.add_argument("-a", action="store_true", dest="flag_a")
    short_parser.add_argument("-b", action="store_true", dest="flag_b")
    assert get_explicit_cli_dests(short_parser, ["-", "-ab"]) == {"flag_a", "flag_b"}

    # Non-dict profile entry in config raises ValueError
    with pytest.raises(ValueError, match="must be a JSON object"):
        resolve_analysis_profile("bad_entry", {"profiles": {"bad_entry": "not_a_dict"}})

    # Relative bam_path resolved against sample_sheet directory
    rel_dir = tmp_path / "subdir"
    rel_dir.mkdir(parents=True, exist_ok=True)
    rel_sheet = rel_dir / "samples.tsv"
    rel_sheet.write_text(
        "sample_id\tbam_path\nS1\tbams/s1.bam\nS2\t/abs/path/s2.bam\n",
        encoding="utf-8",
    )
    rel_loaded = load_sample_sheet(rel_sheet)
    assert rel_loaded.bam_mapping["S1"] == str((rel_dir / "bams/s1.bam").resolve())
    assert rel_loaded.bam_mapping["S2"] == "/abs/path/s2.bam"


@pytest.mark.unit
def test_cli_config_and_sync_coverage(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Cover remaining branches in cli.py, cli_config.py, cli_sync.py, and cli_subcommands.py."""
    from variantcentrifuge import cli_config
    from variantcentrifuge.cli_config import (
        populate_cli_config,
        resolve_presets_and_tumor_normal,
    )
    from variantcentrifuge.cli_sync import sync_cli_and_config
    from variantcentrifuge.config import load_config

    parser = cli.create_parser()
    monkeypatch.setattr(
        cli_config,
        "resolve_tumor_normal_indices",
        lambda *a, **kw: (
            str(kw.get("cli_tumor_idx") or "1"),
            str(kw.get("cli_normal_idx") or "0"),
        ),
    )

    # 1. resolve_presets_and_tumor_normal with empty presets dict and tumor/normal vars
    args = parser.parse_args(
        [
            "-v",
            "in.vcf.gz",
            "--preset",
            "rare",
            "--filters",
            "GEN[{tumor_idx}].DP >= {tumor_dp_min} & GEN[{normal_idx}].DP >= {normal_dp_min}",
            "--final-filter",
            "INFO_AF < {tumor_af_min}",
            "--tumor-sample-index",
            "1",
            "--normal-sample-index",
            "0",
        ]
    )
    cfg_empty_presets: dict[str, Any] = {"presets": {}}
    filt, tn_vars = resolve_presets_and_tumor_normal(args, cfg_empty_presets, args.filters)
    assert filt is not None and "GEN[1].DP >= 20" in filt
    assert tn_vars["tumor_idx"] == "1"
    assert tn_vars["normal_idx"] == "0"

    # Unknown preset raises ValueError
    args_bad_preset = parser.parse_args(["-v", "in.vcf.gz", "--preset", "no_such_preset"])
    with pytest.raises(ValueError, match="not found in config"):
        resolve_presets_and_tumor_normal(args_bad_preset, load_config(), None)

    # 2. populate_cli_config branches: association_tests, coast, igv, exceptions, annotations
    exc_file = tmp_path / "exceptions.tsv"
    exc_file.write_text("chr1\t100\tA\tG\n", encoding="utf-8")
    weights_file = tmp_path / "weights.tsv"
    weights_file.write_text("gene\tweight\nBRCA1\t1.0\n", encoding="utf-8")
    regions_file = tmp_path / "regions.bed"
    regions_file.write_text("chr1\t0\t1000\n", encoding="utf-8")

    args_pop = parser.parse_args(
        [
            "-v",
            "in.vcf.gz",
            "--perform-association",
            "--association-tests",
            "fisher,coast",
            "--coast-weights",
            "1.0,2.0,3.0",
            "--variant-weight-params",
            '{"a": 1}',
            "--gene-prior-weights",
            str(weights_file),
            "--regions-bed",
            str(regions_file),
            "--add-chr",
            "--igv",
            "--bam-mapping-file",
            "bams.tsv",
            "--igv-fasta",
            "ref.fa",
            "--igv-reference",
            "hg38",
            "--filter-exception-file",
            str(exc_file),
            "--append-extra-sample-fields",
            "DP",
            "AD",
            "--extra-sample-field-delimiter",
            "|",
            "--annotate-bed",
            "b.bed",
            "--annotate-gene-list",
            "g.txt",
            "--annotate-json-genes",
            "j.json",
            "--json-gene-mapping",
            "identifier=symbol",
            "--sort-memory-limit",
            "2G",
        ]
    )
    cfg_pop = load_config()
    populate_cli_config(
        args_pop,
        cfg_pop,
        None,
        "GRCh38",
        "PASS",
        "CHROM,POS,REF,ALT",
        tn_vars,
    )
    assert cfg_pop["coast_weights"] == [1.0, 2.0, 3.0]
    assert cfg_pop["variant_weight_params"] == {"a": 1}
    assert cfg_pop["add_chr"] is True
    assert "RAW_FILTER" in str(cfg_pop["fields_to_extract"])
    assert cfg_pop["append_extra_sample_fields"] is True
    assert cfg_pop["extra_sample_field_delimiter"] == "|"
    assert cfg_pop["coast_classification"] is not None

    # COAST with empty fields string
    cfg_coast_empty = load_config()
    populate_cli_config(args_pop, cfg_coast_empty, None, "GRCh38", "PASS", "", tn_vars)
    assert cfg_coast_empty["fields_to_extract"]

    # populate_cli_config validation error branches (with parser=None -> ValueError)
    def _make_base_ns(**overrides: Any) -> argparse.Namespace:
        ns = parser.parse_args(["-v", "in.vcf.gz"])
        for k, v in overrides.items():
            setattr(ns, k, v)
        return ns

    with pytest.raises(ValueError, match="--variant-weight-params"):
        populate_cli_config(
            _make_base_ns(variant_weight_params="[1,2]"),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="--coast-weights must have exactly 3"):
        populate_cli_config(
            _make_base_ns(coast_weights="1.0,2.0"),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="--coast-weights values must be numeric"):
        populate_cli_config(
            _make_base_ns(coast_weights="a,b,c"),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="--gene-prior-weights file not found"):
        populate_cli_config(
            _make_base_ns(gene_prior_weights=str(tmp_path / "missing_w.tsv")),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="--bam-mapping-file must be provided"):
        populate_cli_config(
            _make_base_ns(igv=True, bam_mapping_file=None),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="either --igv-reference or --igv-fasta"):
        populate_cli_config(
            _make_base_ns(igv=True, bam_mapping_file="b.tsv", igv_fasta=None, igv_reference=None),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(FileNotFoundError, match="COAST classification model not found"):
        populate_cli_config(
            _make_base_ns(
                perform_association=True,
                association_tests=["coast"],
                coast_classification="nonexistent_model",
            ),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="--annotate-json-genes requires --json-gene-mapping"):
        populate_cli_config(
            _make_base_ns(annotate_json_genes=["a.json"], json_gene_mapping=None),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="--json-gene-mapping requires --annotate-json-genes"):
        populate_cli_config(
            _make_base_ns(annotate_json_genes=None, json_gene_mapping="id=sym"),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="--regions-bed file not found"):
        populate_cli_config(
            _make_base_ns(regions_bed=str(tmp_path / "no.bed")),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="--association-tests requires --perform-association"):
        populate_cli_config(
            _make_base_ns(perform_association=False, association_tests="fisher"),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="--variant-weight-column requires --perform-association"):
        populate_cli_config(
            _make_base_ns(perform_association=False, variant_weight_column="W"),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="when --variant-weights score_column"):
        populate_cli_config(
            _make_base_ns(
                perform_association=True,
                variant_weights="score_column",
                variant_weight_column=None,
            ),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    with pytest.raises(ValueError, match="requires a non-empty column name"):
        populate_cli_config(
            _make_base_ns(perform_association=True, variant_weights="column:   "),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    for flag_name, flag_val in (
        ("covariate_file", "cov.tsv"),
        ("diagnostics_output", "diag_dir"),
        ("pca", "pca.tsv"),
    ):
        with pytest.raises(ValueError, match="requires --perform-association"):
            populate_cli_config(
                _make_base_ns(perform_association=False, **{flag_name: flag_val}),
                load_config(),
                None,
                "GRCh38",
                "PASS",
                "CHROM",
                tn_vars,
            )

    with pytest.raises(ValueError, match="incompatible with logistic_burden"):
        populate_cli_config(
            _make_base_ns(
                perform_association=True,
                trait_type="quantitative",
                association_tests="logistic_burden",
            ),
            load_config(),
            None,
            "GRCh38",
            "PASS",
            "CHROM",
            tn_vars,
        )

    # 3. sync_cli_and_config remaining branches
    ns_sync = argparse.Namespace(
        vcf_file=None,
        from_tsv=None,
        gene_name=None,
        gene_file=None,
        pseudonymize_table="pseudo.tsv",
        gzip_intermediates=None,
        max_read_depth=None,
        transcript_file="tx.txt",
        annotate_gene_list="genes.txt",
        annotate_json_genes="genes.json",
        json_gene_mapping="id=sym",
        remove_snpeff_sep=True,
        do_not_remove_additional_columns=None,
        extra_columns_to_remove=None,
    )
    cfg_sync: dict[str, Any] = {
        "from_tsv": "input.tsv",
        "max_Read_depth": 50,
        "do_not_remove_additional_columns": True,
    }
    sync_cli_and_config(ns_sync, cfg_sync)
    assert ns_sync.from_tsv == "input.tsv"
    assert ns_sync.pseudonymize is True
    assert ns_sync.max_read_depth == 50
    assert cfg_sync["transcript_file"] == "tx.txt"
    assert cfg_sync["annotate_gene_lists"] == "genes.txt"
    assert cfg_sync["annotate_json_genes"] == "genes.json"
    assert cfg_sync["json_gene_mapping"] == "id=sym"
    assert cfg_sync["remove_snpeff_sep"] is True
    assert ns_sync.do_not_remove_additional_columns is True

    # 4. cli.main debug logging, log_file, and invalid profile return code
    log_file = tmp_path / "logs" / "run.log"
    rc_bad_prof = cli.main(
        [
            "-v",
            "in.vcf.gz",
            "--log-level",
            "DEBUG",
            "--log-file",
            str(log_file),
            "--profile",
            "nonexistent_profile_xyz",
        ]
    )
    assert rc_bad_prof == 1
    assert log_file.is_file()
