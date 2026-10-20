"""Subcommand dispatch and early-exit helpers for VariantCentrifuge CLI."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

logger = logging.getLogger("variantcentrifuge")

RECOGNIZED_SUBCOMMANDS: frozenset[str] = frozenset(
    {
        "run",
        "report",
        "validate-config",
        "build-pm5",
        "download-references",
    }
)


def extract_leading_subcommand(argv: list[str] | None = None) -> tuple[str | None, list[str]]:
    """Extract a recognized leading subcommand if present as the first token.

    Option flags (starting with '-') and option values following flags are never
    treated as subcommands.
    """
    effective = list(sys.argv[1:]) if argv is None else list(argv)
    if effective and not effective[0].startswith("-") and effective[0] in RECOGNIZED_SUBCOMMANDS:
        return effective[0], effective[1:]
    return None, effective


def create_report_subcommand_parser() -> argparse.ArgumentParser:
    """Create argument parser for 'variantcentrifuge report'."""
    parser = argparse.ArgumentParser(
        prog="variantcentrifuge report",
        description="Generate HTML and/or Excel reports directly from an existing TSV file.",
    )
    parser.add_argument(
        "-i",
        "--input-tsv",
        "--tsv",
        dest="input_tsv",
        required=True,
        help="Path to existing variants TSV (or .tsv.gz) file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        default="output",
        help="Output directory for generated reports (default: output).",
    )
    parser.add_argument(
        "--output-file",
        default=None,
        help="Optional explicit path for the Excel workbook output.",
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        help="Optional path to configuration JSON file.",
    )
    parser.add_argument(
        "--html",
        "--html-report",
        dest="html_report",
        action="store_true",
        help="Generate interactive HTML report.",
    )
    parser.add_argument(
        "--xlsx",
        "--excel",
        dest="xlsx",
        action="store_true",
        help="Generate Excel workbook report.",
    )
    return parser


def run_report_subcommand(argv: list[str]) -> int:
    """Execute the 'report' subcommand without starting the VCF extraction pipeline."""
    from .config import load_config
    from .converter import convert_to_excel, produce_report_json
    from .field_profile import resolve_profile
    from .generate_html_report import generate_html_report
    from .links import resolve_links

    parser = create_report_subcommand_parser()
    parsed = parser.parse_args(argv)

    tsv_path = Path(parsed.input_tsv)
    if not tsv_path.is_file() or tsv_path.stat().st_size == 0:
        logger.error("Input TSV file not found or empty: %s", tsv_path)
        return 1

    try:
        cfg = load_config(parsed.config)
        resolve_profile(cfg)
        resolve_links(cfg)
        out_dir = Path(parsed.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        do_html = parsed.html_report or (not parsed.html_report and not parsed.xlsx)
        do_xlsx = parsed.xlsx

        if do_html:
            report_dir = out_dir / "report"
            report_dir.mkdir(parents=True, exist_ok=True)
            produce_report_json(str(tsv_path), str(out_dir))
            variants_json = report_dir / "variants.json"
            summary_json = report_dir / "summary.json"
            generate_html_report(
                variants_json=str(variants_json),
                summary_json=str(summary_json),
                output_dir=str(report_dir),
                cfg=cfg,
            )

        if do_xlsx:
            if parsed.output_file:
                xlsx_target = Path(parsed.output_file)
            else:
                stem = tsv_path.name
                for suffix in (".tsv.gz", ".tsv", ".gz"):
                    if stem.endswith(suffix):
                        stem = stem[: -len(suffix)]
                        break
                xlsx_target = out_dir / f"{stem}.xlsx"
            xlsx_target.parent.mkdir(parents=True, exist_ok=True)
            convert_to_excel(str(tsv_path), cfg, output_file=str(xlsx_target))

        return 0
    except Exception as exc:
        logger.error("Report generation failed: %s", exc)
        return 1


def create_validate_config_subcommand_parser() -> argparse.ArgumentParser:
    """Create argument parser for 'variantcentrifuge validate-config'."""
    parser = argparse.ArgumentParser(
        prog="variantcentrifuge validate-config",
        description="Validate VariantCentrifuge configuration, profiles, scoring, and metadata.",
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        help="Path to configuration JSON file (defaults to package config.json).",
    )
    parser.add_argument(
        "--field-profile",
        default=None,
        help="Field profile name to validate against configuration.",
    )
    parser.add_argument(
        "--profile",
        default=None,
        help="Analysis profile name or file path to validate.",
    )
    parser.add_argument(
        "--scoring-config-path",
        default=None,
        help="Path to scoring model directory to validate.",
    )
    parser.add_argument(
        "-p",
        "--ped",
        "--ped-file",
        dest="ped",
        default=None,
        help="Path to pedigree (.ped) file to validate.",
    )
    parser.add_argument(
        "--phenotype-file",
        default=None,
        help="Path to phenotype file to validate.",
    )
    parser.add_argument(
        "--phenotype-sample-column",
        default="sample_id",
        help="Sample column name in phenotype file (default: sample_id).",
    )
    parser.add_argument(
        "--phenotype-value-column",
        default="phenotype",
        help="Phenotype value column name in phenotype file (default: phenotype).",
    )
    parser.add_argument(
        "--sample-sheet",
        default=None,
        help="Path to sample sheet manifest to validate.",
    )
    parser.add_argument(
        "--preset",
        action="append",
        default=None,
        help="Optional preset name(s) to validate against configuration.",
    )
    return parser


def run_validate_config_subcommand(argv: list[str]) -> int:
    """Validate configuration and companion metadata files without running the pipeline."""
    from .config import load_config
    from .field_profile import resolve_profile
    from .ped_reader import read_pedigree
    from .scoring import read_scoring_config
    from .validators import validate_phenotype_file

    parser = create_validate_config_subcommand_parser()
    parsed = parser.parse_args(argv)

    try:
        cfg = load_config(parsed.config)
        if not isinstance(cfg.get("presets", {}), dict):
            raise ValueError("Configuration 'presets' must be a JSON object.")
        if not isinstance(cfg.get("field_profiles", {}), dict):
            raise ValueError("Configuration 'field_profiles' must be a JSON object.")

        if parsed.profile:
            from .analysis_profiles import resolve_analysis_profile

            prof_data = resolve_analysis_profile(str(parsed.profile), cfg)
            cfg.update(prof_data)

        if parsed.field_profile:
            cfg["field_profile"] = parsed.field_profile
        resolve_profile(cfg)

        presets_dict = cfg.get("presets", {})
        raw_presets = parsed.preset or cfg.get("preset") or cfg.get("presets_to_apply")
        if raw_presets:
            preset_items = [raw_presets] if isinstance(raw_presets, str) else list(raw_presets)
            for item in preset_items:
                for p_name in str(item).split(","):
                    p_clean = p_name.strip()
                    if p_clean and p_clean not in presets_dict:
                        raise ValueError(f"Preset '{p_clean}' not found in configuration.")

        scoring_path = parsed.scoring_config_path or cfg.get("scoring_config_path")
        if scoring_path:
            scoring_cfg = read_scoring_config(str(scoring_path))
            if not isinstance(scoring_cfg.get("variables"), dict) or not scoring_cfg["variables"]:
                raise ValueError(
                    "Scoring configuration 'variables' must be a non-empty JSON object."
                )
            if not isinstance(scoring_cfg.get("formulas"), list) or not scoring_cfg["formulas"]:
                raise ValueError("Scoring configuration 'formulas' must be a non-empty JSON array.")

        sample_sheet_path = parsed.sample_sheet or cfg.get("sample_sheet")
        if sample_sheet_path:
            from .sample_sheet import validate_sample_sheet_conflicts

            validate_sample_sheet_conflicts(parsed, cfg)

        ped_path = parsed.ped or cfg.get("ped_file")
        if ped_path:
            if not Path(ped_path).is_file():
                raise FileNotFoundError(f"Pedigree file not found: {ped_path}")
            read_pedigree(str(ped_path))

        pheno_path = parsed.phenotype_file or cfg.get("phenotype_file")
        if pheno_path:
            s_col = parsed.phenotype_sample_column or cfg.get(
                "phenotype_sample_column", "sample_id"
            )
            v_col = parsed.phenotype_value_column or cfg.get("phenotype_value_column", "phenotype")
            validate_phenotype_file(str(pheno_path), str(s_col), str(v_col), logger)

        print("Configuration valid.")
        return 0
    except Exception as exc:
        logger.error("Configuration validation failed: %s", exc)
        return 1


def dispatch_subcommand(subcommand: str, sub_argv: list[str]) -> int:
    """Dispatch non-'run' subcommands to their dedicated handlers."""
    if subcommand == "report":
        return run_report_subcommand(sub_argv)
    if subcommand == "validate-config":
        return run_validate_config_subcommand(sub_argv)
    if subcommand in ("build-pm5", "download-references"):
        from . import build_pm5_lookup

        return build_pm5_lookup.main(sub_argv)
    raise ValueError(f"Unsupported subcommand: {subcommand}")


def handle_early_utility_options(
    effective_argv: list[str],
    load_config_fn: Any = None,
) -> None:
    """Handle early-exit utility flags for field profiles, VCF annotations, and checkpoints."""
    if load_config_fn is None:
        from .config import load_config as load_config_fn

    if "--list-profiles" in effective_argv:
        from .analysis_profiles import list_analysis_profiles

        ap_parser = argparse.ArgumentParser(description="List analysis profiles")
        ap_parser.add_argument("--list-profiles", action="store_true")
        ap_parser.add_argument("-c", "--config", help="Path to configuration file", default=None)
        ap_args = ap_parser.parse_args(effective_argv)

        ap_cfg = load_config_fn(ap_args.config)
        analysis_profiles = list_analysis_profiles(ap_cfg)
        if analysis_profiles:
            print("Available analysis profiles:")
            for prof in analysis_profiles:
                print(f"  {prof['name']:24s}  {prof['description']}")
        else:
            print("No analysis profiles defined in configuration.")
        sys.exit(0)

    if "--list-field-profiles" in effective_argv:
        from .field_profile import list_profiles

        profile_parser = argparse.ArgumentParser(description="List field profiles")
        profile_parser.add_argument("--list-field-profiles", action="store_true")
        profile_parser.add_argument(
            "-c", "--config", help="Path to configuration file", default=None
        )
        profile_args = profile_parser.parse_args(effective_argv)

        profile_cfg = load_config_fn(profile_args.config)
        profiles = list_profiles(profile_cfg)
        if profiles:
            print("Available field profiles:")
            for prof in profiles:
                print(f"  {prof['name']:12s}  {prof['description']}")
        else:
            print("No field profiles defined in configuration.")
        sys.exit(0)

    if "--show-vcf-annotations" in effective_argv:
        from .vcf_header_parser import (
            format_fields_json,
            format_fields_table,
            parse_vcf_header,
        )

        ann_parser = argparse.ArgumentParser(description="Show VCF annotations")
        ann_parser.add_argument("--show-vcf-annotations", action="store_true")
        ann_parser.add_argument("-v", "--vcf-file", required=True)
        ann_parser.add_argument("--annotation-filter", type=str, default=None)
        ann_parser.add_argument("--annotation-format", choices=["table", "json"], default="table")
        ann_parser.add_argument("--info-only", action="store_true", default=False)
        ann_parser.add_argument("--format-only", action="store_true", default=False)
        ann_args, _ = ann_parser.parse_known_args(effective_argv)

        ann_fields = parse_vcf_header(ann_args.vcf_file)
        if ann_args.annotation_format == "json":
            output = format_fields_json(
                ann_fields,
                pattern=ann_args.annotation_filter,
                info_only=ann_args.info_only,
                format_only=ann_args.format_only,
            )
        else:
            output = format_fields_table(
                ann_fields,
                pattern=ann_args.annotation_filter,
                info_only=ann_args.info_only,
                format_only=ann_args.format_only,
            )
        print(output)
        sys.exit(0)

    if "--show-checkpoint-status" in effective_argv:
        status_parser = argparse.ArgumentParser(description="Show checkpoint status")
        status_parser.add_argument("--show-checkpoint-status", action="store_true")
        status_parser.add_argument("--output-dir", default="output", help="Output directory")
        status_parser.add_argument(
            "--log-level", choices=["DEBUG", "INFO", "WARN", "WARNING", "ERROR"], default="INFO"
        )
        status_parser.add_argument(
            "-c", "--config", help="Path to configuration file", default=None
        )
        status_args = status_parser.parse_args(effective_argv)

        log_level_map = {
            "DEBUG": logging.DEBUG,
            "INFO": logging.INFO,
            "WARN": logging.WARNING,
            "WARNING": logging.WARNING,
            "ERROR": logging.ERROR,
        }
        logging.getLogger("variantcentrifuge").setLevel(log_level_map[status_args.log_level])

        from .checkpoint import PipelineState

        print("Checking checkpoint status for stage-based pipeline...")

        pipeline_state = PipelineState(status_args.output_dir)
        if pipeline_state.load():
            print(pipeline_state.get_summary())
        else:
            print(f"No checkpoint state found in {status_args.output_dir}")
        sys.exit(0)


def handle_listing_and_resume_options(args: argparse.Namespace, cfg: dict[str, Any]) -> None:
    """Handle early-exit listing options and interactive resume."""
    if getattr(args, "list_stages", False):
        from .display_utils import display_available_stages
        from .pipeline import create_stages_from_config
        from .stages.stage_registry import initialize_registry

        initialize_registry()
        try:
            stages = create_stages_from_config(cfg)
            display_available_stages(stages, cfg)
        except Exception as e:
            logger.error(f"Failed to create stages: {e}")
            sys.exit(1)
        sys.exit(0)

    if getattr(args, "list_checkpoints", False):
        from .checkpoint import PipelineState
        from .display_utils import display_enhanced_status

        pipeline_state = PipelineState(args.output_dir)
        if pipeline_state.load():
            display_enhanced_status(pipeline_state, [])
        else:
            print(f"❌ No checkpoint file found in {args.output_dir}")
            print("   Start a new pipeline run with --enable-checkpoint to create checkpoints.")
        sys.exit(0)

    if getattr(args, "interactive_resume", False):
        from .checkpoint import PipelineState
        from .interactive_resume import handle_interactive_resume
        from .pipeline import create_stages_from_config
        from .stages.stage_registry import initialize_registry

        initialize_registry()
        pipeline_state = PipelineState(args.output_dir)
        if not pipeline_state.load():
            print(f"❌ No checkpoint file found in {args.output_dir}")
            print("   Start a new pipeline run with --enable-checkpoint to create checkpoints.")
            sys.exit(1)

        try:
            stages = create_stages_from_config(cfg)
            stage_names = [stage.name for stage in stages]
            selected_stage = handle_interactive_resume(args, pipeline_state, stage_names)
            if selected_stage:
                cfg["resume_from"] = selected_stage
                logger.info(f"Interactive resume: Selected stage '{selected_stage}'")
            else:
                print("❌ No stage selected. Exiting.")
                sys.exit(0)
        except Exception as e:
            logger.error(f"Interactive resume failed: {e}")
            sys.exit(1)
