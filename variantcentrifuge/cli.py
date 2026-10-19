"""Command-line interface for VariantCentrifuge."""

import argparse
import datetime
import logging
import sys
from pathlib import Path
from typing import Any

from .cli_config import (
    _calculate_sort_memory_limit,
    load_fields_from_file,
    populate_cli_config,
    resolve_presets_and_tumor_normal,
)
from .cli_parser_groups import (
    add_analysis_argument_groups,
    add_operational_argument_groups,
    get_curated_examples_epilog,
    populate_advanced_group,
    populate_misc_group,
    populate_performance_group,
)
from .cli_subcommands import handle_listing_and_resume_options
from .cli_sync import sync_cli_and_config
from .config import load_config
from .pipeline import run_pipeline
from .validators import (
    validate_cli_args,
    validate_mandatory_parameters,
    validate_phenotype_file,
    validate_required_tools,
    validate_vcf_file,
)
from .vcf_header_parser import resolve_tumor_normal_indices
from .version import __version__

__all__ = [
    "_calculate_sort_memory_limit",
    "create_parser",
    "main",
    "parse_args",
    "resolve_tumor_normal_indices",
    "validate_cli_args",
]

logger = logging.getLogger("variantcentrifuge")


class _VariantCentrifugeArgumentParser(argparse.ArgumentParser):
    """Parser that enforces -v/--vcf-file for execution while permitting listing flags without it."""

    def parse_args(  # type: ignore[override]
        self,
        args: Any = None,
        namespace: argparse.Namespace | None = None,
    ) -> argparse.Namespace:
        parsed: argparse.Namespace = super().parse_args(args=args, namespace=namespace)
        if not parsed.vcf_file and not (
            getattr(parsed, "list_stages", False)
            or getattr(parsed, "list_checkpoints", False)
            or getattr(parsed, "show_checkpoint_status", False)
        ):
            self.error("the following arguments are required: -v/--vcf-file")
        return parsed


def create_parser() -> argparse.ArgumentParser:
    """Create and return the argument parser for variantcentrifuge CLI."""
    parser = _VariantCentrifugeArgumentParser(
        description="variantcentrifuge: Filter and process VCF files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=get_curated_examples_epilog(),
    )

    # General Options
    general_group = parser.add_argument_group("General Options")
    general_group.add_argument(
        "--version",
        action="version",
        version=f"variantcentrifuge {__version__}",
        help="Show the current version and exit",
    )
    general_group.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARN", "WARNING", "ERROR"],
        default="INFO",
        help="Set the logging level",
    )
    general_group.add_argument(
        "--log-file", help="Path to a file to write logs to (in addition to stderr)."
    )
    general_group.add_argument(
        "-c",
        "--config",
        help="Path to configuration file",
        default=None,
    )

    # Core Input/Output
    io_group = parser.add_argument_group("Core Input/Output")
    io_group.add_argument(
        "-v", "--vcf-file", help="Input VCF file path", required=False, default=None
    )
    io_group.add_argument(
        "-o",
        "--output-file",
        nargs="?",
        const="stdout",
        help="Final output file name or 'stdout'/'-' for stdout",
    )
    io_group.add_argument(
        "--output-dir",
        help="Directory to store intermediate and final output files",
        default="output",
    )
    io_group.add_argument("--xlsx", action="store_true", help="Convert final result to Excel")
    io_group.add_argument(
        "--keep-intermediates",
        action="store_true",
        default=False,
        help="Keep intermediate files (by default, they are deleted after successful pipeline completion).",
    )
    io_group.add_argument(
        "--archive-results",
        action="store_true",
        default=False,
        help="Create a compressed tar.gz archive of the entire results directory after pipeline completion. "
        "Archive will be timestamped and placed in the parent directory of the output folder.",
    )
    io_group.add_argument(
        "--fail-on-empty", action="store_true", help="Exit with code 2 if zero variants remain."
    )

    add_analysis_argument_groups(parser)
    add_operational_argument_groups(parser)

    # Performance & Processing
    performance_group = parser.add_argument_group("Performance & Processing")
    populate_performance_group(performance_group)

    # Advanced & Tuning Options
    advanced_group = parser.add_argument_group("Advanced & Tuning Options")
    populate_advanced_group(advanced_group)

    # Miscellaneous Options
    misc_group = parser.add_argument_group("Miscellaneous Options")
    populate_misc_group(misc_group)

    return parser


def parse_args(args_list: list[str] | None = None) -> argparse.Namespace:
    """Parse command line arguments.

    Parameters
    ----------
    args_list : list, optional
        List of arguments to parse. If None, uses sys.argv

    Returns
    -------
    argparse.Namespace
        Parsed arguments
    """
    parser = create_parser()
    return parser.parse_args(args_list)


def main() -> int:
    """Run main entry point for variantcentrifuge CLI."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )

    start_time: datetime.datetime = datetime.datetime.now()
    logger.info(f"Run started at {start_time.isoformat()}")

    parser = create_parser()

    # For --list-field-profiles, we don't need other required arguments
    if "--list-field-profiles" in sys.argv:
        from .field_profile import list_profiles

        # Parse only the config argument (optional)
        profile_parser = argparse.ArgumentParser(description="List field profiles")
        profile_parser.add_argument("--list-field-profiles", action="store_true")
        profile_parser.add_argument(
            "-c", "--config", help="Path to configuration file", default=None
        )
        profile_args = profile_parser.parse_args()

        profile_cfg = load_config(profile_args.config)
        profiles = list_profiles(profile_cfg)
        if profiles:
            print("Available field profiles:")
            for prof in profiles:
                print(f"  {prof['name']:12s}  {prof['description']}")
        else:
            print("No field profiles defined in configuration.")
        sys.exit(0)

    # For --show-vcf-annotations, we only need the VCF file
    if "--show-vcf-annotations" in sys.argv:
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
        ann_args, _ = ann_parser.parse_known_args()

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

    # For --show-checkpoint-status, we don't need other required arguments
    if "--show-checkpoint-status" in sys.argv:
        # Create a minimal parser just for checkpoint status
        status_parser = argparse.ArgumentParser(description="Show checkpoint status")
        status_parser.add_argument("--show-checkpoint-status", action="store_true")
        status_parser.add_argument("--output-dir", default="output", help="Output directory")
        status_parser.add_argument(
            "--log-level", choices=["DEBUG", "INFO", "WARN", "WARNING", "ERROR"], default="INFO"
        )
        status_parser.add_argument(
            "-c", "--config", help="Path to configuration file", default=None
        )
        status_args = status_parser.parse_args()

        # Configure logging
        log_level_map = {
            "DEBUG": logging.DEBUG,
            "INFO": logging.INFO,
            "WARN": logging.WARNING,
            "WARNING": logging.WARNING,
            "ERROR": logging.ERROR,
        }
        logging.getLogger("variantcentrifuge").setLevel(log_level_map[status_args.log_level])

        # Show checkpoint status
        from .checkpoint import PipelineState

        print("Checking checkpoint status for stage-based pipeline...")

        pipeline_state = PipelineState(status_args.output_dir)
        if pipeline_state.load():
            print(pipeline_state.get_summary())
        else:
            print(f"No checkpoint state found in {status_args.output_dir}")
        sys.exit(0)

    args: argparse.Namespace = parser.parse_args()

    # Configure logging level
    log_level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARN": logging.WARNING,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
    }
    logging.getLogger("variantcentrifuge").setLevel(log_level_map[args.log_level])

    # Print the actual command line that was used to run the program in debug mode
    if args.log_level == "DEBUG":
        import shlex

        # Reconstruct the command line that was used
        command_line = " ".join(shlex.quote(arg) for arg in sys.argv)
        logger.debug(f"Command line invocation: {command_line}")
        logger.debug(f"Python executable: {sys.executable}")
        logger.debug(f"Working directory: {Path.cwd()}")

    # If a log file is specified, add a file handler
    if args.log_file:
        # Ensure the log file directory exists
        log_file_path = Path(args.log_file)
        log_file_path.parent.mkdir(parents=True, exist_ok=True)

        fh = logging.FileHandler(args.log_file)
        fh.setLevel(log_level_map[args.log_level])
        fh.setFormatter(
            logging.Formatter(
                "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
            )
        )
        logger.addHandler(fh)
        logger.debug(f"Logging to file enabled: {args.log_file}")

    logger.debug(f"CLI arguments: {args}")

    # Load configuration
    cfg: dict[str, Any] = load_config(args.config)
    logger.debug(f"Configuration loaded: {cfg}")

    # Resolve field profile (expands {{fragment:param}} templates in presets,
    # merges profile-specific fields and hidden columns into config)
    from .field_profile import resolve_profile

    if args.field_profile:
        cfg["field_profile"] = args.field_profile
    resolve_profile(cfg)

    # Synchronize configuration and CLI options so listing and validation see effective settings
    sync_cli_and_config(args, cfg)

    # Handle information and interactive options that exit early (e.g. --list-stages, --list-checkpoints)
    handle_listing_and_resume_options(args, cfg)

    # Ordinary execution requires VCF input
    if not args.vcf_file:
        parser.error("the following arguments are required: -v/--vcf-file")

    reference: str | None = args.reference or cfg.get("reference")
    filters: str | None = args.filters or cfg.get("filters")
    fields_file = getattr(args, "fields_file", None) or cfg.get("fields_file")

    try:
        validate_cli_args(args, cfg, parser=parser)
        raw_fields: Any = (
            load_fields_from_file(fields_file)
            if (not args.fields and fields_file)
            else (args.fields or cfg.get("fields_to_extract"))
        )
        fields: str | None = " ".join(raw_fields) if isinstance(raw_fields, list) else raw_fields

        filters, tumor_normal_vars = resolve_presets_and_tumor_normal(args, cfg, filters)

        populate_cli_config(args, cfg, parser, reference, filters, fields, tumor_normal_vars)

        # Validate mandatory parameters, external tools, VCF and phenotype inputs
        validate_mandatory_parameters(reference, filters, fields)
        validate_required_tools(args, cfg, logger)
        validate_vcf_file(args.vcf_file, logger)
        if args.phenotype_file:
            validate_phenotype_file(
                args.phenotype_file,
                args.phenotype_sample_column,
                args.phenotype_value_column,
                logger,
            )
    except (ValueError, FileNotFoundError, RuntimeError) as e:
        logger.error(str(e))
        sys.exit(1)

    # Performance configuration
    cfg["max_memory_gb"] = args.max_memory_gb
    cfg["genotype_replacement_method"] = args.genotype_replacement_method

    # Debug logging for performance parameters
    logger.debug(f"Performance config mapping: max_memory_gb={args.max_memory_gb}")
    logger.debug(
        f"Performance config mapping: genotype_replacement_method={args.genotype_replacement_method}"
    )

    # Field extraction configuration
    cfg["add_column"] = args.add_column if args.add_column else cfg.get("add_column")
    cfg["no_replacement"] = bool(
        args.no_replacement if args.no_replacement else cfg.get("no_replacement", False)
    )

    # Phenotype and sample group configuration
    cfg["case_phenotypes"] = (
        args.case_phenotypes if args.case_phenotypes is not None else cfg.get("case_phenotypes")
    )
    cfg["case_phenotypes_file"] = (
        args.case_phenotypes_file
        if args.case_phenotypes_file is not None
        else cfg.get("case_phenotypes_file")
    )
    cfg["case_samples"] = (
        args.case_samples if args.case_samples is not None else cfg.get("case_samples")
    )
    cfg["case_samples_file"] = (
        args.case_samples_file
        if args.case_samples_file is not None
        else cfg.get("case_samples_file")
    )
    cfg["control_phenotypes"] = (
        args.control_phenotypes
        if args.control_phenotypes is not None
        else cfg.get("control_phenotypes")
    )
    cfg["control_phenotypes_file"] = (
        args.control_phenotypes_file
        if args.control_phenotypes_file is not None
        else cfg.get("control_phenotypes_file")
    )
    cfg["control_samples"] = (
        args.control_samples if args.control_samples is not None else cfg.get("control_samples")
    )
    cfg["control_samples_file"] = (
        args.control_samples_file
        if args.control_samples_file is not None
        else cfg.get("control_samples_file")
    )

    # Statistical analysis configuration
    cfg["stats_output_file"] = (
        args.stats_output_file
        if args.stats_output_file is not None
        else cfg.get("stats_output_file")
    )

    # Pseudonymization configuration
    cfg["pseudonymize"] = bool(
        args.pseudonymize if args.pseudonymize else cfg.get("pseudonymize", False)
    )
    cfg["pseudonymize_schema"] = args.pseudonymize_schema
    cfg["pseudonymize_prefix"] = args.pseudonymize_prefix
    cfg["pseudonymize_pattern"] = args.pseudonymize_pattern
    cfg["pseudonymize_category_field"] = args.pseudonymize_category_field
    cfg["pseudonymize_table"] = (
        args.pseudonymize_table
        if args.pseudonymize_table is not None
        else cfg.get("pseudonymize_table")
    )
    cfg["pseudonymize_ped"] = (
        args.pseudonymize_ped if args.pseudonymize_ped is not None else cfg.get("pseudonymize_ped")
    )

    # Validate pseudonymization arguments
    if cfg["pseudonymize"] and not cfg.get("pseudonymize_table"):
        logger.error("--pseudonymize-table is required when using --pseudonymize")
        sys.exit(1)
    if args.pseudonymize_table and not cfg["pseudonymize"]:
        logger.warning("--pseudonymize-table provided without --pseudonymize. It will be ignored.")
    if args.pseudonymize_ped and not cfg["pseudonymize"]:
        logger.warning("--pseudonymize-ped provided without --pseudonymize. It will be ignored.")
    if args.pseudonymize_pattern and args.pseudonymize_schema != "custom":
        logger.warning(
            f"--pseudonymize-pattern is only used with --pseudonymize-schema custom. "
            f"Current schema: {args.pseudonymize_schema}"
        )

    try:
        # Create a new args object with the config properly set
        from types import SimpleNamespace

        refactored_args = SimpleNamespace(**vars(args))
        refactored_args.config = cfg
        refactored_args.ped_file = (
            cfg.get("ped_file") or getattr(args, "ped", None) or getattr(args, "ped_file", None)
        )
        refactored_args.xlsx = cfg.get("xlsx", False)
        refactored_args.keep_intermediates = cfg.get("keep_intermediates", False)
        refactored_args.archive_results = cfg.get("archive_results", False)
        refactored_args.gzip_intermediates = cfg.get("gzip_intermediates", True)
        refactored_args.gene_name = cfg.get("gene_name")
        refactored_args.gene_file = cfg.get("gene_file")
        refactored_args.html_report = cfg.get("html_report", False)
        refactored_args.pseudonymize = cfg.get("pseudonymize", False)
        refactored_args.no_replacement = cfg.get("no_replacement", False)
        refactored_args.late_filtering = cfg.get("late_filtering", False)
        refactored_args.final_filter = cfg.get("final_filter")
        refactored_args.perform_gene_burden = cfg.get("perform_gene_burden", False)
        refactored_args.perform_association = cfg.get("perform_association", False)
        refactored_args.no_stats = cfg.get("no_stats", False)
        refactored_args.fields_file = cfg.get("fields_file")
        refactored_args.igv_flanking = cfg.get("igv_flanking", 50)
        refactored_args.phenotype_file = cfg.get("phenotype_file")
        if not hasattr(refactored_args, "start_time"):
            refactored_args.start_time = start_time

        run_pipeline(refactored_args)  # type: ignore[arg-type]
        return 0
    except SystemExit as e:
        return int(e.code) if e.code is not None else 1
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        return 1


if __name__ == "__main__":
    main()
