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
from .cli_parser_groups import add_analysis_argument_groups, add_operational_argument_groups
from .cli_subcommands import handle_listing_and_resume_options
from .config import load_config
from .pipeline import run_pipeline
from .validators import (
    validate_cli_args,
    validate_mandatory_parameters,
    validate_phenotype_file,
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


def create_parser() -> argparse.ArgumentParser:
    """Create and return the argument parser for variantcentrifuge CLI."""
    parser = argparse.ArgumentParser(description="variantcentrifuge: Filter and process VCF files.")

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
    io_group.add_argument("-v", "--vcf-file", help="Input VCF file path", required=True)
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
    performance_group.add_argument(
        "--threads",
        default="auto",
        help="Number of parallel worker processes to use for variant extraction and filtering. "
        "Use 'auto' to detect available CPU cores (default: auto). "
        "Set to 1 to disable parallelization.",
    )
    performance_group.add_argument(
        "--no-chunked-processing",
        action="store_true",
        help="Disable chunked processing even for large files. "
        "May cause memory issues with very large datasets.",
    )
    performance_group.add_argument(
        "--force-chunked-processing",
        action="store_true",
        help="Force chunked processing even for small files. "
        "Useful for testing or memory-constrained environments.",
    )
    performance_group.add_argument(
        "--sort-memory-limit",
        type=str,
        default=None,
        help="Memory limit for external sort command (default: auto - scales with available memory). "
        "Examples: 500M, 4G, 8G, auto. When 'auto', uses ~10%% of available memory.",
    )
    performance_group.add_argument(
        "--sort-parallel",
        type=int,
        default=4,
        help="Number of parallel threads for sorting (default: 4)",
    )
    performance_group.add_argument(
        "--genotype-replacement-method",
        choices=[
            "auto",
            "sequential",
            "vectorized",
            "chunked-vectorized",
            "parallel",
            "streaming-parallel",
        ],
        default="auto",
        help="Method for genotype replacement processing. "
        "auto: automatically select based on data characteristics and memory availability; "
        "sequential: line-by-line streaming (memory efficient); "
        "vectorized: pandas-based vectorized operations (faster for many samples); "
        "chunked-vectorized: memory-safe vectorized processing for large files; "
        "parallel: multi-threaded chunked processing (for very large files); "
        "streaming-parallel: enhanced streaming pipeline with optimal thread utilization. "
        "Default: auto",
    )
    performance_group.add_argument(
        "--max-memory-gb",
        type=float,
        help="Maximum memory in GB to use for genotype processing (auto-detected if not specified). "
        "Used to intelligently select processing method and chunk sizes.",
    )
    performance_group.add_argument(
        "--force-inheritance-processing",
        action="store_true",
        help="Force inheritance analysis even if it exceeds safe memory limits. "
        "Use with caution on systems with limited memory.",
    )
    performance_group.add_argument(
        "--sample-column-creation-method",
        choices=["iterative", "vectorized", "auto"],
        default="auto",
        help="Method for creating sample columns from GT field during inheritance analysis. "
        "iterative: traditional row-by-row regex parsing (slower but more reliable); "
        "vectorized: pandas vectorized string operations (5-10x faster for large datasets); "
        "auto: automatically select based on dataset size and complexity. "
        "Default: auto",
    )
    performance_group.add_argument(
        "--memory-safety-factor",
        type=float,
        default=0.92,
        help="Fraction of allocated memory to use safely (default: 0.92 = 92%%). "
        "Lower values are more conservative, higher values use more memory.",
    )
    performance_group.add_argument(
        "--inheritance-memory-fraction",
        type=float,
        default=0.85,
        help="Fraction of safe memory to allocate for inheritance analysis (default: 0.85 = 85%%). "
        "Remainder is reserved for other pipeline stages.",
    )

    # Miscellaneous Options
    misc_group = parser.add_argument_group("Miscellaneous Options")
    misc_group.add_argument(
        "--no-gzip-intermediates",
        dest="gzip_intermediates",
        action="store_false",
        default=True,
        help="Disable compression of intermediate TSV files (compression is enabled by default).",
    )
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

        # Validate mandatory parameters
        validate_mandatory_parameters(reference, filters, fields)
        validate_vcf_file(args.vcf_file, logger)
        if args.phenotype_file:
            validate_phenotype_file(
                args.phenotype_file,
                args.phenotype_sample_column,
                args.phenotype_value_column,
                logger,
            )

        populate_cli_config(args, cfg, parser, reference, filters, fields, tumor_normal_vars)
    except (ValueError, FileNotFoundError) as e:
        logger.error(str(e))
        sys.exit(1)

    # Handle information and interactive options that exit early
    handle_listing_and_resume_options(args, cfg)

    # Performance configuration
    cfg["max_memory_gb"] = args.max_memory_gb
    cfg["genotype_replacement_method"] = args.genotype_replacement_method

    # Debug logging for performance parameters
    logger.debug(f"Performance config mapping: max_memory_gb={args.max_memory_gb}")
    logger.debug(
        f"Performance config mapping: genotype_replacement_method={args.genotype_replacement_method}"
    )

    # Core I/O configuration: preserve config values if CLI flag not passed
    cfg["xlsx"] = bool(args.xlsx if args.xlsx else cfg.get("xlsx", False))
    cfg["keep_intermediates"] = bool(
        args.keep_intermediates if args.keep_intermediates else cfg.get("keep_intermediates", False)
    )
    cfg["archive_results"] = bool(
        args.archive_results if args.archive_results else cfg.get("archive_results", False)
    )
    cfg["gzip_intermediates"] = (
        args.gzip_intermediates
        if getattr(args, "gzip_intermediates", None) is not None
        else cfg.get("gzip_intermediates", True)
    )
    logger.debug(f"I/O config mapping: xlsx={cfg['xlsx']}")

    # Gene selection configuration
    cfg["gene_name"] = args.gene_name if args.gene_name is not None else cfg.get("gene_name")
    cfg["gene_file"] = args.gene_file if args.gene_file is not None else cfg.get("gene_file")
    logger.debug(f"Gene config mapping: gene_name={cfg['gene_name']}")

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

    # Reporting configuration
    cfg["html_report"] = bool(
        args.html_report if args.html_report else cfg.get("html_report", False)
    )

    # Pseudonymization configuration
    cfg["pseudonymize"] = bool(
        args.pseudonymize if args.pseudonymize else cfg.get("pseudonymize", False)
    )
    cfg["pseudonymize_schema"] = args.pseudonymize_schema
    cfg["pseudonymize_prefix"] = args.pseudonymize_prefix
    cfg["pseudonymize_pattern"] = args.pseudonymize_pattern
    cfg["pseudonymize_category_field"] = args.pseudonymize_category_field
    cfg["pseudonymize_table"] = args.pseudonymize_table
    cfg["pseudonymize_ped"] = args.pseudonymize_ped

    # Validate pseudonymization arguments
    if cfg["pseudonymize"] and not args.pseudonymize_table:
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
        refactored_args.xlsx = cfg["xlsx"]
        refactored_args.keep_intermediates = cfg["keep_intermediates"]
        refactored_args.archive_results = cfg["archive_results"]
        refactored_args.gzip_intermediates = cfg["gzip_intermediates"]
        refactored_args.gene_name = cfg["gene_name"]
        refactored_args.gene_file = cfg["gene_file"]
        refactored_args.html_report = cfg["html_report"]
        refactored_args.pseudonymize = cfg["pseudonymize"]
        refactored_args.no_replacement = cfg["no_replacement"]
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
