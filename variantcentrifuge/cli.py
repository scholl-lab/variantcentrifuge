"""Command-line interface for VariantCentrifuge."""

import argparse
import datetime
import logging
import sys
from pathlib import Path
from typing import Any

from .cli_config import (
    _calculate_sort_memory_limit,
    populate_cli_config,
    resolve_presets_and_tumor_normal,
)
from .cli_parser_groups import add_analysis_argument_groups
from .config import load_config
from .pipeline import run_pipeline
from .validators import (
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

    # Scoring & Custom Annotations
    scoring_group = parser.add_argument_group("Scoring & Custom Annotations")
    scoring_group.add_argument(
        "--scoring-config-path",
        help="Path to a directory containing a scoring model "
        "(variable_assignment_config.json and formula_config.json).",
    )
    scoring_group.add_argument(
        "--annotate-bed",
        action="append",
        default=[],
        metavar="BED_FILE",
        help="Path to a BED file for annotation. Adds 'Region=<name>' to the Custom_Annotation column for overlapping variants. Can be used multiple times.",
    )
    scoring_group.add_argument(
        "--annotate-gene-list",
        action="append",
        default=[],
        metavar="GENE_LIST_FILE_PATH",
        help="Path to a gene list file (one gene per line, case-insensitive matching). "
        "A new column will be added to the output TSV for each file, "
        "named after a sanitized version of the file's basename (e.g., 'my_gene_set' from 'my_gene_set.txt'), "
        "indicating if the variant's GENE (or one of its comma-separated genes) is in this list ('yes'/'no'). "
        "Specify multiple times for multiple lists.",
    )
    scoring_group.add_argument(
        "--annotate-json-genes",
        action="append",
        default=[],
        metavar="JSON_FILE",
        help="Path to a JSON file with structured gene data. Requires --json-gene-mapping.",
    )
    scoring_group.add_argument(
        "--json-gene-mapping",
        help='JSON string to map fields in JSON gene files. e.g., \'{"identifier":"gene_symbol","dataFields":["panel","inheritance"]}\'',
    )
    scoring_group.add_argument(
        "--json-genes-as-columns",
        action="store_true",
        help="When using --annotate-json-genes, add each specified dataField as a separate column "
        "instead of appending key-value pairs to the Custom_Annotation column.",
    )

    # Reporting & Visualization
    report_group = parser.add_argument_group("Reporting & Visualization")
    report_group.add_argument(
        "--html-report",
        action="store_true",
        help="Generate an interactive HTML report with sortable variant tables and summary plots.",
    )
    report_group.add_argument(
        "--igv",
        action="store_true",
        help="Enable IGV.js integration for genomic visualization.",
    )
    report_group.add_argument(
        "--bam-mapping-file",
        help="Path to a TSV or CSV file mapping sample IDs to BAM files (sample_id,bam_path).",
    )
    report_group.add_argument(
        "--igv-reference",
        help="Genome reference identifier for IGV (e.g., 'hg19' or 'hg38'). Required if --igv is enabled unless --igv-fasta is provided.",
    )
    report_group.add_argument(
        "--igv-fasta",
        help="Path to a local FASTA file for IGV reports. The index file (.fai) must exist in the same location with the same name (e.g., reference.fa.fai). This will be used instead of --igv-reference if both are provided.",
    )
    report_group.add_argument(
        "--igv-ideogram",
        help="Path to an ideogram file for chromosome visualization in IGV reports.",
    )
    report_group.add_argument(
        "--igv-flanking",
        type=int,
        default=None,
        help="Flanking region size in base pairs for IGV reports (default: 50)",
    )
    report_group.add_argument(
        "--igv-max-allele-len-filename",
        type=int,
        default=10,
        help="Maximum length for REF/ALT alleles in IGV report filenames (default: 10). Longer alleles will be truncated and hashed.",
    )
    report_group.add_argument(
        "--igv-hash-len-filename",
        type=int,
        default=6,
        help="Length of hash to append when truncating alleles in filenames (default: 6).",
    )
    report_group.add_argument(
        "--igv-max-variant-part-filename",
        type=int,
        default=50,
        help="Maximum length for the variant part of IGV report filenames",
    )

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
        default="auto",
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

    # Checkpoint & Resume Options
    checkpoint_group = parser.add_argument_group("Checkpoint & Resume Options")
    checkpoint_group.add_argument(
        "--enable-checkpoint",
        action="store_true",
        help="Enable checkpoint system to track pipeline progress and allow resumption",
    )
    checkpoint_group.add_argument(
        "--resume",
        action="store_true",
        help="Resume pipeline from the last successful checkpoint. "
        "Requires --enable-checkpoint and previous run in same output directory",
    )
    checkpoint_group.add_argument(
        "--checkpoint-checksum",
        action="store_true",
        help="Calculate file checksums for checkpoint validation (RECOMMENDED FOR PRODUCTION). "
        "This provides the most reliable validation of file integrity when resuming from checkpoints. "
        "Without this option, interrupted stages will always be re-executed for safety. "
        "The performance overhead is often worth the reliability gain in critical pipelines. "
        "Default is to use file size and modification time only",
    )
    checkpoint_group.add_argument(
        "--show-checkpoint-status",
        action="store_true",
        help="Show checkpoint status for the output directory and exit",
    )
    checkpoint_group.add_argument(
        "--resume-from",
        type=str,
        metavar="STAGE_NAME",
        help="Restart pipeline from a specific stage, re-executing that stage and all subsequent stages. "
        "Stages before the specified stage remain completed. Use --list-stages to see available stages. "
        "Requires --enable-checkpoint",
    )
    checkpoint_group.add_argument(
        "--list-stages",
        action="store_true",
        help="List all available stages for current configuration and exit",
    )
    checkpoint_group.add_argument(
        "--list-checkpoints",
        action="store_true",
        help="List all completed stages from checkpoint file and exit. "
        "Requires existing checkpoint file in output directory",
    )
    checkpoint_group.add_argument(
        "--interactive-resume",
        action="store_true",
        help="Interactive mode to select resume point from completed stages. "
        "Requires --enable-checkpoint and existing checkpoint file",
    )

    # Data Privacy Options
    privacy_group = parser.add_argument_group("Data Privacy Options")
    privacy_group.add_argument(
        "--pseudonymize",
        action="store_true",
        help="Enable sample pseudonymization for privacy-preserving output",
    )
    privacy_group.add_argument(
        "--pseudonymize-schema",
        type=str,
        choices=["sequential", "categorical", "anonymous", "custom"],
        default="sequential",
        help="Pseudonymization naming schema (default: sequential)",
    )
    privacy_group.add_argument(
        "--pseudonymize-prefix",
        type=str,
        default="SAMPLE",
        help="Prefix for sequential schema (default: SAMPLE)",
    )
    privacy_group.add_argument(
        "--pseudonymize-pattern",
        type=str,
        help="Custom pattern for 'custom' schema (e.g., '{prefix}_{category}_{index:03d}')",
    )
    privacy_group.add_argument(
        "--pseudonymize-category-field",
        type=str,
        default="phenotype",
        help="Metadata field for categorical schema (default: phenotype)",
    )
    privacy_group.add_argument(
        "--pseudonymize-table",
        type=str,
        help="Path to save pseudonymization mapping table (required if --pseudonymize)",
    )
    privacy_group.add_argument(
        "--pseudonymize-ped",
        action="store_true",
        help="Also create pseudonymized PED file for sharing",
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


def parse_args(args_list=None):
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
    """Run main entry point for variantcentrifuge CLI.

    Steps:
        1. Parse arguments.
        2. Configure logging and load config.
        3. Validate mandatory parameters (reference, filters, fields).
        4. Validate input files (VCF, phenotype).
        5. Update configuration with CLI parameters.
        6. Run the pipeline.

    Changes for presets:
        - Added a --preset argument which allows specifying one or more predefined filters.
        - If multiple presets are chosen, they are combined with AND (&).
        - If user-specified filters are also provided, they are combined with these presets using AND.
        - Presets must be defined in config.json under "presets".

    Changes for sample substring removal:
        - Added a --remove-sample-substring argument which, if provided, removes the specified substring
          from all sample names extracted from the VCF before any comparisons or mappings.

    Changes for links:
        - Added a --no-links argument to disable adding link columns. By default, links are added.

    Changes for custom columns:
        - Added a --add-column argument which can be given multiple times to append columns with given
          headers (but blank values) to the final output.

    Changes for transcript filtering:
        - Added a --transcript-list and --transcript-file argument which allows specifying one or more
          transcript IDs (by comma-separated list or a file with one transcript ID per line).
        - If provided, these transcripts are used to construct a SnpSift filter expression that
          filters variants by EFF[*].TRID field.

    New changes for extra sample fields in genotype:
        - Replaced the old two-flag logic with a single --append-extra-sample-fields flag that
          optionally takes zero or more fields. If this flag is used with no arguments, it simply
          enables genotype+field appending with no extra fields. If arguments (space-separated) are
          provided, those fields are appended (e.g. DP, AD).
        - We retain --extra-sample-field-delimiter to control the delimiter between genotype and
          extra fields (default ':').
    """
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
    fields: str | None = args.fields or cfg.get("fields_to_extract")

    filters, tumor_normal_vars = resolve_presets_and_tumor_normal(args, cfg, filters)

    # Validate mandatory parameters
    validate_mandatory_parameters(reference, filters, fields)
    validate_vcf_file(args.vcf_file, logger)
    tx_file = args.transcript_file or cfg.get("transcript_file")
    if tx_file and not Path(tx_file).is_file():
        logger.error(f"Transcript file not found: {tx_file}")
        sys.exit(1)

    if args.phenotype_file and (
        args.phenotype_sample_column is None or args.phenotype_value_column is None
    ):
        logger.error(
            "If a phenotype file is provided, both --phenotype-sample-column "
            "and --phenotype-value-column must be specified."
        )
        sys.exit(1)
    if args.phenotype_file:
        validate_phenotype_file(
            args.phenotype_file,
            args.phenotype_sample_column,
            args.phenotype_value_column,
            logger,
        )

    populate_cli_config(args, cfg, parser, reference, filters, fields, tumor_normal_vars)

    # Note: --show-checkpoint-status is handled earlier before argument parsing

    # Validate checkpoint arguments
    if args.resume and not args.enable_checkpoint:
        logger.error("--resume requires --enable-checkpoint to be set")
        sys.exit(1)

    # Validate selective resume arguments
    if args.resume_from and not args.enable_checkpoint:
        logger.error("--resume-from requires --enable-checkpoint to be set")
        sys.exit(1)

    if args.interactive_resume and not args.enable_checkpoint:
        logger.error("--interactive-resume requires --enable-checkpoint to be set")
        sys.exit(1)

    # Mutually exclusive resume options
    resume_options = [args.resume, args.resume_from, args.interactive_resume]
    if sum(bool(x) for x in resume_options) > 1:
        logger.error(
            "Cannot use multiple resume options simultaneously: --resume, --resume-from, --interactive-resume"
        )
        sys.exit(1)

    # Handle information and interactive options that exit early
    if args.list_stages:
        from .display_utils import display_available_stages
        from .stages.stage_registry import initialize_registry

        initialize_registry()

        # Show available stages from registry
        from .pipeline import create_stages_from_config

        try:
            stages = create_stages_from_config(cfg)
            display_available_stages(stages, cfg)
        except Exception as e:
            logger.error(f"Failed to create stages: {e}")
            sys.exit(1)

        sys.exit(0)

    if args.list_checkpoints:
        from .checkpoint import PipelineState
        from .display_utils import display_enhanced_status

        pipeline_state = PipelineState(args.output_dir)
        if pipeline_state.load():
            # For list checkpoints, show enhanced status instead of basic summary
            display_enhanced_status(pipeline_state, [])
        else:
            print(f"❌ No checkpoint file found in {args.output_dir}")
            print("   Start a new pipeline run with --enable-checkpoint to create checkpoints.")

        sys.exit(0)

    if args.interactive_resume:
        from .checkpoint import PipelineState
        from .interactive_resume import handle_interactive_resume
        from .stages.stage_registry import initialize_registry

        initialize_registry()

        pipeline_state = PipelineState(args.output_dir)
        if not pipeline_state.load():
            print(f"❌ No checkpoint file found in {args.output_dir}")
            print("   Start a new pipeline run with --enable-checkpoint to create checkpoints.")
            sys.exit(1)

        # Get available stages for current configuration
        from .pipeline import create_stages_from_config

        try:
            stages = create_stages_from_config(cfg)
            stage_names = [stage.name for stage in stages]

            selected_stage = handle_interactive_resume(args, pipeline_state, stage_names)
            if selected_stage:
                # Update config with selected resume point
                cfg["resume_from"] = selected_stage
                logger.info(f"Interactive resume: Selected stage '{selected_stage}'")
            else:
                print("❌ No stage selected. Exiting.")
                sys.exit(0)
        except Exception as e:
            logger.error(f"Interactive resume failed: {e}")
            sys.exit(1)

    # Performance configuration
    cfg["max_memory_gb"] = args.max_memory_gb
    cfg["genotype_replacement_method"] = args.genotype_replacement_method

    # Debug logging for performance parameters
    logger.debug(f"Performance config mapping: max_memory_gb={args.max_memory_gb}")
    logger.debug(
        f"Performance config mapping: genotype_replacement_method={args.genotype_replacement_method}"
    )

    # Core I/O configuration
    cfg["xlsx"] = args.xlsx
    cfg["keep_intermediates"] = args.keep_intermediates
    cfg["archive_results"] = args.archive_results
    cfg["gzip_intermediates"] = args.gzip_intermediates

    # Debug logging for I/O parameters
    logger.debug(f"I/O config mapping: xlsx={args.xlsx}")
    logger.debug(f"I/O config mapping: keep_intermediates={args.keep_intermediates}")
    logger.debug(f"I/O config mapping: archive_results={args.archive_results}")
    logger.debug(f"I/O config mapping: gzip_intermediates={args.gzip_intermediates}")

    # Gene selection configuration
    cfg["gene_name"] = args.gene_name
    cfg["gene_file"] = args.gene_file

    # Debug logging for gene selection parameters
    logger.debug(f"Gene config mapping: gene_name={args.gene_name}")
    logger.debug(f"Gene config mapping: gene_file={args.gene_file}")

    # Field extraction configuration
    cfg["add_column"] = args.add_column
    cfg["no_replacement"] = args.no_replacement

    # Debug logging for field extraction parameters
    logger.debug(f"Field config mapping: add_column={args.add_column}")
    logger.debug(f"Field config mapping: no_replacement={args.no_replacement}")

    # Phenotype and sample group configuration
    cfg["case_phenotypes"] = args.case_phenotypes
    cfg["case_phenotypes_file"] = args.case_phenotypes_file
    cfg["case_samples"] = args.case_samples
    cfg["case_samples_file"] = args.case_samples_file
    cfg["control_phenotypes"] = args.control_phenotypes
    cfg["control_phenotypes_file"] = args.control_phenotypes_file
    cfg["control_samples"] = args.control_samples
    cfg["control_samples_file"] = args.control_samples_file

    # Debug logging for phenotype and sample group parameters
    logger.debug(f"Phenotype config mapping: case_phenotypes={args.case_phenotypes}")
    logger.debug(f"Phenotype config mapping: case_phenotypes_file={args.case_phenotypes_file}")
    logger.debug(f"Phenotype config mapping: case_samples={args.case_samples}")
    logger.debug(f"Phenotype config mapping: case_samples_file={args.case_samples_file}")
    logger.debug(f"Phenotype config mapping: control_phenotypes={args.control_phenotypes}")
    logger.debug(
        f"Phenotype config mapping: control_phenotypes_file={args.control_phenotypes_file}"
    )
    logger.debug(f"Phenotype config mapping: control_samples={args.control_samples}")
    logger.debug(f"Phenotype config mapping: control_samples_file={args.control_samples_file}")

    # Statistical analysis configuration
    cfg["stats_output_file"] = args.stats_output_file

    # Debug logging for statistical analysis parameters
    logger.debug(f"Stats config mapping: stats_output_file={args.stats_output_file}")

    # Reporting configuration
    cfg["html_report"] = args.html_report

    # Debug logging for reporting parameters
    logger.debug(f"Reporting config mapping: html_report={args.html_report}")

    # Pseudonymization configuration
    cfg["pseudonymize"] = args.pseudonymize
    cfg["pseudonymize_schema"] = args.pseudonymize_schema
    cfg["pseudonymize_prefix"] = args.pseudonymize_prefix
    cfg["pseudonymize_pattern"] = args.pseudonymize_pattern
    cfg["pseudonymize_category_field"] = args.pseudonymize_category_field
    cfg["pseudonymize_table"] = args.pseudonymize_table
    cfg["pseudonymize_ped"] = args.pseudonymize_ped

    # Validate pseudonymization arguments
    if args.pseudonymize and not args.pseudonymize_table:
        logger.error("--pseudonymize-table is required when using --pseudonymize")
        sys.exit(1)
    if args.pseudonymize_table and not args.pseudonymize:
        logger.warning("--pseudonymize-table provided without --pseudonymize. It will be ignored.")
    if args.pseudonymize_ped and not args.pseudonymize:
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
