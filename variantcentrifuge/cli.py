# File: variantcentrifuge/cli.py

import argparse
import datetime
import logging
import sys
from typing import Any, Dict, Optional

from .config import load_config
from .pipeline import run_pipeline
from .validators import validate_mandatory_parameters, validate_phenotype_file, validate_vcf_file
from .version import __version__

logger = logging.getLogger("variantcentrifuge")


def main() -> None:
    """
    Main entry point for variantcentrifuge CLI.

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

    parser = argparse.ArgumentParser(description="variantcentrifuge: Filter and process VCF files.")
    parser.add_argument(
        "--version",
        action="version",
        version=f"variantcentrifuge {__version__}",
        help="Show the current version and exit",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARN", "ERROR"],
        default="INFO",
        help="Set the logging level",
    )
    parser.add_argument(
        "--log-file", help="Path to a file to write logs to (in addition to stderr)."
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Path to configuration file",
        default=None,
    )
    parser.add_argument("-g", "--gene-name", help="Gene or list of genes of interest")
    parser.add_argument("-G", "--gene-file", help="File containing gene names, one per line")
    parser.add_argument("-v", "--vcf-file", help="Input VCF file path", required=True)
    parser.add_argument(
        "-r",
        "--reference",
        help="Reference database for snpEff",
        default=None,
    )
    parser.add_argument("-f", "--filters", help="Filters to apply in SnpSift filter")
    # Preset argument
    parser.add_argument(
        "--preset",
        action="append",
        help=(
            "Apply predefined filtering presets defined in the config file. "
            "Specify multiple times for multiple presets. They are combined with AND. "
            "If custom filters are also given, they are combined with these presets using AND."
        ),
    )
    parser.add_argument("-e", "--fields", help="Fields to extract with SnpSift extractFields")
    parser.add_argument(
        "-o",
        "--output-file",
        nargs="?",
        const="stdout",
        help="Final output file name or 'stdout'/'-' for stdout",
    )
    parser.add_argument("--xlsx", action="store_true", help="Convert final result to Excel")
    parser.add_argument(
        "--no-replacement", action="store_true", help="Skip genotype replacement step"
    )
    parser.add_argument("--stats-output-file", help="File to write analysis statistics")
    parser.add_argument(
        "--perform-gene-burden",
        action="store_true",
        help="Perform gene burden analysis",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory to store intermediate and final output files",
        default="output",
    )
    # MODIFIED: Start of intermediate cleanup feature
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        default=False,
        help="Keep intermediate files (by default, they are deleted after successful pipeline completion).",
    )
    # MODIFIED: End of intermediate cleanup feature
    parser.add_argument("--phenotype-file", help="Path to phenotype file (.csv or .tsv)")
    parser.add_argument(
        "--phenotype-sample-column",
        help="Name of the column containing sample IDs in phenotype file",
    )
    parser.add_argument(
        "--phenotype-value-column",
        help="Name of the column containing phenotype values in phenotype file",
    )
    parser.add_argument(
        "--no-stats",
        action="store_true",
        default=False,
        help="Skip the statistics computation step.",
    )
    parser.add_argument("--case-phenotypes", help="Comma-separated HPO terms defining case group")
    parser.add_argument(
        "--control-phenotypes", help="Comma-separated HPO terms defining control group"
    )
    parser.add_argument("--case-phenotypes-file", help="File with HPO terms for case group")
    parser.add_argument("--control-phenotypes-file", help="File with HPO terms for control group")
    parser.add_argument("--case-samples", help="Comma-separated sample IDs defining the case group")
    parser.add_argument(
        "--control-samples",
        help="Comma-separated sample IDs defining the control group",
    )
    parser.add_argument("--case-samples-file", help="File with sample IDs for case group")
    parser.add_argument("--control-samples-file", help="File with sample IDs for control group")
    parser.add_argument(
        "--gene-burden-mode",
        choices=["samples", "alleles"],
        default="alleles",
        help="Mode for gene burden calculation: 'samples' or 'alleles'",
    )
    parser.add_argument(
        "--correction-method",
        choices=["fdr", "bonferroni"],
        default="fdr",
        help="Multiple testing correction method for gene burden test",
    )
    parser.add_argument(
        "--html-report",
        action="store_true",
        help="Generate an interactive HTML report with sortable variant tables and summary plots.",
    )
    parser.add_argument(
        "--igv",
        action="store_true",
        help="Enable IGV.js integration for genomic visualization.",
    )
    parser.add_argument(
        "--bam-mapping-file",
        help="Path to a TSV or CSV file mapping sample IDs to BAM files (sample_id,bam_path).",
    )
    parser.add_argument(
        "--igv-reference",
        help="Genome reference identifier for IGV (e.g., 'hg19' or 'hg38'). Required if --igv is enabled unless --igv-fasta is provided.",
    )
    # MODIFIED: Start of local IGV FASTA feature
    parser.add_argument(
        "--igv-fasta",
        help="Path to a local FASTA file for IGV reports. This will be used instead of --igv-reference if both are provided.",
    )
    parser.add_argument(
        "--igv-fasta-index",
        help="Path to FASTA index file (.fai). If not provided, igv-reports will attempt to create or locate it.",
    )
    parser.add_argument(
        "--igv-ideogram",
        help="Path to an ideogram file for chromosome visualization in IGV reports.",
    )
    # MODIFIED: End of local IGV FASTA feature

    parser.add_argument(
        "--remove-sample-substring",
        help="If provided, this substring will be removed from all sample names found in the VCF.",
    )

    parser.add_argument(
        "--no-links",
        action="store_true",
        default=False,
        help="Disable adding link columns to the final output (links are added by default).",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for bcftools and related operations.",
    )

    parser.add_argument(
        "--add-column",
        action="append",
        default=[],
        help="Add a blank column with the specified header name to the final TSV output. "
        "Can be specified multiple times.",
    )
    parser.add_argument(
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

    parser.add_argument(
        "--add-chr",
        action="store_true",
        help="Add 'chr' prefix to chromosome names in BED file. Overrides config setting.",
    )

    parser.add_argument(
        "--split-snpeff-lines",
        nargs="?",
        const="before_filters",
        choices=["before_filters", "after_filters"],
        help=(
            "If set, run the SNPeff annotation splitter step in one of two modes:\n"
            "  'before_filters' (default if no value is given) => Split lines right after variant extraction;\n"
            "       recommended if you're also doing transcript filtering, though it can be slower.\n"
            "  'after_filters' => Only split EFF/ANN lines after the main filter step,\n"
            "       produces a smaller intermediate file and can be faster for big multi-annotation VCFs.\n\n"
            "If omitted, splitting is not performed at all."
        ),
    )

    parser.add_argument(
        "--transcript-list",
        help=(
            "Comma-separated list of transcript IDs to filter for. "
            "For example: NM_007294.4,NM_000059.4"
        ),
    )
    parser.add_argument(
        "--transcript-file",
        help=(
            "Path to a file containing transcript IDs, one per line. "
            "For example: a file where each line is a transcript like NM_007294.4"
        ),
    )

    parser.add_argument(
        "--genotype-filter",
        help=(
            "Specify genotype filter(s), e.g. 'het', 'hom', 'comp_het' or "
            "comma-separated combination (e.g. 'het,comp_het'). Only samples/lines "
            "fulfilling these genotypes are kept (unless per-gene file overrides)."
        ),
    )
    parser.add_argument(
        "--gene-genotype-file",
        help=(
            "Path to a file specifying per-gene genotype rules. It must contain at least two columns: "
            "'GENE' and 'GENOTYPES' (a comma-separated list of genotype filters, e.g. het,comp_het). "
            "If a gene is listed here, it overrides any global --genotype-filter."
        ),
    )

    # Unified argument for extra sample fields:
    parser.add_argument(
        "--append-extra-sample-fields",
        nargs="*",
        default=None,
        help=(
            "Optionally append extra sample fields (e.g. DP, AD) next to each genotype in the final output. "
            "If no arguments are provided, this flag is merely enabled without additional fields. If "
            "arguments (space-separated) are provided, those must also appear in the extracted TSV columns "
            "(see --fields). Example usage: --append-extra-sample-fields DP AD"
        ),
    )
    parser.add_argument(
        "--extra-sample-field-delimiter",
        default=":",
        help=(
            "Delimiter to separate genotype from extra sample fields (and from each other). "
            "Default is ':'."
        ),
    )

    args: argparse.Namespace = parser.parse_args()

    # Configure logging level
    log_level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARN": logging.WARNING,
        "ERROR": logging.ERROR,
    }
    logging.getLogger("variantcentrifuge").setLevel(log_level_map[args.log_level])

    # If a log file is specified, add a file handler
    if args.log_file:
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
    cfg: Dict[str, Any] = load_config(args.config)
    logger.debug(f"Configuration loaded: {cfg}")

    reference: Optional[str] = args.reference or cfg.get("reference")
    filters: Optional[str] = args.filters or cfg.get("filters")
    fields: Optional[str] = args.fields or cfg.get("fields_to_extract")

    # Handle presets
    preset_dict = cfg.get("presets", {})
    if args.preset:
        chosen_presets = []
        for p in args.preset:
            p_str = p.strip()
            if p_str in preset_dict:
                chosen_presets.append(f"({preset_dict[p_str]})")
            else:
                logger.error(f"Preset '{p_str}' not found in config.")
                sys.exit(1)
        if chosen_presets:
            combined_preset_filter = " & ".join(chosen_presets)
            if filters:
                filters = f"({combined_preset_filter}) & ({filters})"
            else:
                filters = combined_preset_filter

    # Validate mandatory parameters
    validate_mandatory_parameters(reference, filters, fields)
    validate_vcf_file(args.vcf_file, logger)

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

    cfg["perform_gene_burden"] = args.perform_gene_burden
    cfg["no_stats"] = args.no_stats
    cfg["gene_burden_mode"] = args.gene_burden_mode
    cfg["correction_method"] = args.correction_method

    # Handle add_chr configuration
    if args.add_chr:
        cfg["add_chr"] = True
    # If flag is not provided, use the value from config file (already loaded)

    # IGV parameters
    cfg["igv_enabled"] = args.igv
    cfg["bam_mapping_file"] = args.bam_mapping_file

    # MODIFIED: Start of local IGV FASTA feature
    # Store local FASTA-related parameters
    cfg["igv_fasta"] = args.igv_fasta
    cfg["igv_fasta_index"] = args.igv_fasta_index
    cfg["igv_ideogram"] = args.igv_ideogram

    # Handle precedence and validation
    if args.igv:
        if not args.bam_mapping_file:
            logger.error("For IGV integration, --bam-mapping-file must be provided.")
            sys.exit(1)

        if not args.igv_fasta and not args.igv_reference:
            logger.error(
                "For IGV integration, either --igv-reference or --igv-fasta must be provided."
            )
            sys.exit(1)

        if args.igv_fasta and args.igv_reference:
            logger.warning(
                "Both --igv-fasta and --igv-reference provided. Local FASTA file will take precedence."
            )

        cfg["igv_reference"] = args.igv_reference
    # MODIFIED: End of local IGV FASTA feature

    # Update reference/filters/fields
    cfg["reference"] = reference
    cfg["filters"] = filters
    cfg["fields_to_extract"] = fields

    # Remove sample substring (optional)
    cfg["remove_sample_substring"] = args.remove_sample_substring or None

    # Toggle link columns
    cfg["no_links"] = args.no_links

    # Threads
    cfg["threads"] = args.threads

    # Gene list annotation
    cfg["annotate_gene_list_files"] = args.annotate_gene_list

    # Transcript list/file
    cfg["transcript_list"] = args.transcript_list
    cfg["transcript_file"] = args.transcript_file

    # Genotype filter arguments
    cfg["genotype_filter"] = args.genotype_filter
    cfg["gene_genotype_file"] = args.gene_genotype_file

    # Unified extra sample fields logic:
    if args.append_extra_sample_fields is not None:
        # Flag used => enable genotype+field appending
        cfg["append_extra_sample_fields"] = True
        # Possibly empty list if used with no args
        cfg["extra_sample_fields"] = args.append_extra_sample_fields
    else:
        # Flag not used at all
        cfg["append_extra_sample_fields"] = False
        cfg["extra_sample_fields"] = []

    cfg["extra_sample_field_delimiter"] = args.extra_sample_field_delimiter

    # Manage the new snpEff splitting mode
    cfg["snpeff_splitting_mode"] = (
        args.split_snpeff_lines
    )  # None, 'before_filters', or 'after_filters'

    run_pipeline(args, cfg, start_time)


if __name__ == "__main__":
    main()
