"""Command-line interface for VariantCentrifuge."""

import argparse
import datetime
import logging
import sys
from pathlib import Path
from typing import Any, Dict, Optional

from .config import load_config
from .pipeline import run_refactored_pipeline
from .validators import validate_mandatory_parameters, validate_phenotype_file, validate_vcf_file
from .version import __version__

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
        choices=["DEBUG", "INFO", "WARN", "ERROR"],
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

    # Gene Selection
    gene_group = parser.add_argument_group("Gene Selection")
    gene_group.add_argument("-g", "--gene-name", help="Gene or list of genes of interest")
    gene_group.add_argument("-G", "--gene-file", help="File containing gene names, one per line")
    gene_group.add_argument(
        "--transcript-list",
        help=(
            "Comma-separated list of transcript IDs to filter for. "
            "For example: NM_007294.4,NM_000059.4"
        ),
    )
    gene_group.add_argument(
        "--transcript-file",
        help=(
            "Path to a file containing transcript IDs, one per line. "
            "For example: a file where each line is a transcript like NM_007294.4"
        ),
    )

    # Filtering & Annotation
    filter_group = parser.add_argument_group("Filtering & Annotation")
    filter_group.add_argument(
        "-r",
        "--reference",
        help="Reference database for snpEff",
        default=None,
    )
    filter_group.add_argument("-f", "--filters", help="Filters to apply in SnpSift filter")
    filter_group.add_argument(
        "--bcftools-prefilter",
        type=str,
        default=None,
        help="Optional bcftools expression to pre-filter the VCF file for performance. "
        "This is applied during variant extraction to reduce data early. "
        "Example: 'FILTER=\"PASS\" && INFO/AC<10'",
    )
    filter_group.add_argument(
        "--preset",
        action="append",
        help=(
            "Apply predefined filtering presets defined in the config file. "
            "Specify multiple times for multiple presets. They are combined with AND. "
            "If custom filters are also given, they are combined with these presets using AND."
        ),
    )
    filter_group.add_argument(
        "--late-filtering",
        action="store_true",
        help="Apply filters after scoring and annotation (allows filtering on computed columns like scores). "
        "By default, filters are applied early on VCF data.",
    )
    filter_group.add_argument(
        "--final-filter",
        type=str,
        default=None,
        help="An expression to filter the final results table. Uses pandas query() syntax. "
        "This is applied after all annotations and scores have been calculated. "
        "Example: 'inheritance_score > 0.5 and IMPACT == \"HIGH\"'",
    )
    filter_group.add_argument(
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

    # Field Extraction & Formatting
    format_group = parser.add_argument_group("Field Extraction & Formatting")
    format_group.add_argument("-e", "--fields", help="Fields to extract with SnpSift extractFields")
    format_group.add_argument(
        "--no-replacement", action="store_true", help="Skip genotype replacement step"
    )
    format_group.add_argument(
        "--add-column",
        action="append",
        default=[],
        help="Add a blank column with the specified header name to the final TSV output. "
        "Can be specified multiple times.",
    )
    format_group.add_argument(
        "--no-links",
        action="store_true",
        default=False,
        help="Disable adding link columns to the final output (links are added by default).",
    )
    format_group.add_argument(
        "--add-chr",
        action="store_true",
        help="Add 'chr' prefix to chromosome names in BED file. Overrides config setting.",
    )
    format_group.add_argument(
        "--remove-sample-substring",
        help="If provided, this substring will be removed from all sample names found in the VCF.",
    )
    # Genotype Analysis
    genotype_group = parser.add_argument_group("Genotype Analysis")
    genotype_group.add_argument(
        "--genotype-filter",
        help=(
            "Specify genotype filter(s), e.g. 'het', 'hom', 'comp_het' or "
            "comma-separated combination (e.g. 'het,comp_het'). Only samples/lines "
            "fulfilling these genotypes are kept (unless per-gene file overrides)."
        ),
    )
    genotype_group.add_argument(
        "--gene-genotype-file",
        help=(
            "Path to a file specifying per-gene genotype rules. It must contain at least two columns: "
            "'GENE' and 'GENOTYPES' (a comma-separated list of genotype filters, e.g. het,comp_het). "
            "If a gene is listed here, it overrides any global --genotype-filter."
        ),
    )
    genotype_group.add_argument(
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
    genotype_group.add_argument(
        "--extra-sample-field-delimiter",
        default=":",
        help=(
            "Delimiter to separate genotype from extra sample fields (and from each other). "
            "Default is ':'."
        ),
    )

    # Phenotype & Sample Groups
    phenotype_group = parser.add_argument_group("Phenotype & Sample Groups")
    phenotype_group.add_argument("--phenotype-file", help="Path to phenotype file (.csv or .tsv)")
    phenotype_group.add_argument(
        "--phenotype-sample-column",
        help="Name of the column containing sample IDs in phenotype file",
    )
    phenotype_group.add_argument(
        "--phenotype-value-column",
        help="Name of the column containing phenotype values in phenotype file",
    )
    phenotype_group.add_argument(
        "--case-phenotypes", help="Comma-separated HPO terms defining case group"
    )
    phenotype_group.add_argument(
        "--control-phenotypes", help="Comma-separated HPO terms defining control group"
    )
    phenotype_group.add_argument(
        "--case-phenotypes-file", help="File with HPO terms for case group"
    )
    phenotype_group.add_argument(
        "--control-phenotypes-file", help="File with HPO terms for control group"
    )
    phenotype_group.add_argument(
        "--case-samples", help="Comma-separated sample IDs defining the case group"
    )
    phenotype_group.add_argument(
        "--control-samples",
        help="Comma-separated sample IDs defining the control group",
    )
    phenotype_group.add_argument("--case-samples-file", help="File with sample IDs for case group")
    phenotype_group.add_argument(
        "--control-samples-file", help="File with sample IDs for control group"
    )
    # Statistical Analysis
    stats_group = parser.add_argument_group("Statistical Analysis")
    stats_group.add_argument(
        "--perform-gene-burden",
        action="store_true",
        help="Perform gene burden analysis",
    )
    stats_group.add_argument(
        "--gene-burden-mode",
        choices=["samples", "alleles"],
        default="alleles",
        help="Mode for gene burden calculation: 'samples' or 'alleles'",
    )
    stats_group.add_argument(
        "--correction-method",
        choices=["fdr", "bonferroni"],
        default="fdr",
        help="Multiple testing correction method for gene burden test",
    )
    stats_group.add_argument(
        "--no-stats",
        action="store_true",
        default=False,
        help="Skip the statistics computation step.",
    )
    stats_group.add_argument("--stats-output-file", help="File to write analysis statistics")
    stats_group.add_argument(
        "--stats-config",
        help="Path to custom statistics configuration JSON file. If not provided, uses default statistics.",
    )
    # Inheritance Analysis
    inheritance_group = parser.add_argument_group("Inheritance Analysis")
    inheritance_group.add_argument(
        "--ped",
        help="Path to the PED file defining family structure for inheritance analysis.",
    )
    inheritance_group.add_argument(
        "--inheritance-mode",
        choices=["simple", "columns", "full"],
        default=None,
        help="Define the output format for inheritance analysis. "
        "'simple': Show pattern only (default if analysis is triggered). "
        "'columns': Unpack key details into separate columns. "
        "'full': Show the complete JSON object in a single column. "
        "Note: Analysis is triggered by this flag or --ped. Without --ped, treats all samples as affected singletons.",
    )
    inheritance_group.add_argument(
        "--no-vectorized-comp-het",
        action="store_true",
        help="Disable vectorized compound heterozygous analysis (use original implementation).",
    )

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
        type=int,
        default=1,
        help="Number of parallel worker processes to use for variant extraction and filtering. "
        "Set to 1 to disable parallelization (default: 1).",
    )
    performance_group.add_argument(
        "--chunk-size",
        type=int,
        default=10000,
        help="Number of variants per chunk for memory-efficient processing (default: 10000). "
        "Automatically adjusts to keep genes together.",
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
        default="2G",
        help="Memory limit for external sort command (default: 2G). " "Examples: 500M, 4G, 8G",
    )
    performance_group.add_argument(
        "--sort-parallel",
        type=int,
        default=4,
        help="Number of parallel threads for sorting (default: 4)",
    )
    performance_group.add_argument(
        "--genotype-replacement-method",
        choices=["auto", "sequential", "vectorized", "chunked-vectorized", "parallel"],
        default="auto",
        help="Method for genotype replacement processing. "
        "auto: automatically select based on data characteristics and memory availability; "
        "sequential: line-by-line streaming (memory efficient); "
        "vectorized: pandas-based vectorized operations (faster for many samples); "
        "chunked-vectorized: memory-safe vectorized processing for large files; "
        "parallel: multi-threaded chunked processing (for very large files). "
        "Default: auto",
    )
    performance_group.add_argument(
        "--vectorized-chunk-size",
        type=int,
        default=5000,
        help="Number of rows per chunk for vectorized genotype replacement (default: 5000). "
        "Only used when processing method is vectorized and file is large.",
    )
    performance_group.add_argument(
        "--vectorized-chunk-threshold-mb",
        type=int,
        default=200,
        help="File size threshold in MB for using chunked vectorized processing (default: 200). "
        "Files larger than this will use chunked processing to manage memory usage.",
    )
    performance_group.add_argument(
        "--disable-parallel-genotype-replacement",
        action="store_true",
        help="Disable parallel genotype replacement even when multiple threads are available. "
        "Forces sequential or vectorized processing.",
    )
    performance_group.add_argument(
        "--max-memory-gb",
        type=float,
        help="Maximum memory in GB to use for genotype processing (auto-detected if not specified). "
        "Used to intelligently select processing method and chunk sizes.",
    )
    performance_group.add_argument(
        "--force-genotype-method",
        choices=[
            "sequential",
            "vectorized",
            "chunked-vectorized",
            "parallel-chunked-vectorized",
            "parallel",
        ],
        help="Force specific genotype replacement method, bypassing auto-selection. "
        "Useful for debugging or when you know the optimal method for your data.",
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
        "--gzip-intermediates",
        action="store_true",
        default=True,
        help="Compress intermediate TSV files with gzip to save disk space. "
        "Uses fast compression (level 1) to optimize I/O performance while reducing disk usage. "
        "Use --no-gzip-intermediates to disable.",
    )
    misc_group.add_argument(
        "--no-gzip-intermediates",
        dest="gzip_intermediates",
        action="store_false",
        help="Disable compression of intermediate TSV files.",
    )
    misc_group.add_argument(
        "--genotype-replacement-chunk-size",
        type=int,
        default=10000,
        help="Chunk size for parallel genotype replacement (default: 10000 variants per chunk). "
        "Larger chunks use more memory but may be more efficient for very large datasets.",
    )
    misc_group.add_argument(
        "--disable-parallel-genotype-replacement",
        action="store_true",
        default=False,
        help="Force sequential genotype replacement even for large datasets. "
        "May be useful for debugging or when memory is very limited.",
    )
    misc_group.add_argument(
        "--chunks",
        type=int,
        default=None,
        help="Split file processing into this many chunks. Alias for --chunk-size for backward compatibility.",
    )
    misc_group.add_argument(
        "--max-variants-per-gene",
        type=int,
        default=10000,
        help="Maximum number of variants allowed per gene before forcing chunk split "
        "(default: 10000). Lower values prevent memory exhaustion with large genes in big cohorts.",
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
        choices=["DEBUG", "INFO", "WARN", "ERROR"],
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

    # Gene Selection
    gene_group = parser.add_argument_group("Gene Selection")
    gene_group.add_argument("-g", "--gene-name", help="Gene or list of genes of interest")
    gene_group.add_argument("-G", "--gene-file", help="File containing gene names, one per line")
    gene_group.add_argument(
        "--transcript-list",
        help=(
            "Comma-separated list of transcript IDs to filter for. "
            "For example: NM_007294.4,NM_000059.4"
        ),
    )
    gene_group.add_argument(
        "--transcript-file",
        help=(
            "Path to a file containing transcript IDs, one per line. "
            "For example: a file where each line is a transcript like NM_007294.4"
        ),
    )

    # Filtering & Annotation
    filter_group = parser.add_argument_group("Filtering & Annotation")
    filter_group.add_argument(
        "-r",
        "--reference",
        help="Reference database for snpEff",
        default=None,
    )
    filter_group.add_argument("-f", "--filters", help="Filters to apply in SnpSift filter")
    filter_group.add_argument(
        "--bcftools-prefilter",
        type=str,
        default=None,
        help="Optional bcftools expression to pre-filter the VCF file for performance. "
        "This is applied during variant extraction to reduce data early. "
        "Example: 'FILTER=\"PASS\" && INFO/AC<10'",
    )
    filter_group.add_argument(
        "--preset",
        action="append",
        help=(
            "Apply predefined filtering presets defined in the config file. "
            "Specify multiple times for multiple presets. They are combined with AND. "
            "If custom filters are also given, they are combined with these presets using AND."
        ),
    )
    filter_group.add_argument(
        "--late-filtering",
        action="store_true",
        help="Apply filters after scoring and annotation (allows filtering on computed columns like scores). "
        "By default, filters are applied early on VCF data.",
    )
    filter_group.add_argument(
        "--final-filter",
        type=str,
        default=None,
        help="An expression to filter the final results table. Uses pandas query() syntax. "
        "This is applied after all annotations and scores have been calculated. "
        "Example: 'inheritance_score > 0.5 and IMPACT == \"HIGH\"'",
    )
    filter_group.add_argument(
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

    # Field Extraction & Formatting
    format_group = parser.add_argument_group("Field Extraction & Formatting")
    format_group.add_argument("-e", "--fields", help="Fields to extract with SnpSift extractFields")
    format_group.add_argument(
        "--no-replacement", action="store_true", help="Skip genotype replacement step"
    )
    format_group.add_argument(
        "--add-column",
        action="append",
        default=[],
        help="Add a blank column with the specified header name to the final TSV output. "
        "Can be specified multiple times.",
    )
    format_group.add_argument(
        "--no-links",
        action="store_true",
        default=False,
        help="Disable adding link columns to the final output (links are added by default).",
    )
    format_group.add_argument(
        "--add-chr",
        action="store_true",
        help="Add 'chr' prefix to chromosome names in BED file. Overrides config setting.",
    )
    format_group.add_argument(
        "--remove-sample-substring",
        help="If provided, this substring will be removed from all sample names found in the VCF.",
    )
    # Genotype Analysis
    genotype_group = parser.add_argument_group("Genotype Analysis")
    genotype_group.add_argument(
        "--genotype-filter",
        help=(
            "Specify genotype filter(s), e.g. 'het', 'hom', 'comp_het' or "
            "comma-separated combination (e.g. 'het,comp_het'). Only samples/lines "
            "fulfilling these genotypes are kept (unless per-gene file overrides)."
        ),
    )
    genotype_group.add_argument(
        "--gene-genotype-file",
        help=(
            "Path to a file specifying per-gene genotype rules. It must contain at least two columns: "
            "'GENE' and 'GENOTYPES' (a comma-separated list of genotype filters, e.g. het,comp_het). "
            "If a gene is listed here, it overrides any global --genotype-filter."
        ),
    )
    genotype_group.add_argument(
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
    genotype_group.add_argument(
        "--extra-sample-field-delimiter",
        default=":",
        help=(
            "Delimiter to separate genotype from extra sample fields (and from each other). "
            "Default is ':'."
        ),
    )

    # Phenotype & Sample Groups
    phenotype_group = parser.add_argument_group("Phenotype & Sample Groups")
    phenotype_group.add_argument("--phenotype-file", help="Path to phenotype file (.csv or .tsv)")
    phenotype_group.add_argument(
        "--phenotype-sample-column",
        help="Name of the column containing sample IDs in phenotype file",
    )
    phenotype_group.add_argument(
        "--phenotype-value-column",
        help="Name of the column containing phenotype values in phenotype file",
    )
    phenotype_group.add_argument(
        "--case-phenotypes", help="Comma-separated HPO terms defining case group"
    )
    phenotype_group.add_argument(
        "--control-phenotypes", help="Comma-separated HPO terms defining control group"
    )
    phenotype_group.add_argument(
        "--case-phenotypes-file", help="File with HPO terms for case group"
    )
    phenotype_group.add_argument(
        "--control-phenotypes-file", help="File with HPO terms for control group"
    )
    phenotype_group.add_argument(
        "--case-samples", help="Comma-separated sample IDs defining the case group"
    )
    phenotype_group.add_argument(
        "--control-samples",
        help="Comma-separated sample IDs defining the control group",
    )
    phenotype_group.add_argument("--case-samples-file", help="File with sample IDs for case group")
    phenotype_group.add_argument(
        "--control-samples-file", help="File with sample IDs for control group"
    )
    # Statistical Analysis
    stats_group = parser.add_argument_group("Statistical Analysis")
    stats_group.add_argument(
        "--perform-gene-burden",
        action="store_true",
        help="Perform gene burden analysis",
    )
    stats_group.add_argument(
        "--gene-burden-mode",
        choices=["samples", "alleles"],
        default="alleles",
        help="Mode for gene burden calculation: 'samples' or 'alleles'",
    )
    stats_group.add_argument(
        "--correction-method",
        choices=["fdr", "bonferroni"],
        default="fdr",
        help="Multiple testing correction method for gene burden test",
    )
    stats_group.add_argument(
        "--no-stats",
        action="store_true",
        default=False,
        help="Skip the statistics computation step.",
    )
    stats_group.add_argument("--stats-output-file", help="File to write analysis statistics")
    stats_group.add_argument(
        "--stats-config",
        help="Path to custom statistics configuration JSON file. If not provided, uses default statistics.",
    )
    # Inheritance Analysis
    inheritance_group = parser.add_argument_group("Inheritance Analysis")
    inheritance_group.add_argument(
        "--ped",
        help="Path to the PED file defining family structure for inheritance analysis.",
    )
    inheritance_group.add_argument(
        "--inheritance-mode",
        choices=["simple", "columns", "full"],
        default=None,
        help="Define the output format for inheritance analysis. "
        "'simple': Show pattern only (default if analysis is triggered). "
        "'columns': Unpack key details into separate columns. "
        "'full': Show the complete JSON object in a single column. "
        "Note: Analysis is triggered by this flag or --ped. Without --ped, treats all samples as affected singletons.",
    )
    inheritance_group.add_argument(
        "--no-vectorized-comp-het",
        action="store_true",
        help="Disable vectorized compound heterozygous analysis (use original implementation).",
    )

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
        type=int,
        default=1,
        help="Number of parallel worker processes to use for variant extraction and filtering. "
        "Set to 1 to disable parallelization (default: 1).",
    )
    performance_group.add_argument(
        "--chunk-size",
        type=int,
        default=10000,
        help="Number of variants per chunk for memory-efficient processing (default: 10000). "
        "Automatically adjusts to keep genes together.",
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
        default="2G",
        help="Memory limit for external sort command (default: 2G). " "Examples: 500M, 4G, 8G",
    )
    performance_group.add_argument(
        "--sort-parallel",
        type=int,
        default=4,
        help="Number of parallel threads for sorting (default: 4)",
    )
    performance_group.add_argument(
        "--genotype-replacement-method",
        choices=["auto", "sequential", "vectorized", "chunked-vectorized", "parallel"],
        default="auto",
        help="Method for genotype replacement processing. "
        "auto: automatically select based on data characteristics and memory availability; "
        "sequential: line-by-line streaming (memory efficient); "
        "vectorized: pandas-based vectorized operations (faster for many samples); "
        "chunked-vectorized: memory-safe vectorized processing for large files; "
        "parallel: multi-threaded chunked processing (for very large files). "
        "Default: auto",
    )
    performance_group.add_argument(
        "--vectorized-chunk-size",
        type=int,
        default=5000,
        help="Number of rows per chunk for vectorized genotype replacement (default: 5000). "
        "Only used when processing method is vectorized and file is large.",
    )
    performance_group.add_argument(
        "--vectorized-chunk-threshold-mb",
        type=int,
        default=200,
        help="File size threshold in MB for using chunked vectorized processing (default: 200). "
        "Files larger than this will use chunked processing to manage memory usage.",
    )
    performance_group.add_argument(
        "--disable-parallel-genotype-replacement",
        action="store_true",
        help="Disable parallel genotype replacement even when multiple threads are available. "
        "Forces sequential or vectorized processing.",
    )
    performance_group.add_argument(
        "--max-memory-gb",
        type=float,
        help="Maximum memory in GB to use for genotype processing (auto-detected if not specified). "
        "Used to intelligently select processing method and chunk sizes.",
    )
    performance_group.add_argument(
        "--force-genotype-method",
        choices=[
            "sequential",
            "vectorized",
            "chunked-vectorized",
            "parallel-chunked-vectorized",
            "parallel",
        ],
        help="Force specific genotype replacement method, bypassing auto-selection. "
        "Useful for debugging or when you know the optimal method for your data.",
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

    # For --show-checkpoint-status, we don't need other required arguments
    if "--show-checkpoint-status" in sys.argv:
        # Create a minimal parser just for checkpoint status
        status_parser = argparse.ArgumentParser(description="Show checkpoint status")
        status_parser.add_argument("--show-checkpoint-status", action="store_true")
        status_parser.add_argument("--output-dir", default="output", help="Output directory")
        status_parser.add_argument(
            "--log-level", choices=["DEBUG", "INFO", "WARN", "ERROR"], default="INFO"
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
        "ERROR": logging.ERROR,
    }
    logging.getLogger("variantcentrifuge").setLevel(log_level_map[args.log_level])

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
    cfg["igv_ideogram"] = args.igv_ideogram

    # MODIFIED: IGV filename shortening parameters
    cfg["igv_max_allele_len_filename"] = args.igv_max_allele_len_filename
    cfg["igv_hash_len_filename"] = args.igv_hash_len_filename
    cfg["igv_max_variant_part_filename"] = args.igv_max_variant_part_filename

    # MODIFIED: Start of IGV flanking feature
    # Use CLI argument if provided, otherwise use config value, with a default of 50
    cfg["igv_flanking"] = args.igv_flanking or cfg.get("igv_flanking", 50)
    # MODIFIED: End of IGV flanking feature

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

    # Output file - important to override default config value
    if args.output_file is not None:
        cfg["output_file"] = args.output_file

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

    # Scoring configuration
    cfg["scoring_config_path"] = args.scoring_config_path

    # Phenotype configuration - ensure these are always set
    cfg["phenotype_file"] = args.phenotype_file
    cfg["phenotype_sample_column"] = args.phenotype_sample_column
    cfg["phenotype_value_column"] = args.phenotype_value_column

    # Unified annotation configuration
    cfg["annotate_bed_files"] = args.annotate_bed
    cfg["annotate_json_genes"] = args.annotate_json_genes
    cfg["json_gene_mapping"] = args.json_gene_mapping
    cfg["json_genes_as_columns"] = args.json_genes_as_columns
    # Note: args.annotate_gene_list is already handled above as "annotate_gene_list_files"
    # We'll map it to the new unified system
    cfg["annotate_gene_lists"] = args.annotate_gene_list

    # Validate unified annotation arguments
    if args.annotate_json_genes and not args.json_gene_mapping:
        logger.error("--annotate-json-genes requires --json-gene-mapping to be provided.")
        sys.exit(1)
    if args.json_gene_mapping and not args.annotate_json_genes:
        logger.error("--json-gene-mapping requires --annotate-json-genes to be provided.")
        sys.exit(1)

    # Inheritance analysis configuration
    cfg["ped_file"] = args.ped
    # Enable inheritance analysis if inheritance mode is explicitly specified OR if PED file is provided
    # This allows single-sample analysis without PED
    cfg["calculate_inheritance"] = args.inheritance_mode is not None or args.ped is not None
    # If inheritance is calculated but no mode specified, default to "simple"
    cfg["inheritance_mode"] = args.inheritance_mode or (
        "simple" if cfg["calculate_inheritance"] else None
    )
    cfg["use_vectorized_comp_het"] = (
        not args.no_vectorized_comp_het
    )  # Default to True unless disabled

    # Late filtering configuration
    cfg["late_filtering"] = args.late_filtering

    # Bcftools pre-filtering configuration
    cfg["bcftools_prefilter"] = args.bcftools_prefilter

    # Final filter configuration
    cfg["final_filter"] = args.final_filter

    # Chunked processing configuration
    cfg["chunk_size"] = args.chunk_size
    cfg["no_chunked_processing"] = args.no_chunked_processing
    cfg["force_chunked_processing"] = args.force_chunked_processing
    cfg["sort_memory_limit"] = args.sort_memory_limit
    cfg["sort_parallel"] = args.sort_parallel

    # Statistics configuration
    cfg["stats_config"] = args.stats_config

    # Checkpoint configuration
    cfg["enable_checkpoint"] = args.enable_checkpoint
    cfg["resume"] = args.resume
    cfg["checkpoint_checksum"] = args.checkpoint_checksum
    cfg["resume_from"] = args.resume_from
    cfg["list_stages"] = args.list_stages
    cfg["list_checkpoints"] = args.list_checkpoints
    cfg["interactive_resume"] = args.interactive_resume
    cfg["pipeline_version"] = __version__  # Add pipeline version for compatibility checking

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

        run_refactored_pipeline(refactored_args)
        return 0
    except SystemExit as e:
        return e.code if e.code is not None else 1
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        return 1


if __name__ == "__main__":
    main()
