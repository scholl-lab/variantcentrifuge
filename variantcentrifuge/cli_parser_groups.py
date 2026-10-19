"""CLI argument group definitions for VariantCentrifuge."""

from __future__ import annotations

import argparse

from .cli_analysis_groups import add_analysis_argument_groups

__all__ = [
    "add_analysis_argument_groups",
    "add_operational_argument_groups",
    "add_performance_and_advanced_argument_groups",
    "get_curated_examples_epilog",
    "populate_advanced_group",
    "populate_misc_group",
    "populate_performance_group",
]


def add_operational_argument_groups(parser: argparse.ArgumentParser) -> None:
    """Register scoring, reporting, performance, checkpoint, privacy, and misc argument groups."""
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
        help="Path to a BED file for annotation. Adds 'Region=<name>' to the Custom_Annotation "
        "column for overlapping variants. Can be used multiple times.",
    )
    scoring_group.add_argument(
        "--annotate-gene-list",
        action="append",
        default=[],
        metavar="GENE_LIST_FILE_PATH",
        help="Path to a gene list file (one gene per line, case-insensitive matching). "
        "A new column will be added to the output TSV for each file, "
        "named after a sanitized version of the file's basename (e.g., 'my_gene_set'), "
        "indicating if the variant's GENE is in this list ('yes'/'no'). "
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
        help="JSON string to map fields in JSON gene files. "
        'e.g., \'{"identifier":"gene_symbol","dataFields":["panel","inheritance"]}\'',
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
        help="Genome reference identifier for IGV (e.g., 'hg19' or 'hg38'). "
        "Required if --igv is enabled unless --igv-fasta is provided.",
    )
    report_group.add_argument(
        "--igv-fasta",
        help="Path to a local FASTA file for IGV reports. The index file (.fai) must exist "
        "in the same location with the same name (e.g., reference.fa.fai). "
        "This will be used instead of --igv-reference if both are provided.",
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
        "This provides the most reliable validation of file integrity when resuming checkpoints. "
        "Without this option, interrupted stages will always be re-executed for safety.",
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
        help="Restart pipeline from a specific stage, re-executing that and subsequent stages. "
        "Stages before the specified stage remain completed. "
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


def populate_performance_group(performance_group: argparse._ArgumentGroup) -> None:
    """Populate performance and resource management arguments."""
    performance_group.add_argument(
        "--threads",
        default=None,
        help="Number of parallel worker processes to use for variant extraction and filtering. "
        "Use 'auto' to detect available CPU cores (default: auto). "
        "Set to 1 to disable parallelization.",
    )
    performance_group.add_argument(
        "--sort-memory-limit",
        type=str,
        default=None,
        help="Memory limit for external sort command (default: auto - scales with memory). "
        "Examples: 500M, 4G, 8G, auto. When 'auto', uses ~10%% of available memory.",
    )
    performance_group.add_argument(
        "--sort-parallel",
        type=int,
        default=4,
        help="Number of parallel threads for sorting (default: 4)",
    )
    performance_group.add_argument(
        "--max-memory-gb",
        type=float,
        help="Maximum memory in GB to use for genotype processing (auto-detected if omitted). "
        "Used to intelligently select processing method and chunk sizes.",
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


def populate_advanced_group(advanced_group: argparse._ArgumentGroup) -> None:
    """Populate advanced, developer, and internal tuning arguments."""
    advanced_group.add_argument(
        "--no-chunked-processing",
        action="store_true",
        help="Disable chunked processing even for large files. "
        "May cause memory issues with very large datasets.",
    )
    advanced_group.add_argument(
        "--force-chunked-processing",
        action="store_true",
        help="Force chunked processing even for small files. "
        "Useful for testing or memory-constrained environments.",
    )
    advanced_group.add_argument(
        "--force-inheritance-processing",
        action="store_true",
        help="Force inheritance analysis even if it exceeds safe memory limits. "
        "Use with caution on systems with limited memory.",
    )
    advanced_group.add_argument(
        "--sample-column-creation-method",
        choices=["iterative", "vectorized", "auto"],
        default="auto",
        help="Method for creating sample columns from GT field during inheritance analysis. "
        "iterative: traditional row-by-row regex parsing; "
        "vectorized: pandas vectorized string operations (5-10x faster); "
        "auto: automatically select based on dataset size. Default: auto",
    )
    advanced_group.add_argument(
        "--memory-safety-factor",
        type=float,
        default=0.92,
        help="Fraction of allocated memory to use safely (default: 0.92 = 92%%). "
        "Lower values are more conservative, higher values use more memory.",
    )
    advanced_group.add_argument(
        "--inheritance-memory-fraction",
        type=float,
        default=0.85,
        help="Fraction of safe memory to allocate for inheritance analysis (default: 0.85 = 85%%). "
        "Remainder is reserved for other pipeline stages.",
    )
    advanced_group.add_argument(
        "--no-vectorized-comp-het",
        action="store_true",
        help="Disable vectorized compound heterozygous analysis (use original implementation).",
    )
    advanced_group.add_argument(
        "--igv-max-allele-len-filename",
        type=int,
        default=10,
        help="Maximum length for REF/ALT alleles in IGV report filenames (default: 10). "
        "Longer alleles will be truncated and hashed.",
    )
    advanced_group.add_argument(
        "--igv-hash-len-filename",
        type=int,
        default=6,
        help="Length of hash to append when truncating alleles in filenames (default: 6).",
    )
    advanced_group.add_argument(
        "--igv-max-variant-part-filename",
        type=int,
        default=50,
        help="Maximum length for the variant part of IGV report filenames (default: 50).",
    )


def populate_misc_group(misc_group: argparse._ArgumentGroup) -> None:
    """Populate miscellaneous arguments."""
    misc_group.add_argument(
        "--no-gzip-intermediates",
        dest="gzip_intermediates",
        action="store_false",
        default=True,
        help="Disable compression of intermediate TSV files (compression is enabled by default).",
    )


def add_performance_and_advanced_argument_groups(parser: argparse.ArgumentParser) -> None:
    """Register performance, resource management, and advanced/internal tuning argument groups."""
    performance_group = parser.add_argument_group("Performance & Processing")
    populate_performance_group(performance_group)

    advanced_group = parser.add_argument_group("Advanced & Tuning Options")
    populate_advanced_group(advanced_group)

    misc_group = parser.add_argument_group("Miscellaneous Options")
    populate_misc_group(misc_group)


def get_curated_examples_epilog() -> str:
    """Return curated CLI examples for the help epilog."""
    return """\
Examples:
  # 1. Basic gene-centric variant extraction:
  variantcentrifuge -g BICC1 -v input.vcf.gz -o bicc1_variants.tsv

  # 2. Extract multiple candidate genes with preset filtering and Excel output:
  variantcentrifuge -G candidate_genes.txt -v input.vcf.gz \\
      --preset rare,coding --xlsx -o results.tsv

  # 3. Family inheritance analysis (de novo, recessive, compound het) with PED file:
  variantcentrifuge -G genes.txt -v family.vcf.gz --ped family.ped \\
      --inheritance-mode simple --preset rare,coding

  # 4. Interactive HTML report and summary visualization:
  variantcentrifuge -G genes.txt -v cohort.vcf.gz --preset rare,coding --html-report --xlsx

  # 5. Variant scoring and downstream expression filtering:
  variantcentrifuge -G genes.txt -v input.vcf.gz --scoring-config-path scoring/model \\
      --final-filter 'score > 0.8 and IMPACT == "HIGH"' -o scored.tsv

  # 6. Case-control gene burden and association analysis (SKAT / logistic burden) with covariates:
  variantcentrifuge -G genes.txt -v cohort.vcf.gz --phenotype-file phenotypes.tsv \\
      --phenotype-sample-column sample_id --phenotype-value-column affected \\
      --case-phenotypes affected --control-phenotypes unaffected \\
      --perform-association --association-tests skat,logistic_burden --covariate-file covariates.tsv
"""
