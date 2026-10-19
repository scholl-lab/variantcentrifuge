"""CLI argument group definitions for VariantCentrifuge."""

from __future__ import annotations

import argparse

from .cli_analysis_groups import add_analysis_argument_groups

__all__ = [
    "add_analysis_argument_groups",
    "add_operational_argument_groups",
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
    report_group.add_argument(
        "--igv-max-allele-len-filename",
        type=int,
        default=10,
        help="Maximum length for REF/ALT alleles in IGV report filenames (default: 10). "
        "Longer alleles will be truncated and hashed.",
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
