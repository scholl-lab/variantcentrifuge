"""Domain-specific analysis CLI argument group definitions for VariantCentrifuge."""

from __future__ import annotations

import argparse

from .cli_association_groups import add_association_argument_group


def add_analysis_argument_groups(parser: argparse.ArgumentParser) -> None:
    """Register gene, filtering, formatting, genotype, phenotype, stats, and inheritance groups."""
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
            "For GRCh38, use MANE Select transcript lists matching the VCF namespace "
            "(RefSeq NM_* or Ensembl ENST*). For hg19/GRCh37 RefSeq annotations, use "
            "a GRCh37 RefSeq Select transcript list because MANE is defined on GRCh38."
        ),
    )
    gene_group.add_argument(
        "--transcript-fallback",
        choices=["none", "keep-gene"],
        default=None,
        help=(
            "What to do with a gene none of whose annotated transcripts is in the transcript "
            "list. 'none' (default): drop it (reported as a WARNING and, with "
            "--diagnostics-output, in transcript_filter_dropped_genes.tsv). 'keep-gene': keep "
            "all of that gene's annotations."
        ),
    )
    gene_group.add_argument(
        "--regions-bed",
        help=(
            "BED file to restrict variant extraction to specified regions (e.g., capture kit BED). "
            "The gene BED is intersected with this file before variant extraction; with "
            "--gene-name all, the regions BED is used directly."
        ),
        default=None,
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
        "--field-profile",
        help=(
            "Field profile for annotation-version-specific field names "
            "(e.g., 'dbnsfp4', 'dbnsfp5'). Controls which gnomAD field names "
            "are used in filter presets and field extraction. "
            "Default: value of 'default_field_profile' in config (usually 'dbnsfp4'). "
            "Use --list-field-profiles to see available profiles."
        ),
        default=None,
    )
    filter_group.add_argument(
        "--list-field-profiles",
        action="store_true",
        default=False,
        help="List available field profiles and exit.",
    )
    # Tumor-Normal Filtering
    tumor_group = parser.add_argument_group("Tumor-Normal Filtering")
    tumor_group.add_argument(
        "--tumor-sample-index",
        type=int,
        default=None,
        help="0-based index of the tumor sample in the VCF. "
        "Used to expand {tumor_idx} in somatic/LOH presets. Default: 1",
    )
    tumor_group.add_argument(
        "--normal-sample-index",
        type=int,
        default=None,
        help="0-based index of the normal sample in the VCF. "
        "Used to expand {normal_idx} in somatic/LOH presets. Default: 0",
    )
    tumor_group.add_argument(
        "--tumor-dp-min",
        type=int,
        default=None,
        help="Minimum read depth for tumor sample (default: 20).",
    )
    tumor_group.add_argument(
        "--normal-dp-min",
        type=int,
        default=None,
        help="Minimum read depth for normal sample (default: 20).",
    )
    tumor_group.add_argument(
        "--tumor-af-min",
        type=float,
        default=None,
        help="Minimum allele frequency for tumor sample (default: 0.05).",
    )
    tumor_group.add_argument(
        "--normal-af-max",
        type=float,
        default=None,
        help="Maximum allele frequency for normal sample (default: 0.03).",
    )

    filter_group.add_argument(
        "--show-vcf-annotations",
        action="store_true",
        default=False,
        help="Show available INFO and FORMAT fields from the VCF header and exit. "
        "Useful for discovering field names before configuring extraction. "
        "Requires -v/--vcf-file.",
    )
    filter_group.add_argument(
        "--annotation-filter",
        type=str,
        default=None,
        help="Filter --show-vcf-annotations output to fields matching this substring "
        "(case-insensitive). Example: --annotation-filter gnomAD",
    )
    filter_group.add_argument(
        "--annotation-format",
        choices=["table", "json"],
        default="table",
        help="Output format for --show-vcf-annotations (default: table).",
    )
    ann_source_group = filter_group.add_mutually_exclusive_group()
    ann_source_group.add_argument(
        "--info-only",
        action="store_true",
        default=False,
        help="With --show-vcf-annotations, show only INFO fields.",
    )
    ann_source_group.add_argument(
        "--format-only",
        action="store_true",
        default=False,
        help="With --show-vcf-annotations, show only FORMAT fields.",
    )

    filter_group.add_argument(
        "--late-filtering",
        action="store_true",
        help="Apply filters after scoring and annotation (allows filtering on computed columns). "
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
            "  'before_filters' (default if no value is given) => Split right after extraction;\n"
            "       required for row-level ANN/LOF/NMD consequence qualification when combining\n"
            "       transcript filtering with gene burden or association analysis.\n"
            "  'after_filters' => Only split EFF/ANN lines after the main filter step,\n"
            "       produces a smaller intermediate file and can be faster for big VCFs,\n"
            "       but SnpSift consequence predicates are evaluated at record level.\n\n"
            "If omitted, splitting is not performed at all."
        ),
    )

    filter_group.add_argument(
        "--filter-exception-file",
        help=(
            "5-col TSV of somatic hotspots to rescue to FILTER=PASS during VCF filtering "
            "(does not apply to --bcftools-prefilter or --late-filtering)."
        ),
    )

    # Field Extraction & Formatting
    format_group = parser.add_argument_group("Field Extraction & Formatting")
    format_group.add_argument("-e", "--fields", help="Fields to extract with SnpSift extractFields")
    format_group.add_argument(
        "--fields-file",
        metavar="FIELDS_FILE",
        help="Path to a file containing fields to extract with SnpSift extractFields.",
    )
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
        type=str,
        default=None,
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
            "Path to a file specifying per-gene genotype rules with at least two columns: "
            "'GENE' and 'GENOTYPES' (comma-separated genotype filters, e.g. het,comp_het). "
            "If a gene is listed here, it overrides any global --genotype-filter."
        ),
    )
    genotype_group.add_argument(
        "--append-extra-sample-fields",
        nargs="*",
        default=None,
        help=(
            "Optionally append extra sample fields (e.g. DP, AD) next to each genotype in GT. "
            "If no arguments are provided, this flag is enabled without additional fields. "
            "If arguments (space-separated) are provided, they must appear in extracted columns "
            "(see --fields). Example usage: --append-extra-sample-fields DP AD"
        ),
    )
    genotype_group.add_argument(
        "--extra-sample-field-delimiter",
        default=None,
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
    add_association_argument_group(parser)
    # Inheritance Analysis
    inheritance_group = parser.add_argument_group("Inheritance Analysis")
    inheritance_group.add_argument(
        "--ped",
        help="Path to the PED file defining family structure for inheritance analysis.",
    )
    inheritance_group.add_argument(
        "--ped-file",
        dest="ped",
        help=argparse.SUPPRESS,
        default=None,
    )
    inheritance_group.add_argument(
        "--inheritance-mode",
        choices=["simple", "columns", "full"],
        default=None,
        help="Define the output format for inheritance analysis. "
        "'simple': Show pattern only (default if analysis is triggered). "
        "'columns': Unpack key details into separate columns. "
        "'full': Show the complete JSON object in a single column. "
        "Without --ped, treats all samples as affected singletons.",
    )

    # ClinVar PM5 Annotation
    pm5_group = parser.add_argument_group("ClinVar PM5 Annotation")
    pm5_group.add_argument(
        "--clinvar-pm5-lookup",
        type=str,
        default=None,
        help="Path to a pre-built PM5 lookup table (from variantcentrifuge-build-pm5). "
        "When provided, adds PM5, PM5_evidence_count, and PM5_known_variants columns.",
    )
