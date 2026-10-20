"""Domain-specific analysis CLI argument group definitions for VariantCentrifuge."""

from __future__ import annotations

import argparse


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
    stats_group = parser.add_argument_group("Statistical Analysis")
    stats_group.add_argument(
        "--perform-gene-burden",
        action="store_true",
        help="Perform gene burden analysis",
    )
    stats_group.add_argument(
        "--gene-burden-mode",
        choices=["samples", "alleles"],
        default=None,
        help="Mode for gene burden calculation: 'samples' or 'alleles'",
    )
    stats_group.add_argument(
        "--correction-method",
        choices=["fdr", "bonferroni"],
        default=None,
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
        help="Path to custom statistics configuration JSON file (uses default if omitted).",
    )
    # Association Analysis
    stats_group.add_argument(
        "--perform-association",
        action="store_true",
        help="Perform association analysis using the modular association framework",
    )
    stats_group.add_argument(
        "--association-tests",
        type=str,
        default=None,
        help=(
            "Comma-separated list of association tests to run (default: fisher). "
            "Available: fisher, logistic_burden, linear_burden, skat, skat_python, coast"
        ),
    )
    stats_group.add_argument(
        "--skat-backend",
        choices=["auto", "r", "python"],
        default="python",
        help="SKAT computation backend: python (default), r (deprecated), or auto",
    )
    stats_group.add_argument(
        "--coast-backend",
        choices=["auto", "python"],
        default="python",
        help="COAST computation backend: python (default) or auto",
    )
    stats_group.add_argument(
        "--coast-classification",
        type=str,
        default="sift_polyphen",
        help=(
            "COAST classification model for BMV/DMV/PTV assignment. "
            "Built-in models: 'sift_polyphen' (default), 'cadd'. "
            "Custom models: provide directory name under scoring/coast_classification/."
        ),
    )
    stats_group.add_argument(
        "--covariate-file",
        type=str,
        default=None,
        help=(
            "Path to covariate file (TSV/CSV). First column = sample ID, "
            "remaining columns = covariates. Required when using logistic_burden or linear_burden."
        ),
    )
    stats_group.add_argument(
        "--covariates",
        type=str,
        default=None,
        help=(
            "Comma-separated covariate column names to use from the covariate file. "
            "Default: all columns."
        ),
    )
    stats_group.add_argument(
        "--categorical-covariates",
        type=str,
        default=None,
        help=(
            "Comma-separated column names to force-treat as categorical (one-hot encoded). "
            "Default: auto-detect (non-numeric columns with <=5 unique values)."
        ),
    )
    stats_group.add_argument(
        "--trait-type",
        choices=["binary", "quantitative"],
        default="binary",
        help=(
            "Phenotype type for burden tests. 'binary' uses logistic regression with "
            "Firth fallback; 'quantitative' uses OLS. Default: binary."
        ),
    )
    stats_group.add_argument(
        "--diagnostics-output",
        type=str,
        default=None,
        help=(
            "Path to directory for diagnostics output (lambda_GC, QQ data, summary). "
            "Created automatically if it does not exist."
        ),
    )
    stats_group.add_argument(
        "--variant-weights",
        type=str,
        default="beta:1,25",
        help=(
            "Variant weighting scheme for burden tests. "
            "'beta:a,b' = Beta(MAF;a,b) density (default: 'beta:1,25', SKAT convention). "
            "'uniform' = equal weights. "
            "'cadd' = Beta(MAF) x CADD_phred/40. "
            "'revel' = Beta(MAF) x REVEL_score. "
            "'combined' = Beta(MAF) x functional score (CADD preferred)."
        ),
    )
    stats_group.add_argument(
        "--variant-weight-column",
        type=str,
        default=None,
        help=(
            "Variant-level numeric score column to use when --variant-weights score_column "
            "is set. Ignored when --variant-weights uses inline column:<name>."
        ),
    )
    stats_group.add_argument(
        "--variant-weight-params",
        type=str,
        default=None,
        help=(
            "JSON string of extra weight parameters, e.g. "
            "'{\"cadd_cap\": 30}'. Overrides default weight parameters."
        ),
    )
    # Phase 23: COAST allelic series weights
    stats_group.add_argument(
        "--coast-weights",
        type=str,
        default=None,
        help=(
            "Comma-separated category weights for COAST allelic series test "
            "(BMV, DMV, PTV). Default: '1,2,3'. Example: '1,4,9' for quadratic scaling."
        ),
    )
    # Phase 27: Gene-level parallelization
    stats_group.add_argument(
        "--association-workers",
        type=int,
        default=0,
        help=(
            "Number of parallel worker processes for association analysis. "
            "Default: 0 (auto — uses --threads value). Set to 1 for sequential, "
            "-1 for all cores. "
            "Only effective with Python backends (R backend is not parallel-safe)."
        ),
    )
    # Phase 28: SKAT method and diagnostic thresholds
    stats_group.add_argument(
        "--skat-method",
        choices=["SKAT", "Burden", "SKATO"],
        default="SKAT",
        help=(
            "SKAT method variant: 'SKAT' (default), 'Burden' (burden-only), "
            "or 'SKATO' (SKAT-O omnibus). Requires SKAT test in --association-tests."
        ),
    )
    stats_group.add_argument(
        "--min-cases",
        type=int,
        default=200,
        help=(
            "Minimum number of cases for diagnostics warning. "
            "Default: 200. Warns when n_cases < this value."
        ),
    )
    stats_group.add_argument(
        "--max-case-control-ratio",
        type=float,
        default=20.0,
        help=(
            "Maximum case:control ratio for diagnostics warning. "
            "Default: 20.0. Warns when n_controls/n_cases > this value."
        ),
    )
    stats_group.add_argument(
        "--min-case-carriers",
        type=int,
        default=10,
        help=(
            "Minimum case carrier count per gene for diagnostics warning. "
            "Default: 10. Flags genes where case_carriers < this value."
        ),
    )
    # Phase 32: Unified PCA argument (replaces --pca-file and --pca-tool)
    stats_group.add_argument(
        "--pca",
        dest="pca",
        default=None,
        help=(
            "PCA eigenvectors: path to pre-computed file (PLINK .eigenvec, AKT output, "
            "or generic TSV), or 'akt' to compute via AKT as a pipeline stage."
        ),
    )
    stats_group.add_argument(
        "--pca-components",
        type=int,
        default=10,
        help="Number of principal components to use as covariates (default: 10).",
    )
    stats_group.add_argument(
        "--pca-sites",
        default=None,
        help=(
            "VCF/BCF with sites to restrict AKT PCA computation (passed as -R to akt). "
            "Use the wes.grch37.vcf.gz or wgs.grch37.vcf.gz files shipped with AKT. "
            "Without this, --pca akt uses --force to run on all sites."
        ),
    )
    # Phase 33: Gene-level FDR weighting
    stats_group.add_argument(
        "--gene-prior-weights",
        default=None,
        help=(
            "TSV file with per-gene weights for weighted Benjamini-Hochberg FDR correction. "
            "First column: gene name, weight column (default: 'weight'): positive float. "
            "Genes not in the file receive weight=1.0. Activates weighted BH automatically."
        ),
    )
    stats_group.add_argument(
        "--gene-prior-weight-column",
        default="weight",
        help="Column name in the gene weight file containing weight values (default: 'weight')",
    )
    # Deprecated aliases (backward compatibility — hidden from help)
    stats_group.add_argument("--pca-file", dest="pca", help=argparse.SUPPRESS, default=None)
    stats_group.add_argument("--pca-tool", dest="pca", help=argparse.SUPPRESS, default=None)
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
