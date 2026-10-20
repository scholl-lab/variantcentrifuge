"""Statistical Analysis (association) CLI argument group."""

from __future__ import annotations

import argparse


def _fraction(value: str) -> float:
    """Argparse type: a float in [0, 1] (same rule as the JSON config check)."""
    try:
        f = float(value)
    except ValueError:
        f = float("nan")
    if not 0.0 <= f <= 1.0:  # NaN fails too
        raise argparse.ArgumentTypeError(f"must be a number between 0 and 1, got {value!r}")
    return f


def add_association_argument_group(parser: argparse.ArgumentParser) -> None:
    """Register the Statistical Analysis argument group on *parser*."""
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
        "--score-artifact",
        type=str,
        default=None,
        help=(
            "Also write this cohort's per-gene score/covariance artifact (.npz + .meta.json) "
            "for cross-cohort meta-analysis with `variantcentrifuge meta`. Binary traits only."
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
    # Genotype missing-data policy for matrix-based tests and --score-artifact (#200).
    # None = not set on the CLI, so the JSON "association" section can still provide it.
    stats_group.add_argument(
        "--missing-site-threshold",
        type=_fraction,
        default=None,
        help=(
            "Drop variants with a missing-genotype fraction above this value before "
            "imputation (matrix-based tests and --score-artifact). Default: 0.10."
        ),
    )
    stats_group.add_argument(
        "--missing-sample-threshold",
        type=_fraction,
        default=None,
        help=(
            "Exclude samples with a missing-genotype fraction above this value across the "
            "kept variants of a gene. Default: 0.80."
        ),
    )
    stats_group.add_argument(
        "--missing-as-ref",
        action="store_true",
        default=None,
        help=(
            "Treat missing genotypes (./.) as reference (0/0) in every association test "
            "(Fisher, matrix-based tests, --score-artifact) instead of excluding or imputing "
            "them. Use for merged single-sample VCFs where ./. means 'no variant "
            "called'. Where ./. is a genuine no-call, carriers are miscounted as reference, which "
            "can bias results either way if missingness differs between cases and controls."
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
