# File: variantcentrifuge/association/base.py
# Location: variantcentrifuge/variantcentrifuge/association/base.py
"""
Core abstractions for the association testing framework.

Defines the AssociationTest abstract base class, TestResult dataclass, and
AssociationConfig dataclass that all association tests in the v0.15.0
framework must implement.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger("variantcentrifuge")


@dataclass
class TestResult:
    """
    Result from a single association test on a single gene.

    Fields
    ------
    gene : str
        Gene symbol tested.
    test_name : str
        Short test identifier (e.g. "fisher", "burden", "skat").
    p_value : float | None
        Raw (uncorrected) p-value. None when test is skipped (e.g. zero
        qualifying variants) — None is NOT the same as 1.0 (failure vs skip).
    corrected_p_value : float | None
        Multiple-testing-corrected p-value. Populated by AssociationEngine
        after all genes are tested. None until correction is applied.
    effect_size : float | None
        Primary effect size estimate. For Fisher: odds ratio. For regression
        tests (Phase 19+): beta coefficient.
    ci_lower : float | None
        Lower bound of the confidence interval for effect_size.
    ci_upper : float | None
        Upper bound of the confidence interval for effect_size.
    se : float | None
        Standard error of the effect size estimate. First-class field for
        regression tests (burden, SKAT); None for non-regression tests (Fisher).
    n_cases : int
        Total number of case samples in the analysis.
    n_controls : int
        Total number of control samples in the analysis.
    n_variants : int
        Number of qualifying variants for this gene (after filtering).
    extra : dict
        Test-specific ancillary data (e.g. contingency table, convergence
        flags). Not written to output unless a formatter explicitly accesses it.
    """

    gene: str
    test_name: str
    p_value: float | None
    corrected_p_value: float | None
    effect_size: float | None
    ci_lower: float | None
    ci_upper: float | None
    se: float | None
    n_cases: int
    n_controls: int
    n_variants: int
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class AssociationConfig:
    """
    Configuration for the association analysis framework.

    All fields have defaults that mirror the equivalent keys read from the
    cfg dict in gene_burden.py (lines 373-377), so existing workflows can
    transition without changing semantics.

    Fields
    ------
    correction_method : str
        Multiple-testing correction method. "fdr" (Benjamini-Hochberg) or
        "bonferroni". Default: "fdr".
    gene_burden_mode : str
        Collapsing strategy for carrier counting. "samples" = unique carrier
        samples (CMC/CAST); "alleles" = max allele dosage summed across
        samples (preserves diploid constraint). Default: "samples".
    confidence_interval_method : str
        Method for odds ratio CI computation. Currently "normal_approx" which
        tries score -> normal -> logit methods in sequence. Default:
        "normal_approx".
    confidence_interval_alpha : float
        Significance level for CIs. 0.05 gives 95% CIs. Default: 0.05.
    continuity_correction : float
        Value added to each cell when a zero is present, to stabilise CI
        computation. Default: 0.5 (Haldane-Anscombe correction).
    """

    correction_method: str = "fdr"
    gene_burden_mode: str = "samples"
    confidence_interval_method: str = "normal_approx"
    confidence_interval_alpha: float = 0.05
    continuity_correction: float = 0.5

    # Phase 19: Covariate system + burden test fields
    # All optional; backward compatible with Phase 18 workflows.
    covariate_file: str | None = None
    """Path to tab/CSV covariate file (first column = sample ID). None = no covariates."""

    covariate_columns: list[str] | None = None
    """Subset of covariate columns to use. None = all columns."""

    categorical_covariates: list[str] | None = None
    """Column names to one-hot encode. None = auto-detect (non-numeric with <=5 levels)."""

    trait_type: str = "binary"
    """Phenotype scale: "binary" (logistic, Firth fallback) or "quantitative" (linear OLS)."""

    variant_weights: str = "beta:1,25"
    """Weight scheme: "beta:a,b" (Beta MAF weights, SKAT convention) or "uniform"."""

    variant_weight_params: dict | None = None
    """Extra parameters for weight schemes (e.g. {'cadd_cap': 40.0})."""

    missing_site_threshold: float = 0.10
    """Variants with >threshold fraction missing site-wide are excluded before imputation."""

    missing_sample_threshold: float = 0.80
    """Samples with >threshold fraction missing across kept variants are excluded."""

    firth_max_iter: int = 25
    """Maximum Newton-Raphson iterations for Firth penalized logistic regression fallback."""

    # Phase 20: R SKAT backend fields
    # All optional; backward compatible with Phase 18-19 workflows.
    skat_backend: str = "python"
    """SKAT computation backend: "python" (default), "r" (deprecated, R via rpy2), or "auto"."""

    skat_method: str = "SKAT"
    """SKAT variant to run: "SKAT" (default), "Burden" (burden-only), or "SKATO" (omnibus)."""

    # Phase 22: ACAT-O + diagnostics fields
    min_cases: int = 200
    """Cohort-level warning threshold: warn if n_cases < this value."""

    max_case_control_ratio: float = 20.0
    """Cohort-level warning threshold: warn if n_controls/n_cases > this value."""

    min_case_carriers: int = 10
    """Per-gene warning threshold: flag genes with case_carriers < this value."""

    diagnostics_output: str | None = None
    """Path to diagnostics output directory. None = no diagnostics output."""

    # Phase 23: PCA integration fields
    pca_file: str | None = None
    """Path to pre-computed PCA file (PLINK .eigenvec, AKT output, or generic TSV)."""

    pca_tool: str | None = None
    """PCA computation tool: 'akt' to invoke AKT as subprocess. None = pre-computed file only."""

    pca_components: int = 10
    """Number of principal components to use. Default: 10. Warn if >20."""

    # Phase 23: COAST allelic series fields
    coast_weights: list[float] | None = None
    """Category weights for COAST allelic series (default: [1.0, 2.0, 3.0] for BMV, DMV, PTV)."""

    # Phase 24: Pure Python COAST backend
    coast_backend: str = "python"
    """COAST computation backend: "python" (default), "r" (deprecated, R via rpy2), or "auto"."""

    # Phase 31: Configurable COAST classification model
    coast_classification: str | None = None
    """Absolute path to a scoring/coast_classification/<model>/ directory.
    None = use built-in SIFT/PolyPhen hardcoded logic (backward-compatible default).
    Set by cli.py after resolving --coast-classification model name to a path."""

    # Phase 33: Gene-level FDR weighting
    gene_prior_weights: str | None = None
    """Path to gene-to-weight TSV file for weighted BH FDR correction.
    None = standard (unweighted) BH/Bonferroni (backward-compatible default)."""

    gene_prior_weight_column: str = "weight"
    """Column name in the weight file containing weight values. Default: 'weight'."""

    # Phase 27: Gene-level parallelization
    association_workers: int = 1
    """Number of parallel worker processes for gene-level association analysis.
    Default: 1 (sequential). Set > 1 for parallel execution via ProcessPoolExecutor.
    Set -1 for os.cpu_count(). Only effective when all registered tests have parallel_safe=True."""


class AssociationTest(ABC):
    """
    Abstract base class for all association tests.

    Subclasses implement a specific statistical test (Fisher, burden regression,
    SKAT, etc.) and are registered in AssociationEngine's test registry. Each
    subclass is responsible for a single gene at a time; the engine handles
    iteration and correction.

    Methods
    -------
    name : str (property)
        Short, lowercase identifier used for test registry lookup and output
        column prefixes (e.g. "fisher", "burden", "skat").
    run(gene, contingency_data, config) -> TestResult
        Execute the test for one gene and return a TestResult.
    check_dependencies() -> None
        Raise ImportError if required optional libraries are missing. The
        default implementation is a no-op; subclasses override as needed.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Short lowercase identifier for this test (e.g. 'fisher')."""
        ...

    @abstractmethod
    def run(
        self, gene: str, contingency_data: dict[str, Any], config: AssociationConfig
    ) -> TestResult:
        """
        Run the association test for a single gene.

        Parameters
        ----------
        gene : str
            Gene symbol being tested.
        contingency_data : dict
            Gene-level aggregated data. Keys available from gene_burden
            aggregation:
              - proband_count       : int, total case samples
              - control_count       : int, total control samples
              - proband_carrier_count : int, case carrier samples
              - control_carrier_count : int, control carrier samples
              - proband_allele_count  : int, total alt alleles in cases
              - control_allele_count  : int, total alt alleles in controls
              - n_qualifying_variants : int, variants passing filters
        config : AssociationConfig
            Runtime configuration (correction method, mode, CI params).

        Returns
        -------
        TestResult
            Result with p_value=None when test is skipped (e.g. zero variants).
        """
        ...

    def effect_column_names(self) -> dict[str, str | None]:
        """
        Column name suffixes for this test's effect size output.

        Returns a mapping of semantic role to column suffix. The engine uses
        these to build output column names as ``{test_name}_{suffix}``.

        Default returns OR-based naming (appropriate for Fisher's exact test).
        Regression tests (burden, SKAT) override to return beta/SE naming.
        """
        return {
            "effect": "or",
            "se": None,
            "ci_lower": "or_ci_lower",
            "ci_upper": "or_ci_upper",
        }

    def check_dependencies(self) -> None:  # noqa: B027
        """
        Verify that required optional dependencies are available.

        Raises
        ------
        ImportError
            If a required library is not installed. Called eagerly at engine
            construction so users get a clear error before processing begins.
        """

    def prepare(self, gene_count: int) -> None:  # noqa: B027
        """
        Called by the engine before the per-gene loop.

        Default is a no-op. Subclasses override to set up progress logging,
        emit large-panel warnings, initialize timers, etc.

        Parameters
        ----------
        gene_count : int
            Total number of genes that will be processed.
        """

    def finalize(self) -> None:  # noqa: B027
        """
        Called by the engine after the per-gene loop completes.

        Default is a no-op. Subclasses override to log aggregate timing,
        release resources, or perform post-run cleanup.
        """
