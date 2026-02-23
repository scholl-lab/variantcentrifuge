# File: variantcentrifuge/association/tests/skat_r.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/skat_r.py
"""
RSKATTest — AssociationTest wrapper for the R SKAT backend.

Wraps RSKATBackend to integrate R-based SKAT into the AssociationEngine
framework. Registered in the engine registry as ``"skat"``.

SKAT has no effect size, standard error, or confidence interval — it only
produces a p-value (and optionally rho from SKAT-O). The engine's None-effect
guard handles the all-None ``effect_column_names()`` return without creating
malformed column names like ``skat_None``.

Null model lifecycle
--------------------
The SKAT null model is fit once on the full cohort (all samples, phenotype +
covariates) and then reused for all genes. RSKATTest caches the null model
as ``self._null_model`` after the first gene call. This is safe because the
null model does not depend on the variant data — it depends only on the
phenotype, covariates, and trait type.

Thread safety
-------------
RSKATBackend is not thread-safe (rpy2 restriction). The stage that uses
RSKATTest must declare ``parallel_safe=False``.

Lifecycle hooks
---------------
The engine calls ``prepare(n_genes)`` before the gene loop and ``finalize()``
after. RSKATTest uses these for:
  - prepare(): large panel warning, log interval computation, start message
  - finalize(): total timing summary, backend cleanup (R gc())
"""

from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING, Any

import numpy as np

from variantcentrifuge.association.backends.r_backend import RSKATBackend
from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
from variantcentrifuge.association.tests._utils import (
    parse_weights_beta as _parse_weights_beta_shared,
)

if TYPE_CHECKING:
    from variantcentrifuge.association.backends.base import NullModelResult

logger = logging.getLogger("variantcentrifuge")


class RSKATTest(AssociationTest):
    """
    # DEPRECATED
    SKAT association test using the R SKAT package via rpy2.

    Registered in the engine as ``"skat"``. Uses RSKATBackend for all R
    interaction; this class handles the AssociationTest protocol (lifecycle,
    null model caching, result packing).

    Parameters
    ----------
    None — instantiated by AssociationEngine.from_names() via the registry.

    Raises
    ------
    ImportError
        From check_dependencies() if rpy2 is not importable or the R SKAT
        package is not installed.

    Examples
    --------
    This test is typically used through the engine, not instantiated directly:

    >>> engine = AssociationEngine.from_names(["skat"], config)
    >>> result_df = engine.run_all(gene_burden_data)
    """

    parallel_safe: bool = False  # rpy2 restriction: main thread only

    def __init__(self) -> None:
        import warnings

        warnings.warn(
            "RSKATTest (R SKAT backend) is deprecated and will be removed in v0.17.0. "
            "Use --skat-backend python (default).",
            DeprecationWarning,
            stacklevel=2,
        )
        self._backend: RSKATBackend | None = None
        self._null_model: NullModelResult | None = None

        # Lifecycle tracking (set by prepare())
        self._total_genes: int = 0
        self._log_interval: int = 50
        self._genes_processed: int = 0
        self._start_time: float = 0.0

    @property
    def name(self) -> str:
        """Short identifier for this test used in registry lookup and column prefixes."""
        return "skat"

    def check_dependencies(self) -> None:
        """
        Verify that rpy2 is importable and the R SKAT package is installed.

        Calls RSKATBackend.detect_environment() eagerly so that users get a
        clear error (with R_HOME and install instructions) before any data
        processing begins.

        Raises
        ------
        ImportError
            If rpy2 is not installed, R is not found (R_HOME misconfigured),
            or the SKAT package is not installed in R.
        """
        from variantcentrifuge.association.backends import get_skat_backend

        # RSKATTest is always backed by the R backend — the "r" name is
        # hard-coded here because this class IS the R SKAT test. The engine
        # registry maps "skat" -> RSKATTest; an auto-selected Python backend
        # would be a separate test class (Phase 21).
        backend = get_skat_backend("r")
        backend.detect_environment()
        backend.log_environment()
        self._backend = backend  # type: ignore[assignment]

    def effect_column_names(self) -> dict[str, str | None]:
        """
        Column name suffixes for SKAT output.

        SKAT produces only a p-value (plus SKAT-O rho and other metadata in
        ``extra``). There is no effect size, standard error, or confidence
        interval. All four slots are None; the engine's None-guard skips
        column creation for them.

        Returns
        -------
        dict with all values None.
        """
        return {
            "effect": None,
            "se": None,
            "ci_lower": None,
            "ci_upper": None,
        }

    def prepare(self, gene_count: int) -> None:
        """
        Called by the engine before the gene loop with total gene count.

        Sets up progress logging parameters and emits pre-run messages:
        - WARNING if gene_count > LARGE_PANEL_THRESHOLD (>2000 genes)
        - INFO log announcing start: "R SKAT: beginning analysis of N genes"

        Parameters
        ----------
        gene_count : int
            Total number of genes to be processed.
        """
        self._total_genes = gene_count
        self._genes_processed = 0
        self._log_interval = max(10, min(50, gene_count // 10)) if gene_count > 0 else 50
        self._start_time = time.time()

        if gene_count > RSKATBackend.LARGE_PANEL_THRESHOLD:
            logger.warning(
                f"R SKAT: large gene panel ({gene_count} genes > "
                f"{RSKATBackend.LARGE_PANEL_THRESHOLD} threshold). "
                "Memory usage will be monitored. Consider splitting the panel "
                "if R heap warnings appear."
            )

        logger.info(f"R SKAT: beginning analysis of {gene_count} genes")

    def finalize(self) -> None:
        """
        Called by the engine after the gene loop completes.

        Logs aggregate timing summary and triggers R backend cleanup (R gc()).
        """
        elapsed = time.time() - self._start_time
        logger.info(f"R SKAT complete: {self._genes_processed} genes in {elapsed:.1f}s")
        if self._backend is not None:
            self._backend.cleanup()

    def run(
        self,
        gene: str,
        contingency_data: dict[str, Any],
        config: AssociationConfig,
    ) -> TestResult:
        """
        Run SKAT for a single gene.

        Fits the null model on the first call (lazy), then reuses it for all
        subsequent genes. The genotype matrix, phenotype, and covariate matrix
        are extracted from ``contingency_data`` (same keys as logistic_burden).

        Parameters
        ----------
        gene : str
            Gene symbol being tested.
        contingency_data : dict
            Must contain:
              - ``genotype_matrix`` : np.ndarray (n_samples, n_variants)
              - ``variant_mafs``     : np.ndarray (n_variants,)
              - ``phenotype_vector`` : np.ndarray (n_samples,)
              Optional:
              - ``covariate_matrix`` : np.ndarray (n_samples, k) or None
        config : AssociationConfig
            Runtime configuration (skat_method, trait_type, variant_weights).

        Returns
        -------
        TestResult
            p_value=None when test is skipped (no variants, no genotype matrix).
            effect_size=None, se=None, ci_lower=None, ci_upper=None (SKAT has
            no effect size). extra contains SKAT-specific metadata.

        Notes
        -----
        The backend must be initialised (check_dependencies() called) before
        run() is invoked. The engine calls check_dependencies() at construction
        time (from_names()), so this is guaranteed in normal usage.
        """
        n_cases = int(contingency_data.get("proband_count", 0))
        n_controls = int(contingency_data.get("control_count", 0))
        n_variants = int(contingency_data.get("n_qualifying_variants", 0))

        # Require genotype matrix — absent when only Fisher is requested
        if "genotype_matrix" not in contingency_data:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"skat_warnings": ["NO_GENOTYPE_MATRIX"]},
            )

        geno: np.ndarray = contingency_data["genotype_matrix"]
        phenotype: np.ndarray | None = contingency_data.get("phenotype_vector")
        covariate_matrix: np.ndarray | None = contingency_data.get("covariate_matrix")

        if phenotype is None or geno.shape[0] == 0 or geno.shape[1] == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"skat_warnings": ["NO_PHENOTYPE_OR_EMPTY_MATRIX"]},
            )

        # Backend must have been initialised by check_dependencies()
        if self._backend is None:
            raise RuntimeError(
                "RSKATTest.run() called without check_dependencies(). "
                "Use AssociationEngine.from_names() which calls check_dependencies() "
                "at construction time."
            )

        # Lazy null model: fit once, reuse for all genes
        if self._null_model is None:
            self._null_model = self._backend.fit_null_model(
                phenotype=phenotype,
                covariates=covariate_matrix,
                trait_type=config.trait_type,
            )

        # Parse weight parameters from config (e.g. "beta:1,25" -> (1.0, 25.0))
        weights_beta = _parse_weights_beta(config.variant_weights)

        result = self._backend.test_gene(
            gene=gene,
            genotype_matrix=geno,
            null_model=self._null_model,
            method=config.skat_method,
            weights_beta=weights_beta,
        )

        # Increment local counter (backend also tracks, but we need per-test tracking)
        self._genes_processed += 1

        # Periodic GC + heap monitoring (every GC_INTERVAL genes)
        if self._genes_processed % RSKATBackend.GC_INTERVAL == 0:
            self._backend._run_r_gc()
            heap_mb = self._backend._check_r_heap()
            if heap_mb is not None and heap_mb > RSKATBackend.R_HEAP_WARNING_GB * 1024.0:
                logger.warning(
                    f"R heap usage {heap_mb / 1024.0:.1f} GB exceeds "
                    f"{RSKATBackend.R_HEAP_WARNING_GB} GB threshold"
                )

        # Progress logging at INFO level every log_interval genes
        if self._total_genes > 0 and self._genes_processed % self._log_interval == 0:
            pct = 100.0 * self._genes_processed / self._total_genes
            logger.info(
                f"SKAT progress: {self._genes_processed}/{self._total_genes} genes ({pct:.0f}%)"
            )

        return TestResult(
            gene=gene,
            test_name=self.name,
            p_value=result.get("p_value"),
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            se=None,
            n_cases=n_cases,
            n_controls=n_controls,
            n_variants=result.get("n_variants", n_variants),
            extra={
                "skat_o_rho": result.get("rho"),
                "skat_warnings": "; ".join(result.get("warnings", [])) or None,
                "skat_method": config.skat_method,
                "skat_n_marker_test": result.get("n_marker_test"),
            },
        )


def _parse_weights_beta(variant_weights: str) -> tuple[float, float]:
    """
    Parse weight specification string to Beta distribution parameters.

    Delegates to the shared implementation in
    ``variantcentrifuge.association.tests._utils.parse_weights_beta``.
    Kept here as a module-level name for backward compatibility with any
    code that imports it directly from this module.

    Parameters
    ----------
    variant_weights : str
        Weight scheme, e.g. ``"beta:1,25"``, ``"uniform"``, or any string.

    Returns
    -------
    tuple of (float, float)
        Beta distribution parameters (a1, a2).

    Examples
    --------
    >>> _parse_weights_beta("beta:1,25")
    (1.0, 25.0)
    >>> _parse_weights_beta("uniform")
    (1.0, 1.0)
    """
    return _parse_weights_beta_shared(variant_weights)
