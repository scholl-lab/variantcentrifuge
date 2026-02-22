# File: variantcentrifuge/association/tests/skat_python.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/skat_python.py
"""
PurePythonSKATTest — AssociationTest wrapper for the pure Python SKAT backend.

Wraps PythonSKATBackend to integrate pure-Python SKAT into the AssociationEngine
framework. Registered in the engine registry as ``"skat_python"`` (explicit) and
as a swap target for ``"skat"`` when ``--skat-backend python`` is specified.

This class follows the same structural pattern as RSKATTest (skat_r.py) but:
- Uses PythonSKATBackend instead of RSKATBackend
- Is thread-safe (Python SKAT is pure numpy/scipy)
- Has no R heap monitoring (no R process involved)
- Includes p_method and p_converged in TestResult.extra

SKAT has no effect size, standard error, or confidence interval — it only
produces a p-value (and optionally rho from SKAT-O). The engine's None-effect
guard handles the all-None ``effect_column_names()`` return.

Null model lifecycle
--------------------
The null model is fit once on the first gene call (lazy) and cached as
``self._null_model``. This is safe because the null model does not depend
on the variant data — only on the phenotype, covariates, and trait type.

Thread safety
-------------
PythonSKATBackend is thread-safe. The stage using PurePythonSKATTest does
NOT need to declare ``parallel_safe=False`` (unlike RSKATTest).

Lifecycle hooks
---------------
The engine calls ``prepare(n_genes)`` before the gene loop and ``finalize()``
after. PurePythonSKATTest uses these for:
  - prepare(): log interval computation, start message
  - finalize(): total timing summary, backend cleanup (no-op for Python backend)
"""

from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING, Any

import numpy as np

from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
from variantcentrifuge.association.tests._utils import parse_weights_beta
from variantcentrifuge.association.tests.acat import compute_acat_v
from variantcentrifuge.association.weights import beta_maf_weights

if TYPE_CHECKING:
    from variantcentrifuge.association.backends.base import NullModelResult
    from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

logger = logging.getLogger("variantcentrifuge")


class PurePythonSKATTest(AssociationTest):
    """
    SKAT association test using the pure Python backend (numpy/scipy/statsmodels).

    Registered in the engine as ``"skat_python"`` (explicit) and used as the
    ``"skat"`` implementation when ``--skat-backend python`` is specified.

    Parameters
    ----------
    None — instantiated by AssociationEngine.from_names() via the registry.

    Raises
    ------
    ImportError
        From check_dependencies() if required Python libraries are missing
        (should not occur in normal variantcentrifuge installations).

    Examples
    --------
    This test is typically used through the engine, not instantiated directly:

    >>> engine = AssociationEngine.from_names(["skat_python"], config)
    >>> result_df = engine.run_all(gene_burden_data)
    """

    def __init__(self) -> None:
        self._backend: PythonSKATBackend | None = None
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
        Verify that the Python SKAT backend dependencies are available.

        Calls PythonSKATBackend.detect_environment() eagerly so that users
        get clear diagnostics before any data processing begins.

        The Python backend never raises (numpy/scipy/statsmodels are required
        variantcentrifuge dependencies). This method still calls detect_environment()
        to initialize internal state and log environment info.
        """
        from variantcentrifuge.association.backends import get_skat_backend

        backend = get_skat_backend("python")
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

        Sets up progress logging parameters and emits a start message at INFO level.

        Parameters
        ----------
        gene_count : int
            Total number of genes to be processed.
        """
        self._total_genes = gene_count
        self._genes_processed = 0
        self._log_interval = max(10, min(50, gene_count // 10)) if gene_count > 0 else 50
        self._start_time = time.time()

        logger.info(f"Python SKAT: beginning analysis of {gene_count} genes")

    def finalize(self) -> None:
        """
        Called by the engine after the gene loop completes.

        Logs aggregate timing summary and calls backend cleanup (no-op for Python backend).
        """
        elapsed = time.time() - self._start_time
        logger.info(f"Python SKAT complete: {self._genes_processed} genes in {elapsed:.1f}s")
        if self._backend is not None:
            self._backend.cleanup()

    def run(
        self,
        gene: str,
        contingency_data: dict[str, Any],
        config: AssociationConfig,
    ) -> TestResult:
        """
        Run SKAT for a single gene using the pure Python backend.

        Fits the null model on the first call (lazy), then reuses it for all
        subsequent genes. The genotype matrix, phenotype, and covariate matrix
        are extracted from ``contingency_data``.

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
            no effect size). extra contains SKAT-specific metadata including
            p_method, p_converged, and skip_reason.

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
                extra={"skat_warnings": "NO_GENOTYPE_MATRIX", "acat_v_p": None},
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
                extra={"skat_warnings": "NO_PHENOTYPE_OR_EMPTY_MATRIX", "acat_v_p": None},
            )

        # Backend must have been initialised by check_dependencies()
        if self._backend is None:
            raise RuntimeError(
                "PurePythonSKATTest.run() called without check_dependencies(). "
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
        weights_beta = parse_weights_beta(config.variant_weights)

        result = self._backend.test_gene(
            gene=gene,
            genotype_matrix=geno,
            null_model=self._null_model,
            method=config.skat_method,
            weights_beta=weights_beta,
        )

        # Compute ACAT-V per-variant score test (Phase 25)
        # Uses same Beta(MAF) weights as SKAT for consistency.
        mafs = geno.mean(axis=0) / 2.0
        a1, a2 = weights_beta
        acat_v_weights = beta_maf_weights(mafs, a=a1, b=a2)
        acat_v_p = compute_acat_v(
            geno=geno,
            residuals=self._null_model.extra["residuals"],
            trait_type=config.trait_type,
            sigma2=self._null_model.extra["sigma2"],
            mu_hat=self._null_model.extra.get("mu_hat"),  # None for quantitative
            weights=acat_v_weights,
        )

        # Increment local counter
        self._genes_processed += 1

        # Progress logging at INFO level every log_interval genes
        if self._total_genes > 0 and self._genes_processed % self._log_interval == 0:
            pct = 100.0 * self._genes_processed / self._total_genes
            logger.info(
                f"Python SKAT progress: "
                f"{self._genes_processed}/{self._total_genes} genes ({pct:.0f}%)"
            )

        # Build extra dict with Python backend-specific metadata
        warnings_raw = result.get("warnings", [])
        warnings_str: str | None = "; ".join(warnings_raw) if warnings_raw else None
        skip_reason = result.get("skip_reason")

        extra: dict[str, Any] = {
            "skat_o_rho": result.get("rho"),
            "skat_warnings": warnings_str,
            "skat_method": config.skat_method,
            "skat_n_marker_test": result.get("n_marker_test"),
            "skat_p_method": result.get("p_method"),
            "skat_p_converged": result.get("p_converged"),
            "acat_v_p": acat_v_p,
        }
        if skip_reason is not None:
            extra["skat_skip_reason"] = skip_reason

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
            extra=extra,
        )
