"""
Pure Python SKAT backend using numpy, scipy, and statsmodels.

PythonSKATBackend implements SKATBackend using only Python-native scientific
libraries — no R or rpy2 required. P-values are computed via the three-tier
Davies -> saddlepoint -> Liu chain from davies.py and skat_math.py.

Lifecycle
---------
1. Instantiate PythonSKATBackend()
2. Call detect_environment() — verifies numpy/scipy/statsmodels; never raises
3. Call log_environment() — emits INFO version line
4. Call fit_null_model() — fit once per cohort via statsmodels GLM
5. Call test_gene() for each gene — SKAT, Burden, or SKAT-O
6. Call cleanup() — no-op (no external resources)

Null model
----------
Uses statsmodels GLM with Binomial family (binary) or Gaussian (quantitative).
Extracts resid_response (y - mu_hat) — NOT resid_deviance or resid_pearson,
matching the R SKAT convention.
"""

from __future__ import annotations

import logging
import time
from typing import Any

import numpy as np
import scipy.linalg
import scipy.stats

from variantcentrifuge.association.backends.base import NullModelResult, SKATBackend
from variantcentrifuge.association.backends.davies import (
    _try_load_davies,
    compute_pvalue,
)
from variantcentrifuge.association.backends.skat_math import (
    _GL_A,
    _GL_B,
    _GL_NODES_RAW,
    _GL_W,
    _GL_WEIGHTS_RAW,
    _GL_X,
    _SKATO_RHO_GRID,
    _get_lambda,
    _skato_each_q,
    _skato_get_pvalue,
    _skato_integrate_davies,
    _skato_integrate_liu,
    _skato_optimal_param,
    compute_eigenvalues_filtered,
    compute_skato_q_and_z1,
)

__all__ = [
    "_GL_A",
    "_GL_B",
    "_GL_NODES_RAW",
    "_GL_W",
    "_GL_WEIGHTS_RAW",
    "_GL_X",
    "_SKATO_RHO_GRID",
    "PythonSKATBackend",
    "_get_lambda",
    "_skato_each_q",
    "_skato_get_pvalue",
    "_skato_integrate_davies",
    "_skato_integrate_liu",
    "_skato_optimal_param",
    "compute_eigenvalues_filtered",
    "compute_skato_q_and_z1",
]

logger = logging.getLogger("variantcentrifuge")


class PythonSKATBackend(SKATBackend):
    """SKAT backend using numpy, scipy, and statsmodels (no R required)."""

    def __init__(self) -> None:
        self._numpy_version: str = "<unknown>"
        self._scipy_version: str = "<unknown>"
        self._statsmodels_version: str = "<unknown>"
        self._davies_available: bool = False
        self._start_time: float = time.time()
        self._genes_processed: int = 0

    def detect_environment(self) -> None:
        """
        Verify that numpy, scipy, and statsmodels are available.

        Also probes the Davies C extension. Never raises — Python deps are
        always available (they are declared dependencies of variantcentrifuge).
        If Davies is unavailable, logs INFO and continues with saddlepoint/Liu.
        """
        # Collect numpy version
        try:
            self._numpy_version = np.__version__
        except Exception:
            self._numpy_version = "<unknown>"

        # Collect scipy version
        try:
            import scipy

            self._scipy_version = scipy.__version__
        except Exception:
            self._scipy_version = "<unknown>"

        # Collect statsmodels version
        try:
            import importlib.metadata

            self._statsmodels_version = importlib.metadata.version("statsmodels")
        except Exception:
            try:
                import statsmodels

                self._statsmodels_version = statsmodels.__version__
            except Exception:
                self._statsmodels_version = "<unknown>"

        # Probe Davies C extension (non-fatal)
        self._davies_available = _try_load_davies()
        if not self._davies_available:
            logger.info("Davies C extension unavailable; will use saddlepoint/Liu fallback")

    def log_environment(self) -> None:
        """
        Log versions of numpy, scipy, statsmodels and Davies availability at INFO level.

        Should be called after detect_environment() succeeds.
        """
        davies_status = "available" if self._davies_available else "unavailable (using fallback)"
        logger.info(
            f"Python SKAT backend: numpy={self._numpy_version}, "
            f"scipy={self._scipy_version}, "
            f"statsmodels={self._statsmodels_version}, "
            f"davies_c_ext={davies_status}"
        )

    def fit_null_model(
        self,
        phenotype: np.ndarray,
        covariates: np.ndarray | None,
        trait_type: str,
    ) -> NullModelResult:
        """
        Fit the SKAT null model using statsmodels GLM.

        Uses Binomial family for binary traits, Gaussian for quantitative.
        Extracts resid_response (y - mu_hat), matching R SKAT convention.

        For binary traits, also precomputes the projection half-matrix
        (``proj_half``) needed for the R SKAT-compatible eigenvalue
        computation (see ``_compute_eigenvalues_filtered``).

        Parameters
        ----------
        phenotype : np.ndarray, shape (n_samples,)
            Binary (0/1) or continuous phenotype.
        covariates : np.ndarray or None, shape (n_samples, k)
            Covariate matrix. None = intercept only.
        trait_type : str
            "binary" -> Binomial GLM; "quantitative" -> Gaussian GLM.

        Returns
        -------
        NullModelResult
            Contains fitted model with residuals, sigma2, mu_hat, and
            design_matrix in extra dict.
        """
        import statsmodels.api as sm

        n_obs = len(phenotype)

        # Build design matrix: covariates + intercept
        if covariates is not None and covariates.ndim == 2 and covariates.shape[1] > 0:
            design_matrix = sm.add_constant(covariates, prepend=True, has_constant="add")
        else:
            # Intercept only
            design_matrix = np.ones((n_obs, 1), dtype=np.float64)

        # Select GLM family based on trait type
        family = sm.families.Binomial() if trait_type == "binary" else sm.families.Gaussian()

        # Fit the null model
        glm = sm.GLM(phenotype, design_matrix, family=family)
        fit_result = glm.fit()

        # Extract response residuals: y - mu_hat (NOT deviance or Pearson)
        residuals = np.asarray(fit_result.resid_response, dtype=np.float64)

        # Variance estimate
        sigma2 = (
            float(np.var(residuals, ddof=1))
            if trait_type == "quantitative"
            else 1.0  # Binary trait: sigma2 = 1.0 by SKAT convention
        )

        mu_hat = np.asarray(fit_result.fittedvalues, dtype=np.float64)

        logger.debug(
            f"Null model fit: trait_type={trait_type}, n_samples={n_obs}, sigma2={sigma2:.6f}"
        )

        # Reset tracking for this cohort
        self._genes_processed = 0
        self._start_time = time.time()

        return NullModelResult(
            model=fit_result,
            trait_type=trait_type,
            n_samples=n_obs,
            adjustment=False,
            extra={
                "residuals": residuals,
                "sigma2": sigma2,
                "mu_hat": mu_hat,
                "design_matrix": design_matrix,
            },
        )

    def test_gene(
        self,
        gene: str,
        genotype_matrix: np.ndarray,
        null_model: NullModelResult,
        method: str,
        weights: np.ndarray | None = None,
        weights_beta: tuple[float, float] | None = None,
    ) -> dict[str, Any]:
        """
        Run the SKAT test for a single gene.

        Routes to _test_skat(), _test_burden(), or _test_skato() based on method.

        Parameters
        ----------
        gene : str
            Gene symbol (for logging only).
        genotype_matrix : np.ndarray, shape (n_samples, n_variants)
            Dosage matrix (0/1/2). No NaN.
        null_model : NullModelResult
            Fitted null model from fit_null_model().
        method : str
            "SKAT", "Burden", or "SKATO".
        weights : np.ndarray or None
            Concrete per-variant weights aligned to genotype_matrix columns.
        weights_beta : tuple of (float, float) or None
            Deprecated legacy parameter accepted only for API compatibility.

        Returns
        -------
        dict with keys: p_value, rho, n_variants, n_marker_test, warnings,
                        p_method, p_converged, skip_reason (if applicable)
        """
        self._genes_processed += 1
        n_variants = genotype_matrix.shape[1]
        weights_arr = self._validate_variant_weights(genotype_matrix, weights)

        method_upper = method.upper()
        # Map R SKAT method names to our internal names
        # "optimal.adj" and "optimal" are R's names for SKAT-O
        if method_upper == "SKAT":
            result = self._test_skat(gene, genotype_matrix, null_model, weights_arr)
        elif method_upper == "BURDEN":
            result = self._test_burden(gene, genotype_matrix, null_model, weights_arr)
        elif method_upper in ("SKATO", "OPTIMAL.ADJ", "OPTIMAL"):
            result = self._test_skato(gene, genotype_matrix, null_model, weights_arr)
        else:
            logger.warning(f"Unknown SKAT method '{method}' for gene {gene}; using SKAT")
            result = self._test_skat(gene, genotype_matrix, null_model, weights_arr)

        logger.debug(
            f"Gene {gene}: method={method}, p={result.get('p_value')}, "
            f"n_variants={n_variants}, p_method={result.get('p_method')}"
        )
        return result

    def _validate_variant_weights(
        self, genotype_matrix: np.ndarray, weights: np.ndarray | None
    ) -> np.ndarray:
        if weights is None:
            raise ValueError("PythonSKATBackend requires concrete variant weights")
        weights_arr = np.asarray(weights, dtype=np.float64)
        n_variants = genotype_matrix.shape[1]
        if len(weights_arr) != n_variants:
            raise ValueError(
                f"weights length must match genotype_matrix columns "
                f"(got {len(weights_arr)} and {n_variants})"
            )
        if not np.isfinite(weights_arr).all():
            raise ValueError("weights must be finite")
        return weights_arr

    def _compute_eigenvalues_filtered(
        self,
        geno_weighted: np.ndarray,
        null_model: NullModelResult,
    ) -> np.ndarray:
        """Compute filtered eigenvalues of a SKAT kernel matrix."""
        return compute_eigenvalues_filtered(geno_weighted, null_model)

    def _test_skat(
        self,
        gene: str,
        geno: np.ndarray,
        null_model: NullModelResult,
        weights: np.ndarray,
    ) -> dict[str, Any]:
        """
        SKAT score test for a single gene.

        Parameters
        ----------
        gene : str
            Gene symbol (for logging).
        geno : np.ndarray, shape (n_samples, n_variants)
            Genotype dosage matrix.
        null_model : NullModelResult
            Fitted null model (contains residuals and sigma2 in extra).
        weights : np.ndarray
            Concrete per-variant weights aligned to genotype columns.

        Returns
        -------
        dict with SKAT result fields.
        """
        _n_samples, n_variants = geno.shape
        residuals: np.ndarray = null_model.extra["residuals"]

        # Guard: zero-variant matrix cannot be tested
        if n_variants == 0:
            logger.debug(f"Gene {gene}: n_variants=0; returning p_value=None")
            return {
                "p_value": None,
                "rho": None,
                "n_variants": 0,
                "n_marker_test": 0,
                "warnings": [],
                "p_method": None,
                "p_converged": False,
                "skip_reason": "rank_deficient",
            }

        # Rank check on FULL matrix BEFORE any eigenvalue filtering
        rank = int(np.linalg.matrix_rank(geno))
        if rank < 2:
            logger.debug(f"Gene {gene}: rank={rank} < 2; skipping (rank_deficient)")
            return {
                "p_value": None,
                "rho": None,
                "n_variants": n_variants,
                "n_marker_test": n_variants,
                "warnings": [],
                "p_method": None,
                "p_converged": False,
                "skip_reason": "rank_deficient",
            }

        # Weighted genotype matrix via broadcasting (memory-efficient vs geno @ diag(w))
        geno_weighted = geno * weights[np.newaxis, :]

        # Score statistic: Q = 0.5 * (ZW^T r)^T (ZW^T r)
        score_vec = geno_weighted.T @ residuals  # shape: (n_variants,)
        q_stat = float(score_vec @ score_vec) / 2.0

        # Filtered eigenvalues (projection-adjusted for binary traits)
        lambdas = self._compute_eigenvalues_filtered(geno_weighted, null_model)

        if len(lambdas) == 0:
            logger.debug(f"Gene {gene}: no eigenvalues above threshold; returning p=1.0")
            return {
                "p_value": 1.0,
                "rho": None,
                "n_variants": n_variants,
                "n_marker_test": n_variants,
                "warnings": [],
                "p_method": "liu",
                "p_converged": False,
                "skip_reason": None,
            }

        # Compute p-value via three-tier chain
        p_value, p_method, p_converged = compute_pvalue(q_stat, lambdas)

        # Map NaN/Inf to None (degenerate case)
        p_out: float | None = None if (np.isnan(p_value) or np.isinf(p_value)) else float(p_value)

        return {
            "p_value": p_out,
            "rho": None,
            "n_variants": n_variants,
            "n_marker_test": n_variants,
            "warnings": [],
            "p_method": p_method,
            "p_converged": p_converged,
            "skip_reason": None,
        }

    def _test_burden(
        self,
        gene: str,
        geno: np.ndarray,
        null_model: NullModelResult,
        weights: np.ndarray,
    ) -> dict[str, Any]:
        """
        Burden score test for a single gene.

        Collapses variants to a weighted burden score and performs a score test
        against the null model residuals.

        Parameters
        ----------
        gene : str
            Gene symbol (for logging).
        geno : np.ndarray, shape (n_samples, n_variants)
            Genotype dosage matrix.
        null_model : NullModelResult
            Fitted null model.
        weights : np.ndarray
            Concrete per-variant weights aligned to genotype columns.

        Returns
        -------
        dict with burden result fields.
        """
        _n_samples, n_variants = geno.shape
        residuals: np.ndarray = null_model.extra["residuals"]
        sigma2: float = null_model.extra["sigma2"]

        # Collapse to burden score per sample: b = geno @ weights
        burden = geno @ weights  # shape: (n_samples,)

        # Score test: score = r^T b; Var(score) = sigma2 * ||b||^2
        score = float(residuals @ burden)
        variance = sigma2 * float(burden @ burden)

        if variance <= 0.0:
            return {
                "p_value": None,
                "rho": None,
                "n_variants": n_variants,
                "n_marker_test": n_variants,
                "warnings": ["burden_zero_variance"],
                "p_method": "analytical",
                "p_converged": True,
                "skip_reason": "zero_variance",
            }

        z_stat = score / np.sqrt(variance)
        p_value = float(2.0 * scipy.stats.norm.sf(abs(z_stat)))

        return {
            "p_value": p_value,
            "rho": None,
            "n_variants": n_variants,
            "n_marker_test": n_variants,
            "warnings": [],
            "p_method": "analytical",
            "p_converged": True,
            "skip_reason": None,
        }

    def _test_skato(
        self,
        gene: str,
        geno: np.ndarray,
        null_model: NullModelResult,
        weights: np.ndarray,
    ) -> dict[str, Any]:
        """
        SKAT-O (optimal unified SKAT) for a single gene.

        Implements the full Lee et al. (2012) optimal unified test, matching
        R SKAT's ``SKAT_Optimal_Logistic`` / ``SKAT_Optimal_Linear`` exactly.

        Algorithm:
        1. Compute Q_rho for each rho in the grid
        2. Build projection-adjusted genotype matrix Z1
        3. For each rho, compute eigenvalues via Cholesky mixing of Z1/sqrt(2)
        4. Compute per-rho p-values via Davies -> Liu chain
        5. Compute omnibus p-value by integrating over chi-squared(1)
           (conditions on the shared burden component)

        References
        ----------
        Lee, S., Wu, M.C., Lin, X. (2012). Optimal tests for rare variant
        effects in sequencing association studies. Biostatistics 13:762-775.

        Parameters
        ----------
        gene : str
            Gene symbol (for logging).
        geno : np.ndarray, shape (n_samples, n_variants)
            Genotype dosage matrix.
        null_model : NullModelResult
            Fitted null model.
        weights : np.ndarray
            Concrete per-variant weights aligned to genotype columns.

        Returns
        -------
        dict with SKAT-O result fields including optimal rho.
        """
        _n_samples, n_variants = geno.shape

        # Guard: zero-variant matrix cannot be tested
        if n_variants == 0:
            logger.debug(f"Gene {gene} (SKAT-O): n_variants=0; returning p_value=None")
            return {
                "p_value": None,
                "rho": None,
                "n_variants": 0,
                "n_marker_test": 0,
                "warnings": [],
                "p_method": None,
                "p_converged": False,
                "skip_reason": "rank_deficient",
            }

        # Rank check on full matrix
        rank = int(np.linalg.matrix_rank(geno))
        if rank < 2:
            logger.debug(f"Gene {gene} (SKAT-O): rank={rank} < 2; skipping")
            return {
                "p_value": None,
                "rho": None,
                "n_variants": n_variants,
                "n_marker_test": n_variants,
                "warnings": [],
                "p_method": None,
                "p_converged": False,
                "skip_reason": "rank_deficient",
            }

        # Weighted genotype matrix: Z = geno * weights
        geno_weighted = geno * weights[np.newaxis, :]

        # ── Step 1-2: Compute Q_rho and projection-adjusted Z1/sqrt(2) ──
        q_rho, z1_half = compute_skato_q_and_z1(geno_weighted, null_model, _SKATO_RHO_GRID)

        # ── Step 3-5: Full omnibus p-value ──
        pvalue, pval_each = _skato_get_pvalue(q_rho, z1_half, _SKATO_RHO_GRID)

        # Find optimal rho (rho with minimum per-rho p-value)
        best_idx = int(np.argmin(pval_each))
        best_rho = _SKATO_RHO_GRID[best_idx]
        if best_rho >= 0.999:
            best_rho = 1.0

        logger.debug(
            f"Gene {gene} SKAT-O: p={pvalue:.6e}, rho={best_rho}, "
            f"per-rho: {[f'{p:.2e}' for p in pval_each]}"
        )

        return {
            "p_value": float(pvalue) if not np.isnan(pvalue) else None,
            "rho": best_rho,
            "n_variants": n_variants,
            "n_marker_test": n_variants,
            "warnings": [],
            "p_method": "optimal.adj",
            "p_converged": True,
            "skip_reason": None,
        }

    def cleanup(self) -> None:
        """
        No-op cleanup (no external resources to release).

        Logs timing summary at INFO level.
        """
        elapsed = time.time() - self._start_time
        logger.info(
            f"PythonSKATBackend.cleanup: {self._genes_processed} genes processed in {elapsed:.1f}s"
        )
        logger.debug("Python SKAT cleanup (no-op)")
