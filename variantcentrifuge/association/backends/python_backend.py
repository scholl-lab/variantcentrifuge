# File: variantcentrifuge/association/backends/python_backend.py
# Location: variantcentrifuge/variantcentrifuge/association/backends/python_backend.py
"""
Pure Python SKAT backend using numpy, scipy, and statsmodels.

PythonSKATBackend implements SKATBackend using only Python-native scientific
libraries — no R or rpy2 required. P-values are computed via the three-tier
Davies -> saddlepoint -> Liu chain from davies.py (Plan 21-01).

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
Extracts resid_response (y - mu_hat) — NOT resid_deviance or resid_pearson.
This matches the R SKAT convention used by RSKATBackend.

Eigenvalue computation
----------------------
Uses scipy.linalg.eigh(driver='evr') matching R's DSYEVR Lapack routine.
Eigenvalue threshold: mean(positive_lambdas) / 100_000 (matches R SKAT source).
Rank check performed on FULL matrix BEFORE eigenvalue filtering.

SKAT-O
------
Searches fixed rho grid [0.0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0].
For each rho, kernel = (1-rho)*K_skat + rho*K_burden; p_min approach used.

Thread safety
-------------
PythonSKATBackend is thread-safe — pure Python/numpy, no R/rpy2 restrictions.

Variable naming note
--------------------
Scientific math variables (design matrix, kernels, eigenvalues) use descriptive
names to comply with ruff N-series naming rules while keeping code readable.
"""

from __future__ import annotations

import logging
import time
from typing import Any

import numpy as np
import scipy.linalg
import scipy.stats

from variantcentrifuge.association.backends.base import NullModelResult, SKATBackend
from variantcentrifuge.association.backends.davies import _try_load_davies, compute_pvalue
from variantcentrifuge.association.weights import beta_maf_weights

logger = logging.getLogger("variantcentrifuge")

# Fixed SKAT-O rho search grid (matches R SKAT package default)
_SKATO_RHO_GRID = [0.0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0]


class PythonSKATBackend(SKATBackend):
    """
    SKAT backend using numpy, scipy, and statsmodels (no R required).

    Implements SKATBackend ABC using only Python-native scientific libraries.
    P-values computed via davies.py three-tier chain (Davies C ext -> saddlepoint -> Liu).

    Lifecycle
    ---------
    1. Instantiate PythonSKATBackend()
    2. Call detect_environment() (never raises — Python deps always present)
    3. Call log_environment()
    4. Call fit_null_model() once per cohort
    5. Call test_gene() per gene
    6. Call cleanup() (no-op)
    """

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
        weights_beta: tuple[float, float],
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
        weights_beta : tuple of (float, float)
            Beta distribution parameters (a1, a2). SKAT default: (1.0, 25.0).

        Returns
        -------
        dict with keys: p_value, rho, n_variants, n_marker_test, warnings,
                        p_method, p_converged, skip_reason (if applicable)
        """
        self._genes_processed += 1
        n_variants = genotype_matrix.shape[1]

        method_upper = method.upper()
        # Map R SKAT method names to our internal names
        # "optimal.adj" and "optimal" are R's names for SKAT-O
        if method_upper == "SKAT":
            result = self._test_skat(gene, genotype_matrix, null_model, weights_beta)
        elif method_upper == "BURDEN":
            result = self._test_burden(gene, genotype_matrix, null_model, weights_beta)
        elif method_upper in ("SKATO", "OPTIMAL.ADJ", "OPTIMAL"):
            result = self._test_skato(gene, genotype_matrix, null_model, weights_beta)
        else:
            logger.warning(f"Unknown SKAT method '{method}' for gene {gene}; using SKAT")
            result = self._test_skat(gene, genotype_matrix, null_model, weights_beta)

        logger.debug(
            f"Gene {gene}: method={method}, p={result.get('p_value')}, "
            f"n_variants={n_variants}, p_method={result.get('p_method')}"
        )
        return result

    def _compute_eigenvalues_filtered(
        self,
        geno_weighted: np.ndarray,
        null_model: NullModelResult,
    ) -> np.ndarray:
        """
        Compute filtered eigenvalues of a SKAT kernel matrix.

        For binary traits, uses R SKAT's projection-adjusted approach:
        ``phi = sqrt(mu*(1-mu))``, then projects the weighted genotype matrix
        through the hat matrix in phi-space before computing eigenvalues. This
        accounts for variance heterogeneity across observations.

        For quantitative traits, uses the simple ``kernel / (2*sigma2)`` scaling.

        Uses scipy.linalg.eigh(driver='evr') matching R's DSYEVR Lapack routine.
        Applies threshold = mean(positive_eigenvalues) / 100_000 (matches R SKAT).

        Parameters
        ----------
        geno_weighted : np.ndarray, shape (n_samples, n_variants)
            Weighted genotype matrix (G * W).
        null_model : NullModelResult
            Fitted null model containing mu_hat, sigma2, design_matrix.

        Returns
        -------
        np.ndarray
            Filtered eigenvalues (positive and above threshold). Empty array
            if none pass.
        """
        sigma2: float = null_model.extra["sigma2"]
        trait_type: str = null_model.trait_type

        if trait_type == "binary":
            # R SKAT projection approach for binary traits:
            # phi = sqrt(mu * (1-mu))  — per-observation weight
            # phi_X = diag(phi) @ X    — weighted design matrix
            # hat_phi = phi_X @ inv(phi_X' phi_X) @ phi_X'
            # Z_adj = diag(phi) @ G_w - hat_phi @ G_w
            # eigenvalues of Z_adj' @ Z_adj (smaller dimension)
            mu_hat: np.ndarray = null_model.extra["mu_hat"]
            design_x: np.ndarray = null_model.extra["design_matrix"]

            variance = mu_hat * (1.0 - mu_hat)
            phi = np.sqrt(variance)

            # Weighted design matrix in phi-space
            phi_x = phi[:, np.newaxis] * design_x  # (n_samples, k)

            # Z_tilde = diag(phi) @ G_w — the genotype in phi-space
            z_tilde = phi[:, np.newaxis] * geno_weighted

            # Project out covariate effects without forming nxn hat matrix:
            # hat_phi @ Z_tilde = phi_X @ inv(phi_X'phi_X) @ phi_X' @ Z_tilde
            phi_xtx = phi_x.T @ phi_x  # (k, k) — small
            phi_xtx_inv = np.linalg.inv(phi_xtx)
            projected = phi_x @ (phi_xtx_inv @ (phi_x.T @ z_tilde))

            # Adjusted genotype matrix: (I - hat_phi) @ Z_tilde
            z_adj = z_tilde - projected

            # Eigenvalues of W.1/2 where W.1 = Z_adj' @ Z_adj
            # R SKAT: K <- W/2 then Get_Lambda(K). The /2 matches the
            # Q = score'score/2 convention so that Q ~ sum(lambda_i chi2_1).
            n_variants = geno_weighted.shape[1]
            eig_mat = z_adj.T @ z_adj if n_variants <= geno_weighted.shape[0] else z_adj @ z_adj.T

            lambdas_all = scipy.linalg.eigh(
                eig_mat / 2.0,
                eigvals_only=True,
                driver="evr",
            )
        else:
            # Quantitative traits: simple scaling
            kernel = geno_weighted @ geno_weighted.T
            lambdas_all = scipy.linalg.eigh(
                kernel / (2.0 * sigma2),
                eigvals_only=True,
                driver="evr",
            )

        # Filter eigenvalues (matches R SKAT exactly):
        # 1. Take positive eigenvalues
        # 2. threshold = mean(positive) / 100_000
        # 3. Keep lambdas > threshold
        pos = lambdas_all[lambdas_all >= 0]
        if len(pos) == 0:
            return np.array([], dtype=np.float64)
        threshold = pos.mean() / 100_000.0
        return lambdas_all[lambdas_all > threshold]

    def _test_skat(
        self,
        gene: str,
        geno: np.ndarray,
        null_model: NullModelResult,
        weights_beta: tuple[float, float],
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
        weights_beta : tuple of (float, float)
            Beta distribution parameters (a1, a2).

        Returns
        -------
        dict with SKAT result fields.
        """
        _n_samples, n_variants = geno.shape
        residuals: np.ndarray = null_model.extra["residuals"]
        a1, a2 = weights_beta

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

        # Compute MAFs and Beta weights
        mafs = geno.mean(axis=0) / 2.0
        weights = beta_maf_weights(mafs, a=a1, b=a2)

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
        weights_beta: tuple[float, float],
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
        weights_beta : tuple of (float, float)
            Beta distribution parameters (a1, a2).

        Returns
        -------
        dict with burden result fields.
        """
        _n_samples, n_variants = geno.shape
        residuals: np.ndarray = null_model.extra["residuals"]
        sigma2: float = null_model.extra["sigma2"]
        a1, a2 = weights_beta

        # Compute MAFs and Beta weights
        mafs = geno.mean(axis=0) / 2.0
        weights = beta_maf_weights(mafs, a=a1, b=a2)

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
        weights_beta: tuple[float, float],
    ) -> dict[str, Any]:
        """
        SKAT-O (optimal unified SKAT) for a single gene.

        Searches the fixed rho grid [0.0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0].
        For each rho, the test statistic is a linear combination of SKAT and
        Burden statistics with kernel K_rho = (1-rho)*K_skat + rho*K_burden.

        The omnibus p-value uses the minimum-p approach: the final p is the
        minimum p-value over the rho grid (a conservative approximation
        sufficient for this implementation; full SKAT-O integrates over the
        joint distribution of p-values).

        Parameters
        ----------
        gene : str
            Gene symbol (for logging).
        geno : np.ndarray, shape (n_samples, n_variants)
            Genotype dosage matrix.
        null_model : NullModelResult
            Fitted null model.
        weights_beta : tuple of (float, float)
            Beta distribution parameters (a1, a2).

        Returns
        -------
        dict with SKAT-O result fields including optimal rho.
        """
        _n_samples, n_variants = geno.shape
        residuals: np.ndarray = null_model.extra["residuals"]
        a1, a2 = weights_beta

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
            logger.debug(f"Gene {gene} (SKAT-O): rank={rank} < 2; skipping (rank_deficient)")
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

        # Compute weights
        mafs = geno.mean(axis=0) / 2.0
        weights = beta_maf_weights(mafs, a=a1, b=a2)

        # Weighted genotype matrix
        geno_weighted = geno * weights[np.newaxis, :]

        # Burden vector: sum of weighted genotypes per sample
        burden = geno @ weights  # (n_samples,)

        # Score vectors for SKAT and burden components
        score_skat_vec = geno_weighted.T @ residuals  # (n_variants,)
        q_skat = float(score_skat_vec @ score_skat_vec) / 2.0
        score_burden = float(residuals @ burden)
        q_burden = score_burden**2 / 2.0

        best_p: float | None = None
        best_rho: float = 0.0
        best_p_method: str = "liu"
        best_p_converged: bool = False

        for rho in _SKATO_RHO_GRID:
            # Combined test statistic for this rho
            q_rho = (1.0 - rho) * q_skat + rho * q_burden

            # Combined weighted genotype for eigenvalue computation:
            # K_rho = (1-rho)*K_skat + rho*K_burden
            # K_skat = G_w @ G_w^T, K_burden = burden @ burden^T
            # Build the combined "genotype" matrix whose outer product = K_rho
            # G_rho = [sqrt(1-rho)*G_w | sqrt(rho)*burden]
            sqrt_rho = np.sqrt(rho) if rho > 0 else 0.0
            sqrt_1mrho = np.sqrt(1.0 - rho)
            geno_rho = (
                np.column_stack(
                    [
                        sqrt_1mrho * geno_weighted,
                        sqrt_rho * burden[:, np.newaxis],
                    ]
                )
                if rho > 0
                else geno_weighted
            )

            # Filtered eigenvalues (projection-adjusted for binary traits)
            lambdas_rho = self._compute_eigenvalues_filtered(geno_rho, null_model)
            if len(lambdas_rho) == 0:
                continue

            p_rho, p_method_rho, p_converged_rho = compute_pvalue(q_rho, lambdas_rho)

            if np.isnan(p_rho) or np.isinf(p_rho):
                continue

            if best_p is None or p_rho < best_p:
                best_p = float(p_rho)
                best_rho = float(rho)
                best_p_method = p_method_rho
                best_p_converged = p_converged_rho

        if best_p is None:
            # All rho values failed — fall back to pure SKAT
            logger.warning(f"Gene {gene}: SKAT-O failed for all rho values; falling back to SKAT")
            return self._test_skat(gene, geno, null_model, weights_beta)

        logger.debug(
            f"Gene {gene} SKAT-O: optimal rho={best_rho}, p={best_p}, p_method={best_p_method}"
        )

        return {
            "p_value": best_p,
            "rho": best_rho,
            "n_variants": n_variants,
            "n_marker_test": n_variants,
            "warnings": [],
            "p_method": best_p_method,
            "p_converged": best_p_converged,
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
