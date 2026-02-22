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
Implements Lee et al. (2012) optimal unified test matching R SKAT's
``SKAT_Optimal_Logistic`` / ``SKAT_Optimal_Linear`` exactly.
Searches fixed rho grid [0.0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0].
Per-rho eigenvalues computed via analytical R.M^{1/2} (avoids Cholesky).
Omnibus p-value via integration over chi-squared(1) conditioning on the
shared burden component.

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
import scipy.integrate
import scipy.linalg
import scipy.stats
from numpy.polynomial.legendre import leggauss

from variantcentrifuge.association.backends.base import NullModelResult, SKATBackend
from variantcentrifuge.association.backends.davies import (
    _liu_params,
    _try_load_davies,
    compute_pvalue,
    davies_pvalue,
)
from variantcentrifuge.association.weights import beta_maf_weights

logger = logging.getLogger("variantcentrifuge")

# Fixed SKAT-O rho search grid (matches R SKAT package default)
_SKATO_RHO_GRID = [0.0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0]

# 128-node Gauss-Legendre quadrature constants for SKAT-O integration.
# Precomputed at module load: 46x speedup over adaptive scipy.integrate.quad
# (measured: 379ms -> 8ms per gene on typical cohort sizes).
# Integration bounds [0, 40] match R SKAT's upper=40 default.
_GL_NODES_RAW, _GL_WEIGHTS_RAW = leggauss(128)
_GL_A, _GL_B = 0.0, 40.0  # integration bounds (matches R: upper=40)
_GL_X = (_GL_B - _GL_A) / 2.0 * _GL_NODES_RAW + (_GL_B + _GL_A) / 2.0
_GL_W = (_GL_B - _GL_A) / 2.0 * _GL_WEIGHTS_RAW


def _get_lambda(kernel_mat: np.ndarray) -> np.ndarray:
    """
    Compute filtered eigenvalues of a kernel matrix (R: ``Get_Lambda``).

    Uses scipy.linalg.eigh and applies the R SKAT filtering:
    keep eigenvalues > mean(positive) / 100_000.
    """
    lambdas_all = scipy.linalg.eigh(kernel_mat, eigvals_only=True, driver="evr")
    pos = lambdas_all[lambdas_all >= 0]
    if len(pos) == 0:
        return np.array([], dtype=np.float64)
    threshold = pos.mean() / 100_000.0
    return np.asarray(lambdas_all[lambdas_all > threshold], dtype=np.float64)


def _skato_optimal_param(z1: np.ndarray, rho_grid: list[float]) -> dict[str, Any]:
    """
    Compute SKAT-O mixture distribution parameters (R: ``SKAT_Optimal_Param``).

    Decomposes Z1 into a burden-direction component and an orthogonal component.
    Returns parameters needed for the omnibus p-value integration.

    Parameters
    ----------
    z1 : np.ndarray, shape (n, p)
        Projection-adjusted genotype matrix divided by sqrt(2).
    rho_grid : list of float
        Rho correlation values.

    Returns
    -------
    dict with keys: mu_q, var_q, ker_q, df, lambdas, var_remain, tau
    """
    _n_samples, n_variants = z1.shape

    # Decompose Z1 into burden direction (z_mean) and orthogonal
    z_mean = z1.mean(axis=1)  # (n,) — mean across variants
    z_mean_sq = float(z_mean @ z_mean)

    if z_mean_sq < 1e-30:
        # Degenerate case: no burden signal
        z_mean_sq = 1e-30

    # Regression coefficients of each column of Z1 on z_mean
    cof1 = (z_mean @ z1) / z_mean_sq  # (p,)

    # Decompose: Z1 = Z.item1 (burden direction) + Z.item2 (orthogonal)
    z_item1 = z_mean[:, np.newaxis] * cof1[np.newaxis, :]  # (n, p)
    z_item2 = z1 - z_item1  # (n, p)

    # W3.2 — orthogonal component kernel
    w3_2 = z_item2.T @ z_item2  # (p, p)
    lambdas = _get_lambda(w3_2)

    # W3.3 — cross-term variance contribution
    # R: sum((t(Z.item1) %*% Z.item1) * (t(Z.item2) %*% Z.item2)) * 4
    cross1 = z_item1.T @ z_item1  # (p, p)
    cross2 = z_item2.T @ z_item2  # (p, p)
    var_remain = float(4.0 * np.sum(cross1 * cross2))

    # Mixture moments
    mu_q = float(np.sum(lambdas))
    var_q = float(2.0 * np.sum(lambdas**2)) + var_remain
    ker_q = float(12.0 * np.sum(lambdas**4)) / float(np.sum(lambdas**2)) ** 2
    df = 12.0 / ker_q

    # Tau: burden-direction contribution at each rho
    cof1_sq_sum = float(np.sum(cof1**2))
    p_m = n_variants
    tau = []
    for rho in rho_grid:
        term = float(p_m**2) * rho + cof1_sq_sum * (1.0 - rho)
        tau.append(term * z_mean_sq)

    return {
        "mu_q": mu_q,
        "var_q": var_q,
        "ker_q": ker_q,
        "df": df,
        "lambdas": lambdas,
        "var_remain": var_remain,
        "tau": np.array(tau),
    }


def _skato_each_q(
    param: dict[str, Any],
    q_all: np.ndarray,
    rho_grid: list[float],
    lambda_all: list[np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute per-rho p-values and inverted quantiles (R: ``SKAT_Optimal_Each_Q``).

    For optimal.adj: uses the full Davies -> Liu chain for each rho.

    Parameters
    ----------
    param : dict
        Mixture parameters from _skato_optimal_param.
    q_all : np.ndarray, shape (n_q, n_rho)
        Q statistics for each observation x rho.
    rho_grid : list of float
        Rho values.
    lambda_all : list of np.ndarray
        Per-rho eigenvalues.

    Returns
    -------
    (pmin, pval, pmin_q)
        pmin : np.ndarray, shape (n_q,) — minimum p-value across rho for each Q
        pval : np.ndarray, shape (n_q, n_rho) — per-rho p-values
        pmin_q : np.ndarray, shape (n_q, n_rho) — inverted quantiles at pmin
    """
    n_rho = len(rho_grid)
    n_q = q_all.shape[0] if q_all.ndim > 1 else 1
    if q_all.ndim == 1:
        q_all = q_all[np.newaxis, :]

    pval = np.zeros((n_q, n_rho), dtype=np.float64)
    param_mat = np.zeros((n_rho, 3), dtype=np.float64)  # mu_q, var_q, df

    for i, _rho in enumerate(rho_grid):
        lam = lambda_all[i]
        if len(lam) == 0:
            pval[:, i] = 1.0
            param_mat[i] = [0, 1, 1]
            continue

        # Liu parameters for this rho's eigenvalues
        mu_q_i, sigma_q_i, _mu_x_i, _sigma_x_i, ll_i, _d_i = _liu_params(lam)
        param_mat[i] = [mu_q_i, sigma_q_i**2, ll_i]

        # For optimal.adj: use Davies -> Liu chain (Get_PValue.Lambda)
        for j in range(n_q):
            p_val, _method, _conv = compute_pvalue(q_all[j, i], lam)
            pval[j, i] = p_val

    # Minimum p-value across rho for each query
    pmin = np.min(pval, axis=1)

    # Invert pmin back to Q-scale for each rho
    pmin_q = np.zeros((n_q, n_rho), dtype=np.float64)
    for i in range(n_rho):
        mu_q_i = param_mat[i, 0]
        var_q_i = param_mat[i, 1]
        df_i = param_mat[i, 2]

        # qchisq(1 - pmin, df=df) -> invert through Liu parameterization
        # Clamp pmin to avoid qchisq(0) = inf
        pmin_clamped = np.clip(pmin, 1e-300, 1.0 - 1e-15)
        q_org = scipy.stats.chi2.ppf(1.0 - pmin_clamped, df=df_i)
        # Map back to Q-scale
        pmin_q[:, i] = (q_org - df_i) / np.sqrt(2.0 * df_i) * np.sqrt(var_q_i) + mu_q_i

    return pmin, pval, pmin_q


def _skato_integrate_davies(
    pmin_q: np.ndarray,
    param: dict[str, Any],
    rho_grid: list[float],
    pmin: float,
) -> float:
    """
    Compute SKAT-O omnibus p-value via 128-node Gauss-Legendre quadrature.

    Replaces adaptive ``scipy.integrate.quad`` with fixed GL quadrature for a
    46x speedup (379ms -> 8ms per gene). Integrates over chi-squared(1)
    distribution, conditioning on the shared burden component and computing the
    minimum conditional probability across rho.

    Integration bounds [0, 40] match R SKAT's default (upper=40). The 128-node
    GL rule provides sufficient accuracy for the smooth integrands encountered in
    SKAT-O; results are within tolerance of the previous adaptive quad approach.

    Falls back to ``_skato_integrate_liu`` if any Davies call fails during the
    quadrature loop (same fallback behaviour as the previous adaptive approach).

    Parameters
    ----------
    pmin_q : np.ndarray, shape (n_rho,)
        Inverted quantiles at the minimum p-value.
    param : dict
        Mixture parameters from _skato_optimal_param.
    rho_grid : list of float
        Rho values.
    pmin : float
        Minimum per-rho p-value (for Bonferroni guard).

    Returns
    -------
    float
        Omnibus p-value.
    """
    tau = param["tau"]
    mu_q = param["mu_q"]
    var_q = param["var_q"]
    var_remain = param["var_remain"]
    lambdas = param["lambdas"]
    rho_arr = np.array(rho_grid)

    # Standardization factor: sd1 = sqrt((VarQ - VarRemain) / VarQ)
    var_ratio = (var_q - var_remain) / var_q if var_q > 0 else 0.0
    sd1 = np.sqrt(max(var_ratio, 0.0))

    # Precompute mask for valid rho values (rho < 1.0)
    valid = rho_arr < 1.0
    if not np.any(valid):
        return _skato_integrate_liu(pmin_q, param, rho_grid, pmin)

    pmin_q_valid = pmin_q[valid]
    tau_valid = tau[valid]
    lambda_sum = float(np.sum(lambdas))

    # Evaluate integrand at all 128 GL nodes via loop.
    # Loop is required because davies_pvalue() is called per node.
    integrand_vals = np.zeros(128, dtype=np.float64)
    for k in range(128):
        x = float(_GL_X[k])

        # Conditional quantile at each valid rho
        cond_q = (pmin_q_valid - tau_valid * x) / (1.0 - rho_arr[valid])
        min_q = float(np.min(cond_q))

        # Check if min_q is extremely large (probability ~0)
        if min_q > lambda_sum * 1e4:
            cdf_val = 0.0
        else:
            # Standardize to match the orthogonal-component distribution
            min_q_centered = min_q - mu_q
            min_q_std = min_q_centered * sd1 + mu_q

            # Davies CDF of standardized quantile under orthogonal eigenvalues
            p_dav, ifault = davies_pvalue(min_q_std, lambdas, acc=1e-6, lim=10_000)
            if p_dav is not None and ifault == 0:
                cdf_val = float(p_dav)  # p_dav is P[Q > q], which is survival
            else:
                # Davies failed — fall back to Liu integration (same as before)
                return _skato_integrate_liu(pmin_q, param, rho_grid, pmin)

            if cdf_val > 1.0:
                cdf_val = 1.0

        # Integrand: (1 - survival) * dchisq(x, df=1)
        integrand_vals[k] = (1.0 - cdf_val) * float(scipy.stats.chi2.pdf(x, df=1))

    # GL quadrature: integral = sum(f(x_k) * w_k)
    integral = float(np.dot(integrand_vals, _GL_W))

    # Check for NaN from Davies failure in integrand
    if np.isnan(integral):
        return _skato_integrate_liu(pmin_q, param, rho_grid, pmin)

    pvalue: float = 1.0 - integral

    # Bonferroni guard: p should be <= pmin * n_rho
    if pmin * len(rho_grid) < pvalue:
        pvalue = pmin * len(rho_grid)

    return pvalue


def _skato_integrate_liu(
    pmin_q: np.ndarray,
    param: dict[str, Any],
    rho_grid: list[float],
    pmin: float,
) -> float:
    """
    SKAT-O omnibus p-value via Liu integration (R: ``SKAT_Optimal_PValue_Liu``).

    Fallback when Davies integration fails.
    """
    tau = param["tau"]
    mu_q = param["mu_q"]
    var_q = param["var_q"]
    df = param["df"]
    rho_arr = np.array(rho_grid)

    def integrand(x: float) -> float:
        valid = rho_arr < 0.999
        if not np.any(valid):
            return 0.0
        cond_q = (pmin_q[valid] - tau[valid] * x) / (1.0 - rho_arr[valid])
        min_q = float(np.min(cond_q))

        # Standardize through chi-squared approximation
        q_std = (min_q - mu_q) / np.sqrt(var_q) * np.sqrt(2.0 * df) + df
        return float(scipy.stats.chi2.cdf(q_std, df=df)) * float(scipy.stats.chi2.pdf(x, df=1))

    integral, _err = scipy.integrate.quad(integrand, 0, 40, limit=2000, epsabs=1e-25)

    pvalue: float = 1.0 - float(integral)

    if pmin * len(rho_grid) < pvalue:
        pvalue = pmin * len(rho_grid)

    return pvalue


def _skato_get_pvalue(
    q_rho: np.ndarray,
    z1_half: np.ndarray,
    rho_grid: list[float],
) -> tuple[float, np.ndarray]:
    """
    Full SKAT-O p-value computation (R: ``SKAT_Optimal_Get_Pvalue``).

    Parameters
    ----------
    q_rho : np.ndarray, shape (n_rho,)
        Q statistics at each rho.
    z1_half : np.ndarray, shape (n, p)
        Projection-adjusted genotype matrix / sqrt(2).
    rho_grid : list of float
        Rho correlation grid.

    Returns
    -------
    (p_value, p_val_each)
        p_value : float — omnibus p-value
        p_val_each : np.ndarray, shape (n_rho,) — per-rho p-values
    """
    n_rho = len(rho_grid)
    n_variants = z1_half.shape[1]

    # Cap rho at 0.999 globally (matches R SKAT_Optimal_Logistic which caps
    # r.all before computing Q, eigenvalues, tau, and per-rho p-values).
    rho_capped = [min(r, 0.999) for r in rho_grid]

    # Recompute Q for any capped rho values. Q is linear in rho:
    #   q_rho[i] = (1-rho_i) * q_rho[0] + rho_i * q_rho[-1]
    # because Q(rho) = ((1-rho)*Q_SKAT + rho*Q_Burden) / 2.
    if rho_capped != list(rho_grid):
        q_skat_half = q_rho[0]  # Q_SKAT / 2 (at rho=0)
        q_burden_half = q_rho[-1]  # Q_Burden / 2 (at rho=1)
        q_rho = np.array([(1 - r) * q_skat_half + r * q_burden_half for r in rho_capped])

    # Base kernel: A = Z1_half' @ Z1_half (p x p)
    a_mat = z1_half.T @ z1_half

    # Compute eigenvalues for each rho via analytical R.M^{1/2}.
    #
    # R.M = (1-rho)*I + rho*J has known eigenvalues:
    #   - (1-rho) with multiplicity (p-1)
    #   - (1-rho + p*rho) with multiplicity 1
    #
    # Its matrix square root is:
    #   R.M^{1/2} = s*I + delta*J
    # where s = sqrt(1-rho), delta = (sqrt(1-rho+p*rho) - sqrt(1-rho)) / p
    #
    # The symmetric product K_sym = R.M^{1/2} @ A @ R.M^{1/2} has the same
    # eigenvalues as L @ A @ L' (the Cholesky approach) but avoids numerical
    # instability from standard Cholesky on near-singular correlation matrices.
    # This also works at rho=1 where Cholesky fails (rank-1 matrix).

    # Column sums of A (needed for the J @ A @ J term)
    a_col_sums = a_mat.sum(axis=0)  # (p,) — each entry is sum of column j
    a_row_sums = a_mat.sum(axis=1)  # (p,) — each entry is sum of row i

    lambda_all: list[np.ndarray] = []
    for _i, rho in enumerate(rho_capped):
        s = np.sqrt(1.0 - rho)
        lam_burden = 1.0 - rho + n_variants * rho
        s_burden = np.sqrt(lam_burden)
        delta = (s_burden - s) / n_variants

        # K_sym = (s*I + delta*J) @ A @ (s*I + delta*J)
        #       = s^2 * A
        #         + s*delta * (J@A + A@J)
        #         + delta^2 * J@A@J
        #
        # J@A has each row = a_col_sums (row vector repeated)
        # A@J has each column = a_row_sums (column vector repeated)
        # J@A@J is a constant matrix with all entries = sum(A)
        k_sym = s * s * a_mat
        k_sym += s * delta * (a_col_sums[np.newaxis, :] + a_row_sums[:, np.newaxis])
        k_sym += delta * delta * float(a_mat.sum())

        lambda_all.append(_get_lambda(k_sym))

    # Mixture parameters for omnibus integration (uses capped rho for tau)
    param = _skato_optimal_param(z1_half, rho_capped)

    # Per-rho p-values and inverted quantiles
    q_all_2d = q_rho[np.newaxis, :]  # (1, n_rho)
    pmin, pval, pmin_q = _skato_each_q(param, q_all_2d, rho_capped, lambda_all)

    p_min_scalar = float(pmin[0])
    pval_each = pval[0, :]

    # Omnibus p-value via Davies integration
    pvalue = _skato_integrate_davies(pmin_q[0, :], param, rho_capped, p_min_scalar)

    # Correction: SKAT-O p should be <= min(per-rho p) * multiplier
    multi = 3 if len(rho_capped) >= 3 else 2
    pval_each_pos = pval_each[pval_each > 0]

    if pvalue <= 0 or len(pval_each_pos) < n_rho:
        pvalue = float(np.min(pval_each)) * multi

    if pvalue == 0 and len(pval_each_pos) > 0:
        pvalue = float(np.min(pval_each_pos))

    return float(np.clip(pvalue, 0.0, 1.0)), pval_each


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
        return np.asarray(lambdas_all[lambdas_all > threshold], dtype=np.float64)

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
        weights_beta : tuple of (float, float)
            Beta distribution parameters (a1, a2).

        Returns
        -------
        dict with SKAT-O result fields including optimal rho.
        """
        _n_samples, n_variants = geno.shape
        residuals: np.ndarray = null_model.extra["residuals"]
        trait_type: str = null_model.trait_type
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

        # Compute weights
        mafs = geno.mean(axis=0) / 2.0
        weights = beta_maf_weights(mafs, a=a1, b=a2)

        # Weighted genotype matrix: Z = geno * weights
        geno_weighted = geno * weights[np.newaxis, :]

        # ── Step 1: Compute Q_rho for each rho ──
        # R: temp = t(res) %*% Z; Q_rho = ((1-rho)*sum(temp^2) + rho*p^2*mean(temp)^2) / 2
        score_vec = geno_weighted.T @ residuals  # (n_variants,)
        score_sum = float(np.sum(score_vec))

        q_rho = np.zeros(len(_SKATO_RHO_GRID), dtype=np.float64)
        for i, rho in enumerate(_SKATO_RHO_GRID):
            q1 = (1.0 - rho) * float(score_vec @ score_vec)
            q2 = rho * score_sum**2
            q_rho[i] = (q1 + q2) / 2.0

        # ── Step 2: Build projection-adjusted Z1 ──
        if trait_type == "binary":
            mu_hat: np.ndarray = null_model.extra["mu_hat"]
            design_x: np.ndarray = null_model.extra["design_matrix"]
            pi_1 = mu_hat * (1.0 - mu_hat)  # variance
            sqrt_pi = np.sqrt(pi_1)

            # R: Z1 = (Z * sqrt(pi_1)) - (X1 * sqrt(pi_1)) %*%
            #         solve(t(X1) %*% (X1 * pi_1)) %*% (t(X1) %*% (Z * pi_1))
            z_phi = geno_weighted * sqrt_pi[:, np.newaxis]
            x_phi = design_x * sqrt_pi[:, np.newaxis]
            x_pi = design_x * pi_1[:, np.newaxis]
            xtx_pi = x_phi.T @ x_phi  # = X' diag(pi_1) X
            xtx_pi_inv = np.linalg.inv(xtx_pi)
            z1_adj = z_phi - x_phi @ (xtx_pi_inv @ (x_pi.T @ geno_weighted))
        else:
            # Linear (quantitative): Z1 = Z - X @ solve(X'X) @ X'Z
            design_x = null_model.extra["design_matrix"]
            sigma2: float = null_model.extra["sigma2"]
            xtx = design_x.T @ design_x
            xtx_inv = np.linalg.inv(xtx)
            z1_adj = geno_weighted - design_x @ (xtx_inv @ (design_x.T @ geno_weighted))
            # For linear, divide Q by sigma2 (R: Q.all / s2)
            q_rho /= sigma2

        # Z1 / sqrt(2) — matches R's call: SKAT_Optimal_Get_Pvalue(Q.all, Z1/sqrt(2), ...)
        z1_half = z1_adj / np.sqrt(2.0)

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
