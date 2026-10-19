# File: variantcentrifuge/association/backends/skat_math.py
# Location: variantcentrifuge/association/backends/skat_math.py
"""Mathematical routines and numerical integration for Python SKAT / SKAT-O."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import numpy as np
import scipy.integrate
import scipy.linalg
import scipy.stats
from numpy.polynomial.legendre import leggauss

from variantcentrifuge.association.backends.davies import (
    _liu_params,
    compute_pvalue,
    davies_pvalue,
)

if TYPE_CHECKING:
    from variantcentrifuge.association.backends.base import NullModelResult

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

        pmin_clamped = np.clip(pmin, 1e-300, 1.0 - 1e-15)
        q_org = scipy.stats.chi2.ppf(1.0 - pmin_clamped, df=df_i)
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

    Integrates over chi-squared(1) distribution, conditioning on the shared
    burden component and computing the minimum conditional probability across rho.
    """
    tau = param["tau"]
    mu_q = param["mu_q"]
    var_q = param["var_q"]
    var_remain = param["var_remain"]
    lambdas = param["lambdas"]
    rho_arr = np.array(rho_grid)

    var_ratio = (var_q - var_remain) / var_q if var_q > 0 else 0.0
    sd1 = np.sqrt(max(var_ratio, 0.0))

    valid = rho_arr < 1.0
    if not np.any(valid):
        return _skato_integrate_liu(pmin_q, param, rho_grid, pmin)

    pmin_q_valid = pmin_q[valid]
    tau_valid = tau[valid]
    lambda_sum = float(np.sum(lambdas))

    integrand_vals = np.zeros(128, dtype=np.float64)
    for k in range(128):
        x = float(_GL_X[k])

        cond_q = (pmin_q_valid - tau_valid * x) / (1.0 - rho_arr[valid])
        min_q = float(np.min(cond_q))

        if min_q > lambda_sum * 1e4:
            cdf_val = 0.0
        else:
            min_q_centered = min_q - mu_q
            min_q_std = min_q_centered * sd1 + mu_q

            p_dav, ifault = davies_pvalue(min_q_std, lambdas, acc=1e-6, lim=10_000)
            if p_dav is not None and ifault == 0:
                cdf_val = float(p_dav)
            else:
                return _skato_integrate_liu(pmin_q, param, rho_grid, pmin)

            if cdf_val > 1.0:
                cdf_val = 1.0

        integrand_vals[k] = (1.0 - cdf_val) * float(scipy.stats.chi2.pdf(x, df=1))

    integral = float(np.dot(integrand_vals, _GL_W))

    if np.isnan(integral):
        return _skato_integrate_liu(pmin_q, param, rho_grid, pmin)

    pvalue: float = 1.0 - integral

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

    rho_capped = [min(r, 0.999) for r in rho_grid]

    if rho_capped != list(rho_grid):
        q_skat_half = q_rho[0]
        q_burden_half = q_rho[-1]
        q_rho = np.array([(1 - r) * q_skat_half + r * q_burden_half for r in rho_capped])

    a_mat = z1_half.T @ z1_half
    a_col_sums = a_mat.sum(axis=0)
    a_row_sums = a_mat.sum(axis=1)

    lambda_all: list[np.ndarray] = []
    for _i, rho in enumerate(rho_capped):
        s = np.sqrt(1.0 - rho)
        lam_burden = 1.0 - rho + n_variants * rho
        s_burden = np.sqrt(lam_burden)
        delta = (s_burden - s) / n_variants

        k_sym = s * s * a_mat
        k_sym += s * delta * (a_col_sums[np.newaxis, :] + a_row_sums[:, np.newaxis])
        k_sym += delta * delta * float(a_mat.sum())

        lambda_all.append(_get_lambda(k_sym))

    param = _skato_optimal_param(z1_half, rho_capped)

    q_all_2d = q_rho[np.newaxis, :]
    pmin, pval, pmin_q = _skato_each_q(param, q_all_2d, rho_capped, lambda_all)

    p_min_scalar = float(pmin[0])
    pval_each = pval[0, :]

    pvalue = _skato_integrate_davies(pmin_q[0, :], param, rho_capped, p_min_scalar)

    multi = 3 if len(rho_capped) >= 3 else 2
    pval_each_pos = pval_each[pval_each > 0]

    if pvalue <= 0 or len(pval_each_pos) < n_rho:
        pvalue = float(np.min(pval_each)) * multi

    if pvalue == 0 and len(pval_each_pos) > 0:
        pvalue = float(np.min(pval_each_pos))

    return float(np.clip(pvalue, 0.0, 1.0)), pval_each


def compute_eigenvalues_filtered(
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
        mu_hat: np.ndarray = null_model.extra["mu_hat"]
        design_x: np.ndarray = null_model.extra["design_matrix"]

        variance = mu_hat * (1.0 - mu_hat)
        phi = np.sqrt(variance)

        phi_x = phi[:, np.newaxis] * design_x
        z_tilde = phi[:, np.newaxis] * geno_weighted

        phi_xtx = phi_x.T @ phi_x
        phi_xtx_inv = np.linalg.inv(phi_xtx)
        projected = phi_x @ (phi_xtx_inv @ (phi_x.T @ z_tilde))

        z_adj = z_tilde - projected

        n_variants = geno_weighted.shape[1]
        eig_mat = z_adj.T @ z_adj if n_variants <= geno_weighted.shape[0] else z_adj @ z_adj.T

        lambdas_all = scipy.linalg.eigh(
            eig_mat / 2.0,
            eigvals_only=True,
            driver="evr",
        )
    else:
        kernel = geno_weighted @ geno_weighted.T
        lambdas_all = scipy.linalg.eigh(
            kernel / (2.0 * sigma2),
            eigvals_only=True,
            driver="evr",
        )

    pos = lambdas_all[lambdas_all >= 0]
    if len(pos) == 0:
        return np.array([], dtype=np.float64)
    threshold = pos.mean() / 100_000.0
    return np.asarray(lambdas_all[lambdas_all > threshold], dtype=np.float64)


def compute_skato_q_and_z1(
    geno_weighted: np.ndarray,
    null_model: NullModelResult,
    rho_grid: list[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-rho Q statistics and projection-adjusted Z1/sqrt(2) for SKAT-O."""
    residuals: np.ndarray = null_model.extra["residuals"]
    trait_type: str = null_model.trait_type

    # ── Step 1: Compute Q_rho for each rho ──
    score_vec = geno_weighted.T @ residuals
    score_sum = float(np.sum(score_vec))

    q_rho = np.zeros(len(rho_grid), dtype=np.float64)
    for i, rho in enumerate(rho_grid):
        q1 = (1.0 - rho) * float(score_vec @ score_vec)
        q2 = rho * score_sum**2
        q_rho[i] = (q1 + q2) / 2.0

    # ── Step 2: Build projection-adjusted Z1 ──
    if trait_type == "binary":
        mu_hat: np.ndarray = null_model.extra["mu_hat"]
        design_x: np.ndarray = null_model.extra["design_matrix"]
        pi_1 = mu_hat * (1.0 - mu_hat)
        sqrt_pi = np.sqrt(pi_1)

        z_phi = geno_weighted * sqrt_pi[:, np.newaxis]
        x_phi = design_x * sqrt_pi[:, np.newaxis]
        x_pi = design_x * pi_1[:, np.newaxis]
        xtx_pi = x_phi.T @ x_phi
        xtx_pi_inv = np.linalg.inv(xtx_pi)
        z1_adj = z_phi - x_phi @ (xtx_pi_inv @ (x_pi.T @ geno_weighted))
    else:
        design_x = null_model.extra["design_matrix"]
        sigma2: float = null_model.extra["sigma2"]
        xtx = design_x.T @ design_x
        xtx_inv = np.linalg.inv(xtx)
        z1_adj = geno_weighted - design_x @ (xtx_inv @ (design_x.T @ geno_weighted))
        q_rho /= sigma2

    z1_half = z1_adj / np.sqrt(2.0)
    return q_rho, z1_half
