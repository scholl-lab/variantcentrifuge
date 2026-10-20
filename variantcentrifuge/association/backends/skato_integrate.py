# File: variantcentrifuge/association/backends/skato_integrate.py
# Location: variantcentrifuge/association/backends/skato_integrate.py
"""SKAT-O omnibus p-value: mixture parameters and numerical integration over rho."""

from __future__ import annotations

import itertools
import logging
import threading
from typing import Any

import numpy as np
import scipy.integrate
import scipy.linalg
import scipy.stats
from numpy.polynomial.legendre import leggauss

from variantcentrifuge.association.backends.davies import (
    _DAVIES_RETRY_LIM,
    _liu_params,
    _tilted_sf_or_none,
    compute_pvalue,
    davies_pvalue,
)

logger = logging.getLogger("variantcentrifuge")

# Former fixed 16-node rule; unused since the adaptive integral, kept because
# python_backend re-exports it.
_GL_NODES_RAW, _GL_WEIGHTS_RAW = leggauss(16)
_GL_A, _GL_B = 0.0, 40.0  # integration bounds (matches R: upper=40)
# Conditional-statistic levels (in units of E[Q]) at which the SKAT-O integral is split
_SKATO_Z_LEVELS = np.array([0.0, 1.0, 8.0, 64.0])
# Relative accuracy target of the SKAT-O integral; Davies' absolute acc is set to
# this times pmin, but Davies stops converging below ~1e-9 (two eigenvalues), so
# nodes beyond that resolution use the tilted tail probability (_tilted_sf).
_SKATO_DAVIES_REL_ACC, _SKATO_DAVIES_MIN_ACC = 1e-4, 1e-9
_SKATO_X_MAX = 1500.0  # P[chi2_1 > 1500] underflows to 0
_SKATO_QUAD_LIMIT = 50  # max subintervals per piece of the adaptive SKAT-O integral
_SKATO_TAIL_QUAD_LIMIT = 200  # max subintervals of the x > 40 tail integral
_SQRT_2_OVER_PI = float(np.sqrt(2.0 / np.pi))


class _DaviesFailureError(Exception):
    """Davies or the outer quadrature failed for the SKAT-O integral (R: use Liu)."""


# Whether the last SKAT-O p-value in this thread came from the Davies integral (not Liu);
# reported as p_converged without changing _skato_get_pvalue's return value.
_STATUS = threading.local()


def _skato_last_converged() -> bool:
    return bool(getattr(_STATUS, "converged", True))


def _tilted_sf(z: float, lambdas: np.ndarray) -> float:
    """Tail survival with relative accuracy (davies._tilted_sf); failures go to Liu."""
    p = _tilted_sf_or_none(z, lambdas)
    if p is None:
        raise _DaviesFailureError(f"tilted survival failed at z={z:.6g}")
    return p


def _skato_breakpoints(
    a: np.ndarray, b: np.ndarray, levels: np.ndarray, upper: float = _GL_B
) -> np.ndarray:
    """Points in (0, upper) where the SKAT-O integrand needs a quadrature split.

    m(x) = min_i(a_i - b_i*x) has kinks where the minimising line changes. b >= 0, so m is
    non-increasing and crosses each of ``levels`` at most once, at min_i((a_i - level)/b_i);
    levels mark the support edge of the conditional distribution (its CDF is not smooth
    there) and the scale over which that CDF rises, which can be far narrower than a piece.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        cross = (a[:, None] - a[None, :]) / (b[:, None] - b[None, :])
        roots = np.min((a[:, None] - np.asarray(levels)[None, :]) / b[:, None], axis=0)
    cand = np.unique(cross[np.isfinite(cross) & (cross > _GL_A) & (cross < upper)])
    kinks = [
        x for x in cand if np.argmin(a - b * x * (1 - 1e-9)) != np.argmin(a - b * x * (1 + 1e-9))
    ]
    return np.unique(np.r_[kinks, roots[(roots > _GL_A) & (roots < upper)]])


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
        # isf, not ppf(1 - p): 1 - p rounds to 1 (q = inf) for p < 1e-16 (R: qchisq(1-pmin))
        q_org = scipy.stats.chi2.isf(pmin_clamped, df=df_i)
        pmin_q[:, i] = (q_org - df_i) / np.sqrt(2.0 * df_i) * np.sqrt(var_q_i) + mu_q_i

    return pmin, pval, pmin_q


def _skato_integrate_davies(
    pmin_q: np.ndarray,
    param: dict[str, Any],
    rho_grid: list[float],
    pmin: float,
) -> float:
    """
    Compute SKAT-O omnibus p-value via piecewise Gauss-Legendre quadrature.

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

    # The conditional statistic is compared with Q = sum(lambdas * chi2_1): split where it
    # hits 0 (support edge, levels[0]) and multiples of E[Q] (where its CDF rises)
    levels = mu_q + mu_q * (_SKATO_Z_LEVELS - 1.0) / sd1 if sd1 > 0 else np.empty(0)
    a, b = pmin_q[valid] / (1.0 - rho_arr[valid]), tau[valid] / (1.0 - rho_arr[valid])
    lambda_sum = float(np.sum(lambdas))
    # Davies is absolute-accurate; tie its accuracy to the p-value scale (p >= ~pmin)
    acc_target = _SKATO_DAVIES_REL_ACC * pmin
    acc = float(np.clip(acc_target, _SKATO_DAVIES_MIN_ACC, 1e-6))

    def surv(x: np.ndarray) -> np.ndarray:
        """P[Q > g(x)] for the conditional statistic g(x) = (m(x) - mu_q)*sd1 + mu_q."""
        min_q = np.min(a[None, :] - b[None, :] * x[:, None], axis=1)
        z = (min_q - mu_q) * sd1 + mu_q
        if len(lambdas) == 1:
            # Rank one: exact scaled chi2(1); Davies' series fails near the support edge
            return np.asarray(scipy.stats.chi2.sf(z / lambdas[0], df=1))
        out = np.zeros(x.size)
        for k in np.flatnonzero(min_q <= lambda_sum * 1e4):
            # R's lim=10_000 hits ifault=1 near the support edge and R then discards the
            # whole integral for Liu; the larger cap is identical where R converges.
            p_dav, ifault = davies_pvalue(float(z[k]), lambdas, acc=acc, lim=_DAVIES_RETRY_LIM)
            if ifault != 0 and acc < 1e-6:  # tighter accuracy not reachable here: R's acc
                p_dav, ifault = davies_pvalue(float(z[k]), lambdas, acc=1e-6, lim=_DAVIES_RETRY_LIM)
            if p_dav is None or ifault != 0:
                raise _DaviesFailureError(f"ifault={ifault} at q={z[k]:.6g}")
            out[k] = min(max(float(p_dav), 0.0), 1.0)  # Davies error is absolute
            if acc > acc_target and out[k] * _SKATO_DAVIES_REL_ACC < acc and z[k] > lambda_sum:
                # Davies' absolute error cannot meet the relative budget at this node
                out[k] = _tilted_sf(float(z[k]), lambdas)
        return out

    # Survival form: p = int_0^edge f1(x) S(g(x)) dx + P[chi2_1 > edge], with g(x) <= 0
    # (S = 1) beyond the support edge. R computes 1 - int_0^40 f1 (1 - S), which equals
    # this when edge <= 40 but cancels for small p and drops the mass of x > 40.
    edge = float(np.min((a - levels[0]) / b)) if levels.size else np.inf
    upper = min(edge, _GL_B)
    try:
        # Adaptive Gauss-Kronrod per piece in t = sqrt(x) (dx f1(x) = sqrt(2/pi) e^{-t^2/2} dt),
        # with an absolute budget tied to the p scale (p >= ~pmin) and a subdivision cap.
        cuts = np.sqrt(np.r_[_GL_A, _skato_breakpoints(a, b, levels, upper), upper])
        tol = _SKATO_DAVIES_REL_ACC * pmin / len(cuts)
        pvalue = 0.0
        for lo, hi in itertools.pairwise(cuts):
            v, err, *info = scipy.integrate.quad(
                lambda t: float(np.exp(-t * t / 2) * surv(np.array([t * t]))[0]),
                lo,
                hi,
                epsabs=tol / _SQRT_2_OVER_PI,
                epsrel=_SKATO_DAVIES_REL_ACC,
                limit=_SKATO_QUAD_LIMIT,
                full_output=1,
            )
            if len(info) > 1:
                logger.warning("SKAT-O integral piece [%g, %g]: %s", lo**2, hi**2, info[1])
                raise _DaviesFailureError(f"outer quadrature did not converge (err {err:.2g})")
            pvalue += _SQRT_2_OVER_PI * v
        tail_hi = min(edge, _SKATO_X_MAX)
        # Mass beyond x = 40 is < 2.6e-10; integrate it only where it matters relative to p
        if tail_hi > _GL_B and scipy.stats.chi2.sf(_GL_B, 1) > _SKATO_DAVIES_REL_ACC * pvalue:
            pts = _skato_breakpoints(a, b, levels, tail_hi)
            tail, err, *info = scipy.integrate.quad(
                lambda x: float(surv(np.array([x]))[0] * scipy.stats.chi2.pdf(x, df=1)),
                _GL_B,
                tail_hi,
                points=pts[pts > _GL_B] if np.any(pts > _GL_B) else None,
                epsabs=0.0,
                epsrel=_SKATO_DAVIES_REL_ACC,
                limit=_SKATO_TAIL_QUAD_LIMIT,
                full_output=1,
            )
            if len(info) > 1:
                logger.warning("SKAT-O integral tail x>%g: %s", _GL_B, info[1])
                raise _DaviesFailureError(f"tail quadrature did not converge (err {err:.2g})")
            logger.debug("SKAT-O tail x>%g: %.3g (quad err %.2g)", _GL_B, tail, err)
            pvalue += tail
    except _DaviesFailureError as exc:
        logger.debug("SKAT-O Davies %s; using Liu", exc)
        return _skato_integrate_liu(pmin_q, param, rho_grid, pmin)
    pvalue += float(scipy.stats.chi2.sf(edge, df=1))

    if np.isnan(pvalue):
        return _skato_integrate_liu(pmin_q, param, rho_grid, pmin)

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
    _STATUS.converged = False
    tau = param["tau"]
    mu_q = param["mu_q"]
    var_q = param["var_q"]
    df = param["df"]
    rho_arr = np.array(rho_grid)
    # R integrates over every (capped, < 1) rho; dropping rho=0.999 lost the burden term
    valid = rho_arr < 1.0
    if not np.any(valid):
        return min(1.0, pmin * len(rho_grid))
    a, b = pmin_q[valid] / (1.0 - rho_arr[valid]), tau[valid] / (1.0 - rho_arr[valid])
    scale = np.sqrt(var_q) / np.sqrt(2.0 * df)

    def integrand(x: float) -> float:
        q_std = (float(np.min(a - b * x)) - mu_q) / scale + df
        return float(scipy.stats.chi2.sf(q_std, df=df) * scipy.stats.chi2.pdf(x, df=1))

    # Survival form as in _skato_integrate_davies: q_std <= 0 (survival 1) beyond the edge.
    # Split where q_std hits 0 and multiples of its mean df.
    levels = mu_q + (df * _SKATO_Z_LEVELS - df) * scale
    edge = float(np.min((a - levels[0]) / b))
    upper = min(edge, _SKATO_X_MAX)
    points = _skato_breakpoints(a, b, levels, upper)
    integral, _err = scipy.integrate.quad(
        integrand,
        0,
        upper,
        points=points if points.size else None,
        limit=2000,
        epsabs=0.0,
        epsrel=1e-10,  # cheap integrand
    )

    pvalue: float = float(integral) + float(scipy.stats.chi2.sf(edge, df=1))

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
    _STATUS.converged = True
    n_rho = len(rho_grid)
    n_variants = z1_half.shape[1]

    # Cap rho at 0.999 to match R's SKAT_Optimal_Logistic and avoid singular R_rho at rho=1
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
        # Analytic sqrt R_rho^{1/2} = s*I + delta*J avoids Cholesky failure near rho=1
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
