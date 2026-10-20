"""
Three-tier p-value computation for quadratic forms of chi-squared variables.

This module implements the fallback chain used by GENESIS/GMMAT for computing
P[Q > q] where Q = sum(lambda_j * chi2_1_j):

    Tier 1 — Davies (exact):   CompQuadForm qfc() C++ extension via CFFI.
    Tier 2 — Kuonen (saddle):  Lugannani-Rice saddlepoint (Kuonen 1999).
    Tier 3 — Liu (moments):    Liu moment-matching (Liu et al. 2009).

The Davies tier is only available when the _qfc C extension has been compiled
(via ``pip install`` with a C++ compiler present). When unavailable, the chain
starts at Tier 2.

Davies runs to a ~1e-4 relative accuracy budget. Below p = 1e-5 (where Davies'
absolute accuracy floor cannot meet that budget) it is evaluated on the exponentially
tilted mixture ("davies_tilted"). Saddlepoint, then Liu, are fallbacks only, used when
Davies faults (ifault > 0) or returns p >= 1 / non-finite; they report p_converged=False.

References
----------
Davies, R.B. (1980). The Distribution of a Linear Combination of chi-squared
    Random Variables. JRSS-C 29(3):323-333.
Kuonen, D. (1999). Saddlepoint approximations for distributions of quadratic
    forms in normal variables. Biometrika 86(4):929-935.
Liu, H., Tang, Y., Zhang, H. (2009). A new chi-square approximation to the
    distribution of non-negative definite quadratic forms in non-central normal
    variables. Computational Statistics & Data Analysis 53(2):853-856.
"""

from __future__ import annotations

import logging
import os
import threading

import numpy as np
from numpy.polynomial.laguerre import laggauss
from scipy.optimize import brentq
from scipy.stats import ncx2, norm

logger = logging.getLogger("variantcentrifuge")

# ---------------------------------------------------------------------------
# Davies C extension — lazy loading
# ---------------------------------------------------------------------------

_qfc_lib = None  # The loaded CFFI lib object
_qfc_ffi = None  # The CFFI FFI object
_DAVIES_AVAILABLE = False
_DAVIES_LOAD_ATTEMPTED = False

# Gauss-Laguerre rule for the exponentially tilted tail probability (_tilted_sf_or_none)
_LAG_NODES, _LAG_WEIGHTS = laggauss(16)
# qfc keeps its state in C globals and the CFFI call releases the GIL: serialize calls
_QFC_LOCK = threading.Lock()
# Relative accuracy budget for compute_pvalue. Davies' error is absolute and it stops
# converging below acc ~1e-9, so below p = 1e-9 / 1e-4 = 1e-5 the tilted evaluation is used.
_DAVIES_REL_ACC, _DAVIES_MIN_ACC = 1e-4, 1e-9
_TAIL_P = _DAVIES_MIN_ACC / _DAVIES_REL_ACC
# Davies term cap for a retry after ifault=1 (term limit hit). qfc's result does not depend
# on lim when the smaller cap converges, but lim=10_000 (R) fails for small q with few or
# dominant eigenvalues, where saddlepoint/Liu are off by ~1%.
_DAVIES_RETRY_LIM = 1_000_000


def _try_load_davies() -> bool:
    """
    Attempt to load the compiled _qfc CFFI extension.

    Returns True if the extension is available, False otherwise.
    Results are cached after the first attempt.

    The ``VARIANTCENTRIFUGE_NO_C_EXT`` environment variable disables loading
    (useful for testing pure-Python fallback paths).
    """
    global _qfc_lib, _qfc_ffi, _DAVIES_AVAILABLE, _DAVIES_LOAD_ATTEMPTED

    if _DAVIES_LOAD_ATTEMPTED:
        return _DAVIES_AVAILABLE

    _DAVIES_LOAD_ATTEMPTED = True

    if os.environ.get("VARIANTCENTRIFUGE_NO_C_EXT"):
        logger.debug(
            "VARIANTCENTRIFUGE_NO_C_EXT set: Davies C extension disabled. "
            "Using Kuonen saddlepoint / Liu moment-matching fallback."
        )
        return False

    try:
        from variantcentrifuge import _qfc as _qfc_mod  # type: ignore[attr-defined]

        _qfc_lib = _qfc_mod.lib
        _qfc_ffi = _qfc_mod.ffi
        _DAVIES_AVAILABLE = True
        logger.debug("Davies qfc C extension loaded successfully.")
        return True
    except ImportError:
        logger.debug(
            "Davies qfc C extension not available (_qfc not compiled). "
            "Using Kuonen saddlepoint / Liu moment-matching fallback."
        )
        return False


# ---------------------------------------------------------------------------
# Tier 3: Liu moment-matching (Get_Liu_Params_Mod_Lambda)
# ---------------------------------------------------------------------------


def _liu_params(
    lambdas: np.ndarray,
) -> tuple[float, float, float, float, float, float]:
    """
    Compute Liu moment-matching parameters (R: ``Get_Liu_Params_Mod_Lambda``).

    Parameters
    ----------
    lambdas : np.ndarray
        Positive eigenvalues.

    Returns
    -------
    (mu_q, sigma_q, mu_x, sigma_x, ll, d)
        mu_q, sigma_q : mean/sd of the quadratic form Q
        mu_x, sigma_x : mean/sd of the matched non-central chi-squared
        ll : degrees of freedom
        d : non-centrality parameter
    """
    c1 = float(np.sum(lambdas))
    c2 = float(np.sum(lambdas**2))
    c3 = float(np.sum(lambdas**3))
    c4 = float(np.sum(lambdas**4))

    mu_q = c1
    sigma_q = np.sqrt(2.0 * c2)
    s1 = c3 / c2**1.5
    s2 = c4 / c2**2.0

    if s1**2 > s2:
        a = 1.0 / (s1 - np.sqrt(s1**2 - s2))
        d = s1 * a**3 - a**2
        ll = a**2 - 2.0 * d
    else:
        ll = 1.0 / s2
        a = np.sqrt(ll)
        d = 0.0

    mu_x = ll + d
    sigma_x = np.sqrt(2.0) * a

    return mu_q, sigma_q, mu_x, sigma_x, ll, d


def _liu_pvalue(q: float, lambdas: np.ndarray) -> float:
    """
    Compute P[Q > q] by Liu moment-matching (Liu et al. 2009).

    Matches the R SKAT source ``Get_Liu_Params_Mod_Lambda`` exactly.
    Uses scipy.stats.ncx2.sf for the non-central chi-squared CDF.

    Parameters
    ----------
    q : float
        Observed test statistic.
    lambdas : np.ndarray
        Array of positive eigenvalues (shape: (k,)).

    Returns
    -------
    float
        Estimated p-value in [0, 1].
    """
    if len(lambdas) == 0:
        return float("nan")

    # Cumulants: c[i] = sum(lambda^(i+1)) for i=0..3
    c = np.array(
        [
            np.sum(lambdas),  # c1
            np.sum(lambdas**2),  # c2
            np.sum(lambdas**3),  # c3
            np.sum(lambdas**4),  # c4
        ]
    )

    mu_q = c[0]
    sigma_q = np.sqrt(2.0 * c[1])
    s1 = c[2] / c[1] ** 1.5
    s2 = c[3] / c[1] ** 2.0

    if s1**2 > s2:
        # Noncentrality parameter case
        a = 1.0 / (s1 - np.sqrt(s1**2 - s2))
        d = s1 * a**3 - a**2
        ll = a**2 - 2.0 * d
    else:
        # All-equal eigenvalues: degenerates to central chi-squared
        ll = 1.0 / s2
        a = np.sqrt(ll)
        d = 0.0

    mu_x = ll + d
    sigma_x = np.sqrt(2.0) * a

    # Standardise Q to match scaled non-central chi-squared distribution
    q_norm = (q - mu_q) / sigma_q * sigma_x + mu_x

    # Guard: if q_norm <= 0 the tail probability is essentially 1
    if q_norm <= 0.0:
        return 1.0

    p = float(ncx2.sf(q_norm, df=ll, nc=d))
    # Clamp to [0, 1] for numerical stability
    return float(np.clip(p, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Tier 2: Kuonen Lugannani-Rice saddlepoint (Kuonen 1999)
# ---------------------------------------------------------------------------


def _kuonen_pvalue(q: float, lambdas: np.ndarray) -> float | None:
    """
    Compute P[Q > q] via Lugannani-Rice saddlepoint approximation.

    Uses scipy.optimize.brentq to solve K'(t_hat) = q where
    K(t) = -1/2 * sum log(1 - 2*t*lambda_j) is the CGF.

    Returns None if:
    - q <= mean(lambdas) (below mean; saddlepoint not valid here)
    - brentq fails to converge (ValueError)
    - the computed p-value is outside (0, 1)

    Parameters
    ----------
    q : float
        Observed test statistic.
    lambdas : np.ndarray
        Array of eigenvalues (positive entries only used for K).

    Returns
    -------
    float or None
        Saddlepoint p-value, or None if computation fails.
    """
    lam_pos = lambdas[lambdas > 0]
    if len(lam_pos) == 0:
        return None
    # Scale-free: the brentq bracket/xtol below are absolute
    q, lambdas, lam_pos = q / lam_pos.max(), lambdas / lam_pos.max(), lam_pos / lam_pos.max()

    t_max = 0.4999 / lam_pos.max()  # slightly inside singularity at 1/(2*lam_max)

    # CGF and derivatives for Q = sum(lambda_j * chi2_1)
    def cgf(t: float) -> float:
        """Cumulant generating function K(t)."""
        return float(-0.5 * np.sum(np.log(1.0 - 2.0 * t * lambdas)))

    def cgf_prime(t: float) -> float:
        """First derivative K'(t) = E[Q] at t."""
        return float(np.sum(lambdas / (1.0 - 2.0 * t * lambdas)))

    def cgf_dprime(t: float) -> float:
        """Second derivative K''(t)."""
        return float(2.0 * np.sum(lambdas**2 / (1.0 - 2.0 * t * lambdas) ** 2))

    # Mean of Q: K'(0) = sum(lambdas)
    mu = cgf_prime(0.0)
    if q <= mu:
        # Saddlepoint requires q > mean for upper-tail computation
        return None

    try:
        t_hat = brentq(lambda t: cgf_prime(t) - q, 1e-10, t_max, xtol=1e-12)
    except ValueError:
        return None

    w = float(np.sign(t_hat) * np.sqrt(2.0 * (t_hat * q - cgf(t_hat))))
    u = float(t_hat * np.sqrt(cgf_dprime(t_hat)))

    # Near-singularity guard: when u ~ 0 use simple normal approximation (SIM108)
    p = float(norm.sf(w)) if abs(u) < 1e-8 else float(norm.sf(w + np.log(u / w) / w))

    return p if 0.0 < p < 1.0 else None


# ---------------------------------------------------------------------------
# Tier 1: Davies via CFFI
# ---------------------------------------------------------------------------


def davies_pvalue(
    q: float,
    lambdas: np.ndarray,
    acc: float = 1e-9,
    lim: int = 1_000_000,
) -> tuple[float | None, int]:
    """
    Compute P[Q > q] using the Davies qfc() C++ extension via CFFI.

    Parameters
    ----------
    q : float
        Observed test statistic.
    lambdas : np.ndarray
        Array of positive eigenvalues.
    acc : float
        Maximum absolute error for the Davies integration. Default: 1e-9.
    lim : int
        Maximum number of integration terms. Default: 1_000_000.

    Returns
    -------
    (p_value, ifault)
        p_value is None if the extension is unavailable.
        ifault: 0=success, 1=accuracy not achieved, 2=round-off,
                3=invalid params, 4=integration failed, 5=out of memory.
    """
    if not _try_load_davies():
        return None, -1

    assert _qfc_lib is not None and _qfc_ffi is not None

    ffi = _qfc_ffi
    lib = _qfc_lib

    n = len(lambdas)
    lb1 = ffi.new("double[]", list(lambdas))
    nc1 = ffi.new("double[]", [0.0] * n)  # noncentrality params: all zero
    n1 = ffi.new("int[]", [1] * n)  # degrees of freedom: all 1
    r1 = ffi.new("int[1]", [n])
    sigma = ffi.new("double[1]", [0.0])
    c1 = ffi.new("double[1]", [q])
    lim1 = ffi.new("int[1]", [lim])
    acc_arr = ffi.new("double[1]", [acc])
    trace = ffi.new("double[7]", [0.0] * 7)
    ifault = ffi.new("int[1]", [0])
    res = ffi.new("double[1]", [0.0])

    with _QFC_LOCK:
        lib.qfc(lb1, nc1, n1, r1, sigma, c1, lim1, acc_arr, trace, ifault, res)

    # qfc returns the CDF P[Q <= q]; convert to survival P[Q > q] = 1 - CDF
    # (matches R CompQuadForm::davies which returns Qq = 1 - qfc_result)
    p_value = 1.0 - float(res[0])
    fault = int(ifault[0])
    return p_value, fault


def _tilted_sf_or_none(z: float, lambdas: np.ndarray) -> float | None:
    """P[sum(lambdas * chi2_1) > z] with relative accuracy deep in the upper tail.

    Davies' error is absolute (>= ~1e-9), useless for tail probabilities below that.
    Tilting by exp(t*Q) at the saddlepoint K'(t) = z gives Q' = sum(lambdas' * chi2_1),
    lambdas' = lambdas / (1 - 2 t lambdas), with mean z, and
    P[Q > z] = M(t) e^{-tz} E'[e^{-t(Q'-z)} 1{Q' > z}]
             = M(t) e^{-tz} int_0^inf e^{-s} P'[z < Q' <= z + s/t] ds,
    where Davies only evaluates O(1) probabilities of Q' (Gauss-Laguerre in s).
    Returns None if the tilt or a Davies evaluation fails.
    """
    # ponytail: validated to ~1e-7 relative against exact survival for k <= 3 (down to
    # 1e-16) and to Monte Carlo standard error (<= 1%) for k = 8, 11 (down to 1e-34);
    # k > 3 has no exact reference here. Cost ~17 Davies calls per tail node.
    scale = float(lambdas.max())  # scale-free: brentq's xtol is absolute
    z, lambdas = z / scale, lambdas / scale
    try:
        t = brentq(lambda u: float(np.sum(lambdas / (1 - 2 * u * lambdas))) - z, 0.0, 0.5 - 1e-12)
    except ValueError:  # no sign change: z beyond the representable tilt (or below the mean)
        return None
    if not 0.0 < t < 0.5:
        return None
    lam_t = lambdas / (1 - 2 * t * lambdas)
    log_m = -0.5 * float(np.sum(np.log1p(-2 * t * lambdas)))
    sf = []
    for q in np.r_[z, z + _LAG_NODES / t]:
        p, ifault = davies_pvalue(float(q), lam_t, acc=1e-6, lim=_DAVIES_RETRY_LIM)
        if p is None or ifault != 0 or not np.isfinite(p):
            return None
        sf.append(p)
    out = float(np.exp(log_m - t * z) * np.dot(_LAG_WEIGHTS, sf[0] - np.array(sf[1:])))
    return out if 0.0 <= out <= 1.0 else None


# ---------------------------------------------------------------------------
# Public API: compute_pvalue — three-tier fallback chain
# ---------------------------------------------------------------------------


def compute_pvalue(
    q: float,
    lambdas: np.ndarray,
    acc: float = 1e-6,
    lim: int = 10_000,
) -> tuple[float, str, bool]:
    """
    Compute P[Q > q] using R SKAT-compatible Davies -> Liu fallback chain.

    Matches R SKAT package ``Get_PValue.Lambda`` behavior with enhanced fallback:
    - Davies called with acc=1e-6, lim=10000 (R defaults); on ifault=1 (term limit)
      retried with lim=_DAVIES_RETRY_LIM, which is identical wherever lim=10000 converges
    - Relative accuracy budget 1e-4: where acc is too coarse for p, Davies is refined to
      acc = 1e-4 * p; below p = 1e-5 (acc floor 1e-9) Davies on the exponentially tilted
      mixture is used ("davies_tilted"); if Davies fails, saddlepoint, then Liu
    - Single eigenvalue: always use Liu
    - Three-tier chain when Davies available: Davies -> saddlepoint -> Liu

    Parameters
    ----------
    q : float
        Observed test statistic.
    lambdas : np.ndarray
        Filtered positive eigenvalues of the score kernel matrix.
    acc : float
        Davies integration accuracy. Default: 1e-6 (matches R SKAT).
    lim : int
        Davies integration term limit. Default: 10_000 (matches R SKAT).

    Returns
    -------
    (p_value, p_method, p_converged)
        p_value : float — estimated tail probability (NaN for empty lambdas)
        p_method : str  — "davies", "davies_tilted", "saddlepoint", or "liu"
        p_converged : bool — True only when Davies succeeded (ifault=0)
    """
    # Edge case: empty eigenvalue array
    if len(lambdas) == 0:
        return float("nan"), "liu", False

    # Edge case: non-positive test statistic — tail probability is ~1
    if q <= 0.0:
        return float(_liu_pvalue(q, lambdas)), "liu", False

    # R SKAT: single eigenvalue always uses Liu
    if len(lambdas) == 1:
        return float(_liu_pvalue(q, lambdas)), "liu", False

    # Always compute Liu as fallback (matches R SKAT Get_PValue.Lambda)
    p_liu = float(_liu_pvalue(q, lambdas))

    # --- Tier 1: Davies (requires compiled _qfc extension) ---
    if _try_load_davies():

        def davies(a: float) -> tuple[float | None, int]:
            p, f = davies_pvalue(q, lambdas, acc=a, lim=lim)
            if f == 1 and lim < _DAVIES_RETRY_LIM:
                p, f = davies_pvalue(q, lambdas, acc=a, lim=_DAVIES_RETRY_LIM)
            return p, f

        p_davies, ifault = davies(acc)
        if ifault != 0 and acc < 1e-6:  # requested accuracy unreachable: retry at R's
            acc = 1e-6
            p_davies, ifault = davies(acc)
        p_val = float("nan") if p_davies is None else float(p_davies)
        ok = ifault == 0 and 0.0 < p_val < 1.0
        if ok and p_val * _DAVIES_REL_ACC < acc:
            # absolute accuracy too coarse for p: refine where Davies can, else tail regime
            ok = p_val * _DAVIES_REL_ACC >= _DAVIES_MIN_ACC
            if ok:
                p2, f2 = davies(p_val * _DAVIES_REL_ACC)
                p_val = float(p2) if f2 == 0 and p2 is not None and 0.0 < p2 < 1.0 else p_val
        if ok:
            return p_val, "davies", True

        p_tail = p_val if ifault == 0 and np.isfinite(p_val) else p_liu
        if p_tail < _TAIL_P and q > float(np.sum(lambdas)):
            # Exact identity, Davies on O(1) tilted probabilities: relative accuracy in the
            # tail where plain Davies (absolute error) and the saddlepoint (~5%) are not
            p_t = _tilted_sf_or_none(q, lambdas)
            if p_t is not None and 0.0 < p_t < 1.0:
                return p_t, "davies_tilted", True

        p_sp = _kuonen_pvalue(q, lambdas)
        if p_sp is not None and 0.0 < p_sp < 1.0:
            logger.info(
                "Davies p=%s (ifault=%d); using saddlepoint fallback (p=%s)", p_val, ifault, p_sp
            )
            return float(p_sp), "saddlepoint", False
        logger.info(
            "Davies p=%s (ifault=%d); saddlepoint also failed; using Liu fallback (p=%s)",
            p_val,
            ifault,
            p_liu,
        )
        return p_liu, "liu", False

    # --- Davies unavailable: use saddlepoint -> Liu chain ---
    p_sp = _kuonen_pvalue(q, lambdas)
    if p_sp is not None and 0.0 < p_sp < 1.0:
        return float(p_sp), "saddlepoint", False

    if p_sp is None:
        logger.info(
            "Saddlepoint not applicable (q=%.6g <= mean(lambdas)=%.6g or brentq failed); "
            "using Liu moment-matching.",
            q,
            float(np.sum(lambdas)),
        )
    else:
        logger.info(
            "Saddlepoint returned out-of-range p=%s; using Liu moment-matching.",
            p_sp,
        )

    # --- Tier 3: Liu moment-matching (last resort, always returns a value) ---
    return p_liu, "liu", False
