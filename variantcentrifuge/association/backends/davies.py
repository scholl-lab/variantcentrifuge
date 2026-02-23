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

Fallback triggering (GENESIS+GMMAT hybrid):
    - Davies ifault > 0
    - Davies p >= 1.0
    - Davies p < 1000 * machine_epsilon (~2.22e-13)
    - Davies p <= 1e-5  (proactive saddlepoint to avoid false convergence)

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

import numpy as np
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

_EPS1000 = 1000.0 * np.finfo(float).eps  # ~2.22e-13
_PROACTIVE_THRESHOLD = 1e-5  # GMMAT-style proactive saddlepoint trigger


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

    lib.qfc(lb1, nc1, n1, r1, sigma, c1, lim1, acc_arr, trace, ifault, res)

    # qfc returns the CDF P[Q <= q]; convert to survival P[Q > q] = 1 - CDF
    # (matches R CompQuadForm::davies which returns Qq = 1 - qfc_result)
    p_value = 1.0 - float(res[0])
    fault = int(ifault[0])
    return p_value, fault


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
    - Davies called with acc=1e-6, lim=10000 (R defaults)
    - If Davies ifault != 0 but 0 < p <= 1: keep Davies result (non-converged)
    - If p > 1 or p <= 0: try saddlepoint first, then Liu moment-matching
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
        p_method : str  — "davies", "saddlepoint", or "liu"
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
        p_davies, ifault = davies_pvalue(q, lambdas, acc=acc, lim=lim)

        if p_davies is not None:
            # R SKAT logic: keep Davies result if in valid range,
            # even if ifault != 0 (mark as non-converged)
            is_converged = ifault == 0
            p_val = float(p_davies)

            # Davies out of range: try saddlepoint first, then Liu
            if p_val > 1.0 or p_val <= 0.0:
                p_sp = _kuonen_pvalue(q, lambdas)
                if p_sp is not None and 0.0 < p_sp < 1.0:
                    logger.info(
                        "Davies p=%s out of range (ifault=%d); using saddlepoint fallback (p=%s)",
                        p_val,
                        ifault,
                        p_sp,
                    )
                    return float(p_sp), "saddlepoint", False
                logger.info(
                    "Davies p=%s out of range (ifault=%d); "
                    "saddlepoint also failed; using Liu fallback (p=%s)",
                    p_val,
                    ifault,
                    p_liu,
                )
                return p_liu, "liu", False

            if not is_converged:
                logger.debug(
                    "Davies ifault=%d but p=%s in range; keeping (non-converged)",
                    ifault,
                    p_val,
                )

            return p_val, "davies", is_converged

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
