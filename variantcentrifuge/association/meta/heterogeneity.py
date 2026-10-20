# File: variantcentrifuge/association/meta/heterogeneity.py
"""
Heterogeneity diagnostics for cohort meta-analysis.

- Per-stratum burden effect and standard error from the score/covariance.
- Cochran's Q, its p-value, and I-squared across strata with non-zero information
  (df = usable strata - 1).
"""

from __future__ import annotations

import logging

import numpy as np
import scipy.stats

from variantcentrifuge.association.backends.null_model import _DECISION_ACC

logger = logging.getLogger("variantcentrifuge")


def per_stratum_burden(
    u_list: list[np.ndarray],
    v_list: list[np.ndarray],
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return per-stratum burden effect estimates and standard errors.

    ``beta_c = (wᵀU_c) / (wᵀV_c w)``, ``se_c = 1/sqrt(wᵀV_c w)``. Strata whose
    information is numerically zero get ``beta=nan``, ``se=inf`` (same rule as
    ``combine.meta_burden``: below _DECISION_ACC of its Cauchy-Schwarz scale).
    """
    w = np.asarray(weights, dtype=np.float64)
    betas: list[float] = []
    ses: list[float] = []
    for u, v in zip(u_list, v_list, strict=True):
        score = float(w @ np.asarray(u, dtype=np.float64))
        v = np.asarray(v, dtype=np.float64)
        info = float(w @ (v @ w))
        scale = float((np.abs(w) * np.sqrt(np.clip(np.diag(v), 0.0, None))).sum()) ** 2
        if info <= _DECISION_ACC * scale:
            betas.append(np.nan)
            ses.append(np.inf)
        else:
            betas.append(score / info)
            ses.append(1.0 / np.sqrt(info))
    return np.array(betas), np.array(ses)


def cochran_q(betas: np.ndarray, ses: np.ndarray) -> tuple[float, float, float]:
    """Cochran's Q heterogeneity test. Returns (Q, p_value, I²).

    Uses only strata with finite beta and positive, finite SE. df = usable - 1.
    """
    betas = np.asarray(betas, dtype=np.float64)
    ses = np.asarray(ses, dtype=np.float64)
    usable = np.isfinite(betas) & np.isfinite(ses) & (ses > 0)
    b = betas[usable]
    s = ses[usable]
    k = b.size
    if k < 2:
        return 0.0, 1.0, 0.0
    w = 1.0 / (s**2)
    beta_bar = float(np.sum(w * b) / np.sum(w))
    q = float(np.sum(w * (b - beta_bar) ** 2))
    df = k - 1
    p = float(scipy.stats.chi2.sf(q, df))
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    return q, p, i2


__all__ = [
    "cochran_q",
    "per_stratum_burden",
]
