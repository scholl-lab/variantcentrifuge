# File: variantcentrifuge/association/meta/heterogeneity.py
"""
Heterogeneity diagnostics for cohort meta-analysis.

- Per-stratum burden effect and standard error from the score/covariance.
- Cochran's Q, its p-value, and I-squared across strata with non-zero information
  (df = usable strata - 1).
- Leave-one-stratum-out re-combination.
- An explicit, opt-in ACAT companion combining per-stratum p-values — offered as a
  heterogeneity-robust complement, never a silent substitute for the common-effect
  score test.
"""

from __future__ import annotations

import logging
from collections.abc import Callable

import numpy as np
import scipy.stats

logger = logging.getLogger("variantcentrifuge")


def per_stratum_burden(
    u_list: list[np.ndarray],
    v_list: list[np.ndarray],
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return per-stratum burden effect estimates and standard errors.

    ``beta_c = (wᵀU_c) / (wᵀV_c w)``, ``se_c = 1/sqrt(wᵀV_c w)``. Strata with
    non-positive information get ``beta=nan``, ``se=inf``.
    """
    w = np.asarray(weights, dtype=np.float64)
    betas: list[float] = []
    ses: list[float] = []
    for u, v in zip(u_list, v_list, strict=True):
        score = float(w @ np.asarray(u, dtype=np.float64))
        info = float(w @ (np.asarray(v, dtype=np.float64) @ w))
        if info <= 0.0:
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


def leave_one_out(
    u_list: list[np.ndarray],
    v_list: list[np.ndarray],
    weights: np.ndarray,
    stratum_ids: list[str],
    combine_fn: Callable[[list[np.ndarray], list[np.ndarray], np.ndarray], float | None],
) -> list[tuple[str, float | None]]:
    """Recompute the meta p-value dropping each stratum in turn.

    ``combine_fn`` maps (u_list, v_list, weights) -> p_value (e.g. a lambda over
    ``meta_burden``/``meta_skat``). Returns [(dropped_stratum_id, p_value), ...].
    """
    if len(u_list) < 2:
        return []
    results: list[tuple[str, float | None]] = []
    for i, sid in enumerate(stratum_ids):
        u_sub = [u for j, u in enumerate(u_list) if j != i]
        v_sub = [v for j, v in enumerate(v_list) if j != i]
        results.append((sid, combine_fn(u_sub, v_sub, weights)))
    return results


def acat_companion(p_values: list[float | None]) -> float | None:
    """Explicit ACAT (Cauchy) combination of per-stratum p-values (companion only)."""
    from variantcentrifuge.association.tests.acat import cauchy_combination

    return cauchy_combination(p_values)


__all__ = [
    "acat_companion",
    "cochran_q",
    "leave_one_out",
    "per_stratum_burden",
]
