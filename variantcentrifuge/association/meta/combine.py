# File: variantcentrifuge/association/meta/combine.py
"""
Common-effect meta-analysis of aligned per-stratum score/covariance statistics.

Given per-stratum unweighted score vectors ``U_c`` and projection-adjusted
covariances ``V_c`` (already aligned to a shared union of variants, zero-padded
for variants absent in a stratum) and a single harmonised weight vector ``w``:

    S = Σ_c diag(w) U_c            (weighted meta score)
    Σ = Σ_c diag(w) V_c diag(w)    (weighted meta covariance)

- meta-burden : Z = 1ᵀS / sqrt(1ᵀΣ1)
- meta-SKAT   : Q = ‖S‖²/2, null Q ~ Σ λ_i χ²₁, λ = eig(Σ/2)  (Davies chain)
- meta-SKAT-O : factor Σ/2 = BᵀB exactly, reuse the existing SKAT-O omnibus.

A single-stratum call reproduces ``PythonSKATBackend._test_skat`` /
``_test_skato`` / the projection-adjusted burden score test bit-for-bit (golden
gate), because these formulas are the exact weighted forms of what the backend
computes internally.
"""

from __future__ import annotations

import logging

import numpy as np
import scipy.linalg
import scipy.stats

from variantcentrifuge.association.backends.davies import compute_pvalue
from variantcentrifuge.association.backends.python_backend import (
    _SKATO_RHO_GRID,
    _skato_get_pvalue,
)

logger = logging.getLogger("variantcentrifuge")

# Matches R SKAT / _get_lambda eigenvalue-retention threshold.
_LAMBDA_REL_THRESHOLD = 100_000.0


def _weighted_moments(
    u_list: list[np.ndarray],
    v_list: list[np.ndarray],
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the weighted meta score S and weighted meta covariance Σ."""
    w = np.asarray(weights, dtype=np.float64)
    s = np.zeros_like(w)
    sigma = np.zeros((w.size, w.size), dtype=np.float64)
    for u, v in zip(u_list, v_list, strict=True):
        s += w * np.asarray(u, dtype=np.float64)
        sigma += (w[:, np.newaxis] * np.asarray(v, dtype=np.float64)) * w[np.newaxis, :]
    sigma = 0.5 * (sigma + sigma.T)
    return s, sigma


def _filter_lambdas(mat: np.ndarray) -> np.ndarray:
    """Eigenvalues of a symmetric matrix, kept per R SKAT (mean(pos)/1e5) threshold."""
    lam = scipy.linalg.eigh(mat, eigvals_only=True, driver="evr")
    pos = lam[lam >= 0]
    if pos.size == 0:
        return np.array([], dtype=np.float64)
    return np.asarray(lam[lam > pos.mean() / _LAMBDA_REL_THRESHOLD], dtype=np.float64)


def meta_skat(
    u_list: list[np.ndarray],
    v_list: list[np.ndarray],
    weights: np.ndarray,
) -> tuple[float | None, str]:
    """Meta-SKAT p-value. Returns (p_value, p_method); p_value None if untestable."""
    s, sigma = _weighted_moments(u_list, v_list, weights)
    q = float(s @ s) / 2.0
    lambdas = _filter_lambdas(sigma / 2.0)
    if lambdas.size == 0:
        return 1.0, "degenerate"
    p, method, _converged = compute_pvalue(q, lambdas)
    if not np.isfinite(p):
        return None, method
    return float(np.clip(p, 0.0, 1.0)), method


def meta_burden(
    u_list: list[np.ndarray],
    v_list: list[np.ndarray],
    weights: np.ndarray,
) -> tuple[float | None, float, int]:
    """Common-effect burden meta. Returns (p_value, beta, direction).

    ``beta`` is the score-based effect estimate S_b / Var_b; ``direction`` is its
    sign (+1 risk, -1 protective, 0 undefined).
    """
    s, sigma = _weighted_moments(u_list, v_list, weights)
    score = float(s.sum())
    var = float(sigma.sum())
    if var <= 0.0:
        return None, 0.0, 0
    z = score / np.sqrt(var)
    p = float(2.0 * scipy.stats.norm.sf(abs(z)))
    beta = score / var
    return p, beta, int(np.sign(score))


def meta_skato(
    u_list: list[np.ndarray],
    v_list: list[np.ndarray],
    weights: np.ndarray,
) -> tuple[float | None, float | None]:
    """Meta-SKAT-O omnibus p-value. Returns (p_value, optimal_rho)."""
    s, sigma = _weighted_moments(u_list, v_list, weights)
    a_mat = sigma / 2.0

    # Exact factor B with BᵀB = A: keep ALL non-negative eigenmodes (clip round-off).
    eigvals, eigvecs = scipy.linalg.eigh(a_mat)
    keep = eigvals > 0.0
    if not np.any(keep):
        return None, None
    b = (eigvecs[:, keep] * np.sqrt(eigvals[keep])).T  # (k, p): b.T @ b == a_mat

    score_sq = float(s @ s)
    score_sum_sq = float(s.sum()) ** 2
    q_rho = np.array(
        [((1.0 - rho) * score_sq + rho * score_sum_sq) / 2.0 for rho in _SKATO_RHO_GRID],
        dtype=np.float64,
    )

    p_value, per_rho = _skato_get_pvalue(q_rho, b, _SKATO_RHO_GRID)
    best_idx = int(np.argmin(per_rho))
    best_rho = _SKATO_RHO_GRID[best_idx]
    if best_rho >= 0.999:
        best_rho = 1.0
    if not np.isfinite(p_value):
        return None, best_rho
    return float(np.clip(p_value, 0.0, 1.0)), best_rho


__all__ = ["meta_burden", "meta_skat", "meta_skato"]
