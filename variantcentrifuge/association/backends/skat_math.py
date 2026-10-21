# File: variantcentrifuge/association/backends/skat_math.py
# Location: variantcentrifuge/association/backends/skat_math.py
"""Mathematical routines and numerical integration for Python SKAT / SKAT-O."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scipy.linalg

# SKAT-O integration lives in skato_integrate; re-exported for existing imports.
from variantcentrifuge.association.backends.skato_integrate import (  # noqa: F401
    _GL_A,
    _GL_B,
    _SKATO_DAVIES_MIN_ACC,
    _SKATO_DAVIES_REL_ACC,
    _SKATO_X_MAX,
    _SKATO_Z_LEVELS,
    _DaviesFailureError,
    _get_lambda,
    _skato_breakpoints,
    _skato_each_q,
    _skato_get_pvalue,
    _skato_integrate_davies,
    _skato_integrate_liu,
    _skato_optimal_param,
)

if TYPE_CHECKING:
    from variantcentrifuge.association.backends.base import NullModelResult

# Fixed SKAT-O rho search grid (matches R SKAT package default)
_SKATO_RHO_GRID = [0.0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0]


def compute_eigenvalues_filtered(
    geno_weighted: np.ndarray,
    null_model: NullModelResult,
    projected: tuple[np.ndarray, np.ndarray, float] | None = None,
) -> np.ndarray:
    """
    Compute filtered eigenvalues of the null-projected SKAT kernel ``Z'PZ / 2``.

    ``P`` is the null-model score-covariance projector (see
    ``project_out_nuisance``): Bernoulli-variance weighted for binary traits,
    ``I - X(X'X)^-1 X'`` for quantitative traits.  For quantitative traits the
    matching statistic is ``Q = ||Z'r||^2 / (2 sigma2)`` (R ``SKAT.linear``),
    so the eigenvalues carry no ``sigma2`` and are invariant to outcome units.

    Uses scipy.linalg.eigh(driver='evr') matching R's DSYEVR Lapack routine.
    Applies threshold = mean(positive_eigenvalues) / 100_000 (matches R SKAT),
    floored at ``(bound * sigma_max(V^1/2 Z))**2 / 2`` (round-off of ``Z'PZ``; ``bound`` is
    the per-column round-off bound of ``null_model._project``).

    Parameters
    ----------
    geno_weighted : np.ndarray, shape (n_samples, n_variants)
        Weighted genotype matrix (G * W).
    null_model : NullModelResult
        Fitted null model containing mu_hat (binary) and design_matrix.
    projected : tuple, optional
        ``null_model._project(geno_weighted, null_model)``, if the caller already has it.

    Returns
    -------
    np.ndarray
        Filtered eigenvalues (positive and above threshold). Empty array
        if none pass.
    """
    from variantcentrifuge.association.backends.null_model import _project

    z_adj, z, bound = projected if projected is not None else _project(geno_weighted, null_model)
    n_variants = geno_weighted.shape[1]
    eig_mat = z_adj.T @ z_adj if n_variants <= geno_weighted.shape[0] else z_adj @ z_adj.T
    lambdas_all = scipy.linalg.eigh(eig_mat / 2.0, eigvals_only=True, driver="evr")

    pos = lambdas_all[lambdas_all >= 0]
    if len(pos) == 0:
        return np.array([], dtype=np.float64)
    # R's relative filter, floored at the round-off level of the unprojected kernel
    # (directions explained by covariates must not become noise eigenvalues).
    # Singular-value round-off floor: bound * sigma_max of the UNprojected z (callers zero
    # explained columns first, so they cannot inflate sigma_max).
    threshold = max(pos.mean() / 100_000.0, float(bound * np.linalg.norm(z, 2)) ** 2 / 2.0)
    return np.asarray(lambdas_all[lambdas_all > threshold], dtype=np.float64)


def compute_skato_q_and_z1(
    geno_weighted: np.ndarray,
    null_model: NullModelResult,
    rho_grid: list[float],
    z_adj: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-rho Q statistics and projection-adjusted Z1/sqrt(2) for SKAT-O.

    ``z_adj`` is ``project_out_nuisance(geno_weighted, null_model)`` if the caller has it.
    """
    from variantcentrifuge.association.backends.null_model import project_out_nuisance

    residuals: np.ndarray = null_model.extra["residuals"]
    # Binary sigma2 is 1; quantitative Q is scaled by the residual variance (R SKAT_Optimal_Linear)
    sigma2: float = null_model.extra["sigma2"]

    # ── Step 1: Compute Q_rho for each rho ──
    score_vec = geno_weighted.T @ residuals
    score_sum = float(np.sum(score_vec))

    q_rho = np.zeros(len(rho_grid), dtype=np.float64)
    for i, rho in enumerate(rho_grid):
        q1 = (1.0 - rho) * float(score_vec @ score_vec)
        q2 = rho * score_sum**2
        q_rho[i] = (q1 + q2) / (2.0 * sigma2)

    # ── Step 2: Build projection-adjusted Z1 ──
    if z_adj is None:
        z_adj = project_out_nuisance(geno_weighted, null_model)
    z1_half = z_adj / np.sqrt(2.0)
    return q_rho, z1_half
