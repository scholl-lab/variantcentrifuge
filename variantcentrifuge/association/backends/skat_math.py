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
    _GL_NODES_RAW,
    _GL_WEIGHTS_RAW,
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
