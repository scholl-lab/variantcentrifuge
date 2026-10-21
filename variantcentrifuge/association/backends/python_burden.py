"""Weighted burden score test and the shared skip result of the Python SKAT backend."""

from __future__ import annotations

from typing import Any

import numpy as np
import scipy.stats

from variantcentrifuge.association.backends.base import NullModelResult
from variantcentrifuge.association.backends.null_model import (
    explained_columns,
    project_out_nuisance,
)


def skip_result(n_variants: int, n_marker_test: int, reason: str) -> dict[str, Any]:
    """Untestable-gene result (``p_value=None``) carrying ``reason`` as skip_reason."""
    return {
        "p_value": None,
        "rho": None,
        "n_variants": n_variants,
        "n_marker_test": n_marker_test,
        "warnings": [reason],
        "p_method": None,
        "p_converged": False,
        "skip_reason": reason,
    }


def burden_score_test(
    geno: np.ndarray, null_model: NullModelResult, weights: np.ndarray
) -> dict[str, Any]:
    """
    Burden score test for a single gene.

    Collapses variants (``geno``, n_samples x n_variants) to a weighted burden
    score and performs a score test against the null model residuals.
    """
    _n_samples, n_variants = geno.shape
    residuals: np.ndarray = null_model.extra["residuals"]
    sigma2: float = null_model.extra["sigma2"]

    # Collapse to burden score per sample: b = geno @ weights
    burden = geno @ weights  # shape: (n_samples,)

    # Score test: score = r^T b; Var(score) = sigma2 * b^T P b (null-model score covariance)
    score = float(residuals @ burden)
    burden_adj = project_out_nuisance(burden, null_model)
    variance = sigma2 * float(burden_adj @ burden_adj)

    # Shared three-way rule on the collapsed burden (as in LinearBurdenTest):
    # indeterminate -> skip; in the covariate span -> zero variance.
    burden_status = explained_columns(burden[:, np.newaxis], null_model)
    if burden_status is None:
        return skip_result(n_variants, n_variants, "ill_conditioned_covariates")
    if burden_status[0]:
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
