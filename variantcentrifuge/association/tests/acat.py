# File: variantcentrifuge/association/tests/acat.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/acat.py
"""
ACAT Cauchy combination formula for omnibus rare variant association testing.

Implements the Aggregated Cauchy Association Test (ACAT) as described in:
  Liu et al. (2019). "Acat: A fast and powerful p value combination method
  for rare-variant analysis in sequencing studies." AJHG 104(3):410-421.
  https://doi.org/10.1016/j.ajhg.2019.01.002

  Liu and Xie (2020). "Cauchy combination test: a powerful test with analytic
  p-value calculation under arbitrary dependency structures." JASA 115:393-402.
  https://doi.org/10.1080/01621459.2018.1554485

Numerical stability note (Liu & Xie 2020, Section 2.2):
  For p < 1e-15 (the R ACAT() threshold), the standard tan((0.5 - p) * pi) formula overflows because
  tan approaches +inf. The approximation tan((0.5 - p) * pi) ≈ 1 / (p * pi)
  is used instead, which is accurate to machine precision for tiny p.

Examples
--------
>>> import numpy as np
>>> cauchy_combination(np.array([0.05, 0.05, 0.05]))  # approx 0.05
>>> cauchy_combination(np.array([0.5, 0.5]))           # approx 0.5
>>> cauchy_combination(np.array([0.001, 0.01, 0.05]))  # approx 0.00268
"""

from __future__ import annotations

import logging
from collections.abc import Sequence

import numpy as np
from scipy.stats import norm

from variantcentrifuge.association.backends.base import NullModelResult
from variantcentrifuge.association.backends.null_model import (
    _project,
    classify_projected,
    full_rank_design,
    project_out_nuisance,
)

logger = logging.getLogger("variantcentrifuge")

# Threshold below which the approximation 1/(p*pi) is used instead of tan((0.5-p)*pi).
# This prevents floating-point overflow for extreme p-values (Liu & Xie 2020).
_TINY_P_THRESHOLD = 1e-15  # as R ACAT()


def cauchy_combination(
    p_values: np.ndarray | Sequence[float | None],
    weights: np.ndarray | Sequence[float] | None = None,
) -> float | None:
    """
    Combine p-values via the Cauchy combination formula (Liu & Xie 2020).

    Parameters
    ----------
    p_values : array-like of float | None
        Raw p-values to combine. None and NaN values are filtered out; the rest
        follow R ``ACAT()`` (p = 1 is kept; p outside [0, 1] raises ValueError).
    weights : array-like of float | None, optional
        Non-negative weights for each p-value. Must be the same length as
        p_values (including None entries). Default: equal weights (1.0 each).
        Weights corresponding to filtered (None/NaN) p-values are
        dropped and remaining weights are re-normalized to sum to 1.

    Returns
    -------
    float | None
        Combined p-value in [0, 1], or None if fewer than 1 valid p-value
        exists after filtering or the p-values contain both 0 and 1 (undefined).
        Single valid p-value is returned as-is
        (pass-through, per CONTEXT.md decision: k=1 is informative enough
        to surface without modification).

    Notes
    -----
    The Cauchy statistic is:
        t_stat = sum(w_i * tan((0.5 - p_i) * pi))
    where w_i are normalized weights (sum to 1). The combined p-value is:
        p_combined = 1 - CDF_Cauchy(t_stat) = SF_Cauchy(t_stat)

    For p_i < _TINY_P_THRESHOLD, the approximation tan((0.5 - p) * pi) ≈ 1/(p*pi)
    is used (numerically equivalent but avoids overflow).
    """
    # Convert to numpy array for vectorized operations, handling None/NaN
    p_arr = np.array(
        [np.nan if (v is None or (isinstance(v, float) and np.isnan(v))) else v for v in p_values],
        dtype=float,
    )

    # Build corresponding weights array
    if weights is not None:
        w_arr = np.array(weights, dtype=float)
        if len(w_arr) != len(p_arr):
            raise ValueError(
                f"weights length ({len(w_arr)}) must match p_values length ({len(p_arr)})"
            )
    else:
        w_arr = np.ones(len(p_arr), dtype=float)

    # None/NaN mark tests that did not run; R ACAT() itself never sees them
    valid_mask = ~np.isnan(p_arr)
    p_valid = p_arr[valid_mask]
    w_valid = w_arr[valid_mask]
    n_valid = len(p_valid)

    if n_valid == 0:
        return None
    # Validation as in R ACAT(): p outside [0, 1] is an error, p = 0 with p = 1 is undefined
    if ((p_valid < 0.0) | (p_valid > 1.0)).any():
        raise ValueError(f"P-values must be between 0 and 1, got {p_valid.tolist()}")
    if (p_valid == 0.0).any() and (p_valid == 1.0).any():
        logger.debug("cauchy_combination: p-values of exactly 0 and 1, ACAT undefined")
        return None

    if n_valid == 1:
        # Pass-through: single valid p-value returned as-is (CONTEXT.md decision)
        logger.debug("cauchy_combination: single valid p-value, returning pass-through")
        return float(p_valid[0])

    # Normalize weights to sum to 1
    # All-zero weights: fall back to equal weighting
    w_sum = w_valid.sum()
    w_valid = np.ones(n_valid, dtype=float) / n_valid if w_sum <= 0 else w_valid / w_sum

    # As R ACAT(): 1/(p*pi) below the threshold (p = 0 gives +inf, combined p 0) and
    # tan((0.5 - p) * pi) otherwise (p = 1 gives about -1.6e16, kept like any other p)
    tiny_mask = p_valid < _TINY_P_THRESHOLD
    with np.errstate(divide="ignore"):
        transforms = np.where(tiny_mask, 1.0 / p_valid / np.pi, np.tan((0.5 - p_valid) * np.pi))
    t_stat = float(np.dot(w_valid, transforms))
    if np.isnan(t_stat):  # p = 0 with weight 0: R returns NaN
        return None

    # pcauchy(t, lower.tail = FALSE) without the 0.5 - atan(t)/pi cancellation at large t
    return float(np.arctan2(1.0, t_stat) / np.pi)


def compute_acat_v(
    geno: np.ndarray,
    residuals: np.ndarray,
    trait_type: str,
    sigma2: float,
    mu_hat: np.ndarray | None,
    weights: np.ndarray,
    design_matrix: np.ndarray | None = None,
    null_model: NullModelResult | None = None,
    projection: tuple[np.ndarray, np.ndarray, float] | None = None,
) -> float | None:
    """
    ACAT-V per-variant score test (Liu & Xie 2020, STAARpipeline).

    Computes marginal score test p-values for each variant and combines
    them via Cauchy combination with Beta(MAF) weights. Powerful when
    only a small fraction of variants are causal.

    Parameters
    ----------
    geno : np.ndarray, shape (n_samples, n_variants)
        Genotype dosage matrix (0/1/2). No NaN.
    residuals : np.ndarray, shape (n_samples,)
        Null model residuals (y - mu_hat).
    trait_type : str
        "binary" or "quantitative".
    sigma2 : float
        Variance estimate from null model.
    mu_hat : np.ndarray | None, shape (n_samples,)
        Fitted values from null model. Required for binary traits;
        may be None for quantitative traits (sigma2 is used instead).
    weights : np.ndarray, shape (n_variants,)
        Per-variant weights (typically Beta(MAF; 1, 25)).
    design_matrix : np.ndarray | None, shape (n_samples, k)
        Fitted design matrix (covariates + intercept) from null model; it is
        reduced to a full-rank basis that always contains the intercept.
        When provided, null model covariates are projected out to obtain
        unbiased score test variances (Liu & Xie 2020).
    null_model : NullModelResult | None
        Fitted null model; preferred over ``design_matrix`` because it carries the
        design conditioning used to recognise variants explained by covariates.
    projection : tuple | None
        ``null_model._project`` of the tested (nonzero-weight) columns of ``geno``, if
        the caller already has it (the SKAT backend's ``result["projection"]``).

    Returns
    -------
    float | None
        Combined ACAT-V p-value, or None if no valid per-variant p-values.
    """
    _n_samples, n_variants = geno.shape
    if n_variants == 0:
        return None

    # Per-variant score statistics: S_j = G_j' @ residuals
    score_vec = geno.T @ residuals  # shape: (n_variants,)

    # Per-variant variance: V_j = G_j' @ V @ G_j
    # With null model covariates X projected out (Liu 2020, Get.marginal.pval):
    # Binary: V = diag(G' [W - W X (X' W X)^-1 X' W] G)
    # Quantitative: V = sigma2 * diag(G' [I - X (X' X)^-1 X'] G)
    if null_model is None and design_matrix is not None:
        # Raw design (intercept + covariates): reduce to the shared full-rank basis
        # (duplicate/zero/constant columns dropped, intercept kept) before projecting.
        basis, covariate_r = full_rank_design(
            np.asarray(design_matrix, dtype=np.float64), len(residuals)
        )
        null_model = NullModelResult(
            model=None,
            trait_type="binary"
            if trait_type == "binary" and mu_hat is not None
            else "quantitative",
            n_samples=len(residuals),
            extra={"design_matrix": basis, "covariate_r": covariate_r, "mu_hat": mu_hat},
        )
    if null_model is not None:
        # Shared QR projection; variants in the covariate space get zero variance (skipped)
        # Classify only the variants that are tested (nonzero weight)
        tested = np.asarray(weights) != 0.0
        if projection is None:
            projection = _project(geno[:, tested], null_model)
        explained_tested = classify_projected(*projection)
        if explained_tested is None:
            return None  # ill-conditioned covariates: projections not trustworthy
        explained = ~tested
        explained[np.flatnonzero(tested)[explained_tested]] = True
        # a column subset takes other BLAS kernels (last bits differ): reuse only the full set
        z_adj = projection[0] if tested.all() else project_out_nuisance(geno, null_model)
        g_ss = np.where(explained, 0.0, np.sum(z_adj**2, axis=0))
        var_vec = g_ss if null_model.trait_type == "binary" else sigma2 * g_ss
    else:
        if trait_type == "binary" and mu_hat is not None:
            var_weights = mu_hat * (1.0 - mu_hat)  # (n_samples,)
            g_ss = np.sum(geno**2 * var_weights[:, np.newaxis], axis=0)  # (n_variants,)
            var_vec = g_ss
        else:
            g_ss = np.sum(geno**2, axis=0)  # (n_variants,)
            var_vec = sigma2 * g_ss

    # Single-pass: compute per-variant z-statistics, p-values, and collect
    # valid (p, weight) pairs in one loop
    valid_p: list[float] = []
    valid_w: list[float] = []
    for j in range(n_variants):
        if g_ss[j] <= 1e-12 or var_vec[j] <= 0.0:
            continue  # skip zero-variance or collinear variants
        z_j = score_vec[j] / np.sqrt(var_vec[j])
        # Every marginal p in [0, 1] goes to ACAT, as STAAR; R ACAT() rules handle 0 and 1
        valid_p.append(2.0 * float(norm.sf(abs(z_j))))
        valid_w.append(float(weights[j]))

    if not valid_p:
        return None

    return cauchy_combination(valid_p, valid_w)


def compute_acat_o(test_p_values: dict[str, float | None]) -> float | None:
    """
    Compute ACAT-O (omnibus) p-value by combining test p-values per gene.

    This is the gene-level omnibus test that combines burden, SKAT/SKAT-O and
    ACAT-V p-values via equal-weight Cauchy combination.
    ACAT-O is the primary output for significance ranking.

    Parameters
    ----------
    test_p_values : dict mapping test_name -> p_value
        Dictionary of {test_name: p_value | None} for a single gene.
        Example: {"burden": 0.01, "skat": 0.05, "acat_v": None}

    Returns
    -------
    float | None
        Combined ACAT-O p-value, or None if no valid p-values exist.
        Single valid p-value is returned as-is (pass-through).

    Examples
    --------
    >>> compute_acat_o({"burden": None, "skat": None})
    None
    >>> compute_acat_o({"burden": 0.05, "skat": None})
    0.05
    >>> compute_acat_o({"burden": 0.01, "skat": 0.05})  # < 0.05
    """
    p_list = list(test_p_values.values())
    return cauchy_combination(p_list)
