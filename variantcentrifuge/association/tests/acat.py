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
  For p < 1e-16, the standard tan((0.5 - p) * pi) formula overflows because
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
from typing import Sequence

import numpy as np
from scipy.stats import cauchy

logger = logging.getLogger("variantcentrifuge")

# Threshold below which the approximation 1/(p*pi) is used instead of tan((0.5-p)*pi).
# This prevents floating-point overflow for extreme p-values (Liu & Xie 2020).
_TINY_P_THRESHOLD = 1e-16


def cauchy_combination(
    p_values: np.ndarray | Sequence[float | None],
    weights: np.ndarray | Sequence[float] | None = None,
) -> float | None:
    """
    Combine p-values via the Cauchy combination formula (Liu & Xie 2020).

    Parameters
    ----------
    p_values : array-like of float | None
        Raw p-values to combine. None and NaN values are filtered out.
    weights : array-like of float | None, optional
        Non-negative weights for each p-value. Must be the same length as
        p_values (including None entries). Default: equal weights (1.0 each).
        Weights corresponding to filtered (None/NaN/>=1.0) p-values are
        dropped and remaining weights are re-normalized to sum to 1.

    Returns
    -------
    float | None
        Combined p-value in [0, 1], or None if fewer than 1 valid p-value
        exists after filtering. Single valid p-value is returned as-is
        (pass-through, per CONTEXT.md decision: k=1 is informative enough
        to surface without modification).

    Notes
    -----
    The Cauchy statistic is:
        T = sum(w_i * tan((0.5 - p_i) * pi))
    where w_i are normalized weights (sum to 1). The combined p-value is:
        p_combined = 1 - CDF_Cauchy(T) = SF_Cauchy(T)

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

    # Mask 1: filter None/NaN p-values
    valid_mask = ~np.isnan(p_arr)

    # Mask 2: filter non-informative p-values (p >= 1.0)
    # p=1.0 contributes tan((0.5 - 1.0)*pi) = tan(-pi/2) = -inf which is numerically problematic
    informative_mask = valid_mask & (p_arr < 1.0)

    p_valid = p_arr[informative_mask]
    w_valid = w_arr[informative_mask]

    n_valid = len(p_valid)

    if n_valid == 0:
        return None

    if n_valid == 1:
        # Pass-through: single valid p-value returned as-is (CONTEXT.md decision)
        logger.debug("cauchy_combination: single valid p-value, returning pass-through")
        return float(p_valid[0])

    # Normalize weights to sum to 1
    w_sum = w_valid.sum()
    if w_sum <= 0:
        # All-zero weights: fall back to equal weighting
        w_valid = np.ones(n_valid, dtype=float) / n_valid
    else:
        w_valid = w_valid / w_sum

    # Compute Cauchy transforms with numerical stability guard
    # For p < _TINY_P_THRESHOLD: use approximation 1/(p*pi) (Liu & Xie 2020)
    # For p >= _TINY_P_THRESHOLD: use exact tan((0.5 - p) * pi)
    tiny_mask = p_valid < _TINY_P_THRESHOLD
    transforms = np.empty(n_valid, dtype=float)
    if tiny_mask.any():
        transforms[tiny_mask] = 1.0 / (p_valid[tiny_mask] * np.pi)
    if (~tiny_mask).any():
        transforms[~tiny_mask] = np.tan((0.5 - p_valid[~tiny_mask]) * np.pi)

    # Weighted Cauchy statistic
    T = float(np.dot(w_valid, transforms))

    # Combined p-value: survival function of standard Cauchy distribution
    p_combined = float(np.clip(cauchy.sf(T), 0.0, 1.0))
    return p_combined


def compute_acat_o(test_p_values: dict[str, float | None]) -> float | None:
    """
    Compute ACAT-O (omnibus) p-value by combining test p-values per gene.

    This is the gene-level omnibus test that combines burden test p-values
    (e.g. Fisher, logistic burden) and SKAT p-values via equal-weight Cauchy
    combination. ACAT-O is the primary output for significance ranking.

    Parameters
    ----------
    test_p_values : dict mapping test_name -> p_value
        Dictionary of {test_name: p_value | None} for a single gene.
        Example: {"burden": 0.01, "skat": 0.05, "fisher": None}

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
