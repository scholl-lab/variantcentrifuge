# File: variantcentrifuge/association/weights.py
# Location: variantcentrifuge/variantcentrifuge/association/weights.py
"""
Variant weight schemes for rare variant association burden tests.

Provides public functions:

- ``beta_maf_weights``: Beta(MAF; a, b) weights (default a=1, b=25) matching
  the SKAT R package convention — rare variants receive higher weights.
- ``uniform_weights``: All-ones weight vector for unweighted burden tests.
- ``cadd_weights``: Beta(MAF) x min(CADD_phred / cap, 1.0) weights with
  variant-type-aware fallback for missing CADD scores.
- ``revel_weights``: Beta(MAF) x REVEL_score weights with fallback for missing
  REVEL scores. REVEL is already in [0, 1]; no normalization needed.
- ``combined_weights``: Beta(MAF) x functional score. Uses CADD by default;
  falls back to REVEL if CADD is not provided.
- ``get_weights``: String-spec parser that dispatches to any of the above.

Weight spec string format
-------------------------
``"beta:a,b"``
    Beta distribution weights. ``a`` and ``b`` are float parameters.
    Example: ``"beta:1,25"`` (default, SKAT convention).

``"uniform"``
    Uniform weights; all variants receive weight 1.0.

``"cadd"``
    Beta(MAF) x min(CADD_phred / 40, 1.0). Requires ``cadd_scores`` kwarg.

``"revel"``
    Beta(MAF) x REVEL_score. Requires ``revel_scores`` kwarg.

``"combined"``
    Beta(MAF) x functional score. Uses CADD if provided; falls back to REVEL.

Fallback rules for missing annotation scores
--------------------------------------------
LoF variants (stop_gained, frameshift, splice site) missing CADD/REVEL receive
functional weight 1.0 (maximum — conservative, do not penalize).
Missense and other variants missing scores also receive functional weight 1.0
(Beta(MAF)-only — no up-weighting, no down-weighting).
A warning is logged listing per-category missing counts.
"""

from __future__ import annotations

import logging

import numpy as np
from scipy.stats import beta as _beta_dist

logger = logging.getLogger("variantcentrifuge")

# Variant effect constants for type-aware fallback
LOF_EFFECTS = frozenset(
    {
        "stop_gained",
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
    }
)
MISSENSE_EFFECTS = frozenset({"missense_variant"})


def beta_maf_weights(
    mafs: np.ndarray,
    a: float = 1.0,
    b: float = 25.0,
) -> np.ndarray:
    """
    Compute Beta(MAF; a, b) density weights for each variant.

    Follows the SKAT R package convention (default a=1, b=25): rare variants
    (low MAF) receive substantially higher weights than common variants.

    Parameters
    ----------
    mafs : np.ndarray, shape (n_variants,)
        Minor allele frequencies in [0, 1].
    a : float
        Beta distribution shape parameter alpha. Default: 1.0.
    b : float
        Beta distribution shape parameter beta. Default: 25.0.

    Returns
    -------
    np.ndarray, shape (n_variants,), float64
        Weight for each variant. Values are strictly positive.

    Notes
    -----
    MAF values are clipped to ``[1e-8, 1 - 1e-8]`` before evaluation to
    avoid numerical edge issues at 0 and 1.
    """
    maf_clipped = np.clip(np.asarray(mafs, dtype=np.float64), 1e-8, 1 - 1e-8)
    return _beta_dist.pdf(maf_clipped, a=a, b=b)


def uniform_weights(n_variants: int) -> np.ndarray:
    """
    Return a uniform (all-ones) weight vector.

    Parameters
    ----------
    n_variants : int
        Number of variants.

    Returns
    -------
    np.ndarray, shape (n_variants,), float64
        Vector of ones.
    """
    return np.ones(n_variants, dtype=np.float64)


def _parse_scores_to_float(scores: np.ndarray | None) -> np.ndarray | None:
    """
    Coerce an annotation score array to float64, converting None/'.' to NaN.

    Returns None if input is None.
    """
    if scores is None:
        return None
    arr = np.asarray(scores, dtype=object)
    result = np.empty(len(arr), dtype=np.float64)
    for i, v in enumerate(arr):
        if v is None or (isinstance(v, str) and v.strip() in (".", "", "NA", "nan")):
            result[i] = np.nan
        else:
            try:
                result[i] = float(v)
            except (ValueError, TypeError):
                result[i] = np.nan
    return result


def _log_missing_score_counts(
    nan_mask: np.ndarray,
    variant_effects: np.ndarray | None,
    scheme_name: str,
) -> None:
    """
    Log a warning summarising the number of variants with missing annotation scores
    broken down by variant type (LoF, missense, other).

    Parameters
    ----------
    nan_mask : np.ndarray of bool, shape (n_variants,)
        True where annotation score is missing.
    variant_effects : np.ndarray or None
        Variant effect strings (e.g. "stop_gained", "missense_variant").
    scheme_name : str
        Weight scheme name for log message (e.g. "CADD", "REVEL").
    """
    n_missing = int(nan_mask.sum())
    if n_missing == 0:
        return

    if variant_effects is None:
        logger.warning(
            f"{scheme_name} weights: {n_missing} variant(s) had missing scores "
            "(used Beta(MAF)-only fallback)"
        )
        return

    effects_arr = np.asarray(variant_effects, dtype=object)
    missing_effects = effects_arr[nan_mask]

    n_lof = int(sum(1 for e in missing_effects if e in LOF_EFFECTS))
    n_miss = int(sum(1 for e in missing_effects if e in MISSENSE_EFFECTS))
    n_other = n_missing - n_lof - n_miss

    logger.warning(
        f"{scheme_name} weights: {n_missing} variant(s) had missing scores — "
        f"{n_lof} LoF (functional=1.0), "
        f"{n_miss} missense (fallback), "
        f"{n_other} other (fallback)"
    )


def cadd_weights(
    mafs: np.ndarray,
    cadd_scores: np.ndarray,
    variant_effects: np.ndarray | None = None,
    cap: float = 40.0,
) -> np.ndarray:
    """
    Compute Beta(MAF) x min(CADD_phred / cap, 1.0) weights.

    CADD scores above ``cap`` are capped at 1.0. Variants with missing CADD
    scores (NaN, None, '.') receive functional weight 1.0 (Beta(MAF)-only —
    no up-weighting or down-weighting).

    Parameters
    ----------
    mafs : np.ndarray, shape (n_variants,)
        Minor allele frequencies in [0, 1].
    cadd_scores : np.ndarray, shape (n_variants,)
        CADD Phred scores. May contain NaN, None, or '.' for missing values.
    variant_effects : np.ndarray or None, shape (n_variants,)
        Variant effect strings for type-aware fallback logging. Optional.
    cap : float
        Normalization cap for CADD Phred. Default: 40.0.

    Returns
    -------
    np.ndarray, shape (n_variants,), float64
        Combined Beta(MAF) x functional weights.
    """
    maf_w = beta_maf_weights(np.asarray(mafs, dtype=np.float64))
    scores_f = _parse_scores_to_float(np.asarray(cadd_scores, dtype=object))
    nan_mask = np.isnan(scores_f)

    functional = np.where(nan_mask, 1.0, np.minimum(scores_f / cap, 1.0))

    _log_missing_score_counts(nan_mask, variant_effects, "CADD")

    return maf_w * functional


def revel_weights(
    mafs: np.ndarray,
    revel_scores: np.ndarray,
    variant_effects: np.ndarray | None = None,
) -> np.ndarray:
    """
    Compute Beta(MAF) x REVEL_score weights.

    REVEL scores are already in [0, 1]; no normalization is applied. Variants
    with missing REVEL scores (NaN, None, '.') receive functional weight 1.0
    (Beta(MAF)-only).

    Parameters
    ----------
    mafs : np.ndarray, shape (n_variants,)
        Minor allele frequencies in [0, 1].
    revel_scores : np.ndarray, shape (n_variants,)
        REVEL scores in [0, 1]. May contain NaN, None, or '.' for missing.
    variant_effects : np.ndarray or None, shape (n_variants,)
        Variant effect strings for type-aware fallback logging. Optional.

    Returns
    -------
    np.ndarray, shape (n_variants,), float64
        Combined Beta(MAF) x functional weights.
    """
    maf_w = beta_maf_weights(np.asarray(mafs, dtype=np.float64))
    scores_f = _parse_scores_to_float(np.asarray(revel_scores, dtype=object))
    nan_mask = np.isnan(scores_f)

    functional = np.where(nan_mask, 1.0, scores_f)

    _log_missing_score_counts(nan_mask, variant_effects, "REVEL")

    return maf_w * functional


def combined_weights(
    mafs: np.ndarray,
    cadd_scores: np.ndarray | None = None,
    revel_scores: np.ndarray | None = None,
    variant_effects: np.ndarray | None = None,
    cadd_cap: float = 40.0,
) -> np.ndarray:
    """
    Compute Beta(MAF) x functional score weights using the best available score.

    If both CADD and REVEL scores are provided, CADD is preferred. If only
    REVEL is provided, uses REVEL. If neither is provided, falls back to plain
    Beta(MAF) weights.

    Parameters
    ----------
    mafs : np.ndarray, shape (n_variants,)
        Minor allele frequencies in [0, 1].
    cadd_scores : np.ndarray or None
        CADD Phred scores. Preferred over REVEL when both are provided.
    revel_scores : np.ndarray or None
        REVEL scores in [0, 1]. Used when CADD scores are not provided.
    variant_effects : np.ndarray or None
        Variant effect strings for type-aware fallback logging.
    cadd_cap : float
        Normalization cap for CADD Phred. Default: 40.0.

    Returns
    -------
    np.ndarray, shape (n_variants,), float64
        Combined Beta(MAF) x functional weights.
    """
    if cadd_scores is not None:
        return cadd_weights(mafs, cadd_scores, variant_effects, cap=cadd_cap)
    if revel_scores is not None:
        return revel_weights(mafs, revel_scores, variant_effects)
    # No annotation scores: return plain Beta(MAF) weights
    logger.warning(
        "combined_weights: no CADD or REVEL scores provided; falling back to Beta(MAF)-only weights"
    )
    return beta_maf_weights(np.asarray(mafs, dtype=np.float64))


def get_weights(
    mafs: np.ndarray,
    weight_spec: str,
    *,
    cadd_scores: np.ndarray | None = None,
    revel_scores: np.ndarray | None = None,
    variant_effects: np.ndarray | None = None,
    weight_params: dict | None = None,
) -> np.ndarray:
    """
    Parse a weight specification string and return the corresponding weights.

    Parameters
    ----------
    mafs : np.ndarray, shape (n_variants,)
        Minor allele frequencies, required for Beta weight computation.
        Passed to ``beta_maf_weights``; ignored for ``"uniform"``.
    weight_spec : str
        Specification string. Supported formats:

        - ``"beta:a,b"`` — calls ``beta_maf_weights(mafs, a, b)``
        - ``"uniform"`` — calls ``uniform_weights(len(mafs))``
        - ``"cadd"`` — calls ``cadd_weights(mafs, cadd_scores, ...)``
        - ``"revel"`` — calls ``revel_weights(mafs, revel_scores, ...)``
        - ``"combined"`` — calls ``combined_weights(mafs, cadd_scores, revel_scores, ...)``
    cadd_scores : np.ndarray or None, keyword-only
        CADD Phred scores. Required when ``weight_spec="cadd"``; optional for
        ``"combined"``. Ignored for ``"beta:*"`` and ``"uniform"`` specs.
    revel_scores : np.ndarray or None, keyword-only
        REVEL scores. Required when ``weight_spec="revel"``; optional for
        ``"combined"``. Ignored for ``"beta:*"`` and ``"uniform"`` specs.
    variant_effects : np.ndarray or None, keyword-only
        Variant effect strings for type-aware fallback logging. Optional.
        Ignored for ``"beta:*"`` and ``"uniform"`` specs.
    weight_params : dict or None, keyword-only
        Extra parameters for weight schemes. Currently supports:
        - ``"cadd_cap"`` (float): normalization cap for CADD Phred (default 40.0)

    Returns
    -------
    np.ndarray, shape (n_variants,), float64

    Raises
    ------
    ValueError
        If ``weight_spec`` does not match a known format, or if a functional
        spec (e.g. ``"cadd"``) is used without the required scores array.

    Examples
    --------
    >>> import numpy as np
    >>> mafs = np.array([0.01, 0.05, 0.10])
    >>> get_weights(mafs, "beta:1,25")
    array([...])
    >>> get_weights(mafs, "uniform")
    array([1., 1., 1.])
    """
    mafs_arr = np.asarray(mafs, dtype=np.float64)

    if weight_spec == "uniform":
        return uniform_weights(len(mafs_arr))

    if weight_spec.startswith("beta:"):
        params_str = weight_spec[len("beta:") :]
        try:
            parts = params_str.split(",")
            if len(parts) != 2:
                raise ValueError()
            a = float(parts[0].strip())
            b = float(parts[1].strip())
        except (ValueError, IndexError) as err:
            raise ValueError(
                f"Invalid beta weight spec '{weight_spec}'. "
                "Expected format: 'beta:a,b' where a and b are floats "
                "(e.g. 'beta:1,25')."
            ) from err
        return beta_maf_weights(mafs_arr, a=a, b=b)

    if weight_spec == "cadd":
        if cadd_scores is None:
            raise ValueError(
                "weight_spec='cadd' requires cadd_scores to be provided. "
                "Ensure CADD annotation columns (dbNSFP_CADD_phred) are present in your data."
            )
        cadd_cap = (weight_params or {}).get("cadd_cap", 40.0)
        return cadd_weights(mafs_arr, cadd_scores, variant_effects, cap=cadd_cap)

    if weight_spec == "revel":
        if revel_scores is None:
            raise ValueError(
                "weight_spec='revel' requires revel_scores to be provided. "
                "Ensure REVEL annotation columns (dbNSFP_REVEL_score) are present in your data."
            )
        return revel_weights(mafs_arr, revel_scores, variant_effects)

    if weight_spec == "combined":
        cadd_cap = (weight_params or {}).get("cadd_cap", 40.0)
        return combined_weights(
            mafs_arr,
            cadd_scores=cadd_scores,
            revel_scores=revel_scores,
            variant_effects=variant_effects,
            cadd_cap=cadd_cap,
        )

    raise ValueError(
        f"Unknown weight spec '{weight_spec}'. Supported specs: "
        "'beta:a,b', 'uniform', 'cadd', 'revel', 'combined'."
    )
