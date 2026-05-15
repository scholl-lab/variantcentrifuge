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
- ``score_column_weights``: Beta(MAF) x arbitrary numeric score-column weights,
  or raw normalized score weights when combine_with_beta is false.
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
    return np.asarray(_beta_dist.pdf(maf_clipped, a=a, b=b), dtype=np.float64)


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


def resolve_score_weight_column(
    weight_spec: str,
    variant_weight_column: str | None = None,
) -> str | None:
    """Resolve a score-column weight spec to a DataFrame column name."""
    if weight_spec.startswith("column:"):
        inline_column = weight_spec[len("column:") :].strip()
        if not inline_column:
            raise ValueError("column:<name> requires a non-empty column name")
        if variant_weight_column:
            logger.debug(
                "variant_weights=%r includes inline column %r; ignoring variant_weight_column=%r",
                weight_spec,
                inline_column,
                variant_weight_column,
            )
        return inline_column

    if weight_spec == "score_column":
        if not variant_weight_column:
            raise ValueError("score_column requires variant_weight_column")
        return variant_weight_column

    return None


_SCORE_WEIGHT_DEFAULTS = {
    "score_min": None,
    "score_max": None,
    "floor": 0.0,
    "ceiling": 1.0,
    "combine_with_beta": True,
    "missing": None,
    "beta_a": 1.0,
    "beta_b": 25.0,
}


def _normalize_score_weight_params(weight_params: dict | None) -> dict:
    if weight_params is not None and not isinstance(weight_params, dict):
        raise ValueError(
            f"variant_weight_params must be a dict, got {type(weight_params).__name__}"
        )

    params = dict(_SCORE_WEIGHT_DEFAULTS)
    params.update(weight_params or {})

    score_min = params["score_min"]
    score_max = params["score_max"]
    if (score_min is None) != (score_max is None):
        raise ValueError("score_min and score_max must be provided together")
    if score_min is not None and float(score_min) >= float(score_max):
        raise ValueError("score_min must be less than score_max")

    floor = float(params["floor"])
    ceiling = float(params["ceiling"])
    if floor < 0:
        raise ValueError("floor must be >= 0")
    if ceiling <= 0:
        raise ValueError("ceiling must be > 0")
    if floor > ceiling:
        raise ValueError("floor must be <= ceiling")

    beta_a = float(params["beta_a"])
    beta_b = float(params["beta_b"])
    if beta_a <= 0:
        raise ValueError("beta_a must be > 0")
    if beta_b <= 0:
        raise ValueError("beta_b must be > 0")

    missing = params["missing"]
    if missing not in (None, "neutral", "floor"):
        raise ValueError("missing must be one of None, 'neutral', or 'floor'")

    combine_with_beta = bool(params["combine_with_beta"])
    if missing == "neutral" and not combine_with_beta:
        raise ValueError("missing='neutral' is invalid when combine_with_beta=false")

    params["floor"] = floor
    params["ceiling"] = ceiling
    params["beta_a"] = beta_a
    params["beta_b"] = beta_b
    params["combine_with_beta"] = combine_with_beta
    if score_min is not None:
        params["score_min"] = float(score_min)
        params["score_max"] = float(score_max)
    return params


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
    assert scores_f is not None
    nan_mask = np.isnan(scores_f)

    functional = np.where(nan_mask, 1.0, np.minimum(scores_f / cap, 1.0))

    _log_missing_score_counts(nan_mask, variant_effects, "CADD")

    return np.asarray(maf_w * functional, dtype=np.float64)


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
    assert scores_f is not None
    nan_mask = np.isnan(scores_f)

    functional = np.where(nan_mask, 1.0, scores_f)

    _log_missing_score_counts(nan_mask, variant_effects, "REVEL")

    return np.asarray(maf_w * functional, dtype=np.float64)


def score_column_weights(
    mafs: np.ndarray,
    score_values: np.ndarray,
    *,
    variant_effects: np.ndarray | None = None,
    weight_params: dict | None = None,
) -> np.ndarray:
    """Compute weights from an arbitrary numeric variant score column."""
    mafs_arr = np.asarray(mafs, dtype=np.float64)
    scores_f = _parse_scores_to_float(np.asarray(score_values, dtype=object))
    assert scores_f is not None

    if len(mafs_arr) != len(scores_f):
        raise ValueError(
            f"score_values and mafs must have the same length "
            f"(got {len(scores_f)} and {len(mafs_arr)})"
        )

    params = _normalize_score_weight_params(weight_params)
    finite_mask = np.isfinite(scores_f)
    missing_mask = ~finite_mask

    normalized = np.empty(len(scores_f), dtype=np.float64)
    normalized[:] = np.nan

    if params["score_min"] is not None:
        score_min = float(params["score_min"])
        score_max = float(params["score_max"])
        normalized[finite_mask] = (scores_f[finite_mask] - score_min) / (
            score_max - score_min
        )
    else:
        normalized[finite_mask] = scores_f[finite_mask]
        out_of_unit = finite_mask & ((scores_f < 0.0) | (scores_f > 1.0))
        if bool(out_of_unit.any()):
            logger.warning(
                "score_column weights: %d finite score(s) outside [0, 1] without explicit range",
                int(out_of_unit.sum()),
            )

    functional = np.clip(normalized, params["floor"], params["ceiling"])

    missing_mode = params["missing"]
    if missing_mode is None:
        missing_mode = "neutral" if params["combine_with_beta"] else "floor"

    if missing_mode == "neutral":
        functional[missing_mask] = 1.0
    else:
        functional[missing_mask] = params["floor"]

    _log_missing_score_counts(missing_mask, variant_effects, "score_column")

    if bool(missing_mask.all()) and len(missing_mask) > 0:
        logger.warning("score_column weights: all score values for this gene are missing or invalid")

    if params["combine_with_beta"]:
        maf_w = beta_maf_weights(mafs_arr, a=params["beta_a"], b=params["beta_b"])
        return np.asarray(maf_w * functional, dtype=np.float64)

    return np.asarray(functional, dtype=np.float64)


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
    score_values: np.ndarray | None = None,
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

    if weight_spec == "score_column" or weight_spec.startswith("column:"):
        if weight_spec.startswith("column:") and not weight_spec[len("column:") :].strip():
            raise ValueError("column:<name> requires a non-empty column name")
        if score_values is None:
            raise ValueError(f"weight_spec='{weight_spec}' requires score_values to be provided")
        return score_column_weights(
            mafs_arr,
            score_values,
            variant_effects=variant_effects,
            weight_params=weight_params,
        )

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
        "'beta:a,b', 'uniform', 'cadd', 'revel', 'combined', "
        "'column:<name>', 'score_column'."
    )
