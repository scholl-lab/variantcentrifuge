# File: variantcentrifuge/association/tests/allelic_series_classify.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/allelic_series_classify.py
"""
Variant classification logic for the COAST allelic series test.

Classifies variants into four categories:
  - Code 3 (PTV): HIGH impact AND effect in {stop_gained, frameshift_variant,
    splice_acceptor_variant, splice_donor_variant}
  - Code 2 (DMV): missense_variant AND (SIFT deleterious OR PolyPhen probably/possibly_damaging)
  - Code 1 (BMV): missense_variant AND SIFT tolerated AND PolyPhen benign
  - Code 0 (unclassified): everything else, including missense without predictions

Variants with code 0 are excluded from COAST (include_mask=False) but are still
included in other tests (SKAT, burden) via the standard genotype matrix path.
"""

from __future__ import annotations

import logging
import os
from typing import Any

import numpy as np
import pandas as pd

from variantcentrifuge.scoring import apply_scoring, read_scoring_config

logger = logging.getLogger("variantcentrifuge")

# ---------------------------------------------------------------------------
# Variant classification constants
# ---------------------------------------------------------------------------

PTV_EFFECTS = frozenset(
    {
        "stop_gained",
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
    }
)
MISSENSE_EFFECT = "missense_variant"

# dbNSFP uses single-letter codes: D=damaging/deleterious, T=tolerated, B=benign,
# P=possibly_damaging. Full words used in some annotation pipelines.
SIFT_DAMAGING = frozenset({"deleterious", "D"})
POLYPHEN_DAMAGING = frozenset({"probably_damaging", "possibly_damaging", "D", "P"})
SIFT_BENIGN = frozenset({"tolerated", "T"})
POLYPHEN_BENIGN = frozenset({"benign", "B"})

# Possible column names for SIFT predictions (tried in order)
SIFT_COLUMN_CANDIDATES = ["dbNSFP_SIFT_pred", "SIFT_pred", "sift_pred"]

# Possible column names for PolyPhen predictions (tried in order)
POLYPHEN_COLUMN_CANDIDATES = [
    "dbNSFP_Polyphen2_HDIV_pred",
    "dbNSFP_Polyphen2_HVAR_pred",
    "Polyphen2_HDIV_pred",
    "Polyphen2_HVAR_pred",
    "polyphen2_hdiv_pred",
    "polyphen2_hvar_pred",
]

# Possible EFFECT column names (tried in order)
EFFECT_COLUMN_CANDIDATES = ["EFFECT", "ANN_0__EFFECT", "effect", "ann_0__effect"]

# Possible IMPACT column names (tried in order)
IMPACT_COLUMN_CANDIDATES = ["IMPACT", "ANN_0__IMPACT", "impact", "ann_0__impact"]

# Possible CADD phred column names (tried in order)
CADD_COLUMN_CANDIDATES = ["dbNSFP_CADD_phred", "CADD_phred", "cadd_phred", "CADD_PHRED"]


# ---------------------------------------------------------------------------
# Effect string resolution (COAST-04 / COAST-07)
# ---------------------------------------------------------------------------


def _resolve_effect(effect_str: str) -> str:
    """Resolve '&'-concatenated SnpEff effect string to single highest-priority effect.

    Priority: PTV effects > missense_variant > first part.
    This handles multi-transcript annotations like 'stop_gained&splice_region_variant'.

    Parameters
    ----------
    effect_str : str
        Raw SnpEff effect string, possibly containing '&'-delimited multi-transcript
        annotations (e.g. 'stop_gained&splice_region_variant').

    Returns
    -------
    str
        The single highest-priority effect. Returns the input unchanged if no '&'
        is present (fast path).
    """
    if "&" not in effect_str:
        return effect_str
    parts = [p.strip() for p in effect_str.split("&")]
    for part in parts:
        if part in PTV_EFFECTS:
            return part
    if MISSENSE_EFFECT in parts:
        return MISSENSE_EFFECT
    return parts[0] if parts else effect_str


# ---------------------------------------------------------------------------
# Formula engine classification helper (COAST-06)
# ---------------------------------------------------------------------------


def _classify_via_formula_engine(
    gene_df: pd.DataFrame,
    effect_col: str,
    impact_col: str,
    model_dir: str,
    diagnostics_rows: list[dict[str, Any]] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Classify variants using a scoring formula engine config directory.

    Builds a normalized copy of gene_df with standard column names
    (COAST_EFFECT, COAST_IMPACT, COAST_SIFT, COAST_POLYPHEN, COAST_CADD),
    then calls read_scoring_config() + apply_scoring() to compute coast_category.

    Parameters
    ----------
    gene_df : pd.DataFrame
        Per-variant DataFrame for one gene.
    effect_col : str
        Detected effect column name in gene_df.
    impact_col : str
        Detected impact column name in gene_df.
    model_dir : str
        Absolute path to the scoring config directory.
    diagnostics_rows : list | None
        If not None, per-variant diagnostic dicts are appended here.

    Returns
    -------
    annotation_codes : np.ndarray of int
    include_mask : np.ndarray of bool
    """
    n = len(gene_df)

    # Build normalized working copy with standard column names
    work = pd.DataFrame(index=gene_df.index)

    # Effect column: resolve '&'-concatenated strings (COAST-04)
    work["COAST_EFFECT"] = (
        gene_df[effect_col].astype(object).fillna("").astype(str).map(_resolve_effect)
    )
    work["COAST_IMPACT"] = gene_df[impact_col].astype(object).fillna("").astype(str)

    # SIFT and PolyPhen (for sift_polyphen model)
    cols = set(gene_df.columns)
    for candidate in SIFT_COLUMN_CANDIDATES:
        if candidate in cols:
            work["COAST_SIFT"] = gene_df[candidate].astype(object).fillna("").astype(str)
            break
    if "COAST_SIFT" not in work.columns:
        work["COAST_SIFT"] = ""

    for candidate in POLYPHEN_COLUMN_CANDIDATES:
        if candidate in cols:
            work["COAST_POLYPHEN"] = gene_df[candidate].astype(object).fillna("").astype(str)
            break
    if "COAST_POLYPHEN" not in work.columns:
        work["COAST_POLYPHEN"] = ""

    # CADD (for cadd model)
    for candidate in CADD_COLUMN_CANDIDATES:
        if candidate in cols:
            cadd_raw = gene_df[candidate].astype(object).where(gene_df[candidate].notna(), other=0)
            work["COAST_CADD"] = pd.to_numeric(cadd_raw, errors="coerce").fillna(0.0)
            break
    if "COAST_CADD" not in work.columns:
        work["COAST_CADD"] = 0.0

    # Load and apply scoring config
    try:
        scoring_config = read_scoring_config(model_dir)
    except Exception as exc:
        logger.error(f"classify_variants: failed to load model from '{model_dir}': {exc}")
        raise ValueError(
            f"classify_variants: failed to load model from '{model_dir}': {exc}"
        ) from exc

    try:
        scored = apply_scoring(work, scoring_config, strict=True)
    except Exception as exc:
        logger.error(f"classify_variants: apply_scoring failed for model '{model_dir}': {exc}")
        raise ValueError(
            f"classify_variants: apply_scoring failed for model '{model_dir}': {exc}"
        ) from exc

    if "coast_category" not in scored.columns:
        logger.error(
            f"classify_variants: model '{model_dir}' did not produce 'coast_category' column. "
            f"Available columns: {list(scored.columns)}"
        )
        raise ValueError(
            f"classify_variants: model '{model_dir}' did not produce 'coast_category' column. "
            f"Available columns: {list(scored.columns)}"
        )

    if scored["coast_category"].isna().any():
        logger.error(
            f"classify_variants: model '{model_dir}' produced missing/NaN category outputs."
        )
        raise ValueError(
            f"classify_variants: model '{model_dir}' produced missing/NaN category outputs."
        )

    annotation_codes = scored["coast_category"].astype(int).to_numpy()
    include_mask = annotation_codes > 0

    # Diagnostics (optional)
    if diagnostics_rows is not None:
        model_name = os.path.basename(model_dir)
        original_effects = gene_df[effect_col].astype(object).fillna("").astype(str)
        resolved_effects = work["COAST_EFFECT"]
        for i in range(n):
            idx_val = gene_df.index[i]
            row = gene_df.iloc[i]
            chrom = str(row.get("CHROM", row.get("chrom", idx_val)))
            pos = str(row.get("POS", row.get("pos", "")))
            ref = str(row.get("REF", row.get("ref", "")))
            alt = str(row.get("ALT", row.get("alt", "")))
            variant_id = f"{chrom}:{pos}:{ref}:{alt}" if pos else str(idx_val)
            diagnostics_rows.append(
                {
                    "gene": row.get("GENE", row.get("gene", "")),
                    "variant_id": variant_id,
                    "original_effect": original_effects.iloc[i],
                    "resolved_effect": resolved_effects.iloc[i],
                    "assigned_category": int(annotation_codes[i]),
                    "model_name": model_name,
                }
            )

    return annotation_codes, include_mask


# ---------------------------------------------------------------------------
# Public classification function
# ---------------------------------------------------------------------------


def classify_variants(
    gene_df: pd.DataFrame,
    effect_col: str,
    impact_col: str,
    model_dir: str | None = None,
    diagnostics_rows: list[dict[str, Any]] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Classify variants into COAST categories (BMV/DMV/PTV/unclassified).

    Uses SnpEff EFFECT/IMPACT columns and dbNSFP SIFT/PolyPhen predictions.
    SIFT and PolyPhen column names are auto-detected from gene_df.columns.

    Parameters
    ----------
    gene_df : pd.DataFrame
        Per-variant DataFrame for one gene. Must contain effect_col and
        impact_col. May contain SIFT/PolyPhen columns (auto-detected).
    effect_col : str
        Column name for variant effect (e.g. "EFFECT" or "ANN_0__EFFECT").
    impact_col : str
        Column name for variant impact (e.g. "IMPACT" or "ANN_0__IMPACT").
    model_dir : str | None
        Path to a scoring/coast_classification/<model>/ directory. When provided,
        delegates classification to the scoring formula engine (read_scoring_config +
        apply_scoring). When None, uses the built-in SIFT/PolyPhen hardcoded logic.
    diagnostics_rows : list | None
        If provided (not None), per-variant diagnostic dicts are appended here.
        Each dict contains: gene, variant_id, original_effect, resolved_effect,
        assigned_category, model_name.

    Returns
    -------
    annotation_codes : np.ndarray of int
        Integer code per variant: 3=PTV, 2=DMV, 1=BMV, 0=unclassified.
    include_mask : np.ndarray of bool
        True for variants included in COAST (BMV/DMV/PTV only).
        Unclassified variants (code 0) have include_mask=False.

    Notes
    -----
    Missense variants without SIFT/PolyPhen predictions receive code 0 and
    are excluded from COAST only. They are still used by SKAT/burden via the
    standard genotype matrix path.

    '&'-concatenated SnpEff multi-transcript effects are resolved via
    _resolve_effect() before classification (COAST-04/COAST-07).

    If no SIFT or PolyPhen columns are found at all (default path), logs an error and
    returns all-zero codes with all-False include mask.
    """
    n = len(gene_df)
    annotation_codes = np.zeros(n, dtype=int)

    if model_dir is not None:
        # ── Formula engine path (COAST-06) ──────────────────────────────────────
        return _classify_via_formula_engine(
            gene_df=gene_df,
            effect_col=effect_col,
            impact_col=impact_col,
            model_dir=model_dir,
            diagnostics_rows=diagnostics_rows,
        )

    # ── Default SIFT/PolyPhen hardcoded path ─────────────────────────────────

    # Auto-detect SIFT and PolyPhen column names
    cols = set(gene_df.columns)
    sift_col: str | None = None
    polyphen_col: str | None = None

    for candidate in SIFT_COLUMN_CANDIDATES:
        if candidate in cols:
            sift_col = candidate
            break

    for candidate in POLYPHEN_COLUMN_CANDIDATES:
        if candidate in cols:
            polyphen_col = candidate
            break

    if sift_col is None and polyphen_col is None:
        logger.error(
            "classify_variants: no SIFT or PolyPhen columns found. "
            f"Looked for SIFT candidates {SIFT_COLUMN_CANDIDATES} and "
            f"PolyPhen candidates {POLYPHEN_COLUMN_CANDIDATES}. "
            "All variants will be excluded from COAST."
        )
        return annotation_codes, np.zeros(n, dtype=bool)

    # Extract effect and impact series
    # Convert to object dtype first to handle Categorical columns safely
    # Apply _resolve_effect to handle '&'-concatenated multi-transcript annotations (COAST-04)
    effect_series = gene_df[effect_col].astype(object).fillna("").astype(str).map(_resolve_effect)
    impact_series = gene_df[impact_col].astype(object).fillna("").astype(str)

    # Helper: extract first prediction value from potentially multi-value cell
    # dbNSFP fields can be e.g. "deleterious;tolerated" — take first non-null
    def _get_pred(series: Any, val_set: frozenset[str]) -> np.ndarray:
        """Return boolean array: True where series contains any value in val_set."""
        result = np.zeros(n, dtype=bool)
        for i, val in enumerate(series):
            if pd.isna(val) or val == "" or val == ".":
                continue
            # Multi-value: split on semicolons or commas
            parts = str(val).replace(",", ";").split(";")
            for part in parts:
                part = part.strip()
                if part in val_set:
                    result[i] = True
                    break
        return result

    # Classify PTV: impact == "HIGH" AND effect in PTV_EFFECTS
    is_high_impact = impact_series.str.strip() == "HIGH"
    is_ptv_effect = effect_series.str.strip().isin(PTV_EFFECTS)
    is_ptv = is_high_impact & is_ptv_effect

    # Classify missense
    is_missense = effect_series.str.strip() == MISSENSE_EFFECT

    # Get SIFT and PolyPhen predictions (vectorized)
    sift_damaging = np.zeros(n, dtype=bool)
    sift_benign = np.zeros(n, dtype=bool)
    polyphen_damaging = np.zeros(n, dtype=bool)
    polyphen_benign = np.zeros(n, dtype=bool)

    if sift_col is not None:
        sift_series = gene_df[sift_col].astype(object).fillna("").astype(str)
        sift_damaging = _get_pred(sift_series, SIFT_DAMAGING)
        sift_benign = _get_pred(sift_series, SIFT_BENIGN)

    if polyphen_col is not None:
        polyphen_series = gene_df[polyphen_col].astype(object).fillna("").astype(str)
        polyphen_damaging = _get_pred(polyphen_series, POLYPHEN_DAMAGING)
        polyphen_benign = _get_pred(polyphen_series, POLYPHEN_BENIGN)

    # DMV: missense AND (SIFT damaging OR PolyPhen damaging)
    is_dmv = is_missense & (sift_damaging | polyphen_damaging)

    # BMV: missense AND SIFT benign AND PolyPhen benign
    is_bmv = is_missense & sift_benign & polyphen_benign

    # Assign codes (PTV takes priority over DMV/BMV if a variant is both)
    annotation_codes[is_bmv.to_numpy() if hasattr(is_bmv, "to_numpy") else is_bmv] = 1
    annotation_codes[is_dmv.to_numpy() if hasattr(is_dmv, "to_numpy") else is_dmv] = 2
    annotation_codes[is_ptv.to_numpy() if hasattr(is_ptv, "to_numpy") else is_ptv] = 3

    include_mask = annotation_codes > 0

    # Warn about excluded missense variants
    is_missense_arr = is_missense.to_numpy() if hasattr(is_missense, "to_numpy") else is_missense
    n_excluded_missense = int(np.sum(is_missense_arr & ~include_mask))
    if n_excluded_missense > 0:
        logger.warning(
            f"classify_variants: {n_excluded_missense} missense variant(s) excluded from COAST "
            "due to missing or ambiguous SIFT/PolyPhen predictions."
        )

    return annotation_codes, include_mask
