# File: variantcentrifuge/association/tests/allelic_series.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/allelic_series.py
"""
COASTTest — AssociationTest wrapper for the R AllelicSeries::COAST() test.

Implements the Combined Omnibus Association Test (COAST) for allelic series
testing. COAST has higher power than standard burden/SKAT when effect sizes
are ordered by variant damage severity: BMV (benign missense) < DMV (damaging
missense) < PTV (protein-truncating variants). This ordered alternative
provides substantial power gains for variants with monotone effect size trends.

Reference: McCaw et al. AJHG 2023. "Operating within the context of an allelic
series, COAST maintains type 1 error control and demonstrates superior power
compared to burden and SKAT tests."

Variant Classification
-----------------------
classify_variants() classifies each variant into one of four categories:
  - Code 3 (PTV): HIGH impact AND effect in {stop_gained, frameshift_variant,
    splice_acceptor_variant, splice_donor_variant}
  - Code 2 (DMV): missense_variant AND (SIFT deleterious OR PolyPhen probably/possibly_damaging)
  - Code 1 (BMV): missense_variant AND SIFT tolerated AND PolyPhen benign
  - Code 0 (unclassified): everything else, including missense without predictions

Variants with code 0 are excluded from COAST (include_mask=False) but are still
included in other tests (SKAT, burden) via the standard genotype matrix path.

Thread safety
--------------
COASTTest is not thread-safe (rpy2 restriction). The stage that uses COASTTest
must declare parallel_safe=False.

Lifecycle hooks
----------------
prepare(n_genes): logs start message, sets up progress tracking
finalize(): timing summary, R gc() cleanup
"""

from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING, Any

import numpy as np

from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult

if TYPE_CHECKING:
    pass

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

# COAST-specific R GC interval (same as RSKATTest)
_GC_INTERVAL = 100


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
    gene_df: Any,
    effect_col: str,
    impact_col: str,
    model_dir: str,
    diagnostics_rows: list | None = None,
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
    import os

    import pandas as pd

    from variantcentrifuge.scoring import apply_scoring, read_scoring_config

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
        return np.zeros(n, dtype=int), np.zeros(n, dtype=bool)

    try:
        scored = apply_scoring(work, scoring_config)
    except Exception as exc:
        logger.error(f"classify_variants: apply_scoring failed for model '{model_dir}': {exc}")
        return np.zeros(n, dtype=int), np.zeros(n, dtype=bool)

    if "coast_category" not in scored.columns:
        logger.error(
            f"classify_variants: model '{model_dir}' did not produce 'coast_category' column. "
            f"Available columns: {list(scored.columns)}"
        )
        return np.zeros(n, dtype=int), np.zeros(n, dtype=bool)

    annotation_codes = scored["coast_category"].fillna(0).astype(int).to_numpy()
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
    gene_df: Any,
    effect_col: str,
    impact_col: str,
    model_dir: str | None = None,
    diagnostics_rows: list | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Classify variants into COAST categories (BMV/DMV/PTV/unclassified).

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
    import pandas as pd

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


# ---------------------------------------------------------------------------
# COASTTest class
# ---------------------------------------------------------------------------


class COASTTest(AssociationTest):
    """
    # DEPRECATED
    COAST allelic series test via R AllelicSeries::COAST().

    Registered in the engine as ``"coast"``. Classifies variants into
    BMV/DMV/PTV categories using SnpEff EFFECT/IMPACT and dbNSFP SIFT/PolyPhen
    annotations, then invokes the R AllelicSeries COAST omnibus test.

    COAST combines a burden and a SKAT component that both test the allelic
    series hypothesis. The omnibus p-value from COAST feeds into ACAT-O.
    The burden and SKAT component p-values are stored in extra for diagnostics.

    Parameters
    ----------
    None — instantiated by AssociationEngine.from_names() via the registry.

    Raises
    ------
    ImportError
        From check_dependencies() if rpy2 is not importable or the
        AllelicSeries R package is not installed.

    Examples
    --------
    This test is typically used through the engine, not instantiated directly:

    >>> engine = AssociationEngine.from_names(["coast"], config)
    >>> result_df = engine.run_all(gene_burden_data)
    """

    parallel_safe: bool = False  # rpy2 restriction: main thread only

    def __init__(self) -> None:
        import warnings

        warnings.warn(
            "COASTTest (R COAST backend) is deprecated and will be removed in v0.17.0. "
            "Use --coast-backend python (default).",
            DeprecationWarning,
            stacklevel=2,
        )
        # Lifecycle tracking (set by prepare())
        self._total_genes: int = 0
        self._log_interval: int = 50
        self._genes_processed: int = 0
        self._start_time: float = 0.0

        # Category-level counters (set by prepare(), incremented in run())
        self._n_complete: int = 0
        self._n_partial: int = 0
        self._n_skipped: int = 0

    @property
    def name(self) -> str:
        """Short identifier for this test used in registry lookup and column prefixes."""
        return "coast"

    def check_dependencies(self) -> None:
        """
        Verify that rpy2 is importable and the R AllelicSeries package is installed.

        Raises
        ------
        ImportError
            If rpy2 is not installed, R is not found, or the AllelicSeries
            R package is not installed.
        """
        try:
            import rpy2.robjects  # noqa: F401
        except ImportError as exc:
            raise ImportError(
                "COASTTest requires rpy2. Install with: pip install rpy2\n"
                "Also ensure R is installed and R_HOME is set correctly."
            ) from exc

        try:
            from rpy2.robjects.packages import importr

            importr("AllelicSeries")
        except Exception as exc:
            raise ImportError(
                "COASTTest requires the R AllelicSeries package. "
                "Install in R with: install.packages('AllelicSeries')\n"
                "Or from GitHub: remotes::install_github('zhangz19/AllelicSeries')\n"
                f"Original error: {exc}"
            ) from exc

    def effect_column_names(self) -> dict[str, str | None]:
        """
        Column name suffixes for COAST output.

        COAST produces only p-values (omnibus + burden/SKAT components in extra).
        There is no single effect size. All four effect slots are None; the
        engine's None-guard skips column creation for them.

        Returns
        -------
        dict with all values None.
        """
        return {
            "effect": None,
            "se": None,
            "ci_lower": None,
            "ci_upper": None,
        }

    def prepare(self, gene_count: int) -> None:
        """
        Called by the engine before the per-gene loop.

        Sets up progress tracking and logs the start message.

        Parameters
        ----------
        gene_count : int
            Total number of genes that will be processed.
        """
        self._total_genes = gene_count
        self._genes_processed = 0
        self._log_interval = max(10, min(50, gene_count // 10)) if gene_count > 0 else 50
        self._start_time = time.time()
        self._n_complete = 0
        self._n_partial = 0
        self._n_skipped = 0

        logger.info(f"COAST: beginning allelic series analysis of {gene_count} genes")

    def finalize(self) -> None:
        """
        Called by the engine after the per-gene loop completes.

        Logs aggregate timing summary and triggers R garbage collection.
        """
        elapsed = time.time() - self._start_time
        n_tested = self._n_complete + self._n_partial
        logger.info(f"COAST complete: {self._genes_processed} genes in {elapsed:.1f}s")
        logger.info(
            f"COAST: {n_tested} genes tested "
            f"({self._n_complete} complete, {self._n_partial} partial, "
            f"{self._n_skipped} skipped)"
        )

        # Final R gc() to free heap
        try:
            import rpy2.robjects as ro

            ro.r("gc()")
        except Exception:
            pass  # R cleanup is best-effort

    def run(
        self,
        gene: str,
        contingency_data: dict[str, Any],
        config: AssociationConfig,
    ) -> TestResult:
        """
        Run COAST for a single gene.

        Extracts EFFECT/IMPACT annotations from the gene DataFrame, classifies
        variants into BMV/DMV/PTV categories, then calls R AllelicSeries::COAST().

        Parameters
        ----------
        gene : str
            Gene symbol being tested.
        contingency_data : dict
            Must contain:
              - ``genotype_matrix``  : np.ndarray (n_samples, n_variants)
              - ``phenotype_vector`` : np.ndarray (n_samples,)
              - ``gene_df``          : pd.DataFrame with per-variant annotation rows
            Optional:
              - ``covariate_matrix`` : np.ndarray (n_samples, k) or None
        config : AssociationConfig
            Runtime configuration including coast_weights and trait_type.

        Returns
        -------
        TestResult
            p_value=None when test is skipped (missing variant categories,
            no genotype matrix, or insufficient data).
            extra contains: coast_burden_pvalue, coast_skat_pvalue,
            coast_n_bmv, coast_n_dmv, coast_n_ptv.
        """
        n_cases = int(contingency_data.get("proband_count", 0))
        n_controls = int(contingency_data.get("control_count", 0))
        n_variants = int(contingency_data.get("n_qualifying_variants", 0))

        # Require genotype matrix and phenotype
        if "genotype_matrix" not in contingency_data:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "NO_GENOTYPE_MATRIX"},
            )

        geno: np.ndarray = contingency_data["genotype_matrix"]
        phenotype: np.ndarray | None = contingency_data.get("phenotype_vector")
        covariate_matrix: np.ndarray | None = contingency_data.get("covariate_matrix")
        gene_df = contingency_data.get("gene_df")

        if phenotype is None or geno.shape[0] == 0 or geno.shape[1] == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "NO_PHENOTYPE_OR_EMPTY_MATRIX"},
            )

        if gene_df is None or len(gene_df) == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "NO_GENE_DF"},
            )

        # Find EFFECT and IMPACT columns in gene_df
        effect_col = self._find_column(gene_df, EFFECT_COLUMN_CANDIDATES, "EFFECT")
        impact_col = self._find_column(gene_df, IMPACT_COLUMN_CANDIDATES, "IMPACT")

        if effect_col is None or impact_col is None:
            logger.warning(
                f"COAST [{gene}]: could not find EFFECT ({effect_col}) or "
                f"IMPACT ({impact_col}) columns. Skipping."
            )
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "MISSING_EFFECT_IMPACT_COLUMNS"},
            )

        # Classify variants into COAST categories
        # Pass model_dir when coast_classification config is set
        coast_model_dir = getattr(config, "coast_classification", None)
        anno_codes, include_mask = classify_variants(
            gene_df, effect_col, impact_col, model_dir=coast_model_dir
        )

        # Check if gene_df length aligns with genotype matrix columns
        # (gene_df may be a subset due to site filtering in genotype matrix builder)
        if len(anno_codes) != geno.shape[1]:
            logger.warning(
                f"COAST [{gene}]: annotation codes length ({len(anno_codes)}) "
                f"does not match genotype matrix columns ({geno.shape[1]}). "
                "This may happen when site filtering removed variants. Skipping."
            )
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "ANNOTATION_GENOTYPE_MISMATCH"},
            )

        # Filter to COAST-eligible variants only
        anno_filtered = anno_codes[include_mask]
        geno_filtered = geno[:, include_mask]

        if geno_filtered.shape[1] == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"coast_skip_reason": "NO_CLASSIFIABLE_VARIANTS"},
            )

        # Count variants per category
        n_bmv = int(np.sum(anno_filtered == 1))
        n_dmv = int(np.sum(anno_filtered == 2))
        n_ptv = int(np.sum(anno_filtered == 3))

        # Partial-category logic: skip only when ALL categories are empty
        missing = [cat for cat, n in [("BMV", n_bmv), ("DMV", n_dmv), ("PTV", n_ptv)] if n == 0]
        present = [cat for cat, n in [("BMV", n_bmv), ("DMV", n_dmv), ("PTV", n_ptv)] if n > 0]

        if len(present) == 0:
            self._n_skipped += 1
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={
                    "coast_skip_reason": "ALL_CATEGORIES_EMPTY",
                    "coast_n_bmv": 0,
                    "coast_n_dmv": 0,
                    "coast_n_ptv": 0,
                    "coast_status": "skipped",
                },
            )

        coast_status = "complete" if not missing else "partial"
        if missing:
            logger.debug(
                f"COAST [{gene}]: partial -- missing {', '.join(missing)} "
                f"(BMV={n_bmv}, DMV={n_dmv}, PTV={n_ptv})"
            )

        # Run R AllelicSeries::COAST()
        try:
            burden_p, skat_p, omnibus_p = self._run_r_coast(
                geno_filtered=geno_filtered,
                anno_codes_filtered=anno_filtered,
                phenotype=phenotype,
                covariate_matrix=covariate_matrix,
                config=config,
            )
        except Exception as exc:
            logger.error(
                f"COAST [{gene}]: R AllelicSeries::COAST() failed: {type(exc).__name__}: {exc!r}"
            )
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={
                    "coast_skip_reason": f"R_ERROR:{exc!s}",
                    "coast_n_bmv": n_bmv,
                    "coast_n_dmv": n_dmv,
                    "coast_n_ptv": n_ptv,
                },
            )

        # Increment counter and periodic GC
        self._genes_processed += 1
        if coast_status == "complete":
            self._n_complete += 1
        else:
            self._n_partial += 1
        if self._genes_processed % _GC_INTERVAL == 0:
            try:
                import rpy2.robjects as ro

                ro.r("gc()")
            except Exception:
                pass

        # Progress logging
        if self._total_genes > 0 and self._genes_processed % self._log_interval == 0:
            pct = 100.0 * self._genes_processed / self._total_genes
            logger.info(
                f"COAST progress: {self._genes_processed}/{self._total_genes} genes ({pct:.0f}%)"
            )

        extra_base: dict[str, Any] = {
            "coast_burden_pvalue": burden_p,
            "coast_skat_pvalue": skat_p,
            "coast_n_bmv": n_bmv,
            "coast_n_dmv": n_dmv,
            "coast_n_ptv": n_ptv,
            "coast_status": coast_status,
        }
        if missing:
            extra_base["coast_missing_categories"] = ",".join(missing)

        return TestResult(
            gene=gene,
            test_name=self.name,
            p_value=omnibus_p,
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            se=None,
            n_cases=n_cases,
            n_controls=n_controls,
            n_variants=n_variants,
            extra=extra_base,
        )

    def _find_column(self, gene_df: Any, candidates: list[str], label: str) -> str | None:
        """Find the first available column name from a list of candidates."""
        cols = set(gene_df.columns)
        for candidate in candidates:
            if candidate in cols:
                return candidate
        logger.debug(
            f"COASTTest: could not find {label} column. "
            f"Tried: {candidates}. Available: {sorted(cols)}"
        )
        return None

    def _run_r_coast(
        self,
        geno_filtered: np.ndarray,
        anno_codes_filtered: np.ndarray,
        phenotype: np.ndarray,
        covariate_matrix: np.ndarray | None,
        config: AssociationConfig,
    ) -> tuple[float | None, float | None, float | None]:
        """
        Call R AllelicSeries::COAST() via rpy2 and extract p-values.

        Parameters
        ----------
        geno_filtered : np.ndarray (n_samples, n_variants_filtered)
            Genotype matrix for COAST-eligible variants only.
        anno_codes_filtered : np.ndarray (n_variants_filtered,)
            Integer annotation codes (1=BMV, 2=DMV, 3=PTV).
        phenotype : np.ndarray (n_samples,)
            Phenotype vector.
        covariate_matrix : np.ndarray or None
            Covariate matrix (n_samples, k) or None.
        config : AssociationConfig
            Runtime configuration.

        Returns
        -------
        (burden_p, skat_p, omnibus_p)
            Three p-values from COAST Pvals slot.
        """
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr

        allelic_series = importr("AllelicSeries")

        n_samples, n_variants_filtered = geno_filtered.shape

        # Build R matrix: AllelicSeries expects genotype as (n_samples, n_variants)
        r_geno = ro.r.matrix(
            ro.FloatVector(geno_filtered.T.ravel()),
            nrow=n_samples,
            ncol=n_variants_filtered,
        )
        r_anno = ro.IntVector(anno_codes_filtered.tolist())
        r_pheno = ro.FloatVector(phenotype.tolist())

        coast_weights = getattr(config, "coast_weights", None) or [1.0, 2.0, 3.0]
        r_weights = ro.FloatVector(coast_weights)

        is_binary = getattr(config, "trait_type", "binary") == "binary"
        r_is_binary = ro.BoolVector([is_binary])

        # Build covariate argument
        if covariate_matrix is not None and covariate_matrix.ndim == 2:
            n_cov = covariate_matrix.shape[1]
            r_covar = ro.r.matrix(
                ro.FloatVector(covariate_matrix.T.ravel()),
                nrow=n_samples,
                ncol=n_cov,
            )
        else:
            r_covar = ro.NULL

        result = allelic_series.COAST(
            anno=r_anno,
            geno=r_geno,
            pheno=r_pheno,
            covar=r_covar,
            weights=r_weights,
            is_pheno_binary=r_is_binary,
        )

        # Extract p-values from Pvals slot.
        # AllelicSeries COAST() returns Pvals as a data.frame with columns:
        #   test (str), type (str), pval (numeric)
        # Rows: baseline, ind, max_count, max_ind, sum_count, sum_ind,
        #        allelic_skat, omni
        pvals_slot = result.slots["Pvals"]

        # Convert R data.frame to Python dict keyed by test name
        pvals_dict: dict[str, float | None] = {}
        try:
            # R data.frame columns accessible by name
            test_names = list(pvals_slot.rx2("test"))
            pval_values = list(pvals_slot.rx2("pval"))
            import rpy2.rinterface as ri

            for name, val in zip(test_names, pval_values, strict=False):
                if hasattr(ri, "NA_Real") and val is ri.NA_Real:
                    pvals_dict[str(name)] = None
                elif hasattr(val, "__float__"):
                    pvals_dict[str(name)] = float(val)
                else:
                    pvals_dict[str(name)] = float(val) if val is not None else None
        except Exception:
            # Fallback: try as named vector (older AllelicSeries versions)
            pvals_names = list(pvals_slot.names) if hasattr(pvals_slot, "names") else []
            pvals_values = list(pvals_slot)
            for name, val in zip(pvals_names, pvals_values, strict=False):
                try:
                    pvals_dict[str(name)] = float(val) if val is not None else None
                except (TypeError, ValueError):
                    pvals_dict[str(name)] = None

        logger.debug(f"COAST R Pvals: {pvals_dict}")

        # Map to burden/skat/omnibus using AllelicSeries test names
        burden_p = _extract_pval(
            pvals_dict,
            ["Burden", "Burden_p", "burden", "burden_p", "baseline"],
        )
        skat_p = _extract_pval(
            pvals_dict,
            ["SKAT", "SKAT_p", "skat", "skat_p", "allelic_skat"],
        )
        omnibus_p = _extract_pval(
            pvals_dict,
            ["Omni", "Omni_p", "omni", "omni_p", "Omnibus", "omnibus"],
        )

        # If omnibus not found separately, try O key
        if omnibus_p is None:
            omnibus_p = _extract_pval(pvals_dict, ["O", "o"])

        if omnibus_p is None and pvals_dict:
            # Last resort: take the last value (typically omnibus in AllelicSeries)
            last_key = list(pvals_dict.keys())[-1]
            omnibus_p = pvals_dict[last_key]
            logger.debug(
                f"COAST: omnibus p-value not found by standard keys; "
                f"using last key '{last_key}' = {omnibus_p}. "
                f"Available keys: {list(pvals_dict.keys())}"
            )

        return burden_p, skat_p, omnibus_p


def _extract_pval(pvals_dict: dict[str, float | None], keys: list[str]) -> float | None:
    """Extract p-value from dict by trying multiple key names."""
    for key in keys:
        if key in pvals_dict:
            return pvals_dict[key]
    return None
