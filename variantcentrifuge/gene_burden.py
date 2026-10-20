"""
Gene burden analysis module (CMC/CAST collapsing test for rare variant association).

Supports 'samples' (binary collapse) and 'alleles' (max dosage collapse) modes
with Fisher's exact test and multiple testing correction (FDR/Bonferroni).
"""

import logging
import re
import warnings
from typing import Any

import pandas as pd

try:
    from scipy.stats import fisher_exact
except ImportError:
    fisher_exact = None

import numpy as np

try:
    import statsmodels.stats.multitest as smm
    from statsmodels.stats.contingency_tables import Table2x2
except ImportError:
    smm = None
    Table2x2 = None

try:
    from .association.correction import apply_correction as _apply_correction
except ImportError:
    _apply_correction = None  # type: ignore[assignment]

from .association.aggregation import (  # noqa: F401 - re-exported for existing callers
    _aggregate_gene_burden_from_columns,
    _aggregate_gene_burden_from_gt,
    _aggregate_gene_burden_legacy,
    _called_and_chrom,
    _col_total,
    _gt_to_dosage,
    _legacy_counts,
)
from .gt_reconstruction import gt_columns_match_samples

logger = logging.getLogger("variantcentrifuge")


def warn_legacy_gene_burden(missing_as_ref: bool = False) -> None:
    """Warn that --perform-gene-burden is deprecated in favour of the association engine."""
    msg = (
        "--perform-gene-burden is deprecated and will be removed in the next minor release; "
        "use --perform-association (e.g. --association-tests fisher) instead."
    )
    if missing_as_ref:
        msg += " --missing-as-ref is ignored by the legacy gene burden path."
    logger.warning(msg)
    warnings.warn(msg, DeprecationWarning, stacklevel=2)


_VARIANT_KEY = ["GENE", "CHROM", "POS", "REF", "ALT"]


def drop_transcript_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Keep one row per variant per gene before burden/association tests.

    ``--split-snpeff-lines`` writes one row per transcript, so a variant on two
    transcripts of a gene would enter the genotype matrix twice. snpEff orders
    annotations by severity, so the first row kept is the most severe consequence.
    """
    if df.empty or not set(_VARIANT_KEY).issubset(df.columns):
        return df
    deduped = df.drop_duplicates(subset=_VARIANT_KEY)
    if len(deduped) < len(df):
        logger.info(
            "Collapsed %d transcript-duplicate rows to one row per variant per gene",
            len(df) - len(deduped),
        )
    return deduped


def _find_gt_columns(df: pd.DataFrame) -> list[str]:
    """Find per-sample GT columns in DataFrame, sorted by sample index."""
    cols = [c for c in df.columns if re.match(r"^GEN[_\[](\d+)[_\]\.]+(GT)$", c)]
    return sorted(cols, key=lambda c: int(re.search(r"(\d+)", c).group(1)))  # type: ignore[union-attr]


def _compute_or_confidence_interval(
    table: list, alpha: float, continuity_correction: float = 0.5
) -> tuple[float, float, float]:
    """Return (odds ratio, CI lower, CI upper); Haldane-corrected if any cell is zero."""
    a, b = table[0][0], table[0][1]
    c, d = table[1][0], table[1][1]

    if (a + b == 0) or (c + d == 0) or (a + c == 0) or (b + d == 0) or Table2x2 is None:
        return np.nan, np.nan, np.nan

    table_np = np.array(table, dtype=float)
    table_for_ci = table_np + continuity_correction if 0 in table_np.flatten() else table_np

    try:
        cont_table = Table2x2(table_for_ci)
        ci_lower, ci_upper = cont_table.oddsratio_confint(alpha=alpha, method="normal")
        return float(cont_table.oddsratio), float(ci_lower), float(ci_upper)
    except Exception as e:
        logger.warning(f"Failed to compute OR/CI for table {table_for_ci.tolist()}: {e}")
        return np.nan, np.nan, np.nan


def perform_gene_burden_analysis(
    df: pd.DataFrame,
    cfg: dict[str, Any],
    case_samples: set[str] | None = None,
    control_samples: set[str] | None = None,
    vcf_samples: list[str] | None = None,
) -> pd.DataFrame:
    """
    Perform gene burden analysis using collapsing (CMC/CAST) with Fisher's exact test.

    Supports column-based, packed GT, and legacy per-variant count aggregation.
    """
    logger.debug("Starting gene burden aggregation...")

    mode = cfg.get("gene_burden_mode", "samples")
    correction_method = cfg.get("correction_method", "fdr")
    ci_alpha = cfg.get("confidence_interval_alpha", 0.05)
    continuity_correction = cfg.get("continuity_correction", 0.5)

    attrs = getattr(df, "attrs", {})
    if case_samples is None:
        p_s = cfg.get("proband_samples") or cfg.get("case_samples") or attrs.get("proband_samples")
        if p_s is not None:
            case_samples = set(p_s)
    if control_samples is None:
        c_s = cfg.get("control_samples") or attrs.get("control_samples")
        if c_s is not None:
            control_samples = set(c_s)

    df = drop_transcript_duplicates(df)

    # Determine aggregation strategy (priority: columns > packed GT > legacy)
    has_case_ctrl = case_samples is not None and control_samples is not None
    gt_columns = _find_gt_columns(df) if has_case_ctrl else []
    # Columns only when they are GEN[0..n-1].GT for the n VCF samples (#219).
    use_column_aggregation = has_case_ctrl and gt_columns_match_samples(gt_columns, vcf_samples)
    use_gt_aggregation = has_case_ctrl and not use_column_aggregation and "GT" in df.columns
    if gt_columns and use_gt_aggregation:
        logger.warning(
            "%d per-sample GT columns do not map 1:1 to the %d VCF samples; using packed GT",
            len(gt_columns),
            len(vcf_samples or []),
        )

    if use_column_aggregation and case_samples and control_samples and vcf_samples:
        logger.info(
            "Using column-based aggregation for gene burden (%d GT cols, %d cases, %d ctrls)",
            len(gt_columns),
            len(case_samples),
            len(control_samples),
        )
        gene_burden_data = _aggregate_gene_burden_from_columns(
            df, case_samples, control_samples, vcf_samples, gt_columns
        )
    elif use_gt_aggregation and case_samples and control_samples:
        logger.info(
            "Using packed GT string collapsing for gene burden (%d cases, %d ctrls)",
            len(case_samples),
            len(control_samples),
        )
        gene_burden_data = _aggregate_gene_burden_from_gt(df, case_samples, control_samples)
    else:
        logger.info("Using pre-computed per-variant counts for gene burden (legacy mode)")
        gene_burden_data = _aggregate_gene_burden_legacy(df)

    if not gene_burden_data:
        logger.warning("No genes found with variant data for gene burden analysis.")
        return pd.DataFrame()

    grouped = pd.DataFrame(gene_burden_data)
    grouped = grouped.sort_values("GENE").reset_index(drop=True)

    logger.info(f"Gene burden analysis processing {len(grouped)} genes in deterministic order")

    # Both column-based and GT-based aggregation produce carrier counts
    uses_per_sample_collapsing = use_column_aggregation or use_gt_aggregation

    # Build 2x2 tables and run Fisher's exact test
    results = []
    for row in grouped.itertuples(index=False):
        gene = getattr(row, "GENE", "")
        p_c = getattr(row, "proband_called_count", None)
        p_count = (
            int(p_c)
            if p_c is not None and not pd.isna(p_c)
            else int(getattr(row, "proband_count", 0))
        )
        c_c = getattr(row, "control_called_count", None)
        c_count = (
            int(c_c)
            if c_c is not None and not pd.isna(c_c)
            else int(getattr(row, "control_count", 0))
        )

        if p_count == 0 or c_count == 0:
            continue

        if mode == "samples":
            suffix = "carrier_count" if uses_per_sample_collapsing else "variant_count"
            p_key, c_key = f"proband_{suffix}", f"control_{suffix}"
            p_var, c_var = int(getattr(row, p_key, 0)), int(getattr(row, c_key, 0))
            var_metric = (p_key, c_key)

            p_ref = p_count - p_var
            c_ref = c_count - c_var
            table = [[p_var, c_var], [p_ref, c_ref]]

            if p_ref < 0 or c_ref < 0:
                logger.error(
                    f"Gene {gene} has negative reference counts: p_count={p_count}, "
                    f"p_var={p_var}, p_ref={p_ref}, c_count={c_count}, "
                    f"c_var={c_var}, c_ref={c_ref}"
                )
                continue
        else:
            # Allele-based test: alt alleles vs ref alleles
            p_all = int(getattr(row, "proband_allele_count", 0))
            c_all = int(getattr(row, "control_allele_count", 0))
            p_tot_chrom = int(getattr(row, "proband_called_chromosomes", p_count * 2))
            c_tot_chrom = int(getattr(row, "control_called_chromosomes", c_count * 2))
            p_ref = p_tot_chrom - p_all
            c_ref = c_tot_chrom - c_all
            table = [[p_all, c_all], [p_ref, c_ref]]
            var_metric = ("proband_allele_count", "control_allele_count")

            if p_ref < 0 or c_ref < 0:
                logger.error(
                    f"Gene {gene} has negative reference allele counts: "
                    f"p_count={p_count}, p_all={p_all}, p_ref={p_ref}, "
                    f"c_count={c_count}, c_all={c_all}, c_ref={c_ref}"
                )
                continue

        pval = fisher_exact(table)[1] if fisher_exact is not None else 1.0
        odds_ratio, ci_lower, ci_upper = _compute_or_confidence_interval(
            table, ci_alpha, continuity_correction
        )

        n_vars = int(getattr(row, "n_qualifying_variants", 0)) if uses_per_sample_collapsing else 0

        result_row = {
            "GENE": gene,
            "proband_count": p_count,
            "control_count": c_count,
            var_metric[0]: table[0][0],
            var_metric[1]: table[0][1],
            "raw_p_value": pval,
            "odds_ratio": odds_ratio,
            "or_ci_lower": ci_lower,
            "or_ci_upper": ci_upper,
        }
        if uses_per_sample_collapsing:
            result_row["n_qualifying_variants"] = n_vars

        results.append(result_row)

    if not results:
        logger.warning("No genes found with variant data for gene burden analysis.")
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # Multiple testing correction
    pvals: np.ndarray[Any, Any] = np.asarray(results_df["raw_p_value"].values, dtype=float)
    if _apply_correction is not None:
        corrected_pvals = _apply_correction(pvals, correction_method)
    elif smm is not None:
        method = "bonferroni" if correction_method == "bonferroni" else "fdr_bh"
        corrected_pvals = smm.multipletests(pvals, method=method)[1]
    else:
        logger.warning("statsmodels not available. Skipping multiple testing correction.")
        corrected_pvals = pvals

    results_df["corrected_p_value"] = corrected_pvals

    final_cols = [
        "GENE",
        "proband_count",
        "control_count",
        var_metric[0],
        var_metric[1],
        "raw_p_value",
        "corrected_p_value",
        "odds_ratio",
        "or_ci_lower",
        "or_ci_upper",
    ]

    if uses_per_sample_collapsing and "n_qualifying_variants" in results_df.columns:
        final_cols.append("n_qualifying_variants")

    results_df = results_df[final_cols]
    logger.debug("Gene burden analysis complete.")
    return results_df
