# File: variantcentrifuge/gene_burden.py
# Location: variantcentrifuge/variantcentrifuge/gene_burden.py
"""
Gene burden analysis module.

Provides:
- perform_gene_burden_analysis: Aggregates per-gene counts (samples or alleles),
  performs Fisher's exact test, calculates confidence intervals, and applies multiple
  testing correction.

New Features (Issue #21):
- Adds confidence intervals to the gene burden metrics (e.g., odds ratio).
- Confidence interval calculation method and confidence level can be configured.

Configuration Additions
-----------------------
- "confidence_interval_method": str
    Method for confidence interval calculation. Currently supports:
    - "normal_approx": Uses statsmodels' Table2x2 normal approximation for OR CI.
    Future methods (e.g., "wilson") can be added.
- "confidence_interval_alpha": float
    Significance level for confidence interval calculation. Default: 0.05 for a 95% CI.

Outputs
-------
In addition to existing columns (p-values, odds ratio), the result now includes:
- "or_ci_lower": Lower bound of the odds ratio confidence interval.
- "or_ci_upper": Upper bound of the odds ratio confidence interval.
"""

import logging
from math import isnan
from typing import Dict, Any
import pandas as pd

try:
    from scipy.stats import fisher_exact
except ImportError:
    fisher_exact = None

import statsmodels.stats.multitest as smm
from statsmodels.stats.contingency_tables import Table2x2
import numpy as np

logger = logging.getLogger("variantcentrifuge")


def _compute_or_confidence_interval(table: list, odds_ratio: float, method: str, alpha: float) -> (float, float):
    """
    Compute confidence intervals for the odds ratio using the specified method.

    Parameters
    ----------
    table : list
        2x2 contingency table, e.g. [[a, b], [c, d]].
    odds_ratio : float
        The odds ratio for the given table.
    method : str
        Method for confidence interval calculation.
        Currently supported:
        - "normal_approx": Uses statsmodels Table2x2 normal approximation.
    alpha : float
        Significance level for the confidence interval. 0.05 for 95% CI.

    Returns
    -------
    (float, float)
        Tuple of (ci_lower, ci_upper).

    Notes
    -----
    If odds_ratio is NaN, zero, or infinite, the confidence interval cannot be computed
    and NaN is returned for both bounds.

    normal_approx method (via statsmodels):
    Uses statsmodels.stats.contingency_tables.Table2x2.oddsratio_confint to compute
    the confidence interval. If any cell in the table is zero, or if the OR is invalid,
    NaNs are returned.
    """
    if isnan(odds_ratio) or odds_ratio <= 0 or np.isinf(odds_ratio):
        # Cannot compute a CI if odds_ratio <= 0, NaN, or infinite
        return float('nan'), float('nan')

    a = table[0][0]
    b = table[0][1]
    c = table[1][0]
    d = table[1][1]

    # If any cell is zero, the normal approximation may not be valid
    if a == 0 or b == 0 or c == 0 or d == 0:
        return float('nan'), float('nan')

    # Use statsmodels Table2x2 for CI computation
    if method == "normal_approx":
        try:
            table_np = np.array(table)
            cont_table = Table2x2(table_np)
            ci_lower, ci_upper = cont_table.oddsratio_confint(alpha=alpha, method="normal")
        except Exception:
            # If any unexpected error occurs, fallback to NaN
            ci_lower, ci_upper = float('nan'), float('nan')
    else:
        # Unsupported methods return NaN by default
        ci_lower, ci_upper = float('nan'), float('nan')

    return ci_lower, ci_upper


def perform_gene_burden_analysis(df: pd.DataFrame, cfg: Dict[str, Any]) -> pd.DataFrame:
    """
    Perform gene burden analysis for each gene.

    Steps
    -----
    1. Aggregate variant counts per gene based on the chosen mode (samples or alleles).
    2. Perform Fisher's exact test for each gene.
    3. Apply multiple testing correction (FDR or Bonferroni).
    4. Compute and add confidence intervals for the odds ratio.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with per-variant case/control counts. Must include columns:
        "GENE", "proband_count", "control_count", "proband_variant_count",
        "control_variant_count", "proband_allele_count", "control_allele_count".
    cfg : dict
        Configuration dictionary with keys:
        - "gene_burden_mode": str
            "samples" or "alleles" indicating the aggregation mode.
        - "correction_method": str
            "fdr" or "bonferroni" for multiple testing correction.
        - "confidence_interval_method": str (optional)
            Method for CI calculation ("normal_approx"), defaults to "normal_approx".
        - "confidence_interval_alpha": float (optional)
            Significance level for CI, defaults to 0.05 for a 95% CI.

    Returns
    -------
    pd.DataFrame
        A DataFrame with gene-level burden results, including p-values, odds ratios,
        and confidence intervals.

        The output includes:
        - "GENE"
        - "proband_count", "control_count"
        - Either "proband_variant_count", "control_variant_count" or
          "proband_allele_count", "control_allele_count" depending on mode
        - "raw_p_value"
        - "corrected_p_value"
        - "odds_ratio"
        - "or_ci_lower"
        - "or_ci_upper"

    Notes
    -----
    If no fisher_exact test is available (scipy not installed), p-values default to 1.0.
    Confidence intervals are computed using statsmodels Table2x2 normal approximation
    if "normal_approx" is chosen.
    """

    logger.debug("Starting gene burden aggregation...")

    mode = cfg.get("gene_burden_mode", "alleles")
    correction_method = cfg.get("correction_method", "fdr")
    ci_method = cfg.get("confidence_interval_method", "normal_approx")
    ci_alpha = cfg.get("confidence_interval_alpha", 0.05)

    grouped = df.groupby("GENE").agg(
        {
            "proband_count": "max",
            "control_count": "max",
            "proband_variant_count": "sum",
            "control_variant_count": "sum",
            "proband_allele_count": "sum",
            "control_allele_count": "sum",
        }
    ).reset_index()

    results = []
    for _, row in grouped.iterrows():
        gene = row["GENE"]
        p_count = int(row["proband_count"])
        c_count = int(row["control_count"])

        # Skip if no samples in either group
        if p_count == 0 and c_count == 0:
            continue

        if mode == "samples":
            # Per-sample variant counts
            p_var = int(row["proband_variant_count"])
            c_var = int(row["control_variant_count"])
            p_ref = p_count - p_var
            c_ref = c_count - c_var
            table = [[p_var, c_var], [p_ref, c_ref]]
            var_metric = ("proband_variant_count", "control_variant_count")
        else:
            # Per-allele counts
            p_all = int(row["proband_allele_count"])
            c_all = int(row["control_allele_count"])
            p_ref = p_count * 2 - p_all
            c_ref = c_count * 2 - c_all
            table = [[p_all, c_all], [p_ref, c_ref]]
            var_metric = ("proband_allele_count", "control_allele_count")

        if fisher_exact is not None:
            odds_ratio, pval = fisher_exact(table)
        else:
            odds_ratio = float('nan')
            pval = 1.0

        # Compute confidence interval for the odds ratio
        if isnan(odds_ratio) or odds_ratio is None:
            ci_lower, ci_upper = float('nan'), float('nan')
        else:
            ci_lower, ci_upper = _compute_or_confidence_interval(table, odds_ratio, ci_method, ci_alpha)

        results.append(
            {
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
        )

    if not results:
        logger.warning("No genes found with variant data for gene burden analysis.")
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # Multiple testing correction
    pvals = results_df["raw_p_value"].values
    if correction_method == "bonferroni":
        corrected_pvals = smm.multipletests(pvals, method="bonferroni")[1]
    else:
        # Default to FDR (Benjamini–Hochberg)
        corrected_pvals = smm.multipletests(pvals, method="fdr_bh")[1]

    results_df["corrected_p_value"] = corrected_pvals

    if mode == "samples":
        final_cols = [
            "GENE",
            "proband_count",
            "control_count",
            "proband_variant_count",
            "control_variant_count",
            "raw_p_value",
            "corrected_p_value",
            "odds_ratio",
            "or_ci_lower",
            "or_ci_upper",
        ]
    else:
        final_cols = [
            "GENE",
            "proband_count",
            "control_count",
            "proband_allele_count",
            "control_allele_count",
            "raw_p_value",
            "corrected_p_value",
            "odds_ratio",
            "or_ci_lower",
            "or_ci_upper",
        ]

    results_df = results_df[final_cols]
    logger.debug("Gene burden analysis complete.")
    return results_df
