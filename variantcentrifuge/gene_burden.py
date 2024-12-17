# File: variantcentrifuge/gene_burden.py
# Location: variantcentrifuge/variantcentrifuge/gene_burden.py

"""
Gene burden analysis module.

Provides:
- perform_gene_burden_analysis: Aggregates per-gene counts (samples or alleles),
  performs Fisher's exact test, and applies multiple testing correction.
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

logger = logging.getLogger("variantcentrifuge")


def perform_gene_burden_analysis(df: pd.DataFrame, cfg: Dict[str, Any]) -> pd.DataFrame:
    """
    Perform gene burden analysis for each gene.

    Steps
    -----
    1. Aggregate variant counts per gene based on the chosen mode (samples or alleles).
    2. Perform Fisher's exact test for each gene.
    3. Apply multiple testing correction (FDR or Bonferroni).

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

    Returns
    -------
    pd.DataFrame
        A DataFrame with gene-level burden results, including p-values and odds ratios.
        The output includes:
        - "GENE"
        - "proband_count", "control_count"
        - Either "proband_variant_count", "control_variant_count" or
          "proband_allele_count", "control_allele_count"
        - "raw_p_value"
        - "corrected_p_value"
        - "odds_ratio"

    Notes
    -----
    If no fisher_exact test is available (scipy not installed), p-values default to 1.0.
    """
    logger.debug("Starting gene burden aggregation...")

    mode = cfg.get("gene_burden_mode", "alleles")
    correction_method = cfg.get("correction_method", "fdr")

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
            odds_ratio = None
            pval = 1.0

        results.append(
            {
                "GENE": gene,
                "proband_count": p_count,
                "control_count": c_count,
                var_metric[0]: table[0][0],
                var_metric[1]: table[0][1],
                "raw_p_value": pval,
                "odds_ratio": odds_ratio,
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
        ]

    results_df = results_df[final_cols]
    logger.debug("Gene burden analysis complete.")
    return results_df
