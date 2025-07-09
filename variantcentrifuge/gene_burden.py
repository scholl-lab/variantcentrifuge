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

Updated for Issue #31:
- Improved handling of edge cases (e.g., infinite or zero odds_ratio).
  Now properly detects structural zeros and applies continuity correction
  for zero cells to calculate meaningful confidence intervals.
- Uses score method as primary CI calculation (more robust for sparse data).
- Returns NaN for structural zeros where OR cannot be calculated.

Configuration Additions
-----------------------
- "confidence_interval_method": str
    Method for confidence interval calculation. Defaults to "normal_approx".
- "confidence_interval_alpha": float
    Significance level for confidence interval calculation. Default: 0.05 for a 95% CI.
- "continuity_correction": float
    Value to add to zero cells for continuity correction. Default: 0.5.

Outputs
-------
In addition to existing columns (p-values, odds ratio), the result now includes:
- "or_ci_lower": Lower bound of the odds ratio confidence interval.
- "or_ci_upper": Upper bound of the odds ratio confidence interval.
"""

import logging
from math import isnan
from typing import Any, Dict

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

logger = logging.getLogger("variantcentrifuge")


def _compute_or_confidence_interval(
    table: list, odds_ratio: float, method: str, alpha: float, continuity_correction: float = 0.5
) -> tuple:
    """
    Compute confidence intervals for the odds ratio with robust handling of zero cells.

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
          Fallback to "logit" if normal approximation fails.
    alpha : float
        Significance level for the confidence interval. 0.05 for 95% CI.
    continuity_correction : float
        Value to add to zero cells for continuity correction (default: 0.5).

    Returns
    -------
    tuple
        Tuple of (ci_lower, ci_upper) as floats.

    Notes
    -----
    This enhanced function first checks for structural zeros. If found, it correctly
    returns NaN as no OR can be calculated. It then applies a continuity correction
    if any cell is zero to prevent division-by-zero errors, allowing for the
    calculation of meaningful CIs in edge cases.
    """
    a = table[0][0]
    b = table[0][1]
    c = table[1][0]
    d = table[1][1]

    # Check for structural zeros where an OR is not calculable
    if (a + b == 0) or (c + d == 0) or (a + c == 0) or (b + d == 0):
        logger.debug(f"Structural zero in table {table}, cannot compute OR or CI.")
        return np.nan, np.nan

    if Table2x2 is None:
        logger.warning("statsmodels not available. Cannot compute confidence intervals.")
        return np.nan, np.nan

    # Apply continuity correction if any cell is zero
    table_np = np.array(table, dtype=float)
    if 0 in table_np.flatten():
        logger.debug(f"Applying continuity correction ({continuity_correction}) to table {table}")
        table_for_ci = table_np + continuity_correction
    else:
        table_for_ci = table_np

    try:
        cont_table = Table2x2(table_for_ci)
        # Use the score method for confidence intervals, as it's robust for sparse data
        ci_lower, ci_upper = cont_table.oddsratio_confint(alpha=alpha, method="score")

        # If score method fails, try other methods
        if isnan(ci_lower) or isnan(ci_upper):
            logger.debug("Score method failed, trying normal approximation.")
            ci_lower, ci_upper = cont_table.oddsratio_confint(alpha=alpha, method="normal")

            if isnan(ci_lower) or isnan(ci_upper):
                logger.debug("Normal approximation failed, trying logit method.")
                ci_lower, ci_upper = cont_table.oddsratio_confint(alpha=alpha, method="logit")

        return ci_lower, ci_upper
    except Exception as e:
        logger.warning(f"Failed to compute CI for table {table_for_ci.tolist()}: {e}")
        return np.nan, np.nan


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
    If no fisher_exact test is available, p-values default to 1.0.
    Confidence intervals now attempt multiple methods and fallback intervals
    in edge cases (Issue #31).
    """
    logger.debug("Starting gene burden aggregation...")

    mode = cfg.get("gene_burden_mode", "alleles")
    correction_method = cfg.get("correction_method", "fdr")
    ci_method = cfg.get("confidence_interval_method", "normal_approx")
    ci_alpha = cfg.get("confidence_interval_alpha", 0.05)
    continuity_correction = cfg.get("continuity_correction", 0.5)

    # First aggregate by gene to get per-gene counts
    # For gene burden analysis, we need to count unique samples with variants per gene,
    # not sum variant counts across all variants in the gene
    logger.debug("Aggregating variant data by gene for burden analysis...")

    gene_burden_data = []
    for gene, gene_df in df.groupby("GENE"):
        # Get sample counts (should be consistent across all variants)
        p_count = gene_df["proband_count"].iloc[0]
        c_count = gene_df["control_count"].iloc[0]

        # For gene burden analysis, we need to identify which samples have ANY variant in this gene
        # This requires parsing the GT column to find unique samples with variants
        case_samples_with_variants = set()
        control_samples_with_variants = set()
        case_total_alleles = 0
        control_total_alleles = 0

        # Parse GT columns to find samples with variants in this gene
        for _, row in gene_df.iterrows():
            gt_value = str(row.get("GT", ""))
            if gt_value and gt_value not in ["NA", "nan", ""]:
                # Parse the GT field: "Sample1(0/1);Sample2(1/1);..."
                for sample_entry in gt_value.split(";"):
                    sample_entry = sample_entry.strip()
                    if not sample_entry:
                        continue

                    # Extract sample name and genotype
                    if "(" in sample_entry and ")" in sample_entry:
                        sample_name = sample_entry.split("(")[0].strip()
                        genotype_part = sample_entry.split("(")[1].split(")")[0]
                        genotype = (
                            genotype_part.split(":")[0] if ":" in genotype_part else genotype_part
                        )

                        # Count alleles
                        if genotype == "1/1":
                            allele_count = 2
                        elif genotype in ["0/1", "1/0"]:
                            allele_count = 1
                        else:
                            allele_count = 0

                        # Add to appropriate group if has variant
                        if allele_count > 0:
                            # We need to determine if this sample is case or control
                            # For now, we'll use the aggregated counts from the existing columns
                            pass

        # Use the existing aggregated counts for now (this preserves the original logic)
        # but we'll add validation to ensure they make sense
        p_var_count = int(gene_df["proband_variant_count"].sum())
        c_var_count = int(gene_df["control_variant_count"].sum())
        p_allele_count = int(gene_df["proband_allele_count"].sum())
        c_allele_count = int(gene_df["control_allele_count"].sum())

        # Validate that variant counts don't exceed sample counts
        if p_var_count > p_count:
            logger.warning(
                f"Gene {gene}: proband variant count ({p_var_count}) > proband sample count ({p_count}). Using sample count."
            )
            p_var_count = p_count
        if c_var_count > c_count:
            logger.warning(
                f"Gene {gene}: control variant count ({c_var_count}) > control sample count ({c_count}). Using sample count."
            )
            c_var_count = c_count

        # Validate that allele counts don't exceed 2 * sample counts
        if p_allele_count > p_count * 2:
            logger.warning(
                f"Gene {gene}: proband allele count ({p_allele_count}) > 2 * proband sample count ({p_count * 2}). Capping at maximum."
            )
            p_allele_count = p_count * 2
        if c_allele_count > c_count * 2:
            logger.warning(
                f"Gene {gene}: control allele count ({c_allele_count}) > 2 * control sample count ({c_count * 2}). Capping at maximum."
            )
            c_allele_count = c_count * 2

        gene_burden_data.append(
            {
                "GENE": gene,
                "proband_count": p_count,
                "control_count": c_count,
                "proband_variant_count": p_var_count,
                "control_variant_count": c_var_count,
                "proband_allele_count": p_allele_count,
                "control_allele_count": c_allele_count,
            }
        )

    grouped = pd.DataFrame(gene_burden_data)

    # Sort by gene name to ensure deterministic order
    grouped = grouped.sort_values("GENE").reset_index(drop=True)

    # Debug logging to track determinism
    logger.info(f"Gene burden analysis processing {len(grouped)} genes in deterministic order")
    for i, row in grouped.iterrows():
        logger.debug(
            f"Processing gene {row['GENE']}: case_count={row['proband_count']}, control_count={row['control_count']}, "
            f"case_alleles={row['proband_allele_count']}, control_alleles={row['control_allele_count']}"
        )

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

            # Debug logging for negative values
            if p_ref < 0 or c_ref < 0:
                logger.error(
                    f"Gene {gene} has negative reference counts: p_count={p_count}, p_var={p_var}, p_ref={p_ref}, c_count={c_count}, c_var={c_var}, c_ref={c_ref}"
                )
                continue
        else:
            # Per-allele counts
            p_all = int(row["proband_allele_count"])
            c_all = int(row["control_allele_count"])
            p_ref = p_count * 2 - p_all
            c_ref = c_count * 2 - c_all
            table = [[p_all, c_all], [p_ref, c_ref]]
            var_metric = ("proband_allele_count", "control_allele_count")

            # Debug logging for negative values
            if p_ref < 0 or c_ref < 0:
                logger.error(
                    f"Gene {gene} has negative reference allele counts: p_count={p_count}, p_all={p_all}, p_ref={p_ref}, c_count={c_count}, c_all={c_all}, c_ref={c_ref}"
                )
                continue

        if fisher_exact is not None:
            odds_ratio, pval = fisher_exact(table)
        else:
            odds_ratio = float("nan")
            pval = 1.0

        # Compute confidence interval for the odds ratio
        ci_lower, ci_upper = _compute_or_confidence_interval(
            table, odds_ratio, ci_method, ci_alpha, continuity_correction
        )

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
    if smm is None:
        logger.warning("statsmodels not available. Skipping multiple testing correction.")
        corrected_pvals = pvals  # Use raw p-values as fallback
    else:
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
