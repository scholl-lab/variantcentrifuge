# File: variantcentrifuge/gene_burden.py
# Location: variantcentrifuge/variantcentrifuge/gene_burden.py
"""
Gene burden analysis module.

Implements a collapsing burden test (CMC/CAST) for rare variant association.

Two modes are supported:
- "samples" (default): Counts unique carrier samples per gene (binary collapse).
  Each sample is counted once regardless of how many qualifying variants it carries.
  This is the CMC/CAST collapsing test (Li & Leal 2008, Morgenthaler & Thilly 2007).
- "alleles": For each sample, takes the maximum allele dosage (0/1/2) across all
  qualifying variant sites in the gene, then sums across samples. This preserves
  the diploid constraint (total <= 2*N) required for Fisher's exact test.

Statistical testing uses Fisher's exact test on a 2x2 contingency table with
Benjamini-Hochberg FDR or Bonferroni correction for multiple testing.

References
----------
- Li B, Leal SM. Am J Hum Genet. 2008;83(3):311-321 (CMC method)
- Morgenthaler S, Thilly WG. Mutat Res. 2007;615(1-2):28-56 (CAST)
"""

import logging
from math import isnan
from typing import Any

import pandas as pd

from .stages.output_stages import _find_per_sample_gt_columns as _find_gt_columns

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


def _gt_to_dosage(gt: str) -> int:
    """Fast genotype string to allele dosage conversion.

    Maps common genotype strings to alt allele count:
    '1/1' or '1|1' -> 2, '0/1' or '1/0' or '0|1' or '1|0' -> 1, else -> 0.
    """
    if gt in ("1/1", "1|1"):
        return 2
    if gt in ("0/1", "1/0", "0|1", "1|0"):
        return 1
    return 0


def _aggregate_gene_burden_from_columns(
    df: pd.DataFrame,
    case_samples: set[str],
    control_samples: set[str],
    vcf_samples: list[str],
    gt_columns: list[str],
) -> list[dict]:
    """
    Aggregate gene burden using per-sample GT columns (vectorized, Phase 11).

    Instead of parsing the packed GT string, directly reads pre-split per-sample
    GT columns (GEN_0__GT, GEN_1__GT, ...) for ~10-100x faster aggregation.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with GENE column and per-sample GT columns.
    case_samples : set of str
        Set of case sample IDs.
    control_samples : set of str
        Set of control sample IDs.
    vcf_samples : list of str
        Ordered VCF sample names matching GT column indices.
    gt_columns : list of str
        Per-sample GT column names, sorted by sample index.

    Returns
    -------
    list of dict
        Gene-level burden data with carrier counts and allele counts.
    """
    total_cases = len(case_samples)
    total_controls = len(control_samples)

    # Pre-compute sample classification by column index
    n_cols = min(len(gt_columns), len(vcf_samples))
    case_indices = []
    ctrl_indices = []
    for i in range(n_cols):
        name = vcf_samples[i]
        if name in case_samples:
            case_indices.append(i)
        elif name in control_samples:
            ctrl_indices.append(i)

    case_col_names = [gt_columns[i] for i in case_indices]
    ctrl_col_names = [gt_columns[i] for i in ctrl_indices]

    logger.debug(
        f"Column-based aggregation: {len(case_col_names)} case columns, "
        f"{len(ctrl_col_names)} control columns"
    )

    gene_burden_data = []

    for gene, gene_df in df.groupby("GENE", observed=True):
        n_variants = len(gene_df)

        # For case samples: compute per-sample max dosage across all variants
        p_carrier_count = 0
        p_allele_count = 0
        for col in case_col_names:
            # Convert GT strings to dosages for all variants in this gene
            dosages = gene_df[col].map(_gt_to_dosage)
            max_dosage = int(dosages.max())
            if max_dosage > 0:
                p_carrier_count += 1
                p_allele_count += max_dosage

        # For control samples: same approach
        c_carrier_count = 0
        c_allele_count = 0
        for col in ctrl_col_names:
            dosages = gene_df[col].map(_gt_to_dosage)
            max_dosage = int(dosages.max())
            if max_dosage > 0:
                c_carrier_count += 1
                c_allele_count += max_dosage

        gene_burden_data.append(
            {
                "GENE": gene,
                "proband_count": total_cases,
                "control_count": total_controls,
                "proband_carrier_count": p_carrier_count,
                "control_carrier_count": c_carrier_count,
                "proband_allele_count": p_allele_count,
                "control_allele_count": c_allele_count,
                "n_qualifying_variants": n_variants,
            }
        )

    return gene_burden_data


def _aggregate_gene_burden_from_gt(
    df: pd.DataFrame,
    case_samples: set[str],
    control_samples: set[str],
) -> list[dict]:
    """
    Aggregate gene burden counts from packed GT column using per-sample collapsing.

    Fallback for when per-sample GT columns are not available. Parses the packed
    format "Sample1(0/1);Sample2(1/1)" to extract carrier information.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with GENE and GT columns.
    case_samples : set of str
        Set of case sample IDs.
    control_samples : set of str
        Set of control sample IDs.

    Returns
    -------
    list of dict
        Gene-level burden data with carrier counts and allele counts.
    """
    from .helpers import extract_sample_and_genotype, genotype_to_allele_count

    total_cases = len(case_samples)
    total_controls = len(control_samples)

    gene_burden_data = []

    for gene, gene_df in df.groupby("GENE", observed=True):
        # Track per-sample max allele dosage across all variants in this gene.
        # Using max dosage (not sum) ensures each sample contributes at most
        # 2 alleles to the gene-level count, preserving the diploid constraint.
        case_max_dosage: dict[str, int] = {}
        ctrl_max_dosage: dict[str, int] = {}

        for gt_val in gene_df["GT"]:
            if not isinstance(gt_val, str) or not gt_val.strip():
                continue
            for entry in gt_val.split(";"):
                entry = entry.strip()
                if not entry:
                    continue
                sample_name, genotype = extract_sample_and_genotype(entry)
                if not sample_name:
                    continue
                allele_count = genotype_to_allele_count(genotype)
                if allele_count <= 0:
                    continue

                if sample_name in case_samples:
                    case_max_dosage[sample_name] = max(
                        case_max_dosage.get(sample_name, 0), allele_count
                    )
                elif sample_name in control_samples:
                    ctrl_max_dosage[sample_name] = max(
                        ctrl_max_dosage.get(sample_name, 0), allele_count
                    )

        # Carrier counts: unique samples with any qualifying variant
        p_carrier_count = len(case_max_dosage)
        c_carrier_count = len(ctrl_max_dosage)

        # Allele counts: sum of per-sample max dosages (bounded by 2*N)
        p_allele_count = sum(case_max_dosage.values())
        c_allele_count = sum(ctrl_max_dosage.values())

        gene_burden_data.append(
            {
                "GENE": gene,
                "proband_count": total_cases,
                "control_count": total_controls,
                "proband_carrier_count": p_carrier_count,
                "control_carrier_count": c_carrier_count,
                "proband_allele_count": p_allele_count,
                "control_allele_count": c_allele_count,
                "n_qualifying_variants": len(gene_df),
            }
        )

    return gene_burden_data


def perform_gene_burden_analysis(
    df: pd.DataFrame,
    cfg: dict[str, Any],
    case_samples: set[str] | None = None,
    control_samples: set[str] | None = None,
    vcf_samples: list[str] | None = None,
) -> pd.DataFrame:
    """
    Perform gene burden analysis using a collapsing test with Fisher's exact test.

    When case_samples and control_samples are provided, uses proper per-sample
    collapsing (CMC/CAST method) to avoid double-counting samples with variants
    at multiple sites in the same gene.

    Three aggregation strategies (selected automatically by priority):
    1. Column-based (fastest): Uses per-sample GT columns (GEN_0__GT, etc.)
       when vcf_samples is provided and columns exist in the DataFrame.
    2. Packed GT string: Parses "Sample1(0/1);Sample2(1/1)" format.
    3. Legacy: Sums pre-computed per-variant counts (backward compatibility).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with variant data. Must include "GENE" column.
    cfg : dict
        Configuration dictionary with keys:
        - "gene_burden_mode": "samples" (carrier collapse) or "alleles" (max dosage)
        - "correction_method": "fdr" or "bonferroni"
        - "confidence_interval_method": str (optional, default "normal_approx")
        - "confidence_interval_alpha": float (optional, default 0.05)
    case_samples : set of str, optional
        Case sample IDs. When provided with control_samples, enables proper
        per-sample collapsing.
    control_samples : set of str, optional
        Control sample IDs.
    vcf_samples : list of str, optional
        Ordered VCF sample names. When provided with per-sample GT columns,
        enables fast column-based aggregation.

    Returns
    -------
    pd.DataFrame
        Gene-level burden results with p-values, odds ratios, and CIs.
    """
    logger.debug("Starting gene burden aggregation...")

    mode = cfg.get("gene_burden_mode", "samples")
    correction_method = cfg.get("correction_method", "fdr")
    ci_method = cfg.get("confidence_interval_method", "normal_approx")
    ci_alpha = cfg.get("confidence_interval_alpha", 0.05)
    continuity_correction = cfg.get("continuity_correction", 0.5)

    # Determine aggregation strategy (priority: columns > packed GT > legacy)
    has_case_ctrl = case_samples is not None and control_samples is not None
    gt_columns = _find_gt_columns(df) if has_case_ctrl else []
    use_column_aggregation = bool(
        has_case_ctrl and gt_columns and vcf_samples and len(gt_columns) <= len(vcf_samples)
    )
    use_gt_aggregation = has_case_ctrl and not use_column_aggregation and "GT" in df.columns

    if use_column_aggregation:
        assert case_samples is not None
        assert control_samples is not None
        assert vcf_samples is not None
        logger.info(
            f"Using column-based aggregation for gene burden "
            f"({len(gt_columns)} GT columns, "
            f"{len(case_samples)} cases, {len(control_samples)} controls)"
        )
        gene_burden_data = _aggregate_gene_burden_from_columns(
            df, case_samples, control_samples, vcf_samples, gt_columns
        )
    elif use_gt_aggregation:
        assert case_samples is not None
        assert control_samples is not None
        logger.info(
            "Using packed GT string collapsing for gene burden "
            f"({len(case_samples)} cases, {len(control_samples)} controls)"
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
        p_count = int(getattr(row, "proband_count", 0))
        c_count = int(getattr(row, "control_count", 0))

        if p_count == 0 and c_count == 0:
            continue

        if mode == "samples":
            # Collapsing test (CMC/CAST): carrier vs non-carrier
            if uses_per_sample_collapsing:
                p_var = int(getattr(row, "proband_carrier_count", 0))
                c_var = int(getattr(row, "control_carrier_count", 0))
                var_metric = ("proband_carrier_count", "control_carrier_count")
            else:
                p_var = int(getattr(row, "proband_variant_count", 0))
                c_var = int(getattr(row, "control_variant_count", 0))
                var_metric = ("proband_variant_count", "control_variant_count")

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
            p_ref = p_count * 2 - p_all
            c_ref = c_count * 2 - c_all
            table = [[p_all, c_all], [p_ref, c_ref]]
            var_metric = ("proband_allele_count", "control_allele_count")

            if p_ref < 0 or c_ref < 0:
                logger.error(
                    f"Gene {gene} has negative reference allele counts: "
                    f"p_count={p_count}, p_all={p_all}, p_ref={p_ref}, "
                    f"c_count={c_count}, c_all={c_all}, c_ref={c_ref}"
                )
                continue

        if fisher_exact is not None:
            odds_ratio, pval = fisher_exact(table)
        else:
            odds_ratio = float("nan")
            pval = 1.0

        ci_lower, ci_upper = _compute_or_confidence_interval(
            table, odds_ratio, ci_method, ci_alpha, continuity_correction
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
        if correction_method == "bonferroni":
            corrected_pvals = smm.multipletests(pvals, method="bonferroni")[1]
        else:
            corrected_pvals = smm.multipletests(pvals, method="fdr_bh")[1]
    else:
        logger.warning("statsmodels not available. Skipping multiple testing correction.")
        corrected_pvals = pvals

    results_df["corrected_p_value"] = corrected_pvals

    if mode == "samples":
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

    if uses_per_sample_collapsing and "n_qualifying_variants" in results_df.columns:
        final_cols.append("n_qualifying_variants")

    results_df = results_df[final_cols]
    logger.debug("Gene burden analysis complete.")
    return results_df


def _aggregate_gene_burden_legacy(df: pd.DataFrame) -> list[dict]:
    """
    Legacy aggregation: sum per-variant counts by gene.

    This is the fallback when case/control sample sets are not available.
    Note: this method double-counts samples that carry variants at multiple
    sites in the same gene. Use GT-based aggregation when possible.
    """
    gene_burden_data = []
    for gene, gene_df in df.groupby("GENE", observed=True):
        p_count = gene_df["proband_count"].iloc[0]
        c_count = gene_df["control_count"].iloc[0]

        p_var_count = int(gene_df["proband_variant_count"].sum())
        c_var_count = int(gene_df["control_variant_count"].sum())
        p_allele_count = int(gene_df["proband_allele_count"].sum())
        c_allele_count = int(gene_df["control_allele_count"].sum())

        # Cap to prevent impossible values (band-aid for legacy path)
        if p_var_count > p_count:
            p_var_count = p_count
        if c_var_count > c_count:
            c_var_count = c_count
        if p_allele_count > p_count * 2:
            p_allele_count = p_count * 2
        if c_allele_count > c_count * 2:
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

    return gene_burden_data
