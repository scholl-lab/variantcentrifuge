"""
Gene burden analysis module (CMC/CAST collapsing test for rare variant association).

Supports 'samples' (binary collapse) and 'alleles' (max dosage collapse) modes
with Fisher's exact test and multiple testing correction (FDR/Bonferroni).
"""

import logging
import re
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

from .genotype_utils import is_missing

logger = logging.getLogger("variantcentrifuge")


def _find_gt_columns(df: pd.DataFrame) -> list[str]:
    """Find per-sample GT columns in DataFrame, sorted by sample index."""
    cols = [c for c in df.columns if re.match(r"^GEN[_\[](\d+)[_\]\.]+(GT)$", c)]
    return sorted(cols, key=lambda c: int(re.search(r"(\d+)", c).group(1)))  # type: ignore[union-attr]


def _compute_or_confidence_interval(
    table: list, odds_ratio: float, method: str, alpha: float, continuity_correction: float = 0.5
) -> tuple:
    """Compute confidence intervals for the odds ratio with robust handling of zero cells."""
    a, b = table[0][0], table[0][1]
    c, d = table[1][0], table[1][1]

    if (a + b == 0) or (c + d == 0) or (a + c == 0) or (b + d == 0) or Table2x2 is None:
        return np.nan, np.nan

    table_np = np.array(table, dtype=float)
    table_for_ci = table_np + continuity_correction if 0 in table_np.flatten() else table_np

    try:
        cont_table = Table2x2(table_for_ci)
        ci_lower, ci_upper = cont_table.oddsratio_confint(alpha=alpha, method="normal")
        return ci_lower, ci_upper
    except Exception as e:
        logger.warning(f"Failed to compute CI for table {table_for_ci.tolist()}: {e}")
        return np.nan, np.nan


def _gt_to_dosage(gt: str) -> int:
    """Fast genotype string to allele dosage conversion."""
    if gt in ("1/1", "1|1"):
        return 2
    return 1 if gt in ("0/1", "1/0", "0|1", "1|0", "1") else 0


def _aggregate_gene_burden_from_columns(
    df: pd.DataFrame,
    case_samples: set[str],
    control_samples: set[str],
    vcf_samples: list[str],
    gt_columns: list[str],
) -> list[dict]:
    """Aggregate gene burden using per-sample GT columns (GEN_0__GT, GEN_1__GT, ...)."""
    n_cols = min(len(gt_columns), len(vcf_samples))
    case_col_names = [gt_columns[i] for i in range(n_cols) if vcf_samples[i] in case_samples]
    ctrl_col_names = [gt_columns[i] for i in range(n_cols) if vcf_samples[i] in control_samples]
    gene_burden_data = []

    def _tally(col_names: list[str], gdf: pd.DataFrame) -> tuple[int, int, int, int]:
        carrier, allele, called, called_chrom = 0, 0, 0, 0
        for col in col_names:
            vals = gdf[col]
            if any(
                v is None or (not isinstance(v, str) and pd.isna(v)) or is_missing(str(v).strip())
                for v in vals
            ):
                continue
            called += 1
            s_ploidy = 1 if all("/" not in str(v) and "|" not in str(v) for v in vals) else 2
            called_chrom += s_ploidy
            max_d = int(vals.astype(str).map(_gt_to_dosage).max())
            if max_d > 0:
                carrier += 1
                allele += max_d
        return carrier, allele, called, called_chrom

    for gene, gene_df in df.groupby("GENE", observed=True):
        p_carrier, p_allele, p_called, p_chrom = _tally(case_col_names, gene_df)
        c_carrier, c_allele, c_called, c_chrom = _tally(ctrl_col_names, gene_df)

        gene_burden_data.append(
            {
                "GENE": gene,
                "proband_count": len(case_col_names),
                "control_count": len(ctrl_col_names),
                "proband_called_count": p_called,
                "control_called_count": c_called,
                "proband_called_chromosomes": p_chrom,
                "control_called_chromosomes": c_chrom,
                "proband_carrier_count": p_carrier,
                "control_carrier_count": c_carrier,
                "proband_allele_count": p_allele,
                "control_allele_count": c_allele,
                "n_qualifying_variants": len(gene_df),
            }
        )

    return gene_burden_data


def _aggregate_gene_burden_from_gt(
    df: pd.DataFrame,
    case_samples: set[str],
    control_samples: set[str],
) -> list[dict]:
    """Aggregate gene burden counts from packed GT column using per-sample collapsing."""
    from .helpers import extract_sample_and_genotype, genotype_to_allele_count

    gene_burden_data = []

    for gene, gene_df in df.groupby("GENE", observed=True):
        case_max_dosage: dict[str, int] = {}
        ctrl_max_dosage: dict[str, int] = {}
        missing_cases: set[str] = set()
        missing_ctrls: set[str] = set()
        case_ploidy: dict[str, int] = {}
        ctrl_ploidy: dict[str, int] = {}
        case_seen_rows: dict[str, int] = {}
        ctrl_seen_rows: dict[str, int] = {}

        for gt_val in gene_df["GT"]:
            seen_in_row: set[str] = set()
            if isinstance(gt_val, str) and gt_val.strip():
                for entry in gt_val.split(";"):
                    entry = entry.strip()
                    if not entry:
                        continue
                    sample_name, genotype = extract_sample_and_genotype(entry)
                    if not sample_name:
                        continue
                    seen_in_row.add(sample_name)
                    if is_missing(genotype):
                        if sample_name in case_samples:
                            missing_cases.add(sample_name)
                        elif sample_name in control_samples:
                            missing_ctrls.add(sample_name)
                        continue

                    pl = 2 if ("/" in genotype or "|" in genotype) else 1
                    if sample_name in case_samples:
                        case_ploidy[sample_name] = max(case_ploidy.get(sample_name, 0), pl)
                    elif sample_name in control_samples:
                        ctrl_ploidy[sample_name] = max(ctrl_ploidy.get(sample_name, 0), pl)

                    allele_count = genotype_to_allele_count(genotype)
                    if allele_count > 0:
                        if sample_name in case_samples:
                            cur_d = case_max_dosage.get(sample_name, 0)
                            case_max_dosage[sample_name] = max(cur_d, allele_count)
                        elif sample_name in control_samples:
                            cur_d = ctrl_max_dosage.get(sample_name, 0)
                            ctrl_max_dosage[sample_name] = max(cur_d, allele_count)

            for s in seen_in_row:
                if s in case_samples:
                    case_seen_rows[s] = case_seen_rows.get(s, 0) + 1
                elif s in control_samples:
                    ctrl_seen_rows[s] = ctrl_seen_rows.get(s, 0) + 1

        for s in missing_cases:
            case_max_dosage.pop(s, None)
        for s in missing_ctrls:
            ctrl_max_dosage.pop(s, None)

        n_rows = len(gene_df)
        haploid_cases = sum(
            1
            for s, cnt in case_seen_rows.items()
            if cnt == n_rows and case_ploidy.get(s, 2) == 1 and s not in missing_cases
        )
        haploid_ctrls = sum(
            1
            for s, cnt in ctrl_seen_rows.items()
            if cnt == n_rows and ctrl_ploidy.get(s, 2) == 1 and s not in missing_ctrls
        )
        p_called = len(case_samples - missing_cases)
        c_called = len(control_samples - missing_ctrls)
        p_called_chrom = 2 * p_called - haploid_cases
        c_called_chrom = 2 * c_called - haploid_ctrls
        p_carrier_count = len(case_max_dosage)
        c_carrier_count = len(ctrl_max_dosage)
        p_allele_count = sum(case_max_dosage.values())
        c_allele_count = sum(ctrl_max_dosage.values())

        gene_burden_data.append(
            {
                "GENE": gene,
                "proband_count": len(case_samples),
                "control_count": len(control_samples),
                "proband_called_count": p_called,
                "control_called_count": c_called,
                "proband_called_chromosomes": p_called_chrom,
                "control_called_chromosomes": c_called_chrom,
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
    Perform gene burden analysis using collapsing (CMC/CAST) with Fisher's exact test.

    Supports column-based, packed GT, and legacy per-variant count aggregation.
    """
    logger.debug("Starting gene burden aggregation...")

    mode = cfg.get("gene_burden_mode", "samples")
    correction_method = cfg.get("correction_method", "fdr")
    ci_method = cfg.get("confidence_interval_method", "normal_approx")
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

    # Determine aggregation strategy (priority: columns > packed GT > legacy)
    has_case_ctrl = case_samples is not None and control_samples is not None
    gt_columns = _find_gt_columns(df) if has_case_ctrl else []
    use_column_aggregation = bool(
        has_case_ctrl and gt_columns and vcf_samples and len(gt_columns) <= len(vcf_samples)
    )
    use_gt_aggregation = has_case_ctrl and not use_column_aggregation and "GT" in df.columns

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

        if fisher_exact is not None:
            odds_ratio, pval = fisher_exact(table)
        else:
            odds_ratio = float("nan")
            pval = 1.0

        has_zero_cell = any(cell == 0 for row in table for cell in row)
        if has_zero_cell and Table2x2 is not None:
            a, b, c, d = table[0][0], table[0][1], table[1][0], table[1][1]
            if (a + b == 0) or (c + d == 0) or (a + c == 0) or (b + d == 0):
                odds_ratio, ci_lower, ci_upper = np.nan, np.nan, np.nan
            else:
                table_for_or_ci = np.array(table, dtype=float) + continuity_correction
                try:
                    cont = Table2x2(table_for_or_ci)
                    odds_ratio = float(cont.oddsratio)
                    ci_lower, ci_upper = cont.oddsratio_confint(alpha=ci_alpha, method="normal")
                except Exception as exc:
                    logger.warning("Failed to compute OR/CI: %s", exc)
                    odds_ratio, ci_lower, ci_upper = np.nan, np.nan, np.nan
        else:
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


def _col_total(gdf: pd.DataFrame, col: str, fallback_col: str) -> int:
    if col in gdf.columns:
        return int(gdf[col].iloc[0])
    return int(gdf[fallback_col].max()) if fallback_col in gdf.columns else 0


def _parse_carrier_token(token: str) -> tuple[str, int, int]:
    parts = token.strip().rsplit(":", 2)
    if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
        return parts[0], int(parts[1]), int(parts[2])
    if len(parts) >= 2 and parts[-1].isdigit():
        return token.strip().rsplit(":", 1)[0], int(parts[-1]), 2
    return token.strip(), 1, 2


def _called_and_chrom(gdf: pd.DataFrame, tot: int, prefix: str) -> tuple[int, int]:
    """Compute callable sample count and callable chromosomes for legacy aggregation."""
    m_col, c_col = f"{prefix}_missing_samples", f"{prefix}_called_count"
    chrom_col, car_col = f"{prefix}_called_chromosomes", f"{prefix}_carrier_samples"

    if car_col not in gdf.columns and len(gdf) > 1:
        has_mis = (c_col in gdf.columns and int((tot - gdf[c_col] > 0).sum()) > 0) or (
            m_col in gdf.columns
            and any(
                x.strip() for s in gdf[m_col].dropna() if s for x in str(s).split(";") if x.strip()
            )
        )
        if has_mis:
            return 0, 0

    if m_col in gdf.columns:
        missing = {
            x.strip() for s in gdf[m_col].dropna() if s for x in str(s).split(";") if x.strip()
        }
        called = max(0, tot - len(missing))
    elif c_col in gdf.columns:
        called = (
            int(gdf[c_col].min())
            if (len(gdf) == 1 or int((tot - gdf[c_col] > 0).sum()) == 0)
            else 0
        )
    else:
        called = tot

    if called == 0:
        chrom = 0
    elif car_col in gdf.columns and m_col in gdf.columns:
        haploid_rows: dict[str, int] = {}
        for s_str in gdf[car_col].dropna():
            if not s_str or not str(s_str).strip():
                continue
            seen_in_row: set[str] = set()
            for item in str(s_str).split(";"):
                if not item.strip():
                    continue
                sname, _, pl = _parse_carrier_token(item)
                if sname not in missing and pl == 1 and sname not in seen_in_row:
                    seen_in_row.add(sname)
                    haploid_rows[sname] = haploid_rows.get(sname, 0) + 1
        n_all_haploid = sum(1 for cnt in haploid_rows.values() if cnt == len(gdf))
        chrom = 2 * called - n_all_haploid
    elif chrom_col in gdf.columns:
        chrom = max(0, min(called * 2, int(gdf[chrom_col].min())))
    else:
        chrom = called * 2
    return called, chrom


def _legacy_counts(gdf: pd.DataFrame, called: int, chrom: int, prefix: str) -> tuple[int, int]:
    if called == 0:
        return 0, 0
    m_col, car_col = f"{prefix}_missing_samples", f"{prefix}_carrier_samples"
    if car_col in gdf.columns and m_col in gdf.columns:
        missing = {
            x.strip() for s in gdf[m_col].dropna() if s for x in str(s).split(";") if x.strip()
        }
        dosages: dict[str, int] = {}
        for s_str in gdf[car_col].dropna():
            if not s_str or not str(s_str).strip():
                continue
            for item in str(s_str).split(";"):
                if not item.strip():
                    continue
                sname, d, _ = _parse_carrier_token(item)
                if sname not in missing and d > 0:
                    dosages[sname] = max(dosages.get(sname, 0), d)
        return min(called, len(dosages)), min(chrom, sum(dosages.values()))

    v_col, a_col = f"{prefix}_variant_count", f"{prefix}_allele_count"
    v_cnt = min(called, int(gdf[v_col].sum())) if v_col in gdf.columns else 0
    a_cnt = min(chrom, int(gdf[a_col].sum())) if a_col in gdf.columns else 0
    return v_cnt, a_cnt


def _aggregate_gene_burden_legacy(df: pd.DataFrame) -> list[dict]:
    """Legacy aggregation: sum per-variant counts by gene."""
    gene_burden_data = []
    for gene, gene_df in df.groupby("GENE", observed=True):
        tot_p = _col_total(gene_df, "proband_count", "proband_called_count")
        tot_c = _col_total(gene_df, "control_count", "control_called_count")
        p_called, p_chrom = _called_and_chrom(gene_df, tot_p, "proband")
        c_called, c_chrom = _called_and_chrom(gene_df, tot_c, "control")
        p_var_count, p_allele_count = _legacy_counts(gene_df, p_called, p_chrom, "proband")
        c_var_count, c_allele_count = _legacy_counts(gene_df, c_called, c_chrom, "control")

        gene_burden_data.append(
            {
                "GENE": gene,
                "proband_count": tot_p,
                "control_count": tot_c,
                "proband_called_count": p_called,
                "control_called_count": c_called,
                "proband_called_chromosomes": p_chrom,
                "control_called_chromosomes": c_chrom,
                "proband_variant_count": p_var_count,
                "control_variant_count": c_var_count,
                "proband_allele_count": p_allele_count,
                "control_allele_count": c_allele_count,
            }
        )

    return gene_burden_data
