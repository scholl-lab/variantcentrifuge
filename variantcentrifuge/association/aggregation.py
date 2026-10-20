"""Per-gene case/control burden aggregation shared by association and gene burden."""

import logging
from typing import Any

import pandas as pd

from ..genotype_utils import is_missing

logger = logging.getLogger("variantcentrifuge")


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
    missing_as_ref: bool = False,
) -> list[dict]:
    """Aggregate gene burden using per-sample GT columns (GEN_0__GT, GEN_1__GT, ...).

    Complete-case by default: a sample with any missing call in the gene is not counted.
    ``missing_as_ref`` counts missing calls as reference instead (#200).
    ``gt_columns`` must be ``GEN[0..n-1].GT`` for the n ``vcf_samples`` (#219).
    """
    from ..gt_reconstruction import gt_columns_match_samples

    if not gt_columns_match_samples(gt_columns, vcf_samples):
        raise ValueError(
            f"{len(gt_columns)} per-sample GT columns do not map 1:1 (GEN[0..n-1].GT, in "
            f"order) to the {len(vcf_samples)} VCF samples"
        )
    n_cols = len(gt_columns)
    case_col_names = [gt_columns[i] for i in range(n_cols) if vcf_samples[i] in case_samples]
    ctrl_col_names = [gt_columns[i] for i in range(n_cols) if vcf_samples[i] in control_samples]
    gene_burden_data = []

    def _tally(col_names: list[str], gdf: pd.DataFrame) -> tuple[int, int, int, int]:
        carrier, allele, called, called_chrom = 0, 0, 0, 0
        for col in col_names:
            vals = gdf[col]
            miss = [
                v is None or (not isinstance(v, str) and pd.isna(v)) or is_missing(str(v).strip())
                for v in vals
            ]
            if any(miss) and not missing_as_ref:
                continue
            called += 1
            # Ploidy from the explicit calls; "." (haploid missing) counts as haploid too.
            shown = [str(v) for v, m in zip(vals, miss, strict=True) if not m or str(v) == "."]
            s_ploidy = 1 if shown and all("/" not in v and "|" not in v for v in shown) else 2
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
    missing_as_ref: bool = False,
) -> list[dict]:
    """Aggregate gene burden counts from packed GT column using per-sample collapsing.

    Complete-case by default; ``missing_as_ref`` counts missing calls as reference (#200).
    """
    from ..helpers import extract_sample_and_genotype, genotype_to_allele_count

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
                    if is_missing(genotype) and missing_as_ref:
                        genotype = "0" if genotype.strip() == "." else "0/0"
                    elif is_missing(genotype):
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


def _col_total(gdf: pd.DataFrame, col: str, fallback_col: str) -> int:
    if col in gdf.columns:
        return int(gdf[col].iloc[0])
    return int(gdf[fallback_col].max()) if fallback_col in gdf.columns else 0


def _called_and_chrom(gdf: pd.DataFrame, tot: int, prefix: str, gene: Any) -> tuple[int, int]:
    """Compute callable sample count and callable chromosomes for legacy aggregation."""
    c_col, chrom_col = f"{prefix}_called_count", f"{prefix}_called_chromosomes"
    if c_col not in gdf.columns:
        called = tot
    elif len(gdf) > 1 and int((tot - gdf[c_col] > 0).sum()) > 0:
        logger.warning("Gene %s skipped: legacy counts cannot tell which samples are missing", gene)
        return 0, 0
    else:
        called = int(gdf[c_col].min())

    if called == 0:
        return 0, 0
    if chrom_col in gdf.columns:
        return called, max(0, min(called * 2, int(gdf[chrom_col].min())))
    return called, called * 2


def _legacy_counts(gdf: pd.DataFrame, called: int, chrom: int, prefix: str) -> tuple[int, int]:
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
        p_called, p_chrom = _called_and_chrom(gene_df, tot_p, "proband", gene)
        c_called, c_chrom = _called_and_chrom(gene_df, tot_c, "control", gene)
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
