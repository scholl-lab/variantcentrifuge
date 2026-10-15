"""
Per-sample genotype column reconstruction and output collapsing utilities.

Reconstructs packed GT columns from per-sample GEN[i].GT columns and removes
per-sample GEN columns before materializing final TSV/Excel outputs.
"""

import logging
import re

import pandas as pd

logger = logging.getLogger(__name__)


def _find_per_sample_gt_columns(df: pd.DataFrame) -> list[str]:
    """Find per-sample GT columns in DataFrame.

    Matches columns like GEN_0__GT (sanitized) or GEN[0].GT (original).
    Returns them sorted by sample index.
    """
    gt_cols = []
    for col in df.columns:
        # Match sanitized: GEN_0__GT, GEN_123__GT
        # Match original: GEN[0].GT, GEN[123].GT
        if re.match(r"^GEN[_\[](\d+)[_\]\.]+(GT)$", col):
            gt_cols.append(col)
    # Sort by numeric index
    gt_cols.sort(key=lambda c: int(re.search(r"(\d+)", c).group(1)))  # type: ignore[union-attr]
    return gt_cols


def reconstruct_gt_column(
    df: pd.DataFrame,
    vcf_samples: list[str],
    extra_sample_fields: list[str] | None = None,
    delimiter: str = ":",
    sample_mapping: dict[str, str] | None = None,
) -> pd.DataFrame:
    """
    Reconstruct packed GT column from per-sample GT columns (Phase 11).

    Takes per-sample GT columns (GEN[0].GT / GEN_0__GT) from bcftools query output
    and reconstructs the packed format "Sample1(0/1);Sample2(1/1)" for TSV/Excel output.

    Only includes samples with variants (skips 0/0, ./., NA, empty) unless
    extra_sample_fields are explicitly requested, in which case 0/0 samples are
    retained alongside their extra sample fields.
    Drops the per-sample columns after reconstruction.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with per-sample GT columns
    vcf_samples : list of str
        List of sample IDs in VCF order
    extra_sample_fields : list of str, optional
        Explicit per-sample FORMAT fields (e.g. ["GEN[*].DP", "GEN[*].AD"]) to pack
        after GT using ``delimiter``.
    delimiter : str, default ":"
        Delimiter used to join GT and extra_sample_fields inside sample parentheses.
    sample_mapping : dict of str to str, optional
        Pseudonymization mapping from original sample IDs to pseudonymized IDs.

    Returns
    -------
    pd.DataFrame
        DataFrame with reconstructed GT column, per-sample columns dropped
    """
    gt_cols = _find_per_sample_gt_columns(df)
    if not gt_cols:
        logger.warning("No per-sample GT columns found for reconstruction")
        return df

    n_samples = min(len(gt_cols), len(vcf_samples))
    if len(gt_cols) != len(vcf_samples):
        logger.warning(
            f"GT column count ({len(gt_cols)}) != sample count ({len(vcf_samples)}), "
            f"using first {n_samples}"
        )

    s_cols: dict[tuple[int, str], str] = {}
    for col in df.columns:
        m = re.match(r"^GEN[_\[](?P<idx>\d+)[_\]\.]+(?P<field>[A-Za-z0-9_]+)$", col)
        if m and m.group("field") != "GT":
            s_cols[(int(m.group("idx")), m.group("field"))] = col

    if extra_sample_fields:
        norm = [re.sub(r"^GEN[\[_][^\].]+[\]_.]+", "", str(f)).strip() for f in extra_sample_fields]
        extra_fields = [
            f for f in norm if f and f != "GT" and any((i, f) in s_cols for i in range(n_samples))
        ]
    else:
        extra_fields = []

    logger.info(f"Reconstructing GT column from {n_samples} per-sample columns")

    # Vectorized approach: build GT strings using numpy for speed
    base_skip = {"./.", ".|.", ".", "NA", "", "nan", "None"}
    skip_values = base_skip if extra_fields else (base_skip | {"0/0", "0|0"})

    # Extract GT values as a 2D array for fast access
    gt_data = df[gt_cols[:n_samples]].values  # shape: (n_rows, n_samples)
    sample_ids = [
        sample_mapping.get(s, s) if sample_mapping else s for s in vcf_samples[:n_samples]
    ]

    gt_strings = []
    for row_idx in range(len(df)):
        entries = []
        for col_idx in range(n_samples):
            val = gt_data[row_idx, col_idx]
            if val is None:
                continue
            try:
                if pd.isna(val):
                    continue
            except (TypeError, ValueError):
                pass
            gt_str = str(val).strip()
            if gt_str in skip_values:
                continue
            if extra_fields:
                vals = [
                    ""
                    if (ecol := s_cols.get((col_idx, ef))) is None
                    or (rv := df[ecol].iat[row_idx]) is None
                    or (not isinstance(rv, str) and pd.isna(rv))
                    else str(rv).strip()
                    for ef in extra_fields
                ]
                entries.append(f"{sample_ids[col_idx]}({delimiter.join([gt_str, *vals])})")
            else:
                entries.append(f"{sample_ids[col_idx]}({gt_str})")
        gt_strings.append(";".join(entries))

    df["GT"] = gt_strings

    # Drop per-sample GT and requested extra GEN columns
    extra_cols = [col for (_, f), col in s_cols.items() if f in extra_fields]
    df = df.drop(columns=gt_cols + extra_cols)
    logger.info(f"Dropped {len(gt_cols)} per-sample GT columns from output")

    return df


def _collapse_per_sample_gt_columns_for_output(
    df: pd.DataFrame,
    vcf_samples: list[str] | None,
    *,
    no_replacement: bool = False,
    extra_sample_fields: list[str] | None = None,
    delimiter: str = ":",
    sample_mapping: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Remove raw per-sample GT columns from final outputs unless explicitly requested.

    Analysis stages may keep GEN[N].GT columns for burden/association logic. Final
    TSV/Excel output should expose only the packed GT representation by default.
    """
    if no_replacement:
        return df

    gt_cols = _find_per_sample_gt_columns(df)
    if not gt_cols:
        return df

    if ("GT" not in df.columns or extra_sample_fields) and vcf_samples:
        eff_mapping = dict(sample_mapping) if sample_mapping else {}
        if not eff_mapping and "GT" in df.columns:
            n_s = min(len(gt_cols), len(vcf_samples))
            std_skip = {"./.", ".|.", ".", "NA", "", "nan", "None", "0/0", "0|0"}
            for r_idx in range(len(df)):
                raw_gt = df["GT"].iat[r_idx]
                if not isinstance(raw_gt, str) or not raw_gt:
                    continue
                names = [m.group(1) for m in re.finditer(r"([^;()\s]+)\([^)]*\)", raw_gt)]
                active = [
                    i
                    for i in range(n_s)
                    if (v := df[gt_cols[i]].iat[r_idx]) is not None
                    and str(v).strip() not in std_skip
                ]
                if len(names) == len(active):
                    for i, name in zip(active, names, strict=False):
                        eff_mapping.setdefault(vcf_samples[i], name)
        df = reconstruct_gt_column(
            df,
            vcf_samples,
            extra_sample_fields,
            delimiter,
            sample_mapping=eff_mapping or None,
        )
    elif "GT" not in df.columns:
        logger.warning(
            "Per-sample GT columns found in final output but VCF sample names are "
            "unavailable; leaving columns unchanged"
        )
        return df

    gen_cols = [c for c in df.columns if re.match(r"^GEN[_\[](\d+)[_\]\.]+([A-Za-z0-9_]+)$", c)]
    if gen_cols:
        df = df.drop(columns=gen_cols)
        logger.info(f"Dropped {len(gen_cols)} per-sample GT columns from final output")
    return df
