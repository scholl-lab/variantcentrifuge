"""
Per-sample genotype column reconstruction and output collapsing utilities.

Reconstructs packed GT columns from per-sample GEN[i].GT columns and removes
per-sample GEN columns before materializing final TSV/Excel outputs.
"""

import logging
import re

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
_GT_COLUMN = re.compile(r"^GEN[_\[](\d+)[_\]\.]+(GT)$")


def _gt_column_index(col: str) -> int:
    """GEN index of a per-sample GT column name (``GEN_3__GT`` -> 3), else -1."""
    m = _GT_COLUMN.match(col)
    return int(m.group(1)) if m else -1


def gt_columns_match_samples(gt_columns: list[str], vcf_samples: list[str] | None) -> bool:
    """True when *gt_columns* are exactly ``GEN[0..n-1].GT``, in order, for n *vcf_samples*.

    Per-sample columns carry the VCF sample index, not the ID; only a 1:1 set can be
    paired with the sample list by position (#219).
    """
    return bool(vcf_samples) and [_gt_column_index(c) for c in gt_columns] == list(
        range(len(vcf_samples or []))
    )


def _copy_gt_columns_to_samples(
    df: pd.DataFrame, gt_cols: list[str], vcf_samples: list[str]
) -> int:
    """Copy each ``GEN[i].GT`` column to a column named ``vcf_samples[i]``, in place.

    Columns whose index has no VCF sample are skipped. Returns the number copied.
    """
    n = 0
    for col in gt_cols:
        i = _gt_column_index(col)
        if 0 <= i < len(vcf_samples):
            df[vcf_samples[i]] = df[col]
            n += 1
    return n


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
    df = df.copy(deep=False)
    gt_cols = _find_per_sample_gt_columns(df)
    if not gt_cols:
        logger.warning("No per-sample GT columns found for reconstruction")
        return df

    if not gt_columns_match_samples(gt_cols, vcf_samples):
        # A missing GEN index would read that sample as reference in the packed GT (#219).
        seen = {_gt_column_index(c) for c in gt_cols}
        missing = [i for i in range(len(vcf_samples)) if i not in seen]
        raise ValueError(
            f"Cannot reconstruct GT: {len(gt_cols)} per-sample GT columns for "
            f"{len(vcf_samples)} VCF samples; GEN[0..n-1].GT are required, missing indices "
            f"{missing[:10]}"
        )
    indexed = [(i, c) for i, c in enumerate(gt_cols)]
    n_samples = len(indexed)

    s_cols: dict[tuple[int, str], str] = {}
    for col in df.columns:
        m = re.match(r"^GEN[_\[](?P<idx>\d+)[_\]\.]+(?P<field>[A-Za-z0-9_]+)$", col)
        if m and m.group("field") != "GT":
            s_cols[(int(m.group("idx")), m.group("field"))] = col

    if extra_sample_fields:
        norm = [re.sub(r"^GEN[\[_][^\].]+[\]_.]+", "", str(f)).strip() for f in extra_sample_fields]
        extra_fields = [
            f for f in norm if f and f != "GT" and any((i, f) in s_cols for i, _ in indexed)
        ]
    else:
        extra_fields = []

    logger.info(f"Reconstructing GT column from {n_samples} per-sample columns")

    skip_values = set() if extra_fields else {"0/0", "0|0"}
    missing_sentinels = {"NA", "", "nan", "None"}

    n_rows = len(df)
    sample_cols = [c for _, c in indexed]
    sample_ids = [
        sample_mapping.get(vcf_samples[i], vcf_samples[i]) if sample_mapping else vcf_samples[i]
        for i, _ in indexed
    ]

    # Per extra field: (column names present, sample index -> position in that list)
    extra_layout: dict[str, tuple[list[str], dict[int, int]]] = {}
    for ef in extra_fields:
        present = [k for k, (i, _) in enumerate(indexed) if (i, ef) in s_cols]
        extra_layout[ef] = (
            [s_cols[(indexed[k][0], ef)] for k in present],
            {k: p for p, k in enumerate(present)},
        )

    def _clean(rv: object) -> str:
        return "" if rv is None or (not isinstance(rv, str) and pd.isna(rv)) else str(rv).strip()

    # ~1M cells per batch bounds peak memory (~100 bytes per cell in flight). Within a batch
    # every cell is normalised at once and Python only touches the cells that are kept.
    batch_size = max(1, 1_000_000 // max(1, n_samples))
    gt_strings: list[str] = []
    for start in range(0, n_rows, batch_size):
        batch = df.iloc[start : start + batch_size]
        flat = pd.Series(batch[sample_cols].to_numpy(dtype=object).ravel(), dtype=object)
        gt_flat = flat.astype(str).str.strip()
        gt_flat = gt_flat.where(~(flat.isna() | gt_flat.isin(missing_sentinels)), "./.")
        gt_vals = gt_flat.to_numpy(dtype=object)
        kept = np.flatnonzero(~gt_flat.isin(skip_values).to_numpy())
        extra_vals = {
            ef: (batch[cols].to_numpy(dtype=object), pos)
            for ef, (cols, pos) in extra_layout.items()
        }

        entries_by_row: list[list[str]] = [[] for _ in range(len(batch))]
        for flat_idx in kept:
            r_idx, col_idx = divmod(int(flat_idx), n_samples)
            gt_str = gt_vals[flat_idx]
            if extra_fields:
                vals = []
                for ef in extra_fields:  # keeps requested order, duplicates included
                    arr, pos = extra_vals[ef]
                    vals.append(_clean(arr[r_idx, pos[col_idx]]) if col_idx in pos else "")
                entry = f"{sample_ids[col_idx]}({delimiter.join([gt_str, *vals])})"
            else:
                entry = f"{sample_ids[col_idx]}({gt_str})"
            entries_by_row[r_idx].append(entry)
        gt_strings.extend(";".join(entries) for entries in entries_by_row)

    df["GT"] = gt_strings

    # Drop per-sample GT and requested extra GEN columns
    # Remove per-sample GT and requested extra GEN columns with ONE delete: deleting them
    # one by one re-splits the remaining blocks each time (quadratic, #173). Relabel them
    # to a shared sentinel so a single ``del`` keeps the other blocks shared with the input;
    # DataFrame.drop would copy every kept column without Copy-on-Write (pandas < 3).
    drop = set(gt_cols) | {col for (_, f), col in s_cols.items() if f in extra_fields}
    sentinel = object()
    df.columns = pd.Index([sentinel if c in drop else c for c in df.columns], name=df.columns.name)
    del df[sentinel]
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
            idx = [_gt_column_index(c) for c in gt_cols]
            std_skip = {"0/0", "0|0"}
            for r_idx in range(len(df)):
                raw_gt = df["GT"].iat[r_idx]
                if not isinstance(raw_gt, str) or not raw_gt:
                    continue
                names = [m.group(1) for m in re.finditer(r"([^;()\s]+)\([^)]*\)", raw_gt)]
                active = [
                    i
                    for c, i in zip(gt_cols, idx, strict=True)
                    if 0 <= i < len(vcf_samples)
                    and (v := df[c].iat[r_idx]) is not None
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
