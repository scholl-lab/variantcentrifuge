# File: variantcentrifuge/association/pca.py
# Location: variantcentrifuge/variantcentrifuge/association/pca.py
"""
PCA file loading, format detection, sample alignment, and covariate merge.

Supports three PCA file formats:
- PLINK .eigenvec with header (``#FID``/``FID`` first line)
- PLINK .eigenvec without header (FID + IID + numeric columns; two sample ID columns)
- AKT stdout / generic TSV (one sample ID column, remaining columns numeric)

The returned matrix is always float64 with shape ``(n_samples, n_components)`` aligned
to ``vcf_samples`` order.

Public API
----------
load_pca_file(filepath, vcf_samples, n_components=10)
    Load a PCA file and return aligned (matrix, column_names) tuple.
merge_pca_covariates(pca_matrix, pca_col_names, covariate_matrix, covariate_col_names)
    Horizontally stack PCA columns onto an existing covariate matrix.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger("variantcentrifuge")


def _is_numeric(s: str) -> bool:
    """Return True if *s* can be parsed as a float."""
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


def _detect_pca_format(lines: list[str], sep: str | None = None) -> str:
    """Detect PCA file format from the first few non-empty lines.

    Parameters
    ----------
    lines : list[str]
        First 5 (or fewer) non-empty lines from the file.
    sep : str or None
        Field separator; ``None`` splits on any whitespace.

    Returns
    -------
    str
        One of: ``"plink_header"``, ``"plink_nohdr"``, ``"akt_or_generic"``,
        ``"generic"``.
    """
    if not lines:
        return "generic"

    parts = [p.strip() for p in lines[0].strip().split(sep)]

    # PLINK .eigenvec with header: first token is exactly #FID or FID (a sample "fido" is data)
    if parts[0].lower() in ("#fid", "fid"):
        return "plink_header"

    # Inspect first non-empty data line (could be line 0 itself if no header)
    if len(parts) < 2:
        return "generic"

    # PLINK without header: two non-numeric ID columns then numeric columns
    if (
        not _is_numeric(parts[0])
        and not _is_numeric(parts[1])
        and len(parts) > 2
        and all(_is_numeric(p) for p in parts[2:])
    ):
        return "plink_nohdr"

    # AKT / generic: first column non-numeric, rest numeric
    if not _is_numeric(parts[0]) and all(_is_numeric(p) for p in parts[1:]):
        return "akt_or_generic"

    return "generic"


def _clean_ids(ids: pd.Index | pd.Series, remove_sample_substring: str | None) -> pd.Index:
    """Strip IDs and apply the VCF sample-name substring removal."""
    ids = pd.Index(ids).astype(str).str.strip()
    if remove_sample_substring:
        ids = ids.str.replace(remove_sample_substring, "", regex=False)
    return ids


def load_pca_file(
    filepath: str,
    vcf_samples: Sequence[str],
    n_components: int = 10,
    remove_sample_substring: str | None = None,
) -> tuple[np.ndarray, list[str]]:
    """Load a PCA file and return a sample-aligned PC matrix.

    Parameters
    ----------
    filepath : str
        Path to a PCA file.  Auto-detects one of three formats:

        1. **PLINK .eigenvec with header** — first line starts with ``#FID`` or
           ``FID``; uses ``IID`` column (second column) as sample ID.
        2. **PLINK .eigenvec without header** — first column is FID, second is
           IID (both non-numeric), remaining columns are PCs.  Uses IID.
        3. **AKT stdout / generic TSV** — first column is sample ID, remaining
           columns are numeric PCs.  Columns may be unnamed (AKT) or named
           (generic TSV header).

    vcf_samples : sequence of str
        Ordered list of sample IDs as they appear in the VCF.  The returned
        matrix rows will match this order exactly.
    n_components : int
        Number of principal components to extract.  If the file contains fewer
        PCs, a warning is logged and all available PCs are returned.  If
        ``n_components > 20`` a warning is also logged.

    Returns
    -------
    pc_matrix : np.ndarray, shape (n_samples, k), float64
        Principal component matrix aligned to ``vcf_samples`` order.
    col_names : list[str]
        Column names for each PC (``["PC1", "PC2", ...]``).

    Raises
    ------
    ValueError
        If any sample in ``vcf_samples`` is absent from the PCA file.
        All missing samples are reported in the error message.
        Also if ``n_components < 1``, the file has no PC columns, or a used PC
        cell is NaN/non-numeric (usually a delimiter or format problem).
    """
    if n_components > 20:
        logger.warning(
            "pca: requested %d principal components (>20). "
            "High PC counts risk overfitting and may not improve stratification correction.",
            n_components,
        )

    if n_components < 1:
        raise ValueError(f"pca: n_components must be >= 1, got {n_components}")

    # 1. Delimiter: .csv is comma; otherwise tab only for a clean TSV (same tab count on
    #    every line, no whitespace inside data fields), else any whitespace. Never sniff:
    #    csv.Sniffer misreads akt's "ID<TAB>  padded values" layout (#221).
    with open(filepath) as fh:
        lines = [ln.rstrip("\r\n") for ln in fh if ln.strip()]
    tab_rows = [ln.split("\t") for ln in lines]
    clean_tsv = (
        bool(tab_rows)
        and len(tab_rows[0]) >= 2
        and len({len(r) for r in tab_rows}) == 1
        # header names may contain spaces; data fields after the ID must not (akt pads with them)
        and not any(len(f.split()) > 1 for r in (tab_rows[1:] or tab_rows) for f in r[1:])
    )
    is_csv = os.path.splitext(filepath)[1].lower() == ".csv"
    sep = "," if is_csv else ("\t" if clean_tsv else r"\s+")

    # 2. Detect format and header from the first non-empty line
    split_sep = None if sep == r"\s+" else sep
    parts = [p.strip() for p in lines[0].split(split_sep)] if lines else []
    fmt = _detect_pca_format(lines[:1], split_sep)
    n_id = 2 if fmt in ("plink_header", "plink_nohdr") else 1
    # Header iff the first non-ID field is not a number
    has_header = fmt == "plink_header" or (len(parts) > n_id and not _is_numeric(parts[n_id]))

    # 3. Read all cells as strings; IID (PLINK) or the single ID column becomes the index
    df = pd.read_csv(filepath, sep=sep, header=0 if has_header else None, dtype=str)
    if not has_header:
        if n_id == 1 and df.shape[1] > 2:
            # Headerless "ID PC1 ..." vs PLINK "FID IID PC1 ..." with a numeric IID:
            # pick the layout whose ID column matches the VCF samples
            vcf_set = set(vcf_samples)
            n_col0, n_col1 = (
                len(vcf_set.intersection(_clean_ids(df.iloc[:, c], remove_sample_substring)))
                for c in (0, 1)
            )
            # exactly one layout may match; a partial match in the other is ambiguous
            if min(n_col0, n_col1) > 0 or n_col0 == n_col1:
                raise ValueError(
                    f"PCA file {filepath}: cannot tell whether column 1 (sample ID) or "
                    f"column 2 (PLINK FID/IID) holds the sample IDs ({n_col0} vs {n_col1} VCF "
                    "samples match). Add a header line, e.g. '#FID IID PC1 PC2 ...' or "
                    "'sample_id PC1 PC2 ...'."
                )
            n_id = 2 if n_col1 > n_col0 else 1
        df.columns = ["FID", "IID"][2 - n_id :] + [
            f"PC{i}" for i in range(1, df.shape[1] - n_id + 1)
        ]
    df.columns = [str(c).lstrip("#") for c in df.columns]
    df = df.set_index(df.columns[n_id - 1]).iloc[:, n_id - 1 :]
    df.index = _clean_ids(df.index, remove_sample_substring)

    # 4. Validate PC count (requested vs available)
    available = df.shape[1]
    if available == 0:
        raise ValueError(
            f"PCA file {filepath}: no principal component columns after the sample ID."
        )
    if available < n_components:
        logger.warning(
            "pca: requested %d components but file only has %d — using %d.",
            n_components,
            available,
            available,
        )
        n_components = available

    # 5. Coerce the used PCs to numbers; any NaN/non-numeric cell is a hard error
    df = df.iloc[:, :n_components].apply(pd.to_numeric, errors="coerce")
    bad = df.isna()
    if bad.any().any():
        bad_rows = bad.any(axis=1)
        raise ValueError(
            f"PCA file {filepath}: {int(bad.values.sum())} NaN/non-numeric cell(s) in "
            f"column(s) {list(bad.columns[bad.any()])}, e.g. sample(s) "
            f"{list(df.index[bad_rows][:3])} ({int(bad_rows.sum())} affected). Likely a "
            "delimiter/format problem: expected whitespace-separated (PLINK .eigenvec, "
            "akt pca stdout, TSV) or comma-separated (.csv) rows of sample ID(s) followed "
            "by numeric PCs, with all rows the same length."
        )

    dup = df.index[df.index.duplicated()].unique()
    if len(dup):
        raise ValueError(
            f"PCA file {filepath}: {len(dup)} sample ID(s) appear on more than one row, "
            f"e.g. {list(dup[:3])}. Each sample must have exactly one row."
        )

    # 6. Sample alignment — VCF samples missing from PCA file are an error
    vcf_samples_list = list(vcf_samples)
    missing = set(vcf_samples_list) - set(df.index)
    if missing:
        sorted_missing = sorted(missing)
        preview = sorted_missing[:10]
        suffix = f" ... and {len(sorted_missing) - 10} more" if len(sorted_missing) > 10 else ""
        raise ValueError(
            f"{len(sorted_missing)} VCF sample(s) missing from PCA file: {preview}{suffix}"
        )

    # Warn about extra PCA samples not in VCF (silently dropped)
    extra = set(df.index) - set(vcf_samples_list)
    if extra:
        extra_preview = sorted(extra)[:5]
        logger.warning(
            "pca: %d extra sample(s) in PCA file not in VCF (ignored): %s%s",
            len(extra),
            extra_preview,
            " ..." if len(extra) > 5 else "",
        )

    # 7. Reindex to VCF sample order
    df_aligned = df.reindex(vcf_samples_list)

    pc_matrix = df_aligned.values.astype(np.float64)
    col_names = [f"PC{i + 1}" for i in range(n_components)]

    return pc_matrix, col_names


def merge_pca_covariates(
    pca_matrix: np.ndarray,
    pca_col_names: list[str],
    covariate_matrix: np.ndarray | None,
    covariate_col_names: list[str],
) -> tuple[np.ndarray, list[str]]:
    """Merge PCA columns into an existing covariate matrix.

    Parameters
    ----------
    pca_matrix : np.ndarray, shape (n_samples, k_pca)
        Principal component matrix from ``load_pca_file``.
    pca_col_names : list[str]
        Column names for *pca_matrix* (e.g. ``["PC1", "PC2"]``).
    covariate_matrix : np.ndarray or None
        Existing covariate matrix, shape ``(n_samples, k_cov)``, or ``None``
        when no covariates have been loaded.
    covariate_col_names : list[str]
        Column names for *covariate_matrix*.  Ignored when
        *covariate_matrix* is ``None``.

    Returns
    -------
    merged_matrix : np.ndarray, shape (n_samples, k_pca + k_cov)
        Horizontally stacked matrix (covariates first, PCs appended).
    merged_col_names : list[str]
        Concatenated column names.
    """
    if covariate_matrix is None:
        return pca_matrix, list(pca_col_names)

    assert covariate_matrix.shape[0] == pca_matrix.shape[0], (
        f"Shape mismatch: covariate_matrix has {covariate_matrix.shape[0]} samples "
        f"but pca_matrix has {pca_matrix.shape[0]} samples"
    )

    merged = np.hstack([covariate_matrix, pca_matrix])
    merged_names = list(covariate_col_names) + list(pca_col_names)
    return merged, merged_names
