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

import csv
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


def _detect_pca_format(lines: list[str]) -> str:
    """Detect PCA file format from the first few non-empty lines.

    Parameters
    ----------
    lines : list[str]
        First 5 (or fewer) non-empty lines from the file.

    Returns
    -------
    str
        One of: ``"plink_header"``, ``"plink_nohdr"``, ``"akt_or_generic"``,
        ``"generic"``.
    """
    if not lines:
        return "generic"

    first = lines[0].strip()

    # PLINK .eigenvec with header: first line starts with #FID or FID
    first_lower = first.lower()
    if first_lower.startswith("#fid") or first_lower.startswith("fid"):
        return "plink_header"

    # Inspect first non-empty data line (could be line 0 itself if no header)
    parts = first.split()
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
    """
    if n_components > 20:
        logger.warning(
            "pca: requested %d principal components (>20). "
            "High PC counts risk overfitting and may not improve stratification correction.",
            n_components,
        )

    # 1. Auto-detect delimiter from extension, fallback to csv.Sniffer
    ext = os.path.splitext(filepath)[1].lower()
    if ext in (".tsv", ".tab"):
        sep: str | None = "\t"
    elif ext == ".csv":
        sep = ","
    elif ext == ".eigenvec":
        sep = None  # whitespace — PLINK and AKT both use mixed tab+space
    else:
        try:
            with open(filepath) as fh:
                sample_text = fh.read(2048)
            sep = csv.Sniffer().sniff(sample_text).delimiter
        except Exception:
            sep = None  # let pandas decide (whitespace for PLINK .eigenvec)

    # 2. Read first few lines to detect format (before pandas parsing)
    with open(filepath) as fh:
        raw_lines = [fh.readline() for _ in range(5)]
    non_empty_lines = [ln.strip() for ln in raw_lines if ln.strip()]
    fmt = _detect_pca_format(non_empty_lines)

    # 3. Read file with pandas
    # PLINK .eigenvec files are whitespace-delimited
    read_sep = sep if sep is not None else r"\s+"

    if fmt == "plink_header":
        # Header row present; strip leading # from column names
        df = pd.read_csv(filepath, sep=read_sep, header=0, engine="python")
        df.columns = [c.lstrip("#") for c in df.columns]
        # IID is the second column after FID
        id_col = df.columns[1]  # IID
        df = df.set_index(id_col)
        # Drop FID column (first column) if still present
        fid_col = df.columns[0] if len(df.columns) > 0 else None
        if fid_col is not None and fid_col.upper() == "FID":
            df = df.drop(columns=[fid_col])

    elif fmt == "plink_nohdr":
        # No header; first column = FID, second = IID, rest = PCs
        df = pd.read_csv(filepath, sep=read_sep, header=None, engine="python")
        # Assign names: FID, IID, PC1, PC2, ...
        n_cols = df.shape[1]
        pc_col_names_raw = [f"PC{i}" for i in range(1, n_cols - 1)]
        df.columns = ["FID", "IID", *pc_col_names_raw]
        df = df.set_index("IID")
        df = df.drop(columns=["FID"])

    else:
        # AKT stdout or generic TSV: first column = sample ID
        # Check if first line looks like a header (non-numeric first cell)
        has_header = False
        if non_empty_lines:
            parts = non_empty_lines[0].split()
            if parts and not _is_numeric(parts[0]):
                # Could be a header; check if second token is also non-numeric
                # (which would suggest it's a data row with non-numeric sample ID)
                # Generic TSV with header: all remaining tokens would also be string labels
                if len(parts) > 1 and not _is_numeric(parts[1]):
                    has_header = True
                elif len(parts) > 1 and _is_numeric(parts[1]):
                    # First column is sample ID, rest are numeric values — no header
                    has_header = False

        if has_header:
            df = pd.read_csv(filepath, sep=read_sep, header=0, index_col=0, engine="python")
        else:
            # AKT format: no header; first column = sample ID
            df = pd.read_csv(filepath, sep=read_sep, header=None, engine="python")
            n_cols = df.shape[1]
            pc_col_names_raw = [f"PC{i}" for i in range(1, n_cols)]
            df.columns = ["sample_id", *pc_col_names_raw]
            df = df.set_index("sample_id")

    # Ensure index is string
    df.index = df.index.astype(str)

    # Apply sample substring removal to match VCF sample name processing
    if remove_sample_substring:
        df.index = df.index.str.replace(remove_sample_substring, "", regex=False)

    # 4. Select only numeric columns (drop any remaining string columns)
    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if not numeric_cols:
        # Try converting all columns to numeric
        df = df.apply(pd.to_numeric, errors="coerce")
        numeric_cols = list(df.columns)
    else:
        df = df[numeric_cols]

    # 5. Slice to first n_components columns
    available = len(df.columns)
    if available < n_components:
        logger.warning(
            "pca: requested %d components but file only has %d — using %d.",
            n_components,
            available,
            available,
        )
        n_components = available

    df = df.iloc[:, :n_components]

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
    assert not df_aligned.isnull().any().any(), (
        "NaN in PCA matrix after reindex — unreachable if step 6 passed"
    )

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
