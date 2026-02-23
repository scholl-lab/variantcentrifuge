# File: variantcentrifuge/association/covariates.py
# Location: variantcentrifuge/variantcentrifuge/association/covariates.py
"""
Covariate loading, sample alignment, and one-hot encoding for association tests.

Provides ``load_covariates()`` which reads a delimited covariate file, aligns
rows to VCF sample order, one-hot encodes categorical columns, and checks for
high multicollinearity.

The returned array is always float64 with shape ``(n_samples, k)`` and
guaranteed to contain no NaN (ValueError raised before any stats run if
alignment would produce NaN).
"""

from __future__ import annotations

import csv
import logging
import os
from collections.abc import Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger("variantcentrifuge")


def load_covariates(
    filepath: str,
    vcf_samples: Sequence[str],
    covariate_columns: list[str] | None = None,
    categorical_columns: list[str] | None = None,
) -> tuple[np.ndarray, list[str]]:
    """
    Load, align, and encode a covariate file for association analysis.

    Parameters
    ----------
    filepath : str
        Path to a delimited file where the first column contains sample IDs
        and the first row is a header. Delimiter is auto-detected from file
        extension (.tsv/.tab -> tab, .csv -> comma) with csv.Sniffer fallback.
    vcf_samples : sequence of str
        Ordered list of sample IDs as they appear in the VCF. The returned
        matrix rows will match this order exactly.
    covariate_columns : list[str] | None
        If provided, only these columns are retained after loading. None
        means all columns are kept.
    categorical_columns : list[str] | None
        Columns to one-hot encode with ``pd.get_dummies(drop_first=True)``.
        None means auto-detect: non-numeric columns with ``nunique() <= 5``.

    Returns
    -------
    X : np.ndarray, shape (n_samples, k), float64
        Covariate matrix aligned to ``vcf_samples`` order. No NaN present.
    column_names : list[str]
        Names of columns in ``X`` after encoding (useful for diagnostics).

    Raises
    ------
    ValueError
        If any sample in ``vcf_samples`` is absent from the covariate file.
        All missing samples are reported in the error message.

    Notes
    -----
    Extra samples in the covariate file that are not in ``vcf_samples`` are
    silently dropped after a logged warning (up to 5 IDs shown).

    High multicollinearity (condition number > 1000) triggers a WARNING but
    does not abort — the caller decides whether to proceed.
    """
    # 1. Auto-detect delimiter from extension, fallback to csv.Sniffer
    ext = os.path.splitext(filepath)[1].lower()
    if ext in (".tsv", ".tab"):
        sep = "\t"
    elif ext == ".csv":
        sep = ","
    else:
        try:
            with open(filepath) as fh:
                sample_text = fh.read(2048)
            sep = csv.Sniffer().sniff(sample_text).delimiter
        except Exception:
            sep = "\t"  # bioinformatics default

    # 2. Load file (first column = sample ID, header required)
    df = pd.read_csv(filepath, sep=sep, index_col=0)
    # Ensure index is string (sample IDs from VCF are always strings)
    df.index = df.index.astype(str)

    # 3. Column selection
    if covariate_columns is not None:
        df = df[covariate_columns]

    # 4. VCF samples missing from covariate file -> abort before any stats
    vcf_samples_list = list(vcf_samples)
    missing = set(vcf_samples_list) - set(df.index)
    if missing:
        sorted_missing = sorted(missing)
        preview = sorted_missing[:10]
        suffix = f" ... and {len(sorted_missing) - 10} more" if len(sorted_missing) > 10 else ""
        raise ValueError(
            f"{len(sorted_missing)} VCF sample(s) missing from covariate file: {preview}{suffix}"
        )

    # 5. Warn about extra covariate samples not in VCF (then ignore them)
    extra = set(df.index) - set(vcf_samples_list)
    if extra:
        extra_preview = sorted(extra)[:5]
        logger.warning(
            "Covariate file has %d extra sample(s) not in VCF: %s%s",
            len(extra),
            extra_preview,
            " ..." if len(extra) > 5 else "",
        )

    # 6. Align to VCF sample order (CRITICAL)
    df_aligned = df.reindex(vcf_samples_list)
    # This assertion is unreachable if step 4 passed, but guards against bugs
    assert not df_aligned.isnull().any().any(), (
        "NaN in covariate matrix after reindex — unreachable if step 4 passed"
    )

    # 7. One-hot encode categorical columns
    if categorical_columns is None:
        # Auto-detect: non-numeric columns with <=5 unique values
        categorical_columns = [
            col
            for col in df_aligned.columns
            if not pd.api.types.is_numeric_dtype(df_aligned[col]) and df_aligned[col].nunique() <= 5
        ]

    if categorical_columns:
        df_aligned = pd.get_dummies(
            df_aligned,
            columns=categorical_columns,
            drop_first=True,
            dtype=float,
        )

    # 8. Multicollinearity check (warn only, never abort)
    covariate_matrix = df_aligned.values.astype(np.float64)
    if covariate_matrix.shape[1] >= 2:
        cond_num = np.linalg.cond(covariate_matrix)
        if cond_num > 1000:
            logger.warning(
                "High multicollinearity in covariate matrix "
                "(condition number: %.1f). Results may be unreliable.",
                cond_num,
            )

    column_names = list(df_aligned.columns)
    return covariate_matrix, column_names
