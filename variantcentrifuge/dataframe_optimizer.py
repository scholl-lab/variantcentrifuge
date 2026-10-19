"""DataFrame optimization utilities for the stage-based pipeline.

Provides optimizations for pandas DataFrame loading and processing:
- PyArrow engine for 5-15x CSV read speedup
- Auto-detection of low-cardinality columns for categorical dtype (50-75% memory reduction)
- Column name sanitization for itertuples compatibility
- Memory threshold-based pass-through decision making

All optimizations preserve byte-identical output.
"""

import csv
import gc
import keyword
import logging
import re
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
import psutil

logger = logging.getLogger(__name__)


def detect_categorical_columns(
    csv_path: str | Path,
    sep: str = "\t",
    cardinality_threshold: float = 0.5,
    max_sample_rows: int = 10000,
) -> dict[str, str]:
    """Detect low-cardinality columns suitable for categorical dtype.

    Samples the first N rows of a CSV file and identifies columns with
    low cardinality (few unique values relative to row count). These
    columns benefit from categorical dtype conversion, which provides
    significant memory reduction (typically 50-75%).

    Parameters
    ----------
    csv_path : str or Path
        Path to the CSV/TSV file to analyze
    sep : str, default "\t"
        Field separator
    cardinality_threshold : float, default 0.5
        If (unique_values / row_count) < threshold, column is categorical
    max_sample_rows : int, default 10000
        Number of rows to sample for detection

    Returns
    -------
    dict[str, str]
        Dictionary mapping column names to dtype strings ("category" or "str")
        Suitable for passing to pd.read_csv(dtype=...)

    Examples
    --------
    >>> dtype_map = detect_categorical_columns("variants.tsv")
    >>> df = pd.read_csv("variants.tsv", sep="\t", dtype=dtype_map)
    """
    # Sample first N rows for analysis
    try:
        sample_df = pd.read_csv(
            csv_path,
            sep=sep,
            nrows=max_sample_rows,
            dtype=str,
            keep_default_na=False,
            na_values=[""],
        )
    except Exception as e:
        logger.warning(f"Failed to sample {csv_path} for categorical detection: {e}")
        return {}

    if len(sample_df) == 0:
        logger.debug("Empty DataFrame, no categorical columns detected")
        return {}

    dtype_map = {}
    categorical_cols = []

    for col in sample_df.columns:
        nunique = sample_df[col].nunique()
        row_count = len(sample_df)
        cardinality = nunique / row_count if row_count > 0 else 1.0

        if cardinality < cardinality_threshold:
            dtype_map[col] = "category"
            categorical_cols.append(col)
        else:
            dtype_map[col] = "str"

    if categorical_cols:
        logger.debug(
            f"Detected {len(categorical_cols)} categorical columns "
            f"(threshold={cardinality_threshold}): {categorical_cols[:10]}"
            + ("..." if len(categorical_cols) > 10 else "")
        )
    else:
        logger.debug("No categorical columns detected")

    return dtype_map


def rename_invalid_identifiers(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, str]]:
    """Rename columns with invalid Python identifiers for itertuples compatibility.

    Columns with names like 'GEN[0].GT', 'ANN[0].EFFECT', or starting with digits
    cannot be accessed as attributes in itertuples namedtuples. This function
    sanitizes them to valid Python identifiers.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with potentially invalid column names

    Returns
    -------
    tuple[pd.DataFrame, dict[str, str]]
        (renamed_df, rename_map) where rename_map maps old_name -> new_name
        Allows reversing the mapping for output if needed

    Examples
    --------
    >>> df = pd.DataFrame({"GEN[0].GT": [1, 2], "CHROM": [3, 4]})
    >>> renamed_df, rename_map = rename_invalid_identifiers(df)
    >>> rename_map
    {'GEN[0].GT': 'GEN_0__GT'}
    >>> renamed_df.columns.tolist()
    ['GEN_0__GT', 'CHROM']
    """
    rename_map = {}
    new_columns = []
    seen_names = set()

    for col in df.columns:
        new_name = col

        # Check if valid identifier and not a keyword
        if not col.isidentifier() or keyword.iskeyword(col):
            # Sanitize: replace non-alphanumeric chars with underscore
            new_name = re.sub(r"[^a-zA-Z0-9_]", "_", col)

            # Prefix with 'col_' if starts with digit
            if new_name and new_name[0].isdigit():
                new_name = f"col_{new_name}"

            # Handle empty strings after sanitization
            if not new_name:
                new_name = "col_unnamed"

            # Handle duplicates by appending _N suffix
            original_new_name = new_name
            counter = 1
            while new_name in seen_names:
                new_name = f"{original_new_name}_{counter}"
                counter += 1

            rename_map[col] = new_name
            logger.info(f"Renamed column '{col}' -> '{new_name}' for itertuples compatibility")

        seen_names.add(new_name)
        new_columns.append(new_name)

    # Rename the DataFrame
    df_renamed = df.copy()
    df_renamed.columns = new_columns

    if rename_map:
        logger.info(f"Renamed {len(rename_map)} columns with invalid identifiers")
    else:
        logger.debug("No columns needed renaming")

    return df_renamed, rename_map


def should_use_memory_passthrough(
    df: pd.DataFrame,
    threshold_ratio: float = 0.25,
) -> bool:
    """Determine if DataFrame is small enough for in-memory pass-through.

    Checks if DataFrame memory footprint is under threshold relative to
    available system RAM. Used to decide between in-memory vs disk-based
    stage pass-through.

    Targets normal desktops with 8-16GB RAM. Conservative threshold
    prevents OOM while still benefiting from in-memory optimization.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to check
    threshold_ratio : float, default 0.25
        Max ratio of DataFrame memory to available RAM (0.0-1.0)
        Default 0.25 = use at most 25% of available memory

    Returns
    -------
    bool
        True if DataFrame fits within memory threshold

    Examples
    --------
    >>> df = pd.DataFrame(...)
    >>> if should_use_memory_passthrough(df):
    ...     context.variants_df = df  # In-memory pass-through
    ... else:
    ...     # Fall back to disk-based pass-through
    """
    try:
        available_memory = psutil.virtual_memory().available
        df_memory = df.memory_usage(deep=True).sum()

        threshold_memory = available_memory * threshold_ratio
        use_passthrough = bool(df_memory <= threshold_memory)

        logger.info(
            f"DataFrame memory: {df_memory / 1024**2:.1f} MB, "
            f"Available: {available_memory / 1024**2:.1f} MB, "
            f"Threshold: {threshold_memory / 1024**2:.1f} MB "
            f"({threshold_ratio * 100:.0f}%) -> "
            f"{'Using' if use_passthrough else 'Not using'} in-memory pass-through"
        )

        return use_passthrough
    except Exception as e:
        logger.warning(f"Failed to check memory for pass-through decision: {e}")
        # Conservative fallback: don't use pass-through if can't determine memory
        return False


def load_optimized_dataframe(
    file_path: str | Path,
    sep: str = "\t",
    compression: str | None = None,
    use_pyarrow: bool = True,
    optimize_dtypes: bool = True,
    sanitize_columns: bool = True,
    quoting: int = csv.QUOTE_NONE,
    **extra_kwargs: Any,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Load CSV/TSV with all DataFrame optimizations applied.

    Main entry point for optimized DataFrame loading. Combines:
    - PyArrow engine for 5-15x read speedup
    - Categorical dtype auto-detection for 50-75% memory reduction
    - Column name sanitization for itertuples compatibility

    All optimizations preserve byte-identical pipeline output.

    Parameters
    ----------
    file_path : str or Path
        Path to CSV/TSV file
    sep : str, default "\t"
        Field separator
    compression : str or None, default None
        Compression type ('gzip', 'bz2', etc.)
    use_pyarrow : bool, default True
        Use PyArrow engine if possible (falls back to C engine if incompatible params)
    optimize_dtypes : bool, default True
        Auto-detect categorical columns
    sanitize_columns : bool, default True
        Rename invalid Python identifiers
    quoting : int, default csv.QUOTE_NONE
        Control CSV quoting behavior (e.g. csv.QUOTE_MINIMAL for TSV artifacts with JSON)
    **extra_kwargs
        Additional arguments passed to pd.read_csv

    Returns
    -------
    tuple[pd.DataFrame, dict[str, str]]
        (df, column_rename_map) where rename_map maps old_name -> new_name

    Examples
    --------
    >>> df, rename_map = load_optimized_dataframe("variants.tsv")
    >>> # Use df for analysis with itertuples
    >>> for row in df.itertuples(index=False):
    ...     print(row.CHROM, row.POS, row.GEN_0__GT)
    >>> # Store rename_map for output reversal if needed
    >>> context.column_rename_map = rename_map
    """
    file_path = Path(file_path)

    # Log memory before loading
    mem_before = psutil.Process().memory_info().rss / 1024**2

    # Step 1: Detect categorical columns if requested
    dtype_map = {}
    if optimize_dtypes:
        dtype_map = detect_categorical_columns(file_path, sep=sep)

    # Step 2: Determine engine and load DataFrame
    can_use_pyarrow = (
        use_pyarrow
        and quoting in (csv.QUOTE_NONE, csv.QUOTE_MINIMAL)
        and not (set(extra_kwargs) - {"on_bad_lines"})
        and extra_kwargs.get("on_bad_lines", "error") in ("error", "warn", "skip")
    )

    df: pd.DataFrame | None = None
    if can_use_pyarrow:
        try:
            import gzip
            import warnings

            import pyarrow as pa  # type: ignore[import-untyped]
            import pyarrow.csv as pv  # type: ignore[import-untyped]

            opener = (
                gzip.open if (compression == "gzip" or str(file_path).endswith(".gz")) else open
            )
            with opener(file_path, "rt", encoding="utf-8", errors="replace") as f:
                first_line = f.readline()

            header_text = first_line.lstrip("\ufeff").rstrip("\r\n")
            if header_text:
                if quoting != csv.QUOTE_NONE:
                    header_cols = next(
                        csv.reader([header_text], delimiter=sep, quoting=cast(Any, quoting))
                    )
                else:
                    header_cols = header_text.split(sep)
                if len(header_cols) != len(set(header_cols)):
                    raise ValueError("Duplicate column headers; falling back to C parser")
                on_bad_lines = extra_kwargs.get("on_bad_lines", "error")

                def invalid_handler(row: Any) -> str:
                    if row.actual_columns < row.expected_columns:
                        return "error"
                    if on_bad_lines == "warn":
                        warnings.warn(
                            f"Skipping line: expected {row.expected_columns} fields, "
                            f"saw {row.actual_columns}",
                            pd.errors.ParserWarning,
                            stacklevel=2,
                        )
                        return "skip"
                    elif on_bad_lines == "skip":
                        return "skip"
                    return "error"

                if compression:
                    source = pa.input_stream(str(file_path), compression=compression)
                elif str(file_path).endswith(".gz"):
                    source = pa.input_stream(str(file_path), compression="gzip")
                else:
                    source = pa.input_stream(str(file_path))

                quote_char: str | bool = '"' if quoting == csv.QUOTE_MINIMAL else False
                double_quote: bool = quoting == csv.QUOTE_MINIMAL

                table = pv.read_csv(
                    source,
                    parse_options=pv.ParseOptions(
                        delimiter=sep,
                        quote_char=quote_char,
                        double_quote=double_quote,
                        invalid_row_handler=invalid_handler
                        if on_bad_lines in ("warn", "skip")
                        else None,
                    ),
                    convert_options=pv.ConvertOptions(
                        column_types={c: pa.string() for c in header_cols},
                        strings_can_be_null=True,
                        null_values=[""],
                    ),
                )
                df = table.to_pandas()
                df = df.fillna(np.nan)
                if dtype_map:
                    for col, dt in dtype_map.items():
                        if col in df.columns and dt == "category":
                            df[col] = df[col].astype("category")
                logger.debug("Successfully loaded DataFrame using PyArrow engine")
        except Exception as e:
            logger.debug(f"PyArrow loading failed, falling back to C engine: {e}")
            df = None

    if df is None:
        read_kwargs = {
            "sep": sep,
            "dtype": dtype_map if dtype_map else str,
            "keep_default_na": False,
            "na_values": [""],
            "compression": compression,
            "quoting": quoting,
            "low_memory": False,
            **extra_kwargs,
        }
        df = pd.read_csv(file_path, **read_kwargs)

    # Step 4: Sanitize columns if requested
    column_rename_map: dict[str, str] = {}
    if sanitize_columns:
        df, column_rename_map = rename_invalid_identifiers(df)

    # Log memory after loading
    gc.collect()  # Force GC to get accurate post-load memory
    mem_after = psutil.Process().memory_info().rss / 1024**2
    mem_delta = mem_after - mem_before

    logger.info(
        f"Loaded {len(df)} rows x {len(df.columns)} columns "
        f"(Memory: {df.memory_usage(deep=True).sum() / 1024**2:.1f} MB, "
        f"Process RSS delta: {mem_delta:+.1f} MB)"
    )

    return df, column_rename_map


def get_column_rename_map(context: Any) -> dict[str, str]:
    """Retrieve stored column rename mapping from PipelineContext.

    Helper to get the rename map stored by DataFrameLoadingStage.
    Useful for output stages that need to reverse column renames.

    Parameters
    ----------
    context : PipelineContext
        Pipeline context containing column_rename_map

    Returns
    -------
    dict[str, str]
        Column rename mapping (old_name -> new_name), or empty dict if not found

    Examples
    --------
    >>> rename_map = get_column_rename_map(context)
    >>> if rename_map:
    ...     # Reverse mapping for output
    ...     reverse_map = {v: k for k, v in rename_map.items()}
    ...     df.rename(columns=reverse_map, inplace=True)
    """
    if hasattr(context, "column_rename_map"):
        return dict(context.column_rename_map)
    else:
        logger.debug("No column_rename_map found in context")
        return {}
