"""Unit tests for DataFrame optimizer utilities."""

import keyword
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from variantcentrifuge.dataframe_optimizer import (
    detect_categorical_columns,
    get_column_rename_map,
    load_optimized_dataframe,
    rename_invalid_identifiers,
    should_use_memory_passthrough,
)


@pytest.mark.unit
def test_detect_categorical_columns_low_cardinality():
    """Test that low-cardinality columns are detected as categorical."""
    # Create temp TSV with low and high cardinality columns
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        # CHROM: 3 unique values in 100 rows = 3% cardinality -> categorical
        # POS: 100 unique values in 100 rows = 100% cardinality -> str
        # REF: 4 unique values = 4% cardinality -> categorical
        f.write("CHROM\tPOS\tREF\n")
        for i in range(100):
            chrom = f"chr{i % 3 + 1}"  # chr1, chr2, chr3
            pos = str(1000 + i)  # All unique
            ref = ["A", "C", "G", "T"][i % 4]
            f.write(f"{chrom}\t{pos}\t{ref}\n")
        temp_path = f.name

    try:
        dtype_map = detect_categorical_columns(temp_path, sep="\t", cardinality_threshold=0.5)

        # Check results
        assert "CHROM" in dtype_map
        assert dtype_map["CHROM"] == "category"
        assert "POS" in dtype_map
        assert dtype_map["POS"] == "str"
        assert "REF" in dtype_map
        assert dtype_map["REF"] == "category"
    finally:
        Path(temp_path).unlink()


@pytest.mark.unit
def test_detect_categorical_columns_empty_file():
    """Test categorical detection on empty file."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("CHROM\tPOS\n")  # Header only
        temp_path = f.name

    try:
        dtype_map = detect_categorical_columns(temp_path)
        assert dtype_map == {}
    finally:
        Path(temp_path).unlink()


@pytest.mark.unit
def test_detect_categorical_columns_custom_threshold():
    """Test categorical detection with custom threshold."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("EFFECT\n")
        for i in range(100):
            # 10 unique values = 10% cardinality
            f.write(f"effect_{i % 10}\n")
        temp_path = f.name

    try:
        # With threshold=0.5, 10% < 50% -> categorical
        dtype_map = detect_categorical_columns(temp_path, cardinality_threshold=0.5)
        assert dtype_map["EFFECT"] == "category"

        # With threshold=0.05, 10% > 5% -> str
        dtype_map = detect_categorical_columns(temp_path, cardinality_threshold=0.05)
        assert dtype_map["EFFECT"] == "str"
    finally:
        Path(temp_path).unlink()


@pytest.mark.unit
def test_rename_invalid_identifiers_bracket_notation():
    """Test renaming columns with bracket notation."""
    df = pd.DataFrame(
        {
            "GEN[0].GT": [1, 2, 3],
            "ANN[0].EFFECT": [4, 5, 6],
            "CHROM": [7, 8, 9],  # Valid, should not change
        }
    )

    renamed_df, rename_map = rename_invalid_identifiers(df)

    # Check rename map
    assert "GEN[0].GT" in rename_map
    assert rename_map["GEN[0].GT"] == "GEN_0__GT"
    assert "ANN[0].EFFECT" in rename_map
    assert rename_map["ANN[0].EFFECT"] == "ANN_0__EFFECT"
    assert "CHROM" not in rename_map  # Valid identifier, no rename

    # Check DataFrame columns
    assert "GEN_0__GT" in renamed_df.columns
    assert "ANN_0__EFFECT" in renamed_df.columns
    assert "CHROM" in renamed_df.columns
    assert "GEN[0].GT" not in renamed_df.columns

    # Check data preserved
    pd.testing.assert_series_equal(renamed_df["GEN_0__GT"], df["GEN[0].GT"], check_names=False)
    pd.testing.assert_series_equal(renamed_df["CHROM"], df["CHROM"])


@pytest.mark.unit
def test_rename_invalid_identifiers_digit_start():
    """Test renaming columns starting with digits."""
    df = pd.DataFrame(
        {
            "0_foo": [1, 2],
            "123": [3, 4],
        }
    )

    renamed_df, rename_map = rename_invalid_identifiers(df)

    # Columns starting with digit should be prefixed
    assert rename_map["0_foo"] == "col_0_foo"
    assert rename_map["123"] == "col_123"
    assert "col_0_foo" in renamed_df.columns
    assert "col_123" in renamed_df.columns


@pytest.mark.unit
def test_rename_invalid_identifiers_keywords():
    """Test renaming Python keyword columns."""
    df = pd.DataFrame(
        {
            "class": [1, 2],
            "return": [3, 4],
            "if": [5, 6],
        }
    )

    _renamed_df, rename_map = rename_invalid_identifiers(df)

    # Keywords should be renamed
    assert "class" in rename_map
    assert "return" in rename_map
    assert "if" in rename_map
    assert all(keyword.iskeyword(col) or not col.isidentifier() for col in rename_map)


@pytest.mark.unit
def test_rename_invalid_identifiers_duplicates():
    """Test handling of duplicate column names after sanitization."""
    df = pd.DataFrame(
        {
            "foo.bar": [1],
            "foo-bar": [2],
            "foo bar": [3],
        }
    )

    _renamed_df, rename_map = rename_invalid_identifiers(df)

    # All would sanitize to 'foo_bar', so should get _1, _2 suffixes
    new_names = list(rename_map.values())
    assert len(new_names) == 3
    assert len(set(new_names)) == 3  # All unique
    assert "foo_bar" in new_names
    assert any("foo_bar_1" in name or "foo_bar_2" in name for name in new_names)


@pytest.mark.unit
def test_rename_invalid_identifiers_no_changes_needed():
    """Test that valid identifiers are not renamed."""
    df = pd.DataFrame(
        {
            "CHROM": [1, 2],
            "POS": [3, 4],
            "valid_name": [5, 6],
        }
    )

    renamed_df, rename_map = rename_invalid_identifiers(df)

    assert rename_map == {}
    assert renamed_df.columns.tolist() == df.columns.tolist()


@pytest.mark.unit
@patch("variantcentrifuge.dataframe_optimizer.psutil.virtual_memory")
def test_should_use_memory_passthrough_under_threshold(mock_memory):
    """Test pass-through decision when DataFrame is under threshold."""
    # Mock 4GB available memory
    mock_memory.return_value = MagicMock(available=4 * 1024**3)

    # Create small DataFrame (~1MB)
    df = pd.DataFrame({"col": ["x" * 100] * 1000})

    # With threshold=0.25, 1MB < 1GB -> should use pass-through
    result = should_use_memory_passthrough(df, threshold_ratio=0.25)
    assert result is True


@pytest.mark.unit
@patch("variantcentrifuge.dataframe_optimizer.psutil.virtual_memory")
def test_should_use_memory_passthrough_over_threshold(mock_memory):
    """Test pass-through decision when DataFrame exceeds threshold."""
    # Mock 1GB available memory
    mock_memory.return_value = MagicMock(available=1 * 1024**3)

    # Create large DataFrame (~500MB)
    df = pd.DataFrame({"col": ["x" * 1000] * 500000})

    # With threshold=0.25, 500MB > 250MB -> should NOT use pass-through
    result = should_use_memory_passthrough(df, threshold_ratio=0.25)
    assert result is False


@pytest.mark.unit
@patch("variantcentrifuge.dataframe_optimizer.psutil.virtual_memory")
def test_should_use_memory_passthrough_exception_handling(mock_memory):
    """Test that memory check exceptions return False (conservative fallback)."""
    # Mock exception
    mock_memory.side_effect = Exception("Memory check failed")

    df = pd.DataFrame({"col": [1, 2, 3]})

    result = should_use_memory_passthrough(df)
    assert result is False


@pytest.mark.unit
def test_load_optimized_dataframe_basic():
    """Test basic optimized DataFrame loading."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("CHROM\tPOS\tGEN[0].GT\n")
        for i in range(10):
            f.write(f"chr1\t{1000 + i}\t0/1\n")
        temp_path = f.name

    try:
        df, rename_map = load_optimized_dataframe(temp_path, sep="\t")

        # Check DataFrame loaded
        assert len(df) == 10
        assert len(df.columns) == 3

        # Check column sanitization happened
        assert "GEN_0__GT" in df.columns
        assert "GEN[0].GT" not in df.columns
        assert rename_map["GEN[0].GT"] == "GEN_0__GT"

        # Check categorical dtype applied
        assert df["CHROM"].dtype.name == "category"  # Low cardinality
    finally:
        Path(temp_path).unlink()


@pytest.mark.unit
def test_load_optimized_dataframe_pyarrow_engine():
    """Test that PyArrow engine is used when possible."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("CHROM\tPOS\n")
        f.write("chr1\t1000\n")
        temp_path = f.name

    try:
        # Mock pd.read_csv to verify engine parameter
        with patch("variantcentrifuge.dataframe_optimizer.pd.read_csv") as mock_read_csv:
            mock_read_csv.return_value = pd.DataFrame({"CHROM": ["chr1"], "POS": ["1000"]})

            _df, _rename_map = load_optimized_dataframe(temp_path, use_pyarrow=True)

            # Verify PyArrow engine was used
            call_kwargs = mock_read_csv.call_args[1]
            assert call_kwargs.get("engine") == "pyarrow"
    finally:
        Path(temp_path).unlink()


@pytest.mark.unit
def test_load_optimized_dataframe_pyarrow_fallback():
    """Test PyArrow fallback to C engine with unsupported params."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("CHROM\tPOS\n")
        f.write("chr1\t1000\n")
        temp_path = f.name

    try:
        # Mock pd.read_csv to verify engine parameter
        with patch("variantcentrifuge.dataframe_optimizer.pd.read_csv") as mock_read_csv:
            mock_read_csv.return_value = pd.DataFrame({"CHROM": ["chr1"], "POS": ["1000"]})

            # Pass unsupported param
            _df, _rename_map = load_optimized_dataframe(
                temp_path, use_pyarrow=True, on_bad_lines="warn"
            )

            # Verify PyArrow was NOT used (unsupported param)
            call_kwargs = mock_read_csv.call_args[1]
            assert call_kwargs.get("engine") != "pyarrow"
            assert call_kwargs.get("low_memory") is False  # C engine fallback
    finally:
        Path(temp_path).unlink()


@pytest.mark.unit
def test_load_optimized_dataframe_no_optimization():
    """Test loading with optimizations disabled."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("GEN[0].GT\tCHROM\n")
        f.write("0/1\tchr1\n")
        temp_path = f.name

    try:
        df, rename_map = load_optimized_dataframe(
            temp_path,
            optimize_dtypes=False,
            sanitize_columns=False,
        )

        # No categorical dtype — should be string-like (object or StringDtype)
        assert pd.api.types.is_string_dtype(df["CHROM"])

        # No column renaming
        assert "GEN[0].GT" in df.columns
        assert rename_map == {}
    finally:
        Path(temp_path).unlink()


@pytest.mark.unit
def test_load_optimized_dataframe_compression():
    """Test loading compressed file."""
    import gzip

    with tempfile.NamedTemporaryFile(mode="wb", suffix=".tsv.gz", delete=False) as f:
        with gzip.open(f, "wt") as gz:
            gz.write("CHROM\tPOS\n")
            gz.write("chr1\t1000\n")
        temp_path = f.name

    try:
        df, _rename_map = load_optimized_dataframe(temp_path, compression="gzip")

        assert len(df) == 1
        assert df["CHROM"].iloc[0] == "chr1"
    finally:
        Path(temp_path).unlink()


@pytest.mark.unit
def test_get_column_rename_map_present():
    """Test retrieving column rename map from context."""
    # Mock context
    context = MagicMock()
    context.column_rename_map = {"GEN[0].GT": "GEN_0__GT"}

    rename_map = get_column_rename_map(context)
    assert rename_map == {"GEN[0].GT": "GEN_0__GT"}


@pytest.mark.unit
def test_get_column_rename_map_missing():
    """Test retrieving column rename map when not present."""
    # Mock context without rename map
    context = MagicMock(spec=[])  # No attributes

    rename_map = get_column_rename_map(context)
    assert rename_map == {}


@pytest.mark.unit
def test_rename_invalid_identifiers_empty_after_sanitization():
    """Test handling of columns that become empty after sanitization."""
    df = pd.DataFrame(
        {
            "!!!": [1, 2],
            "...": [3, 4],
        }
    )

    renamed_df, _rename_map = rename_invalid_identifiers(df)

    # Should create valid names
    assert all(col.isidentifier() and not keyword.iskeyword(col) for col in renamed_df.columns)
    assert len(renamed_df.columns) == 2


@pytest.mark.unit
def test_load_optimized_dataframe_preserves_data():
    """Test that optimization preserves DataFrame data exactly."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("CHROM\tPOS\tREF\tALT\tGEN[0].GT\n")
        f.write("chr1\t1000\tA\tC\t0/1\n")
        f.write("chr2\t2000\tG\tT\t1/1\n")
        f.write("chr1\t3000\tC\tA\t0/0\n")
        temp_path = f.name

    try:
        # Load with optimizations
        df_opt, _rename_map_opt = load_optimized_dataframe(temp_path)

        # Load without optimizations
        df_plain, _rename_map_plain = load_optimized_dataframe(
            temp_path,
            use_pyarrow=False,
            optimize_dtypes=False,
            sanitize_columns=False,
        )

        # Data should be identical (accounting for column renames)
        assert len(df_opt) == len(df_plain)
        assert df_opt["CHROM"].tolist() == df_plain["CHROM"].tolist()
        assert df_opt["POS"].tolist() == df_plain["POS"].tolist()
        assert df_opt["GEN_0__GT"].tolist() == df_plain["GEN[0].GT"].tolist()
    finally:
        Path(temp_path).unlink()
