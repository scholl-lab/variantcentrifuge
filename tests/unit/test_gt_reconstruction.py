"""
Tests for GT column reconstruction from per-sample columns (Phase 11).

These tests prove that reconstruct_gt_column() produces identical output format
to the old genotype replacement pipeline. Per-sample columns use the naming
convention from bcftools query: GEN[0].GT, GEN[1].GT, etc. (or sanitized
GEN_0__GT, GEN_1__GT after DataFrame loading).
"""

import pandas as pd
import pytest

from variantcentrifuge.stages.output_stages import (
    _find_per_sample_gt_columns,
    reconstruct_gt_column,
)


def _make_gt_col(idx: int, sanitized: bool = True) -> str:
    """Helper to create GT column name."""
    return f"GEN_{idx}__GT" if sanitized else f"GEN[{idx}].GT"


@pytest.mark.unit
def test_find_per_sample_gt_columns_sanitized():
    """Test detection of sanitized GT column names."""
    df = pd.DataFrame({"CHROM": [], "GEN_0__GT": [], "GEN_1__GT": [], "GEN_2__GT": []})
    cols = _find_per_sample_gt_columns(df)
    assert cols == ["GEN_0__GT", "GEN_1__GT", "GEN_2__GT"]


@pytest.mark.unit
def test_find_per_sample_gt_columns_original():
    """Test detection of original (unsanitized) GT column names."""
    df = pd.DataFrame({"CHROM": [], "GEN[0].GT": [], "GEN[1].GT": []})
    cols = _find_per_sample_gt_columns(df)
    assert cols == ["GEN[0].GT", "GEN[1].GT"]


@pytest.mark.unit
def test_find_per_sample_gt_columns_sorted():
    """Test that columns are sorted by index, not lexicographically."""
    df = pd.DataFrame({"GEN_10__GT": [], "GEN_2__GT": [], "GEN_0__GT": [], "GEN_1__GT": []})
    cols = _find_per_sample_gt_columns(df)
    assert cols == ["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_10__GT"]


@pytest.mark.unit
def test_find_per_sample_gt_columns_no_match():
    """Test that non-GT columns are not matched."""
    df = pd.DataFrame({"CHROM": [], "POS": [], "GENE": [], "GT": []})
    cols = _find_per_sample_gt_columns(df)
    assert cols == []


@pytest.mark.unit
def test_basic_reconstruction():
    """Test basic GT reconstruction with mixed genotypes."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [12345, 67890],
            "GEN_0__GT": ["0/1", "1/1"],
            "GEN_1__GT": ["0/0", "0/1"],
            "GEN_2__GT": ["1/1", "./."],
        }
    )
    vcf_samples = ["Sample1", "Sample2", "Sample3"]

    result = reconstruct_gt_column(df, vcf_samples)

    # Row 0: Sample1(0/1) and Sample3(1/1) - skip Sample2 (0/0)
    # Row 1: Sample1(1/1) and Sample2(0/1) - skip Sample3 (./.)
    assert "GT" in result.columns
    assert result["GT"].iloc[0] == "Sample1(0/1);Sample3(1/1)"
    assert result["GT"].iloc[1] == "Sample1(1/1);Sample2(0/1)"

    # Per-sample columns should be dropped
    assert "GEN_0__GT" not in result.columns
    assert "GEN_1__GT" not in result.columns
    assert "GEN_2__GT" not in result.columns


@pytest.mark.unit
def test_all_reference_genotypes():
    """Test GT reconstruction when all samples are reference."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [12345, 67890],
            "GEN_0__GT": ["0/0", "./."],
            "GEN_1__GT": ["0/0", "0/0"],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert "GT" in result.columns
    assert result["GT"].iloc[0] == ""
    assert result["GT"].iloc[1] == ""


@pytest.mark.unit
def test_missing_values():
    """Test GT reconstruction with NA/NaN/empty values."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2", "chr3"],
            "POS": [111, 222, 333],
            "GEN_0__GT": ["0/1", pd.NA, ""],
            "GEN_1__GT": [None, "1/1", "."],
            "GEN_2__GT": ["1/0", "0/0", "0/1"],
        }
    )
    vcf_samples = ["Sample1", "Sample2", "Sample3"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert result["GT"].iloc[0] == "Sample1(0/1);Sample3(1/0)"
    assert result["GT"].iloc[1] == "Sample2(1/1)"
    assert result["GT"].iloc[2] == "Sample3(0/1)"


@pytest.mark.unit
def test_single_sample():
    """Test GT reconstruction with only one sample."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [12345, 67890],
            "GEN_0__GT": ["0/1", "1/1"],
        }
    )
    vcf_samples = ["SampleA"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert result["GT"].iloc[0] == "SampleA(0/1)"
    assert result["GT"].iloc[1] == "SampleA(1/1)"


@pytest.mark.unit
def test_many_samples():
    """Test GT reconstruction with 12 samples."""
    samples = [f"Sample{i:02d}" for i in range(1, 13)]
    data = {
        "CHROM": ["chr1"],
        "POS": [12345],
    }
    # Samples 1, 3, 5, 7, 9, 11 have variants; others are 0/0 or ./.
    gts = ["0/1", "0/0", "1/1", "0/0", "0/1", "./.", "1/0", "0/0", "1/1", "0/0", "0/1", "0/0"]
    for i, gt in enumerate(gts):
        data[f"GEN_{i}__GT"] = [gt]

    df = pd.DataFrame(data)
    result = reconstruct_gt_column(df, samples)

    expected = "Sample01(0/1);Sample03(1/1);Sample05(0/1);Sample07(1/0);Sample09(1/1);Sample11(0/1)"
    assert result["GT"].iloc[0] == expected


@pytest.mark.unit
def test_phased_genotypes():
    """Test GT reconstruction with phased genotypes (|)."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2", "chr3"],
            "POS": [12345, 67890, 11111],
            "GEN_0__GT": ["0|1", "1|0", "0|0"],
            "GEN_1__GT": ["0|0", "1|1", ".|."],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert result["GT"].iloc[0] == "Sample1(0|1)"
    assert result["GT"].iloc[1] == "Sample1(1|0);Sample2(1|1)"
    assert result["GT"].iloc[2] == ""


@pytest.mark.unit
def test_hemizygous_genotypes():
    """Test GT reconstruction with hemizygous genotypes (single allele)."""
    df = pd.DataFrame(
        {
            "CHROM": ["chrX", "chrX"],
            "POS": [12345, 67890],
            "GEN_0__GT": ["1", "0"],
            "GEN_1__GT": ["0/1", "0/0"],
        }
    )
    vcf_samples = ["MaleSample", "FemaleSample"]

    result = reconstruct_gt_column(df, vcf_samples)

    # "0" is not in skip list, so hemizygous ref is included
    assert result["GT"].iloc[0] == "MaleSample(1);FemaleSample(0/1)"
    assert result["GT"].iloc[1] == "MaleSample(0)"


@pytest.mark.unit
def test_per_sample_columns_dropped():
    """Test that per-sample columns are removed after reconstruction."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [12345],
            "REF": ["A"],
            "ALT": ["G"],
            "GEN_0__GT": ["0/1"],
            "GEN_1__GT": ["1/1"],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert "CHROM" in result.columns
    assert "POS" in result.columns
    assert "REF" in result.columns
    assert "ALT" in result.columns
    assert "GT" in result.columns
    assert "GEN_0__GT" not in result.columns
    assert "GEN_1__GT" not in result.columns


@pytest.mark.unit
def test_column_ordering():
    """Test that GT column is added at the end."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [12345],
            "REF": ["A"],
            "ALT": ["G"],
            "GEN_0__GT": ["0/1"],
        }
    )
    vcf_samples = ["Sample1"]

    result = reconstruct_gt_column(df, vcf_samples)

    columns = result.columns.tolist()
    assert columns == ["CHROM", "POS", "REF", "ALT", "GT"]


@pytest.mark.unit
def test_empty_dataframe():
    """Test GT reconstruction with zero rows."""
    df = pd.DataFrame(
        {
            "CHROM": pd.Series(dtype=str),
            "POS": pd.Series(dtype=int),
            "GEN_0__GT": pd.Series(dtype=str),
            "GEN_1__GT": pd.Series(dtype=str),
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert "GT" in result.columns
    assert len(result) == 0
    assert "GEN_0__GT" not in result.columns
    assert "GEN_1__GT" not in result.columns


@pytest.mark.unit
def test_original_column_names():
    """Test GT reconstruction with unsanitized GEN[N].GT column names."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "GEN[0].GT": ["0/1"],
            "GEN[1].GT": ["1/1"],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert result["GT"].iloc[0] == "Sample1(0/1);Sample2(1/1)"
    assert "GEN[0].GT" not in result.columns
    assert "GEN[1].GT" not in result.columns


@pytest.mark.unit
def test_sample_count_mismatch_warning(caplog):
    """Test warning when GT column count differs from sample count."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "GEN_0__GT": ["0/1"],
            "GEN_1__GT": ["1/1"],
        }
    )
    # More samples than GT columns
    vcf_samples = ["Sample1", "Sample2", "Sample3"]

    result = reconstruct_gt_column(df, vcf_samples)

    assert "GT column count (2) != sample count (3)" in caplog.text
    # Should still work with available columns
    assert result["GT"].iloc[0] == "Sample1(0/1);Sample2(1/1)"
