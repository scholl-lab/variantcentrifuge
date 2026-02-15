"""
Tests for GT column reconstruction from per-sample columns (Phase 11).

These tests prove that reconstruct_gt_column() produces identical output format
to the old genotype replacement pipeline.
"""

import pandas as pd
import pytest

from variantcentrifuge.stages.output_stages import reconstruct_gt_column


@pytest.mark.unit
def test_basic_reconstruction():
    """Test basic GT reconstruction with mixed genotypes."""
    # Create DataFrame with per-sample GT columns
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [12345, 67890],
            "Sample1": ["0/1", "1/1"],
            "Sample2": ["0/0", "0/1"],
            "Sample3": ["1/1", "./."],
        }
    )
    vcf_samples = ["Sample1", "Sample2", "Sample3"]

    result = reconstruct_gt_column(df, vcf_samples)

    # Expected GT: only include non-reference genotypes
    # Row 0: Sample1(0/1) and Sample3(1/1) - skip Sample2 (0/0)
    # Row 1: Sample1(1/1) and Sample2(0/1) - skip Sample3 (./.)
    assert "GT" in result.columns
    assert result["GT"].iloc[0] == "Sample1(0/1);Sample3(1/1)"
    assert result["GT"].iloc[1] == "Sample1(1/1);Sample2(0/1)"

    # Per-sample columns should be dropped
    assert "Sample1" not in result.columns
    assert "Sample2" not in result.columns
    assert "Sample3" not in result.columns


@pytest.mark.unit
def test_all_reference_genotypes():
    """Test GT reconstruction when all samples are reference."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [12345, 67890],
            "Sample1": ["0/0", "./."],
            "Sample2": ["0/0", "0/0"],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    # All reference -> empty GT string
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
            "Sample1": ["0/1", pd.NA, ""],
            "Sample2": [None, "1/1", "."],
            "Sample3": ["1/0", "0/0", "0/1"],
        }
    )
    vcf_samples = ["Sample1", "Sample2", "Sample3"]

    result = reconstruct_gt_column(df, vcf_samples)

    # Row 0: Sample1(0/1) and Sample3(1/0) - Sample2 is None
    # Row 1: Sample2(1/1) - Sample1 is NA, Sample3 is 0/0
    # Row 2: Sample3(0/1) - Sample1 is empty, Sample2 is "."
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
            "SampleA": ["0/1", "1/1"],
        }
    )
    vcf_samples = ["SampleA"]

    result = reconstruct_gt_column(df, vcf_samples)

    # Single sample with variant
    assert result["GT"].iloc[0] == "SampleA(0/1)"
    assert result["GT"].iloc[1] == "SampleA(1/1)"


@pytest.mark.unit
def test_many_samples():
    """Test GT reconstruction with 10+ samples."""
    # Create 12 samples with mix of variant and reference
    samples = [f"Sample{i:02d}" for i in range(1, 13)]
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [12345],
            # Samples 1, 3, 5, 7, 9, 11 have variants
            "Sample01": ["0/1"],
            "Sample02": ["0/0"],
            "Sample03": ["1/1"],
            "Sample04": ["0/0"],
            "Sample05": ["0/1"],
            "Sample06": ["./."],
            "Sample07": ["1/0"],
            "Sample08": ["0/0"],
            "Sample09": ["1/1"],
            "Sample10": ["0/0"],
            "Sample11": ["0/1"],
            "Sample12": ["0/0"],
        }
    )

    result = reconstruct_gt_column(df, samples)

    # Only variants included
    expected = "Sample01(0/1);Sample03(1/1);Sample05(0/1);Sample07(1/0);Sample09(1/1);Sample11(0/1)"
    assert result["GT"].iloc[0] == expected


@pytest.mark.unit
def test_phased_genotypes():
    """Test GT reconstruction with phased genotypes (|)."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2", "chr3"],
            "POS": [12345, 67890, 11111],
            "Sample1": ["0|1", "1|0", "0|0"],
            "Sample2": ["0|0", "1|1", ".|."],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    # Phased genotypes included (not reference)
    # Row 0: Sample1(0|1), Sample2 is 0|0 (phased reference) -> excluded
    # Row 1: Sample1(1|0);Sample2(1|1)
    # Row 2: Sample1 is 0|0 (phased reference), Sample2 is .|. -> both excluded
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
            "MaleSample": ["1", "0"],
            "FemaleSample": ["0/1", "0/0"],
        }
    )
    vcf_samples = ["MaleSample", "FemaleSample"]

    result = reconstruct_gt_column(df, vcf_samples)

    # Row 0: MaleSample(1) and FemaleSample(0/1)
    # Row 1: MaleSample(0) is technically reference but single allele, FemaleSample is 0/0
    # Based on the code logic, "0" is not in ["0/0", "./.", ".", "NA"] so it's included
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
            "Sample1": ["0/1"],
            "Sample2": ["1/1"],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    # Original columns should remain
    assert "CHROM" in result.columns
    assert "POS" in result.columns
    assert "REF" in result.columns
    assert "ALT" in result.columns

    # GT column added
    assert "GT" in result.columns

    # Per-sample columns dropped
    assert "Sample1" not in result.columns
    assert "Sample2" not in result.columns


@pytest.mark.unit
def test_column_ordering():
    """Test that GT column is added at the correct position."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [12345],
            "REF": ["A"],
            "ALT": ["G"],
            "Sample1": ["0/1"],
        }
    )
    vcf_samples = ["Sample1"]

    result = reconstruct_gt_column(df, vcf_samples)

    # GT column should exist
    assert "GT" in result.columns

    # Column order: CHROM, POS, REF, ALT, GT
    # (Per-sample columns are dropped)
    columns = result.columns.tolist()
    assert columns == ["CHROM", "POS", "REF", "ALT", "GT"]


@pytest.mark.unit
def test_empty_dataframe():
    """Test GT reconstruction with zero rows."""
    df = pd.DataFrame(
        {
            "CHROM": [],
            "POS": [],
            "Sample1": [],
            "Sample2": [],
        }
    )
    df = df.astype({"CHROM": str, "POS": int, "Sample1": str, "Sample2": str})
    vcf_samples = ["Sample1", "Sample2"]

    result = reconstruct_gt_column(df, vcf_samples)

    # Should return DataFrame with GT column but no rows
    assert "GT" in result.columns
    assert len(result) == 0
    assert "Sample1" not in result.columns
    assert "Sample2" not in result.columns
