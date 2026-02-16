"""
Tests for full Phase 11 pipeline data flow (bcftools query -> analysis -> GT reconstruction).

These tests prove the new pipeline path produces identical output to the old path.
"""

import pandas as pd
import pytest

from variantcentrifuge.stages.analysis_stages import create_sample_columns_from_gt_intelligent
from variantcentrifuge.stages.output_stages import reconstruct_gt_column
from variantcentrifuge.stages.processing_stages import GenotypeReplacementStage


@pytest.mark.unit
def test_genotype_replacement_stage_noop():
    """Test that GenotypeReplacementStage is now a no-op (Phase 11)."""
    # GenotypeReplacementStage._process just returns context immediately
    # We can verify this by checking the code path
    import inspect

    # Read the _process method source
    source = inspect.getsource(GenotypeReplacementStage._process)

    # Verify it contains the Phase 11 skip message at the beginning
    assert "Phase 11" in source
    skip_msg = "genotype replacement skipped"
    assert skip_msg.title() in source or skip_msg in source.lower()

    # Verify it returns context immediately (check that return statement is near the top)
    lines = source.split("\n")
    # Find the line with "return context"
    return_line_idx = None
    for i, line in enumerate(lines):
        if "return context" in line and not line.strip().startswith("#"):
            return_line_idx = i
            break

    # The return should be within the first ~10 lines (after docstring and log)
    assert return_line_idx is not None, "Should have 'return context' statement"
    err_msg = f"return context should be near top, found at line {return_line_idx}"
    assert return_line_idx < 10, err_msg


@pytest.mark.unit
def test_full_data_flow_simulation():
    """
    Test full Phase 11 data flow: bcftools query -> analysis -> GT reconstruction.

    This simulates the complete pipeline with per-sample columns from bcftools.
    """
    # 1. Simulate bcftools query output (per-sample GT columns, sanitized names)
    # After DataFrame loading, GEN[0].GT -> GEN_0__GT (sanitized for itertuples)
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2", "chr3"],
            "POS": [12345, 67890, 11111],
            "REF": ["A", "C", "G"],
            "ALT": ["G", "T", "A"],
            "ANN_0__GENE": ["GENE1", "GENE2", "GENE3"],
            "ANN_0__EFFECT": ["missense_variant", "frameshift_variant", "stop_gained"],
            "GEN_0__GT": ["0/1", "1/1", "0/0"],  # Per-sample column from bcftools
            "GEN_1__GT": ["1/1", "0/1", "0/1"],  # Per-sample column from bcftools
            "GEN_2__GT": ["0/0", "./.", "1/1"],  # Per-sample column from bcftools
        }
    )

    vcf_samples = ["Sample1", "Sample2", "Sample3"]

    # 2. In Phase 11, per-sample columns are GEN_N__GT (not sample-named).
    # create_sample_columns_from_gt_intelligent is not called because
    # InheritanceAnalysisStage works directly with per-sample GT columns.

    # 3. Simulate analysis (inheritance analysis would work with per-sample columns)
    # Just add a mock analysis column
    df["Inheritance_Pattern"] = ["AD", "AR", "de_novo"]

    # 4. Reconstruct GT column for output
    output_df = reconstruct_gt_column(df, vcf_samples)

    # 5. Verify GT column is correctly reconstructed
    assert "GT" in output_df.columns
    # Row 0: Sample1=0/1, Sample2=1/1 (skip Sample3=0/0)
    assert output_df["GT"].iloc[0] == "Sample1(0/1);Sample2(1/1)"
    # Row 1: Sample1=1/1, Sample2=0/1 (skip Sample3=./.)
    assert output_df["GT"].iloc[1] == "Sample1(1/1);Sample2(0/1)"
    # Row 2: Sample1=0/0 (skip), Sample2=0/1, Sample3=1/1
    assert output_df["GT"].iloc[2] == "Sample2(0/1);Sample3(1/1)"

    # 6. Verify per-sample columns are dropped
    assert "GEN_0__GT" not in output_df.columns
    assert "GEN_1__GT" not in output_df.columns
    assert "GEN_2__GT" not in output_df.columns

    # 7. Verify output has expected column structure (same as old pipeline)
    expected_columns = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "ANN_0__GENE",
        "ANN_0__EFFECT",
        "Inheritance_Pattern",
        "GT",
    ]
    assert output_df.columns.tolist() == expected_columns


@pytest.mark.unit
def test_backwards_compatibility_packed_gt():
    """
    Test backwards compatibility with old-style packed GT column.

    If a DataFrame with old-style packed GT is loaded (e.g., from checkpoint),
    the pipeline should still work.
    """
    # Old-style DataFrame with packed GT column
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [12345, 67890],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
            "GT": ["Sample1(0/1);Sample2(1/1)", "Sample1(1/1)"],
        }
    )

    vcf_samples = ["Sample1", "Sample2", "Sample3"]

    # create_sample_columns_from_gt_intelligent should create per-sample columns
    result_df = create_sample_columns_from_gt_intelligent(df, vcf_samples)

    # Should have created per-sample columns
    assert "Sample1" in result_df.columns
    assert "Sample2" in result_df.columns
    assert "Sample3" in result_df.columns

    # Verify values (samples not in GT string get "./.")
    assert result_df["Sample1"].iloc[0] == "0/1"
    assert result_df["Sample2"].iloc[0] == "1/1"
    assert result_df["Sample3"].iloc[0] == "./."  # Not in GT string

    assert result_df["Sample1"].iloc[1] == "1/1"
    assert result_df["Sample2"].iloc[1] == "./."  # Not in GT string
    assert result_df["Sample3"].iloc[1] == "./."  # Not in GT string


@pytest.mark.unit
def test_sample_columns_detection():
    """Test that sample column detection works correctly."""
    # Test case 1: Per-sample columns exist (Phase 11 new pipeline)
    # In this case, we already have per-sample columns from bcftools query
    df1 = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [12345],
            "Sample1": ["0/1"],
            "Sample2": ["1/1"],
            "GT": ["Sample1(0/1);Sample2(1/1)"],  # May or may not exist
        }
    )
    vcf_samples1 = ["Sample1", "Sample2"]

    result1 = create_sample_columns_from_gt_intelligent(df1, vcf_samples1)
    # Should return unchanged because sample columns exist
    assert "Sample1" in result1.columns
    assert "Sample2" in result1.columns

    # Test case 2: Per-sample columns don't exist, but GT exists (backwards compatibility)
    df2 = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [12345],
            "GT": ["Sample1(0/1);Sample2(1/1)"],
        }
    )
    vcf_samples2 = ["Sample1", "Sample2"]

    result2 = create_sample_columns_from_gt_intelligent(df2, vcf_samples2)
    # Should have created sample columns from GT
    assert "Sample1" in result2.columns
    assert "Sample2" in result2.columns


@pytest.mark.unit
def test_output_column_structure_equivalence():
    """
    Test that output column structure matches old pipeline.

    Old pipeline: many intermediate columns + GT column
    New pipeline: same columns, GT reconstructed at output time
    """
    # Create DataFrame with per-sample GT columns (new pipeline)
    df_new = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [12345],
            "REF": ["A"],
            "ALT": ["G"],
            "ANN_0__GENE": ["BRCA1"],
            "GEN_0__GT": ["0/1"],
            "GEN_1__GT": ["1/1"],
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    # Reconstruct GT
    output_new = reconstruct_gt_column(df_new, vcf_samples)

    # Create DataFrame with packed GT (old pipeline)
    df_old = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [12345],
            "REF": ["A"],
            "ALT": ["G"],
            "ANN_0__GENE": ["BRCA1"],
            "GT": ["Sample1(0/1);Sample2(1/1)"],
        }
    )

    # Column sets should match (order may differ)
    assert set(output_new.columns) == set(df_old.columns)

    # GT column content should match
    assert output_new["GT"].iloc[0] == df_old["GT"].iloc[0]


@pytest.mark.unit
def test_empty_dataframe_handling():
    """Test that pipeline handles empty DataFrames gracefully."""
    df = pd.DataFrame(
        {
            "CHROM": pd.Series(dtype=str),
            "POS": pd.Series(dtype=int),
            "GEN_0__GT": pd.Series(dtype=str),
            "GEN_1__GT": pd.Series(dtype=str),
        }
    )
    vcf_samples = ["Sample1", "Sample2"]

    # Should handle empty DataFrame
    result = reconstruct_gt_column(df, vcf_samples)

    assert "GT" in result.columns
    assert len(result) == 0
    assert "GEN_0__GT" not in result.columns
    assert "GEN_1__GT" not in result.columns
