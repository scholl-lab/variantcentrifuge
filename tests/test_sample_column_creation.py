"""Tests for the unified sample column creation function."""

import pandas as pd
import pytest

from variantcentrifuge.stages.analysis_stages import (
    create_sample_columns_from_gt,
    create_sample_columns_from_gt_vectorized,
    create_sample_columns_from_gt_intelligent,
)


class TestSampleColumnCreation:
    """Test the unified sample column creation function."""

    def test_replaced_genotype_format(self):
        """Test sample column creation from replaced genotype format."""
        # Test data with replaced genotype format
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1"],
                "POS": [100, 200],
                "GT": [
                    "Sample1(0/1);Sample2(0/0);Sample3(1/1)",
                    "Sample1(1/0);Sample2(1/1);Sample3(0/1)",
                ],
            }
        )

        vcf_samples = ["Sample1", "Sample2", "Sample3"]

        result = create_sample_columns_from_gt(df, vcf_samples, separator=";")

        # Check that sample columns were created
        assert "Sample1" in result.columns
        assert "Sample2" in result.columns
        assert "Sample3" in result.columns

        # Check genotype values
        assert result.iloc[0]["Sample1"] == "0/1"
        assert result.iloc[0]["Sample2"] == "0/0"
        assert result.iloc[0]["Sample3"] == "1/1"

        assert result.iloc[1]["Sample1"] == "1/0"
        assert result.iloc[1]["Sample2"] == "1/1"
        assert result.iloc[1]["Sample3"] == "0/1"

    def test_snpsift_format(self):
        """Test sample column creation from SnpSift format."""
        # Test data with SnpSift format
        df = pd.DataFrame(
            {"CHROM": ["chr1", "chr1"], "POS": [100, 200], "GT": ["0/1,0/0,1/1", "1/0,1/1,0/1"]}
        )

        vcf_samples = ["Sample1", "Sample2", "Sample3"]

        result = create_sample_columns_from_gt(df, vcf_samples, snpsift_sep=",")

        # Check that sample columns were created
        assert "Sample1" in result.columns
        assert "Sample2" in result.columns
        assert "Sample3" in result.columns

        # Check genotype values
        assert result.iloc[0]["Sample1"] == "0/1"
        assert result.iloc[0]["Sample2"] == "0/0"
        assert result.iloc[0]["Sample3"] == "1/1"

        assert result.iloc[1]["Sample1"] == "1/0"
        assert result.iloc[1]["Sample2"] == "1/1"
        assert result.iloc[1]["Sample3"] == "0/1"

    def test_empty_dataframe(self):
        """Test handling of empty DataFrame."""
        df = pd.DataFrame({"GT": []})
        vcf_samples = ["Sample1", "Sample2"]

        result = create_sample_columns_from_gt(df, vcf_samples)

        # Should return original DataFrame without error
        assert len(result) == 0
        assert "GT" in result.columns

    def test_missing_gt_column(self):
        """Test error handling when GT column is missing."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": [100]})
        vcf_samples = ["Sample1"]

        with pytest.raises(ValueError, match="GT column not found"):
            create_sample_columns_from_gt(df, vcf_samples)

    def test_sample_columns_already_exist(self):
        """Test that existing sample columns are not overwritten."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": [100],
                "GT": ["Sample1(0/1);Sample2(0/0)"],
                "Sample1": ["existing_value"],
                "Sample2": ["existing_value"],
            }
        )
        vcf_samples = ["Sample1", "Sample2"]

        result = create_sample_columns_from_gt(df, vcf_samples)

        # Sample columns should remain unchanged
        assert result.iloc[0]["Sample1"] == "existing_value"
        assert result.iloc[0]["Sample2"] == "existing_value"

    def test_missing_genotype_data(self):
        """Test handling of missing genotype data."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1"],
                "POS": [100, 200],
                "GT": ["Sample1(0/1);Sample2(0/0)", "NA"],
            }
        )
        vcf_samples = ["Sample1", "Sample2"]

        result = create_sample_columns_from_gt(df, vcf_samples)

        # First row should have genotypes
        assert result.iloc[0]["Sample1"] == "0/1"
        assert result.iloc[0]["Sample2"] == "0/0"

        # Second row should have missing genotypes
        assert result.iloc[1]["Sample1"] == "./."
        assert result.iloc[1]["Sample2"] == "./."

    def test_custom_separator(self):
        """Test with custom separator."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": [100], "GT": ["Sample1(0/1)|Sample2(0/0)"]})
        vcf_samples = ["Sample1", "Sample2"]

        result = create_sample_columns_from_gt(df, vcf_samples, separator="|")

        assert result.iloc[0]["Sample1"] == "0/1"
        assert result.iloc[0]["Sample2"] == "0/0"

    def test_unrecognized_format(self):
        """Test handling of unrecognized GT format."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": [100], "GT": ["some_unrecognized_format"]})
        vcf_samples = ["Sample1", "Sample2"]

        # Should not raise an error, just log a warning
        result = create_sample_columns_from_gt(df, vcf_samples)

        # Sample columns should not be created
        assert "Sample1" not in result.columns
        assert "Sample2" not in result.columns


class TestVectorizedSampleColumnCreation:
    """Test the vectorized sample column creation function."""

    def test_vectorized_vs_iterative_equivalence_replaced_format(self):
        """Test that vectorized and iterative methods produce identical results for replaced format."""
        # Test data with replaced genotype format
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1", "chr2"],
                "POS": [100, 200, 300],
                "GT": [
                    "Sample1(0/1);Sample2(0/0);Sample3(1/1)",
                    "Sample1(1/0);Sample2(1/1);Sample3(0/1)",
                    "Sample1(./.);Sample2(0/1);Sample3(1/0)",
                ],
            }
        )

        vcf_samples = ["Sample1", "Sample2", "Sample3"]

        # Process with both methods
        result_iterative = create_sample_columns_from_gt(df.copy(), vcf_samples, separator=";")
        result_vectorized = create_sample_columns_from_gt_vectorized(
            df.copy(), vcf_samples, separator=";"
        )

        # Results should be identical
        assert list(result_iterative.columns) == list(result_vectorized.columns)

        for sample in vcf_samples:
            pd.testing.assert_series_equal(
                result_iterative[sample], result_vectorized[sample], check_names=False
            )

    def test_vectorized_vs_iterative_equivalence_snpsift_format(self):
        """Test that vectorized and iterative methods produce identical results for SnpSift format."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2"],
                "POS": [100, 200],
                "GT": ["0/1,0/0,1/1", "1/0,1/1,0/1"],
            }
        )

        vcf_samples = ["Sample1", "Sample2", "Sample3"]

        # Process with both methods
        result_iterative = create_sample_columns_from_gt(df.copy(), vcf_samples, snpsift_sep=",")
        result_vectorized = create_sample_columns_from_gt_vectorized(
            df.copy(), vcf_samples, snpsift_sep=","
        )

        # Results should be identical
        assert list(result_iterative.columns) == list(result_vectorized.columns)

        for sample in vcf_samples:
            pd.testing.assert_series_equal(
                result_iterative[sample], result_vectorized[sample], check_names=False
            )

    def test_vectorized_performance_large_dataset(self):
        """Test vectorized method performance on larger dataset."""
        import time

        # Create larger test dataset
        num_variants = 1000
        num_samples = 100
        sample_names = [f"Sample{i}" for i in range(1, num_samples + 1)]

        # Generate realistic GT data
        gt_data = []
        for v in range(num_variants):
            sample_entries = []
            for s in range(num_samples):
                genotypes = ["0/0", "0/1", "1/1", "./."]
                gt = genotypes[v % len(genotypes)]  # Cycle through genotypes
                sample_entries.append(f"Sample{s+1}({gt})")
            gt_data.append(";".join(sample_entries))

        df = pd.DataFrame(
            {
                "CHROM": [f"chr{(i % 22) + 1}" for i in range(num_variants)],
                "POS": list(range(1000, 1000 + num_variants)),
                "GT": gt_data,
            }
        )

        # Time vectorized method
        start_time = time.time()
        result_vectorized = create_sample_columns_from_gt_vectorized(
            df.copy(), sample_names, separator=";"
        )
        vectorized_time = time.time() - start_time

        # Verify results
        assert len(result_vectorized) == num_variants
        for sample in sample_names:
            assert sample in result_vectorized.columns

        # Check some specific values
        assert result_vectorized.iloc[0]["Sample1"] == "0/0"  # First variant, first sample
        assert result_vectorized.iloc[1]["Sample1"] == "0/1"  # Second variant, first sample

        print(
            f"Vectorized method processed {num_variants:,} variants × {num_samples:,} samples in {vectorized_time:.3f}s"
        )

        # Should be reasonably fast (this is a performance indicator, not a strict requirement)
        ops_per_second = (
            (num_variants * num_samples) / vectorized_time if vectorized_time > 0 else 0
        )
        print(f"Performance: {ops_per_second:,.0f} operations/second")


class TestIntelligentSampleColumnCreation:
    """Test the intelligent sample column creation function."""

    def test_auto_selection_small_dataset(self):
        """Test that auto selection chooses iterative for small datasets."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": [100],
                "GT": ["Sample1(0/1);Sample2(0/0)"],
            }
        )
        vcf_samples = ["Sample1", "Sample2"]

        # Should work with auto selection
        result = create_sample_columns_from_gt_intelligent(df, vcf_samples, method="auto")

        assert "Sample1" in result.columns
        assert "Sample2" in result.columns
        assert result.iloc[0]["Sample1"] == "0/1"
        assert result.iloc[0]["Sample2"] == "0/0"

    def test_explicit_method_selection(self):
        """Test explicit method selection works."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": [100],
                "GT": ["Sample1(0/1);Sample2(0/0)"],
            }
        )
        vcf_samples = ["Sample1", "Sample2"]

        # Test explicit iterative selection
        result_iterative = create_sample_columns_from_gt_intelligent(
            df.copy(), vcf_samples, method="iterative"
        )

        # Test explicit vectorized selection
        result_vectorized = create_sample_columns_from_gt_intelligent(
            df.copy(), vcf_samples, method="vectorized"
        )

        # Both should produce same results
        for sample in vcf_samples:
            assert result_iterative.iloc[0][sample] == result_vectorized.iloc[0][sample]


if __name__ == "__main__":
    # Run tests with logging
    import logging

    logging.basicConfig(level=logging.DEBUG)
    pytest.main([__file__, "-v"])
