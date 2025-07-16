"""
Tests for the unified sample column creation function.
"""

import pandas as pd
import pytest

from variantcentrifuge.stages.analysis_stages import create_sample_columns_from_gt


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


if __name__ == "__main__":
    # Run tests with logging
    import logging

    logging.basicConfig(level=logging.DEBUG)
    pytest.main([__file__, "-v"])
