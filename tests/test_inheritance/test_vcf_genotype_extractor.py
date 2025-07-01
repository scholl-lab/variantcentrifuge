"""Test cases for VCF genotype extraction module."""

import pandas as pd
import pytest
from unittest.mock import patch, MagicMock
from variantcentrifuge.vcf_genotype_extractor import (
    extract_genotypes_from_vcf,
    create_variant_genotype_map,
    augment_dataframe_with_vcf_genotypes,
)


class TestGenotypeExtraction:
    """Test genotype extraction from VCF files."""

    @patch("subprocess.run")
    def test_extract_genotypes_basic(self, mock_run):
        """Test basic genotype extraction from VCF."""
        # Mock bcftools output
        mock_run.return_value = MagicMock(
            stdout="1\t1000\tA\tT\t0/1\t0/0\t0/1\n2\t2000\tC\tG\t1/1\t0/1\t0/1\n",
            stderr="",
            returncode=0,
        )

        variant_positions = [("1", "1000", "A", "T"), ("2", "2000", "C", "G")]
        sample_list = ["Sample1", "Sample2", "Sample3"]

        result = extract_genotypes_from_vcf("test.vcf", variant_positions, sample_list)

        # Check results
        assert "1:1000:A>T" in result
        assert result["1:1000:A>T"]["Sample1"] == "0/1"
        assert result["1:1000:A>T"]["Sample2"] == "0/0"
        assert result["1:1000:A>T"]["Sample3"] == "0/1"

        assert "2:2000:C>G" in result
        assert result["2:2000:C>G"]["Sample1"] == "1/1"
        assert result["2:2000:C>G"]["Sample2"] == "0/1"
        assert result["2:2000:C>G"]["Sample3"] == "0/1"

    @patch("subprocess.run")
    def test_extract_genotypes_empty(self, mock_run):
        """Test extraction with no variants."""
        variant_positions = []
        sample_list = ["Sample1"]

        result = extract_genotypes_from_vcf("test.vcf", variant_positions, sample_list)

        assert result == {}
        # Should not call subprocess if no variants
        mock_run.assert_not_called()

    @patch("subprocess.run")
    def test_extract_genotypes_error(self, mock_run):
        """Test handling of bcftools error."""
        mock_run.side_effect = Exception("bcftools not found")

        variant_positions = [("1", "1000", "A", "T")]
        sample_list = ["Sample1"]

        # The function will raise the exception, not return empty dict
        with pytest.raises(Exception, match="bcftools not found"):
            extract_genotypes_from_vcf("test.vcf", variant_positions, sample_list)

    @patch("subprocess.run")
    def test_extract_genotypes_missing_samples(self, mock_run):
        """Test extraction with fewer genotypes than samples."""
        # Only 2 genotypes but 3 samples
        mock_run.return_value = MagicMock(
            stdout="1\t1000\tA\tT\t0/1\t0/0\n", stderr="", returncode=0
        )

        variant_positions = [("1", "1000", "A", "T")]
        sample_list = ["Sample1", "Sample2", "Sample3"]

        result = extract_genotypes_from_vcf("test.vcf", variant_positions, sample_list)

        assert result["1:1000:A>T"]["Sample1"] == "0/1"
        assert result["1:1000:A>T"]["Sample2"] == "0/0"
        # Sample3 should not be in the result (no genotype data)
        assert "Sample3" not in result["1:1000:A>T"]


class TestVariantGenotypeMapping:
    """Test creating variant to genotype mappings."""

    @patch("variantcentrifuge.vcf_genotype_extractor.extract_genotypes_from_vcf")
    def test_create_variant_genotype_map(self, mock_extract):
        """Test creating genotype map from DataFrame."""
        # Mock genotype extraction
        mock_extract.return_value = {
            "1:1000:A>T": {"Sample1": "0/1", "Sample2": "0/0"},
            "2:2000:C>G": {"Sample1": "1/1", "Sample2": "0/1"},
        }

        df = pd.DataFrame(
            [
                {"CHROM": "1", "POS": "1000", "REF": "A", "ALT": "T"},
                {"CHROM": "2", "POS": "2000", "REF": "C", "ALT": "G"},
                {"CHROM": "3", "POS": "3000", "REF": "G", "ALT": "A"},  # Not in VCF
            ]
        )

        sample_list = ["Sample1", "Sample2"]

        result = create_variant_genotype_map(df, "test.vcf", sample_list)

        # Check mapping
        assert 0 in result  # First row
        assert result[0] == {"Sample1": "0/1", "Sample2": "0/0"}

        assert 1 in result  # Second row
        assert result[1] == {"Sample1": "1/1", "Sample2": "0/1"}

        assert 2 in result  # Third row
        assert result[2] == {}  # No genotype data

    @patch("variantcentrifuge.vcf_genotype_extractor.extract_genotypes_from_vcf")
    def test_create_variant_genotype_map_missing_columns(self, mock_extract):
        """Test handling of missing required columns."""
        mock_extract.return_value = {}

        # Missing POS column
        df = pd.DataFrame(
            [
                {"CHROM": "1", "REF": "A", "ALT": "T"},
            ]
        )

        result = create_variant_genotype_map(df, "test.vcf", ["Sample1"])

        # Should handle gracefully
        assert 0 in result
        assert result[0] == {}


class TestDataFrameAugmentation:
    """Test augmenting DataFrame with VCF genotypes."""

    @patch("variantcentrifuge.vcf_genotype_extractor.create_variant_genotype_map")
    def test_augment_dataframe_basic(self, mock_create_map):
        """Test basic DataFrame augmentation."""
        mock_create_map.return_value = {
            0: {"Sample1": "0/1", "Sample2": "0/0"},
            1: {"Sample1": "1/1", "Sample2": "0/1"},
            2: {"Sample1": "./.", "Sample2": "./."},
        }

        df = pd.DataFrame(
            [
                {"CHROM": "1", "POS": "1000", "REF": "A", "ALT": "T"},
                {"CHROM": "2", "POS": "2000", "REF": "C", "ALT": "G"},
                {"CHROM": "3", "POS": "3000", "REF": "G", "ALT": "A"},
            ]
        )

        sample_list = ["Sample1", "Sample2"]

        result_df = augment_dataframe_with_vcf_genotypes(df, "test.vcf", sample_list)

        # Check that sample columns were added
        assert "Sample1" in result_df.columns
        assert "Sample2" in result_df.columns

        # Check genotype values
        assert result_df.loc[0, "Sample1"] == "0/1"
        assert result_df.loc[0, "Sample2"] == "0/0"
        assert result_df.loc[1, "Sample1"] == "1/1"
        assert result_df.loc[1, "Sample2"] == "0/1"
        assert result_df.loc[2, "Sample1"] == "./."
        assert result_df.loc[2, "Sample2"] == "./."

    @patch("variantcentrifuge.vcf_genotype_extractor.create_variant_genotype_map")
    def test_augment_dataframe_existing_columns(self, mock_create_map):
        """Test augmentation when sample columns already exist."""
        mock_create_map.return_value = {
            0: {"Sample1": "0/1"},
        }

        # DataFrame already has Sample1 column
        df = pd.DataFrame(
            [
                {"CHROM": "1", "POS": "1000", "Sample1": "OLD_VALUE"},
            ]
        )

        result_df = augment_dataframe_with_vcf_genotypes(df, "test.vcf", ["Sample1"])

        # Should overwrite existing column
        assert result_df.loc[0, "Sample1"] == "0/1"

    @patch("variantcentrifuge.vcf_genotype_extractor.create_variant_genotype_map")
    def test_augment_dataframe_partial_data(self, mock_create_map):
        """Test augmentation with partial genotype data."""
        # Only Sample1 has data for first variant
        mock_create_map.return_value = {
            0: {"Sample1": "0/1"},  # No Sample2 data
            1: {"Sample2": "0/1"},  # No Sample1 data
        }

        df = pd.DataFrame(
            [
                {"CHROM": "1", "POS": "1000"},
                {"CHROM": "2", "POS": "2000"},
            ]
        )

        result_df = augment_dataframe_with_vcf_genotypes(df, "test.vcf", ["Sample1", "Sample2"])

        # Check missing data is filled with ./.
        assert result_df.loc[0, "Sample1"] == "0/1"
        assert result_df.loc[0, "Sample2"] == "./."
        assert result_df.loc[1, "Sample1"] == "./."
        assert result_df.loc[1, "Sample2"] == "0/1"


class TestIntegration:
    """Integration tests with real-like data."""

    def test_full_workflow_simulation(self):
        """Simulate full workflow with mocked subprocess."""
        with patch("subprocess.run") as mock_run:
            # Simulate bcftools output
            mock_run.return_value = MagicMock(
                stdout="1\t1000\tA\tT\t0/1\n1\t2000\tC\tG\t0/1\n", stderr="", returncode=0
            )

            # Create DataFrame
            df = pd.DataFrame(
                [
                    {"CHROM": "1", "POS": "1000", "REF": "A", "ALT": "T", "GENE": "GENE1"},
                    {"CHROM": "1", "POS": "2000", "REF": "C", "ALT": "G", "GENE": "GENE1"},
                ]
            )

            # Augment with genotypes
            result_df = augment_dataframe_with_vcf_genotypes(df, "test.vcf", ["Sample1"])

            # Verify
            assert "Sample1" in result_df.columns
            assert all(result_df["Sample1"] == "0/1")

            # Both variants should have genotypes for compound het detection
            assert len(result_df[result_df["Sample1"] == "0/1"]) == 2
