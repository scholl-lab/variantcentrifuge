"""
Test module for utility functions in variantcentrifuge/utils.py.

This module contains tests for various utility functions in the
variantcentrifuge.utils module, focusing on header normalization
for SnpEff and SnpSift output.
"""

from variantcentrifuge.utils import normalize_snpeff_headers, normalize_vcf_headers


class TestHeaderNormalization:
    """Tests for header normalization functions."""

    def test_indexed_field_renaming(self):
        """Test that indexed fields are properly renamed from GEN[index].FIELD to FIELD_index."""
        # Test header with SnpSift genotype fields that use indices
        test_header = "#CHROM\tPOS\tREF\tALT\tGEN[0].AF\tGEN[1].AF\tGEN[2].DP\tGEN[0].DP"
        expected = "#CHROM\tPOS\tREF\tALT\tAF_0\tAF_1\tDP_2\tDP_0"

        # Apply the normalization
        result = normalize_vcf_headers([test_header])

        # Verify the transformation
        assert result[0] == expected, f"Expected '{expected}', got '{result[0]}'"

    def test_standard_snpeff_prefix_removal(self):
        """Test that standard SnpEff prefixes are properly removed."""
        # Test header with standard SnpEff prefixes
        test_header = "#CHROM\tPOS\tREF\tALT\tANN[*].IMPACT\tANN[0].GENE\tGEN[*].AF"
        expected = "#CHROM\tPOS\tREF\tALT\tIMPACT\tGENE\tAF"

        # Apply the normalization
        result = normalize_vcf_headers([test_header])

        # Verify the transformation
        assert result[0] == expected, f"Expected '{expected}', got '{result[0]}'"

    def test_lof_nmd_prefix_transformation(self):
        """Test that LOF and NMD prefixes are properly transformed."""
        # Test header with LOF and NMD prefixes
        test_header = "#CHROM\tPOS\tREF\tALT\tLOF[*].GENE\tNMD[0].GENE"
        expected = "#CHROM\tPOS\tREF\tALT\tLOF_GENE\tNMD_GENE"

        # Apply the normalization
        result = normalize_vcf_headers([test_header])

        # Verify the transformation
        assert result[0] == expected, f"Expected '{expected}', got '{result[0]}'"

    def test_mixed_transformations(self):
        """Test a combination of different transformations."""
        # Test header with a mix of transformation types
        test_header = "#CHROM\tPOS\tREF\tALT\tGEN[0].AF\tANN[*].IMPACT\tLOF[*].GENE\tGEN[2].DP"
        expected = "#CHROM\tPOS\tREF\tALT\tAF_0\tIMPACT\tLOF_GENE\tDP_2"

        # Apply the normalization
        result = normalize_vcf_headers([test_header])

        # Verify the transformation
        assert result[0] == expected, f"Expected '{expected}', got '{result[0]}'"

    def test_empty_lines(self):
        """Test that empty input is handled gracefully."""
        # Test with empty input
        assert normalize_vcf_headers([]) == []

    def test_no_transformation_needed(self):
        """Test that headers without prefixes remain unchanged."""
        # Test with a header that doesn't need transformation
        test_header = "#CHROM\tPOS\tREF\tALT\tIMPACT\tGENE\tAF"

        # Apply the normalization
        result = normalize_vcf_headers([test_header])

        # Verify nothing changed
        assert result[0] == test_header

    def test_alias_function(self):
        """Test that the deprecated alias function produces the same result."""
        # Test header with various prefixes
        test_header = "#CHROM\tPOS\tREF\tALT\tGEN[0].AF\tANN[*].IMPACT"

        # Apply both normalization functions
        result1 = normalize_vcf_headers([test_header])
        result2 = normalize_snpeff_headers([test_header])

        # Verify they produce the same result
        assert result1 == result2, "Alias function should produce the same result"
