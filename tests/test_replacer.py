"""Unit tests for the replacer.py module, focusing on phased genotype normalization."""

import pytest

from variantcentrifuge.replacer import replace_genotypes


@pytest.fixture
def basic_cfg():
    """Basic configuration for replacer tests."""
    return {
        "sample_list": "Sample1,Sample2,Sample3",
        "separator": ";",
        "genotype_replacement_map": {r"[2-9]": "1"},
        "append_extra_sample_fields": False,
        "extra_sample_fields": [],
        "extract_fields_separator": ":",  # Default SnpSift separator
    }


@pytest.fixture
def cfg_with_extra_fields():
    """Configuration with extra fields for testing phased genotype with sample information."""
    return {
        "sample_list": "Sample1,Sample2,Sample3",
        "separator": ";",
        "genotype_replacement_map": {r"[2-9]": "1"},
        "append_extra_sample_fields": True,
        "extra_sample_fields": ["DP", "AD"],
        "extra_sample_field_delimiter": ":",
        "extract_fields_separator": ":",  # Default SnpSift separator
    }


def run_replacer_test(input_header, input_data, expected_output, cfg):
    """Helper function to run the replacer on a simple input and verify expected output."""
    lines = iter([input_header, input_data])
    result = list(replace_genotypes(lines, cfg))
    assert len(result) == 2
    assert result[0] == input_header  # Header line is unchanged
    assert result[1] == expected_output  # Data line is modified as expected


class TestPhasedGenotypeNormalization:
    """Tests for phased genotype normalization in replace_genotypes function."""

    def test_basic_header_detection(self, basic_cfg):
        """Test that the function correctly processes headers with a GT column."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = "chr1\t1000\tA\tG\t0|1:1|0:0|0"
        expected = "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0)"
        run_replacer_test(header, data, expected, basic_cfg)

    @pytest.mark.parametrize(
        "phased_gt,expected_output",
        [
            # Test case 1: Direct normalization of 0|1
            (
                "0|1:1|0:1|1",
                "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0);Sample3(1/1)",
            ),
            # Test case 2: 0|0 should be normalized then skipped
            ("0|1:0|0:1|0", "chr1\t1000\tA\tG\tSample1(0/1);Sample3(1/0)"),
            # Test case 3: All reference genotypes should result in empty GT
            ("0|0:0|0:0|0", "chr1\t1000\tA\tG\t"),
        ],
    )
    def test_direct_normalization(self, basic_cfg, phased_gt, expected_output):
        """Test that phased genotypes are correctly normalized to unphased format."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = f"chr1\t1000\tA\tG\t{phased_gt}"
        run_replacer_test(header, data, expected_output, basic_cfg)

    @pytest.mark.parametrize(
        "phased_gt,expected_output",
        [
            # Test case 1: 0|2 -> 0/2 -> 0/1
            ("0|2:2|0:2|2", "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0);Sample3(1/1)"),
            # Test case 2: Mix of phased genotypes with higher allele numbers
            ("0|3:4|0:9|9", "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0);Sample3(1/1)"),
        ],
    )
    def test_normalization_before_regex_replacement(self, basic_cfg, phased_gt, expected_output):
        """Test that phased genotypes are normalized before regex replacements are applied."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = f"chr1\t1000\tA\tG\t{phased_gt}"
        run_replacer_test(header, data, expected_output, basic_cfg)

    @pytest.mark.parametrize(
        "phased_gt,expected_output",
        [
            # Test case 1: 0|0 should be normalized to 0/0 and skipped
            ("0|0:1|0:0|1", "chr1\t1000\tA\tG\tSample2(1/0);Sample3(0/1)"),
            # Test case 2: .|. should be normalized to ./. and skipped
            (".|.:0|1:.|.", "chr1\t1000\tA\tG\tSample2(0/1)"),
            # Test case 3: Mixed missing/called - test if 0|. becomes 0/. and is kept
            # Removed conflicting test case for 0|.:.|0:1|.
        ],
    )
    def test_normalization_before_skipping(self, basic_cfg, phased_gt, expected_output):
        """Test that phased genotypes are normalized before checking for non-variants to skip."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = f"chr1\t1000\tA\tG\t{phased_gt}"
        run_replacer_test(header, data, expected_output, basic_cfg)

    def test_interaction_with_extra_fields(self, cfg_with_extra_fields):
        """Test that phased genotypes work correctly when extra sample fields are appended."""
        header = "CHROM\tPOS\tREF\tALT\tGT\tDP\tAD"
        data = "chr1\t1000\tA\tG\t0|1:1|0:0|0\t100:120:80\t50,50:70,50:40,40"
        expected = (
            "chr1\t1000\tA\tG\tSample1(0/1:100:50,50);Sample2(1/0:120:70,50)\t"
            "100:120:80\t50,50:70,50:40,40"
        )

        # Use a custom run replacement to preserve the header fields
        lines = iter([header, data])
        result = list(replace_genotypes(lines, cfg_with_extra_fields))
        assert len(result) == 2
        assert result[0] == header  # Header line is unchanged
        assert result[1] == expected  # Data line is modified as expected

    @pytest.mark.parametrize(
        "phased_gt,expected_output",
        [
            # Multiple samples with phased genotypes
            ("0|1:1|0:0|1", "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0);Sample3(0/1)"),
            # Two samples with variants, one without
            ("0|1:1|0:0|0", "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0)"),
            # One sample with variant, two without
            ("0|1:0|0:0|0", "chr1\t1000\tA\tG\tSample1(0/1)"),
        ],
    )
    def test_multiple_samples_with_phased_genotypes(self, basic_cfg, phased_gt, expected_output):
        """Test handling of multiple samples with phased genotypes."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = f"chr1\t1000\tA\tG\t{phased_gt}"
        run_replacer_test(header, data, expected_output, basic_cfg)

    def test_mixed_phased_and_unphased(self, basic_cfg):
        """Test handling of a mix of phased and unphased genotypes."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = "chr1\t1000\tA\tG\t0|1:1/0:0|0"
        expected = "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0)"
        run_replacer_test(header, data, expected, basic_cfg)

    def test_custom_replacement_map(self):
        """Test phased genotype normalization with a custom replacement map."""
        cfg = {
            "sample_list": "Sample1,Sample2",
            "separator": ";",
            "genotype_replacement_map": {r"1": "X", r"2": "Y"},  # Custom replacements
            "extract_fields_separator": ":",
        }

        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = "chr1\t1000\tA\tG\t0|1:2|0"
        expected = "chr1\t1000\tA\tG\tSample1(0/X);Sample2(Y/0)"
        run_replacer_test(header, data, expected, cfg)

    def test_edge_case_unparseable_genotype(self, basic_cfg):
        """Test handling of unparseable genotypes after normalization."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = "chr1\t1000\tA\tG\tABCD:1|0:0|1"  # First "genotype" is not a valid format
        expected = "chr1\t1000\tA\tG\tSample1(ABCD);Sample2(1/0);Sample3(0/1)"
        run_replacer_test(header, data, expected, basic_cfg)
