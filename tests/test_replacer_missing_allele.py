"""Unit tests for the replacer.py module, focusing on missing allele normalization."""

import pytest

from variantcentrifuge.replacer import replace_genotypes


@pytest.fixture
def basic_cfg():
    """Provide basic configuration for replacer tests."""
    return {
        "sample_list": "Sample1,Sample2,Sample3",
        "separator": ";",
        "genotype_replacement_map": {r"[2-9]": "1"},
        "append_extra_sample_fields": False,
        "extra_sample_fields": [],
        # Default SnpSift separator
        "extract_fields_separator": ":",
    }


@pytest.fixture
def cfg_with_extra_fields():
    """Provide configuration with extra fields for testing genotype."""
    return {
        "sample_list": "Sample1,Sample2,Sample3",
        "separator": ";",
        "genotype_replacement_map": {r"[2-9]": "1"},
        "append_extra_sample_fields": True,
        "extra_sample_fields": ["DP", "AD"],
        "extra_sample_field_delimiter": ":",
        # Default SnpSift separator
        "extract_fields_separator": ":",
    }


def run_replacer_test(input_header, input_data, expected_output, cfg):
    """Run the replacer on a simple input and verify expected output."""
    lines = iter([input_header, input_data])
    result = list(replace_genotypes(lines, cfg))
    assert len(result) == 2
    assert result[0] == input_header  # Header line is unchanged
    assert result[1] == expected_output  # Data line is modified as expected


class TestMissingAlleleNormalization:
    """Tests for missing allele normalization in replace_genotypes function."""

    def test_basic_header_detection(self, basic_cfg):
        """Test that the function correctly processes headers with a GT column."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = "chr1\t1000\tA\tG\t./1:1/.:0/0"
        expected = "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0)"
        run_replacer_test(header, data, expected, basic_cfg)

    @pytest.mark.parametrize(
        "missing_gt,expected_output",
        [
            # Test case 1: Basic ./1 normalization
            (
                "./1:1/.:1/1",
                "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0);Sample3(1/1)",
            ),
            # Test case 2: 0/0 should be skipped
            ("./1:0/0:1/.", "chr1\t1000\tA\tG\tSample1(0/1);Sample3(1/0)"),
            # Test case 3: ./. should be skipped
            ("./1:./.:1/.", "chr1\t1000\tA\tG\tSample1(0/1);Sample3(1/0)"),
        ],
    )
    def test_basic_missing_allele_normalization(self, basic_cfg, missing_gt, expected_output):
        """Test that genotypes with one missing allele are correctly normalized."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = f"chr1\t1000\tA\tG\t{missing_gt}"
        run_replacer_test(header, data, expected_output, basic_cfg)

    @pytest.mark.parametrize(
        "missing_gt,expected_output",
        [
            # Test case 1: ./2 -> 0/2 -> 0/1
            ("./2:2/.:3/4", "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0);Sample3(1/1)"),
            # Test case 2: Higher allele numbers
            ("./9:7/.:8/9", "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0);Sample3(1/1)"),
        ],
    )
    def test_normalization_with_higher_alleles(self, basic_cfg, missing_gt, expected_output):
        """Test that missing alleles are normalized before regex replacements are applied."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = f"chr1\t1000\tA\tG\t{missing_gt}"
        run_replacer_test(header, data, expected_output, basic_cfg)

    @pytest.mark.parametrize(
        "phased_missing_gt,expected_output",
        [
            # Test case 1: .|1 -> ./1 -> 0/1
            (".|1:1|.:.|.", "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0)"),
            # Test case 2: Mixed phased notation
            (".|2:3|.:2|3", "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0);Sample3(1/1)"),
        ],
    )
    def test_interaction_with_phasing(self, basic_cfg, phased_missing_gt, expected_output):
        """Test interaction between phased notation and missing allele normalization."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = f"chr1\t1000\tA\tG\t{phased_missing_gt}"
        run_replacer_test(header, data, expected_output, basic_cfg)

    @pytest.mark.parametrize(
        "standard_gt,expected_output",
        [
            # Test case 1: Standard genotypes shouldn't be affected
            ("0/1:1/1:0/0", "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/1)"),
            # Test case 2: Mixed standard and non-standard
            ("0/1:1/.:./1", "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0);Sample3(0/1)"),
        ],
    )
    def test_unaffected_genotypes(self, basic_cfg, standard_gt, expected_output):
        """Test that standard genotypes are unaffected by missing allele normalization."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = f"chr1\t1000\tA\tG\t{standard_gt}"
        run_replacer_test(header, data, expected_output, basic_cfg)

    def test_interaction_with_extra_fields(self, cfg_with_extra_fields):
        """Test that missing allele normalization works with extra sample fields."""
        header = "CHROM\tPOS\tREF\tALT\tGT\tDP\tAD"
        data = "chr1\t1000\tA\tG\t./1:1/.:0/0\t100:120:80\t50,50:70,50:40,40"
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

    def test_multiple_samples_with_mixed_genotypes(self, basic_cfg):
        """Test handling of multiple samples with mixed genotypes including missing alleles."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = "chr1\t1000\tA\tG\t./1:1/.:0/0"
        expected = "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0)"
        run_replacer_test(header, data, expected, basic_cfg)

    def test_edge_cases(self, basic_cfg):
        """Test edge cases for missing allele normalization."""
        header = "CHROM\tPOS\tREF\tALT\tGT"
        # Test with multi-digit allele numbers
        data = "chr1\t1000\tA\tG\t./12:12/.:./."
        # After normalization and replacement
        expected = "chr1\t1000\tA\tG\tSample1(0/1);Sample2(1/0)"
        run_replacer_test(header, data, expected, basic_cfg)

    def test_custom_replacement_with_missing(self):
        """Test missing allele normalization with custom replacement map."""
        cfg = {
            "sample_list": "Sample1,Sample2",
            "separator": ";",
            "genotype_replacement_map": {r"1": "X", r"2": "Y"},  # Custom replacements
            "extract_fields_separator": ":",
        }

        header = "CHROM\tPOS\tREF\tALT\tGT"
        data = "chr1\t1000\tA\tG\t./1:2/."
        expected = "chr1\t1000\tA\tG\tSample1(0/X);Sample2(Y/0)"
        run_replacer_test(header, data, expected, cfg)
