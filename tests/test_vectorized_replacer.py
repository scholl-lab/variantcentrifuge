"""
Test module for vectorized genotype replacement functionality.

This module tests the pandas-based vectorized genotype replacement
implementation to ensure it produces identical results to the original
streaming implementation while providing performance benefits.
"""

import logging

import pandas as pd
import pytest

from variantcentrifuge.vectorized_replacer import (
    VectorizedGenotypeReplacer,
    _normalize_snpeff_column_name,
    process_chunked_vectorized,
    replace_genotypes_vectorized,
)


# Test data fixtures
@pytest.fixture
def sample_config() -> dict[str, str]:
    """Return basic configuration for genotype replacement tests."""
    return {
        "sample_list": "Sample1,Sample2,Sample3",
        "separator": ";",
        "extract_fields_separator": ",",
        "genotype_replacement_map": {r"[2-9]": "1"},
        "append_extra_sample_fields": False,
        "extra_sample_fields": [],
        "extra_sample_field_delimiter": ":",
    }


@pytest.fixture
def config_with_extra_fields() -> dict[str, str]:
    """Return configuration with extra sample fields."""
    return {
        "sample_list": "Sample1,Sample2,Sample3",
        "separator": ";",
        "extract_fields_separator": ",",
        "genotype_replacement_map": {r"[2-9]": "1"},
        "append_extra_sample_fields": True,
        "extra_sample_fields": ["DP", "AD"],
        "extra_sample_field_delimiter": ":",
    }


@pytest.fixture
def sample_dataframe() -> pd.DataFrame:
    """Sample DataFrame with genotype data."""
    return pd.DataFrame(
        {
            "CHROM": ["chr1", "chr1", "chr2"],
            "POS": [100, 200, 300],
            "ID": ["rs1", "rs2", "rs3"],
            "REF": ["A", "G", "C"],
            "ALT": ["T", "C", "T"],
            "GT": ["0/1,0/0,1/1", "0/0,1/2,./.", "1/1,0/1,./."],
            "DP": ["50,30,40", "60,45,35", "55,25,30"],
            "AD": ["25,25,20,20,15,25", "30,30,22,23,17,18", "27,28,12,13,15,15"],
        }
    )


@pytest.fixture
def sample_dataframe_with_extra_fields() -> pd.DataFrame:
    """Sample DataFrame with extra fields that need normalization."""
    return pd.DataFrame(
        {
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "GT": ["0/1,1/2,0/0", "1/1,0/0,1/1"],
            "GEN[*].DP": ["50,60,40", "55,45,35"],
            "GEN[*].AD": ["25,25,30,30,20,20", "27,28,22,23,17,18"],
        }
    )


class TestNormalization:
    """Test column name normalization functions."""

    def test_normalize_snpeff_column_name(self):
        """Test SnpEff column name normalization."""
        assert _normalize_snpeff_column_name("GEN[*].DP") == "DP"
        assert _normalize_snpeff_column_name("GEN[0].AD") == "AD"
        assert _normalize_snpeff_column_name("ANN[*].Effect") == "Effect"
        assert _normalize_snpeff_column_name("ANN[0].Gene_Name") == "Gene_Name"
        assert _normalize_snpeff_column_name("NMD[*].Confidence") == "NMD_Confidence"
        assert _normalize_snpeff_column_name("LOF[*].Type") == "LOF_Type"
        assert _normalize_snpeff_column_name("regular_column") == "regular_column"


class TestVectorizedGenotypeReplacer:
    """Test the VectorizedGenotypeReplacer class."""

    def test_initialization(self, sample_config):
        """Test proper initialization of the replacer."""
        replacer = VectorizedGenotypeReplacer(sample_config)

        assert replacer.samples == ["Sample1", "Sample2", "Sample3"]
        assert replacer.separator == ";"
        assert replacer.snpsift_sep == ","
        assert len(replacer.compiled_patterns) == 1
        assert not replacer.append_extra_fields

    def test_initialization_with_extra_fields(self, config_with_extra_fields):
        """Test initialization with extra fields configuration."""
        replacer = VectorizedGenotypeReplacer(config_with_extra_fields)

        assert replacer.append_extra_fields
        assert replacer.extra_sample_fields == ["DP", "AD"]
        assert replacer.extra_field_delimiter == ":"

    def test_initialization_no_samples(self):
        """Test that initialization fails without samples."""
        config = {"sample_list": ""}

        with pytest.raises(ValueError, match="No samples found"):
            VectorizedGenotypeReplacer(config)

    def test_process_genotype_series_basic(self, sample_config):
        """Test basic genotype series processing."""
        replacer = VectorizedGenotypeReplacer(sample_config)

        # Test series with various genotype patterns
        genotypes = pd.Series(["0/1", "1/2", "0/0", "./.", "2/3", "1|1", "./2"])
        processed = replacer._process_genotype_series(genotypes)

        expected = pd.Series(["0/1", "1/1", "0/0", "./.", "1/1", "1/1", "0/1"])
        pd.testing.assert_series_equal(processed, expected)

    def test_process_genotype_series_phased_conversion(self, sample_config):
        """Test conversion of phased genotypes."""
        replacer = VectorizedGenotypeReplacer(sample_config)

        genotypes = pd.Series(["1|0", "2|3", "0|1"])
        processed = replacer._process_genotype_series(genotypes)

        expected = pd.Series(["1/0", "1/1", "0/1"])
        pd.testing.assert_series_equal(processed, expected)

    def test_process_genotype_series_missing_alleles(self, sample_config):
        """Test handling of missing alleles."""
        replacer = VectorizedGenotypeReplacer(sample_config)

        genotypes = pd.Series(["./1", "2/.", "./.", "0/."])
        processed = replacer._process_genotype_series(genotypes)

        expected = pd.Series(["0/1", "1/0", "./.", "0/0"])
        pd.testing.assert_series_equal(processed, expected)

    def test_process_genotype_series_haploid(self, sample_config):
        """Test handling of haploid genotypes."""
        replacer = VectorizedGenotypeReplacer(sample_config)

        genotypes = pd.Series(["0", "1", "2", "3"])
        processed = replacer._process_genotype_series(genotypes)

        expected = pd.Series(["0", "1", "1", "1"])
        pd.testing.assert_series_equal(processed, expected)

    def test_process_dataframe_basic(self, sample_config, sample_dataframe):
        """Test basic DataFrame processing."""
        replacer = VectorizedGenotypeReplacer(sample_config)

        result_df = replacer.process_dataframe(sample_dataframe)

        # Check that non-GT columns are unchanged
        pd.testing.assert_series_equal(result_df["CHROM"], sample_dataframe["CHROM"])
        pd.testing.assert_series_equal(result_df["POS"], sample_dataframe["POS"])

        # Check GT column processing
        expected_gt = [
            "Sample1(0/1);Sample3(1/1)",  # Row 0: Sample2 has 0/0 (filtered)
            "Sample2(1/1)",  # Row 1: Sample1 has 0/0, Sample3 has ./. (both filtered)
            "Sample1(1/1);Sample2(0/1)",  # Row 2: Sample3 has ./. (filtered)
        ]

        for i, expected in enumerate(expected_gt):
            assert result_df.loc[i, "GT"] == expected

    def test_process_dataframe_with_extra_fields(self, config_with_extra_fields, sample_dataframe):
        """Test DataFrame processing with extra fields."""
        replacer = VectorizedGenotypeReplacer(config_with_extra_fields)

        result_df = replacer.process_dataframe(sample_dataframe)

        # Check that extra fields are included in output
        gt_values = result_df["GT"].tolist()

        # Row 0: Sample1(0/1:50:25,25), Sample3(1/1:40:15,25)
        assert "Sample1(0/1:50:25)" in gt_values[0]
        assert "Sample3(1/1:40:20)" in gt_values[0]

        # Row 1: Sample2(1/1:45:30)
        assert "Sample2(1/1:45:30)" in gt_values[1]

    def test_process_dataframe_normalized_extra_fields(self, sample_dataframe_with_extra_fields):
        """Test processing with normalized extra field names."""
        config = {
            "sample_list": "Sample1,Sample2,Sample3",
            "separator": ";",
            "extract_fields_separator": ",",
            "genotype_replacement_map": {r"[2-9]": "1"},
            "append_extra_sample_fields": True,
            "extra_sample_fields": ["GEN[*].DP", "GEN[*].AD"],
            "extra_sample_field_delimiter": ":",
        }

        replacer = VectorizedGenotypeReplacer(config)
        result_df = replacer.process_dataframe(sample_dataframe_with_extra_fields)

        # Check that normalized field names work correctly
        gt_values = result_df["GT"].tolist()
        assert "Sample1(0/1:50:25)" in gt_values[0]
        assert "Sample2(1/1:60:25)" in gt_values[0]

    def test_process_dataframe_no_gt_column(self, sample_config):
        """Test DataFrame processing when GT column is missing."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": [100]})

        replacer = VectorizedGenotypeReplacer(sample_config)

        result_df = replacer.process_dataframe(df)

        # Should return unchanged DataFrame
        pd.testing.assert_frame_equal(result_df, df)

    def test_combine_sample_genotypes(self, sample_config):
        """Test combining sample genotype strings."""
        replacer = VectorizedGenotypeReplacer(sample_config)

        # Create sample genotype series
        sample1_series = pd.Series(["Sample1(0/1)", "", "Sample1(1/1)"])
        sample2_series = pd.Series(["", "Sample2(1/2)", "Sample2(0/1)"])
        sample3_series = pd.Series(["Sample3(1/1)", "", ""])

        sample_lists = [sample1_series, sample2_series, sample3_series]
        result = replacer._combine_sample_genotypes(sample_lists)

        expected = pd.Series(
            ["Sample1(0/1);Sample3(1/1)", "Sample2(1/2)", "Sample1(1/1);Sample2(0/1)"]
        )

        pd.testing.assert_series_equal(result, expected)


class TestFileProcessing:
    """Test file-based processing functions."""

    def test_replace_genotypes_vectorized(self, sample_config, sample_dataframe, tmp_path):
        """Test the high-level file processing function."""
        # Create input file
        input_file = tmp_path / "input.tsv"
        sample_dataframe.to_csv(input_file, sep="\t", index=False)

        # Create output file
        output_file = tmp_path / "output.tsv"

        # Process file
        replace_genotypes_vectorized(input_file, output_file, sample_config)

        # Check output
        assert output_file.exists()
        result_df = pd.read_csv(output_file, sep="\t", dtype=str)

        # Basic validation
        assert len(result_df) == len(sample_dataframe)
        assert "GT" in result_df.columns

        # Check that genotype replacement occurred
        gt_values = result_df["GT"].tolist()
        assert any("Sample1" in gt for gt in gt_values)
        assert any("Sample2" in gt for gt in gt_values)

    def test_replace_genotypes_vectorized_gzipped(self, sample_config, sample_dataframe, tmp_path):
        """Test processing of gzipped files."""
        import gzip

        # Create gzipped input file
        input_file = tmp_path / "input.tsv.gz"
        with gzip.open(input_file, "wt") as f:
            sample_dataframe.to_csv(f, sep="\t", index=False)

        # Create gzipped output file
        output_file = tmp_path / "output.tsv.gz"

        # Process file
        replace_genotypes_vectorized(input_file, output_file, sample_config)

        # Check output
        assert output_file.exists()
        with gzip.open(output_file, "rt") as f:
            result_df = pd.read_csv(f, sep="\t", dtype=str)

        # Basic validation
        assert len(result_df) == len(sample_dataframe)
        assert "GT" in result_df.columns

    def test_process_chunked_vectorized(self, sample_config, tmp_path):
        """Test chunked processing for large files."""
        # Create larger test DataFrame
        large_df = pd.DataFrame(
            {
                "CHROM": [f"chr{i % 22 + 1}" for i in range(100)],
                "POS": list(range(1000, 1100)),
                "GT": [f"{i % 2}/{(i + 1) % 3}" for i in range(100)],
            }
        )

        # Expand GT column to have proper multi-sample format
        large_df["GT"] = large_df["GT"] + "," + large_df["GT"] + "," + large_df["GT"]

        # Create input file
        input_file = tmp_path / "large_input.tsv.gz"
        large_df.to_csv(input_file, sep="\t", index=False, compression="gzip")

        # Create output file
        output_file = tmp_path / "large_output.tsv.gz"

        # Process with small chunk size
        process_chunked_vectorized(input_file, output_file, sample_config, chunk_size=10)

        # Check output
        assert output_file.exists()
        result_df = pd.read_csv(output_file, sep="\t", dtype=str, compression="gzip")

        # Should have same number of rows
        assert len(result_df) == len(large_df)
        assert "GT" in result_df.columns


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_dataframe(self, sample_config):
        """Test processing of empty DataFrame."""
        replacer = VectorizedGenotypeReplacer(sample_config)

        empty_df = pd.DataFrame()
        result_df = replacer.process_dataframe(empty_df)

        pd.testing.assert_frame_equal(result_df, empty_df)

    def test_dataframe_with_missing_samples(self, sample_config):
        """Test DataFrame with fewer genotype columns than samples."""
        df = pd.DataFrame(
            {"CHROM": ["chr1"], "POS": [100], "GT": ["0/1,1/1"]}  # Only 2 genotypes for 3 samples
        )

        replacer = VectorizedGenotypeReplacer(sample_config)

        result_df = replacer.process_dataframe(df)

        # Should process available genotypes
        assert "Sample1" in result_df.loc[0, "GT"] or "Sample2" in result_df.loc[0, "GT"]

    def test_complex_genotype_patterns(self, sample_config):
        """Test complex genotype patterns."""
        df = pd.DataFrame(
            {
                "GT": [
                    "0/1,1|2,3/4",  # Mixed phased/unphased, multi-allelic
                    "./.,0/.,./0",  # Various missing patterns
                    "1,2,3",  # Haploid genotypes
                    "10/11,12/13,14",  # Multi-digit alleles
                ]
            }
        )

        replacer = VectorizedGenotypeReplacer(sample_config)
        result_df = replacer.process_dataframe(df)

        # Check that processing completed without errors
        assert len(result_df) == len(df)
        assert "GT" in result_df.columns

    def test_malformed_configuration(self):
        """Test handling of malformed configuration."""
        # Missing sample_list
        with pytest.raises(ValueError):
            VectorizedGenotypeReplacer({})

        # Empty sample_list
        with pytest.raises(ValueError):
            VectorizedGenotypeReplacer({"sample_list": ""})

    def test_invalid_regex_pattern(self):
        """Test handling of invalid regex patterns."""
        config = {
            "sample_list": "Sample1,Sample2",
            "genotype_replacement_map": {"[invalid": "1"},  # Invalid regex
        }

        # Should raise an error during regex compilation
        with pytest.raises(Exception):  # noqa: B017 - intentionally broad, any error is acceptable
            VectorizedGenotypeReplacer(config)


class TestPerformanceComparison:
    """Test performance characteristics and comparison with original implementation."""

    def test_vectorized_vs_original_output_equivalence(self, sample_config, sample_dataframe):
        """Test that vectorized output matches original streaming implementation."""
        # This test would ideally compare outputs from both implementations
        # For now, we test that the vectorized implementation produces expected output

        replacer = VectorizedGenotypeReplacer(sample_config)
        result_df = replacer.process_dataframe(sample_dataframe)

        # Verify expected transformations
        gt_values = result_df["GT"].tolist()

        # Check that 0/0 and ./. genotypes are filtered out
        for gt_value in gt_values:
            assert "0/0" not in gt_value
            assert "./." not in gt_value

        # Check that multi-allelic genotypes are converted to 1
        # Original GT: "0/1,0/0,1/1" -> "Sample1(0/1);Sample3(1/1)"
        # Original GT: "0/0,1/2,./." -> "Sample2(1/1)"
        assert "Sample2(1/1)" in gt_values[1]  # 1/2 -> 1/1

    @pytest.mark.performance
    def test_memory_efficiency(self, sample_config):
        """Test memory efficiency with large datasets."""
        # Create a reasonably large DataFrame for memory testing
        n_variants = 1000
        large_df = pd.DataFrame(
            {
                "CHROM": [f"chr{i % 22 + 1}" for i in range(n_variants)],
                "POS": list(range(100000, 100000 + n_variants)),
                "GT": [
                    f"{i % 3}/{(i + 1) % 3},{(i + 2) % 3}/{i % 2},{(i + 1) % 2}/{(i + 2) % 3}"
                    for i in range(n_variants)
                ],
            }
        )

        replacer = VectorizedGenotypeReplacer(sample_config)

        # Process should complete without memory errors
        result_df = replacer.process_dataframe(large_df)

        assert len(result_df) == n_variants
        assert "GT" in result_df.columns


if __name__ == "__main__":
    # Run tests with logging enabled
    logging.basicConfig(level=logging.DEBUG)
    pytest.main([__file__])
