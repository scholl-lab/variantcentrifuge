# File: tests/test_filters.py
# Location: variantcentrifuge/tests/test_filters.py

"""
Tests for filters module.

This file contains tests for the filtering functions to ensure
they handle variant streams as expected.
"""

import pytest
import pandas as pd
from variantcentrifuge.filters import (
    extract_variants,
    filter_dataframe_with_query,
)


def test_extract_variants_no_bed(monkeypatch):
    """Test extract_variants when run_command is mocked."""

    called_commands = []

    def mock_run_command(cmd):
        # Record the command
        called_commands.append(cmd)
        return 0

    monkeypatch.setattr("variantcentrifuge.filters.run_command", mock_run_command)
    cfg = {}
    output_file = "test_output.vcf.gz"
    result = extract_variants("fake.vcf", "fake.bed", cfg, output_file)
    # Check that the function returns the expected output file path
    assert result == output_file
    # Check that bcftools was called with the correct basic arguments
    assert len(called_commands) == 1
    cmd = called_commands[0]
    assert cmd[0] == "bcftools"
    assert cmd[1] == "view"
    assert "-R" in cmd
    assert "fake.bed" in cmd
    assert "-o" in cmd
    assert output_file in cmd
    assert "fake.vcf" in cmd
    # Should NOT have -i flag when no prefilter is specified
    assert "-i" not in cmd


def test_extract_variants_with_prefilter(monkeypatch):
    """Test extract_variants with bcftools prefilter."""

    called_commands = []

    def mock_run_command(cmd):
        called_commands.append(cmd)
        return 0

    monkeypatch.setattr("variantcentrifuge.filters.run_command", mock_run_command)

    cfg = {"threads": 4, "bcftools_prefilter": 'FILTER="PASS" && INFO/AC<10'}
    output_file = "test_output.vcf.gz"
    result = extract_variants("fake.vcf", "fake.bed", cfg, output_file)

    assert result == output_file
    assert len(called_commands) == 1

    cmd = called_commands[0]
    assert cmd[0] == "bcftools"
    assert cmd[1] == "view"
    assert "--threads" in cmd
    assert "4" in cmd
    assert "-R" in cmd
    assert "fake.bed" in cmd
    assert "-i" in cmd
    # Find the position of -i and check the next element is the filter
    i_idx = cmd.index("-i")
    assert cmd[i_idx + 1] == 'FILTER="PASS" && INFO/AC<10'
    assert "-o" in cmd
    assert output_file in cmd
    assert "fake.vcf" in cmd


def test_apply_bcftools_prefilter(monkeypatch):
    """Test apply_bcftools_prefilter with mocked run_command."""

    called_commands = []

    def mock_run_command(cmd, output_file=None):
        # Record the command for verification
        called_commands.append(cmd)
        return 0

    monkeypatch.setattr("variantcentrifuge.filters.run_command", mock_run_command)

    cfg = {"threads": 4}
    input_vcf = "input.vcf.gz"
    output_vcf = "output.vcf.gz"
    filter_expr = 'FILTER=="PASS" & INFO/AC < 10'

    result = apply_bcftools_prefilter(input_vcf, output_vcf, filter_expr, cfg)

    # Check return value
    assert result == output_vcf

    # Check that bcftools view was called with correct arguments
    assert len(called_commands) == 2  # view and index commands

    # Check the view command
    view_cmd = called_commands[0]
    assert view_cmd[0] == "bcftools"
    assert view_cmd[1] == "view"
    assert "--threads" in view_cmd
    assert "4" in view_cmd
    assert "-i" in view_cmd
    assert filter_expr in view_cmd
    assert "-Oz" in view_cmd
    assert "-o" in view_cmd
    assert output_vcf in view_cmd
    assert input_vcf in view_cmd

    # Check the index command
    index_cmd = called_commands[1]
    assert index_cmd[0] == "bcftools"
    assert index_cmd[1] == "index"
    assert "--threads" in index_cmd
    assert "4" in index_cmd
    assert output_vcf in index_cmd


def test_apply_bcftools_prefilter_default_threads(monkeypatch):
    """Test apply_bcftools_prefilter with default thread count."""

    called_commands = []

    def mock_run_command(cmd, output_file=None):
        called_commands.append(cmd)
        return 0

    monkeypatch.setattr("variantcentrifuge.filters.run_command", mock_run_command)

    # No threads specified in config
    cfg = {}
    input_vcf = "input.vcf.gz"
    output_vcf = "output.vcf.gz"
    filter_expr = "INFO/DP > 20"

    result = apply_bcftools_prefilter(input_vcf, output_vcf, filter_expr, cfg)

    assert result == output_vcf

    # Check that default thread count (1) was used
    view_cmd = called_commands[0]
    assert "--threads" in view_cmd
    thread_idx = view_cmd.index("--threads")
    assert view_cmd[thread_idx + 1] == "1"


class TestFilterDataframeWithQuery:
    """Tests for filter_dataframe_with_query function."""

    def test_empty_filter_expression(self):
        """Test that empty filter expression returns unmodified DataFrame."""
        df = pd.DataFrame(
            {"CHROM": ["chr1", "chr2", "chr3"], "POS": [100, 200, 300], "score": [0.5, 0.8, 0.3]}
        )

        result = filter_dataframe_with_query(df, "")
        pd.testing.assert_frame_equal(result, df)

    def test_none_filter_expression(self):
        """Test that None filter expression returns unmodified DataFrame."""
        df = pd.DataFrame(
            {"CHROM": ["chr1", "chr2", "chr3"], "POS": [100, 200, 300], "score": [0.5, 0.8, 0.3]}
        )

        result = filter_dataframe_with_query(df, None)
        pd.testing.assert_frame_equal(result, df)

    def test_simple_numeric_filter(self):
        """Test simple numeric comparison filter."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr3"],
                "POS": [100, 200, 300],
                "score": ["0.5", "0.8", "0.3"],  # String values that should be converted
            }
        )

        result = filter_dataframe_with_query(df, "score > 0.4")
        assert len(result) == 2
        assert result["CHROM"].tolist() == ["chr1", "chr2"]

    def test_string_equality_filter(self):
        """Test string equality filter."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr3"],
                "IMPACT": ["HIGH", "MODERATE", "HIGH"],
                "score": [0.5, 0.8, 0.3],
            }
        )

        result = filter_dataframe_with_query(df, 'IMPACT == "HIGH"')
        assert len(result) == 2
        assert result["CHROM"].tolist() == ["chr1", "chr3"]

    def test_complex_filter_with_and(self):
        """Test complex filter with AND condition."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr3"],
                "IMPACT": ["HIGH", "HIGH", "MODERATE"],
                "score": ["0.5", "0.8", "0.3"],
            }
        )

        result = filter_dataframe_with_query(df, 'IMPACT == "HIGH" and score > 0.6')
        assert len(result) == 1
        assert result["CHROM"].tolist() == ["chr2"]

    def test_complex_filter_with_or(self):
        """Test complex filter with OR condition."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr3"],
                "IMPACT": ["HIGH", "MODERATE", "LOW"],
                "score": ["0.9", "0.4", "0.2"],
            }
        )

        result = filter_dataframe_with_query(df, 'IMPACT == "HIGH" or score > 0.8')
        assert len(result) == 1
        assert result["CHROM"].tolist() == ["chr1"]

    def test_filter_with_in_operator(self):
        """Test filter using 'in' operator."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr3"],
                "Inheritance_Pattern": ["de_novo", "autosomal_recessive", "compound_heterozygous"],
            }
        )

        result = filter_dataframe_with_query(
            df, 'Inheritance_Pattern in ["de_novo", "compound_heterozygous"]'
        )
        assert len(result) == 2
        assert set(result["Inheritance_Pattern"]) == {"de_novo", "compound_heterozygous"}

    def test_filter_with_string_contains(self):
        """Test filter using string contains method."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr3"],
                "Custom_Annotation": [
                    "InGeneList=cancer_panel",
                    "Region=exon",
                    "InGeneList=cardiac_panel",
                ],
            }
        )

        result = filter_dataframe_with_query(df, 'Custom_Annotation.str.contains("cancer")')
        assert len(result) == 1
        assert result["CHROM"].tolist() == ["chr1"]

    def test_numeric_conversion_with_mixed_types(self):
        """Test that numeric conversion handles mixed types correctly."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr3"],
                "dbNSFP_CADD_phred": ["25.5", "NA", "30.2"],
                "gene": ["BRCA1", "TP53", "EGFR"],
            }
        )

        result = filter_dataframe_with_query(df, "dbNSFP_CADD_phred > 26")
        assert len(result) == 1
        assert result["gene"].tolist() == ["EGFR"]

    def test_invalid_filter_expression(self, caplog):
        """Test that invalid filter expression returns unfiltered data with error log."""
        df = pd.DataFrame({"CHROM": ["chr1", "chr2", "chr3"], "POS": [100, 200, 300]})

        # Invalid syntax
        result = filter_dataframe_with_query(df, 'CHROM === "chr1"')  # Triple equals invalid

        # Should return original dataframe
        pd.testing.assert_frame_equal(result, df)

        # Check error was logged
        assert "Failed to apply final filter expression" in caplog.text

    def test_filter_on_computed_columns(self):
        """Test filtering on columns that might be computed during analysis."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr3"],
                "gene": ["BRCA1", "TP53", "EGFR"],
                "inheritance_score": ["0.95", "0.45", "0.78"],
                "Inheritance_Confidence": ["0.9", "0.3", "0.8"],
            }
        )

        result = filter_dataframe_with_query(
            df, "inheritance_score > 0.7 and Inheritance_Confidence > 0.75"
        )
        assert len(result) == 2
        assert set(result["gene"]) == {"BRCA1", "EGFR"}
