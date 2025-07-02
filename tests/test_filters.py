# File: tests/test_filters.py
# Location: variantcentrifuge/tests/test_filters.py

"""
Tests for filters module.

This file contains tests for the filtering functions to ensure
they handle variant streams as expected.
"""

import pytest
from variantcentrifuge.filters import extract_variants, apply_bcftools_prefilter


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
    
    cfg = {
        "threads": 4,
        "bcftools_prefilter": 'FILTER="PASS" && INFO/AC<10'
    }
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
    filter_expr = 'INFO/DP > 20'
    
    result = apply_bcftools_prefilter(input_vcf, output_vcf, filter_expr, cfg)
    
    assert result == output_vcf
    
    # Check that default thread count (1) was used
    view_cmd = called_commands[0]
    assert "--threads" in view_cmd
    thread_idx = view_cmd.index("--threads")
    assert view_cmd[thread_idx + 1] == "1"
