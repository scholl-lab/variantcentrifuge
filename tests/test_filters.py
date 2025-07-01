# File: tests/test_filters.py
# Location: variantcentrifuge/tests/test_filters.py

"""
Tests for filters module.

This file contains tests for the filtering functions to ensure
they handle variant streams as expected.
"""

from variantcentrifuge.filters import extract_variants


def test_extract_variants_no_bed(monkeypatch):
    """Test extract_variants when run_command is mocked."""

    def mock_run_command(cmd):
        # Instead of yielding content, just simulate success
        return 0

    monkeypatch.setattr("variantcentrifuge.filters.run_command", mock_run_command)
    cfg = {}
    output_file = "test_output.vcf.gz"
    result = extract_variants("fake.vcf", "fake.bed", cfg, output_file)
    # Check that the function returns the expected output file path
    assert result == output_file
