# File: tests/test_filters.py
# Location: variantcentrifuge/tests/test_filters.py

"""
Tests for filters module.

This file contains tests for the filtering functions to ensure
they handle variant streams as expected.
"""

from variantcentrifuge.filters import extract_variants


def test_extract_variants_no_bed(monkeypatch):
    """
    Test extract_variants when run_command_stream is mocked.
    """

    def mock_run_command(cmd, input_stream=None):
        yield "##fileformat=VCFv4.2"
        yield "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        yield "1\t1000\t.\tA\tG\t.\t.\t."

    monkeypatch.setattr("variantcentrifuge.filters.run_command_stream", mock_run_command)
    cfg = {}
    lines = list(extract_variants("fake.vcf", "fake.bed", cfg))
    assert len(lines) > 0
    assert lines[0].startswith("##fileformat")
