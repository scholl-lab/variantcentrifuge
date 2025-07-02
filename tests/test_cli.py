# File: tests/test_cli.py
# Location: variantcentrifuge/tests/test_cli.py

"""
Tests for CLI module.

This file contains tests ensuring the CLI runs and shows help correctly.
"""

import subprocess


def test_cli_help():
    """Test that the CLI help message can be displayed."""
    cmd = ["python", "-m", "variantcentrifuge.cli", "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    assert "usage:" in result.stdout


def test_bcftools_prefilter_in_help():
    """Test that the bcftools-prefilter argument is in the help text."""
    cmd = ["python", "-m", "variantcentrifuge.cli", "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    assert "--bcftools-prefilter" in result.stdout
    # Check that the help text mentions bcftools expression
    assert "bcftools expression" in result.stdout.lower()
