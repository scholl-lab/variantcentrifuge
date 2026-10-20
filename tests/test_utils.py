"""
Test module for utility functions in variantcentrifuge/utils.py.

This module contains tests for various utility functions in the
variantcentrifuge.utils module, focusing on header normalization
for SnpEff and SnpSift output.
"""

from variantcentrifuge.utils import normalize_vcf_headers


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


def test_get_tool_version_stderr_parsing(monkeypatch):
    """Test that get_tool_version parses snpEff and SnpSift versions and skips echo/error noise."""
    import subprocess

    from variantcentrifuge.utils import get_tool_version

    calls = []

    def fake_run(cmd, **_kwargs):
        calls.append(cmd)
        if cmd[0] == "snpEff":
            return subprocess.CompletedProcess(
                cmd, 0, stdout="", stderr="SnpEff\t5.2\t2023-09-29\n"
            )
        return subprocess.CompletedProcess(
            cmd,
            255,
            stdout="",
            stderr=(
                "Error: Unknown command '-version'\n"
                "SnpSift -version -version\n"
                "SnpSift version 5.2 (build 2023-09-29), by Pablo Cingolani\n"
            ),
        )

    monkeypatch.setattr("variantcentrifuge.utils.subprocess.run", fake_run)
    assert get_tool_version("snpEff") == "SnpEff\t5.2\t2023-09-29"
    assert (
        get_tool_version("SnpSift") == "SnpSift version 5.2 (build 2023-09-29), by Pablo Cingolani"
    )
    assert ["snpEff", "-version"] in calls
    assert ["SnpSift", "annotate", "-version"] in calls


def test_get_tool_version_parses_bare_snpsift_version(monkeypatch):
    """SnpSift 5.2 `annotate -version` prints only the number (seen on 1000G e2e)."""
    import subprocess

    from variantcentrifuge.utils import get_tool_version

    monkeypatch.setattr(
        "variantcentrifuge.utils.subprocess.run",
        lambda cmd, **_kw: subprocess.CompletedProcess(cmd, 0, stdout="5.2\n", stderr=""),
    )
    assert get_tool_version("SnpSift") == "SnpSift 5.2"


def test_get_tool_version_errors_and_git_checkout(monkeypatch, tmp_path):
    """Test that invalid tool/git outputs return N/A and git queries VariantCentrifuge's repo."""
    import subprocess
    from pathlib import Path

    import variantcentrifuge.utils as utils_mod
    from variantcentrifuge.utils import get_tool_version

    expected_repo_dir = str(Path(utils_mod.__file__).resolve().parent)
    monkeypatch.chdir(tmp_path)

    # Real git invocation from a non-repo working directory still resolves VariantCentrifuge SHA
    real_sha = get_tool_version("git")
    assert len(real_sha) >= 7 and all(c in "0123456789abcdefABCDEF" for c in real_sha)

    def fake_error_run(cmd, **_kwargs):
        if cmd[0] == "git":
            assert cmd[:3] == ["git", "-C", expected_repo_dir]
            return subprocess.CompletedProcess(
                cmd,
                128,
                stdout="",
                stderr="fatal: not a git repository (or any of the parent directories): .git\n",
            )
        return subprocess.CompletedProcess(
            cmd, 1, stdout="", stderr="Error: Unable to access jarfile SnpSift.jar\n"
        )

    monkeypatch.setattr("variantcentrifuge.utils.subprocess.run", fake_error_run)
    assert get_tool_version("git") == "N/A"
    assert get_tool_version("SnpSift") == "N/A"
    assert get_tool_version("snpEff") == "N/A"
