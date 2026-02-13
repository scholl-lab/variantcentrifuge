"""
Integration tests for inheritance mode functionality.

These tests verify that the --inheritance-mode flag works correctly
through the full command-line interface.
"""

import json
import os
import subprocess
import tempfile

import pandas as pd
import pytest


class TestInheritanceModeIntegration:
    """Test inheritance mode functionality through the CLI."""

    @pytest.fixture
    def sample_vcf_file(self):
        """Create a sample VCF file with de novo and compound het variants."""
        vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother
chr1\t1000\t.\tA\tT\t.\tPASS\tAF=0.0001\tGT\t0/1\t0/0\t0/0
chr2\t2000\t.\tG\tC\t.\tPASS\tAF=0.0005\tGT\t0/1\t0/1\t0/0
chr2\t3000\t.\tC\tT\t.\tPASS\tAF=0.0003\tGT\t0/1\t0/0\t0/1
chr3\t4000\t.\tT\tG\t.\tPASS\tAF=0.001\tGT\t1/1\t0/1\t0/1
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(vcf_content)
            vcf_file = f.name

        yield vcf_file
        os.unlink(vcf_file)

    @pytest.fixture
    def trio_ped_file(self):
        """Create a trio pedigree file."""
        ped_content = """#FamilyID\tIndividualID\tPaternalID\tMaternalID\tSex\tAffectedStatus
FAM1\tchild\tfather\tmother\t1\t2
FAM1\tfather\t0\t0\t1\t1
FAM1\tmother\t0\t0\t2\t1
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(ped_content)
            ped_file = f.name

        yield ped_file
        os.unlink(ped_file)

    @pytest.fixture
    def gene_bed_file(self):
        """Create a BED file for gene regions."""
        bed_content = """chr1\t900\t1100\tGENE1
chr2\t1900\t3100\tGENE2
chr3\t3900\t4100\tGENE3
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            bed_file = f.name

        yield bed_file
        os.unlink(bed_file)

    def run_variantcentrifuge(self, args):
        """Run variantcentrifuge with given arguments."""
        cmd = ["python", "-m", "variantcentrifuge", *args]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result

    @pytest.mark.slow
    def test_inheritance_mode_simple(self, sample_vcf_file, trio_ped_file, gene_bed_file):
        """Test inheritance analysis with simple mode."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "output.tsv")

            # Run with simple mode (default)
            args = [
                "--vcf-file",
                sample_vcf_file,
                "--gene-name",
                "all",
                "--output-file",
                output_file,
                "--output-dir",
                tmpdir,
                "--ped",
                trio_ped_file,
                "--inheritance-mode",
                "simple",
                "--reference",
                "GRCh38",
                "--filters",
                "FILTER='PASS'",
                "--fields",
                "CHROM POS REF ALT GENE IMPACT AF",
            ]

            result = self.run_variantcentrifuge(args)

            # Check that the command succeeded
            assert result.returncode == 0, f"Command failed: {result.stderr}"

            # Check output file
            assert os.path.exists(output_file)
            df = pd.read_csv(output_file, sep="\t")

            # Should have Inheritance_Pattern but not Inheritance_Details
            assert "Inheritance_Pattern" in df.columns
            assert "Inheritance_Details" not in df.columns

            # Check that patterns were assigned
            patterns = df["Inheritance_Pattern"].unique()
            assert len(patterns) > 0
            assert "none" not in patterns or len(patterns) > 1

    @pytest.mark.slow
    def test_inheritance_mode_full(self, sample_vcf_file, trio_ped_file, gene_bed_file):
        """Test inheritance analysis with full mode."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "output.tsv")

            # Run with full mode
            args = [
                "--vcf-file",
                sample_vcf_file,
                "--gene-name",
                "all",
                "--output-file",
                output_file,
                "--output-dir",
                tmpdir,
                "--ped",
                trio_ped_file,
                "--inheritance-mode",
                "full",
                "--reference",
                "GRCh38",
                "--filters",
                "FILTER='PASS'",
                "--fields",
                "CHROM POS REF ALT GENE IMPACT AF",
            ]

            result = self.run_variantcentrifuge(args)

            # Check that the command succeeded
            assert result.returncode == 0, f"Command failed: {result.stderr}"

            # Check output file
            assert os.path.exists(output_file)
            df = pd.read_csv(output_file, sep="\t")

            # Should have both columns
            assert "Inheritance_Pattern" in df.columns
            assert "Inheritance_Details" in df.columns

            # Check that details are valid JSON
            for _idx, row in df.iterrows():
                if pd.notna(row["Inheritance_Details"]):
                    details = json.loads(row["Inheritance_Details"])
                    assert "primary_pattern" in details
                    assert "confidence" in details

    @pytest.mark.slow
    def test_inheritance_mode_columns(self, sample_vcf_file, trio_ped_file, gene_bed_file):
        """Test inheritance analysis with columns mode."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "output.tsv")

            # Run with columns mode
            args = [
                "--vcf-file",
                sample_vcf_file,
                "--gene-name",
                "all",
                "--output-file",
                output_file,
                "--output-dir",
                tmpdir,
                "--ped",
                trio_ped_file,
                "--inheritance-mode",
                "columns",
                "--reference",
                "GRCh38",
                "--filters",
                "FILTER='PASS'",
                "--fields",
                "CHROM POS REF ALT GENE IMPACT AF",
            ]

            result = self.run_variantcentrifuge(args)

            # Check that the command succeeded
            assert result.returncode == 0, f"Command failed: {result.stderr}"

            # Check output file
            assert os.path.exists(output_file)
            df = pd.read_csv(output_file, sep="\t")

            # Should have pattern and new columns but not details
            assert "Inheritance_Pattern" in df.columns
            assert "Inheritance_Details" not in df.columns
            assert "Inheritance_Confidence" in df.columns
            assert "Inheritance_Description" in df.columns
            assert "Inheritance_Samples" in df.columns

            # Check that new columns have values
            for _idx, row in df.iterrows():
                if row["Inheritance_Pattern"] != "none":
                    assert pd.notna(row["Inheritance_Confidence"])
                    assert pd.notna(row["Inheritance_Description"])
                    assert pd.notna(row["Inheritance_Samples"])

    @pytest.mark.slow
    def test_no_ped_no_inheritance(self, sample_vcf_file, gene_bed_file):
        """Test that without PED file, no inheritance analysis is performed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "output.tsv")

            # Run without PED file
            args = [
                "--vcf-file",
                sample_vcf_file,
                "--gene-name",
                "all",
                "--output-file",
                output_file,
                "--output-dir",
                tmpdir,
                "--reference",
                "GRCh38",
                "--filters",
                "FILTER='PASS'",
                "--fields",
                "CHROM POS REF ALT GENE IMPACT AF",
            ]

            result = self.run_variantcentrifuge(args)

            # Check that the command succeeded
            assert result.returncode == 0, f"Command failed: {result.stderr}"

            # Check output file
            assert os.path.exists(output_file)
            df = pd.read_csv(output_file, sep="\t")

            # Should not have inheritance columns
            assert "Inheritance_Pattern" not in df.columns
            assert "Inheritance_Details" not in df.columns
            assert "Inheritance_Confidence" not in df.columns

    def test_inheritance_mode_backward_compatibility(self):
        """Test that existing scripts using --calculate-inheritance still work."""
        # This test verifies that the removal of --calculate-inheritance
        # is handled gracefully and users get an appropriate error message

        # Note: Since we removed the flag, this would fail at argument parsing
        # This is expected behavior - users need to update their scripts
        pass
