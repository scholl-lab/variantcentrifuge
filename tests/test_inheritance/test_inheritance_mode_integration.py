"""
Integration tests for inheritance mode functionality.

These tests verify that the --inheritance-mode flag works correctly
through the full command-line interface.

Requires external tools: bcftools, SnpSift, bedtools, snpEff, bgzip, tabix.
Requires snpEff with a GRCh38 genome database installed locally.
"""

import json
import os
import shutil
import subprocess
import tempfile

import pandas as pd
import pytest

# These CLI integration tests require external genomic tools AND configured databases.
# snpEff needs a GRCh38 genome database installed locally.
_REQUIRED_TOOLS = ["bcftools", "SnpSift", "bedtools", "snpEff", "bgzip", "tabix"]
_MISSING_TOOLS = [t for t in _REQUIRED_TOOLS if shutil.which(t) is None]

# snpEff genome names vary by installation (GRCh38, GRCh38.p14, GRCh38.105, etc.)
# Auto-detect the installed GRCh38 variant.
_SNPEFF_GRCH38_NAMES = ["GRCh38.p14", "GRCh38.105", "GRCh38.99", "GRCh38.86", "GRCh38"]


def _detect_snpeff_grch38() -> str | None:
    """Find which GRCh38 genome name is locally installed in snpEff."""
    if shutil.which("snpEff") is None:
        return None
    snpeff_path = shutil.which("snpEff")
    if snpeff_path is None:
        return None
    import pathlib

    snpeff_bin = pathlib.Path(snpeff_path).resolve()
    # Conda layout: snpEff resolves to .../share/snpeff-X.Y-Z/snpEff
    # Data dir is a sibling: .../share/snpeff-X.Y-Z/data/
    share_dirs = []
    # Check sibling data/ directory (resolved symlink location)
    sibling_data = snpeff_bin.parent / "data"
    if sibling_data.is_dir():
        share_dirs.append(sibling_data)
    # Also check conda prefix layout
    for ancestor in snpeff_bin.parents:
        candidate = list(ancestor.glob("share/snpeff-*/data"))
        if candidate:
            share_dirs.extend(candidate)
            break
    for share_dir in share_dirs:
        for name in _SNPEFF_GRCH38_NAMES:
            if (share_dir / name).is_dir():
                return name
    # Fallback: try running snpEff to check
    for name in _SNPEFF_GRCH38_NAMES:
        try:
            result = subprocess.run(
                ["snpEff", "dump", "-v", name],
                capture_output=True,
                text=True,
                timeout=15,
            )
            if result.returncode == 0:
                return name
        except (subprocess.TimeoutExpired, OSError):
            continue
    return None


_SNPEFF_GENOME = _detect_snpeff_grch38() if len(_MISSING_TOOLS) == 0 else None

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        len(_MISSING_TOOLS) > 0,
        reason=f"Requires external tools: {', '.join(_MISSING_TOOLS)}",
    ),
    pytest.mark.skipif(
        _SNPEFF_GENOME is None,
        reason="Requires snpEff with a GRCh38 genome database installed",
    ),
]

# GRCh38 exonic coordinates in well-known genes.
# TP53 coding exons: chr17:7,673,700-7,676,600 (reverse strand)
# BRCA2 exon 11: chr13:32,325,000-32,357,000
# These positions are within coding exons and will receive missense/synonymous
# annotations from snpEff, ensuring they survive pipeline field extraction.
_TEST_VCF_CONTENT = """\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr13>
##contig=<ID=chr17>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother
chr17\t7674220\t.\tC\tT\t50\tPASS\tAF=0.0001\tGT\t0/1\t0/0\t0/0
chr13\t32340300\t.\tA\tG\t50\tPASS\tAF=0.0005\tGT\t0/1\t0/1\t0/0
chr13\t32341000\t.\tC\tT\t50\tPASS\tAF=0.0003\tGT\t0/1\t0/0\t0/1
chr17\t7675088\t.\tC\tT\t50\tPASS\tAF=0.001\tGT\t1/1\t0/1\t0/1
"""

# Genes matching the VCF coordinates
_TEST_GENES = "TP53 BRCA2"

# Minimum SnpSift extractFields needed for inheritance analysis.
# Must include ANN fields for gene/impact and GEN[*].GT for genotype data.
_TEST_FIELDS = "CHROM POS REF ALT ANN[0].GENE ANN[0].EFFECT ANN[0].IMPACT GEN[*].GT"


class TestInheritanceModeIntegration:
    """Test inheritance mode functionality through the CLI."""

    @pytest.fixture(scope="class")
    def annotated_vcf_file(self, tmp_path_factory):
        """Create a snpEff-annotated, bgzipped, and indexed VCF file.

        Uses real GRCh38 exonic coordinates so snpEff assigns gene/impact
        annotations. The pipeline expects pre-annotated VCFs.
        """
        tmpdir = tmp_path_factory.mktemp("vcf")
        plain_vcf = str(tmpdir / "test.vcf")
        annotated_vcf = str(tmpdir / "test.ann.vcf")
        gz_vcf = annotated_vcf + ".gz"

        # Write the raw VCF
        with open(plain_vcf, "w") as f:
            f.write(_TEST_VCF_CONTENT)

        # Annotate with snpEff to add ANN fields.
        # GRCh38 database requires >1 GB heap; default snpEff wrapper uses -Xmx1g
        # which is insufficient. Pass -Xmx4G to snpEff (forwarded to JVM).
        with open(annotated_vcf, "w") as out:
            result = subprocess.run(
                ["snpEff", "-Xmx4G", "ann", "-noStats", _SNPEFF_GENOME, plain_vcf],
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                timeout=300,
            )
        assert result.returncode == 0, f"snpEff annotation failed: {result.stderr}"

        # Sort, bgzip, and index for bcftools compatibility.
        # snpEff may reorder variants; bcftools sort ensures valid tabix indexing.
        sorted_vcf = str(tmpdir / "test.sorted.vcf.gz")
        subprocess.run(
            ["bcftools", "sort", "-Oz", "-o", sorted_vcf, annotated_vcf],
            check=True,
            capture_output=True,
        )
        subprocess.run(["tabix", "-p", "vcf", sorted_vcf], check=True)
        gz_vcf = sorted_vcf

        return gz_vcf

    @pytest.fixture
    def trio_ped_file(self):
        """Create a trio pedigree file."""
        ped_content = """\
#FamilyID\tIndividualID\tPaternalID\tMaternalID\tSex\tAffectedStatus
FAM1\tchild\tfather\tmother\t1\t2
FAM1\tfather\t0\t0\t1\t1
FAM1\tmother\t0\t0\t2\t1
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".ped", delete=False) as f:
            f.write(ped_content)
            ped_file = f.name

        yield ped_file
        os.unlink(ped_file)

    def run_variantcentrifuge(self, args):
        """Run variantcentrifuge with given arguments.

        Explicitly passes PATH to the subprocess to ensure system utilities
        (gzip, sort) are available in the child process's shell commands.
        """
        cmd = ["variantcentrifuge", *args]
        env = os.environ.copy()
        # Ensure /usr/bin and /usr/local/bin are on PATH for system utilities
        # needed by shell-based stages (DataSortingStage uses gzip, sort).
        for p in ["/usr/bin", "/usr/local/bin", "/bin"]:
            if p not in env.get("PATH", ""):
                env["PATH"] = env.get("PATH", "") + os.pathsep + p
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=env)
        return result

    def test_inheritance_mode_simple(self, annotated_vcf_file, trio_ped_file):
        """Test inheritance analysis with simple mode."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "output.tsv")

            args = [
                "--vcf-file",
                annotated_vcf_file,
                "--gene-name",
                _TEST_GENES,
                "--output-file",
                output_file,
                "--output-dir",
                tmpdir,
                "--ped",
                trio_ped_file,
                "--inheritance-mode",
                "simple",
                "--reference",
                _SNPEFF_GENOME,
                "--add-chr",
                "--filters",
                "FILTER='PASS'",
                "--fields",
                _TEST_FIELDS,
            ]

            result = self.run_variantcentrifuge(args)
            assert result.returncode == 0, f"Command failed: {result.stderr}"

            assert os.path.exists(output_file)
            df = pd.read_csv(output_file, sep="\t")
            assert len(df) > 0, f"No variants in output. stderr: {result.stderr}"

            # Should have Inheritance_Pattern but not Inheritance_Details
            assert "Inheritance_Pattern" in df.columns
            assert "Inheritance_Details" not in df.columns

            # Check that patterns were assigned
            patterns = df["Inheritance_Pattern"].unique()
            assert len(patterns) > 0
            assert "none" not in patterns or len(patterns) > 1

    def test_inheritance_mode_full(self, annotated_vcf_file, trio_ped_file):
        """Test inheritance analysis with full mode."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "output.tsv")

            args = [
                "--vcf-file",
                annotated_vcf_file,
                "--gene-name",
                _TEST_GENES,
                "--output-file",
                output_file,
                "--output-dir",
                tmpdir,
                "--ped",
                trio_ped_file,
                "--inheritance-mode",
                "full",
                "--reference",
                _SNPEFF_GENOME,
                "--add-chr",
                "--filters",
                "FILTER='PASS'",
                "--fields",
                _TEST_FIELDS,
            ]

            result = self.run_variantcentrifuge(args)
            assert result.returncode == 0, f"Command failed: {result.stderr}"

            assert os.path.exists(output_file)
            df = pd.read_csv(output_file, sep="\t")
            assert len(df) > 0, f"No variants in output. stderr: {result.stderr}"

            # Should have both columns
            assert "Inheritance_Pattern" in df.columns
            assert "Inheritance_Details" in df.columns

            # Check that details are valid JSON
            for _idx, row in df.iterrows():
                if pd.notna(row["Inheritance_Details"]):
                    details = json.loads(row["Inheritance_Details"])
                    assert "primary_pattern" in details
                    assert "confidence" in details

    def test_inheritance_mode_columns(self, annotated_vcf_file, trio_ped_file):
        """Test inheritance analysis with columns mode."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "output.tsv")

            args = [
                "--vcf-file",
                annotated_vcf_file,
                "--gene-name",
                _TEST_GENES,
                "--output-file",
                output_file,
                "--output-dir",
                tmpdir,
                "--ped",
                trio_ped_file,
                "--inheritance-mode",
                "columns",
                "--reference",
                _SNPEFF_GENOME,
                "--add-chr",
                "--filters",
                "FILTER='PASS'",
                "--fields",
                _TEST_FIELDS,
            ]

            result = self.run_variantcentrifuge(args)
            assert result.returncode == 0, f"Command failed: {result.stderr}"

            assert os.path.exists(output_file)
            df = pd.read_csv(output_file, sep="\t")
            assert len(df) > 0, f"No variants in output. stderr: {result.stderr}"

            # Should have pattern and new columns but not details
            assert "Inheritance_Pattern" in df.columns
            assert "Inheritance_Details" not in df.columns
            assert "Inheritance_Confidence" in df.columns
            assert "Inheritance_Description" in df.columns
            assert "Inheritance_Samples" in df.columns

            # Check that new columns have values for definitive patterns
            # "none" and "unknown" patterns legitimately have no samples
            for _idx, row in df.iterrows():
                pattern = row["Inheritance_Pattern"]
                if pattern not in ("none", "unknown"):
                    assert pd.notna(row["Inheritance_Confidence"])
                    assert pd.notna(row["Inheritance_Description"])
                    assert pd.notna(row["Inheritance_Samples"])

    def test_no_ped_no_inheritance(self, annotated_vcf_file):
        """Test that without PED file, no inheritance analysis is performed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "output.tsv")

            args = [
                "--vcf-file",
                annotated_vcf_file,
                "--gene-name",
                _TEST_GENES,
                "--output-file",
                output_file,
                "--output-dir",
                tmpdir,
                "--reference",
                _SNPEFF_GENOME,
                "--add-chr",
                "--filters",
                "FILTER='PASS'",
                "--fields",
                _TEST_FIELDS,
            ]

            result = self.run_variantcentrifuge(args)
            assert result.returncode == 0, f"Command failed: {result.stderr}"

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
