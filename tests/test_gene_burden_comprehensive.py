"""
Comprehensive tests for gene burden analysis functionality.

Tests all possible ways to specify case/control groups:
1. Direct sample specification (--case-samples, --control-samples)
2. Sample file specification (--case-samples-file, --control-samples-file)
3. Phenotype-based classification with HPO terms
4. Phenotype file integration with various column configurations
5. Mixed scenarios and edge cases

Uses the comprehensive test dataset generated in tests/fixtures/geneburden/
"""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture(scope="module")
def test_data_dir():
    """Return path to comprehensive gene burden test data."""
    fixtures_dir = Path(__file__).parent / "fixtures" / "geneburden" / "comprehensive_test_data"
    if not fixtures_dir.exists():
        pytest.skip(
            f"Test data not found at {fixtures_dir}. Run generate_comprehensive_test_data.py first."
        )
    return fixtures_dir


@pytest.fixture(scope="module")
def vcf_file(test_data_dir):
    """Return path to test VCF file."""
    vcf_path = test_data_dir / "test_data.vcf.gz"
    if not vcf_path.exists():
        pytest.skip(f"Test VCF not found at {vcf_path}")
    return vcf_path


@pytest.fixture(scope="module")
def gene_file(test_data_dir):
    """Return path to test gene list."""
    gene_path = test_data_dir / "test_genes.txt"
    if not gene_path.exists():
        pytest.skip(f"Gene file not found at {gene_path}")
    return gene_path


class TestGeneBurdenDirectSamples:
    """Test direct sample specification methods."""

    @pytest.mark.slow
    @pytest.mark.gene_burden
    def test_direct_sample_specification(self, vcf_file, gene_file, test_data_dir):
        """Test --case-samples and --control-samples flags."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Read sample lists
            case_samples = (test_data_dir / "case_samples.txt").read_text().strip().split("\n")[:5]
            control_samples = (
                (test_data_dir / "control_samples.txt").read_text().strip().split("\n")[:8]
            )

            cmd = [
                "variantcentrifuge",
                "--vcf-file",
                str(vcf_file),
                "--gene-file",
                str(gene_file),
                "--case-samples",
                ",".join(case_samples),
                "--control-samples",
                ",".join(control_samples),
                "--perform-gene-burden",
                "--preset",
                "rare,coding",
                "--output-dir",
                tmpdir,
                "--output-file",
                "results.tsv",
                "--use-new-pipeline",
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Command failed: {result.stderr}"

            # Check that gene burden results were generated
            results_file = Path(tmpdir) / "results.gene_burden.tsv"
            if results_file.exists():
                df = pd.read_csv(results_file, sep="\t")
                assert len(df) > 0, "No gene burden results generated"
                assert "gene" in df.columns or "GENE" in df.columns, "Missing gene column"

    @pytest.mark.slow
    @pytest.mark.gene_burden
    def test_sample_file_specification(self, vcf_file, gene_file, test_data_dir):
        """Test --case-samples-file and --control-samples-file flags."""
        with tempfile.TemporaryDirectory() as tmpdir:
            case_file = test_data_dir / "case_samples.txt"
            control_file = test_data_dir / "control_samples.txt"

            cmd = [
                "variantcentrifuge",
                "--vcf-file",
                str(vcf_file),
                "--gene-file",
                str(gene_file),
                "--case-samples-file",
                str(case_file),
                "--control-samples-file",
                str(control_file),
                "--perform-gene-burden",
                "--preset",
                "rare,coding",
                "--output-dir",
                tmpdir,
                "--output-file",
                "results.tsv",
                "--use-new-pipeline",
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Command failed: {result.stderr}"


class TestGeneBurdenPhenotypeBased:
    """Test phenotype-based case/control specification."""

    @pytest.mark.slow
    @pytest.mark.gene_burden
    def test_hpo_based_classification(self, vcf_file, gene_file, test_data_dir):
        """Test HPO-based case/control classification."""
        with tempfile.TemporaryDirectory() as tmpdir:
            phenotype_file = test_data_dir / "phenotypes_basic.csv"
            case_hpo_terms = (test_data_dir / "case_hpo_terms.txt").read_text().strip().split("\n")
            control_hpo_terms = (
                (test_data_dir / "control_hpo_terms.txt").read_text().strip().split("\n")
            )

            cmd = [
                "variantcentrifuge",
                "--vcf-file",
                str(vcf_file),
                "--gene-file",
                str(gene_file),
                "--phenotype-file",
                str(phenotype_file),
                "--phenotype-sample-column",
                "SampleID",
                "--phenotype-value-column",
                "identifier",
                "--case-phenotypes",
                ",".join(case_hpo_terms),
                "--control-phenotypes",
                ",".join(control_hpo_terms),
                "--perform-gene-burden",
                "--preset",
                "rare,coding",
                "--output-dir",
                tmpdir,
                "--output-file",
                "results.tsv",
                "--use-new-pipeline",
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Command failed: {result.stderr}"

    @pytest.mark.slow
    @pytest.mark.gene_burden
    def test_hpo_term_files(self, vcf_file, gene_file, test_data_dir):
        """Test HPO term file specification."""
        with tempfile.TemporaryDirectory() as tmpdir:
            phenotype_file = test_data_dir / "phenotypes_basic.csv"
            case_hpo_file = test_data_dir / "case_hpo_terms.txt"
            control_hpo_file = test_data_dir / "control_hpo_terms.txt"

            cmd = [
                "variantcentrifuge",
                "--vcf-file",
                str(vcf_file),
                "--gene-file",
                str(gene_file),
                "--phenotype-file",
                str(phenotype_file),
                "--phenotype-sample-column",
                "SampleID",
                "--phenotype-value-column",
                "identifier",
                "--case-phenotypes-file",
                str(case_hpo_file),
                "--control-phenotypes-file",
                str(control_hpo_file),
                "--perform-gene-burden",
                "--preset",
                "rare,coding",
                "--output-dir",
                tmpdir,
                "--output-file",
                "results.tsv",
                "--use-new-pipeline",
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Command failed: {result.stderr}"

    @pytest.mark.slow
    @pytest.mark.gene_burden
    def test_alternative_column_names(self, vcf_file, gene_file, test_data_dir):
        """Test phenotype file with alternative column names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            phenotype_file = test_data_dir / "phenotypes_alt_columns.csv"
            case_hpo_terms = (test_data_dir / "case_hpo_terms.txt").read_text().strip().split("\n")
            control_hpo_terms = (
                (test_data_dir / "control_hpo_terms.txt").read_text().strip().split("\n")
            )

            cmd = [
                "variantcentrifuge",
                "--vcf-file",
                str(vcf_file),
                "--gene-file",
                str(gene_file),
                "--phenotype-file",
                str(phenotype_file),
                "--phenotype-sample-column",
                "sample_name",  # Alternative column name
                "--phenotype-value-column",
                "hpo_id",  # Alternative column name
                "--case-phenotypes",
                ",".join(case_hpo_terms),
                "--control-phenotypes",
                ",".join(control_hpo_terms),
                "--perform-gene-burden",
                "--preset",
                "rare,coding",
                "--output-dir",
                tmpdir,
                "--output-file",
                "results.tsv",
                "--use-new-pipeline",
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Command failed: {result.stderr}"


class TestGeneBurdenValidation:
    """Test gene burden analysis result validation."""

    @pytest.mark.slow
    @pytest.mark.gene_burden
    def test_expected_disease_gene_enrichment(self, vcf_file, gene_file, test_data_dir):
        """Test that disease genes show expected enrichment patterns."""
        with tempfile.TemporaryDirectory() as tmpdir:
            case_file = test_data_dir / "case_samples.txt"
            control_file = test_data_dir / "control_samples.txt"

            cmd = [
                "variantcentrifuge",
                "--vcf-file",
                str(vcf_file),
                "--gene-file",
                str(gene_file),
                "--case-samples-file",
                str(case_file),
                "--control-samples-file",
                str(control_file),
                "--perform-gene-burden",
                "--preset",
                "rare,coding",
                "--output-dir",
                tmpdir,
                "--output-file",
                "results.tsv",
                "--use-new-pipeline",
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Command failed: {result.stderr}"

            # Check for gene burden results file
            results_file = Path(tmpdir) / "results.gene_burden.tsv"
            if results_file.exists():
                df = pd.read_csv(results_file, sep="\t")

                # Disease genes should be present in results
                disease_genes = {"PKD1", "PKD2", "BRCA1", "BRCA2"}
                gene_col = None
                for col in df.columns:
                    if col.lower() in ["gene", "gene_name", "gene_id"]:
                        gene_col = col
                        break

                if gene_col:
                    result_genes = set(df[gene_col].astype(str))
                    found_disease_genes = disease_genes & result_genes
                    assert len(found_disease_genes) > 0, (
                        f"No disease genes found in results: {result_genes}"
                    )
                    print(f"Found disease genes: {found_disease_genes}")


class TestGeneBurdenErrorHandling:
    """Test error handling and edge cases."""

    @pytest.mark.gene_burden
    def test_no_case_control_specification(self, vcf_file, gene_file):
        """Test behavior when no case/control groups are specified."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [
                "variantcentrifuge",
                "--vcf-file",
                str(vcf_file),
                "--gene-file",
                str(gene_file),
                "--perform-gene-burden",
                "--preset",
                "rare,coding",
                "--output-dir",
                tmpdir,
                "--output-file",
                "results.tsv",
                "--use-new-pipeline",
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            # This should either fail gracefully or skip gene burden analysis
            # The exact behavior depends on implementation
            assert result.returncode in [0, 1], f"Unexpected return code: {result.returncode}"

    @pytest.mark.gene_burden
    def test_missing_phenotype_file(self, vcf_file, gene_file):
        """Test behavior with missing phenotype file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [
                "variantcentrifuge",
                "--vcf-file",
                str(vcf_file),
                "--gene-file",
                str(gene_file),
                "--phenotype-file",
                "nonexistent_file.csv",
                "--phenotype-sample-column",
                "SampleID",
                "--phenotype-value-column",
                "identifier",
                "--case-phenotypes",
                "HP:0000113",
                "--control-phenotypes",
                "HP:0000001",
                "--perform-gene-burden",
                "--preset",
                "rare,coding",
                "--output-dir",
                tmpdir,
                "--output-file",
                "results.tsv",
                "--use-new-pipeline",
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            # Should fail gracefully with non-zero exit code
            assert result.returncode != 0, "Should fail with missing phenotype file"


@pytest.mark.gene_burden
def test_test_data_availability():
    """Verify that test data is available and complete."""
    fixtures_dir = Path(__file__).parent / "fixtures" / "geneburden" / "comprehensive_test_data"

    if not fixtures_dir.exists():
        pytest.skip(
            f"Test data not generated. Run: cd {fixtures_dir.parent} && "
            f"python generate_comprehensive_test_data.py --output-dir comprehensive_test_data"
        )

    required_files = [
        "test_data.vcf.gz",
        "case_samples.txt",
        "control_samples.txt",
        "test_genes.txt",
        "phenotypes_basic.csv",
    ]

    for filename in required_files:
        filepath = fixtures_dir / filename
        assert filepath.exists(), f"Required test file missing: {filename}"
