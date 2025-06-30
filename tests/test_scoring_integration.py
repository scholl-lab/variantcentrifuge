# File: tests/test_scoring_integration.py

import os
import tempfile
import shutil
import subprocess
import pandas as pd
import pytest

# Check if we can run the full pipeline
try:
    import statsmodels  # noqa: F401

    STATSMODELS_AVAILABLE = True
except ImportError:
    STATSMODELS_AVAILABLE = False

# Skip all tests in this module if statsmodels is not available
pytestmark = pytest.mark.skipif(
    not STATSMODELS_AVAILABLE, reason="statsmodels is required for integration tests"
)

# Get the directory where test files are located
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
FIXTURES_DIR = os.path.join(TEST_DIR, "fixtures")
PROJECT_ROOT = os.path.dirname(TEST_DIR)


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for test outputs."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Cleanup after test
    shutil.rmtree(temp_dir)


@pytest.fixture
def annotated_vcf():
    """Path to the annotated test VCF file."""
    vcf_path = os.path.join(FIXTURES_DIR, "test_variants_scoring.GRCh37.annotated.vcf.gz")
    if not os.path.exists(vcf_path):
        pytest.skip(f"Annotated VCF fixture not found: {vcf_path}")
    return vcf_path


@pytest.fixture
def scoring_config_dir():
    """Path to the scoring configuration directory."""
    return os.path.join(PROJECT_ROOT, "scoring", "nephro_variant_score")


def test_scoring_integration_with_annotated_vcf(annotated_vcf, scoring_config_dir, temp_output_dir):
    """Test the full scoring pipeline with an annotated VCF file."""

    # Define output file
    output_file = os.path.join(temp_output_dir, "scored_variants.tsv")

    # Build the command
    cmd = [
        "variantcentrifuge",
        "--vcf-file",
        annotated_vcf,
        "--gene-name",
        "all",  # Required argument
        "--reference",
        "GRCh37.75",
        "--fields",
        (
            "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT ANN[0].EFFECT "
            "dbNSFP_CADD_phred dbNSFP_gnomAD_exomes_AF dbNSFP_gnomAD_genomes_AF GT"
        ),
        "--filters",
        "ANN[0].IMPACT != 'MODIFIER'",  # Use proper annotation syntax
        "--output-file",
        output_file,
        "--output-dir",
        temp_output_dir,
        "--scoring-config-path",
        scoring_config_dir,
        "--no-replacement",  # Skip genotype replacement for simpler test
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command succeeded
    assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"

    # Check that output file was created
    assert os.path.exists(output_file), "Output file was not created"

    # Read the output and verify scoring columns were added
    df = pd.read_csv(output_file, sep="\t")

    # Check that the scoring column exists
    assert (
        "nephro_variant_score" in df.columns
    ), "Scoring column 'nephro_variant_score' not found in output"

    # Check that scores were calculated (not all NaN)
    assert not df["nephro_variant_score"].isna().all(), "All scores are NaN"

    # Check that scores are in expected range (0 to 1 for logistic model)
    valid_scores = df["nephro_variant_score"].dropna()
    assert (valid_scores >= 0).all() and (
        valid_scores <= 1
    ).all(), "Scores outside expected range [0, 1]"

    # Verify that the formula used the expected columns
    # Note: The scoring module renames columns according to the variable mapping
    # IMPACT -> impact_variant, EFFECT -> consequence_terms_variant
    expected_column_mappings = {
        "GENE": "GENE",  # GENE is not renamed
        "IMPACT": "impact_variant",  # Renamed by scoring config
        "EFFECT": "consequence_terms_variant",  # Renamed by scoring config
    }

    for original, expected in expected_column_mappings.items():
        assert (
            expected in df.columns
        ), f"Expected column '{expected}' (from {original}) not found in output"


def test_scoring_with_missing_columns(annotated_vcf, scoring_config_dir, temp_output_dir):
    """Test scoring behavior when some expected columns are missing."""

    output_file = os.path.join(temp_output_dir, "scored_variants_minimal.tsv")

    # Use minimal fields that don't include all scoring variables
    cmd = [
        "variantcentrifuge",
        "--vcf-file",
        annotated_vcf,
        "--gene-name",
        "all",  # Required argument
        "--reference",
        "GRCh37.75",
        "--fields",
        "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT GT",  # Missing CADD and gnomAD fields
        "--filters",
        "ANN[0].IMPACT != 'MODIFIER'",
        "--output-file",
        output_file,
        "--output-dir",
        temp_output_dir,
        "--scoring-config-path",
        scoring_config_dir,
        "--no-replacement",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Command should still succeed (scoring uses defaults for missing columns)
    assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"

    # Check output
    assert os.path.exists(output_file), "Output file was not created"

    df = pd.read_csv(output_file, sep="\t")
    assert "nephro_variant_score" in df.columns, "Scoring column not found"

    # Scores should still be calculated using default values
    assert not df["nephro_variant_score"].isna().all(), "All scores are NaN when using defaults"


def test_scoring_with_all_genes(annotated_vcf, scoring_config_dir, temp_output_dir):
    """Test scoring with gene set to 'all' to process all variants."""

    output_file = os.path.join(temp_output_dir, "scored_all_genes.tsv")

    cmd = [
        "variantcentrifuge",
        "--vcf-file",
        annotated_vcf,
        "--gene-name",
        "all",  # Process all genes
        "--reference",
        "GRCh37.75",
        "--fields",
        (
            "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT ANN[0].EFFECT "
            "dbNSFP_CADD_phred dbNSFP_gnomAD_exomes_AF GT"
        ),
        "--filters",
        "ANN[0].IMPACT != 'MODIFIER'",
        "--output-file",
        output_file,
        "--output-dir",
        temp_output_dir,
        "--scoring-config-path",
        scoring_config_dir,
        "--no-replacement",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"

    # Verify output
    assert os.path.exists(output_file), "Output file was not created"
    df = pd.read_csv(output_file, sep="\t")

    # Should have multiple genes (might have fewer if some are filtered out)
    assert len(df) > 0, "Expected at least some variants in output"
    # Check for the GENE column - it will be present even if some genes are filtered
    assert "GENE" in df.columns, "GENE column not found in output"

    # All variants should have scores
    assert "nephro_variant_score" in df.columns, "Scoring column not found"
    assert len(df) > 0, "No variants in output"


@pytest.mark.parametrize(
    "impact_filter,expected_min_variants",
    [
        ("ANN[0].IMPACT = 'HIGH'", 1),  # At least some HIGH impact variants
        ("ANN[0].IMPACT IN ('HIGH', 'MODERATE')", 2),  # More variants with HIGH or MODERATE
        ("1", 5),  # No filter (always true), should have many variants
    ],
)
def test_scoring_with_different_filters(
    annotated_vcf, scoring_config_dir, temp_output_dir, impact_filter, expected_min_variants
):
    """Test scoring with different filter conditions."""

    output_file = os.path.join(
        temp_output_dir, f"scored_filtered_{impact_filter.replace(' ', '_')}.tsv"
    )

    cmd = [
        "variantcentrifuge",
        "--vcf-file",
        annotated_vcf,
        "--gene-name",
        "all",  # Required argument
        "--reference",
        "GRCh37.75",
        "--fields",
        (
            "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT ANN[0].EFFECT "
            "dbNSFP_CADD_phred dbNSFP_gnomAD_exomes_AF GT"
        ),
        "--filters",
        impact_filter,
        "--output-file",
        output_file,
        "--output-dir",
        temp_output_dir,
        "--scoring-config-path",
        scoring_config_dir,
        "--no-replacement",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check success
    assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"
    assert os.path.exists(output_file), "Output file was not created"

    # Verify filtering worked and scoring was applied
    df = pd.read_csv(output_file, sep="\t")
    assert (
        len(df) >= expected_min_variants
    ), f"Expected at least {expected_min_variants} variants, got {len(df)}"
    assert "nephro_variant_score" in df.columns, "Scoring column not found"

    # Verify scores are calculated
    valid_scores = df["nephro_variant_score"].dropna()
    assert len(valid_scores) > 0, "No valid scores calculated"
