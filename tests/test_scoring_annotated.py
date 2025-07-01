"""Tests for scoring with annotated data."""

import os
import gzip
import pandas as pd
import pytest
from variantcentrifuge.scoring import apply_scoring, read_scoring_config

# Get the directory where test files are located
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
FIXTURES_DIR = os.path.join(TEST_DIR, "fixtures")
PROJECT_ROOT = os.path.dirname(TEST_DIR)


@pytest.fixture
def scoring_config():
    """Load the nephro_variant_score configuration."""
    config_path = os.path.join(PROJECT_ROOT, "scoring", "nephro_variant_score")
    if not os.path.exists(config_path):
        pytest.skip(f"Scoring config not found: {config_path}")
    return read_scoring_config(config_path)


@pytest.fixture
def sample_annotated_data():
    """Create a sample DataFrame mimicking annotated VCF data."""
    # This mimics what would be extracted from the annotated VCF
    data = {
        "CHROM": ["1", "2", "4", "6", "9"],
        "POS": [45508799, 169139571, 15502427, 10400506, 14789040],
        "REF": ["A", "C", "A", "G", "T"],
        "ALT": ["AT", "T", "T", "A", "C"],
        "GENE": ["MUTYH", "ABCB11", "CC2D2A", "TFAP2A", "NFKB1"],
        "IMPACT": ["MODERATE", "HIGH", "MODERATE", "HIGH", "LOW"],
        "EFFECT": [
            "missense_variant",
            "stop_gained",
            "missense_variant",
            "frameshift_variant",
            "synonymous_variant",
        ],
        "dbNSFP_CADD_phred": [23.5, 35.0, 22.1, 28.9, 5.2],
        "dbNSFP_gnomAD_exomes_AF": [0.001, 0.0, 0.0005, 0.0001, 0.1],
        "dbNSFP_gnomAD_genomes_AF": [0.0008, 0.0, 0.0003, 0.0001, 0.09],
        "GT": ["Sample1:0/1", "Sample1:0/1", "Sample1:0/1", "Sample1:0/1", "Sample1:0/1"],
    }
    return pd.DataFrame(data)


def test_scoring_with_annotated_data(sample_annotated_data, scoring_config):
    """Test scoring on data that mimics annotated VCF output."""
    # Apply scoring
    scored_df = apply_scoring(sample_annotated_data, scoring_config)

    # Check that the score column was added
    assert "nephro_variant_score" in scored_df.columns

    # Check that scores were calculated
    assert not scored_df["nephro_variant_score"].isna().all()

    # Verify scores are in expected range (0 to 1 for logistic regression)
    scores = scored_df["nephro_variant_score"]
    assert (scores >= 0).all() and (scores <= 1).all()

    # Check that HIGH impact variants generally have higher scores
    # Note: The scoring module renames columns to match the config
    high_impact_scores = scored_df[scored_df["impact_variant"] == "HIGH"][
        "nephro_variant_score"
    ].mean()
    low_impact_scores = scored_df[scored_df["impact_variant"] == "LOW"][
        "nephro_variant_score"
    ].mean()
    assert high_impact_scores > low_impact_scores

    # Check that variants with higher CADD scores have higher pathogenicity scores
    # Use the renamed column 'cadd_phred_variant'
    high_cadd = scored_df[scored_df["cadd_phred_variant"] > 25]
    low_cadd = scored_df[scored_df["cadd_phred_variant"] < 10]
    if len(high_cadd) > 0 and len(low_cadd) > 0:
        assert high_cadd["nephro_variant_score"].mean() > low_cadd["nephro_variant_score"].mean()


def test_scoring_with_missing_annotations(sample_annotated_data, scoring_config):
    """Test scoring when some annotations are missing."""
    # Remove some annotation columns
    df_missing = sample_annotated_data.drop(columns=["dbNSFP_gnomAD_genomes_AF"])

    # Apply scoring - should still work with defaults
    scored_df = apply_scoring(df_missing, scoring_config)

    assert "nephro_variant_score" in scored_df.columns
    assert not scored_df["nephro_variant_score"].isna().all()


def test_scoring_formula_components(sample_annotated_data, scoring_config):
    """Test that the scoring formula components work correctly."""
    # Apply scoring
    scored_df = apply_scoring(sample_annotated_data, scoring_config)

    # Verify the formula was applied
    assert "nephro_variant_score" in scored_df.columns

    # Test specific cases
    # HIGH impact with high CADD and low frequency should have high score
    # Use renamed columns matching the config
    high_risk_mask = (
        (scored_df["impact_variant"] == "HIGH")
        & (scored_df["cadd_phred_variant"] > 30)
        & (scored_df["gnomade_variant"] < 0.001)
    )
    if high_risk_mask.any():
        high_risk_scores = scored_df[high_risk_mask]["nephro_variant_score"]
        assert (high_risk_scores > 0.5).all(), "High risk variants should have high scores"


def test_read_scoring_config():
    """Test that the scoring configuration can be read correctly."""
    config_path = os.path.join(PROJECT_ROOT, "scoring", "nephro_variant_score")

    if not os.path.exists(config_path):
        pytest.skip("Scoring config directory not found")

    config = read_scoring_config(config_path)

    # Check structure
    assert "variables" in config
    assert "formulas" in config

    # Check that we have the expected variables
    variables = config["variables"]
    assert "dbNSFP_CADD_phred" in variables
    assert "dbNSFP_gnomAD_exomes_AF" in variables
    # Updated config uses ANN[0].IMPACT or IMPACT depending on the source
    assert any(key in variables for key in ["IMPACT", "ANN[0].IMPACT"])

    # Check that we have at least one formula
    assert len(config["formulas"]) > 0
    assert "nephro_variant_score" in config["formulas"][0]


def test_scoring_with_real_vcf_sample():
    """Test scoring with a sample from the real annotated VCF."""
    vcf_path = os.path.join(FIXTURES_DIR, "test_variants_scoring.GRCh37.annotated.vcf.gz")

    if not os.path.exists(vcf_path):
        pytest.skip("Annotated VCF not found")

    # Read a few lines from the VCF to create a test DataFrame
    # This is a simplified version - in reality you'd parse the VCF properly
    # For this test, we'll just verify the file exists and is readable
    try:
        with gzip.open(vcf_path, "rt") as f:
            lines = []
            for line in f:
                if not line.startswith("#"):
                    lines.append(line.strip())
                    if len(lines) >= 5:  # Just get a few variants
                        break

        assert len(lines) > 0, "No variant lines found in VCF"

    except Exception as e:
        pytest.fail(f"Failed to read annotated VCF: {e}")
