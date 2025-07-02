"""Test the inheritance scoring functionality."""

import json
import pandas as pd
import pytest
from variantcentrifuge.scoring import apply_scoring, read_scoring_config
import os

# Get the project root directory
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(TEST_DIR)


@pytest.fixture
def inheritance_score_config():
    """Load the inheritance scoring configuration."""
    config_path = os.path.join(PROJECT_ROOT, "scoring", "inheritance_score")
    if not os.path.exists(config_path):
        pytest.skip(f"Inheritance scoring config not found: {config_path}")
    return read_scoring_config(config_path)


@pytest.fixture
def sample_inheritance_data():
    """Create sample DataFrame with inheritance analysis results."""
    data = {
        "CHROM": ["1", "2", "3", "4", "5"],
        "POS": [1000, 2000, 3000, 4000, 5000],
        "GENE": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
        "Inheritance_Pattern": [
            "de_novo", 
            "autosomal_recessive", 
            "compound_heterozygous",
            "autosomal_dominant",
            "unknown"
        ],
        "Inheritance_Details": [
            # De novo with segregation p-value
            json.dumps({"segregation_p_value": 0.0001}),
            # Autosomal recessive without segregation (should get penalty)
            json.dumps({}),
            # Compound het with good segregation
            json.dumps({"segregation_p_value": 0.005}),
            # Autosomal dominant with moderate segregation
            json.dumps({"segregation_p_value": 0.05}),
            # Unknown pattern (no penalty for missing segregation)
            json.dumps({})
        ]
    }
    return pd.DataFrame(data)


def test_inheritance_scoring_basic(sample_inheritance_data, inheritance_score_config):
    """Test basic inheritance scoring functionality."""
    scored_df = apply_scoring(sample_inheritance_data, inheritance_score_config)
    
    # Check that only the final score column is present (not intermediates)
    assert "inheritance_score" in scored_df.columns
    assert "base_score" not in scored_df.columns
    assert "segregation_p_value" not in scored_df.columns
    assert "enhancement_factor" not in scored_df.columns
    assert "enhanced_score" not in scored_df.columns
    
    # Verify scores are in expected range
    scores = scored_df["inheritance_score"]
    assert (scores >= 0).all() and (scores <= 1).all()


def test_inheritance_scoring_values(sample_inheritance_data, inheritance_score_config):
    """Test specific inheritance score calculations."""
    scored_df = apply_scoring(sample_inheritance_data, inheritance_score_config)
    scores = scored_df["inheritance_score"]
    
    # De novo should have high score (0.95 base, no penalty since it has segregation)
    assert abs(scores[0] - 0.95) < 0.01  # Base score with no penalty
    
    # Autosomal recessive without segregation should be penalized
    # Base is 0.8, with 0.2 penalty factor: 0.8 * (1 - 0.2) = 0.64
    assert abs(scores[1] - 0.64) < 0.01
    
    # Compound het with segregation should have base score (no penalty)
    assert abs(scores[2] - 0.8) < 0.01  # Base 0.8, no penalty
    
    # Autosomal dominant with segregation has base score (no penalty)
    assert abs(scores[3] - 0.4) < 0.01  # Base score
    
    # Unknown pattern has low score but no penalty (exempt from penalty)
    assert abs(scores[4] - 0.1) < 0.01  # Base score of 0.1


def test_inheritance_scoring_with_missing_details():
    """Test inheritance scoring when Inheritance_Details is missing or malformed."""
    data = {
        "CHROM": ["1", "2"],
        "POS": [1000, 2000],
        "GENE": ["GENE1", "GENE2"],
        "Inheritance_Pattern": ["de_novo", "autosomal_recessive"],
        "Inheritance_Details": [
            "not valid json",  # Invalid JSON
            ""  # Empty string
        ]
    }
    df = pd.DataFrame(data)
    
    config_path = os.path.join(PROJECT_ROOT, "scoring", "inheritance_score")
    if not os.path.exists(config_path):
        pytest.skip("Inheritance scoring config not found")
    
    config = read_scoring_config(config_path)
    
    # Should not crash on invalid data
    scored_df = apply_scoring(df, config)
    assert "inheritance_score" in scored_df.columns
    
    # De novo without segregation data should still get base score
    assert abs(scored_df["inheritance_score"][0] - 0.95) < 0.01


def test_compound_het_possible_patterns():
    """Test scoring of compound_heterozygous_possible patterns."""
    data = {
        "CHROM": ["1", "2"],
        "POS": [1000, 2000],
        "GENE": ["GENE1", "GENE2"],
        "Inheritance_Pattern": [
            "compound_heterozygous_possible",
            "compound_heterozygous_possible_no_pedigree"
        ],
        "Inheritance_Details": [
            json.dumps({}),  # No segregation data
            json.dumps({})   # No segregation data
        ]
    }
    df = pd.DataFrame(data)
    
    config_path = os.path.join(PROJECT_ROOT, "scoring", "inheritance_score")
    if not os.path.exists(config_path):
        pytest.skip("Inheritance scoring config not found")
    
    config = read_scoring_config(config_path)
    scored_df = apply_scoring(df, config)
    
    # Both should get base score of 0.4
    # First one gets penalty (0.4 * 0.8 = 0.32)
    assert abs(scored_df["inheritance_score"][0] - 0.32) < 0.01
    
    # Second one is exempt from penalty (keeps 0.4)
    assert abs(scored_df["inheritance_score"][1] - 0.4) < 0.01