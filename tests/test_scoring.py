"""Tests for scoring module."""

import pandas as pd
import pytest
from variantcentrifuge.scoring import apply_scoring


@pytest.fixture
def sample_df():
    """Provide a sample DataFrame for testing."""
    data = {
        "dbNSFP_CADD_phred": [10.0, 25.0, 5.0],
        "dbNSFP_gnomAD_exomes_AF": [0.1, 0.001, 0.0],
        "IMPACT": ["MODERATE", "HIGH", "LOW"],
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_scoring_config():
    """Provide a sample scoring configuration."""
    return {
        "variables": {
            "dbNSFP_CADD_phred": "cadd",
            "dbNSFP_gnomAD_exomes_AF": "freq",
            "IMPACT": "impact_val",
        },
        "formulas": [
            {"pathogenicity_score": "cadd + (1 - freq) * 10"},
            {
                "impact_score": (
                    "((impact_val == 'HIGH') * 10) + ((impact_val == 'MODERATE') * 5) + "
                    "((impact_val == 'LOW') * 1)"
                )
            },
        ],
    }


def test_apply_scoring_calculates_correct_scores(sample_df, sample_scoring_config):
    """Test that apply_scoring adds correct score columns to the DataFrame."""
    scored_df = apply_scoring(sample_df, sample_scoring_config)

    assert "pathogenicity_score" in scored_df.columns
    assert "impact_score" in scored_df.columns

    # Verify calculations
    # Row 1: 10.0 + (1 - 0.1) * 10 = 10 + 9 = 19.0
    # Row 2: 25.0 + (1 - 0.001) * 10 = 25 + 9.99 = 34.99
    # Row 3: 5.0 + (1 - 0.0) * 10 = 5 + 10 = 15.0
    pd.testing.assert_series_equal(
        scored_df["pathogenicity_score"], pd.Series([19.0, 34.99, 15.0]), check_names=False
    )

    # Verify impact scores
    pd.testing.assert_series_equal(
        scored_df["impact_score"], pd.Series([5, 10, 1]), check_names=False
    )


def test_apply_scoring_handles_missing_columns(sample_df, sample_scoring_config):
    """Test that scoring handles missing columns by using defaults."""
    # Remove a column that the scoring config expects
    df_missing_col = sample_df.drop(columns=["dbNSFP_gnomAD_exomes_AF"])

    # Modify config to provide a default
    config_with_default = sample_scoring_config.copy()
    config_with_default["variables"]["dbNSFP_gnomAD_exomes_AF"] = "freq|default:0.5"

    scored_df = apply_scoring(df_missing_col, config_with_default)

    assert "pathogenicity_score" in scored_df.columns
    # All rows should use the default freq of 0.5
    # Row 1: 10 + (1 - 0.5) * 10 = 15
    # Row 2: 25 + (1 - 0.5) * 10 = 30
    # Row 3: 5 + (1 - 0.5) * 10 = 10
    pd.testing.assert_series_equal(
        scored_df["pathogenicity_score"], pd.Series([15.0, 30.0, 10.0]), check_names=False
    )
