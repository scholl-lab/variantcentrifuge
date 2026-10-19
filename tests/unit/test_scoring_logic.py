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


def test_scoring_hides_intermediate_columns():
    """Test that scoring correctly hides intermediate columns when output_scores is specified."""
    # Create a simple DataFrame
    df = pd.DataFrame(
        {
            "value1": [1.0, 2.0, 3.0],
            "value2": [4.0, 5.0, 6.0],
        }
    )

    # Multi-step config with output_scores
    config = {
        "variables": {
            "value1": "v1",
            "value2": "v2",
        },
        "output_scores": ["final_score"],  # Only keep final_score
        "formulas": [
            {"intermediate1": "v1 * 2"},
            {"intermediate2": "v2 * 3"},
            {"final_score": "intermediate1 + intermediate2"},
        ],
    }

    scored_df = apply_scoring(df, config)

    # Check that only final_score was kept
    assert "final_score" in scored_df.columns
    assert "intermediate1" not in scored_df.columns
    assert "intermediate2" not in scored_df.columns

    # Verify the final score is correct
    # Row 1: (1 * 2) + (4 * 3) = 2 + 12 = 14
    # Row 2: (2 * 2) + (5 * 3) = 4 + 15 = 19
    # Row 3: (3 * 2) + (6 * 3) = 6 + 18 = 24
    pd.testing.assert_series_equal(
        scored_df["final_score"], pd.Series([14.0, 19.0, 24.0]), check_names=False
    )

    # Original columns should still be present
    assert "value1" in scored_df.columns
    assert "value2" in scored_df.columns


def test_scoring_backward_compatibility_without_output_scores():
    """Test that scoring keeps all columns when output_scores is not specified."""
    # Create a simple DataFrame
    df = pd.DataFrame(
        {
            "value1": [1.0, 2.0, 3.0],
            "value2": [4.0, 5.0, 6.0],
        }
    )

    # Multi-step config WITHOUT output_scores
    config = {
        "variables": {
            "value1": "v1",
            "value2": "v2",
        },
        "formulas": [
            {"intermediate1": "v1 * 2"},
            {"intermediate2": "v2 * 3"},
            {"final_score": "intermediate1 + intermediate2"},
        ],
    }

    scored_df = apply_scoring(df, config)

    # Check that all calculated columns are kept
    assert "final_score" in scored_df.columns
    assert "intermediate1" in scored_df.columns
    assert "intermediate2" in scored_df.columns

    # Verify all calculations
    pd.testing.assert_series_equal(
        scored_df["intermediate1"], pd.Series([2.0, 4.0, 6.0]), check_names=False
    )
    pd.testing.assert_series_equal(
        scored_df["intermediate2"], pd.Series([12.0, 15.0, 18.0]), check_names=False
    )
    pd.testing.assert_series_equal(
        scored_df["final_score"], pd.Series([14.0, 19.0, 24.0]), check_names=False
    )


def test_nephro_score_equal_af_reproducibility():
    """Verify legacy and v5 produce identical scores when exome AF == genome AF == joint AF."""
    from variantcentrifuge.scoring import read_scoring_config

    cfg_legacy = read_scoring_config("scoring/nephro_candidate_score")
    cfg_v5 = read_scoring_config("scoring/nephro_candidate_score_v5")

    df_legacy = pd.DataFrame(
        [
            {
                "dbNSFP_gnomAD_exomes_AF": 0.0005,
                "dbNSFP_gnomAD_genomes_AF": 0.0005,
                "dbNSFP_CADD_phred": 20.0,
                "EFFECT": "missense_variant",
                "IMPACT": "MODERATE",
                "ngs": 0.0,
                "Inheritance_Pattern": "unknown",
                "Inheritance_Details": "",
            }
        ]
    )
    df_v5 = pd.DataFrame(
        [
            {
                "dbNSFP_gnomAD4.1_joint_AF": 0.0005,
                "dbNSFP_CADD_phred": 20.0,
                "EFFECT": "missense_variant",
                "IMPACT": "MODERATE",
                "ngs": 0.0,
                "Inheritance_Pattern": "unknown",
                "Inheritance_Details": "",
            }
        ]
    )

    res_legacy = apply_scoring(df_legacy, cfg_legacy)
    res_v5 = apply_scoring(df_v5, cfg_v5)

    score_legacy = res_legacy["nephro_variant_score"].iloc[0]
    score_v5 = res_v5["nephro_variant_score"].iloc[0]

    assert pytest.approx(score_legacy, rel=1e-6) == 0.056695344
    assert pytest.approx(score_v5, rel=1e-6) == 0.056695344
    assert pytest.approx(score_legacy, rel=1e-9) == score_v5


def test_nephro_score_unequal_af_reproducibility():
    """Verify legacy and v5 behavior under unequal AF inputs (Issue #144)."""
    from variantcentrifuge.scoring import read_scoring_config

    cfg_legacy = read_scoring_config("scoring/nephro_candidate_score")
    cfg_v5 = read_scoring_config("scoring/nephro_candidate_score_v5")

    # Exome AF = 0, Genome AF = 0.001, Joint AF = 0.0005
    df_legacy = pd.DataFrame(
        [
            {
                "dbNSFP_gnomAD_exomes_AF": 0.0,
                "dbNSFP_gnomAD_genomes_AF": 0.001,
                "dbNSFP_CADD_phred": 20.0,
                "EFFECT": "missense_variant",
                "IMPACT": "MODERATE",
                "ngs": 0.0,
                "Inheritance_Pattern": "unknown",
                "Inheritance_Details": "",
            }
        ]
    )
    df_v5 = pd.DataFrame(
        [
            {
                "dbNSFP_gnomAD4.1_joint_AF": 0.0005,
                "dbNSFP_CADD_phred": 20.0,
                "EFFECT": "missense_variant",
                "IMPACT": "MODERATE",
                "ngs": 0.0,
                "Inheritance_Pattern": "unknown",
                "Inheritance_Details": "",
            }
        ]
    )

    res_legacy = apply_scoring(df_legacy, cfg_legacy)
    res_v5 = apply_scoring(df_v5, cfg_v5)

    score_legacy = res_legacy["nephro_variant_score"].iloc[0]
    score_v5 = res_v5["nephro_variant_score"].iloc[0]

    assert pytest.approx(score_legacy, rel=1e-6) == 0.443311077
    assert pytest.approx(score_v5, rel=1e-6) == 0.056695344
