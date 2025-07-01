"""
Test that scoring module doesn't leave temporary variable columns in output.

This test verifies that gnomAD variable columns (gnomade_variant, gnomadg_variant)
created during formula evaluation are removed from the final output.
"""

import json
import os
import tempfile

import pandas as pd
import pytest

from variantcentrifuge.scoring import apply_scoring, read_scoring_config


class TestScoringColumns:
    """Test that scoring module properly manages columns."""

    @pytest.fixture
    def sample_df(self):
        """Create a sample DataFrame with gnomAD columns."""
        return pd.DataFrame(
            {
                "CHROM": ["1", "2"],
                "POS": [1000, 2000],
                "REF": ["A", "G"],
                "ALT": ["T", "C"],
                "GENE": ["GENE1", "GENE2"],
                "dbNSFP_gnomAD_exomes_AF": ["0.001", "0.002"],
                "dbNSFP_gnomAD_genomes_AF": ["0.003", "0.004"],
                "dbNSFP_CADD_phred": ["15.5", "20.1"],
                "EFFECT": ["missense_variant", "synonymous_variant"],
                "IMPACT": ["MODERATE", "LOW"],
            }
        )

    @pytest.fixture
    def scoring_config_dir(self):
        """Create a temporary scoring configuration directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Variable assignment config
            var_config = {
                "variables": {
                    "dbNSFP_gnomAD_exomes_AF": "gnomade_variant|default:0.0",
                    "dbNSFP_gnomAD_genomes_AF": "gnomadg_variant|default:0.0",
                    "dbNSFP_CADD_phred": "cadd_phred_variant|default:0.0",
                    "EFFECT": "consequence_terms_variant|default:''",
                    "IMPACT": "impact_variant|default:''",
                }
            }

            # Formula config - simplified version
            formula_config = {
                "formulas": [
                    {
                        "test_score": (
                            "gnomade_variant * 1000 + gnomadg_variant * 100 + cadd_phred_variant"
                        )
                    }
                ]
            }

            # Write configs
            with open(os.path.join(tmpdir, "variable_assignment_config.json"), "w") as f:
                json.dump(var_config, f)

            with open(os.path.join(tmpdir, "formula_config.json"), "w") as f:
                json.dump(formula_config, f)

            yield tmpdir

    def test_temporary_columns_removed(self, sample_df, scoring_config_dir):
        """Test that temporary variable columns are removed after scoring."""
        # Read scoring config
        scoring_config = read_scoring_config(scoring_config_dir)

        # Apply scoring
        result_df = apply_scoring(sample_df, scoring_config)

        # Check that original columns are preserved
        original_cols = list(sample_df.columns)
        for col in original_cols:
            assert col in result_df.columns, f"Original column {col} should be preserved"

        # Check that score column was added
        assert "test_score" in result_df.columns, "Score column should be added"

        # Check that temporary variable columns were removed
        temp_var_cols = [
            "gnomade_variant", "gnomadg_variant", "cadd_phred_variant",
            "consequence_terms_variant", "impact_variant"
        ]
        for col in temp_var_cols:
            assert col not in result_df.columns, \
                f"Temporary variable column {col} should be removed"

        # Verify the score was calculated correctly
        # Score = gnomade_variant * 1000 + gnomadg_variant * 100 + cadd_phred_variant
        # Row 1: 0.001 * 1000 + 0.003 * 100 + 15.5 = 1 + 0.3 + 15.5 = 16.8
        # Row 2: 0.002 * 1000 + 0.004 * 100 + 20.1 = 2 + 0.4 + 20.1 = 22.5
        assert abs(result_df["test_score"].iloc[0] - 16.8) < 0.01
        assert abs(result_df["test_score"].iloc[1] - 22.5) < 0.01

    def test_scoring_preserves_column_order(self, sample_df, scoring_config_dir):
        """Test that original column order is preserved."""
        # Read scoring config
        scoring_config = read_scoring_config(scoring_config_dir)

        # Apply scoring
        result_df = apply_scoring(sample_df, scoring_config)

        # Check that original columns appear in the same order
        original_cols = list(sample_df.columns)
        result_cols = list(result_df.columns)

        # The first n columns should match the original
        for i, col in enumerate(original_cols):
            assert result_cols[i] == col, f"Column order mismatch at position {i}"

        # Score columns should be at the end
        assert result_cols[-1] == "test_score"

    def test_missing_columns_with_defaults(self):
        """Test that missing columns are handled with defaults but not kept in output."""
        # DataFrame missing some columns
        df = pd.DataFrame(
            {
                "CHROM": ["1"],
                "POS": [1000],
                "REF": ["A"],
                "ALT": ["T"],
                "GENE": ["GENE1"],
                # Missing: dbNSFP_gnomAD_exomes_AF, dbNSFP_gnomAD_genomes_AF, etc.
            }
        )

        scoring_config = {
            "variables": {
                "dbNSFP_gnomAD_exomes_AF": "gnomade_variant|default:0.0",
                "dbNSFP_gnomAD_genomes_AF": "gnomadg_variant|default:0.0",
            },
            "formulas": [{"simple_score": "gnomade_variant + gnomadg_variant"}],
        }

        result_df = apply_scoring(df, scoring_config)

        # Original columns should be preserved
        for col in df.columns:
            assert col in result_df.columns

        # Score should be calculated with defaults (0.0 + 0.0 = 0.0)
        assert "simple_score" in result_df.columns
        assert result_df["simple_score"].iloc[0] == 0.0

        # Temporary columns should not be in output
        assert "gnomade_variant" not in result_df.columns
        assert "gnomadg_variant" not in result_df.columns
