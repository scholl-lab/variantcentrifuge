# File: test_gene_burden.py
# Location: variantcentrifuge/tests/test_gene_burden.py
"""
Test suite for gene burden analysis with focus on edge cases.

Tests cover:
- Zero cells in contingency tables
- Structural zeros (entire row/column zero)
- Infinite and zero odds ratios
- Confidence interval calculation methods
- Integration with the full pipeline
"""

import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch

from variantcentrifuge.gene_burden import (
    perform_gene_burden_analysis,
    _compute_or_confidence_interval,
)


class TestOddsRatioEdgeCases:
    """Test odds ratio calculation with various edge cases through gene burden analysis."""

    @pytest.fixture
    def sample_config(self):
        """Provide basic configuration for gene burden analysis."""
        return {
            "gene_burden_mode": "samples",
            "correction_method": "fdr",
            "confidence_interval_method": "normal_approx",
            "confidence_interval_alpha": 0.05,
            "continuity_correction": 0.5,
        }


class TestGeneBurdenEdgeCases:
    """Test gene burden analysis with edge case scenarios."""

    @pytest.fixture
    def edge_case_data(self):
        """Create test data with various edge cases."""
        data = [
            # Gene with zero variants in controls (infinite OR)
            {
                "GENE": "GENE1",
                "proband_count": 100,
                "control_count": 100,
                "proband_variant_count": 10,
                "control_variant_count": 0,
                "proband_allele_count": 15,
                "control_allele_count": 0,
            },
            # Gene with zero variants in cases (zero OR)
            {
                "GENE": "GENE2",
                "proband_count": 100,
                "control_count": 100,
                "proband_variant_count": 0,
                "control_variant_count": 10,
                "proband_allele_count": 0,
                "control_allele_count": 15,
            },
            # Gene with no variants at all (structural zero)
            {
                "GENE": "GENE3",
                "proband_count": 100,
                "control_count": 100,
                "proband_variant_count": 0,
                "control_variant_count": 0,
                "proband_allele_count": 0,
                "control_allele_count": 0,
            },
            # Gene with only one sample total
            {
                "GENE": "GENE4",
                "proband_count": 1,
                "control_count": 0,
                "proband_variant_count": 1,
                "control_variant_count": 0,
                "proband_allele_count": 2,
                "control_allele_count": 0,
            },
            # Normal gene for comparison
            {
                "GENE": "GENE5",
                "proband_count": 100,
                "control_count": 100,
                "proband_variant_count": 10,
                "control_variant_count": 5,
                "proband_allele_count": 15,
                "control_allele_count": 7,
            },
        ]
        return pd.DataFrame(data)

    def test_original_implementation_edge_cases(self, edge_case_data):
        """Test current implementation with edge cases."""
        config = {
            "gene_burden_mode": "samples",
            "correction_method": "fdr",
            "confidence_interval_method": "normal_approx",
            "confidence_interval_alpha": 0.05,
        }

        result = perform_gene_burden_analysis(edge_case_data, config)

        # Check that all genes are in results
        assert len(result) == 5

        # Check GENE1 (infinite OR case)
        gene1 = result[result["GENE"] == "GENE1"].iloc[0]
        assert not np.isnan(gene1["odds_ratio"])  # Should not be NaN
        assert gene1["or_ci_lower"] > 0  # Should have valid CI
        assert gene1["or_ci_upper"] > gene1["or_ci_lower"]

        # Check GENE2 (zero OR case)
        gene2 = result[result["GENE"] == "GENE2"].iloc[0]
        assert gene2["odds_ratio"] == 0  # Fisher exact returns 0
        assert gene2["or_ci_lower"] >= 0
        assert gene2["or_ci_upper"] > 0

        # Check GENE3 (no variants)
        gene3 = result[result["GENE"] == "GENE3"].iloc[0]
        assert np.isnan(gene3["odds_ratio"])  # Structural zero

        # Check GENE4 (no controls - structural zero)
        gene4 = result[result["GENE"] == "GENE4"].iloc[0]
        assert np.isnan(gene4["odds_ratio"])  # Structural zero
        assert np.isnan(gene4["or_ci_lower"])
        assert np.isnan(gene4["or_ci_upper"])

    def test_allele_mode_edge_cases(self, edge_case_data):
        """Test edge cases in allele counting mode."""
        config = {
            "gene_burden_mode": "alleles",
            "correction_method": "bonferroni",
            "confidence_interval_method": "normal_approx",
            "confidence_interval_alpha": 0.05,
        }

        result = perform_gene_burden_analysis(edge_case_data, config)

        # Should produce valid results for genes with samples
        assert len(result) >= 3
        assert all(
            col in result.columns
            for col in [
                "proband_allele_count",
                "control_allele_count",
                "odds_ratio",
                "or_ci_lower",
                "or_ci_upper",
            ]
        )


class TestConfidenceIntervalMethods:
    """Test different confidence interval calculation methods."""

    def test_compute_or_ci_fallback_behavior(self):
        """Test the fallback behavior of _compute_or_confidence_interval."""
        # Test with invalid odds ratio (NaN) - structural zero
        table = [[0, 0], [0, 0]]
        ci_lower, ci_upper = _compute_or_confidence_interval(
            table, float("nan"), "normal_approx", 0.05
        )

        # Should return NaN for structural zeros
        assert np.isnan(ci_lower)
        assert np.isnan(ci_upper)

    def test_compute_or_ci_with_zero_or(self):
        """Test CI computation when OR is zero."""
        table = [[0, 10], [100, 90]]
        ci_lower, ci_upper = _compute_or_confidence_interval(table, 0.0, "normal_approx", 0.05)

        # With continuity correction, should return valid CI
        assert ci_lower >= 0
        assert ci_upper > ci_lower
        assert not np.isinf(ci_upper)

    def test_compute_or_ci_with_infinite_or(self):
        """Test CI computation when OR is infinite."""
        table = [[10, 0], [90, 100]]
        ci_lower, ci_upper = _compute_or_confidence_interval(
            table, float("inf"), "normal_approx", 0.05
        )

        # With continuity correction, should return valid CI
        assert ci_lower > 0
        assert not np.isnan(ci_lower)
        assert not np.isnan(ci_upper)

    @patch("variantcentrifuge.gene_burden.Table2x2", None)
    def test_compute_or_ci_without_statsmodels(self):
        """Test fallback when statsmodels is not available."""
        table = [[10, 5], [20, 30]]
        ci_lower, ci_upper = _compute_or_confidence_interval(table, 3.0, "normal_approx", 0.05)

        # Should return NaN when statsmodels not available
        assert np.isnan(ci_lower)
        assert np.isnan(ci_upper)


class TestIntegrationEdgeCases:
    """Test edge cases in the context of the full pipeline."""

    @pytest.fixture
    def mock_variant_data(self):
        """Create mock variant data that will trigger edge cases."""
        # Multiple variants per gene to test aggregation
        data = []

        # GENE1: Multiple variants, but no controls have any
        for i in range(5):
            data.append(
                {
                    "GENE": "GENE1",
                    "proband_count": 50,
                    "control_count": 50,
                    "proband_variant_count": 1,
                    "control_variant_count": 0,
                    "proband_allele_count": 1,
                    "control_allele_count": 0,
                }
            )

        # GENE2: Sparse data with some zeros
        data.extend(
            [
                {
                    "GENE": "GENE2",
                    "proband_count": 10,
                    "control_count": 10,
                    "proband_variant_count": 1,
                    "control_variant_count": 0,
                    "proband_allele_count": 1,
                    "control_allele_count": 0,
                },
                {
                    "GENE": "GENE2",
                    "proband_count": 10,
                    "control_count": 10,
                    "proband_variant_count": 0,
                    "control_variant_count": 1,
                    "proband_allele_count": 0,
                    "control_allele_count": 2,
                },
            ]
        )

        return pd.DataFrame(data)

    def test_aggregation_with_edge_cases(self, mock_variant_data):
        """Test that aggregation handles edge cases correctly."""
        config = {"gene_burden_mode": "samples", "correction_method": "fdr"}

        result = perform_gene_burden_analysis(mock_variant_data, config)

        # Check GENE1 aggregation (should sum to 5 variants in cases, 0 in controls)
        gene1 = result[result["GENE"] == "GENE1"].iloc[0]
        assert gene1["proband_variant_count"] == 5
        assert gene1["control_variant_count"] == 0

        # Check GENE2 aggregation
        gene2 = result[result["GENE"] == "GENE2"].iloc[0]
        assert gene2["proband_variant_count"] == 1
        assert gene2["control_variant_count"] == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
