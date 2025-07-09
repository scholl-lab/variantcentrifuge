"""Integration tests for custom statistics with the existing stats module."""

import json
import os
import tempfile

import pandas as pd
import pytest

from variantcentrifuge import stats


@pytest.fixture
def sample_df():
    """Create a sample DataFrame for testing."""
    data = {
        "GENE": ["BRCA1", "BRCA1", "BRCA2", "BRCA2", "TP53"],
        "IMPACT": ["HIGH", "MODERATE", "HIGH", "LOW", "HIGH"],
        "EFFECT": [
            "stop_gained",
            "missense_variant",
            "frameshift_variant",
            "synonymous_variant",
            "stop_gained",
        ],
        "dbNSFP_CADD_phred": [35.0, 22.5, 40.0, 5.2, 38.0],
        "gnomAD_exomes_AF": [0.0001, 0.01, 0.0002, 0.1, 0.00005],
        "GT": [
            "SAMPLE1(0/1)",
            "SAMPLE1(1/1);SAMPLE2(0/1)",
            "SAMPLE2(1/1)",
            "SAMPLE1(0/1)",
            "SAMPLE1(1/1);SAMPLE2(0/1)",
        ],
        "proband_count": [1, 2, 1, 1, 2],
        "control_count": [0, 1, 1, 0, 1],
        "proband_allele_count": [1, 3, 2, 1, 3],
        "control_allele_count": [0, 1, 2, 0, 1],
    }
    return pd.DataFrame(data)


@pytest.fixture
def custom_stats_config():
    """Create custom statistics configuration for testing."""
    config = {
        "stats_version": "1.0",
        "dataset_stats": [
            {
                "name": "ultra_rare_count",
                "expression": "(df['gnomAD_exomes_AF'] < 0.0001).sum()",
                "required_columns": ["gnomAD_exomes_AF"],
            }
        ],
        "gene_stats": [
            {
                "name": "max_cadd_score",
                "expression": "group_df['dbNSFP_CADD_phred'].max()",
                "groupby": "GENE",
                "required_columns": ["dbNSFP_CADD_phred"],
            },
            {
                "name": "high_impact_fraction",
                "expression": "(group_df['IMPACT'] == 'HIGH').mean()",
                "groupby": "GENE",
                "required_columns": ["IMPACT"],
            },
        ],
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(config, f)
        return f.name


class TestStatsIntegration:
    """Test integration of custom stats with existing stats module."""

    def test_compute_custom_stats_with_config(self, sample_df, custom_stats_config):
        """Test compute_custom_stats with a configuration file."""
        try:
            results = stats.compute_custom_stats(sample_df, custom_stats_config)

            assert "dataset" in results
            assert "genes" in results

            # Check dataset stats
            dataset_stats = results["dataset"]
            ultra_rare = dataset_stats[dataset_stats["metric"] == "ultra_rare_count"]["value"].iloc[
                0
            ]
            assert ultra_rare == 1  # Only TP53 has AF < 0.0001 (BRCA1 has AF = 0.0001)

            # Check gene stats
            gene_stats = results["genes"]
            assert set(gene_stats["GENE"].values) == {"BRCA1", "BRCA2", "TP53"}

            # BRCA1 should have max CADD of 35.0
            brca1_stats = gene_stats[gene_stats["GENE"] == "BRCA1"]
            assert brca1_stats["max_cadd_score"].iloc[0] == 35.0
            assert brca1_stats["high_impact_fraction"].iloc[0] == 0.5  # 1 HIGH out of 2

        finally:
            os.unlink(custom_stats_config)

    def test_compute_custom_stats_with_default(self, sample_df):
        """Test compute_custom_stats with default configuration."""
        # This should use the default_stats_config.json
        results = stats.compute_custom_stats(sample_df)

        assert isinstance(results, dict)
        # Default config should compute basic stats
        if "dataset" in results and not results["dataset"].empty:
            assert "total_variants" in results["dataset"]["metric"].values

    def test_merge_with_custom_stats(self, sample_df, custom_stats_config):
        """Test merging custom stats with traditional stats."""
        try:
            # Compute traditional stats
            gene_stats = stats.compute_gene_stats(sample_df)
            impact_summary = stats.compute_impact_summary(sample_df)
            variant_type_summary = stats.compute_variant_type_summary(sample_df)

            # Compute custom stats
            custom_stats = stats.compute_custom_stats(sample_df, custom_stats_config)

            # Merge all stats
            merged = stats.merge_and_format_stats(
                gene_stats, impact_summary, variant_type_summary, custom_stats
            )

            # Check that all stats are present
            assert "GENE" in merged.columns
            assert "proband_count" in merged.columns  # Traditional stat
            assert "HIGH" in merged.columns  # Impact summary
            assert "max_cadd_score" in merged.columns  # Custom stat
            assert "high_impact_fraction" in merged.columns  # Custom stat

            # Verify values
            brca1_row = merged[merged["GENE"] == "BRCA1"].iloc[0]
            assert brca1_row["proband_count"] == 3  # 1 + 2
            assert brca1_row["HIGH"] == 1  # One HIGH impact variant
            assert brca1_row["max_cadd_score"] == 35.0
            assert brca1_row["high_impact_fraction"] == 0.5

        finally:
            os.unlink(custom_stats_config)

    def test_merge_without_custom_stats(self, sample_df):
        """Test that merge still works without custom stats (backward compatibility)."""
        # Compute traditional stats
        gene_stats = stats.compute_gene_stats(sample_df)
        impact_summary = stats.compute_impact_summary(sample_df)
        variant_type_summary = stats.compute_variant_type_summary(sample_df)

        # Merge without custom stats
        merged = stats.merge_and_format_stats(gene_stats, impact_summary, variant_type_summary)

        # Should work as before
        assert "GENE" in merged.columns
        assert "proband_count" in merged.columns
        assert "HIGH" in merged.columns
        assert "max_cadd_score" not in merged.columns  # Custom stat not present

    def test_error_handling(self, sample_df):
        """Test error handling for invalid config path."""
        results = stats.compute_custom_stats(sample_df, "/nonexistent/path.json")
        assert results == {}  # Should return empty dict on error

    def test_custom_stats_with_missing_columns(self, sample_df):
        """Test custom stats when DataFrame is missing expected columns."""
        # Create config that expects columns not in our DataFrame
        config = {
            "gene_stats": [
                {
                    "name": "missing_column_stat",
                    "expression": "group_df['NONEXISTENT'].mean()",
                    "groupby": "GENE",
                    "required_columns": ["NONEXISTENT"],
                }
            ]
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config, f)
            config_path = f.name

        try:
            results = stats.compute_custom_stats(sample_df, config_path)
            # Should handle gracefully - stat should be skipped
            assert "genes" in results
            assert results["genes"].empty or "missing_column_stat" not in results["genes"].columns
        finally:
            os.unlink(config_path)
