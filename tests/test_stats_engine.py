"""Tests for the configurable statistics engine."""

import json
import os
import tempfile

import pandas as pd
import pytest

from variantcentrifuge.stats_engine import StatsEngine, _validate_expression


@pytest.fixture
def sample_df():
    """Create a sample DataFrame for testing."""
    data = {
        "GENE": ["BRCA1", "BRCA1", "BRCA2", "BRCA2", "TP53", "TP53"],
        "IMPACT": ["HIGH", "MODERATE", "HIGH", "LOW", "HIGH", "MODERATE"],
        "Consequence": [
            "stop_gained",
            "missense_variant",
            "frameshift_variant",
            "synonymous_variant",
            "stop_gained",
            "missense_variant",
        ],
        "dbNSFP_CADD_phred": [35.0, 22.5, 40.0, 5.2, 38.0, 25.0],
        "gnomAD_exomes_AF": [0.0001, 0.01, 0.0002, 0.1, 0.00005, 0.005],
        "SAMPLE1_GT": ["0/1", "1/1", "0/1", "0/0", "1/1", "0/1"],
        "SAMPLE2_GT": ["0/0", "0/1", "1/1", "0/1", "0/1", "0/0"],
    }
    return pd.DataFrame(data)


@pytest.fixture
def basic_config():
    """Create basic statistics configuration."""
    return {
        "stats_version": "1.0",
        "dataset_stats": [
            {"name": "total_variants", "expression": "len(df)"},
            {
                "name": "high_impact_ratio",
                "expression": "(df['IMPACT'] == 'HIGH').mean()",
                "required_columns": ["IMPACT"],
            },
        ],
        "gene_stats": [
            {"name": "variant_count", "expression": "size()", "groupby": "GENE"},
            {
                "name": "mean_cadd",
                "expression": "group_df['dbNSFP_CADD_phred'].mean()",
                "groupby": "GENE",
                "required_columns": ["dbNSFP_CADD_phred"],
            },
        ],
    }


class TestStatsEngine:
    """Test the StatsEngine class."""

    def test_init_with_dict(self, basic_config):
        """Test initialization with a dictionary configuration."""
        engine = StatsEngine(basic_config)
        assert engine.config == basic_config
        assert engine.results == {}

    def test_init_with_file(self, basic_config):
        """Test initialization with a JSON file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(basic_config, f)
            config_path = f.name

        try:
            engine = StatsEngine(config_path)
            assert engine.config == basic_config
        finally:
            os.unlink(config_path)

    def test_dataset_stats(self, sample_df, basic_config):
        """Test dataset-level statistics computation."""
        engine = StatsEngine(basic_config)
        results = engine.compute(sample_df)

        assert "dataset" in results
        dataset_stats = results["dataset"]

        # Check total_variants
        assert len(dataset_stats[dataset_stats["metric"] == "total_variants"]) == 1
        total_variants = dataset_stats[dataset_stats["metric"] == "total_variants"]["value"].iloc[0]
        assert total_variants == 6

        # Check high_impact_ratio
        assert len(dataset_stats[dataset_stats["metric"] == "high_impact_ratio"]) == 1
        high_impact_ratio = dataset_stats[dataset_stats["metric"] == "high_impact_ratio"][
            "value"
        ].iloc[0]
        assert abs(high_impact_ratio - 0.5) < 0.001  # 3 out of 6 are HIGH

    def test_gene_stats(self, sample_df, basic_config):
        """Test gene-level statistics computation."""
        engine = StatsEngine(basic_config)
        results = engine.compute(sample_df)

        assert "genes" in results
        gene_stats = results["genes"]

        # Check variant counts
        assert set(gene_stats["GENE"].values) == {"BRCA1", "BRCA2", "TP53"}
        assert gene_stats[gene_stats["GENE"] == "BRCA1"]["variant_count"].iloc[0] == 2
        assert gene_stats[gene_stats["GENE"] == "BRCA2"]["variant_count"].iloc[0] == 2
        assert gene_stats[gene_stats["GENE"] == "TP53"]["variant_count"].iloc[0] == 2

        # Check mean CADD scores
        brca1_mean_cadd = gene_stats[gene_stats["GENE"] == "BRCA1"]["mean_cadd"].iloc[0]
        assert abs(brca1_mean_cadd - 28.75) < 0.001  # (35.0 + 22.5) / 2

    def test_grouped_stats(self, sample_df):
        """Test custom grouped statistics."""
        config = {
            "grouped_stats": [
                {
                    "name": "impact_consequence_counts",
                    "expression": "size()",
                    "groupby": ["IMPACT", "Consequence"],
                    "output_format": "long",
                }
            ]
        }

        engine = StatsEngine(config)
        results = engine.compute(sample_df)

        assert "groups" in results
        assert "impact_consequence_counts" in results["groups"]

        grouped_df = results["groups"]["impact_consequence_counts"]
        assert len(grouped_df) == 4  # 4 unique IMPACT/Consequence combinations

    def test_pivot_output_format(self, sample_df):
        """Test pivot output format for grouped stats."""
        config = {
            "grouped_stats": [
                {
                    "name": "gene_impact_matrix",
                    "expression": "size()",
                    "groupby": ["GENE", "IMPACT"],
                    "output_format": "pivot",
                }
            ]
        }

        engine = StatsEngine(config)
        results = engine.compute(sample_df)

        matrix = results["groups"]["gene_impact_matrix"]
        assert matrix.loc["BRCA1", "HIGH"] == 1
        assert matrix.loc["BRCA1", "MODERATE"] == 1
        assert matrix.loc["BRCA2", "LOW"] == 1

    def test_missing_required_columns(self, sample_df):
        """Test handling of missing required columns."""
        config = {
            "dataset_stats": [
                {
                    "name": "missing_column_stat",
                    "expression": "df['NONEXISTENT'].mean()",
                    "required_columns": ["NONEXISTENT"],
                }
            ]
        }

        engine = StatsEngine(config)
        results = engine.compute(sample_df)

        # Should skip the stat with missing column
        assert "dataset" in results
        assert results["dataset"].empty

    def test_invalid_expression(self, sample_df):
        """Test handling of invalid expressions."""
        config = {
            "dataset_stats": [{"name": "invalid_stat", "expression": "this is not valid python"}]
        }

        engine = StatsEngine(config)
        results = engine.compute(sample_df)

        # Should handle error gracefully
        assert "dataset" in results
        assert results["dataset"].empty

    def test_complex_expressions(self, sample_df):
        """Test complex expressions with multiple operations."""
        config = {
            "dataset_stats": [
                {
                    "name": "rare_high_impact_count",
                    "expression": (
                        "((df['gnomAD_exomes_AF'] < 0.001) & (df['IMPACT'] == 'HIGH')).sum()"
                    ),
                    "required_columns": ["gnomAD_exomes_AF", "IMPACT"],
                }
            ],
            "gene_stats": [
                {
                    "name": "max_cadd_rare_variants",
                    "expression": (
                        "group_df[group_df['gnomAD_exomes_AF'] < 0.001]"
                        "['dbNSFP_CADD_phred'].max() "
                        "if len(group_df[group_df['gnomAD_exomes_AF'] < 0.001]) > 0 else 0"
                    ),
                    "groupby": "GENE",
                    "required_columns": ["gnomAD_exomes_AF", "dbNSFP_CADD_phred"],
                }
            ],
        }

        engine = StatsEngine(config)
        results = engine.compute(sample_df)

        # Check rare high impact count
        rare_high_count = results["dataset"][
            results["dataset"]["metric"] == "rare_high_impact_count"
        ]["value"].iloc[0]
        assert rare_high_count == 3  # BRCA1 HIGH, BRCA2 HIGH, TP53 HIGH all have AF < 0.001

        # Check max CADD for rare variants per gene
        gene_stats = results["genes"]
        brca1_max = gene_stats[gene_stats["GENE"] == "BRCA1"]["max_cadd_rare_variants"].iloc[0]
        assert brca1_max == 35.0  # Only the HIGH impact variant has AF < 0.001

    def test_formatted_output(self, sample_df, basic_config):
        """Test the formatted output method."""
        engine = StatsEngine(basic_config)
        engine.compute(sample_df)

        formatted = engine.get_formatted_output()
        assert "=== Dataset Statistics ===" in formatted
        assert "=== Gene Statistics ===" in formatted
        assert "total_variants: 6" in formatted
        assert "BRCA1" in formatted


class TestExpressionValidation:
    """Test AST-based expression validation blocks code injection."""

    def test_valid_expressions_pass(self):
        """Verify legitimate stats expressions are allowed."""
        valid_expressions = [
            "len(df)",
            "(df['IMPACT'] == 'HIGH').mean()",
            "group_df['dbNSFP_CADD_phred'].mean()",
            "pd.to_numeric(df['AF'], errors='coerce').mean()",
            "((df['AF'] < 0.001) & (df['IMPACT'] == 'HIGH')).sum()",
            "len([col for col in df.columns if col.endswith('_GT')])",
            "group_df['col'].str.contains('Pathogenic', case=False, na=False).sum()",
            "group_df['col'].nunique()",
            "df['col'].max() if len(df) > 0 else 0",
        ]
        for expr in valid_expressions:
            _validate_expression(expr)  # should not raise

    def test_import_blocked(self):
        """Verify __import__ is rejected (unknown name)."""
        with pytest.raises(ValueError, match="Unknown name"):
            _validate_expression("__import__('os').system('rm -rf /')")

    def test_dunder_attribute_blocked(self):
        """Verify dunder attribute access is rejected."""
        with pytest.raises(ValueError, match="private/dunder attribute"):
            _validate_expression("df.__class__.__bases__[0].__subclasses__()")

    def test_dunder_init_blocked(self):
        """Verify __init__ access is rejected."""
        with pytest.raises(ValueError, match="private/dunder attribute"):
            _validate_expression("().__class__.__init__.__globals__")

    def test_unknown_name_blocked(self):
        """Verify unknown variable names are rejected."""
        with pytest.raises(ValueError, match="Unknown name"):
            _validate_expression("os.system('whoami')")

    def test_open_blocked(self):
        """Verify open() is rejected (not in allowed names)."""
        with pytest.raises(ValueError, match="Unknown name"):
            _validate_expression("open('/etc/passwd').read()")

    def test_eval_blocked(self):
        """Verify nested eval is rejected."""
        with pytest.raises(ValueError, match="Unknown name"):
            _validate_expression("eval('__import__(\"os\")')")

    def test_exec_blocked(self):
        """Verify exec is rejected."""
        with pytest.raises(ValueError, match="Unknown name"):
            _validate_expression("exec('import os')")

    def test_injection_via_config_blocked(self, sample_df):
        """Verify malicious expressions in config are blocked at runtime."""
        config = {
            "dataset_stats": [
                {
                    "name": "malicious",
                    "expression": "().__class__.__bases__[0].__subclasses__()",
                }
            ]
        }
        engine = StatsEngine(config)
        results = engine.compute(sample_df)
        # Should be blocked, returning empty results
        assert results["dataset"].empty

    def test_syntax_error_rejected(self):
        """Verify syntax errors are caught during validation."""
        with pytest.raises(ValueError, match="Syntax error"):
            _validate_expression("this is not valid python")
