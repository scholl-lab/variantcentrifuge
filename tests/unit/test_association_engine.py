"""
Unit tests for AssociationEngine orchestration.

Tests from_names() construction, error handling, run_all() gene sort order,
output schema, zero-variant gene skipping, and correction application.
"""

from __future__ import annotations

import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.engine import AssociationEngine


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def default_config():
    """Standard AssociationConfig for tests."""
    return AssociationConfig(
        correction_method="fdr",
        gene_burden_mode="samples",
    )


def _make_gene_data(
    gene: str,
    p_carriers: int,
    c_carriers: int,
    p_total: int,
    c_total: int,
    n_variants: int,
    p_alleles: int = 0,
    c_alleles: int = 0,
) -> dict:
    """Build a gene burden data dict matching _aggregate_gene_burden_from_* output."""
    return {
        "GENE": gene,
        "proband_count": p_total,
        "control_count": c_total,
        "proband_carrier_count": p_carriers,
        "control_carrier_count": c_carriers,
        "proband_allele_count": p_alleles or p_carriers,
        "control_allele_count": c_alleles or c_carriers,
        "n_qualifying_variants": n_variants,
    }


# ---------------------------------------------------------------------------
# from_names() construction tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestAssociationEngineConstruction:
    """Tests for AssociationEngine.from_names() factory."""

    def test_from_names_fisher_succeeds(self, default_config):
        """from_names(['fisher'], config) creates engine with Fisher test."""
        engine = AssociationEngine.from_names(["fisher"], default_config)
        assert engine is not None
        assert "fisher" in engine._tests

    def test_from_names_unknown_test_raises_value_error(self, default_config):
        """from_names(['skat'], config) raises ValueError mentioning available tests."""
        with pytest.raises(ValueError, match="Available tests: fisher"):
            AssociationEngine.from_names(["skat"], default_config)

    def test_from_names_unknown_test_error_includes_test_name(self, default_config):
        """ValueError message includes the unknown test name."""
        with pytest.raises(ValueError, match="skat"):
            AssociationEngine.from_names(["skat"], default_config)

    def test_from_names_empty_list_creates_engine_no_tests(self, default_config):
        """from_names([], config) creates engine with empty test set."""
        engine = AssociationEngine.from_names([], default_config)
        assert engine is not None
        assert len(engine._tests) == 0

    def test_from_names_fisher_check_dependencies_passes(self, default_config):
        """Fisher test check_dependencies() does not raise (scipy/statsmodels installed)."""
        # If this raises, scipy or statsmodels is missing from the test environment
        engine = AssociationEngine.from_names(["fisher"], default_config)
        assert engine is not None


# ---------------------------------------------------------------------------
# run_all() output schema tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestAssociationEngineRunAll:
    """Tests for AssociationEngine.run_all() output."""

    def test_run_all_single_gene_returns_dataframe(self, default_config):
        """run_all with one gene returns a DataFrame with expected columns."""
        import pandas as pd

        gene_data = [
            _make_gene_data("BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3)
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1

    def test_run_all_output_columns(self, default_config):
        """run_all output contains all expected columns."""
        gene_data = [
            _make_gene_data("BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3)
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        expected_cols = {
            "gene",
            "n_cases",
            "n_controls",
            "n_variants",
            "fisher_p_value",
            "fisher_corrected_p_value",
            "fisher_or",
            "fisher_or_ci_lower",
            "fisher_or_ci_upper",
        }
        assert expected_cols.issubset(set(result.columns))

    def test_run_all_multiple_genes_returns_sorted(self, default_config):
        """run_all with multiple genes returns them sorted alphabetically by gene."""
        gene_data = [
            _make_gene_data("TP53", p_carriers=3, c_carriers=0, p_total=50, c_total=50, n_variants=2),
            _make_gene_data("BRCA1", p_carriers=5, c_carriers=1, p_total=50, c_total=50, n_variants=3),
            _make_gene_data("MYH7", p_carriers=2, c_carriers=1, p_total=50, c_total=50, n_variants=1),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert list(result["gene"]) == ["BRCA1", "MYH7", "TP53"]

    def test_run_all_zero_variant_gene_excluded(self, default_config):
        """Gene with n_qualifying_variants=0 is excluded from results (p_value=None)."""
        gene_data = [
            _make_gene_data("BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3),
            # Zero-variant gene — should be excluded
            _make_gene_data("EMPTY", p_carriers=0, c_carriers=0, p_total=100, c_total=200, n_variants=0),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert "EMPTY" not in result["gene"].values
        assert len(result) == 1

    def test_run_all_all_zero_variant_genes_returns_empty_df(self, default_config):
        """run_all with all zero-variant genes returns empty DataFrame."""
        gene_data = [
            _make_gene_data("GENE1", p_carriers=0, c_carriers=0, p_total=50, c_total=50, n_variants=0),
            _make_gene_data("GENE2", p_carriers=0, c_carriers=0, p_total=50, c_total=50, n_variants=0),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert result.empty

    def test_run_all_empty_input_returns_empty_df(self, default_config):
        """run_all with empty gene_burden_data returns empty DataFrame."""
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all([])
        assert result.empty

    def test_run_all_corrected_p_values_populated(self, default_config):
        """run_all populates corrected p-values for all genes."""
        gene_data = [
            _make_gene_data("BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3),
            _make_gene_data("TP53", p_carriers=3, c_carriers=0, p_total=100, c_total=200, n_variants=2),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        # All corrected p-values must be non-None (correction was applied)
        assert result["fisher_corrected_p_value"].notna().all()

    def test_run_all_sort_order_reverse_input_still_alphabetical(self, default_config):
        """Input in reverse alpha order -> output still alphabetical."""
        # GENE names in reverse order to verify sort is applied
        gene_data = [
            _make_gene_data("ZZGENE", p_carriers=3, c_carriers=0, p_total=50, c_total=50, n_variants=1),
            _make_gene_data("MMGENE", p_carriers=2, c_carriers=1, p_total=50, c_total=50, n_variants=1),
            _make_gene_data("AAGENE", p_carriers=5, c_carriers=1, p_total=50, c_total=50, n_variants=1),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert list(result["gene"]) == ["AAGENE", "MMGENE", "ZZGENE"]

    def test_run_all_n_cases_n_controls_correct(self, default_config):
        """run_all output has correct n_cases and n_controls from contingency data."""
        gene_data = [
            _make_gene_data("BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3)
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        row = result.iloc[0]
        assert row["n_cases"] == 100
        assert row["n_controls"] == 200
        assert row["n_variants"] == 3

    def test_run_all_no_tests_returns_empty_df(self, default_config):
        """Engine with no registered tests returns empty DataFrame."""
        gene_data = [
            _make_gene_data("BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3)
        ]
        engine = AssociationEngine.from_names([], default_config)
        result = engine.run_all(gene_data)

        assert result.empty
