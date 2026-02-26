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
        """from_names(['nonexistent'], config) raises ValueError mentioning available tests."""
        with pytest.raises(ValueError, match="Available tests:"):
            AssociationEngine.from_names(["nonexistent_test_xyz"], default_config)

    def test_from_names_unknown_test_error_includes_test_name(self, default_config):
        """ValueError message includes the unknown test name."""
        with pytest.raises(ValueError, match="nonexistent_test_xyz"):
            AssociationEngine.from_names(["nonexistent_test_xyz"], default_config)

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
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            )
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1

    def test_run_all_output_columns(self, default_config):
        """run_all output contains all expected columns.

        ARCH-03: Primary tests no longer have corrected_p_value columns.
        FDR is applied only to ACAT-O, which appears as acat_o_pvalue and
        acat_o_qvalue.
        """
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            )
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        expected_cols = {
            "gene",
            "n_cases",
            "n_controls",
            "n_variants",
            "fisher_pvalue",
            # ARCH-03: no fisher_corrected_pvalue — FDR applied to ACAT-O only
            "acat_o_pvalue",
            "acat_o_qvalue",
            "fisher_or",
            "fisher_or_ci_lower",
            "fisher_or_ci_upper",
        }
        assert expected_cols.issubset(set(result.columns))
        # Primary test must NOT have a corrected_p_value column (ARCH-03)
        assert "fisher_corrected_pvalue" not in result.columns

    def test_run_all_multiple_genes_returns_sorted(self, default_config):
        """run_all with multiple genes returns them sorted alphabetically by gene."""
        gene_data = [
            _make_gene_data(
                "TP53", p_carriers=3, c_carriers=0, p_total=50, c_total=50, n_variants=2
            ),
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=50, c_total=50, n_variants=3
            ),
            _make_gene_data(
                "MYH7", p_carriers=2, c_carriers=1, p_total=50, c_total=50, n_variants=1
            ),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert list(result["gene"]) == ["BRCA1", "MYH7", "TP53"]

    def test_run_all_zero_variant_gene_excluded(self, default_config):
        """Gene with n_qualifying_variants=0 is excluded from results (p_value=None)."""
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            ),
            # Zero-variant gene — should be excluded
            _make_gene_data(
                "EMPTY", p_carriers=0, c_carriers=0, p_total=100, c_total=200, n_variants=0
            ),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert "EMPTY" not in result["gene"].values
        assert len(result) == 1

    def test_run_all_all_zero_variant_genes_returns_empty_df(self, default_config):
        """run_all with all zero-variant genes returns empty DataFrame."""
        gene_data = [
            _make_gene_data(
                "GENE1", p_carriers=0, c_carriers=0, p_total=50, c_total=50, n_variants=0
            ),
            _make_gene_data(
                "GENE2", p_carriers=0, c_carriers=0, p_total=50, c_total=50, n_variants=0
            ),
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
        """run_all populates ACAT-O corrected p-values for all genes.

        ARCH-03: FDR is applied only to ACAT-O p-values (not per-test).
        Primary test p-values (fisher_pvalue) are raw/uncorrected.
        acat_o_qvalue is populated after correction.
        """
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            ),
            _make_gene_data(
                "TP53", p_carriers=3, c_carriers=0, p_total=100, c_total=200, n_variants=2
            ),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        # ACAT-O corrected p-values must be non-None (correction was applied)
        assert result["acat_o_qvalue"].notna().all()
        # Primary test p-values are raw (no per-test FDR per ARCH-03)
        assert "fisher_corrected_pvalue" not in result.columns
        # Fisher raw p-values are present
        assert result["fisher_pvalue"].notna().all()

    def test_run_all_sort_order_reverse_input_still_alphabetical(self, default_config):
        """Input in reverse alpha order -> output still alphabetical."""
        # GENE names in reverse order to verify sort is applied
        gene_data = [
            _make_gene_data(
                "ZZGENE", p_carriers=3, c_carriers=0, p_total=50, c_total=50, n_variants=1
            ),
            _make_gene_data(
                "MMGENE", p_carriers=2, c_carriers=1, p_total=50, c_total=50, n_variants=1
            ),
            _make_gene_data(
                "AAGENE", p_carriers=5, c_carriers=1, p_total=50, c_total=50, n_variants=1
            ),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert list(result["gene"]) == ["AAGENE", "MMGENE", "ZZGENE"]

    def test_run_all_n_cases_n_controls_correct(self, default_config):
        """run_all output has correct n_cases and n_controls from contingency data."""
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            )
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
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            )
        ]
        engine = AssociationEngine.from_names([], default_config)
        result = engine.run_all(gene_data)

        assert result.empty


# ---------------------------------------------------------------------------
# Test-aware column naming
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestAssociationEngineColumnNaming:
    """Tests that column naming is test-aware (OR for Fisher, beta/SE for burden)."""

    def test_fisher_uses_or_columns(self, default_config):
        """Fisher test produces _or, _or_ci_lower, _or_ci_upper columns."""
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            )
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert "fisher_or" in result.columns
        assert "fisher_or_ci_lower" in result.columns
        assert "fisher_or_ci_upper" in result.columns
        # Fisher should NOT have _se column
        assert "fisher_se" not in result.columns

    def test_burden_uses_beta_se_columns(self, default_config):
        """Burden tests produce _beta, _se, _beta_ci_lower, _beta_ci_upper columns."""
        import numpy as np

        # Need genotype matrix for burden tests
        np.random.seed(42)
        n = 60
        geno = np.zeros((n, 3))
        for v in range(3):
            carriers = np.random.choice(n, size=int(n * 0.1), replace=False)
            geno[carriers, v] = 1
        phenotype = np.zeros(n)
        phenotype[:30] = 1
        mafs = geno.mean(axis=0) / 2

        gene_data = [
            {
                "GENE": "TEST_GENE",
                "proband_count": 30,
                "control_count": 30,
                "proband_carrier_count": 5,
                "control_carrier_count": 3,
                "proband_allele_count": 5,
                "control_allele_count": 3,
                "n_qualifying_variants": 3,
                "genotype_matrix": geno,
                "variant_mafs": mafs,
                "phenotype_vector": phenotype,
                "covariate_matrix": None,
            }
        ]

        engine = AssociationEngine.from_names(["logistic_burden"], default_config)
        result = engine.run_all(gene_data)

        assert "logistic_burden_beta" in result.columns
        assert "logistic_burden_se" in result.columns
        assert "logistic_burden_beta_ci_lower" in result.columns
        assert "logistic_burden_beta_ci_upper" in result.columns
        # Should NOT have _or columns
        assert "logistic_burden_or" not in result.columns


# ---------------------------------------------------------------------------
# ACAT-O integration tests (Phase 22 / ARCH-03)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestAssociationEngineAcatO:
    """Tests for ACAT-O omnibus column presence and FDR strategy (ARCH-03).

    These tests verify that:
    1. acat_o_pvalue and acat_o_qvalue columns are always present
    2. Primary test corrected_p_value columns are never present (ARCH-03)
    3. Single-test pass-through: acat_o_pvalue == fisher_pvalue when k=1
    4. FDR is applied only to ACAT-O (corrected p-values populated for valid genes)
    5. Zero-variant genes have acat_o_pvalue=None
    """

    def test_engine_acat_o_columns_present(self, default_config):
        """run_all with fisher test produces acat_o_pvalue and acat_o_qvalue columns."""
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            )
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        assert "acat_o_pvalue" in result.columns, "acat_o_pvalue column missing from output"
        assert "acat_o_qvalue" in result.columns, "acat_o_qvalue column missing from output"

    def test_engine_acat_o_single_test_passthrough(self, default_config):
        """With only fisher test, acat_o_pvalue equals fisher_pvalue (k=1 pass-through).

        IMPL-35: Single valid p-value passes through cauchy_combination unchanged.
        When there is only one primary test, ACAT-O is identical to that test's p-value.
        """
        import pytest

        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            )
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        row = result.iloc[0]
        assert row["acat_o_pvalue"] is not None
        assert row["fisher_pvalue"] is not None
        assert row["acat_o_pvalue"] == pytest.approx(row["fisher_pvalue"], rel=1e-9), (
            "With k=1 test, acat_o_pvalue should equal fisher_pvalue (pass-through)"
        )

    def test_engine_acat_o_fdr_applied(self, default_config):
        """acat_o_qvalue is not None for genes with valid ACAT-O p-values.

        ARCH-03: FDR is applied to ACAT-O across all genes. After correction,
        acat_o_qvalue should be populated for every gene that had a
        valid acat_o_pvalue.
        """
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            ),
            _make_gene_data(
                "TP53", p_carriers=3, c_carriers=0, p_total=100, c_total=200, n_variants=2
            ),
            _make_gene_data(
                "MYH7", p_carriers=8, c_carriers=2, p_total=100, c_total=200, n_variants=5
            ),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        # All genes that have acat_o_pvalue should also have acat_o_qvalue
        genes_with_acat_p = result[result["acat_o_pvalue"].notna()]
        assert len(genes_with_acat_p) > 0, "No genes had valid acat_o_pvalue"
        assert genes_with_acat_p["acat_o_qvalue"].notna().all(), (
            "Some genes with valid acat_o_pvalue have None acat_o_qvalue"
        )

    def test_engine_primary_tests_no_fdr(self, default_config):
        """Primary test corrected_p_value columns are absent (ARCH-03).

        FDR is applied only to ACAT-O. Primary test raw p-values exist for
        diagnostic signal decomposition, but they are NOT independently corrected.
        """
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            ),
            _make_gene_data(
                "TP53", p_carriers=3, c_carriers=0, p_total=100, c_total=200, n_variants=2
            ),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        # Fisher must not have a corrected p-value column (ARCH-03)
        assert "fisher_corrected_pvalue" not in result.columns, (
            "ARCH-03 violated: fisher_corrected_pvalue should not be in output columns"
        )

    def test_engine_acat_o_none_for_zero_variant_gene(self, default_config):
        """Gene with 0 variants is excluded from results; no row means no acat_o_pvalue.

        IMPL-02: p_value=None for zero-variant genes (not 1.0). These genes are
        skipped in the engine output entirely (not included as rows with None values).
        """
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            ),
            _make_gene_data(
                "EMPTY", p_carriers=0, c_carriers=0, p_total=100, c_total=200, n_variants=0
            ),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        # EMPTY gene (zero variants) is excluded entirely
        assert "EMPTY" not in result["gene"].values, (
            "Zero-variant gene EMPTY should be excluded from output"
        )
        # Only BRCA1 should be present
        assert len(result) == 1
        assert result.iloc[0]["gene"] == "BRCA1"
        # BRCA1 has valid acat_o_pvalue
        assert result.iloc[0]["acat_o_pvalue"] is not None

    def test_engine_acat_o_corrected_p_in_range(self, default_config):
        """acat_o_qvalue values are in [0, 1]."""
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            ),
            _make_gene_data(
                "TP53", p_carriers=3, c_carriers=0, p_total=100, c_total=200, n_variants=2
            ),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        for val in result["acat_o_qvalue"].dropna():
            assert 0.0 <= val <= 1.0, f"acat_o_qvalue {val} outside [0, 1]"

    def test_engine_acat_o_raw_p_in_range(self, default_config):
        """acat_o_pvalue raw p-values are in [0, 1]."""
        gene_data = [
            _make_gene_data(
                "BRCA1", p_carriers=5, c_carriers=1, p_total=100, c_total=200, n_variants=3
            ),
        ]
        engine = AssociationEngine.from_names(["fisher"], default_config)
        result = engine.run_all(gene_data)

        for val in result["acat_o_pvalue"].dropna():
            assert 0.0 <= val <= 1.0, f"acat_o_pvalue {val} outside [0, 1]"
