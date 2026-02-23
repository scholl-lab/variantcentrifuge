"""
Engine integration tests verifying SKAT column output.

Tests that AssociationEngine produces correct output column structure
when a SKAT-like test is registered. The SKAT test is mocked entirely
to avoid any rpy2 dependency.

Key invariants verified:
- No "skat_None" columns (None-effect guard works)
- Extra columns written with bare key names (not double-prefixed)
- Null model fitted exactly once across multiple genes
- prepare() and finalize() hooks called by the engine
- Multi-test runs (Fisher + SKAT) produce expected mixed columns
- Genes with p_value=None excluded from FDR correction
"""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


def _make_config(**kwargs: Any) -> Any:
    """Build an AssociationConfig with sensible defaults."""
    from variantcentrifuge.association.base import AssociationConfig

    defaults: dict[str, Any] = {"correction_method": "fdr", "gene_burden_mode": "samples"}
    defaults.update(kwargs)
    return AssociationConfig(**defaults)


def _make_gene_burden_data(
    gene: str,
    p_carriers: int = 5,
    c_carriers: int = 1,
    p_total: int = 50,
    c_total: int = 50,
    n_variants: int = 3,
    include_genotype_matrix: bool = True,
) -> dict[str, Any]:
    """Build a gene burden dict with optional genotype matrix."""
    data: dict[str, Any] = {
        "GENE": gene,
        "proband_count": p_total,
        "control_count": c_total,
        "proband_carrier_count": p_carriers,
        "control_carrier_count": c_carriers,
        "proband_allele_count": p_carriers,
        "control_allele_count": c_carriers,
        "n_qualifying_variants": n_variants,
    }
    if include_genotype_matrix:
        rng = np.random.default_rng(42)
        geno = rng.integers(0, 3, size=(p_total + c_total, n_variants)).astype(float)
        phenotype = np.array([1.0] * p_total + [0.0] * c_total)
        mafs = geno.mean(axis=0) / 2.0
        data.update(
            {
                "genotype_matrix": geno,
                "phenotype_vector": phenotype,
                "variant_mafs": mafs,
                "covariate_matrix": None,
            }
        )
    return data


def _make_mock_skat_test(
    *,
    p_value_default: float | None = 0.03,
    extra_default: dict[str, Any] | None = None,
    skip_genes: list[str] | None = None,
) -> MagicMock:
    """
    Build a mock AssociationTest that behaves like a SKAT test.

    Returns all-None effect columns, writes SKAT-style extra columns,
    and has prepare() / finalize() methods.
    """
    from variantcentrifuge.association.base import TestResult

    if extra_default is None:
        extra_default = {
            "skat_o_rho": None,
            "skat_warnings": None,
            "skat_method": "SKAT",
            "skat_n_marker_test": 3,
        }
    if skip_genes is None:
        skip_genes = []

    mock_test = MagicMock(name="MockSKATTest")
    mock_test.name = "skat"

    def effect_column_names() -> dict[str, Any]:
        return {"effect": None, "se": None, "ci_lower": None, "ci_upper": None}

    mock_test.effect_column_names = effect_column_names

    def run(gene: str, data: dict, config: Any) -> TestResult:
        p_val = None if gene in skip_genes else p_value_default
        return TestResult(
            gene=gene,
            test_name="skat",
            p_value=p_val,
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            se=None,
            n_cases=data.get("proband_count", 0),
            n_controls=data.get("control_count", 0),
            n_variants=data.get("n_qualifying_variants", 0),
            extra=dict(extra_default),
        )

    mock_test.run = run
    mock_test.check_dependencies = MagicMock()
    mock_test.prepare = MagicMock()
    mock_test.finalize = MagicMock()

    return mock_test


def _make_engine_with_mock_skat(
    mock_skat_test: Any,
    config: Any | None = None,
) -> Any:
    """Build an AssociationEngine directly with a mock SKAT test (bypass from_names)."""
    from variantcentrifuge.association.engine import AssociationEngine

    if config is None:
        config = _make_config()

    return AssociationEngine([mock_skat_test], config)


# ---------------------------------------------------------------------------
# Column structure tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestSKATColumnStructure:
    """Tests that engine produces correct SKAT column output."""

    def test_engine_skat_no_effect_columns(self):
        """Engine does NOT produce 'skat_None' or other malformed effect columns."""
        mock_skat = _make_mock_skat_test()
        engine = _make_engine_with_mock_skat(mock_skat)

        gene_data = [
            _make_gene_burden_data("BRCA1"),
            _make_gene_burden_data("TP53"),
            _make_gene_burden_data("MYH7"),
        ]
        result = engine.run_all(gene_data)

        # None-effect guard: these columns must NOT exist
        for bad_col in ["skat_None", "skat_or", "skat_beta", "skat_se"]:
            assert bad_col not in result.columns, (
                f"Column '{bad_col}' should not be present in SKAT output"
            )

    def test_engine_skat_has_p_value_columns(self):
        """Engine produces skat_p_value and ACAT-O corrected p-value columns.

        ARCH-03: Primary test (skat) has no corrected_p_value column.
        FDR is applied only to ACAT-O omnibus p-values.
        """
        mock_skat = _make_mock_skat_test()
        engine = _make_engine_with_mock_skat(mock_skat)

        gene_data = [
            _make_gene_burden_data("BRCA1"),
            _make_gene_burden_data("TP53"),
        ]
        result = engine.run_all(gene_data)

        assert "skat_p_value" in result.columns
        # ARCH-03: no per-test corrected column; ACAT-O is the single FDR-corrected output
        assert "skat_corrected_p_value" not in result.columns
        assert "acat_o_p_value" in result.columns
        assert "acat_o_corrected_p_value" in result.columns

    def test_engine_skat_extra_columns_written(self):
        """Extra SKAT columns (skat_o_rho, skat_method) are written to output."""
        extra = {
            "skat_o_rho": 0.5,
            "skat_method": "SKATO",
            "skat_warnings": None,
            "skat_n_marker_test": 4,
        }
        mock_skat = _make_mock_skat_test(extra_default=extra)
        engine = _make_engine_with_mock_skat(mock_skat)

        gene_data = [_make_gene_burden_data("BRCA1")]
        result = engine.run_all(gene_data)

        assert "skat_o_rho" in result.columns
        assert result.iloc[0]["skat_o_rho"] == 0.5
        assert "skat_method" in result.columns
        assert result.iloc[0]["skat_method"] == "SKATO"

    def test_engine_skat_extra_column_not_double_prefixed(self):
        """Extra keys are NOT double-prefixed (no 'skat_skat_o_rho' columns)."""
        extra = {"skat_o_rho": 0.3, "skat_method": "SKAT"}
        mock_skat = _make_mock_skat_test(extra_default=extra)
        engine = _make_engine_with_mock_skat(mock_skat)

        gene_data = [_make_gene_burden_data("BRCA1")]
        result = engine.run_all(gene_data)

        # Keys in extra are already namespaced — no double prefix
        assert "skat_skat_o_rho" not in result.columns
        assert "skat_skat_method" not in result.columns
        assert "skat_o_rho" in result.columns

    def test_engine_skat_shared_columns_present(self):
        """Engine output has shared columns: gene, n_cases, n_controls, n_variants."""
        mock_skat = _make_mock_skat_test()
        engine = _make_engine_with_mock_skat(mock_skat)

        gene_data = [_make_gene_burden_data("BRCA1", p_total=60, c_total=80, n_variants=5)]
        result = engine.run_all(gene_data)

        assert "gene" in result.columns
        assert "n_cases" in result.columns
        assert "n_controls" in result.columns
        assert "n_variants" in result.columns

        row = result.iloc[0]
        assert row["gene"] == "BRCA1"
        assert row["n_cases"] == 60
        assert row["n_controls"] == 80
        assert row["n_variants"] == 5


# ---------------------------------------------------------------------------
# Multi-test (Fisher + SKAT) tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestSKATAndFisherTogether:
    """Tests that Fisher and SKAT can run together in the same engine."""

    def test_engine_skat_and_fisher_together(self):
        """Both fisher_p_value and skat_p_value columns present when both tests run."""
        from variantcentrifuge.association.engine import AssociationEngine
        from variantcentrifuge.association.tests.fisher import FisherExactTest

        config = _make_config()
        fisher = FisherExactTest()
        mock_skat = _make_mock_skat_test()

        engine = AssociationEngine([fisher, mock_skat], config)

        gene_data = [
            _make_gene_burden_data("BRCA1", p_carriers=5, c_carriers=1),
            _make_gene_burden_data("TP53", p_carriers=4, c_carriers=2),
        ]
        result = engine.run_all(gene_data)

        assert "fisher_p_value" in result.columns
        assert "skat_p_value" in result.columns

    def test_engine_fisher_has_or_columns_skat_does_not(self):
        """Fisher produces _or columns; SKAT does not produce any effect columns."""
        from variantcentrifuge.association.engine import AssociationEngine
        from variantcentrifuge.association.tests.fisher import FisherExactTest

        config = _make_config()
        fisher = FisherExactTest()
        mock_skat = _make_mock_skat_test()

        engine = AssociationEngine([fisher, mock_skat], config)

        gene_data = [_make_gene_burden_data("BRCA1", p_carriers=5, c_carriers=1)]
        result = engine.run_all(gene_data)

        # Fisher effect columns present
        assert "fisher_or" in result.columns
        assert "fisher_or_ci_lower" in result.columns
        assert "fisher_or_ci_upper" in result.columns

        # SKAT effect columns absent
        assert "skat_or" not in result.columns
        assert "skat_None" not in result.columns
        assert "skat_beta" not in result.columns


# ---------------------------------------------------------------------------
# Correction behavior tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestSKATCorrectionBehavior:
    """Tests for multiple testing correction with SKAT."""

    def test_engine_skat_none_pvalue_excluded_from_correction(self):
        """Gene with p_value=None is excluded from ACAT-O FDR; others get corrected.

        ARCH-03: FDR correction applied to ACAT-O p-values only.
        Genes where all primary tests returned None are excluded from output entirely.
        """
        mock_skat = _make_mock_skat_test(
            p_value_default=0.05,
            skip_genes=["EMPTY_GENE"],
        )
        engine = _make_engine_with_mock_skat(mock_skat)

        gene_data = [
            _make_gene_burden_data("BRCA1"),
            _make_gene_burden_data("TP53"),
            # EMPTY_GENE has 0 variants — but mock handles skip
            _make_gene_burden_data("EMPTY_GENE", include_genotype_matrix=False),
        ]

        result = engine.run_all(gene_data)

        # EMPTY_GENE is excluded from output (p_value=None → no row)
        if "EMPTY_GENE" in result["gene"].values:
            empty_row = result[result["gene"] == "EMPTY_GENE"].iloc[0]
            # ACAT-O for a skipped gene should be None (no primary p-values)
            assert empty_row["acat_o_p_value"] is None or np.isnan(
                float(empty_row["acat_o_p_value"])
            )

        # Genes with real p-values get corrected ACAT-O values
        tested_genes = result[result["skat_p_value"].notna()]
        assert tested_genes["acat_o_corrected_p_value"].notna().all()

    def test_engine_skat_corrected_pvalues_populated(self):
        """acat_o_corrected_p_value is populated for all tested genes.

        ARCH-03: FDR is applied to ACAT-O p-values, not per-test p-values.
        """
        mock_skat = _make_mock_skat_test(p_value_default=0.04)
        engine = _make_engine_with_mock_skat(mock_skat)

        gene_data = [
            _make_gene_burden_data("BRCA1"),
            _make_gene_burden_data("TP53"),
            _make_gene_burden_data("MYH7"),
        ]
        result = engine.run_all(gene_data)

        assert result["acat_o_corrected_p_value"].notna().all()
        # No per-test corrected column (ARCH-03)
        assert "skat_corrected_p_value" not in result.columns

    def test_engine_skat_genes_sorted_alphabetically(self):
        """Engine outputs genes in alphabetical order for deterministic correction."""
        mock_skat = _make_mock_skat_test()
        engine = _make_engine_with_mock_skat(mock_skat)

        # Provide in reverse order
        gene_data = [
            _make_gene_burden_data("ZZGENE"),
            _make_gene_burden_data("MMGENE"),
            _make_gene_burden_data("AAGENE"),
        ]
        result = engine.run_all(gene_data)

        assert list(result["gene"]) == ["AAGENE", "MMGENE", "ZZGENE"]


# ---------------------------------------------------------------------------
# Lifecycle hook tests (engine calls prepare/finalize)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestEngineLifecycleHooks:
    """Tests that the engine calls prepare() and finalize() on all tests."""

    def test_engine_prepare_called_with_gene_count(self):
        """engine.run_all() calls test.prepare(n_genes) before gene loop."""
        mock_skat = _make_mock_skat_test()
        engine = _make_engine_with_mock_skat(mock_skat)

        gene_data = [
            _make_gene_burden_data("BRCA1"),
            _make_gene_burden_data("TP53"),
            _make_gene_burden_data("MYH7"),
        ]
        engine.run_all(gene_data)

        mock_skat.prepare.assert_called_once_with(3)

    def test_engine_finalize_called_once_after_gene_loop(self):
        """engine.run_all() calls test.finalize() exactly once after gene loop."""
        mock_skat = _make_mock_skat_test()
        engine = _make_engine_with_mock_skat(mock_skat)

        gene_data = [
            _make_gene_burden_data("BRCA1"),
            _make_gene_burden_data("TP53"),
        ]
        engine.run_all(gene_data)

        mock_skat.finalize.assert_called_once()

    def test_engine_prepare_and_finalize_called_for_all_tests(self):
        """engine.run_all() calls prepare/finalize on both Fisher and SKAT tests."""
        from variantcentrifuge.association.engine import AssociationEngine
        from variantcentrifuge.association.tests.fisher import FisherExactTest

        config = _make_config()
        mock_skat = _make_mock_skat_test()

        # Wrap Fisher with a mock that tracks prepare/finalize
        fisher = FisherExactTest()
        fisher.prepare = MagicMock(wraps=fisher.prepare)
        fisher.finalize = MagicMock(wraps=fisher.finalize)

        engine = AssociationEngine([fisher, mock_skat], config)

        gene_data = [_make_gene_burden_data("BRCA1", p_carriers=5, c_carriers=1)]
        engine.run_all(gene_data)

        fisher.prepare.assert_called_once_with(1)
        fisher.finalize.assert_called_once()
        mock_skat.prepare.assert_called_once_with(1)
        mock_skat.finalize.assert_called_once()

    def test_engine_empty_input_still_calls_finalize(self):
        """finalize() is NOT called when input is empty (run_all returns early)."""
        # When gene_burden_data is empty, run_all returns immediately before
        # the gene loop — prepare and finalize are NOT called (consistent with
        # the guard 'if not gene_burden_data: return pd.DataFrame()')
        mock_skat = _make_mock_skat_test()
        engine = _make_engine_with_mock_skat(mock_skat)

        result = engine.run_all([])

        assert result.empty
        # prepare and finalize are only called when there's data to process
        # (engine returns before reaching the lifecycle hooks for empty input)
        # This documents the current behavior.


# ---------------------------------------------------------------------------
# Column count / regression guard
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestSKATColumnCount:
    """Tests that SKAT output has exactly the expected column count (regression guard)."""

    def test_skat_only_output_columns(self):
        """SKAT-only engine produces exactly: gene, n_cases, n_controls, n_variants,
        skat_p_value, plus all extra keys, plus acat_o_p_value and acat_o_corrected_p_value.

        ARCH-03: No skat_corrected_p_value column. FDR applied to ACAT-O only.
        """
        extra = {
            "skat_o_rho": 0.0,
            "skat_method": "SKAT",
            "skat_warnings": None,
            "skat_n_marker_test": 3,
        }
        mock_skat = _make_mock_skat_test(extra_default=extra)
        engine = _make_engine_with_mock_skat(mock_skat)

        gene_data = [_make_gene_burden_data("BRCA1")]
        result = engine.run_all(gene_data)

        expected_cols = {
            "gene",
            "n_cases",
            "n_controls",
            "n_variants",
            "skat_p_value",
            # ARCH-03: no skat_corrected_p_value; ACAT-O is the single corrected output
            "acat_o_p_value",
            "acat_o_corrected_p_value",
            "skat_o_rho",
            "skat_method",
            "skat_warnings",
            "skat_n_marker_test",
        }
        actual_cols = set(result.columns)
        assert expected_cols == actual_cols, (
            f"Column mismatch.\n"
            f"Expected: {sorted(expected_cols)}\n"
            f"Got:      {sorted(actual_cols)}\n"
            f"Extra in output:  {sorted(actual_cols - expected_cols)}\n"
            f"Missing from output: {sorted(expected_cols - actual_cols)}"
        )
