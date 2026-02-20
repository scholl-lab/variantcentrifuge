"""
Unit tests for RSKATTest with mocked RSKATBackend.

Tests the RSKATTest wrapper logic without invoking rpy2 at all.
RSKATBackend is fully mocked — tests verify the AssociationTest
protocol implementation (null model caching, skip conditions,
extra columns, lifecycle hooks, GC triggering).

Covers requirements: SKAT-01 through SKAT-09
"""

from __future__ import annotations

import logging
import time
from typing import Any
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------


def _make_config(**kwargs: Any) -> Any:
    """Build an AssociationConfig with sensible defaults."""
    from variantcentrifuge.association.base import AssociationConfig

    defaults = {
        "trait_type": "binary",
        "skat_method": "SKAT",
        "variant_weights": "beta:1,25",
    }
    defaults.update(kwargs)
    return AssociationConfig(**defaults)


def _make_geno_data(
    gene: str = "BRCA1",
    n_samples: int = 60,
    n_variants: int = 4,
    n_cases: int = 30,
    n_controls: int = 30,
    include_genotype_matrix: bool = True,
) -> dict[str, Any]:
    """Build a contingency_data dict including genotype matrix."""
    rng = np.random.default_rng(42)
    geno = rng.integers(0, 3, size=(n_samples, n_variants)).astype(float)
    phenotype = np.array([1.0] * n_cases + [0.0] * n_controls)
    mafs = geno.mean(axis=0) / 2.0

    data: dict[str, Any] = {
        "GENE": gene,
        "proband_count": n_cases,
        "control_count": n_controls,
        "proband_carrier_count": 5,
        "control_carrier_count": 2,
        "proband_allele_count": 5,
        "control_allele_count": 2,
        "n_qualifying_variants": n_variants,
        "phenotype_vector": phenotype,
        "variant_mafs": mafs,
        "covariate_matrix": None,
    }
    if include_genotype_matrix:
        data["genotype_matrix"] = geno
    return data


def _make_skat_test_with_mock_backend() -> tuple[Any, MagicMock]:
    """
    Create an RSKATTest instance with a mocked RSKATBackend.

    Returns
    -------
    (test, mock_backend)
        test: RSKATTest ready to call run()
        mock_backend: The injected mock backend
    """
    from variantcentrifuge.association.tests.skat_r import RSKATTest

    test = RSKATTest()

    mock_backend = MagicMock(name="RSKATBackend")
    mock_backend.GC_INTERVAL = 100  # match real constant

    # Default test_gene return value
    mock_backend.test_gene.return_value = {
        "p_value": 0.03,
        "rho": None,
        "n_variants": 4,
        "n_marker_test": 4,
        "warnings": [],
    }

    # fit_null_model returns a NullModelResult
    from variantcentrifuge.association.backends.base import NullModelResult

    mock_backend.fit_null_model.return_value = NullModelResult(
        model=MagicMock(name="null_model"),
        trait_type="binary",
        n_samples=60,
        adjustment=True,
    )

    # Inject backend directly (bypass check_dependencies)
    test._backend = mock_backend

    return test, mock_backend


# ---------------------------------------------------------------------------
# Basic property tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestRSKATTestBasicProperties:
    """Tests for RSKATTest name and effect_column_names."""

    def test_name_is_skat(self):
        """RSKATTest.name returns 'skat'."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        test = RSKATTest()
        assert test.name == "skat"

    def test_effect_column_names_all_none(self):
        """effect_column_names() returns all None values (SKAT has no effect size)."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        test = RSKATTest()
        col_names = test.effect_column_names()

        assert col_names["effect"] is None
        assert col_names["se"] is None
        assert col_names["ci_lower"] is None
        assert col_names["ci_upper"] is None

    def test_effect_column_names_returns_all_four_keys(self):
        """effect_column_names() returns a dict with exactly the four expected keys."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        test = RSKATTest()
        col_names = test.effect_column_names()

        assert set(col_names.keys()) == {"effect", "se", "ci_lower", "ci_upper"}


# ---------------------------------------------------------------------------
# Null model caching tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestNullModelCaching:
    """Tests for null model fitted exactly once across multiple gene calls."""

    def test_null_model_fitted_once(self):
        """fit_null_model is called exactly once across 3 different gene calls."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        config = _make_config()

        for gene in ["BRCA1", "TP53", "MYH7"]:
            data = _make_geno_data(gene=gene)
            test.run(gene, data, config)

        mock_backend.fit_null_model.assert_called_once()

    def test_null_model_uses_trait_type_from_config(self):
        """fit_null_model is called with trait_type from AssociationConfig."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        config = _make_config(trait_type="quantitative")

        data = _make_geno_data()
        test.run("GENE1", data, config)

        call_kwargs = mock_backend.fit_null_model.call_args.kwargs
        assert call_kwargs.get("trait_type") == "quantitative"

    def test_null_model_reused_second_call(self):
        """After first run, _null_model is set and reused."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        config = _make_config()
        data = _make_geno_data()

        assert test._null_model is None
        test.run("GENE1", data, config)
        null_model_after_first = test._null_model
        test.run("GENE2", data, config)
        null_model_after_second = test._null_model

        # Same object, not refitted
        assert null_model_after_first is null_model_after_second
        assert mock_backend.fit_null_model.call_count == 1


# ---------------------------------------------------------------------------
# Skip condition tests (p_value=None)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestSkipConditions:
    """Tests for conditions where run() returns p_value=None."""

    def test_run_no_genotype_matrix_returns_none_pvalue(self):
        """run() with no genotype_matrix in contingency_data returns p_value=None."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        config = _make_config()

        data = _make_geno_data(include_genotype_matrix=False)
        result = test.run("GENE1", data, config)

        assert result.p_value is None
        mock_backend.test_gene.assert_not_called()

    def test_run_zero_variants_returns_none_pvalue(self):
        """run() with genotype_matrix of shape (100, 0) returns p_value=None."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        config = _make_config()

        data = _make_geno_data(n_variants=0)
        data["genotype_matrix"] = np.zeros((60, 0))
        result = test.run("GENE1", data, config)

        assert result.p_value is None
        mock_backend.test_gene.assert_not_called()

    def test_run_no_phenotype_returns_none_pvalue(self):
        """run() without phenotype_vector returns p_value=None."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        config = _make_config()

        data = _make_geno_data()
        data.pop("phenotype_vector", None)  # no key at all
        result = test.run("GENE1", data, config)

        assert result.p_value is None
        mock_backend.test_gene.assert_not_called()

    def test_run_backend_none_raises_runtime_error(self):
        """run() raises RuntimeError when _backend is None (check_dependencies not called)."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        test = RSKATTest()
        # _backend is None — no check_dependencies called
        config = _make_config()
        data = _make_geno_data()

        with pytest.raises(RuntimeError, match="check_dependencies"):
            test.run("GENE1", data, config)


# ---------------------------------------------------------------------------
# Extra column output tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestExtraColumns:
    """Tests for extra columns written by RSKATTest."""

    def test_run_returns_skat_extra_columns(self):
        """run() packs SKAT-specific data into result.extra."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        config = _make_config(skat_method="SKATO")

        mock_backend.test_gene.return_value = {
            "p_value": 0.03,
            "rho": 0.5,
            "warnings": ["w1"],
            "n_variants": 4,
            "n_marker_test": 5,
        }

        data = _make_geno_data()
        result = test.run("GENE1", data, config)

        assert "skat_o_rho" in result.extra
        assert result.extra["skat_o_rho"] == 0.5
        assert "skat_method" in result.extra
        assert result.extra["skat_method"] == "SKATO"
        assert "skat_n_marker_test" in result.extra
        assert result.extra["skat_n_marker_test"] == 5

    def test_run_skat_warnings_joined_in_extra(self):
        """R warnings are joined with '; ' in skat_warnings extra column."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        config = _make_config()

        mock_backend.test_gene.return_value = {
            "p_value": 0.04,
            "rho": None,
            "warnings": ["moment adjustment applied", "convergence warning"],
            "n_variants": 3,
            "n_marker_test": 3,
        }

        data = _make_geno_data()
        result = test.run("GENE1", data, config)

        assert "skat_warnings" in result.extra
        assert "moment adjustment applied" in result.extra["skat_warnings"]
        assert "convergence warning" in result.extra["skat_warnings"]

    def test_run_skat_warnings_none_when_empty(self):
        """skat_warnings is None when R returns no warnings."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        config = _make_config()

        mock_backend.test_gene.return_value = {
            "p_value": 0.05,
            "rho": None,
            "warnings": [],
            "n_variants": 4,
            "n_marker_test": 4,
        }

        data = _make_geno_data()
        result = test.run("GENE1", data, config)

        assert result.extra["skat_warnings"] is None

    def test_run_effect_fields_all_none(self):
        """run() result has effect_size, se, ci_lower, ci_upper all None."""
        test, _ = _make_skat_test_with_mock_backend()
        config = _make_config()
        data = _make_geno_data()

        result = test.run("GENE1", data, config)

        assert result.effect_size is None
        assert result.se is None
        assert result.ci_lower is None
        assert result.ci_upper is None


# ---------------------------------------------------------------------------
# Lifecycle hook tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLifecycleHooks:
    """Tests for prepare() and finalize() lifecycle hooks."""

    def test_prepare_large_panel_warning(self, caplog):
        """prepare(3000) emits WARNING about large gene panel."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        test = RSKATTest()

        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            test.prepare(3000)

        warning_msgs = [r for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warning_msgs) >= 1
        assert any("3000" in m.message for m in warning_msgs)

    def test_prepare_small_panel_no_warning(self, caplog):
        """prepare(100) does NOT emit large panel WARNING."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        test = RSKATTest()

        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            test.prepare(100)

        warning_msgs = [r for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warning_msgs) == 0

    def test_prepare_sets_log_interval_for_large_panel(self):
        """prepare(500) sets log_interval to max(10, min(50, 500//10)) = 50."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        test = RSKATTest()
        test.prepare(500)

        assert test._log_interval == 50

    def test_prepare_sets_log_interval_for_small_panel(self):
        """prepare(30) sets log_interval to max(10, min(50, 30//10)) = 10."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        test = RSKATTest()
        test.prepare(30)

        assert test._log_interval == 10

    def test_finalize_logs_timing(self, caplog):
        """finalize() logs INFO with 'R SKAT complete' and gene count."""
        test, _ = _make_skat_test_with_mock_backend()
        test._genes_processed = 42
        test._start_time = time.time() - 5.0  # 5 seconds ago

        with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
            test.finalize()

        info_msgs = [r for r in caplog.records if r.levelno == logging.INFO]
        assert any("R SKAT complete" in m.message for m in info_msgs)
        assert any("42" in m.message for m in info_msgs)

    def test_finalize_calls_backend_cleanup(self):
        """finalize() calls backend.cleanup()."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        test._start_time = time.time()

        test.finalize()

        mock_backend.cleanup.assert_called_once()

    def test_finalize_no_backend_no_error(self):
        """finalize() does not raise when _backend is None."""
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        test = RSKATTest()
        test._start_time = time.time()
        # _backend is None
        test.finalize()  # Should not raise


# ---------------------------------------------------------------------------
# GC interval tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestGCInterval:
    """Tests for periodic GC triggering every 100 genes in RSKATTest."""

    def test_gc_triggered_every_100_genes(self):
        """backend._run_r_gc is called once after 150 gene calls (at gene 100)."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        # Set GC_INTERVAL on mock backend class
        mock_backend.GC_INTERVAL = 100

        # Mock _check_r_heap to return a safe value
        mock_backend._check_r_heap.return_value = 100.0  # MB, below threshold

        config = _make_config()

        # Patch RSKATBackend.GC_INTERVAL used in RSKATTest.run()
        with patch("variantcentrifuge.association.tests.skat_r.RSKATBackend.GC_INTERVAL", 100):
            for i in range(150):
                data = _make_geno_data(gene=f"GENE_{i:04d}")
                test.run(f"GENE_{i:04d}", data, config)

        assert mock_backend._run_r_gc.call_count == 1

    def test_gc_not_triggered_before_interval(self):
        """backend._run_r_gc is not called when fewer than 100 genes processed."""
        test, mock_backend = _make_skat_test_with_mock_backend()
        config = _make_config()

        with patch("variantcentrifuge.association.tests.skat_r.RSKATBackend.GC_INTERVAL", 100):
            for i in range(50):
                data = _make_geno_data(gene=f"GENE_{i:04d}")
                test.run(f"GENE_{i:04d}", data, config)

        assert mock_backend._run_r_gc.call_count == 0


# ---------------------------------------------------------------------------
# Progress tracking tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestProgressTracking:
    """Tests for gene counting and progress logging."""

    def test_genes_processed_increments(self):
        """_genes_processed increments on each successful run() call."""
        test, _ = _make_skat_test_with_mock_backend()
        config = _make_config()

        assert test._genes_processed == 0
        test.run("GENE1", _make_geno_data("GENE1"), config)
        assert test._genes_processed == 1
        test.run("GENE2", _make_geno_data("GENE2"), config)
        assert test._genes_processed == 2

    def test_genes_processed_not_incremented_for_skipped(self):
        """_genes_processed is NOT incremented when run() skips (no genotype_matrix)."""
        test, _ = _make_skat_test_with_mock_backend()
        config = _make_config()

        data = _make_geno_data(include_genotype_matrix=False)
        test.run("GENE1", data, config)

        assert test._genes_processed == 0

    def test_run_result_has_correct_gene_name(self):
        """TestResult.gene matches the gene argument passed to run()."""
        test, _ = _make_skat_test_with_mock_backend()
        config = _make_config()

        result = test.run("MY_GENE", _make_geno_data("MY_GENE"), config)
        assert result.gene == "MY_GENE"

    def test_run_result_has_correct_test_name(self):
        """TestResult.test_name matches 'skat'."""
        test, _ = _make_skat_test_with_mock_backend()
        config = _make_config()

        result = test.run("GENE1", _make_geno_data("GENE1"), config)
        assert result.test_name == "skat"
