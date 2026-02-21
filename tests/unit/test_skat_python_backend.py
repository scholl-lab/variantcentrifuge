"""
Unit tests for PythonSKATBackend.

Tests the pure Python SKAT backend using numpy, scipy, and statsmodels.
No R or rpy2 required — all computations are native Python.

Covers requirements: SKAT-05, SKAT-06, SKAT-07, SKAT-10

Test categories:
- Environment detection and version logging
- Null model fitting (binary and quantitative)
- SKAT test_gene: basic, rank-deficient, edge cases
- SKAT-O test_gene: rho grid search
- Burden test_gene: analytical p-value
- Cleanup (no-op)
"""

from __future__ import annotations

import logging
from unittest.mock import patch

import numpy as np
import pytest

from variantcentrifuge.association.backends.python_backend import _SKATO_RHO_GRID, PythonSKATBackend

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def backend():
    """Create and initialise a PythonSKATBackend."""
    b = PythonSKATBackend()
    b.detect_environment()
    return b


@pytest.fixture
def binary_null_model(backend):
    """Null model for a 50-sample binary trait (no covariates)."""
    rng = np.random.default_rng(42)
    n = 50
    phenotype = rng.integers(0, 2, n).astype(np.float64)
    return backend.fit_null_model(phenotype, None, "binary")


@pytest.fixture
def quantitative_null_model(backend):
    """Null model for a 100-sample quantitative trait (no covariates)."""
    rng = np.random.default_rng(43)
    n = 100
    phenotype = rng.standard_normal(n)
    return backend.fit_null_model(phenotype, None, "quantitative")


@pytest.fixture
def standard_geno():
    """50-sample, 5-variant genotype matrix with rank >= 2."""
    rng = np.random.default_rng(100)
    return rng.integers(0, 3, (50, 5)).astype(np.float64)


@pytest.fixture
def rank2_geno():
    """Larger genotype matrix guaranteed to have rank >= 2 for SKAT-O tests."""
    rng = np.random.default_rng(200)
    # Ensure rank >= 2 by construction
    geno = rng.integers(0, 3, (100, 10)).astype(np.float64)
    return geno


# ---------------------------------------------------------------------------
# Environment detection tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestDetectEnvironment:
    """Tests for PythonSKATBackend.detect_environment()."""

    def test_detect_environment_succeeds(self):
        """detect_environment() completes without raising."""
        b = PythonSKATBackend()
        # Must not raise
        b.detect_environment()
        # Version strings should be populated
        assert b._numpy_version != ""
        assert b._scipy_version != ""
        assert b._statsmodels_version != ""

    def test_detect_environment_sets_version_strings(self):
        """detect_environment() sets non-empty version strings for all deps."""
        b = PythonSKATBackend()
        b.detect_environment()
        assert b._numpy_version not in ("", "<unknown>"), (
            f"numpy version not set: {b._numpy_version}"
        )
        assert b._scipy_version not in ("", "<unknown>"), (
            f"scipy version not set: {b._scipy_version}"
        )
        assert b._statsmodels_version not in ("", "<unknown>"), (
            f"statsmodels version not set: {b._statsmodels_version}"
        )

    def test_log_environment_logs_versions(self, caplog):
        """log_environment() emits INFO with numpy/scipy/statsmodels versions."""
        b = PythonSKATBackend()
        b.detect_environment()
        with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
            b.log_environment()

        info_messages = "\n".join(r.message for r in caplog.records if r.levelno == logging.INFO)
        assert b._numpy_version in info_messages, "numpy version not found in logs"
        assert b._scipy_version in info_messages, "scipy version not found in logs"
        assert b._statsmodels_version in info_messages, "statsmodels version not found in logs"


# ---------------------------------------------------------------------------
# Null model fitting tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestFitNullModel:
    """Tests for PythonSKATBackend.fit_null_model()."""

    def test_fit_null_model_binary_trait_type(self, backend):
        """Binary phenotype -> NullModelResult with trait_type='binary'."""
        rng = np.random.default_rng(1)
        phenotype = rng.integers(0, 2, 50).astype(np.float64)
        result = backend.fit_null_model(phenotype, None, "binary")
        assert result.trait_type == "binary"

    def test_fit_null_model_binary_sigma2_is_one(self, backend):
        """Binary trait: sigma2 = 1.0 by SKAT convention."""
        rng = np.random.default_rng(2)
        phenotype = rng.integers(0, 2, 50).astype(np.float64)
        result = backend.fit_null_model(phenotype, None, "binary")
        assert result.extra["sigma2"] == pytest.approx(1.0)

    def test_fit_null_model_binary_residuals_in_extra(self, backend):
        """Binary: extra dict contains 'residuals' key with correct shape."""
        rng = np.random.default_rng(3)
        phenotype = rng.integers(0, 2, 50).astype(np.float64)
        result = backend.fit_null_model(phenotype, None, "binary")
        assert "residuals" in result.extra
        assert result.extra["residuals"].shape == (50,)

    def test_fit_null_model_quantitative_trait_type(self, backend):
        """Quantitative phenotype -> NullModelResult with trait_type='quantitative'."""
        rng = np.random.default_rng(4)
        phenotype = rng.standard_normal(100)
        result = backend.fit_null_model(phenotype, None, "quantitative")
        assert result.trait_type == "quantitative"

    def test_fit_null_model_quantitative_sigma2_positive(self, backend):
        """Quantitative trait: sigma2 > 0."""
        rng = np.random.default_rng(5)
        phenotype = rng.standard_normal(100)
        result = backend.fit_null_model(phenotype, None, "quantitative")
        assert result.extra["sigma2"] > 0.0

    def test_fit_null_model_n_samples_correct(self, backend):
        """n_samples in NullModelResult matches phenotype length."""
        rng = np.random.default_rng(6)
        phenotype = rng.integers(0, 2, 75).astype(np.float64)
        result = backend.fit_null_model(phenotype, None, "binary")
        assert result.n_samples == 75

    def test_fit_null_model_with_covariates_changes_residuals(self, backend):
        """Covariates affect residuals — result differs from no-covariate model."""
        rng = np.random.default_rng(7)
        n = 100
        phenotype = rng.integers(0, 2, n).astype(np.float64)
        covariates = rng.standard_normal((n, 2))

        result_no_cov = backend.fit_null_model(phenotype, None, "binary")
        result_with_cov = backend.fit_null_model(phenotype, covariates, "binary")

        # Residuals must differ when covariates are included
        resid_diff = np.abs(
            result_no_cov.extra["residuals"] - result_with_cov.extra["residuals"]
        ).max()
        assert resid_diff > 1e-8, (
            f"Residuals should differ with covariates, max diff = {resid_diff:.2e}"
        )

    def test_fit_null_model_residuals_are_response(self, backend):
        """
        Residuals are y - mu_hat (response residuals), not deviance or Pearson.

        R SKAT convention uses resid_response. Verify by computing expected residuals
        directly: extra['residuals'] == phenotype - extra['mu_hat'].
        """
        rng = np.random.default_rng(8)
        phenotype = rng.integers(0, 2, 60).astype(np.float64)
        result = backend.fit_null_model(phenotype, None, "binary")
        mu_hat = result.extra["mu_hat"]
        expected_residuals = phenotype - mu_hat
        np.testing.assert_allclose(
            result.extra["residuals"],
            expected_residuals,
            atol=1e-10,
            err_msg="Residuals should be y - mu_hat (response residuals)",
        )

    def test_fit_null_model_resets_genes_processed(self, backend):
        """fit_null_model() resets the internal gene counter to 0."""
        backend._genes_processed = 999
        rng = np.random.default_rng(9)
        phenotype = rng.integers(0, 2, 50).astype(np.float64)
        backend.fit_null_model(phenotype, None, "binary")
        assert backend._genes_processed == 0


# ---------------------------------------------------------------------------
# SKAT test_gene tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestSKATTestGene:
    """Tests for PythonSKATBackend.test_gene() with method='SKAT'."""

    def test_skat_basic_returns_dict_with_required_keys(
        self, backend, binary_null_model, standard_geno
    ):
        """Basic SKAT call returns dict with p_value, p_method, n_variants."""
        result = backend.test_gene("GENE1", standard_geno, binary_null_model, "SKAT", (1.0, 25.0))
        assert "p_value" in result
        assert "p_method" in result
        assert "p_converged" in result
        assert result["n_variants"] == 5

    def test_skat_basic_p_value_float_or_none(self, backend, binary_null_model, standard_geno):
        """SKAT p_value is a float in [0, 1] or None."""
        result = backend.test_gene("GENE1", standard_geno, binary_null_model, "SKAT", (1.0, 25.0))
        p = result["p_value"]
        if p is not None:
            assert isinstance(p, float)
            assert 0.0 <= p <= 1.0, f"p_value {p:.6e} not in [0, 1]"

    def test_skat_p_method_in_expected_set(self, backend, binary_null_model, standard_geno):
        """p_method must be 'davies', 'saddlepoint', 'liu', or None (skipped)."""
        result = backend.test_gene("GENE1", standard_geno, binary_null_model, "SKAT", (1.0, 25.0))
        assert result["p_method"] in {"davies", "saddlepoint", "liu", None}

    def test_skat_p_converged_is_bool(self, backend, binary_null_model, standard_geno):
        """p_converged must be a bool."""
        result = backend.test_gene("GENE1", standard_geno, binary_null_model, "SKAT", (1.0, 25.0))
        assert isinstance(result["p_converged"], bool)

    def test_skat_rank_deficient_returns_none_pvalue(self, backend, binary_null_model):
        """Rank-deficient matrix: p_value=None, skip_reason='rank_deficient'."""
        rng = np.random.default_rng(11)
        n = 50
        # All columns identical => rank = 1 < 2
        base_col = rng.integers(0, 3, n).astype(np.float64)
        geno = np.column_stack([base_col] * 5)
        result = backend.test_gene("GENE_RD", geno, binary_null_model, "SKAT", (1.0, 25.0))
        assert result["p_value"] is None
        assert result["skip_reason"] == "rank_deficient"

    def test_skat_zero_variants_returns_none(self, backend, binary_null_model):
        """Zero-variant matrix (0 columns): p_value=None."""
        geno = np.zeros((50, 0))
        result = backend.test_gene("GENE_EMPTY", geno, binary_null_model, "SKAT", (1.0, 25.0))
        assert result["p_value"] is None
        assert result["n_variants"] == 0

    def test_skat_single_variant_rank_deficient(self, backend, binary_null_model):
        """Single-variant matrix: rank < 2, should be skipped."""
        rng = np.random.default_rng(12)
        geno = rng.integers(0, 3, (50, 1)).astype(np.float64)
        result = backend.test_gene("GENE_1V", geno, binary_null_model, "SKAT", (1.0, 25.0))
        assert result["p_value"] is None
        assert result["skip_reason"] == "rank_deficient"

    def test_skat_all_zero_genotypes_graceful(self, backend, binary_null_model):
        """All-zero genotype matrix: handled gracefully (p_value=None or 1.0)."""
        geno = np.zeros((50, 5))
        result = backend.test_gene("GENE_ZEROS", geno, binary_null_model, "SKAT", (1.0, 25.0))
        # Should not raise; p_value is None (rank_deficient) or 1.0
        assert result["p_value"] is None or result["p_value"] == pytest.approx(1.0), (
            f"Expected None or 1.0 for all-zero genotypes, got {result['p_value']}"
        )

    def test_skat_eigenvalue_driver_evr(self, backend, binary_null_model, standard_geno):
        """Eigenvalue computation uses scipy.linalg.eigh with driver='evr'."""
        import scipy.linalg

        calls: list[dict] = []
        original_eigh = scipy.linalg.eigh

        def spy_eigh(*args, **kwargs):
            calls.append(kwargs)
            return original_eigh(*args, **kwargs)

        patch_target = "variantcentrifuge.association.backends.python_backend.scipy.linalg.eigh"
        with patch(patch_target, spy_eigh):
            backend.test_gene("GENE_EVR", standard_geno, binary_null_model, "SKAT", (1.0, 25.0))

        evr_calls = [c for c in calls if c.get("driver") == "evr"]
        assert len(evr_calls) >= 1, f"Expected eigh with driver='evr', got kwargs: {calls}"

    def test_skat_increments_genes_processed(self, backend, binary_null_model, standard_geno):
        """Each test_gene call increments the internal gene counter."""
        initial = backend._genes_processed
        backend.test_gene("G1", standard_geno, binary_null_model, "SKAT", (1.0, 25.0))
        backend.test_gene("G2", standard_geno, binary_null_model, "SKAT", (1.0, 25.0))
        assert backend._genes_processed == initial + 2

    def test_skat_unknown_method_defaults_to_skat(self, backend, binary_null_model, standard_geno):
        """Unknown method falls back to SKAT without raising."""
        result = backend.test_gene(
            "G_UNK", standard_geno, binary_null_model, "UNKNOWN", (1.0, 25.0)
        )
        # Should return a result (not crash)
        assert "p_value" in result
        assert result["n_variants"] == 5


# ---------------------------------------------------------------------------
# SKAT-O test_gene tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestSKATOTestGene:
    """Tests for PythonSKATBackend.test_gene() with method='SKATO'."""

    def test_skato_returns_rho_in_unit_interval(self, backend, rank2_geno):
        """SKAT-O returns optimal rho in [0, 1]."""
        rng = np.random.default_rng(20)
        n = rank2_geno.shape[0]
        phenotype = rng.integers(0, 2, n).astype(np.float64)
        null = backend.fit_null_model(phenotype, None, "binary")
        result = backend.test_gene("GENE_O", rank2_geno, null, "SKATO", (1.0, 25.0))
        rho = result["rho"]
        if rho is not None:
            assert 0.0 <= rho <= 1.0, f"rho={rho} not in [0, 1]"

    def test_skato_rho_from_grid(self, backend, rank2_geno):
        """SKAT-O optimal rho must be from the 7-element fixed grid."""
        rng = np.random.default_rng(21)
        n = rank2_geno.shape[0]
        phenotype = rng.integers(0, 2, n).astype(np.float64)
        null = backend.fit_null_model(phenotype, None, "binary")
        result = backend.test_gene("GENE_O", rank2_geno, null, "SKATO", (1.0, 25.0))
        rho = result["rho"]
        if rho is not None:
            assert any(abs(rho - r) < 1e-10 for r in _SKATO_RHO_GRID), (
                f"rho={rho} not in grid {_SKATO_RHO_GRID}"
            )

    def test_skato_rho_grid_has_7_elements(self):
        """SKAT-O rho grid should have exactly 7 elements [0, 0.01, ..., 1.0]."""
        assert len(_SKATO_RHO_GRID) == 7
        assert _SKATO_RHO_GRID[0] == 0.0
        assert _SKATO_RHO_GRID[-1] == 1.0

    def test_skato_returns_valid_p_value(self, backend, rank2_geno):
        """SKAT-O p_value is a valid float in [0, 1] or None."""
        rng = np.random.default_rng(22)
        n = rank2_geno.shape[0]
        phenotype = rng.integers(0, 2, n).astype(np.float64)
        null = backend.fit_null_model(phenotype, None, "binary")
        result = backend.test_gene("GENE_O", rank2_geno, null, "SKATO", (1.0, 25.0))
        p = result["p_value"]
        if p is not None:
            assert 0.0 <= p <= 1.0, f"SKAT-O p={p:.6e} not in [0, 1]"

    def test_skato_rank_deficient_returns_none(self, backend):
        """Rank-deficient matrix: SKAT-O returns p_value=None."""
        rng = np.random.default_rng(23)
        n = 100
        col = rng.integers(0, 3, n).astype(np.float64)
        geno = np.column_stack([col] * 5)  # rank = 1
        phenotype = rng.integers(0, 2, n).astype(np.float64)
        null = backend.fit_null_model(phenotype, None, "binary")
        result = backend.test_gene("GENE_O_RD", geno, null, "SKATO", (1.0, 25.0))
        assert result["p_value"] is None
        assert result["skip_reason"] == "rank_deficient"


# ---------------------------------------------------------------------------
# Burden test_gene tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestBurdenTestGene:
    """Tests for PythonSKATBackend.test_gene() with method='Burden'."""

    def test_burden_returns_p_value(self, backend, binary_null_model, standard_geno):
        """Burden test returns a valid p_value."""
        result = backend.test_gene(
            "GENE_B", standard_geno, binary_null_model, "Burden", (1.0, 25.0)
        )
        p = result["p_value"]
        if p is not None:
            assert 0.0 <= p <= 1.0, f"Burden p={p:.6e} not in [0, 1]"

    def test_burden_p_method_is_analytical(self, backend, binary_null_model, standard_geno):
        """Burden p_method should be 'analytical'."""
        result = backend.test_gene(
            "GENE_B", standard_geno, binary_null_model, "Burden", (1.0, 25.0)
        )
        if result["p_value"] is not None:
            assert result["p_method"] == "analytical", (
                f"Expected 'analytical', got {result['p_method']}"
            )

    def test_burden_no_rho(self, backend, binary_null_model, standard_geno):
        """Burden test does not produce a rho value."""
        result = backend.test_gene(
            "GENE_B", standard_geno, binary_null_model, "Burden", (1.0, 25.0)
        )
        assert result["rho"] is None

    def test_burden_all_zero_geno_handles_gracefully(self, backend, binary_null_model):
        """Burden on all-zero genotype handles gracefully (no crash)."""
        geno = np.zeros((50, 5))
        result = backend.test_gene("GENE_BZ", geno, binary_null_model, "Burden", (1.0, 25.0))
        # Should not raise; p_value may be None due to zero variance
        assert "p_value" in result


# ---------------------------------------------------------------------------
# Cleanup test
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestCleanup:
    """Tests for PythonSKATBackend.cleanup()."""

    def test_cleanup_noop(self, backend):
        """cleanup() does not raise — no external resources to release."""
        # Should not raise
        backend.cleanup()

    def test_cleanup_logs_timing(self, backend, caplog):
        """cleanup() logs INFO with gene count and timing."""
        with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
            backend.cleanup()
        info_msgs = [r.message for r in caplog.records if r.levelno == logging.INFO]
        assert any("cleanup" in m.lower() or "genes" in m.lower() for m in info_msgs), (
            f"Expected cleanup log, got: {info_msgs}"
        )
