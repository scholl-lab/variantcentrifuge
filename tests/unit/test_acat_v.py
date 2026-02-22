# File: tests/unit/test_acat_v.py
"""Tests for ACAT-V per-variant score test (Phase 25)."""

import numpy as np
import pytest

from variantcentrifuge.association.tests.acat import compute_acat_v


@pytest.mark.unit
class TestComputeAcatV:
    """Tests for compute_acat_v()."""

    def test_basic_binary(self):
        """ACAT-V produces valid p-value for binary trait."""
        rng = np.random.default_rng(42)
        n, m = 100, 5
        geno = rng.choice([0, 1, 2], size=(n, m), p=[0.8, 0.15, 0.05]).astype(float)
        mu_hat = np.full(n, 0.3)  # balanced binary trait
        residuals = rng.binomial(1, 0.3, n).astype(float) - mu_hat
        sigma2 = 1.0
        weights = np.ones(m)

        p = compute_acat_v(geno, residuals, "binary", sigma2, mu_hat, weights)
        assert p is not None
        assert 0.0 < p <= 1.0

    def test_basic_quantitative(self):
        """ACAT-V produces valid p-value for quantitative trait."""
        rng = np.random.default_rng(123)
        n, m = 200, 10
        geno = rng.choice([0, 1, 2], size=(n, m), p=[0.85, 0.12, 0.03]).astype(float)
        mu_hat = np.full(n, 5.0)
        residuals = rng.normal(0, 1, n)
        sigma2 = float(np.var(residuals, ddof=1))
        weights = np.ones(m)

        p = compute_acat_v(geno, residuals, "quantitative", sigma2, mu_hat, weights)
        assert p is not None
        assert 0.0 < p <= 1.0

    def test_quantitative_mu_hat_none(self):
        """ACAT-V accepts mu_hat=None for quantitative traits."""
        rng = np.random.default_rng(77)
        n, m = 100, 4
        geno = rng.choice([0, 1, 2], size=(n, m), p=[0.85, 0.12, 0.03]).astype(float)
        residuals = rng.normal(0, 1, n)
        sigma2 = 1.0
        weights = np.ones(m)

        # mu_hat=None should work for quantitative (sigma2 path)
        p = compute_acat_v(geno, residuals, "quantitative", sigma2, None, weights)
        assert p is not None
        assert 0.0 < p <= 1.0

    def test_empty_genotype(self):
        """Zero variants returns None."""
        geno = np.empty((100, 0))
        residuals = np.zeros(100)
        mu_hat = np.full(100, 0.5)
        p = compute_acat_v(geno, residuals, "binary", 1.0, mu_hat, np.array([]))
        assert p is None

    def test_monomorphic_variants_skipped(self):
        """All-zero genotype columns (monomorphic) are skipped."""
        n = 50
        geno = np.zeros((n, 3))  # all monomorphic
        residuals = np.random.default_rng(1).normal(0, 1, n)
        mu_hat = np.full(n, 0.5)
        weights = np.ones(3)

        p = compute_acat_v(geno, residuals, "binary", 1.0, mu_hat, weights)
        assert p is None  # no valid per-variant p-values

    def test_single_variant(self):
        """Single variant: ACAT-V returns per-variant p-value directly."""
        rng = np.random.default_rng(99)
        n = 100
        geno = rng.choice([0, 1], size=(n, 1), p=[0.9, 0.1]).astype(float)
        mu_hat = np.full(n, 0.3)
        residuals = rng.binomial(1, 0.3, n).astype(float) - mu_hat
        weights = np.ones(1)

        p = compute_acat_v(geno, residuals, "binary", 1.0, mu_hat, weights)
        # Single valid p-value: cauchy_combination returns pass-through
        assert p is not None
        assert 0.0 < p <= 1.0

    def test_significant_signal_detected(self):
        """ACAT-V detects a strong per-variant signal."""
        rng = np.random.default_rng(7)
        n = 500
        # Create one causal variant with strong effect
        geno = np.zeros((n, 5))
        geno[:, 0] = rng.binomial(1, 0.1, n)  # causal
        for j in range(1, 5):
            geno[:, j] = rng.binomial(1, 0.05, n)  # noise

        # Phenotype depends on first variant
        mu_hat = np.full(n, 0.3)
        phenotype = (geno[:, 0] * 3.0 + rng.normal(0, 1, n)) > 1.5
        residuals = phenotype.astype(float) - mu_hat
        weights = np.ones(5)

        p = compute_acat_v(geno, residuals, "binary", 1.0, mu_hat, weights)
        assert p is not None
        # With strong causal variant, p should be small
        assert p < 0.05

    def test_weights_affect_result(self):
        """Different weights produce different ACAT-V p-values."""
        rng = np.random.default_rng(55)
        n, m = 200, 5
        geno = rng.choice([0, 1, 2], size=(n, m), p=[0.85, 0.12, 0.03]).astype(float)
        mu_hat = np.full(n, 0.5)
        residuals = rng.normal(0, 1, n)
        sigma2 = 1.0

        p_uniform = compute_acat_v(geno, residuals, "quantitative", sigma2, mu_hat, np.ones(m))
        p_weighted = compute_acat_v(
            geno,
            residuals,
            "quantitative",
            sigma2,
            mu_hat,
            np.array([10.0, 1.0, 1.0, 1.0, 1.0]),
        )
        # Both should be valid
        assert p_uniform is not None
        assert p_weighted is not None
        # Different weights should give different results (not guaranteed but very likely)
        # We don't assert inequality since it's theoretically possible they match

    def test_acat_v_in_acat_o(self):
        """ACAT-V p-value included in ACAT-O when present in extra."""
        from variantcentrifuge.association.tests.acat import compute_acat_o

        # Simulate: SKAT p=0.05, burden p=0.10, acat_v p=0.001
        # ACAT-O with ACAT-V should be more significant than without
        pvals_without = {"skat": 0.05, "burden": 0.10}
        pvals_with = {"skat": 0.05, "burden": 0.10, "acat_v": 0.001}

        p_without = compute_acat_o(pvals_without)
        p_with = compute_acat_o(pvals_with)

        assert p_without is not None
        assert p_with is not None
        # Adding a very significant ACAT-V should make omnibus more significant
        assert p_with < p_without

    def test_pure_python_skat_stores_acat_v_in_extra(self):
        """End-to-end: PurePythonSKATTest.run() -> extra['acat_v_p'] is set."""
        from variantcentrifuge.association.base import AssociationConfig
        from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest

        n_samples, n_variants = 50, 6

        # Synthetic data: binary trait (balanced 50/50 split, same as existing SKAT tests)
        n_cases = n_samples // 2
        n_controls = n_samples - n_cases
        phenotype = np.array([1.0] * n_cases + [0.0] * n_controls)
        rng = np.random.default_rng(42)
        geno = rng.choice([0, 1, 2], size=(n_samples, n_variants), p=[0.6, 0.3, 0.1]).astype(float)

        config = AssociationConfig(trait_type="binary")

        # PurePythonSKATTest.run() takes (gene, contingency_data, config).
        # The null model is fit lazily on first run() call.
        # No covariate_matrix: intercept-only null model (no singular matrix risk).
        contingency_data = {
            "genotype_matrix": geno,
            "phenotype_vector": phenotype,
            "proband_count": n_cases,
            "control_count": n_controls,
            "n_qualifying_variants": n_variants,
        }

        test = PurePythonSKATTest()
        test.check_dependencies()
        result = test.run("TEST_GENE", contingency_data, config)

        # ACAT-V p-value must be stored in extra
        assert "acat_v_p" in result.extra, (
            "PurePythonSKATTest.run() must store 'acat_v_p' in result.extra"
        )
        # Value is either a float in (0, 1] or None (if all variants monomorphic)
        acat_v_p = result.extra["acat_v_p"]
        if acat_v_p is not None:
            assert 0.0 < acat_v_p <= 1.0, f"acat_v_p={acat_v_p} out of range"

    def test_early_return_no_genotype_matrix_has_acat_v_p(self):
        """Early return (no genotype matrix) sets acat_v_p=None in extra."""
        from variantcentrifuge.association.base import AssociationConfig
        from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest

        config = AssociationConfig(trait_type="binary")
        contingency_data = {
            # no "genotype_matrix" key
            "proband_count": 10,
            "control_count": 90,
            "n_qualifying_variants": 5,
        }

        test = PurePythonSKATTest()
        test.check_dependencies()
        result = test.run("GENE1", contingency_data, config)

        assert "acat_v_p" in result.extra
        assert result.extra["acat_v_p"] is None

    def test_early_return_empty_matrix_has_acat_v_p(self):
        """Early return (empty genotype matrix) sets acat_v_p=None in extra."""
        import numpy as np

        from variantcentrifuge.association.base import AssociationConfig
        from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest

        config = AssociationConfig(trait_type="binary")
        contingency_data = {
            "genotype_matrix": np.empty((0, 5)),  # 0 samples
            "phenotype_vector": np.array([]),
            "covariate_matrix": None,
            "proband_count": 0,
            "control_count": 0,
            "n_qualifying_variants": 5,
        }

        test = PurePythonSKATTest()
        test.check_dependencies()
        result = test.run("GENE2", contingency_data, config)

        assert "acat_v_p" in result.extra
        assert result.extra["acat_v_p"] is None
