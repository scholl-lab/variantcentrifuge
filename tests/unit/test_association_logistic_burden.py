"""
Unit tests for LogisticBurdenTest with manual statsmodels.Logit validation.

The critical acceptance criterion for BURDEN-01: LogisticBurdenTest must produce
the same p-value, beta, and SE as calling statsmodels.Logit directly on the same
inputs (intercept + burden score + covariates).

Tests also cover:
- effect_size = beta (log-odds), not OR
- SE as first-class field (result.se)
- Firth fallback on perfect separation (PERFECT_SEPARATION warning)
- Warning codes: FIRTH_CONVERGE_FAIL, ZERO_CARRIERS_ONE_GROUP, LOW_CARRIER_COUNT
- NO_GENOTYPE_MATRIX early return
- Carrier count statistics in extra dict
- Never silently omits gene (always returns TestResult)
"""

from __future__ import annotations

import numpy as np
import pytest
import statsmodels.api as sm

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.tests.logistic_burden import LogisticBurdenTest

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def synthetic_binary_data() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Synthetic data: 100 samples, 5 variants, binary phenotype with known signal.

    Returns
    -------
    geno : np.ndarray, shape (100, 5)
        Sparse rare variant genotype matrix.
    phenotype : np.ndarray, shape (100,)
        Binary outcome correlated with burden.
    mafs : np.ndarray, shape (5,)
        Per-variant MAFs.
    """
    np.random.seed(42)
    n = 100
    geno = np.zeros((n, 5))
    for v in range(5):
        carriers = np.random.choice(n, size=int(n * 0.05), replace=False)
        geno[carriers, v] = 1
    # Phenotype: correlated with burden
    burden_score = geno.sum(axis=1)
    prob = 1 / (1 + np.exp(-(burden_score * 1.5 - 0.5)))
    phenotype = (np.random.random(n) < prob).astype(float)
    mafs = geno.mean(axis=0) / 2
    return geno, phenotype, mafs


@pytest.fixture
def default_config() -> AssociationConfig:
    """Default AssociationConfig with uniform weights for predictable burden scores."""
    return AssociationConfig(variant_weights="uniform")


# ---------------------------------------------------------------------------
# Helper to build contingency_data
# ---------------------------------------------------------------------------


def _make_contingency(
    geno: np.ndarray,
    phenotype: np.ndarray,
    mafs: np.ndarray,
    covariate_matrix: np.ndarray | None = None,
) -> dict:
    """Build contingency_data dict for LogisticBurdenTest.run()."""
    n_cases = int(phenotype.sum())
    n_controls = int((phenotype == 0).sum())
    return {
        "genotype_matrix": geno,
        "variant_mafs": mafs,
        "phenotype_vector": phenotype,
        "covariate_matrix": covariate_matrix,
        "proband_count": n_cases,
        "control_count": n_controls,
        "n_qualifying_variants": geno.shape[1],
    }


def test_logistic_burden_passes_score_values_to_weight_resolver(
    monkeypatch, synthetic_binary_data
) -> None:
    from variantcentrifuge.association import weights as weights_module

    geno, phenotype, mafs = synthetic_binary_data
    score_values = np.array([0.2, 0.4, 0.6, 0.8, 1.0], dtype=np.float64)
    seen_kwargs: dict[str, object] = {}

    def fake_get_weights(mafs_arg, weight_spec, **kwargs):
        seen_kwargs.update(kwargs)
        return np.ones(len(mafs_arg), dtype=np.float64)

    monkeypatch.setattr(weights_module, "get_weights", fake_get_weights)
    config = AssociationConfig(variant_weights="column:nephro_candidate_score")
    contingency = _make_contingency(geno, phenotype, mafs)
    contingency["score_values"] = score_values

    result = LogisticBurdenTest().run("GENE1", contingency, config)

    assert result.p_value is not None
    assert seen_kwargs["score_values"] is score_values


# ---------------------------------------------------------------------------
# Manual statsmodels.Logit validation (BURDEN-01)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLogisticBurdenMatchesStatsmodels:
    """LogisticBurdenTest output matches manual statsmodels.Logit — BURDEN-01."""

    def test_logistic_burden_matches_manual_statsmodels(
        self, synthetic_binary_data, default_config
    ) -> None:
        """p-value, beta, se match direct sm.Logit call on same inputs (rel=1e-6)."""
        geno, phenotype, mafs = synthetic_binary_data
        config = default_config

        # Build burden manually (uniform weights = sum of dosages)
        from variantcentrifuge.association.weights import get_weights

        weights = get_weights(mafs, config.variant_weights)
        burden = geno @ weights

        # Manual statsmodels reference
        design_manual = sm.add_constant(burden.reshape(-1, 1))
        fit_manual = sm.Logit(phenotype, design_manual).fit(disp=False, maxiter=100)
        manual_p = float(fit_manual.pvalues[1])
        manual_beta = float(fit_manual.params[1])
        manual_se = float(fit_manual.bse[1])

        # LogisticBurdenTest
        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, config)

        assert result.p_value is not None, "Expected valid p-value"
        assert result.p_value == pytest.approx(manual_p, rel=1e-6), (
            f"p-value mismatch: got {result.p_value}, expected {manual_p}"
        )
        assert result.extra["beta"] == pytest.approx(manual_beta, rel=1e-6), (
            f"beta mismatch: got {result.extra['beta']}, expected {manual_beta}"
        )
        assert result.extra["se"] == pytest.approx(manual_se, rel=1e-6), (
            f"se mismatch: got {result.extra['se']}, expected {manual_se}"
        )

    def test_logistic_burden_with_covariates_matches_manual(
        self, synthetic_binary_data, default_config
    ) -> None:
        """With covariates, burden coefficient matches direct sm.Logit on same design matrix."""
        geno, phenotype, mafs = synthetic_binary_data
        config = default_config
        np.random.seed(7)
        n = len(phenotype)

        # 2-column covariate matrix (age-like, PC1-like)
        covariate_matrix = np.column_stack(
            [
                np.random.normal(50, 10, n),
                np.random.normal(0, 1, n),
            ]
        )

        from variantcentrifuge.association.weights import get_weights

        weights = get_weights(mafs, config.variant_weights)
        burden = geno @ weights

        # Manual statsmodels: intercept + burden + 2 covariates
        design_manual = sm.add_constant(burden.reshape(-1, 1))
        design_manual = np.column_stack([design_manual, covariate_matrix])
        fit_manual = sm.Logit(phenotype, design_manual).fit(disp=False, maxiter=100)
        # Burden coefficient is at index 1 (0=intercept, 1=burden, 2..=covariates)
        manual_p = float(fit_manual.pvalues[1])
        manual_beta = float(fit_manual.params[1])

        # LogisticBurdenTest with covariates
        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs, covariate_matrix)
        result = test.run("GENE1", contingency, config)

        assert result.p_value is not None
        assert result.p_value == pytest.approx(manual_p, rel=1e-6), (
            f"p-value with covariates mismatch: got {result.p_value}, expected {manual_p}"
        )
        assert result.extra["beta"] == pytest.approx(manual_beta, rel=1e-6)

    def test_logistic_burden_effect_size_is_beta(
        self, synthetic_binary_data, default_config
    ) -> None:
        """effect_size = beta (log-odds), not exp(beta). CI on log-odds scale."""
        geno, phenotype, mafs = synthetic_binary_data
        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, default_config)

        assert result.p_value is not None
        beta = result.extra["beta"]

        # effect_size IS the raw beta (log-odds), not OR
        assert result.effect_size == pytest.approx(beta, rel=1e-6), (
            f"effect_size should be beta={beta}, got {result.effect_size}"
        )
        # CI on log-odds scale can be negative
        if result.ci_lower is not None:
            assert result.ci_lower < result.ci_upper


# ---------------------------------------------------------------------------
# Warning codes (separation, convergence, carrier counts)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLogisticBurdenWarnings:
    """Warning code tests for LogisticBurdenTest."""

    def test_logistic_burden_perfect_separation_firth(
        self, default_config: AssociationConfig
    ) -> None:
        """Perfect separation triggers PERFECT_SEPARATION; Firth produces finite estimates."""
        np.random.seed(0)
        n = 40
        # All carriers are cases (perfect separation)
        geno = np.zeros((n, 3))
        # Assign carriers only to cases (S0-S9)
        geno[:10, 0] = 1  # first 10 = cases, all carriers
        phenotype = np.zeros(n)
        phenotype[:20] = 1  # cases are first 20
        mafs = np.array([0.125, 0.0, 0.0])

        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("SEP_GENE", contingency, default_config)

        # PERFECT_SEPARATION should be flagged
        warnings = result.extra.get("warnings", [])
        assert "PERFECT_SEPARATION" in warnings, (
            f"Expected PERFECT_SEPARATION warning, got: {warnings}"
        )
        # Firth should produce a result (p_value is not None)
        assert result.p_value is not None, (
            "Firth fallback should produce a finite p-value for perfect separation"
        )
        # Effect size should be finite (not inf or nan)
        if result.effect_size is not None:
            assert np.isfinite(result.effect_size), "Effect size should be finite after Firth"

    def test_logistic_burden_low_carrier_count(self, default_config: AssociationConfig) -> None:
        """Fewer than 3 carriers total -> LOW_CARRIER_COUNT warning."""
        np.random.seed(1)
        n = 50
        geno = np.zeros((n, 2))
        # Only 2 carriers total
        geno[0, 0] = 1
        geno[1, 0] = 1
        phenotype = np.zeros(n)
        phenotype[:25] = 1  # 25 cases, 25 controls
        mafs = np.array([0.02, 0.0])

        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("LOW_CARRIER_GENE", contingency, default_config)

        warnings = result.extra.get("warnings", [])
        assert "LOW_CARRIER_COUNT" in warnings, (
            f"Expected LOW_CARRIER_COUNT warning, got: {warnings}"
        )

    def test_logistic_burden_zero_carriers_one_group(
        self, default_config: AssociationConfig
    ) -> None:
        """All carriers in cases, none in controls -> separation-related warning."""
        np.random.seed(2)
        n = 60
        geno = np.zeros((n, 2))
        # 5 carriers, all in cases (S0-S4)
        geno[:5, 0] = 1
        phenotype = np.zeros(n)
        phenotype[:30] = 1  # 30 cases, 30 controls

        mafs = np.array([5.0 / (2 * n), 0.0])

        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("ZERO_CTRL_GENE", contingency, default_config)

        warnings = result.extra.get("warnings", [])
        # Either PERFECT_SEPARATION or ZERO_CARRIERS_ONE_GROUP expected
        assert "PERFECT_SEPARATION" in warnings or "ZERO_CARRIERS_ONE_GROUP" in warnings, (
            f"Expected separation/zero-carriers warning, got: {warnings}"
        )

    def test_logistic_burden_no_genotype_matrix(self, default_config: AssociationConfig) -> None:
        """Missing genotype_matrix key -> p_value=None and NO_GENOTYPE_MATRIX warning."""
        contingency = {
            "proband_count": 50,
            "control_count": 50,
            "n_qualifying_variants": 3,
            # No genotype_matrix key!
        }

        test = LogisticBurdenTest()
        result = test.run("FISHER_ONLY_GENE", contingency, default_config)

        assert result.p_value is None, "Expected p_value=None when genotype_matrix absent"
        warnings = result.extra.get("warnings", [])
        assert "NO_GENOTYPE_MATRIX" in warnings, (
            f"Expected NO_GENOTYPE_MATRIX warning, got: {warnings}"
        )

    def test_logistic_burden_result_is_always_test_result(
        self, default_config: AssociationConfig
    ) -> None:
        """Any valid contingency data always returns TestResult (never raises exception)."""
        from variantcentrifuge.association.base import TestResult

        # Minimal valid data
        n = 30
        geno = np.zeros((n, 2))
        geno[:5, 0] = 1
        phenotype = np.zeros(n)
        phenotype[:15] = 1
        mafs = np.array([5.0 / (2 * n), 0.0])

        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("ROBUST_GENE", contingency, default_config)

        # Must always return TestResult, never raise
        assert isinstance(result, TestResult)
        assert result.gene == "ROBUST_GENE"
        assert result.test_name == "logistic_burden"


# ---------------------------------------------------------------------------
# Carrier count statistics
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLogisticBurdenCarrierCounts:
    """Carrier count statistics in extra dict."""

    def test_logistic_burden_carrier_counts(self, synthetic_binary_data, default_config) -> None:
        """n_carriers, n_carriers_cases, n_carriers_controls match manual computation."""
        geno, phenotype, mafs = synthetic_binary_data

        from variantcentrifuge.association.weights import get_weights

        weights = get_weights(mafs, default_config.variant_weights)
        burden = geno @ weights
        carriers = burden > 0

        expected_n_carriers = int(carriers.sum())
        expected_cases = int((carriers & (phenotype == 1)).sum())
        expected_controls = int((carriers & (phenotype == 0)).sum())

        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, default_config)

        assert result.extra["n_carriers"] == expected_n_carriers, (
            f"n_carriers mismatch: got {result.extra['n_carriers']}, expected {expected_n_carriers}"
        )
        assert result.extra["n_carriers_cases"] == expected_cases
        assert result.extra["n_carriers_controls"] == expected_controls

    def test_logistic_burden_never_omits_gene(self, default_config: AssociationConfig) -> None:
        """Gene with no carriers returns TestResult (not exception)."""
        n = 40
        geno = np.zeros((n, 3))  # all ref — no carriers
        phenotype = np.zeros(n)
        phenotype[:20] = 1
        mafs = np.zeros(3)

        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("ZERO_CARRIER_GENE", contingency, default_config)

        # Should return TestResult, not raise exception
        from variantcentrifuge.association.base import TestResult

        assert isinstance(result, TestResult), "Expected TestResult, never an exception"
        # p_value is None (no carriers -> no signal) or valid
        assert result.p_value is None or isinstance(result.p_value, float)

    def test_logistic_burden_result_metadata(self, synthetic_binary_data, default_config) -> None:
        """TestResult has correct gene name and test_name."""
        geno, phenotype, mafs = synthetic_binary_data
        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("MY_GENE", contingency, default_config)

        assert result.gene == "MY_GENE"
        assert result.test_name == "logistic_burden"

    def test_logistic_burden_effect_size_is_beta_not_or(
        self, synthetic_binary_data, default_config
    ) -> None:
        """effect_size is raw beta (log-odds), not exp(beta)."""
        geno, phenotype, mafs = synthetic_binary_data
        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, default_config)

        if result.p_value is not None:
            beta = result.extra["beta"]
            assert result.effect_size == pytest.approx(beta, rel=1e-6)
            # Verify it's NOT exp(beta)
            assert result.effect_size != pytest.approx(float(np.exp(beta)), rel=1e-2)

    def test_logistic_burden_ci_bounds_consistent_with_effect_size(
        self, synthetic_binary_data, default_config
    ) -> None:
        """ci_lower <= effect_size <= ci_upper (log-odds scale)."""
        geno, phenotype, mafs = synthetic_binary_data
        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, default_config)

        if result.p_value is not None and result.ci_lower is not None:
            assert result.ci_lower <= result.effect_size <= result.ci_upper, (
                f"Beta CI inconsistent: [{result.ci_lower}, {result.ci_upper}] "
                f"does not contain beta={result.effect_size}"
            )

    def test_logistic_burden_se_field(self, synthetic_binary_data, default_config) -> None:
        """result.se is populated and matches extra['se']."""
        geno, phenotype, mafs = synthetic_binary_data
        test = LogisticBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, default_config)

        if result.p_value is not None:
            assert result.se is not None, "SE should be populated for successful fit"
            assert result.se > 0, "SE must be positive"
            assert result.se == pytest.approx(result.extra["se"], rel=1e-6)

    def test_logistic_burden_check_dependencies(self) -> None:
        """check_dependencies does not raise (statsmodels is available)."""
        test = LogisticBurdenTest()
        # Should not raise ImportError
        test.check_dependencies()

    def test_logistic_burden_name_property(self) -> None:
        """LogisticBurdenTest.name == 'logistic_burden'."""
        test = LogisticBurdenTest()
        assert test.name == "logistic_burden"


# ---------------------------------------------------------------------------
# Issue #143 Regression Tests: Firth step-halving and convergence safeguards
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestFirthConvergenceSafeguards:
    """Regression tests for Firth convergence and step-halving safeguards (issue #143)."""

    def test_unconverged_max_iter_returns_none_and_propagates_warning(self):
        """Unconverged max_iter returns None and forced fallback reports FIRTH_CONVERGE_FAIL."""
        from variantcentrifuge.association.tests.logistic_burden import (
            _firth_logistic,
        )

        y = np.array([0.0, 0.0, 1.0, 1.0])
        x = np.column_stack([np.ones(4), y])

        # 1. max_iter=0 and max_iter=1 must return None
        assert _firth_logistic(y, x, max_iter=0) is None
        assert _firth_logistic(y, x, max_iter=1) is None

        # 2. Forced fallback through LogisticBurdenTest
        config = AssociationConfig(
            trait_type="binary",
            firth_max_iter=1,  # insufficient iterations to converge
        )
        test = LogisticBurdenTest()
        contingency_data = {
            "genotype_matrix": y[:, np.newaxis],  # perfect separation: cases=1, controls=0
            "phenotype_vector": y,
            "variant_mafs": np.array([0.5]),
            "proband_count": 2,
            "control_count": 2,
            "n_qualifying_variants": 1,
        }
        result = test.run("TEST_SEPARATED", contingency_data, config)
        assert result.p_value is None
        assert result.effect_size is None
        assert result.se is None
        assert "FIRTH_CONVERGE_FAIL" in result.extra.get("warnings", [])

    def test_singular_trial_returns_negative_infinity(self):
        """Singular trial point returns -np.inf and cannot outrank a finite objective."""
        from variantcentrifuge.association.tests.logistic_burden import (
            _firth_loglik,
        )

        y = np.array([0.0, 0.0, 1.0, 1.0])
        x = np.column_stack([np.ones(4), y])

        # Singular point with det(X'WX) == 0
        singular_ll = _firth_loglik(np.array([-100.0, 200.0]), y, x)
        assert singular_ll == -np.inf

        # Valid initial objective is strictly greater than -np.inf
        initial_ll = _firth_loglik(np.zeros(2), y, x)
        assert np.isfinite(initial_ll)
        assert singular_ll < initial_ll

    def test_exhausted_step_halving_returns_none(self, monkeypatch):
        """Exhausted step-halving cannot commit an unchecked candidate
        or manufacture convergence.
        """
        from variantcentrifuge.association.tests.logistic_burden import (
            _firth_logistic,
        )

        y = np.array([0.0, 0.0, 1.0, 1.0])
        x = np.column_stack([np.ones(4), y])

        # Mock _firth_loglik so candidate steps always produce -inf (failing step-halving)
        import variantcentrifuge.association.tests.logistic_burden as lb_mod

        calls = 0

        def mock_loglik(beta, y_arr, x_arr):
            nonlocal calls
            calls += 1
            if calls == 1:
                return 0.0  # old_ll
            return -np.inf  # candidate always worse

        monkeypatch.setattr(lb_mod, "_firth_loglik", mock_loglik)
        result = _firth_logistic(y, x, max_iter=10)
        assert result is None

    def test_convergent_separated_firth_matches_analytic(self):
        """Convergent separated design matches analytic Firth estimate."""
        from variantcentrifuge.association.tests.logistic_burden import (
            _firth_logistic,
        )

        # 2x2 table: a=2, b=0, c=0, d=2 => beta = 2*log(5) ~ 3.2188758
        y = np.array([1.0, 1.0, 0.0, 0.0])
        x = np.column_stack([np.ones(4), y])
        res = _firth_logistic(y, x, max_iter=50, tol=1e-5)
        assert res is not None
        analytic_beta = 2.0 * float(np.log(5.0))
        np.testing.assert_allclose(res.params[1], analytic_beta, rtol=1e-4)

    def test_worsening_step_with_large_loglik_not_accepted(self, monkeypatch):
        """Worsening step with large loglik (-1000 vs -1000.001) is rejected by rtol=0.0."""
        import variantcentrifuge.association.tests.logistic_burden as lb_mod
        from variantcentrifuge.association.tests.logistic_burden import (
            _firth_logistic,
        )

        y = np.r_[np.ones(251), np.zeros(249), np.ones(249), np.zeros(251)]
        x = np.column_stack([np.ones(1000), np.r_[np.ones(500), np.zeros(500)]])

        def mock_loglik(beta, y_arr, x_arr):
            if np.all(beta == 0):
                return -1000.0
            return -1000.001

        monkeypatch.setattr(lb_mod, "_firth_loglik", mock_loglik)
        result = _firth_logistic(y, x, max_iter=25, tol=0.01)
        assert result is None

    def test_firth_tp53_two_group_predictor_scale_reconstruction(self):
        """TP53 two-group predictor-scale reconstruction matches analytic Firth estimate."""
        from variantcentrifuge.association.tests.logistic_burden import (
            _firth_logistic,
        )

        weight = 1.874434234570756e-8
        # 12 case homozygotes, 2 control homozygotes, 10 control wildtypes
        y = np.array([1.0] * 12 + [0.0] * 2 + [0.0] * 10)
        x_pred = np.array([2.0 * weight] * 12 + [2.0 * weight] * 2 + [0.0] * 10)
        x = np.column_stack([np.ones(len(y)), x_pred])

        analytic_beta = (float(np.log(12.5 / 2.5)) - float(np.log(0.5 / 10.5))) / (2.0 * weight)
        res = _firth_logistic(y, x, max_iter=50, tol=0.1)
        assert res is not None
        np.testing.assert_allclose(res.params[1], analytic_beta, rtol=1e-4)


def test_firth_converges_when_step_halving_hits_loglik_rounding_noise():
    """An uncentred covariate (birth-year-like) leaves a flat ridge where the Newton step
    exceeds tol but no step changes the penalised log-likelihood beyond rounding noise.
    That is the optimum (Newton decrement ~0), not a failure (#198 review)."""
    from variantcentrifuge.association.tests.logistic_burden import _firth_logistic

    rng = np.random.default_rng(1)
    n = 2000
    for _rep in range(25):  # replicate 24 of this stream stalled step-halving at tol 1e-7
        y = (rng.uniform(size=n) < 0.3).astype(float)
        year = rng.normal(1970, 10, n)
        burden = rng.binomial(2, 0.002, (n, 6)).sum(1) * 23.0
    x = np.column_stack([np.ones(n), burden, year])
    fit = _firth_logistic(y, x, max_iter=50, tol=1e-7)
    assert fit is not None
    # same optimum as the centred (well-conditioned) parametrisation
    xc = np.column_stack([np.ones(n), burden, year - 1970])
    ref = _firth_logistic(y, xc, max_iter=50, tol=1e-10)
    assert fit.params[1] == pytest.approx(ref.params[1], rel=1e-6)
    assert fit.bse[1] == pytest.approx(ref.bse[1], rel=1e-6)
