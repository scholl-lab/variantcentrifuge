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
