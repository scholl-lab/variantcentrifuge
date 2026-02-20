"""
Unit tests for LinearBurdenTest with manual statsmodels.OLS validation.

The critical acceptance criterion for BURDEN-02: LinearBurdenTest must produce
the same beta, SE, and p-value as calling statsmodels.OLS directly on the same
inputs (intercept + burden score + covariates).

Tests also cover:
- Beta = effect_size (not OR)
- CI coverage at 95% level for synthetic data with known signal
- NO_GENOTYPE_MATRIX early return
- Carrier count statistics in extra dict
"""

from __future__ import annotations

import numpy as np
import pytest
import statsmodels.api as sm

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.tests.linear_burden import LinearBurdenTest

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def synthetic_quantitative_data() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Synthetic data: 100 samples, 5 variants, continuous phenotype with known signal.

    Returns
    -------
    geno : np.ndarray, shape (100, 5)
        Sparse rare variant genotype matrix.
    phenotype : np.ndarray, shape (100,)
        Continuous outcome = 5 + 2*burden + noise.
    mafs : np.ndarray, shape (5,)
        Per-variant MAFs.
    """
    np.random.seed(42)
    n = 100
    geno = np.zeros((n, 5))
    for v in range(5):
        carriers = np.random.choice(n, size=int(n * 0.1), replace=False)
        geno[carriers, v] = 1
    burden = geno.sum(axis=1)
    phenotype = 5.0 + 2.0 * burden + np.random.normal(0, 1, n)
    mafs = geno.mean(axis=0) / 2
    return geno, phenotype, mafs


@pytest.fixture
def default_config() -> AssociationConfig:
    """Default AssociationConfig with uniform weights for predictable burden scores."""
    return AssociationConfig(variant_weights="uniform")


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _make_contingency(
    geno: np.ndarray,
    phenotype: np.ndarray,
    mafs: np.ndarray,
    covariate_matrix: np.ndarray | None = None,
) -> dict:
    """Build contingency_data dict for LinearBurdenTest.run()."""
    return {
        "genotype_matrix": geno,
        "variant_mafs": mafs,
        "phenotype_vector": phenotype,
        "covariate_matrix": covariate_matrix,
        "proband_count": 0,  # N/A for quantitative trait
        "control_count": 0,
        "n_qualifying_variants": geno.shape[1],
    }


# ---------------------------------------------------------------------------
# Manual statsmodels.OLS validation (BURDEN-02)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLinearBurdenMatchesStatsmodels:
    """LinearBurdenTest output matches manual statsmodels.OLS — BURDEN-02."""

    def test_linear_burden_matches_manual_statsmodels(
        self, synthetic_quantitative_data, default_config
    ) -> None:
        """beta, se, p-value match direct sm.OLS call on same inputs (rel=1e-6)."""
        geno, phenotype, mafs = synthetic_quantitative_data
        config = default_config

        from variantcentrifuge.association.weights import get_weights

        weights = get_weights(mafs, config.variant_weights)
        burden = geno @ weights

        # Manual statsmodels reference
        design_manual = sm.add_constant(burden.reshape(-1, 1))
        fit_manual = sm.OLS(phenotype, design_manual).fit()
        manual_p = float(fit_manual.pvalues[1])
        manual_beta = float(fit_manual.params[1])
        manual_se = float(fit_manual.bse[1])

        # LinearBurdenTest
        test = LinearBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, config)

        assert result.p_value is not None, "Expected valid p-value"
        assert result.p_value == pytest.approx(manual_p, rel=1e-6), (
            f"p-value mismatch: got {result.p_value}, expected {manual_p}"
        )
        assert result.effect_size == pytest.approx(manual_beta, rel=1e-6), (
            f"beta mismatch: got {result.effect_size}, expected {manual_beta}"
        )
        assert result.extra["se"] == pytest.approx(manual_se, rel=1e-6), (
            f"se mismatch: got {result.extra['se']}, expected {manual_se}"
        )

    def test_linear_burden_with_covariates(
        self, synthetic_quantitative_data, default_config
    ) -> None:
        """Burden coefficient still extracted at index 1 when covariates present."""
        geno, phenotype, mafs = synthetic_quantitative_data
        config = default_config
        np.random.seed(13)
        n = len(phenotype)

        # 2-column covariate matrix
        covariate_matrix = np.column_stack(
            [
                np.random.normal(50, 10, n),
                np.random.normal(0, 1, n),
            ]
        )

        from variantcentrifuge.association.weights import get_weights

        weights = get_weights(mafs, config.variant_weights)
        burden = geno @ weights

        # Manual: intercept + burden + 2 covariates
        design_manual = sm.add_constant(burden.reshape(-1, 1))
        design_manual = np.column_stack([design_manual, covariate_matrix])
        fit_manual = sm.OLS(phenotype, design_manual).fit()
        # Burden is at index 1
        manual_beta = float(fit_manual.params[1])
        manual_p = float(fit_manual.pvalues[1])

        # LinearBurdenTest with covariates
        test = LinearBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs, covariate_matrix)
        result = test.run("GENE1", contingency, config)

        assert result.p_value is not None
        assert result.effect_size == pytest.approx(manual_beta, rel=1e-6), (
            f"beta mismatch with covariates: got {result.effect_size}, expected {manual_beta}"
        )
        assert result.p_value == pytest.approx(manual_p, rel=1e-6)

    def test_linear_burden_ci_correct(self, synthetic_quantitative_data, default_config) -> None:
        """95% CI matches manual statsmodels.OLS conf_int at index 1."""
        geno, phenotype, mafs = synthetic_quantitative_data
        config = default_config

        from variantcentrifuge.association.weights import get_weights

        weights = get_weights(mafs, config.variant_weights)
        burden = geno @ weights
        design = sm.add_constant(burden.reshape(-1, 1))
        fit = sm.OLS(phenotype, design).fit()
        ci = fit.conf_int()
        expected_ci_lower = float(ci[1, 0])
        expected_ci_upper = float(ci[1, 1])

        test = LinearBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, config)

        assert result.ci_lower is not None and result.ci_upper is not None
        assert result.ci_lower == pytest.approx(expected_ci_lower, rel=1e-6), (
            f"ci_lower mismatch: got {result.ci_lower}, expected {expected_ci_lower}"
        )
        assert result.ci_upper == pytest.approx(expected_ci_upper, rel=1e-6)

    def test_linear_burden_effect_size_is_beta_not_or(
        self, synthetic_quantitative_data, default_config
    ) -> None:
        """effect_size is beta (raw coefficient), not exp(beta) — linear test reports beta."""
        geno, phenotype, mafs = synthetic_quantitative_data

        from variantcentrifuge.association.weights import get_weights

        config = default_config
        weights = get_weights(mafs, config.variant_weights)
        burden = geno @ weights

        # Reference beta from manual OLS
        design = sm.add_constant(burden.reshape(-1, 1))
        fit = sm.OLS(phenotype, design).fit()
        reference_beta = float(fit.params[1])

        test = LinearBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, config)

        assert result.effect_size is not None
        # effect_size should match raw beta (not exp(beta))
        assert result.effect_size == pytest.approx(reference_beta, rel=1e-6), (
            "effect_size for linear test should be raw beta"
        )
        # Verify it's NOT exp(beta) — they should differ for large betas
        exp_beta = float(np.exp(reference_beta))
        assert abs(result.effect_size - exp_beta) > 1.0, (
            "effect_size appears to be exp(beta) rather than raw beta"
        )


# ---------------------------------------------------------------------------
# Warning codes and edge cases
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLinearBurdenWarnings:
    """Warning codes and edge cases for LinearBurdenTest."""

    def test_linear_burden_no_genotype_matrix(self, default_config: AssociationConfig) -> None:
        """Missing genotype_matrix -> p_value=None and NO_GENOTYPE_MATRIX warning."""
        contingency = {
            "proband_count": 0,
            "control_count": 0,
            "n_qualifying_variants": 3,
            # No genotype_matrix key!
        }

        test = LinearBurdenTest()
        result = test.run("FISHER_ONLY_GENE", contingency, default_config)

        assert result.p_value is None, "Expected p_value=None when genotype_matrix absent"
        warnings = result.extra.get("warnings", [])
        assert "NO_GENOTYPE_MATRIX" in warnings, (
            f"Expected NO_GENOTYPE_MATRIX warning, got: {warnings}"
        )

    def test_linear_burden_no_phenotype_returns_none(
        self, default_config: AssociationConfig
    ) -> None:
        """Missing phenotype_vector -> p_value=None."""
        np.random.seed(0)
        n = 30
        geno = np.random.randint(0, 2, size=(n, 3)).astype(float)
        mafs = np.array([0.1, 0.05, 0.02])

        contingency = {
            "genotype_matrix": geno,
            "variant_mafs": mafs,
            "phenotype_vector": None,  # missing
            "proband_count": 0,
            "control_count": 0,
            "n_qualifying_variants": 3,
        }

        test = LinearBurdenTest()
        result = test.run("NO_PHENO_GENE", contingency, default_config)

        assert result.p_value is None

    def test_linear_burden_carrier_count(self, synthetic_quantitative_data, default_config) -> None:
        """n_carriers in extra matches sum(burden > 0)."""
        geno, phenotype, mafs = synthetic_quantitative_data

        from variantcentrifuge.association.weights import get_weights

        config = default_config
        weights = get_weights(mafs, config.variant_weights)
        burden = geno @ weights
        expected_carriers = int((burden > 0).sum())

        test = LinearBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, config)

        assert result.extra["n_carriers"] == expected_carriers, (
            f"n_carriers mismatch: got {result.extra['n_carriers']}, expected {expected_carriers}"
        )

    def test_linear_burden_result_metadata(
        self, synthetic_quantitative_data, default_config
    ) -> None:
        """TestResult has correct gene name and test_name."""
        geno, phenotype, mafs = synthetic_quantitative_data
        test = LinearBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("MY_QUANT_GENE", contingency, default_config)

        assert result.gene == "MY_QUANT_GENE"
        assert result.test_name == "linear_burden"

    def test_linear_burden_check_dependencies(self) -> None:
        """check_dependencies does not raise (statsmodels is available)."""
        test = LinearBurdenTest()
        test.check_dependencies()

    def test_linear_burden_name_property(self) -> None:
        """LinearBurdenTest.name == 'linear_burden'."""
        test = LinearBurdenTest()
        assert test.name == "linear_burden"

    def test_linear_burden_ci_bounds_ordered(
        self, synthetic_quantitative_data, default_config
    ) -> None:
        """ci_lower < ci_upper always."""
        geno, phenotype, mafs = synthetic_quantitative_data
        test = LinearBurdenTest()
        contingency = _make_contingency(geno, phenotype, mafs)
        result = test.run("GENE1", contingency, default_config)

        if result.ci_lower is not None and result.ci_upper is not None:
            assert result.ci_lower < result.ci_upper, (
                f"ci_lower={result.ci_lower} should be < ci_upper={result.ci_upper}"
            )
