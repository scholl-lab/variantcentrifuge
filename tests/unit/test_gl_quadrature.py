"""
Unit tests for the SKAT-O omnibus integration (piecewise adaptive quadrature).

Covers:
- SKAT-O integration accuracy (smoke test: valid p-value from SKAT-O integration)
- parallel_safe attributes on all 6 Python/R backend test classes

Requirements verified: Plan 27-01 Task 2
"""

from __future__ import annotations

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# GL accuracy smoke test
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestGLQuadratureAccuracy:
    """Smoke test that GL quadrature produces valid SKAT-O p-values."""

    def test_gl_quadrature_accuracy(self) -> None:
        """
        _skato_integrate_davies (using GL) produces a valid p-value in (0, 1].

        Creates a synthetic SKAT-O scenario (n=200, p=20 binary trait) using
        PythonSKATBackend and verifies the omnibus p-value is a valid float in
        the range (0, 1]. This is a smoke test ensuring GL integration produces
        sensible results comparable to the previous adaptive quad approach.
        """
        from variantcentrifuge.association.backends.python_backend import PythonSKATBackend
        from variantcentrifuge.association.weights import beta_maf_weights

        rng = np.random.default_rng(42)
        n, p = 200, 20

        # Synthetic binary phenotype: ~50% prevalence
        phenotype = rng.integers(0, 2, n).astype(np.float64)

        # Synthetic genotype matrix: rare variants (MAF ~1-5%)
        mafs = rng.uniform(0.01, 0.05, p)
        geno = np.zeros((n, p), dtype=np.float64)
        for j in range(p):
            geno[:, j] = rng.binomial(2, mafs[j], n).astype(np.float64)

        # Ensure rank >= 2
        assert np.linalg.matrix_rank(geno) >= 2, "Test genotype matrix has rank < 2"

        backend = PythonSKATBackend()
        backend.detect_environment()

        null_model = backend.fit_null_model(phenotype, None, "binary")
        weights = beta_maf_weights(geno.mean(axis=0) / 2.0, a=1.0, b=25.0)
        result = backend.test_gene(
            gene="TEST_GENE",
            genotype_matrix=geno,
            null_model=null_model,
            method="SKATO",
            weights=weights,
            weights_beta=(1.0, 25.0),
        )

        p_value = result.get("p_value")
        assert p_value is not None, "SKAT-O returned p_value=None for valid input"
        assert isinstance(p_value, float), f"p_value is not float: {type(p_value)}"
        assert 0.0 < p_value <= 1.0, f"p_value={p_value} out of range (0, 1]"

    def test_skato_integral_matches_adaptive_quadrature(self) -> None:
        """Regression: a single 128-node rule over x was ~5% off R on null genes."""
        import scipy.integrate
        import scipy.stats

        from variantcentrifuge.association.backends import skat_math as sm
        from variantcentrifuge.association.backends.davies import davies_pvalue

        rng = np.random.default_rng(7)
        n, p = 700, 12
        geno = rng.binomial(2, rng.uniform(0.003, 0.03, p), size=(n, p)).astype(np.float64)
        res = rng.integers(0, 2, n) - 0.5
        z1_half = (geno - geno.mean(axis=0)) * 0.5 / np.sqrt(2.0)
        score = geno.T @ res
        rho = [0.0, 0.01, 0.04, 0.09, 0.25, 0.5, 0.999]
        q = np.array([((1 - r) * score @ score + r * score.sum() ** 2) / 8.0 for r in rho])
        param = sm._skato_optimal_param(z1_half, rho)
        a_mat = z1_half.T @ z1_half
        lam = []
        for r in rho:
            s = np.sqrt(1 - r)
            d = (np.sqrt(1 - r + p * r) - s) / p
            k = s * s * a_mat + s * d * (a_mat.sum(0)[None, :] + a_mat.sum(1)[:, None])
            lam.append(sm._get_lambda(k + d * d * a_mat.sum()))
        pmin, _pval, pmin_q = sm._skato_each_q(param, q[None, :], rho, lam)
        got = sm._skato_integrate_davies(pmin_q[0], param, rho, float(pmin[0]))

        r_arr = np.array(rho)
        sd1 = np.sqrt((param["var_q"] - param["var_remain"]) / param["var_q"])

        def f(x: float) -> float:
            m = float(np.min((pmin_q[0] - param["tau"] * x) / (1 - r_arr)))
            c = davies_pvalue((m - param["mu_q"]) * sd1 + param["mu_q"], param["lambdas"])[0]
            return (1 - min(c, 1.0)) * float(scipy.stats.chi2.pdf(x, 1))

        ref = 1 - scipy.integrate.quad(f, 0, 40, limit=1000, epsabs=1e-25)[0]
        assert 1e-4 < ref < 0.9
        assert got == pytest.approx(min(ref, float(pmin[0]) * len(rho)), rel=2e-3)


# ---------------------------------------------------------------------------
# parallel_safe attribute tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestParallelSafeAttributes:
    """Tests for parallel_safe class attributes on all 6 test classes."""

    def test_fisher_exact_test_parallel_safe_true(self) -> None:
        """FisherExactTest.parallel_safe is True (pure scipy, thread-safe)."""
        from variantcentrifuge.association.tests.fisher import FisherExactTest

        assert FisherExactTest.parallel_safe is True, "FisherExactTest.parallel_safe should be True"

    def test_logistic_burden_test_parallel_safe_true(self) -> None:
        """LogisticBurdenTest.parallel_safe is True (pure statsmodels, thread-safe)."""
        from variantcentrifuge.association.tests.logistic_burden import LogisticBurdenTest

        assert LogisticBurdenTest.parallel_safe is True, (
            "LogisticBurdenTest.parallel_safe should be True"
        )

    def test_linear_burden_test_parallel_safe_true(self) -> None:
        """LinearBurdenTest.parallel_safe is True (pure statsmodels, thread-safe)."""
        from variantcentrifuge.association.tests.linear_burden import LinearBurdenTest

        assert LinearBurdenTest.parallel_safe is True, (
            "LinearBurdenTest.parallel_safe should be True"
        )

    def test_pure_python_skat_test_parallel_safe_true(self) -> None:
        """PurePythonSKATTest.parallel_safe is True (pure numpy/scipy, thread-safe)."""
        from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest

        assert PurePythonSKATTest.parallel_safe is True, (
            "PurePythonSKATTest.parallel_safe should be True"
        )

    def test_pure_python_coast_test_parallel_safe_true(self) -> None:
        """PurePythonCOASTTest.parallel_safe is True (no rpy2 dependency)."""
        from variantcentrifuge.association.tests.allelic_series_python import PurePythonCOASTTest

        assert PurePythonCOASTTest.parallel_safe is True, (
            "PurePythonCOASTTest.parallel_safe should be True"
        )

    def test_coast_test_parallel_safe_false(self) -> None:
        """COASTTest.parallel_safe is False (rpy2 restriction: main thread only)."""
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            from variantcentrifuge.association.tests.allelic_series import COASTTest

        assert COASTTest.parallel_safe is False, (
            "COASTTest.parallel_safe should be False (rpy2 main-thread restriction)"
        )

    def test_parallel_safe_is_class_attribute_not_instance(self) -> None:
        """parallel_safe is accessible as a class attribute (not just on instances)."""
        from variantcentrifuge.association.tests.fisher import FisherExactTest
        from variantcentrifuge.association.tests.linear_burden import LinearBurdenTest
        from variantcentrifuge.association.tests.logistic_burden import LogisticBurdenTest
        from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest

        # Access via class directly (not via instantiation)
        assert hasattr(FisherExactTest, "parallel_safe")
        assert hasattr(LogisticBurdenTest, "parallel_safe")
        assert hasattr(LinearBurdenTest, "parallel_safe")
        assert hasattr(PurePythonSKATTest, "parallel_safe")
