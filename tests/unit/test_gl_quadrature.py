"""
Unit tests for 128-node Gauss-Legendre quadrature in SKAT-O integration.

Covers:
- GL constants shape and domain validation
- GL quadrature accuracy (smoke test: valid p-value from SKAT-O integration)
- parallel_safe attributes on all 6 Python/R backend test classes

Requirements verified: Plan 27-01 Task 2
"""

from __future__ import annotations

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# GL constants tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestGLConstants:
    """Tests for module-level GL quadrature constants."""

    def test_gl_constants_shape(self) -> None:
        """_GL_X and _GL_W each have exactly 128 elements."""
        from variantcentrifuge.association.backends.python_backend import _GL_W, _GL_X

        assert len(_GL_X) == 128, f"Expected 128 GL nodes, got {len(_GL_X)}"
        assert len(_GL_W) == 128, f"Expected 128 GL weights, got {len(_GL_W)}"

    def test_gl_nodes_in_bounds(self) -> None:
        """All GL nodes fall within [0, 40] (integration bounds matching R upper=40)."""
        from variantcentrifuge.association.backends.python_backend import _GL_X

        assert float(np.min(_GL_X)) >= 0.0, f"GL nodes below 0: min={np.min(_GL_X)}"
        assert float(np.max(_GL_X)) <= 40.0, f"GL nodes above 40: max={np.max(_GL_X)}"

    def test_gl_weights_positive(self) -> None:
        """All GL weights are strictly positive (property of GL quadrature)."""
        from variantcentrifuge.association.backends.python_backend import _GL_W

        assert np.all(_GL_W > 0), "Some GL weights are non-positive"

    def test_gl_weights_sum_matches_interval(self) -> None:
        """Sum of GL weights equals interval length 40 = (b - a)."""
        from variantcentrifuge.association.backends.python_backend import _GL_W

        # For a [a, b] transformed GL rule: sum(w_k) = b - a
        assert abs(float(np.sum(_GL_W)) - 40.0) < 1e-10, (
            f"Sum of GL weights {np.sum(_GL_W):.6f} != 40.0"
        )


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
        result = backend.test_gene(
            gene="TEST_GENE",
            genotype_matrix=geno,
            null_model=null_model,
            method="SKATO",
            weights_beta=(1.0, 25.0),
        )

        p_value = result.get("p_value")
        assert p_value is not None, "SKAT-O returned p_value=None for valid input"
        assert isinstance(p_value, float), f"p_value is not float: {type(p_value)}"
        assert 0.0 < p_value <= 1.0, f"p_value={p_value} out of range (0, 1]"

    def test_gl_vs_known_smooth_integral(self) -> None:
        """
        GL quadrature accurately integrates a smooth function over [0, 40].

        Integrates exp(-x/20) from 0 to 40, which has the exact analytical value
        20*(1 - exp(-2)) = 17.29329... This validates the GL constants independently
        of SKAT-O logic using a function without singularities that GL handles exactly.
        Note: chi2(1) has a singularity at x=0 (pdf->inf) so is not used here;
        the SKAT-O integrand's (1-cdf) factor cancels that singularity in practice.
        """
        import math

        from variantcentrifuge.association.backends.python_backend import _GL_W, _GL_X

        # Integrate exp(-x/20) from 0 to 40 using GL quadrature
        f_vals = np.array([math.exp(-float(x) / 20.0) for x in _GL_X])
        gl_integral = float(np.dot(f_vals, _GL_W))

        # Exact value: integral(exp(-x/20), 0, 40) = 20*(1 - exp(-2))
        true_val = 20.0 * (1.0 - math.exp(-2.0))  # ~17.29329...

        assert abs(gl_integral - true_val) < 1e-10, (
            f"GL integral of exp(-x/20) over [0,40] = {gl_integral:.14f}, "
            f"expected {true_val:.14f}, diff={abs(gl_integral - true_val):.2e}"
        )


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
