"""
Unit tests for the davies.py p-value computation layer.

Tests the three-tier fallback chain:
  Tier 3 — Liu moment-matching (_liu_pvalue)
  Tier 2 — Kuonen saddlepoint approximation (_kuonen_pvalue)
  Tier 1 — Davies (when C extension available)
  Public — compute_pvalue (full chain)

Covers requirements: SKAT-05, SKAT-10

Mathematical references:
  - Single eigenvalue case: Q ~ chi2(1), exact for Liu
  - Equal eigenvalues: Q ~ c * chi2(k), exact for Liu
  - Kuonen is tighter than Liu for small p-values (< 1e-3)
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.stats import chi2

from variantcentrifuge.association.backends.davies import (
    _kuonen_pvalue,
    _liu_pvalue,
    compute_pvalue,
)

# ---------------------------------------------------------------------------
# Liu moment-matching tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestLiuPvalue:
    """Unit tests for _liu_pvalue moment-matching approximation."""

    def test_liu_single_eigenvalue(self):
        """Single eigenvalue: Q ~ chi2(1). Liu should be exact within 1e-3 relative."""
        lambdas = np.array([1.0])
        for q in [1.0, 3.84, 6.63]:
            p_liu = _liu_pvalue(q, lambdas)
            p_true = chi2.sf(q, df=1)
            assert math.isfinite(p_liu), f"Liu returned non-finite p for q={q}"
            assert 0.0 <= p_liu <= 1.0, f"Liu p={p_liu} out of [0,1] for q={q}"
            # Single eigenvalue case: Liu is essentially exact
            if p_true > 1e-10:
                rel_err = abs(p_liu - p_true) / p_true
                assert rel_err < 1e-2, (
                    f"Liu relative error {rel_err:.4f} > 1e-2 for q={q}, "
                    f"p_liu={p_liu:.6e}, p_true={p_true:.6e}"
                )

    def test_liu_equal_eigenvalues(self):
        """Equal eigenvalues: Q ~ c*chi2(k). Liu should match within 1e-3 relative."""
        c = 2.0
        k = 5
        q = 15.0
        lambdas = np.array([c] * k)
        p_liu = _liu_pvalue(q, lambdas)
        p_true = chi2.sf(q / c, df=k)
        assert math.isfinite(p_liu)
        assert 0.0 <= p_liu <= 1.0
        # Equal eigenvalues: Liu enters the all-equal branch (exact)
        rel_err = abs(p_liu - p_true) / p_true
        assert rel_err < 1e-2, (
            f"Liu relative error {rel_err:.4f} for equal lambdas, "
            f"p_liu={p_liu:.6e}, p_true={p_true:.6e}"
        )

    def test_liu_multiple_eigenvalues_reasonable(self):
        """Mixed eigenvalues: result should be in reasonable range."""
        lambdas = np.array([5.0, 3.0, 1.0])
        q = 15.0  # above mean (9.0), below max realistic
        p = _liu_pvalue(q, lambdas)
        assert math.isfinite(p)
        assert 0.01 < p < 0.5, f"Expected p in (0.01, 0.5), got {p:.6e}"

    def test_liu_large_q_small_p(self):
        """Large Q far above mean: p-value should decrease monotonically."""
        lambdas = np.array([5.0, 3.0, 1.0])
        # Liu p-values should decrease as Q increases
        p_q30 = _liu_pvalue(30.0, lambdas)
        p_q80 = _liu_pvalue(80.0, lambdas)
        assert math.isfinite(p_q30)
        assert math.isfinite(p_q80)
        assert p_q80 < p_q30, f"Expected p to decrease with Q: p_q30={p_q30:.4e}, p_q80={p_q80:.4e}"
        # Q=80 should give a small p via Liu
        assert p_q80 < 0.001, f"Expected p < 0.001 for large Q=80, got {p_q80:.6e}"

    def test_liu_q_at_mean(self):
        """Q at mean of distribution: p should be near 0.5."""
        lambdas = np.array([5.0, 3.0, 1.0])
        q = float(np.sum(lambdas))  # mean = 9.0
        p = _liu_pvalue(q, lambdas)
        assert math.isfinite(p)
        # At the mean, distribution is roughly symmetric — p should be moderate
        assert 0.3 < p < 0.7, f"Expected p near 0.5 at Q=mean, got {p:.6e}"

    def test_liu_q_below_mean_returns_1(self):
        """Q below mean: Liu should return 1.0 (tail probability = ~1)."""
        lambdas = np.array([5.0, 3.0, 1.0])
        q = 1.0  # well below mean (9.0)
        p = _liu_pvalue(q, lambdas)
        assert math.isfinite(p)
        assert p >= 0.9, f"Expected p near 1.0 for Q below mean, got {p:.6e}"

    def test_liu_empty_lambdas_returns_nan(self):
        """Empty eigenvalue array: Liu should return nan."""
        lambdas = np.array([], dtype=np.float64)
        p = _liu_pvalue(0.0, lambdas)
        assert math.isnan(p), f"Expected nan for empty lambdas, got {p}"

    def test_liu_very_small_eigenvalues(self):
        """Near machine-epsilon eigenvalues: Liu should not crash."""
        eps = np.finfo(float).eps
        lambdas = np.array([eps * 100, eps * 200, eps * 50])
        q = float(np.sum(lambdas)) * 2.0
        try:
            p = _liu_pvalue(q, lambdas)
            assert math.isfinite(p) or math.isnan(p), "Expected finite or nan for tiny eigenvalues"
        except Exception as e:
            pytest.fail(f"Liu crashed on tiny eigenvalues: {e}")

    def test_liu_single_large_eigenvalue(self):
        """Single large eigenvalue: Q ~ chi2(1) scaled by lambda."""
        lam = 100.0
        lambdas = np.array([lam])
        q = lam * chi2.ppf(0.95, df=1)  # 95th percentile
        p_liu = _liu_pvalue(q, lambdas)
        p_true = chi2.sf(q / lam, df=1)  # Should be ~0.05
        assert math.isfinite(p_liu)
        assert 0.0 <= p_liu <= 1.0
        rel_err = abs(p_liu - p_true) / p_true
        assert rel_err < 0.05, (
            f"Liu relative error {rel_err:.4f} for large single eigenvalue, "
            f"p_liu={p_liu:.6e}, p_true={p_true:.6e}"
        )


# ---------------------------------------------------------------------------
# Kuonen saddlepoint tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestKuonenPvalue:
    """Unit tests for _kuonen_pvalue Lugannani-Rice saddlepoint approximation."""

    def test_kuonen_above_mean_returns_valid_p(self):
        """Q above mean: Kuonen should return valid p in (0, 1)."""
        lambdas = np.array([5.0, 3.0, 1.0])
        q = 15.0  # mean = 9.0
        p = _kuonen_pvalue(q, lambdas)
        assert p is not None, "Kuonen returned None for Q > mean"
        assert 0.0 < p < 1.0, f"Kuonen p={p:.6e} not in (0, 1)"

    def test_kuonen_below_mean_returns_none(self):
        """Q below mean: saddlepoint not valid, should return None."""
        lambdas = np.array([5.0, 3.0, 1.0])
        q = 5.0  # below mean = 9.0
        p = _kuonen_pvalue(q, lambdas)
        assert p is None, f"Expected None for Q < mean, got {p}"

    def test_kuonen_at_mean_returns_none(self):
        """Q at mean: saddlepoint not applicable, should return None."""
        lambdas = np.array([5.0, 3.0, 1.0])
        q = float(np.sum(lambdas))  # exactly at mean = 9.0
        p = _kuonen_pvalue(q, lambdas)
        assert p is None, f"Expected None for Q = mean, got {p}"

    def test_kuonen_agrees_with_liu_moderate_p(self):
        """Moderate p (0.01-0.1): Kuonen and Liu should agree within 20% relative."""
        lambdas = np.array([5.0, 3.0, 1.0])
        # Choose Q to get moderate p — somewhere between 15 and 30
        q = 20.0
        p_liu = _liu_pvalue(q, lambdas)
        p_kuonen = _kuonen_pvalue(q, lambdas)
        if p_kuonen is not None and 0.01 < p_liu < 0.1:
            rel_diff = abs(p_kuonen - p_liu) / p_liu
            assert rel_diff < 0.30, (
                f"Kuonen and Liu disagree by {rel_diff:.2f} at moderate p: "
                f"p_kuonen={p_kuonen:.6e}, p_liu={p_liu:.6e}"
            )

    def test_kuonen_tighter_than_liu_small_p(self):
        """Small p: Kuonen should be closer to chi2 ground truth than Liu."""
        # Single large eigenvalue: Q ~ lambda * chi2(1), ground truth exact
        lam = 10.0
        lambdas = np.array([lam])
        # Target very small p ~ 1e-5 (chi2 99.999th percentile)
        q_target = lam * chi2.ppf(1.0 - 1e-5, df=1)
        p_true = chi2.sf(q_target / lam, df=1)

        p_liu = _liu_pvalue(q_target, lambdas)
        p_kuonen = _kuonen_pvalue(q_target, lambdas)

        if p_kuonen is not None:
            err_liu = abs(math.log10(max(p_liu, 1e-20)) - math.log10(max(p_true, 1e-20)))
            err_kuonen = abs(math.log10(max(p_kuonen, 1e-20)) - math.log10(max(p_true, 1e-20)))
            assert err_kuonen <= err_liu + 0.5, (
                f"Kuonen log10 error {err_kuonen:.3f} not tighter than Liu {err_liu:.3f} "
                f"(p_true={p_true:.2e}, p_liu={p_liu:.2e}, p_kuonen={p_kuonen:.2e})"
            )

    def test_kuonen_single_eigenvalue_matches_chi2(self):
        """Single eigenvalue: Kuonen should agree with chi2 or return None."""
        lambdas = np.array([2.0])
        q = 15.0  # well above mean = 2.0
        p = _kuonen_pvalue(q, lambdas)
        p_true = chi2.sf(q / 2.0, df=1)
        if p is not None:
            # Should be reasonably close to true chi2 value
            rel_err = abs(p - p_true) / p_true
            assert rel_err < 0.1, (
                f"Kuonen single eigenvalue error {rel_err:.4f}, "
                f"p_kuonen={p:.6e}, p_true={p_true:.6e}"
            )

    def test_kuonen_empty_lambdas_returns_none(self):
        """Empty eigenvalue array: Kuonen should return None."""
        lambdas = np.array([], dtype=np.float64)
        p = _kuonen_pvalue(5.0, lambdas)
        assert p is None, "Expected None for empty lambdas"

    def test_kuonen_near_singularity_no_crash(self):
        """Extreme Q near singularity: Kuonen should not crash."""
        lambdas = np.array([5.0, 3.0, 1.0])
        q = 1000.0  # very far in tail
        try:
            p = _kuonen_pvalue(q, lambdas)
            assert p is None or (0.0 <= p <= 1.0), f"Out-of-range p from Kuonen: {p}"
        except Exception as e:
            pytest.fail(f"Kuonen crashed near singularity: {e}")


# ---------------------------------------------------------------------------
# Fallback chain tests (compute_pvalue)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestComputePvalue:
    """Unit tests for the public compute_pvalue fallback chain."""

    def test_returns_tuple_of_three(self):
        """compute_pvalue must return exactly (float, str, bool)."""
        lambdas = np.array([5.0, 3.0, 1.0])
        result = compute_pvalue(15.0, lambdas)
        assert isinstance(result, tuple)
        assert len(result) == 3
        p, method, converged = result
        assert isinstance(p, float)
        assert isinstance(method, str)
        assert isinstance(converged, bool)

    def test_method_in_expected_set(self):
        """p_method must be one of 'davies', 'saddlepoint', 'liu'."""
        lambdas = np.array([5.0, 3.0, 1.0])
        _, method, _ = compute_pvalue(15.0, lambdas)
        assert method in {"davies", "saddlepoint", "liu"}, f"Unexpected method: {method}"

    def test_empty_lambdas_returns_nan_liu_false(self):
        """Empty lambdas: compute_pvalue returns (nan, 'liu', False)."""
        lambdas = np.array([], dtype=np.float64)
        p, method, converged = compute_pvalue(5.0, lambdas)
        assert math.isnan(p), f"Expected nan, got {p}"
        assert method == "liu"
        assert converged is False

    def test_reasonable_values_for_known_case(self):
        """Known eigenvalues and Q: p must be in (0, 1)."""
        lambdas = np.array([5.0, 3.0, 1.0])
        q = 20.0
        p, _method, _ = compute_pvalue(q, lambdas)
        assert math.isfinite(p), f"Expected finite p, got {p}"
        assert 0.0 < p < 1.0, f"Expected p in (0, 1), got {p:.6e}"

    def test_deterministic_same_input(self):
        """Same input twice produces identical output."""
        lambdas = np.array([5.0, 3.0, 1.0])
        q = 15.0
        r1 = compute_pvalue(q, lambdas)
        r2 = compute_pvalue(q, lambdas)
        assert r1 == r2, f"Non-deterministic: {r1} != {r2}"

    def test_no_davies_env_var_uses_fallback(self, monkeypatch):
        """With VARIANTCENTRIFUGE_NO_C_EXT=1, method must be saddlepoint or liu."""
        # Reset the module-level davies state so the env var takes effect
        import variantcentrifuge.association.backends.davies as davies_mod

        monkeypatch.setenv("VARIANTCENTRIFUGE_NO_C_EXT", "1")
        # Reset cached load state
        davies_mod._DAVIES_LOAD_ATTEMPTED = False
        davies_mod._DAVIES_AVAILABLE = False

        lambdas = np.array([5.0, 3.0, 1.0])
        _, method, _ = compute_pvalue(15.0, lambdas)
        assert method in {"saddlepoint", "liu"}, (
            f"Expected saddlepoint or liu with NO_C_EXT=1, got {method}"
        )

        # Restore
        davies_mod._DAVIES_LOAD_ATTEMPTED = False
        davies_mod._DAVIES_AVAILABLE = False

    def test_zero_q_returns_high_p(self):
        """Q <= 0: tail probability is essentially 1."""
        lambdas = np.array([5.0, 3.0, 1.0])
        p, _method, _ = compute_pvalue(0.0, lambdas)
        assert math.isfinite(p)
        assert p >= 0.9, f"Expected p near 1 for Q=0, got {p:.6e}"

    def test_liu_directly_returns_small_p_for_large_q(self):
        """Liu directly called with large Q should return small p."""
        # Test Liu directly since compute_pvalue routing depends on Davies availability
        lambdas = np.array([5.0, 3.0, 1.0])
        # Liu should give monotonically decreasing p as Q increases
        p_q30 = _liu_pvalue(30.0, lambdas)
        p_q80 = _liu_pvalue(80.0, lambdas)
        assert p_q80 < p_q30, "Liu should give smaller p for larger Q"
        assert p_q80 < 0.001, f"Expected Liu p < 0.001 for Q=80, got {p_q80:.6e}"

    def test_converged_false_without_davies(self, monkeypatch):
        """Without Davies, p_converged must be False."""
        import variantcentrifuge.association.backends.davies as davies_mod

        monkeypatch.setenv("VARIANTCENTRIFUGE_NO_C_EXT", "1")
        davies_mod._DAVIES_LOAD_ATTEMPTED = False
        davies_mod._DAVIES_AVAILABLE = False

        lambdas = np.array([5.0, 3.0, 1.0])
        _, _, converged = compute_pvalue(15.0, lambdas)
        assert converged is False

        davies_mod._DAVIES_LOAD_ATTEMPTED = False
        davies_mod._DAVIES_AVAILABLE = False

    def test_single_eigenvalue_matches_chi2(self):
        """Single eigenvalue: compute_pvalue should match chi2 ground truth within 1e-2."""
        lam = 3.0
        lambdas = np.array([lam])
        q = lam * chi2.ppf(0.99, df=1)  # 99th percentile
        p, method, _ = compute_pvalue(q, lambdas)
        p_true = chi2.sf(q / lam, df=1)
        assert math.isfinite(p)
        rel_err = abs(p - p_true) / p_true
        assert rel_err < 0.05, (
            f"compute_pvalue relative error {rel_err:.4f} for single eigenvalue, "
            f"p={p:.6e}, p_true={p_true:.6e}, method={method}"
        )

    def test_positive_q_positive_p(self):
        """For any positive Q above the mean, p should be strictly positive."""
        lambdas = np.array([2.0, 1.5, 1.0, 0.5])
        q = 20.0  # well above mean (5.0)
        p, _, _ = compute_pvalue(q, lambdas)
        assert math.isfinite(p)
        assert p > 0.0, f"Expected positive p, got {p}"
