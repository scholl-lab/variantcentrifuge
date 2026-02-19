"""
Unit tests for association multiple testing correction module.

Verifies that apply_correction() produces bit-identical output to direct
statsmodels.stats.multitest.multipletests() calls for the same inputs.
"""

from __future__ import annotations

import numpy as np
import pytest
import statsmodels.stats.multitest as smm

from variantcentrifuge.association.correction import apply_correction


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _smm_fdr(pvals):
    """Direct smm FDR call — the ground truth reference."""
    return smm.multipletests(np.asarray(pvals, dtype=float), method="fdr_bh")[1]


def _smm_bonferroni(pvals):
    """Direct smm Bonferroni call — the ground truth reference."""
    return smm.multipletests(np.asarray(pvals, dtype=float), method="bonferroni")[1]


# ---------------------------------------------------------------------------
# FDR parity tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestApplyCorrectionFDR:
    """apply_correction(method='fdr') matches direct smm.multipletests(method='fdr_bh')."""

    def test_fdr_parity_typical_values(self):
        """FDR correction matches smm for typical p-values."""
        pvals = [0.01, 0.05, 0.1, 0.2, 0.5, 0.9]
        expected = _smm_fdr(pvals)
        result = apply_correction(pvals, method="fdr")
        np.testing.assert_array_almost_equal(result, expected, decimal=15)

    def test_fdr_parity_small_p_values(self):
        """FDR correction matches smm for very small p-values (1e-10)."""
        pvals = [1e-10, 1e-8, 1e-6, 0.001, 0.01]
        expected = _smm_fdr(pvals)
        result = apply_correction(pvals, method="fdr")
        np.testing.assert_array_almost_equal(result, expected, decimal=15)

    def test_fdr_parity_all_ones(self):
        """FDR correction with all p-values=1.0 matches smm."""
        pvals = [1.0, 1.0, 1.0, 1.0]
        expected = _smm_fdr(pvals)
        result = apply_correction(pvals, method="fdr")
        np.testing.assert_array_almost_equal(result, expected, decimal=15)

    def test_fdr_parity_single_value(self):
        """FDR correction with single p-value matches smm."""
        pvals = [0.04]
        expected = _smm_fdr(pvals)
        result = apply_correction(pvals, method="fdr")
        np.testing.assert_array_almost_equal(result, expected, decimal=15)

    def test_fdr_default_method_is_fdr(self):
        """Default method parameter is 'fdr' (Benjamini-Hochberg)."""
        pvals = [0.01, 0.05, 0.2, 0.5]
        expected = apply_correction(pvals, method="fdr")
        default = apply_correction(pvals)  # no method arg
        np.testing.assert_array_equal(default, expected)

    def test_fdr_unknown_method_treated_as_fdr(self):
        """Unknown method falls back to FDR (matches spec: any unknown -> fdr)."""
        pvals = [0.01, 0.05, 0.2]
        fdr_result = apply_correction(pvals, method="fdr")
        unknown_result = apply_correction(pvals, method="unknown_method")
        np.testing.assert_array_equal(unknown_result, fdr_result)


# ---------------------------------------------------------------------------
# Bonferroni parity tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestApplyCorrectionBonferroni:
    """apply_correction(method='bonferroni') matches direct smm.multipletests(method='bonferroni')."""

    def test_bonferroni_parity_typical_values(self):
        """Bonferroni correction matches smm for typical p-values."""
        pvals = [0.01, 0.05, 0.1, 0.5, 0.9]
        expected = _smm_bonferroni(pvals)
        result = apply_correction(pvals, method="bonferroni")
        np.testing.assert_array_almost_equal(result, expected, decimal=15)

    def test_bonferroni_parity_small_p_values(self):
        """Bonferroni correction matches smm for very small p-values."""
        pvals = [1e-10, 1e-8, 1e-5, 0.001]
        expected = _smm_bonferroni(pvals)
        result = apply_correction(pvals, method="bonferroni")
        np.testing.assert_array_almost_equal(result, expected, decimal=15)

    def test_bonferroni_parity_all_ones(self):
        """Bonferroni with all p=1.0 returns all 1.0 (capped at 1)."""
        pvals = [1.0, 1.0, 1.0]
        expected = _smm_bonferroni(pvals)
        result = apply_correction(pvals, method="bonferroni")
        np.testing.assert_array_almost_equal(result, expected, decimal=15)

    def test_bonferroni_parity_single_value(self):
        """Bonferroni with single p-value: p * 1 = p (capped at 1.0)."""
        pvals = [0.03]
        expected = _smm_bonferroni(pvals)
        result = apply_correction(pvals, method="bonferroni")
        np.testing.assert_array_almost_equal(result, expected, decimal=15)


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestApplyCorrectionEdgeCases:
    """Edge cases for apply_correction."""

    def test_empty_list_returns_empty_array(self):
        """Empty input returns empty numpy array."""
        result = apply_correction([])
        assert isinstance(result, np.ndarray)
        assert len(result) == 0

    def test_numpy_array_input_accepted(self):
        """apply_correction accepts numpy array input as well as list."""
        pvals_list = [0.01, 0.05, 0.2]
        pvals_np = np.array(pvals_list)
        result_list = apply_correction(pvals_list)
        result_np = apply_correction(pvals_np)
        np.testing.assert_array_equal(result_list, result_np)

    def test_returns_numpy_array(self):
        """Return type is always np.ndarray."""
        result = apply_correction([0.1, 0.2])
        assert isinstance(result, np.ndarray)

    def test_order_preserved(self):
        """Output p-values maintain input order (not sorted)."""
        # With FDR BH, result depends on rank; verify order is preserved by
        # checking the smm output directly matches apply_correction order
        pvals = [0.5, 0.01, 0.2, 0.001]
        expected = _smm_fdr(pvals)
        result = apply_correction(pvals, method="fdr")
        np.testing.assert_array_equal(result, expected)

    def test_large_input_matches_smm(self):
        """apply_correction handles 100+ p-values and matches smm exactly."""
        rng = np.random.default_rng(42)
        pvals = rng.uniform(0, 1, 100).tolist()
        expected = _smm_fdr(pvals)
        result = apply_correction(pvals, method="fdr")
        np.testing.assert_array_almost_equal(result, expected, decimal=15)
