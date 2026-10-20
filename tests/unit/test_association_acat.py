"""
Unit tests for ACAT-O Cauchy combination formula.

Tests numerical accuracy against published values (Liu & Xie 2020), edge cases,
and compute_acat_o() integration. All tests use only the public API from
variantcentrifuge.association.tests.acat.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from variantcentrifuge.association.tests.acat import cauchy_combination, compute_acat_o

# ---------------------------------------------------------------------------
# Numerical accuracy tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestCauchyCombinationNumerical:
    """Numerical accuracy tests for the Cauchy combination formula."""

    def test_cauchy_equal_inputs_001(self):
        """cauchy_combination([p, p, p]) should return approximately p for p=0.01."""
        result = cauchy_combination([0.01, 0.01, 0.01])
        assert result == pytest.approx(0.01, rel=1e-4)

    def test_cauchy_equal_inputs_005(self):
        """cauchy_combination([p, p, p]) should return approximately p for p=0.05."""
        result = cauchy_combination([0.05, 0.05, 0.05])
        assert result == pytest.approx(0.05, rel=1e-4)

    def test_cauchy_equal_inputs_01(self):
        """cauchy_combination([p, p, p]) should return approximately p for p=0.1."""
        result = cauchy_combination([0.1, 0.1, 0.1])
        assert result == pytest.approx(0.1, rel=1e-4)

    def test_cauchy_equal_inputs_05(self):
        """cauchy_combination([p, p, p]) should return approximately p for p=0.5."""
        result = cauchy_combination([0.5, 0.5, 0.5])
        assert result == pytest.approx(0.5, rel=1e-4)

    def test_cauchy_published_values(self):
        """cauchy_combination([0.001, 0.01, 0.05]) matches Liu & Xie 2020 Table 1.

        Published reference value: ~0.00268 from the ACAT paper.
        This is the key validation that confirms our Cauchy formula matches the
        original publication.
        """
        result = cauchy_combination([0.001, 0.01, 0.05])
        assert result is not None
        assert result == pytest.approx(0.00268, abs=0.001)

    def test_cauchy_null_inputs(self):
        """cauchy_combination([0.5, 0.5]) should return approximately 0.5 (null p-values)."""
        result = cauchy_combination([0.5, 0.5])
        assert result is not None
        assert result == pytest.approx(0.5, abs=0.01)

    def test_cauchy_very_small_p_no_nan_or_inf(self):
        """cauchy_combination([1e-20, 0.5]) returns a valid float, not NaN or Inf.

        Tests the 1/(p*pi) branch for tiny p-values (Liu & Xie 2020 Section 2.2).
        The result must be a finite number in [0, 1].
        """
        result = cauchy_combination([1e-20, 0.5])
        assert result is not None
        assert not math.isnan(result)
        assert not math.isinf(result)
        assert 0.0 <= result <= 1.0

    def test_cauchy_very_small_p_is_small(self):
        """cauchy_combination([1e-20, 0.5]) should yield a small combined p-value.

        One test has a very significant result (p=1e-20); combining with p=0.5
        should drive the combined p-value well below 0.5.
        """
        result = cauchy_combination([1e-20, 0.5])
        assert result is not None
        assert result < 0.1

    def test_cauchy_weighted_closer_to_small_p(self):
        """cauchy_combination([0.01, 0.5], weights=[0.9, 0.1]) is closer to 0.01 than equal weight.

        With 90% weight on p=0.01 and 10% on p=0.5, the combined p should be
        smaller (more significant) than equal-weighted combination.
        """
        weighted = cauchy_combination([0.01, 0.5], weights=[0.9, 0.1])
        equal_weight = cauchy_combination([0.01, 0.5])

        assert weighted is not None
        assert equal_weight is not None
        # Weighted result should be more significant (smaller p-value) since
        # we put 90% of weight on the more significant input
        assert weighted < equal_weight

    def test_cauchy_combined_between_inputs(self):
        """Combining p=0.01 and p=0.05 yields a combined p more significant than the worst.

        For equal weights, the Cauchy combination of [0.01, 0.05] should yield a value
        smaller than the larger input (0.05) but not necessarily smaller than the best input.
        """
        result = cauchy_combination([0.01, 0.05])
        assert result is not None
        # Combined is more significant (smaller) than the least significant input
        assert result < 0.05
        # Combined is a valid p-value in [0, 1]
        assert 0.0 <= result <= 1.0

    def test_cauchy_result_is_valid_p_value(self):
        """All valid inputs should produce a result in [0, 1]."""
        result = cauchy_combination([0.001, 0.01, 0.05, 0.1, 0.5])
        assert result is not None
        assert 0.0 <= result <= 1.0


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestCauchyCombinationEdgeCases:
    """Edge case tests for cauchy_combination."""

    def test_cauchy_all_none(self):
        """cauchy_combination([None, None]) returns None (no valid p-values)."""
        result = cauchy_combination([None, None])
        assert result is None

    def test_cauchy_single_value(self):
        """cauchy_combination([0.05]) returns 0.05 pass-through (k=1 case).

        Single valid p-value is informative enough to surface without modification
        (CONTEXT.md decision: k=1 pass-through, per IMPL-35).
        """
        result = cauchy_combination([0.05])
        assert result == pytest.approx(0.05)

    def test_cauchy_single_value_identity(self):
        """Pass-through for single value preserves exact float value."""
        result = cauchy_combination([0.123456789])
        assert result == pytest.approx(0.123456789, rel=1e-9)

    def test_cauchy_with_nones_filters_correctly(self):
        """cauchy_combination([0.01, None, 0.05]) combines only the two non-None values."""
        result_with_none = cauchy_combination([0.01, None, 0.05])
        result_without_none = cauchy_combination([0.01, 0.05])
        assert result_with_none is not None
        assert result_with_none == pytest.approx(result_without_none, rel=1e-9)

    def test_cauchy_with_all_nones_except_one(self):
        """cauchy_combination([None, 0.05, None]) returns pass-through 0.05."""
        result = cauchy_combination([None, 0.05, None])
        assert result == pytest.approx(0.05)

    def test_cauchy_all_ones_returns_one(self):
        """cauchy_combination([1.0, 1.0]) is 1.0, as R ACAT() (#216; was None)."""
        assert cauchy_combination([1.0, 1.0]) == 1.0

    def test_cauchy_empty(self):
        """cauchy_combination([]) returns None (empty input)."""
        result = cauchy_combination([])
        assert result is None

    def test_cauchy_with_nan_filtered(self):
        """NaN values are filtered the same as None values."""
        result_with_nan = cauchy_combination([0.01, float("nan"), 0.05])
        result_clean = cauchy_combination([0.01, 0.05])
        assert result_with_nan is not None
        assert not math.isnan(result_with_nan)
        assert result_with_nan == pytest.approx(result_clean, rel=1e-9)

    def test_cauchy_weights_length_mismatch_raises_value_error(self):
        """weights of wrong length raises ValueError."""
        with pytest.raises(ValueError, match="weights length"):
            cauchy_combination([0.01, 0.05], weights=[0.5])


# ---------------------------------------------------------------------------
# compute_acat_o tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestComputeAcatO:
    """Tests for compute_acat_o() which wraps cauchy_combination for named test results."""

    def test_acat_o_two_tests(self):
        """compute_acat_o with two tests returns a valid combined p-value."""
        result = compute_acat_o({"burden": 0.01, "skat": 0.05})
        assert result is not None
        assert 0.0 <= result <= 1.0
        # Combined should be more significant than either input
        assert result < 0.05

    def test_acat_o_single_test_passthrough(self):
        """compute_acat_o with one valid test returns that p-value as pass-through."""
        result = compute_acat_o({"burden": 0.05, "skat": None})
        assert result == pytest.approx(0.05)

    def test_acat_o_no_tests_returns_none(self):
        """compute_acat_o with all None returns None."""
        result = compute_acat_o({"burden": None, "skat": None})
        assert result is None

    def test_acat_o_three_tests(self):
        """compute_acat_o with three tests returns a combined p-value."""
        result = compute_acat_o({"fisher": 0.01, "burden": 0.05, "skat": 0.1})
        assert result is not None
        assert 0.0 <= result <= 1.0
        # Should be more significant than the worst input
        assert result < 0.1

    def test_acat_o_empty_dict_returns_none(self):
        """compute_acat_o with empty dict returns None."""
        result = compute_acat_o({})
        assert result is None

    def test_acat_o_all_ones_returns_one(self):
        """compute_acat_o with all p=1.0 returns 1.0, as R ACAT() (#216)."""
        assert compute_acat_o({"burden": 1.0, "skat": 1.0}) == 1.0

    def test_acat_o_matches_cauchy_combination(self):
        """compute_acat_o should produce the same result as cauchy_combination."""
        test_pvals = {"fisher": 0.01, "burden": 0.05}
        acat_result = compute_acat_o(test_pvals)
        direct_result = cauchy_combination(list(test_pvals.values()))
        assert acat_result == pytest.approx(direct_result, rel=1e-9)


# ---------------------------------------------------------------------------
# Extreme p-value underflow tests
# ---------------------------------------------------------------------------


class TestCauchyCombinationUnderflow:
    """Verify cauchy_combination handles extreme p-values without NaN/inf."""

    def test_all_components_at_1e_300(self):
        """All p-values at 1e-300 should produce a finite, very small result."""
        result = cauchy_combination([1e-300, 1e-300, 1e-300])
        assert result is not None
        assert np.isfinite(result)
        assert 0.0 <= result <= 1.0
        assert result < 1e-10  # must be very significant

    def test_mix_of_extreme_and_moderate(self):
        """Mix of 1e-300 and 0.5 should produce a finite small result."""
        result = cauchy_combination([1e-300, 0.5, 0.5])
        assert result is not None
        assert np.isfinite(result)
        assert 0.0 <= result <= 1.0
        assert result < 0.1

    def test_p_zero_exact_does_not_crash(self):
        """p=0.0 exactly should not produce NaN — clamped to smallest float."""
        result = cauchy_combination([0.0, 0.5])
        assert result is not None
        assert np.isfinite(result)
        assert 0.0 <= result <= 1.0

    def test_p_zero_all_does_not_crash(self):
        """All p=0.0 should not produce NaN."""
        result = cauchy_combination([0.0, 0.0, 0.0])
        assert result is not None
        assert np.isfinite(result)
        assert 0.0 <= result <= 1.0

    def test_weighted_extreme_p_values(self):
        """Weighted Cauchy with extreme p-values should be finite."""
        result = cauchy_combination(
            [1e-300, 1e-200, 0.01, 0.5],
            weights=[1.0, 1.0, 1.0, 1.0],
        )
        assert result is not None
        assert np.isfinite(result)
        assert 0.0 <= result <= 1.0

    def test_coast_7_component_extreme(self):
        """Simulate COAST 7-component omnibus with some extreme p-values."""
        # Weights matching COAST: [1,1,1,1,1,1,6] (50% burden, 50% SKAT)
        weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 6.0]
        p_values = [0.01, 0.05, 0.001, 0.1, 0.005, 0.02, 1e-200]
        result = cauchy_combination(p_values, weights=weights)
        assert result is not None
        assert np.isfinite(result)
        assert 0.0 <= result <= 1.0
        assert result < 0.01  # SKAT at 1e-200 should dominate


# ---------------------------------------------------------------------------
# R ACAT() parity (#216)
# ---------------------------------------------------------------------------

# R 4.5.3, `conda run -n metaskat-check Rscript scripts/generate_acat_parity.R`, ACAT() copied
# verbatim from yaowuliu/ACAT R/ACAT.R @ e737d7d (Liu et al. 2019); "# weighted" section.
# Output: _ledger/evidence23/C1.md. Weights 1 x 6, 6 are the COAST omnibus weights.
_W7 = (1, 1, 1, 1, 1, 1, 6)
R_ACAT_WEIGHTED = [
    ((0.01, 0.05, 0.001, 0.1, 0.005, 0.02, 1e-200), _W7, 2e-200),
    ((0.01, 1.0, 0.3, 1.0, 0.5, 0.2, 0.04), _W7, 0.99999999999999989),
    ((1e-6, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), _W7, 1.0),
    ((0.0, 0.5), (1, 1), 0.0),
    ((1.0, 1.0), (1, 1), 1.0),
    ((5e-16, 0.5), (1, 3), 2.0000000000000005e-15),
]


def _r_acat_unweighted():
    from tests.unit.test_meta_multimask import R_ACAT  # same R script, unweighted section

    return R_ACAT


@pytest.mark.parametrize(("ps", "expected"), _r_acat_unweighted())
def test_cauchy_matches_r_acat_unweighted(ps, expected):
    assert cauchy_combination(list(ps)) == pytest.approx(expected, rel=1e-12, abs=1e-15)


@pytest.mark.parametrize(("ps", "weights", "expected"), R_ACAT_WEIGHTED)
def test_cauchy_matches_r_acat_weighted(ps, weights, expected):
    got = cauchy_combination(list(ps), weights=list(weights))
    assert got == pytest.approx(expected, rel=1e-12, abs=1e-15)


@pytest.mark.parametrize("bad", [1.0000001, -1e-12, 2.0])
def test_cauchy_rejects_p_outside_unit_interval(bad):
    """R ACAT() stops on p outside [0, 1]; so does cauchy_combination."""
    with pytest.raises(ValueError, match="between 0 and 1"):
        cauchy_combination([0.01, bad])


def test_cauchy_exact_zero_and_one_is_undefined():
    """R ACAT() stops on p = 0 together with p = 1: the combination is undefined (None)."""
    assert cauchy_combination([0.0, 1.0, 0.5]) is None
    assert cauchy_combination([None, 0.0, float("nan")]) == 0.0  # single p passes through
    assert cauchy_combination([0.0, 0.5], weights=[0.0, 1.0]) is None  # R: 0/0 -> NaN
