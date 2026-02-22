# File: tests/unit/test_davies_fallback.py
# Location: tests/unit/test_davies_fallback.py
"""
Unit tests for saddlepoint-before-Liu fallback in compute_pvalue (Phase 25).

Verifies:
- When Davies returns out-of-range p-value, saddlepoint is tried before Liu
- When saddlepoint is valid (0 < p < 1), it is used and method is "saddlepoint"
- When saddlepoint fails, Liu is used as final fallback
- When Davies is in-range, behavior is unchanged (no saddlepoint inserted)
- Saddlepoint produces valid p-values for large q (extreme tail)
"""

from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest

from variantcentrifuge.association.backends.davies import (
    _kuonen_pvalue,
    _liu_pvalue,
    compute_pvalue,
)


@pytest.mark.unit
class TestDaviesOutOfRangeFallback:
    """Tests for the Davies out-of-range saddlepoint-before-Liu fallback."""

    def test_davies_out_of_range_tries_saddlepoint(self):
        """When Davies returns p > 1, saddlepoint is tried before Liu."""
        lambdas = np.array([2.0, 1.5, 1.0, 0.5])
        q = 20.0  # large q where saddlepoint should be valid

        # Mock davies_pvalue to return out-of-range result
        with (
            patch(
                "variantcentrifuge.association.backends.davies.davies_pvalue",
                return_value=(2.0, 0),  # p=2.0 is out of range
            ),
            patch(
                "variantcentrifuge.association.backends.davies._try_load_davies",
                return_value=True,
            ),
        ):
            p, method, converged = compute_pvalue(q, lambdas)

        # Must use saddlepoint (not Liu), since q=20.0 > mean(lambdas)=5.0
        assert method == "saddlepoint"
        assert 0.0 < p < 1.0
        assert converged is False

    def test_davies_out_of_range_zero_uses_saddlepoint(self):
        """When Davies returns p <= 0, saddlepoint is tried before Liu."""
        lambdas = np.array([2.0, 1.5, 1.0, 0.5])
        q = 20.0

        with (
            patch(
                "variantcentrifuge.association.backends.davies.davies_pvalue",
                return_value=(0.0, 0),  # p=0.0 is out of range (<=0)
            ),
            patch(
                "variantcentrifuge.association.backends.davies._try_load_davies",
                return_value=True,
            ),
        ):
            p, method, converged = compute_pvalue(q, lambdas)

        assert method == "saddlepoint"
        assert 0.0 < p < 1.0
        assert converged is False

    def test_davies_out_of_range_saddlepoint_fails_falls_to_liu(self):
        """When Davies is out-of-range AND saddlepoint fails, Liu is the fallback."""
        # q <= mean(lambdas) causes saddlepoint to return None
        lambdas = np.array([1.0, 1.0, 1.0, 1.0])
        q = 2.0  # q <= mean(lambdas)=4.0? No, 2.0 < 4.0, saddlepoint returns None

        with (
            patch(
                "variantcentrifuge.association.backends.davies.davies_pvalue",
                return_value=(2.0, 0),  # p=2.0 is out of range
            ),
            patch(
                "variantcentrifuge.association.backends.davies._try_load_davies",
                return_value=True,
            ),
        ):
            p, method, converged = compute_pvalue(q, lambdas)

        # q=2.0 < mean(lambdas)=4.0, so saddlepoint returns None -> fall to Liu
        assert method == "liu"
        assert 0.0 <= p <= 1.0
        assert converged is False

    def test_davies_out_of_range_saddlepoint_also_out_of_range_falls_to_liu(self):
        """When both Davies and saddlepoint are out of range, Liu is used."""
        lambdas = np.array([2.0, 1.5, 1.0, 0.5])
        q = 20.0

        # Mock _kuonen_pvalue to return None (out of range / not applicable)
        with (
            patch(
                "variantcentrifuge.association.backends.davies.davies_pvalue",
                return_value=(2.0, 0),  # Davies out of range
            ),
            patch(
                "variantcentrifuge.association.backends.davies._try_load_davies",
                return_value=True,
            ),
            patch(
                "variantcentrifuge.association.backends.davies._kuonen_pvalue",
                return_value=None,  # saddlepoint fails
            ),
        ):
            p, method, converged = compute_pvalue(q, lambdas)

        assert method == "liu"
        assert 0.0 <= p <= 1.0
        assert converged is False


@pytest.mark.unit
class TestDaviesInRangeBehaviorUnchanged:
    """Tests that in-range Davies behavior is unchanged by the new fallback."""

    def test_davies_in_range_returns_davies_method(self):
        """When Davies returns a valid p-value (0 < p <= 1), use Davies directly."""
        lambdas = np.array([2.0, 1.5, 1.0, 0.5])
        q = 10.0

        with (
            patch(
                "variantcentrifuge.association.backends.davies.davies_pvalue",
                return_value=(0.05, 0),  # p=0.05 is in range
            ),
            patch(
                "variantcentrifuge.association.backends.davies._try_load_davies",
                return_value=True,
            ),
        ):
            p, method, converged = compute_pvalue(q, lambdas)

        assert method == "davies"
        assert p == pytest.approx(0.05)
        assert converged is True  # ifault=0

    def test_davies_in_range_non_converged_kept(self):
        """When Davies returns p in range with ifault!=0, p is kept (non-converged)."""
        lambdas = np.array([2.0, 1.5, 1.0, 0.5])
        q = 10.0

        with (
            patch(
                "variantcentrifuge.association.backends.davies.davies_pvalue",
                return_value=(0.05, 1),  # p=0.05 in range but ifault=1
            ),
            patch(
                "variantcentrifuge.association.backends.davies._try_load_davies",
                return_value=True,
            ),
        ):
            p, method, converged = compute_pvalue(q, lambdas)

        assert method == "davies"
        assert p == pytest.approx(0.05)
        assert converged is False  # ifault != 0

    def test_no_davies_uses_saddlepoint_then_liu(self):
        """Without Davies C extension, saddlepoint -> Liu chain unchanged."""
        lambdas = np.array([2.0, 1.5, 1.0, 0.5])
        q = 20.0  # q > mean(lambdas)=5.0, saddlepoint valid

        with patch(
            "variantcentrifuge.association.backends.davies._try_load_davies",
            return_value=False,  # Davies not available
        ):
            p, method, converged = compute_pvalue(q, lambdas)

        # Should use saddlepoint when Davies unavailable and q > mean
        assert method == "saddlepoint"
        assert 0.0 < p < 1.0
        assert converged is False


@pytest.mark.unit
class TestSaddlepointValidForLargeQ:
    """Tests that saddlepoint produces valid p-values in the extreme tail."""

    def test_kuonen_pvalue_valid_for_large_q(self):
        """_kuonen_pvalue returns valid p-value in (0, 1) for q >> mean(lambdas)."""
        lambdas = np.array([2.0, 1.5, 1.0, 0.5])  # mean = 5.0
        q = 30.0  # large q, should produce small p-value

        p = _kuonen_pvalue(q, lambdas)

        assert p is not None
        assert 0.0 < p < 1.0
        # For large q, p should be small
        assert p < 0.01

    def test_kuonen_pvalue_returns_none_below_mean(self):
        """_kuonen_pvalue returns None when q <= mean(lambdas) (not valid for upper tail)."""
        lambdas = np.array([2.0, 1.5, 1.0, 0.5])  # mean = 5.0
        q = 3.0  # q < mean(lambdas) = 5.0

        p = _kuonen_pvalue(q, lambdas)

        assert p is None

    def test_saddlepoint_agrees_with_liu_approximately(self):
        """Saddlepoint and Liu should agree to within an order of magnitude."""
        lambdas = np.array([2.0, 1.5, 1.0, 0.5])
        q = 15.0  # q > mean(lambdas)=5.0

        p_sp = _kuonen_pvalue(q, lambdas)
        p_liu = _liu_pvalue(q, lambdas)

        assert p_sp is not None
        assert 0.0 < p_sp < 1.0
        assert 0.0 < p_liu < 1.0
        # Both should agree within a factor of 10x for this case
        log_ratio = abs(np.log10(p_sp) - np.log10(p_liu))
        assert log_ratio < 1.0, f"p_sp={p_sp}, p_liu={p_liu}: too different"
