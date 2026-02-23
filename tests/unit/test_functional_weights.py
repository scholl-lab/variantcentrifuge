# File: tests/unit/test_functional_weights.py
"""
Unit tests for functional variant weight schemes: CADD, REVEL, and combined.

Covers:
- cadd_weights(): Beta(MAF) x min(CADD_phred / cap, 1.0), missing fallback
- revel_weights(): Beta(MAF) x REVEL_score, missing fallback
- combined_weights(): CADD preferred, falls back to REVEL or Beta-only
- get_weights() dispatch for "cadd", "revel", "combined", "beta:*", "uniform"
- Type-aware fallback logging (LoF / missense / other counts in warning)
- Functional weights produce DIFFERENT numerical results from uniform weights
- Existing "beta:1,25" and "uniform" specs unchanged (regression check)
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from variantcentrifuge.association.weights import (
    LOF_EFFECTS,
    MISSENSE_EFFECTS,
    beta_maf_weights,
    cadd_weights,
    combined_weights,
    get_weights,
    revel_weights,
    uniform_weights,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def simple_mafs() -> np.ndarray:
    """Three variants with MAFs: very rare, rare, less rare."""
    return np.array([0.001, 0.01, 0.05], dtype=np.float64)


@pytest.fixture
def simple_cadd() -> np.ndarray:
    """CADD Phred scores for three variants: high, moderate, low."""
    return np.array([35.0, 20.0, 5.0], dtype=np.float64)


@pytest.fixture
def simple_revel() -> np.ndarray:
    """REVEL scores for three variants (already in [0, 1])."""
    return np.array([0.9, 0.5, 0.1], dtype=np.float64)


# ---------------------------------------------------------------------------
# cadd_weights() tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestCaddWeights:
    """Tests for cadd_weights() functional weight scheme."""

    def test_normal_case_formula(self, simple_mafs: np.ndarray, simple_cadd: np.ndarray) -> None:
        """cadd_weights = Beta(MAF) * min(CADD/40, 1.0)."""
        result = cadd_weights(simple_mafs, simple_cadd)
        maf_w = beta_maf_weights(simple_mafs)
        expected = maf_w * np.array([35.0 / 40.0, 20.0 / 40.0, 5.0 / 40.0])
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_cadd_above_cap_capped_at_1(self, simple_mafs: np.ndarray) -> None:
        """CADD > 40 is capped at 1.0 (functional weight = 1.0)."""
        cadd = np.array([50.0, 40.0, 45.0])
        result = cadd_weights(simple_mafs, cadd)
        maf_w = beta_maf_weights(simple_mafs)
        expected = maf_w * np.array([1.0, 1.0, 1.0])
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_missing_cadd_lof_gets_max_functional_weight(self, simple_mafs: np.ndarray) -> None:
        """LoF variant missing CADD -> functional = 1.0 (max weight, no penalty)."""
        cadd = np.array([np.nan, 20.0, 5.0])
        effects = np.array(["stop_gained", "missense_variant", "synonymous_variant"])
        result = cadd_weights(simple_mafs, cadd, variant_effects=effects)
        maf_w = beta_maf_weights(simple_mafs)
        # LoF missing -> functional=1.0
        assert result[0] == pytest.approx(maf_w[0] * 1.0)

    def test_missing_cadd_missense_fallback(self, simple_mafs: np.ndarray) -> None:
        """Missense variant missing CADD -> functional = 1.0 (Beta(MAF)-only)."""
        cadd = np.array([20.0, np.nan, 5.0])
        effects = np.array(["stop_gained", "missense_variant", "synonymous_variant"])
        result = cadd_weights(simple_mafs, cadd, variant_effects=effects)
        maf_w = beta_maf_weights(simple_mafs)
        # Missense missing -> functional=1.0
        assert result[1] == pytest.approx(maf_w[1] * 1.0)

    def test_missing_cadd_no_effects_uses_fallback(self, simple_mafs: np.ndarray) -> None:
        """When variant_effects not provided, missing CADD still gets functional=1.0."""
        cadd = np.array([np.nan, 20.0, np.nan])
        result = cadd_weights(simple_mafs, cadd)
        maf_w = beta_maf_weights(simple_mafs)
        assert result[0] == pytest.approx(maf_w[0] * 1.0)
        assert result[2] == pytest.approx(maf_w[2] * 1.0)

    def test_custom_cap_applied(self, simple_mafs: np.ndarray) -> None:
        """Custom cap of 30.0 normalizes CADD differently."""
        cadd = np.array([30.0, 15.0, 5.0])
        result = cadd_weights(simple_mafs, cadd, cap=30.0)
        maf_w = beta_maf_weights(simple_mafs)
        expected = maf_w * np.array([1.0, 0.5, 5.0 / 30.0])
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_dot_string_treated_as_missing(self, simple_mafs: np.ndarray) -> None:
        """'.' string in CADD scores is treated as missing (functional=1.0)."""
        cadd_with_dot = np.array([".", "20.0", "5.0"], dtype=object)
        result = cadd_weights(simple_mafs, cadd_with_dot)
        maf_w = beta_maf_weights(simple_mafs)
        assert result[0] == pytest.approx(maf_w[0] * 1.0)
        assert result[1] == pytest.approx(maf_w[1] * 20.0 / 40.0)

    def test_none_string_treated_as_missing(self, simple_mafs: np.ndarray) -> None:
        """None in CADD scores is treated as missing (functional=1.0)."""
        cadd_with_none = np.array([None, 20.0, 5.0], dtype=object)
        result = cadd_weights(simple_mafs, cadd_with_none)
        maf_w = beta_maf_weights(simple_mafs)
        assert result[0] == pytest.approx(maf_w[0] * 1.0)


# ---------------------------------------------------------------------------
# revel_weights() tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestRevelWeights:
    """Tests for revel_weights() functional weight scheme."""

    def test_normal_case_formula(self, simple_mafs: np.ndarray, simple_revel: np.ndarray) -> None:
        """revel_weights = Beta(MAF) * REVEL_score (no normalization)."""
        result = revel_weights(simple_mafs, simple_revel)
        maf_w = beta_maf_weights(simple_mafs)
        expected = maf_w * simple_revel
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_missing_revel_lof_gets_max_weight(self, simple_mafs: np.ndarray) -> None:
        """LoF variant missing REVEL -> functional = 1.0."""
        revel = np.array([np.nan, 0.5, 0.1])
        effects = np.array(["frameshift_variant", "missense_variant", "synonymous_variant"])
        result = revel_weights(simple_mafs, revel, variant_effects=effects)
        maf_w = beta_maf_weights(simple_mafs)
        assert result[0] == pytest.approx(maf_w[0] * 1.0)

    def test_revel_is_already_zero_to_one(self, simple_mafs: np.ndarray) -> None:
        """REVEL scores in [0,1] are used directly without normalization."""
        revel = np.array([1.0, 0.5, 0.0])
        result = revel_weights(simple_mafs, revel)
        maf_w = beta_maf_weights(simple_mafs)
        # Score=1.0 -> full Beta(MAF) weight; score=0.0 -> zero weight
        assert result[0] == pytest.approx(maf_w[0] * 1.0)
        assert result[1] == pytest.approx(maf_w[1] * 0.5)
        assert result[2] == pytest.approx(maf_w[2] * 0.0)

    def test_missing_revel_no_effects_fallback(self, simple_mafs: np.ndarray) -> None:
        """Missing REVEL without variant_effects -> functional=1.0."""
        revel = np.array([np.nan, 0.7, np.nan])
        result = revel_weights(simple_mafs, revel)
        maf_w = beta_maf_weights(simple_mafs)
        assert result[0] == pytest.approx(maf_w[0] * 1.0)
        assert result[1] == pytest.approx(maf_w[1] * 0.7)
        assert result[2] == pytest.approx(maf_w[2] * 1.0)

    def test_dot_string_treated_as_missing(self, simple_mafs: np.ndarray) -> None:
        """'.' in REVEL treated as missing."""
        revel = np.array([".", "0.8", "0.2"], dtype=object)
        result = revel_weights(simple_mafs, revel)
        maf_w = beta_maf_weights(simple_mafs)
        assert result[0] == pytest.approx(maf_w[0] * 1.0)
        assert result[1] == pytest.approx(maf_w[1] * 0.8)


# ---------------------------------------------------------------------------
# combined_weights() tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestCombinedWeights:
    """Tests for combined_weights() functional weight scheme."""

    def test_uses_cadd_when_provided(
        self, simple_mafs: np.ndarray, simple_cadd: np.ndarray
    ) -> None:
        """combined_weights uses CADD when cadd_scores is provided."""
        result = combined_weights(simple_mafs, cadd_scores=simple_cadd)
        expected = cadd_weights(simple_mafs, simple_cadd)
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_uses_revel_when_no_cadd(
        self, simple_mafs: np.ndarray, simple_revel: np.ndarray
    ) -> None:
        """combined_weights falls back to REVEL when cadd_scores is None."""
        result = combined_weights(simple_mafs, revel_scores=simple_revel)
        expected = revel_weights(simple_mafs, simple_revel)
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_cadd_preferred_over_revel(
        self, simple_mafs: np.ndarray, simple_cadd: np.ndarray, simple_revel: np.ndarray
    ) -> None:
        """When both CADD and REVEL provided, CADD is preferred."""
        result = combined_weights(simple_mafs, cadd_scores=simple_cadd, revel_scores=simple_revel)
        expected = cadd_weights(simple_mafs, simple_cadd)
        np.testing.assert_allclose(result, expected, rtol=1e-10)
        # Verify it's different from REVEL-only
        revel_only = revel_weights(simple_mafs, simple_revel)
        assert not np.allclose(result, revel_only)

    def test_no_scores_falls_back_to_beta_maf(self, simple_mafs: np.ndarray) -> None:
        """Without any scores, combined_weights returns plain Beta(MAF)."""
        result = combined_weights(simple_mafs)
        expected = beta_maf_weights(simple_mafs)
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_custom_cadd_cap_via_param(self, simple_mafs: np.ndarray) -> None:
        """Custom cadd_cap is forwarded to cadd_weights."""
        cadd = np.array([30.0, 15.0, 5.0])
        result = combined_weights(simple_mafs, cadd_scores=cadd, cadd_cap=30.0)
        expected = cadd_weights(simple_mafs, cadd, cap=30.0)
        np.testing.assert_allclose(result, expected, rtol=1e-10)


# ---------------------------------------------------------------------------
# get_weights() dispatch tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestGetWeightsDispatch:
    """Tests for get_weights() string-spec dispatch including new functional specs."""

    def test_cadd_spec_dispatches_to_cadd_weights(
        self, simple_mafs: np.ndarray, simple_cadd: np.ndarray
    ) -> None:
        """get_weights('cadd') dispatches to cadd_weights."""
        result = get_weights(simple_mafs, "cadd", cadd_scores=simple_cadd)
        expected = cadd_weights(simple_mafs, simple_cadd)
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_revel_spec_dispatches_to_revel_weights(
        self, simple_mafs: np.ndarray, simple_revel: np.ndarray
    ) -> None:
        """get_weights('revel') dispatches to revel_weights."""
        result = get_weights(simple_mafs, "revel", revel_scores=simple_revel)
        expected = revel_weights(simple_mafs, simple_revel)
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_combined_spec_dispatches_to_combined_weights(
        self, simple_mafs: np.ndarray, simple_cadd: np.ndarray
    ) -> None:
        """get_weights('combined') dispatches to combined_weights."""
        result = get_weights(simple_mafs, "combined", cadd_scores=simple_cadd)
        expected = combined_weights(simple_mafs, cadd_scores=simple_cadd)
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_cadd_without_scores_raises_value_error(self, simple_mafs: np.ndarray) -> None:
        """get_weights('cadd') with no cadd_scores raises ValueError."""
        with pytest.raises(ValueError, match="cadd_scores"):
            get_weights(simple_mafs, "cadd")

    def test_revel_without_scores_raises_value_error(self, simple_mafs: np.ndarray) -> None:
        """get_weights('revel') with no revel_scores raises ValueError."""
        with pytest.raises(ValueError, match="revel_scores"):
            get_weights(simple_mafs, "revel")

    def test_cadd_with_weight_params_uses_custom_cap(self, simple_mafs: np.ndarray) -> None:
        """get_weights('cadd') with weight_params={'cadd_cap': 30} uses cap=30."""
        cadd = np.array([30.0, 15.0, 5.0])
        result = get_weights(
            simple_mafs, "cadd", cadd_scores=cadd, weight_params={"cadd_cap": 30.0}
        )
        expected = cadd_weights(simple_mafs, cadd, cap=30.0)
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_combined_with_weight_params_cadd_cap(self, simple_mafs: np.ndarray) -> None:
        """get_weights('combined') passes cadd_cap from weight_params."""
        cadd = np.array([30.0, 15.0, 5.0])
        result = get_weights(
            simple_mafs, "combined", cadd_scores=cadd, weight_params={"cadd_cap": 30.0}
        )
        expected = combined_weights(simple_mafs, cadd_scores=cadd, cadd_cap=30.0)
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    # Regression tests: existing specs must still work
    def test_beta_spec_regression(self, simple_mafs: np.ndarray) -> None:
        """get_weights('beta:1,25') still works (regression)."""
        result = get_weights(simple_mafs, "beta:1,25")
        expected = beta_maf_weights(simple_mafs, a=1.0, b=25.0)
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_uniform_spec_regression(self, simple_mafs: np.ndarray) -> None:
        """get_weights('uniform') still works (regression)."""
        result = get_weights(simple_mafs, "uniform")
        expected = uniform_weights(len(simple_mafs))
        np.testing.assert_allclose(result, expected, rtol=1e-10)

    def test_unknown_spec_raises_value_error(self, simple_mafs: np.ndarray) -> None:
        """Unknown weight spec raises ValueError."""
        with pytest.raises(ValueError, match="Unknown weight spec"):
            get_weights(simple_mafs, "invalid_spec")

    def test_existing_callers_unchanged_signature(self, simple_mafs: np.ndarray) -> None:
        """get_weights(mafs, spec) with no kwargs still works (backward compat)."""
        # This is the existing call signature — must not break
        result_beta = get_weights(simple_mafs, "beta:1,25")
        result_uniform = get_weights(simple_mafs, "uniform")
        assert result_beta.shape == (3,)
        assert result_uniform.shape == (3,)
        np.testing.assert_array_equal(result_uniform, [1.0, 1.0, 1.0])


# ---------------------------------------------------------------------------
# Type-aware fallback logging tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestFallbackLogging:
    """Tests that missing score counts are logged per variant category."""

    def test_cadd_missing_lof_logged(
        self, simple_mafs: np.ndarray, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Missing CADD on LoF variant -> log includes LoF count."""
        cadd = np.array([np.nan, 20.0, 5.0])
        effects = np.array(["stop_gained", "missense_variant", "synonymous_variant"])
        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            cadd_weights(simple_mafs, cadd, variant_effects=effects)
        # Should mention LoF count > 0
        assert any("LoF" in record.message for record in caplog.records)
        assert any("1 LoF" in record.message for record in caplog.records)

    def test_cadd_missing_missense_logged(
        self, simple_mafs: np.ndarray, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Missing CADD on missense variant -> log includes missense count."""
        cadd = np.array([20.0, np.nan, 5.0])
        effects = np.array(["stop_gained", "missense_variant", "synonymous_variant"])
        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            cadd_weights(simple_mafs, cadd, variant_effects=effects)
        assert any("missense" in record.message for record in caplog.records)
        assert any("1 missense" in record.message for record in caplog.records)

    def test_cadd_no_missing_no_warning(
        self, simple_mafs: np.ndarray, caplog: pytest.LogCaptureFixture
    ) -> None:
        """No missing CADD scores -> no warning logged."""
        cadd = np.array([35.0, 20.0, 5.0])
        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            cadd_weights(simple_mafs, cadd)
        assert not any("CADD" in record.message for record in caplog.records)

    def test_revel_missing_lof_logged(
        self, simple_mafs: np.ndarray, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Missing REVEL on LoF variant -> log includes LoF count."""
        revel = np.array([np.nan, 0.5, 0.1])
        effects = np.array(["splice_donor_variant", "missense_variant", "synonymous_variant"])
        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            revel_weights(simple_mafs, revel, variant_effects=effects)
        assert any("LoF" in record.message for record in caplog.records)

    def test_cadd_missing_without_effects_logged(
        self, simple_mafs: np.ndarray, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Missing CADD without variant_effects -> generic warning logged."""
        cadd = np.array([np.nan, 20.0, np.nan])
        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            cadd_weights(simple_mafs, cadd)
        assert any("missing" in record.message.lower() for record in caplog.records)

    def test_cadd_all_three_categories_in_one_warning(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """One warning message covers LoF, missense, and other counts."""
        mafs = np.array([0.001, 0.01, 0.05, 0.02])
        cadd = np.array([np.nan, np.nan, np.nan, np.nan])
        effects = np.array(
            ["stop_gained", "missense_variant", "synonymous_variant", "intron_variant"]
        )
        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            cadd_weights(mafs, cadd, variant_effects=effects)
        warning_msgs = [r.message for r in caplog.records if "CADD" in r.message]
        assert len(warning_msgs) == 1
        msg = warning_msgs[0]
        assert "1 LoF" in msg
        assert "1 missense" in msg
        assert "2 other" in msg


# ---------------------------------------------------------------------------
# Numerical difference from uniform (confirms weights are actually applied)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestWeightsDifferFromUniform:
    """Confirm functional weights produce different results from uniform."""

    def test_cadd_weights_differ_from_uniform(
        self, simple_mafs: np.ndarray, simple_cadd: np.ndarray
    ) -> None:
        """cadd_weights must produce different values than uniform weights."""
        cadd_w = cadd_weights(simple_mafs, simple_cadd)
        unif_w = uniform_weights(len(simple_mafs))
        assert not np.allclose(cadd_w, unif_w), (
            "CADD weights should differ from uniform for non-trivial scores"
        )

    def test_revel_weights_differ_from_uniform(
        self, simple_mafs: np.ndarray, simple_revel: np.ndarray
    ) -> None:
        """revel_weights must produce different values than uniform weights."""
        revel_w = revel_weights(simple_mafs, simple_revel)
        unif_w = uniform_weights(len(simple_mafs))
        assert not np.allclose(revel_w, unif_w), (
            "REVEL weights should differ from uniform for non-trivial scores"
        )

    def test_cadd_weights_differ_from_beta_maf_only(
        self, simple_mafs: np.ndarray, simple_cadd: np.ndarray
    ) -> None:
        """cadd_weights (Beta x CADD) differs from Beta(MAF)-only weights."""
        cadd_w = cadd_weights(simple_mafs, simple_cadd)
        beta_w = beta_maf_weights(simple_mafs)
        assert not np.allclose(cadd_w, beta_w), (
            "CADD weights should differ from Beta(MAF)-only when CADD != 1.0"
        )

    def test_revel_weights_differ_from_beta_maf_only(
        self, simple_mafs: np.ndarray, simple_revel: np.ndarray
    ) -> None:
        """revel_weights (Beta x REVEL) differs from Beta(MAF)-only weights."""
        revel_w = revel_weights(simple_mafs, simple_revel)
        beta_w = beta_maf_weights(simple_mafs)
        assert not np.allclose(revel_w, beta_w), (
            "REVEL weights should differ from Beta(MAF)-only when REVEL != 1.0"
        )

    def test_get_weights_cadd_differs_from_get_weights_uniform(
        self, simple_mafs: np.ndarray, simple_cadd: np.ndarray
    ) -> None:
        """get_weights('cadd') produces different results from get_weights('uniform')."""
        cadd_w = get_weights(simple_mafs, "cadd", cadd_scores=simple_cadd)
        unif_w = get_weights(simple_mafs, "uniform")
        assert not np.allclose(cadd_w, unif_w)

    def test_get_weights_revel_differs_from_get_weights_uniform(
        self, simple_mafs: np.ndarray, simple_revel: np.ndarray
    ) -> None:
        """get_weights('revel') produces different results from get_weights('uniform')."""
        revel_w = get_weights(simple_mafs, "revel", revel_scores=simple_revel)
        unif_w = get_weights(simple_mafs, "uniform")
        assert not np.allclose(revel_w, unif_w)


# ---------------------------------------------------------------------------
# Constants tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestVariantEffectConstants:
    """Tests for LOF_EFFECTS and MISSENSE_EFFECTS constants."""

    def test_lof_effects_are_frozenset(self) -> None:
        """LOF_EFFECTS is a frozenset."""
        assert isinstance(LOF_EFFECTS, frozenset)

    def test_lof_effects_contains_canonical_effects(self) -> None:
        """LOF_EFFECTS includes standard LoF variant types."""
        assert "stop_gained" in LOF_EFFECTS
        assert "frameshift_variant" in LOF_EFFECTS
        assert "splice_acceptor_variant" in LOF_EFFECTS
        assert "splice_donor_variant" in LOF_EFFECTS

    def test_missense_effects_contains_missense(self) -> None:
        """MISSENSE_EFFECTS includes missense_variant."""
        assert "missense_variant" in MISSENSE_EFFECTS

    def test_lof_and_missense_are_disjoint(self) -> None:
        """LOF_EFFECTS and MISSENSE_EFFECTS must not overlap."""
        assert LOF_EFFECTS.isdisjoint(MISSENSE_EFFECTS)
