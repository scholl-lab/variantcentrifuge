"""
Unit tests for variantcentrifuge.association.weights.

Tests cover:
- beta_maf_weights: rare variants upweighted, monotonic decrease, custom parameters (WEIGHT-01)
- uniform_weights: all-ones, correct length (WEIGHT-02)
- get_weights: string spec dispatch, invalid spec raises ValueError
"""

from __future__ import annotations

import numpy as np
import pytest

from variantcentrifuge.association.weights import beta_maf_weights, get_weights, uniform_weights

# ---------------------------------------------------------------------------
# beta_maf_weights tests (WEIGHT-01)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestBetaMafWeights:
    """Beta(MAF; a, b) weight scheme — WEIGHT-01."""

    def test_beta_weights_rare_upweighted(self) -> None:
        """Rare variants (low MAF) receive higher weights than common variants."""
        mafs = np.array([0.001, 0.01, 0.1])
        weights = beta_maf_weights(mafs)

        # Rare variant (MAF=0.001) should have higher weight than common (MAF=0.1)
        assert weights[0] > weights[1], "MAF=0.001 should have higher weight than MAF=0.01"
        assert weights[1] > weights[2], "MAF=0.01 should have higher weight than MAF=0.1"

    def test_beta_weights_monotonic_decreasing(self) -> None:
        """Beta weights are monotonically non-increasing with increasing MAF."""
        # Dense MAF grid from rare to common
        mafs = np.linspace(0.001, 0.5, 50)
        weights = beta_maf_weights(mafs)

        # Each weight should be >= the next (monotonically non-increasing)
        for i in range(len(weights) - 1):
            assert weights[i] >= weights[i + 1], (
                f"Weights not monotonically decreasing at index {i}: "
                f"w[{i}]={weights[i]:.4f} < w[{i + 1}]={weights[i + 1]:.4f}"
            )

    def test_beta_weights_all_positive(self) -> None:
        """All weights are strictly positive for MAFs in (0, 1)."""
        mafs = np.array([0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5])
        weights = beta_maf_weights(mafs)

        assert (weights > 0).all(), f"All weights must be positive, got {weights}"

    def test_beta_weights_default_params(self) -> None:
        """Default parameters (a=1, b=25) match SKAT R convention."""
        mafs = np.array([0.01])
        w1 = beta_maf_weights(mafs)
        w2 = beta_maf_weights(mafs, a=1.0, b=25.0)

        assert w1[0] == pytest.approx(w2[0]), "Default params should match a=1, b=25"

    def test_beta_weights_custom_params(self) -> None:
        """Custom parameters (a=0.5, b=0.5) produce different weights."""
        mafs = np.array([0.01, 0.05, 0.1])
        w_default = beta_maf_weights(mafs)
        w_custom = beta_maf_weights(mafs, a=0.5, b=0.5)

        # With a=0.5, b=0.5 the Beta dist is U-shaped — different from default
        assert not np.allclose(w_default, w_custom), (
            "Custom a=0.5, b=0.5 should produce different weights than default a=1, b=25"
        )

    def test_beta_weights_output_shape(self) -> None:
        """Output shape matches input shape."""
        mafs = np.array([0.01, 0.05, 0.1, 0.2])
        weights = beta_maf_weights(mafs)

        assert weights.shape == mafs.shape

    def test_beta_weights_clipping_at_boundaries(self) -> None:
        """MAF=0 and MAF=1 are clipped; no inf or nan returned."""
        mafs = np.array([0.0, 1.0])
        weights = beta_maf_weights(mafs)

        assert not np.isnan(weights).any(), "No NaN expected for boundary MAFs"
        assert not np.isinf(weights).any(), "No inf expected for boundary MAFs"

    def test_beta_weights_output_dtype_float64(self) -> None:
        """Output is float64."""
        mafs = np.array([0.01, 0.1])
        weights = beta_maf_weights(mafs)

        assert weights.dtype == np.float64


# ---------------------------------------------------------------------------
# uniform_weights tests (WEIGHT-02)
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestUniformWeights:
    """Uniform weight scheme — WEIGHT-02."""

    def test_uniform_weights_all_ones(self) -> None:
        """All weights are 1.0."""
        weights = uniform_weights(5)

        np.testing.assert_array_equal(weights, np.ones(5))

    def test_uniform_weights_correct_length(self) -> None:
        """Length matches n_variants."""
        for n in [1, 10, 100]:
            weights = uniform_weights(n)
            assert len(weights) == n, f"Expected length {n}, got {len(weights)}"

    def test_uniform_weights_dtype(self) -> None:
        """Output is float64."""
        weights = uniform_weights(3)
        assert weights.dtype == np.float64

    def test_uniform_weights_equivalent_to_equal_weights(self) -> None:
        """Uniform weights produce same burden score as manually summing genotype rows."""
        np.random.seed(99)
        n_samples, n_variants = 10, 5
        gmat = np.random.randint(0, 3, size=(n_samples, n_variants)).astype(float)

        weights = uniform_weights(n_variants)
        burden_weighted = gmat @ weights  # each variant gets weight 1.0

        # Manual: row sum (sum each variant dosage per sample)
        burden_manual = gmat.sum(axis=1)

        (
            np.testing.assert_allclose(burden_weighted, burden_manual, rtol=1e-10),
            ("Uniform weights burden should equal manual row sum"),
        )


# ---------------------------------------------------------------------------
# get_weights dispatch tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestGetWeights:
    """String spec parsing and dispatch in get_weights."""

    def test_get_weights_beta_spec_default(self) -> None:
        """'beta:1,25' dispatches to beta_maf_weights with a=1, b=25."""
        mafs = np.array([0.01, 0.05, 0.1])
        w_get = get_weights(mafs, "beta:1,25")
        w_direct = beta_maf_weights(mafs, a=1.0, b=25.0)

        np.testing.assert_allclose(w_get, w_direct, rtol=1e-10)

    def test_get_weights_beta_spec_custom(self) -> None:
        """'beta:2,10' dispatches to beta_maf_weights with a=2, b=10."""
        mafs = np.array([0.01, 0.05, 0.1])
        w_get = get_weights(mafs, "beta:2,10")
        w_direct = beta_maf_weights(mafs, a=2.0, b=10.0)

        np.testing.assert_allclose(w_get, w_direct, rtol=1e-10)

    def test_get_weights_uniform_spec(self) -> None:
        """'uniform' returns all-ones vector."""
        mafs = np.array([0.01, 0.05, 0.1])
        weights = get_weights(mafs, "uniform")

        np.testing.assert_array_equal(weights, np.ones(3))

    def test_get_weights_invalid_spec_raises_value_error(self) -> None:
        """Unknown spec raises ValueError."""
        mafs = np.array([0.01, 0.05])

        with pytest.raises(ValueError, match=r"[Uu]nknown"):
            get_weights(mafs, "unknown_scheme")

    def test_get_weights_invalid_beta_format_raises(self) -> None:
        """Malformed beta spec raises ValueError."""
        mafs = np.array([0.01, 0.05])

        with pytest.raises(ValueError):
            get_weights(mafs, "beta:1")  # missing second parameter

    def test_get_weights_uniform_ignores_maf_values(self) -> None:
        """Uniform weights don't depend on MAF values."""
        mafs_a = np.array([0.001, 0.01, 0.1])
        mafs_b = np.array([0.5, 0.3, 0.2])

        w_a = get_weights(mafs_a, "uniform")
        w_b = get_weights(mafs_b, "uniform")

        # Both should be all-ones of the same length
        np.testing.assert_array_equal(w_a, w_b)

    def test_get_weights_output_shape_matches_input(self) -> None:
        """Output length matches number of variants for both beta and uniform."""
        mafs = np.array([0.01, 0.05, 0.1, 0.2])

        for spec in ["beta:1,25", "uniform"]:
            weights = get_weights(mafs, spec)
            assert len(weights) == len(mafs), f"Length mismatch for spec '{spec}'"
