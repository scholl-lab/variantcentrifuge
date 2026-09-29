"""Variant alignment and pooled-frequency weight harmonisation."""

from __future__ import annotations

import numpy as np
import pytest

from variantcentrifuge.association.meta.align import (
    align_strata,
    canonical_key,
    pooled_frequency_weights,
)
from variantcentrifuge.association.meta.artifact import GeneScore


def _gs(keys, u, ac, an):
    p = len(keys)
    return GeneScore(
        keys, np.array(u, float), np.eye(p), np.array(ac, float), np.array(an, float), True
    )


def test_canonical_key_strips_chr_and_uppercases():
    assert canonical_key("chr1", 100, "a", "g") == "1:100:A:G"
    assert canonical_key("1", "100", "A", "G") == "1:100:A:G"


def test_align_union_zero_pads_absent_variants():
    a = _gs(["1:1:A:G", "1:2:A:G"], [3.0, 5.0], [3, 5], [200, 200])
    b = _gs(["1:2:A:G", "1:3:A:G"], [7.0, 11.0], [7, 11], [300, 300])
    u_list, v_list, ac, an, keys = align_strata([a, b])

    assert keys == ["1:1:A:G", "1:2:A:G", "1:3:A:G"]
    # Stratum A has no data for 1:3 -> zero score there; B none for 1:1.
    np.testing.assert_allclose(u_list[0], [3.0, 5.0, 0.0])
    np.testing.assert_allclose(u_list[1], [0.0, 7.0, 11.0])
    # Pooled AC/AN summed on the shared variant only.
    np.testing.assert_allclose(ac, [3.0, 12.0, 11.0])
    np.testing.assert_allclose(an, [200.0, 500.0, 300.0])
    # Covariance zero-block for absent variants.
    assert v_list[0][2, 2] == 0.0
    assert v_list[1][0, 0] == 0.0


def test_pooled_beta_weight_uses_minor_allele_frequency():
    # af = 0.9 -> maf = 0.1; weight equals Beta(0.1;1,25) pdf.
    w = pooled_frequency_weights(np.array([180.0]), np.array([200.0]), "beta:1,25")
    import scipy.stats

    assert w[0] == pytest.approx(scipy.stats.beta.pdf(0.1, 1, 25))


def test_uniform_weights_and_rejected_scheme():
    w = pooled_frequency_weights(np.array([1.0, 2.0]), np.array([100.0, 100.0]), "uniform")
    np.testing.assert_allclose(w, [1.0, 1.0])
    with pytest.raises(ValueError, match="Unsupported meta weight scheme"):
        pooled_frequency_weights(np.array([1.0]), np.array([100.0]), "cadd")
