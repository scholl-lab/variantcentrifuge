"""Heterogeneity + power acceptance criteria (#1 concordant, #2 opposite)."""

from __future__ import annotations

import numpy as np

from tests.unit.meta_helpers import make_stratum
from variantcentrifuge.association.meta.align import align_strata, pooled_frequency_weights
from variantcentrifuge.association.meta.combine import meta_burden
from variantcentrifuge.association.meta.heterogeneity import (
    cochran_q,
    leave_one_out,
    per_stratum_burden,
)

GENE = "GENEA"
KEYS = {GENE: [f"1:{100 + j}:A:G" for j in range(12)]}


def _combine_two(art_a, art_b, scheme="beta:1,25"):
    gs = [art_a.genes[GENE], art_b.genes[GENE]]
    u_list, v_list, ac, an, _keys = align_strata(gs)
    w = pooled_frequency_weights(ac, an, scheme)
    return u_list, v_list, w


def test_concordant_effects_gain_evidence():
    """Acceptance #1: same-direction strata → meta beats either stratum alone."""
    a = make_stratum("A", seed=101, gene_keys=KEYS, effect_gene=GENE, effect_sign=+1.0)
    b = make_stratum("B", seed=202, gene_keys=KEYS, effect_gene=GENE, effect_sign=+1.0)

    u_list, v_list, w = _combine_two(a, b)
    meta_p, _beta, direction = meta_burden(u_list, v_list, w)
    p_a, _, _ = meta_burden([u_list[0]], [v_list[0]], w)
    p_b, _, _ = meta_burden([u_list[1]], [v_list[1]], w)

    assert meta_p < min(p_a, p_b)
    assert direction == 1


def test_opposite_directions_are_distinguishable():
    """Acceptance #2: opposite strata → common-effect attenuated, heterogeneity flagged."""
    a = make_stratum(
        "A", seed=303, gene_keys=KEYS, effect_gene=GENE, effect_sign=+1.0, effect_size=1.6
    )
    b = make_stratum(
        "B", seed=404, gene_keys=KEYS, effect_gene=GENE, effect_sign=-1.0, effect_size=1.6
    )

    u_list, v_list, w = _combine_two(a, b)
    burden_p, _beta, _direction = meta_burden(u_list, v_list, w)
    betas, ses = per_stratum_burden(u_list, v_list, w)
    _q, het_p, i2 = cochran_q(betas, ses)

    # Opposite per-stratum directions and heterogeneity strictly more significant
    # than the (cancelled) common-effect burden.
    assert np.sign(betas[0]) != np.sign(betas[1])
    assert het_p < burden_p
    assert i2 > 0.5


def test_leave_one_out_returns_per_stratum_dropped_pvalues():
    a = make_stratum("A", seed=505, gene_keys=KEYS, effect_gene=GENE, effect_sign=+1.0)
    b = make_stratum("B", seed=606, gene_keys=KEYS, effect_gene=GENE, effect_sign=+1.0)
    u_list, v_list, w = _combine_two(a, b)
    loo = leave_one_out(u_list, v_list, w, ["A", "B"], lambda u, v, ww: meta_burden(u, v, ww)[0])
    assert [sid for sid, _p in loo] == ["A", "B"]
    assert all(p is not None for _sid, p in loo)


def test_cochran_q_df_uses_usable_strata_only():
    # One stratum has zero information (SE inf) -> ignored; usable df collapses.
    q, p, i2 = cochran_q(np.array([0.5, np.nan]), np.array([0.1, np.inf]))
    assert (q, p, i2) == (0.0, 1.0, 0.0)
