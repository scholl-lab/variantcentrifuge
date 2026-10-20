"""Acceptance #3: null calibration under imbalance + stable degenerate-input handling.

The common-effect score meta uses a normal/Davies approximation. It is well
calibrated under MODERATE imbalance; SEVERE imbalance requires the deferred
genotype-level saddlepoint (SPA) backend. These tests assert the honest,
documented behaviour — gross calibration at moderate imbalance, and graceful
failure (never a crash) for degenerate genes — not exact severe-imbalance tails.
"""

from __future__ import annotations

import numpy as np

from variantcentrifuge.association.meta.align import align_strata, pooled_frequency_weights
from variantcentrifuge.association.meta.artifact import build_stratum_artifact
from variantcentrifuge.association.meta.combine import meta_burden, meta_skat, meta_skato


def _null_stratum(stratum_id, seed, n, prevalence, n_genes, p):
    rng = np.random.default_rng(seed)
    genes = {}
    for gi in range(n_genes):
        mafs = rng.uniform(0.02, 0.06, size=p)
        g = rng.binomial(2, mafs[np.newaxis, :], size=(n, p)).astype(float)
        genes[f"G{gi}"] = (g, g.copy(), [f"1:{gi * 100 + j}:A:G" for j in range(p)])
    intercept = np.log(prevalence / (1 - prevalence))
    y = (rng.uniform(size=n) < prevalence).astype(float)
    _ = intercept  # prevalence-driven; no genotype effect (null)
    return build_stratum_artifact(
        stratum_id=stratum_id,
        genes=genes,
        phenotype=y,
        covariates=None,
        trait_type="binary",
        genome_build="GRCh38",
        annotation_version="v1",
        mask_id="pLoF",
        mask_hash="h1",
        n_cases=int(y.sum()),
        n_controls=int((1 - y).sum()),
        sample_id_hashes=[f"{stratum_id}_{i}" for i in range(n)],
    )


def test_null_calibration_moderate_imbalance():
    n_genes = 150
    a = _null_stratum("A", seed=1, n=400, prevalence=0.12, n_genes=n_genes, p=8)
    b = _null_stratum("B", seed=2, n=400, prevalence=0.12, n_genes=n_genes, p=8)

    burden_ps, skat_ps = [], []
    for gi in range(n_genes):
        gene = f"G{gi}"
        u_list, v_list, ac, an, _keys = align_strata([a.genes[gene], b.genes[gene]])
        w = pooled_frequency_weights(ac, an, "beta:1,25")
        bp, _beta, _d = meta_burden(u_list, v_list, w)
        sp, _m = meta_skat(u_list, v_list, w)
        if bp is not None:
            burden_ps.append(bp)
        if sp is not None:
            skat_ps.append(sp)

    burden_ps = np.array(burden_ps)
    skat_ps = np.array(skat_ps)
    # Gross calibration: mean near 0.5 and no severe type-I inflation at alpha=0.05.
    assert 0.35 < burden_ps.mean() < 0.65
    assert 0.35 < skat_ps.mean() < 0.65
    assert (burden_ps < 0.05).mean() < 0.15
    assert (skat_ps < 0.05).mean() < 0.15


def test_monomorphic_gene_is_stable():
    """A gene monomorphic in every stratum must not crash — it yields no signal."""
    rng = np.random.default_rng(0)
    n, p = 200, 5
    zero = np.zeros((n, p))
    y = (rng.uniform(size=n) < 0.2).astype(float)
    keys = [f"1:{j}:A:G" for j in range(p)]
    art = build_stratum_artifact(
        stratum_id="A",
        genes={"MONO": (zero, zero.copy(), keys)},
        phenotype=y,
        covariates=None,
        trait_type="binary",
        genome_build="GRCh38",
        annotation_version="v1",
        mask_id="pLoF",
        mask_hash="h1",
        n_cases=int(y.sum()),
        n_controls=int((1 - y).sum()),
        sample_id_hashes=[f"A_{i}" for i in range(n)],
    )
    u_list, v_list, ac, an, _keys = align_strata([art.genes["MONO"], art.genes["MONO"]])
    w = pooled_frequency_weights(ac, an, "beta:1,25")

    burden_p, _beta, direction = meta_burden(u_list, v_list, w)
    skat_p, _method = meta_skat(u_list, v_list, w)
    skato_p, _rho = meta_skato(u_list, v_list, w)

    assert burden_p is None and direction == 0  # zero variance → undefined, not a crash
    assert skat_p in (None, 1.0)
    assert skato_p is None
