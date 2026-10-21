"""Calibrated variances in score artifacts and meta-analysis (#223, #198).

``v_robust`` = R ``SKATBinary_Robust``'s calibrated score variance (``VarS``), checked
against the R references of tests/fixtures/skat_robust_parity via the kernel.
"""

from __future__ import annotations

import json
import logging
import pathlib

import numpy as np
import pytest
import scipy.stats

from tests.unit.meta_helpers import make_stratum
from tests.unit.test_skat_robust_parity import MOMENTS, _load
from variantcentrifuge.association.meta.align import align_strata
from variantcentrifuge.association.meta.artifact import (
    GeneScore,
    compute_gene_score,
    fit_stratum_null,
    load_artifact,
    save_artifact,
)
from variantcentrifuge.association.meta.combine import meta_burden, meta_skat
from variantcentrifuge.association.meta.manifest import validate_comparability

V0230 = pathlib.Path(__file__).parent.parent / "fixtures" / "score_artifact_v0230" / "old"
KEYS = {"G1": [f"1:{100 + j}:A:G" for j in range(3)], "G2": [f"2:{200 + j}:C:T" for j in range(4)]}


def _with_v_robust(art, scale: float = 0.5):
    for gs in art.genes.values():
        gs.v_robust = scale * np.diag(gs.v)
    return art


def test_v0230_artifact_loads_without_v_robust():
    art = load_artifact(str(V0230))  # written by v0.23.0 (schema 2)
    assert art.schema_version == 2
    assert set(art.genes) == {"G1", "G2"}
    assert all(gs.v_robust is None for gs in art.genes.values())


def test_v_robust_round_trips_and_default_artifacts_have_no_key(tmp_path):
    art = _with_v_robust(make_stratum("a", seed=1, gene_keys=KEYS))
    save_artifact(art, str(tmp_path / "a"))
    back = load_artifact(str(tmp_path / "a"))
    assert back.schema_version == 3  # v_robust present
    for gene, gs in art.genes.items():
        np.testing.assert_array_equal(back.genes[gene].v_robust, gs.v_robust)
    plain = make_stratum("b", seed=2, gene_keys=KEYS)
    save_artifact(plain, str(tmp_path / "b"))
    with np.load(tmp_path / "b.npz") as data:
        assert not [k for k in data.files if k.endswith("::v_robust")]
    # Without v_robust the artifact is a v0.23.0 (schema 2) artifact, readable there
    assert json.loads((tmp_path / "b.meta.json").read_text())["schema_version"] == 2
    assert load_artifact(str(tmp_path / "b")).schema_version == 2


@pytest.mark.parametrize(("cohort", "gene"), [("r19", "mixed"), ("r81", "le10")])
def test_compute_gene_score_v_robust_matches_r(cohort, gene):
    # R's Phi diagonal is w^2 VarS / r; with r = 1 (no burden shrinkage) diag(Phi)/w^2 = VarS.
    geno, w, y, _mu, design = _load(cohort, gene)
    null = fit_stratum_null(y, design[:, 1:], "binary")
    keys = [f"1:{i}:A:G" for i in range(geno.shape[1])]
    gs = compute_gene_score(geno, geno, null, keys, robust=True)
    assert gs.v_robust is not None and compute_gene_score(geno, geno, null, keys).v_robust is None
    r_phi = np.array([row for _, row in MOMENTS[(cohort, gene)]])
    r_var = np.diag(r_phi) / w**2
    ratio = gs.v_robust / r_var  # = r (one number per gene) up to the GLM fit tolerance
    np.testing.assert_allclose(ratio, ratio[0], rtol=1e-6)
    assert 0 < ratio[0] <= 1.0 + 1e-6


def test_explained_columns_get_zero_v_robust():
    geno, _w, y, _mu, design = _load("r19", "mixed")
    null = fit_stratum_null(y, design[:, 1:], "binary")
    geno = np.column_stack([geno, np.ones(y.size)])  # in the intercept span
    keys = [f"1:{i}:A:G" for i in range(geno.shape[1])]
    gs = compute_gene_score(geno, geno, null, keys, robust=True)
    assert gs.v_robust is not None and gs.v_robust[-1] == 0.0 and (gs.v_robust[:-1] > 0).all()


def test_align_uses_calibrated_covariance_only_when_every_stratum_has_it():
    a = make_stratum("a", seed=1, gene_keys=KEYS)
    b = make_stratum("b", seed=2, gene_keys=KEYS)
    raw = align_strata([a.genes["G2"], b.genes["G2"]])[1]
    _with_v_robust(a, 0.5)
    mixed = align_strata([a.genes["G2"], b.genes["G2"]])[1]
    for got, want in zip(mixed, raw, strict=True):
        np.testing.assert_array_equal(got, want)
    _with_v_robust(b, 0.25)
    cal = align_strata([a.genes["G2"], b.genes["G2"]])[1]
    np.testing.assert_allclose(cal[0], 0.5 * raw[0], rtol=1e-12)  # D = sqrt(0.5) I
    np.testing.assert_allclose(cal[1], 0.25 * raw[1], rtol=1e-12)


def test_meta_with_calibrated_variance_is_less_significant():
    a = make_stratum("a", seed=1, gene_keys=KEYS, effect_gene="G2", effect_sign=1)
    b = make_stratum("b", seed=2, gene_keys=KEYS, effect_gene="G2", effect_sign=1)
    _with_v_robust(a, 2.0)
    _with_v_robust(b, 2.0)
    w = np.ones(4)
    u, v_cal, *_ = align_strata([a.genes["G2"], b.genes["G2"]])
    raw = [gs.v for gs in (a.genes["G2"], b.genes["G2"])]
    # Inflated variance (V~ = 2 V): larger burden and SKAT p-values
    assert meta_burden(u, raw, w)[0] < meta_burden(u, v_cal, w)[0]  # type: ignore[operator]
    assert meta_skat(u, raw, w)[0] < meta_skat(u, v_cal, w)[0]  # type: ignore[operator]


def test_mixed_strata_warn_once_and_fall_back(caplog):
    old = load_artifact(str(V0230))
    new = _with_v_robust(make_stratum("new", seed=3, gene_keys=KEYS))
    with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
        validate_comparability([old, new])
    msgs = [r.getMessage() for r in caplog.records if "--skat-robust" in r.getMessage()]
    assert len(msgs) == 1 and "old=default" in msgs[0] and "new=robust" in msgs[0]
    assert all(gs.v_robust is None for gs in new.genes.values())


def test_all_robust_strata_keep_v_robust(caplog):
    a = _with_v_robust(make_stratum("a", seed=1, gene_keys=KEYS))
    b = _with_v_robust(make_stratum("b", seed=2, gene_keys=KEYS))
    with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
        validate_comparability([a, b])
    assert not [r for r in caplog.records if "Strata mix" in r.getMessage()]
    info = [r for r in caplog.records if "calibrated variances" in r.getMessage()]
    assert len(info) == 1 and info[0].levelno == logging.INFO
    assert all(isinstance(gs, GeneScore) and gs.v_robust is not None for gs in b.genes.values())


def test_stage_hook_and_emit_store_v_robust_only_when_robust(tmp_path):
    from tests.unit.test_missing_policy import _builder, _sparse_gene
    from variantcentrifuge.association.meta.artifact import emit_stratum_artifact
    from variantcentrifuge.stages.association_stage_helpers import _score_artifact_hook

    df, samples = _sparse_gene(carriers=3)  # unbalanced carriers: nonzero scores
    builder = _builder(df, samples, missing_as_ref=True)
    kw = {"stratum_id": "s", "genome_build": "b", "annotation_version": "a"}
    kw |= {"mask_id": "m", "mask_hash": "h"}
    plain = emit_stratum_artifact(str(tmp_path / "p"), {"G1": builder}, **kw)
    robust = emit_stratum_artifact(str(tmp_path / "r"), {"G1": builder}, robust=True, **kw)
    assert plain.genes["G1"].v_robust is None
    v_rob = robust.genes["G1"].v_robust
    assert v_rob is not None and v_rob.shape == (4,) and (v_rob > 0).all()
    np.testing.assert_array_equal(robust.genes["G1"].u, plain.genes["G1"].u)
    np.testing.assert_array_equal(robust.genes["G1"].v, plain.genes["G1"].v)
    hook = _score_artifact_hook({"G1": builder}, robust=True)
    from variantcentrifuge.association.meta.artifact import _build_for_artifact

    entry = hook("G1", builder, _build_for_artifact(builder))[2]
    np.testing.assert_array_equal(entry.v_robust, v_rob)


def test_meta_effects_and_heterogeneity_use_raw_covariance():
    """U = 20, V = 40 in both strata, V~ = 160 / 40 (Codex probe): equal effects stay equal.

    v0.24 dev used V~ in beta = U / V and in Cochran's Q: 0.125 / 0.5, het p 0.034.
    """
    from types import SimpleNamespace

    from variantcentrifuge.association.meta.runner import _gene_row

    def stratum(sid, v_robust):
        gs = GeneScore(
            variant_keys=["1:100:A:G"],
            u=np.array([20.0]),
            v=np.array([[40.0]]),
            ac=np.array([30.0]),
            an=np.array([2000.0]),
            null_converged=True,
            v_robust=np.array([v_robust]),
        )
        return SimpleNamespace(stratum_id=sid, n_cases=100, n_controls=900), gs

    row = _gene_row("G1", "m", [stratum("a", 160.0), stratum("b", 40.0)], "uniform")
    assert row["meta_burden_beta"] == pytest.approx(0.5)  # (20 + 20) / (40 + 40)
    assert row["het_p"] == pytest.approx(1.0) and row["het_q"] == pytest.approx(0.0)
    # the association p-value is still calibrated: z = 40 / sqrt(160 + 40)
    assert row["meta_burden_p"] == pytest.approx(2 * scipy.stats.norm.sf(40 / np.sqrt(200)))
