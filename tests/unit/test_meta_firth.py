"""Inverse-variance meta of per-stratum Firth burden log-ORs (#198, step 1).

Parity fixtures and references come from ``scripts/generate_firth_meta_parity.R``
(R 4.5.3, logistf 1.26.1, metafor 5.0.1; provenance line at the top of
``r_reference.csv``). R read the same CSVs, weighted each stratum's burden with
``dbeta(stratum MAF, 1, 25)``, fitted ``logistf(y ~ burden + c1 + c2)`` and took the
classic Firth Wald SE, the inverse Fisher information at logistf's estimate (logistf
1.26's ``fit$var`` uses the augmented-data information instead and its default CI is
profile likelihood; see the R script), then ``metafor::rma(yi, sei, method="FE")``. Stratum s2 has
788 cases : 101 controls and every G2 carrier there is a case (complete separation).

Tolerance ``rel=1e-6`` as in ``test_skat_null_parity.py``.
"""

from __future__ import annotations

import csv
import json
import pathlib

import numpy as np
import pandas as pd
import pytest

from tests.unit.meta_helpers import make_stratum
from variantcentrifuge.association.meta.artifact import load_artifact, save_artifact
from variantcentrifuge.association.meta.firth_meta import (
    FIRTH_COLUMNS,
    firth_ivw_columns,
    ivw_meta,
    stratum_firth_burden,
)
from variantcentrifuge.association.meta.runner import OUTPUT_COLUMNS, meta_analyze_artifacts

FIXTURES = pathlib.Path(__file__).parent.parent / "fixtures" / "firth_meta_parity"
STRATA = ("s1", "s2", "s3")


def _r_reference() -> dict[tuple[str, str, str], float]:
    with open(FIXTURES / "r_reference.csv") as fh:
        rows = csv.DictReader(line for line in fh if not line.startswith("#"))
        return {(r["gene"], r["stratum"], r["quantity"]): float(r["value"]) for r in rows}


R_REF = _r_reference()


def _stratum_fit(stratum: str, gene: str) -> dict | None:
    d = pd.read_csv(FIXTURES / f"{stratum}.csv")
    geno = d.filter(regex=rf"^{gene}_\d+$").to_numpy(dtype=float)
    ac, an = geno.sum(axis=0), np.full(geno.shape[1], 2.0 * len(d))
    cov = d[["c1", "c2"]].to_numpy(dtype=float)
    return stratum_firth_burden(geno, ac, an, d["y"].to_numpy(dtype=float), cov, "beta:1,25")


@pytest.mark.parametrize("gene", ["G1", "G2"])
@pytest.mark.parametrize("stratum", STRATA)
def test_stratum_firth_matches_logistf(stratum, gene):
    fit = _stratum_fit(stratum, gene)
    assert fit is not None
    assert fit["beta"] == pytest.approx(R_REF[(gene, stratum, "beta")], rel=1e-6)
    assert fit["se"] == pytest.approx(R_REF[(gene, stratum, "se")], rel=1e-6)
    assert fit["weights"] == "beta:1,25"


@pytest.mark.parametrize("gene", ["G1", "G2"])
def test_ivw_meta_matches_metafor(gene):
    fits = [_stratum_fit(s, gene) for s in STRATA]
    beta, se, p = ivw_meta([f["beta"] for f in fits], [f["se"] for f in fits])
    assert beta == pytest.approx(R_REF[(gene, "meta", "beta")], rel=1e-6)
    assert se == pytest.approx(R_REF[(gene, "meta", "se")], rel=1e-6)
    assert p == pytest.approx(R_REF[(gene, "meta", "p")], rel=1e-6)


def test_ivw_meta_skips_non_finite_and_handles_empty():
    assert ivw_meta([1.0, np.nan], [0.5, 0.2]) == pytest.approx((1.0, 0.5, 0.04550026389635842))
    assert ivw_meta([np.nan], [np.nan]) == (None, None, None)
    assert ivw_meta([], []) == (None, None, None)


def test_stratum_firth_none_without_burden_variance_or_supported_weights():
    y = np.array([0.0, 1.0] * 10)
    geno = np.zeros((20, 2))
    assert stratum_firth_burden(geno, np.zeros(2), np.full(2, 40.0), y, None, "beta:1,25") is None
    geno[:3, 0] = 1
    assert stratum_firth_burden(geno, np.array([3.0, 0]), np.full(2, 40.0), y, None, "cadd") is None
    fit = stratum_firth_burden(geno, np.array([3.0, 0]), np.full(2, 40.0), y, None, "uniform")
    assert fit is not None and np.isfinite(fit["beta"]) and fit["weights"] == "uniform"


def _with_firth(art, beta, se, weights="beta:1,25"):
    for gs in art.genes.values():
        gs.firth = {"beta": beta, "se": se, "weights": weights}
    return art


KEYS = {"GENE_A": ["1:100:A:G", "1:200:C:T", "1:300:G:A"]}


def test_runner_columns_from_artifact_firth_fields():
    arts = [
        _with_firth(make_stratum("s1", 1, gene_keys=KEYS), 0.5, 0.25),
        _with_firth(make_stratum("s2", 2, gene_keys=KEYS), 0.2, 0.5),
        _with_firth(make_stratum("s3", 3, gene_keys=KEYS), 0.9, 0.3, weights="uniform"),
    ]
    row = meta_analyze_artifacts(arts).iloc[0]
    beta, se, p = ivw_meta([0.5, 0.2], [0.25, 0.5])  # s3 used other weights: excluded
    assert row["burden_firth_ivw_beta"] == pytest.approx(beta)
    assert row["burden_firth_ivw_se"] == pytest.approx(se)
    assert row["burden_firth_ivw_p"] == pytest.approx(p)
    assert row["burden_firth_ivw_n_strata"] == 2
    assert set(FIRTH_COLUMNS) <= set(OUTPUT_COLUMNS)
    # alternative p only: primary unchanged by the Firth fields
    plain = meta_analyze_artifacts(
        [make_stratum(s, i, gene_keys=KEYS) for i, s in enumerate(STRATA, 1)]
    )
    assert row["primary_p"] == plain.iloc[0]["primary_p"]


def test_runner_columns_empty_without_firth_fields():
    row = meta_analyze_artifacts(
        [make_stratum(s, i, gene_keys=KEYS) for i, s in enumerate(STRATA, 1)]
    ).iloc[0]
    assert row[["burden_firth_ivw_beta", "burden_firth_ivw_se", "burden_firth_ivw_p"]].isna().all()
    assert row["burden_firth_ivw_n_strata"] == 0


def test_firth_ivw_columns_scheme_match_is_normalised():
    from types import SimpleNamespace

    gs = [SimpleNamespace(firth={"beta": 0.1, "se": 0.1, "weights": "Beta:1,25"})]
    assert firth_ivw_columns(gs, " beta:1,25 ")["burden_firth_ivw_n_strata"] == 1


def test_firth_field_roundtrip_and_v1_compat(tmp_path):
    art = _with_firth(make_stratum("s1", 1, gene_keys=KEYS), 0.5, 0.25)
    path = str(tmp_path / "s1")
    save_artifact(art, path)
    assert load_artifact(path).genes["GENE_A"].firth == {
        "beta": 0.5,
        "se": 0.25,
        "weights": "beta:1,25",
    }
    # artifacts without the field (v1, or v2 written before #198) load with firth=None
    meta_path = pathlib.Path(f"{path}.meta.json")
    meta = json.loads(meta_path.read_text())
    del meta["genes"]["GENE_A"]["firth"]
    meta["schema_version"] = 1
    del meta["missing_policy"]
    meta_path.write_text(json.dumps(meta))
    assert load_artifact(path).genes["GENE_A"].firth is None


def test_emit_stores_firth_fit_with_run_weights(tmp_path):
    """Emission fits Firth on the ALT-coded, QC-retained matrix with the run's weights."""
    from tests.unit.test_meta_emit import _builder
    from variantcentrifuge.association.meta.artifact import emit_stratum_artifact

    b, alt, pheno = _builder()
    rng = np.random.default_rng(7)
    b.covariate_matrix = rng.normal(size=(len(pheno), 1))
    kw = {
        "stratum_id": "s1",
        "genome_build": "GRCh38",
        "annotation_version": "a",
        "mask_id": "m",
        "mask_hash": "h",
    }
    art = emit_stratum_artifact(str(tmp_path / "s1"), {"G1": b}, firth_weights="beta:1,25", **kw)
    expected = stratum_firth_burden(
        alt, alt.sum(axis=0), np.full(2, 2.0 * len(pheno)), pheno, b.covariate_matrix, "beta:1,25"
    )
    assert art.genes["G1"].firth == pytest.approx(expected)
    assert load_artifact(str(tmp_path / "s1")).genes["G1"].firth == pytest.approx(expected)
    # without firth_weights (direct API default) nothing is fitted
    art = emit_stratum_artifact(str(tmp_path / "s2"), {"G1": b}, **kw)
    assert art.genes["G1"].firth is None


def test_stratum_firth_none_when_fit_fails(monkeypatch):
    import variantcentrifuge.association.meta.firth_meta as fm

    monkeypatch.setattr(fm, "_firth_logistic", lambda *a, **k: None)
    geno = np.zeros((20, 1))
    geno[:3, 0] = 1
    y = np.array([0.0, 1.0] * 10)
    assert (
        fm.stratum_firth_burden(geno, np.array([3.0]), np.array([40.0]), y, None, "uniform") is None
    )


def test_stratum_firth_uses_tol_1e7(monkeypatch):
    """tol 1e-8 can make step-halving fail on rounding noise in large strata (review G)."""
    import variantcentrifuge.association.meta.firth_meta as fm

    seen = {}
    monkeypatch.setattr(fm, "_firth_logistic", lambda *a, **k: seen.update(k))
    geno = np.zeros((20, 1))
    geno[:3, 0] = 1
    fm.stratum_firth_burden(geno, np.array([3.0]), np.array([40.0]), np.ones(20), None, "uniform")
    assert seen["tol"] == 1e-7  # 1e-6 misses R parity at rel 1e-6 (s3/G2)


def test_runner_warns_once_per_stratum_for_excluded_firth_fits(caplog):
    keys = {"GENE_A": KEYS["GENE_A"], "GENE_B": ["2:100:A:G", "2:200:C:T"]}
    ok = _with_firth(make_stratum("s1", 1, gene_keys=keys), 0.5, 0.25)
    other = _with_firth(make_stratum("s2", 2, gene_keys=keys), 0.2, 0.5, weights="uniform")
    partial = _with_firth(make_stratum("s3", 3, gene_keys=keys), 0.9, 0.3)
    partial.genes["GENE_B"].firth = None
    with caplog.at_level("WARNING", logger="variantcentrifuge"):
        meta_analyze_artifacts([ok, other, partial])
    msgs = [r.getMessage() for r in caplog.records if "Firth" in r.getMessage()]
    assert len(msgs) == 2
    assert any("s2" in m and "2 of 2" in m and "uniform" in m for m in msgs)
    assert any("s3" in m and "1 of 2" in m and "no Firth fit" in m for m in msgs)
