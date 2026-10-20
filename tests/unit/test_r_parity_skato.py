"""SKAT-O and cohort meta-analysis parity with R SKAT 2.2.5 / MetaSKAT 0.90.

Reference values were computed once in R on the SYNTH-like COAST fixtures
(tests/fixtures/coast_golden), so CI needs no R:

- SKAT-O: ``SKAT(Z * w, SKAT_Null_Model(y ~ 1, out_type="D", Adjustment=FALSE),
  method="SKATO", weights=rep(1, k), is_check_genotype=FALSE)`` with the
  allelic-series weights ``w`` from ``_compute_allelic_skat_weights``.
- Meta: samples split into two cohorts (first/second half), covariate
  ``(i mod 5) - 2``; ``MetaSKAT_wZ(Z, Meta_Null_Model(..., out_type="D"),
  combined.weight=TRUE, weights.beta=c(1, 25))`` with r.corr=1 (burden),
  r.corr=0 (SKAT) and method="optimal" (SKAT-O).
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import pytest

from tests.unit.meta_helpers import NO_MISSING_POLICY
from variantcentrifuge.association.backends.coast_python import (
    PythonCOASTBackend,
    _compute_allelic_skat_weights,
)
from variantcentrifuge.association.backends.python_backend import PythonSKATBackend
from variantcentrifuge.association.meta.align import align_strata, pooled_frequency_weights
from variantcentrifuge.association.meta.artifact import build_stratum_artifact
from variantcentrifuge.association.meta.combine import meta_burden, meta_skat, meta_skato

FIXTURES = pathlib.Path(__file__).parent.parent / "fixtures" / "coast_golden"

R_SKATO = {
    "synth_g1_all_categories": 0.09329995,
    "synth_g2_partial_bmv_ptv": 0.8489676,
    "synth_g4_all_bmv": 0.6787618,
    "synth_g5_tied_scores": 0.6141558,
}

# (burden, SKAT, SKAT-O) from MetaSKAT_wZ
R_META = {
    "synth_g1_all_categories": (1.3761737103e-01, 3.1161925411e-01, 2.3617350562e-01),
    "synth_g2_partial_bmv_ptv": (7.6297853584e-01, 5.4366614891e-01, 7.3596155732e-01),
    "synth_g4_all_bmv": (4.8314675831e-01, 6.9478618897e-01, 6.9090889947e-01),
    "synth_g5_tied_scores": (5.6578349890e-01, 9.0507617940e-01, 7.7532596445e-01),
}


def _load(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    geno = pd.read_csv(FIXTURES / f"{name}_geno.csv").to_numpy(dtype=np.float64)
    pheno = pd.read_csv(FIXTURES / f"{name}_pheno.csv")["pheno"].to_numpy(dtype=np.float64)
    anno = pd.read_csv(FIXTURES / f"{name}_anno.csv")["anno"].to_numpy(dtype=np.int32)
    return geno, pheno, anno


@pytest.mark.parametrize("name", list(R_SKATO))
def test_skato_matches_r_skat(name: str) -> None:
    geno, pheno, anno = _load(name)
    backend = PythonSKATBackend()
    backend.detect_environment()
    coast = PythonCOASTBackend()
    coast._skat_backend = backend
    weights = _compute_allelic_skat_weights(geno, anno, [1.0, 2.0, 3.0])
    null = backend.fit_null_model(pheno, None, "binary")
    p = coast._run_allelic_skat(geno, weights, null)
    assert p == pytest.approx(R_SKATO[name], rel=2e-5)  # R printed to 7 significant digits


@pytest.mark.parametrize("name", list(R_META))
def test_meta_matches_metaskat(name: str) -> None:
    geno, pheno, _anno = _load(name)
    n = len(pheno)
    halves = (slice(0, n // 2), slice(n // 2, n))
    cov = ((np.arange(n) % 5) - 2.0)[:, None]
    keys = [f"1:{100 + j}:A:G" for j in range(geno.shape[1])]
    scores = [
        build_stratum_artifact(
            f"s{i}",
            {"G": (geno[h], geno[h].copy(), keys)},
            pheno[h],
            cov[h],
            "binary",
            genome_build="GRCh37",
            annotation_version="x",
            mask_id="m",
            mask_hash="m",
            n_cases=int(pheno[h].sum()),
            n_controls=int((1 - pheno[h]).sum()),
            missing_policy=NO_MISSING_POLICY,
        ).genes["G"]
        for i, h in enumerate(halves)
    ]
    u, v, ac, an, _keys = align_strata(scores)
    w = pooled_frequency_weights(ac, an)
    burden, skat, skato = R_META[name]
    assert meta_burden(u, v, w)[0] == pytest.approx(burden, rel=1e-6)
    assert meta_skat(u, v, w)[0] == pytest.approx(skat, rel=1e-4)
    assert meta_skato(u, v, w)[0] == pytest.approx(skato, rel=1e-3)


# Leave-one-stratum-out vs MetaSKAT 0.90 (#193). Samples split into three contiguous
# thirds [0, n//3), [n//3, 2n//3), [2n//3, n); covariate (i mod 5) - 2; dropping
# stratum d = MetaSKAT_wZ on the other two cohorts (Meta_Null_Model(split(y, coh),
# split(x, coh), out_type="D"), combined.weight=TRUE, weights.beta=c(1, 25)) with
# r.corr=1 (burden) and method="optimal" (SKAT-O). R script:
# scripts/generate_meta_loo_parity.R (not needed by CI).
# R SKAT 2.2.5 / MetaSKAT 0.90, values printed with %.10e. Full meta of all three
# strata in R: burden 1.4694995906e-01, SKAT-O 2.5123531457e-01 (g1).
R_META_LOO = {  # drop -> (burden, SKAT-O)
    "synth_g1_all_categories": {
        "s0": (2.3895780219e-01, 3.6888167586e-01),
        "s1": (5.5628031673e-01, 5.9624053551e-01),
        "s2": (1.9540758946e-02, 2.1222307002e-02),
    },
    "synth_g2_partial_bmv_ptv": {
        "s0": (5.8783323555e-01, 2.0415400546e-01),
        "s1": (8.1646649749e-01, 8.5148618302e-01),
        "s2": (3.4676349023e-01, 5.1331509752e-01),
    },
    "synth_g4_all_bmv": {
        "s0": (9.3927842992e-01, 6.3667919535e-01),
        "s1": (2.4664680070e-01, 3.9314885468e-01),
        "s2": (4.7527030293e-01, 6.8166110432e-01),
    },
}


@pytest.mark.parametrize("name", list(R_META_LOO))
def test_named_loo_matches_metaskat(name: str) -> None:
    from variantcentrifuge.association.meta.runner import meta_analyze_artifacts

    geno, pheno, _anno = _load(name)
    n = len(pheno)
    thirds = (slice(0, n // 3), slice(n // 3, 2 * n // 3), slice(2 * n // 3, n))
    cov = ((np.arange(n) % 5) - 2.0)[:, None]
    keys = [f"1:{100 + j}:A:G" for j in range(geno.shape[1])]
    arts = [
        build_stratum_artifact(
            f"s{i}",
            {"G": (geno[h], geno[h].copy(), keys)},
            pheno[h],
            cov[h],
            "binary",
            genome_build="GRCh37",
            annotation_version="x",
            mask_id="m",
            mask_hash="m",
            n_cases=int(pheno[h].sum()),
            n_controls=int((1 - pheno[h]).sum()),
            missing_policy=NO_MISSING_POLICY,
        )
        for i, h in enumerate(thirds)
    ]
    row = meta_analyze_artifacts(arts).iloc[0]
    burden = {
        s.split(":")[0]: float(s.split(":")[1]) for s in row["loo_burden_by_stratum"].split(";")
    }
    skato = {
        s.split(":")[0]: float(s.split(":")[1]) for s in row["loo_skato_by_stratum"].split(";")
    }
    for drop, (r_burden, r_skato) in R_META_LOO[name].items():
        assert burden[drop] == pytest.approx(r_burden, rel=1e-6)
        assert skato[drop] == pytest.approx(r_skato, rel=1e-3)
    worst = max(R_META_LOO[name], key=lambda d: R_META_LOO[name][d][0])
    assert row["loo_max_burden_stratum"] == worst
    assert row["loo_max_burden_p"] == pytest.approx(R_META_LOO[name][worst][0], rel=1e-6)
