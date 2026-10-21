"""Imbalance-robust binary SKAT kernel (association/spa.py) vs R SKATBinary_Robust (#223).

Fixtures and references come from ``scripts/generate_skat_robust_parity.R`` (R 4.5.3,
SKAT 2.2.5, SPAtest 3.1.2; provenance line at the top of each reference CSV). R fitted
``SKAT_Null_Model(y ~ x1 + x2, out_type = "D", Adjustment = FALSE)``; its fitted mean
``mu`` is stored in the cohort CSVs so both sides see identical inputs. Cohorts: 2000
samples at 1:9 and 900 samples at 8:1; genes with MAC <= 10, > 10, mixed and a single
variant, for ``SKATBinary_Robust(Z, obj, method = "SKAT" | "SKATO")``.

Tolerances:
- per-variant SPA/ER p-values and the calibrated (Score, Phi): deterministic, compared
  at rel 1e-9 (observed max 2e-12 / 2e-13);
- gene SKAT p-values: against R's Davies at acc = 1e-12 on R's own (Score, Phi), at the
  backend's Davies relative budget (rel 1e-4) or that absolute accuracy (2e-12);
- gene p-values as printed by SKATBinary_Robust: R runs qfc at acc = 1e-6 absolute (and
  falls back to Liu below it), so they are matched at abs 1e-6 / rel 1e-3 (the
  repository's SKAT-O convention; observed here: SKAT-O rel 5.7e-3 at abs 5.6e-7). The
  bound holds on these fixtures only: SKAT-O is
  evaluated by the repository's SKAT-O omnibus (equal to the default-path SKAT-O to 1e-12
  on identical (Score, Phi)), and R's ``SKAT_META_Optimal`` integration runs Davies at
  acc = 1e-6 inside the integrand. On other synthetic genes the review measured SKAT-O
  differences from R's printed p up to rel 2.9e-3 (abs 4.7e-4 at p = 0.33), the same on
  R's own captured (Score, Phi), i.e. from R's integration, not from the port.
"""

from __future__ import annotations

import csv
import pathlib

import numpy as np
import pytest
import scipy.stats

from variantcentrifuge.association.spa import (
    calibrated_variance,
    exact_pvalue,
    robust_moments,
    robust_skat,
    saddle_pvalues,
    variant_pvalues,
)

FIXTURES = pathlib.Path(__file__).parent.parent / "fixtures" / "skat_robust_parity"
GENES = ["le10", "gt10", "mixed", "single"]
COHORTS = ["r19", "r81"]


def _rows(name: str) -> list[dict[str, str]]:
    with open(FIXTURES / name) as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def _load(cohort: str, gene: str):
    data = np.genfromtxt(FIXTURES / f"{cohort}.csv", delimiter=",", names=True)
    cols = data.dtype.names
    assert cols is not None
    geno = np.column_stack([data[c] for c in cols if c.startswith(f"{gene}_")])
    y, mu = data["y"], data["mu"]
    design = np.column_stack([np.ones_like(y), data["x1"], data["x2"]])
    maf = geno.mean(axis=0) / 2.0
    weights = scipy.stats.beta.pdf(np.minimum(maf, 1.0 - maf), 1, 25)  # R weights.beta = c(1, 25)
    return geno, weights, y, mu, design


REF = {(r["cohort"], r["gene"], r["method"]): r for r in _rows("r_reference.csv")}
SINGLE: dict[tuple[str, str], list[float]] = {}
for _r in _rows("r_single.csv"):
    SINGLE.setdefault((_r["cohort"], _r["gene"]), []).append(float(_r["p_value"]))
MOMENTS: dict[tuple[str, str], list[tuple[float, list[float]]]] = {}
for _r in _rows("r_moments.csv"):
    MOMENTS.setdefault((_r["cohort"], _r["gene"]), []).append(
        (float(_r["score"]), [float(v) for v in _r["phi"].split()])
    )


@pytest.mark.parametrize("gene", GENES)
@pytest.mark.parametrize("cohort", COHORTS)
def test_variant_pvalues_match_r(cohort, gene):
    geno, _w, y, mu, design = _load(cohort, gene)
    _score, p, _g1, _v = variant_pvalues(geno, y, mu, design)
    np.testing.assert_allclose(p, SINGLE[(cohort, gene)], rtol=1e-9, atol=0)


@pytest.mark.parametrize("gene", GENES)
@pytest.mark.parametrize("cohort", COHORTS)
def test_calibrated_moments_match_r(cohort, gene):
    geno, w, y, mu, design = _load(cohort, gene)
    score, cov, _p_burden = robust_moments(geno, w, y, mu, design)
    assert cov is not None
    r_score = np.array([s for s, _ in MOMENTS[(cohort, gene)]])
    r_phi = np.array([row for _, row in MOMENTS[(cohort, gene)]])
    np.testing.assert_allclose(score * w, r_score, rtol=0, atol=1e-9 * np.abs(r_score).max())
    np.testing.assert_allclose(cov * np.outer(w, w), r_phi, rtol=0, atol=1e-9 * np.abs(r_phi).max())


@pytest.mark.parametrize("gene", GENES)
@pytest.mark.parametrize("cohort", COHORTS)
def test_robust_skat_matches_r_tight_davies(cohort, gene):
    geno, w, y, mu, design = _load(cohort, gene)
    p, rho = robust_skat(geno, w, y, mu, design, "SKAT")
    assert rho is None
    ref = float(REF[(cohort, gene, "SKAT")]["p_davies_1e12"])
    assert p == pytest.approx(ref, rel=1e-4, abs=2e-12)


@pytest.mark.parametrize("method", ["SKAT", "SKATO"])
@pytest.mark.parametrize("gene", GENES)
@pytest.mark.parametrize("cohort", COHORTS)
def test_robust_skat_matches_r_skatbinary_robust(cohort, gene, method):
    geno, w, y, mu, design = _load(cohort, gene)
    p, _rho = robust_skat(geno, w, y, mu, design, method)
    ref = float(REF[(cohort, gene, method)]["p_value"])
    assert p == pytest.approx(ref, rel=1e-3, abs=1e-6)


def test_single_variant_skato_is_skat():
    geno, w, y, mu, design = _load("r19", "single")
    assert robust_skat(geno, w, y, mu, design, "SKATO") == (
        robust_skat(geno, w, y, mu, design, "SKAT")[0],
        0.0,
    )


def test_exact_pvalue_is_hypergeometric_mid_p_without_covariates():
    # Constant mu: carrier case counts are hypergeometric; heterozygous carriers make
    # the statistic |#case carriers - k mu|, so the ER mid-p has a closed form.
    n, ncase, k = 50, 10, 6
    y = np.zeros(n)
    y[:ncase] = 1.0
    mu = np.full(n, ncase / n)
    geno = np.zeros(n)
    geno[[0, 1, 2, 20, 30, 40]] = 1.0  # 3 case carriers
    hyper = scipy.stats.hypergeom(n, ncase, k)
    dev = np.abs(np.arange(k + 1) - k * ncase / n)
    obs = dev[3]
    expected = hyper.pmf(np.arange(k + 1))[dev > obs + 1e-12].sum()
    expected += hyper.pmf(np.arange(k + 1))[np.isclose(dev, obs)].sum() / 2.0
    assert exact_pvalue(geno, y, mu) == pytest.approx(expected, rel=1e-12)


def test_saddle_pvalues_below_cutoff_is_normal_and_above_tracks_binomial():
    # g = carrier indicator, constant mu: q ~ Binomial(k, mu) exactly.
    n, k, mu0 = 400, 60, 0.05
    mu = np.full(n, mu0)
    g = np.zeros((n, 2))
    g[:k, :] = 1.0
    q = np.array([3.0, 11.0])  # z ~ 0.0 and ~ 4.7
    p = saddle_pvalues(q, mu, g, 2.0)
    assert p[0] == pytest.approx(scipy.stats.chi2.sf((3.0 - k * mu0) ** 2 / (k * mu0 * 0.95), 1))
    binom = scipy.stats.binom(k, mu0)
    mid_p = binom.sf(11) + binom.pmf(11) / 2.0  # the mirrored lower tail (q < 0) is empty
    normal = scipy.stats.chi2.sf((11.0 - k * mu0) ** 2 / (k * mu0 * 0.95), 1)
    assert normal < mid_p / 40  # the normal approximation is badly anti-conservative here
    assert p[1] == pytest.approx(mid_p, rel=0.2)


def test_saddle_pvalue_at_support_edge_uses_one_tail():
    # q at the largest attainable score: that tail's root is infinite (probability 0 by SPA)
    n = 200
    mu = np.full(n, 0.1)
    g = np.zeros((n, 1))
    g[:5, 0] = 1.0
    p = saddle_pvalues(np.array([5.0]), mu, g, 2.0)
    assert 0.0 <= p[0] < 1e-3


def test_calibrated_variance_edges():
    v = calibrated_variance(np.array([2.0, 2.0, 2.0]), np.array([0.0, 1.0, 0.05]))
    assert v[0] == pytest.approx(4.0 / 500.0)
    assert v[1] == 0.0
    assert v[2] == pytest.approx(4.0 / scipy.stats.chi2.isf(0.05, 1))


def test_rank_one_skato_falls_back_to_skat():
    # R Met_SKAT_Get_Pvalue: rank(Phi) <= 1 -> r.corr = 0 (SKAT)
    geno, w, y, mu, design = _load("r19", "single")
    geno2 = np.column_stack([geno[:, 0], geno[:, 0]])
    w2 = np.array([w[0], w[0]])
    p_skat = robust_skat(geno2, w2, y, mu, design, "SKAT")[0]
    assert robust_skat(geno2, w2, y, mu, design, "SKATO") == (p_skat, 0.0)


def test_cancelling_burden_returns_robust_burden_p():
    # Opposite single-carrier variants: burden score exactly 0, so the burden
    # recalibration is undefined (R errors there); the robust burden p (= 1) is returned.
    n = 40
    y = np.zeros(n)
    y[:20] = 1.0
    mu = np.full(n, 0.5)
    geno = np.zeros((n, 2))
    geno[0, 0] = 1.0  # one case carrier
    geno[30, 1] = 1.0  # one control carrier
    p, rho = robust_skat(geno, np.ones(2), y, mu, np.ones((n, 1)), "SKATO")
    assert (p, rho) == (1.0, None)


@pytest.mark.parametrize("seed", [None, 1, 2])
def test_opposing_variants_keep_their_signal(seed):
    """500:500, uniform weights, variant 1 in 50 cases only, variant 2 in 50 controls only.

    The burden score cancels to round-off (S ~ 1e-13, p_burden ~ 1 - 2e-14), so
    S^2 / qchisq(1 - p, 1) is noise. R SKAT 2.2.5 SKATBinary_Robust (intercept-only null,
    weights.beta = c(1, 1)) prints SKAT p 2.5e-21 in sample order and 1.3e-16 after a
    permutation (its own round-off: r = VarQ / VarQ_2 = 0.91 / 0.70); with the exact
    normal-branch variance (r = 1) the port gives 5.1e-23 in every order (v0.24 dev: 0.905
    and 0.968). The exact-zero case keeps R's robust burden p (test above).
    """
    n = 1000
    y = np.r_[np.ones(500), np.zeros(500)]
    geno = np.zeros((n, 2))
    geno[:50, 0] = 1.0
    geno[500:550, 1] = 1.0
    if seed is not None:
        order = np.random.default_rng(seed).permutation(n)
        y, geno = y[order], geno[order]
    mu = np.full(n, 0.5)
    p_skat, _ = robust_skat(geno, np.ones(2), y, mu, np.ones((n, 1)), "SKAT")
    p_skato, _ = robust_skat(geno, np.ones(2), y, mu, np.ones((n, 1)), "SKATO")
    assert p_skat == pytest.approx(5.1057684528e-23, rel=1e-6)
    assert p_skato < 1e-20


def test_unadjusted_variants_keep_the_exact_score_variance():
    # |z| <= 2: p~ is the normal p, so V~ = S^2 / qchisq(p~, 1) is the variance itself;
    # it is returned directly, without the qchisq round-off (large at a near-zero score).
    n = 400
    y = np.zeros(n)
    y[:40] = 1.0
    mu = np.full(n, 0.1)
    geno = np.zeros((n, 2))
    geno[[0, 50, 60, 70, 80, 90, 100, 110, 120, 130], 0] = 1.0  # score ~ 0
    geno[[1, 2, 3, 200, 210, 220, 230, 240, 250, 260, 270, 280], 1] = 1.0
    score, p, g1, v = variant_pvalues(geno, y, mu, np.ones((n, 1)))
    var = (g1**2 * (mu * (1 - mu))[:, None]).sum(axis=0)
    assert abs(score[0]) < 1e-12 and p[0] > 0.999  # 1 case carrier of 10 at mu = 0.1
    np.testing.assert_allclose(v, var, rtol=1e-13)  # old S^2/qchisq form: 19% off here
