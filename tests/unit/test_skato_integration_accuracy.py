"""SKAT-O omnibus integral accuracy (issue #181).

``tests/fixtures/skato_integration/params.json`` holds SKAT-O integration inputs
(pmin_q, tau, mixture moments, lambdas, pmin) captured from seeded synthetic
binary data sets. Each case is checked against an independent adaptive integral
(scipy.integrate.quad in t = sqrt(x) with explicit kink/support breakpoints and
tight tolerances; exact chi2(1) CDF for rank one, Davies at acc=1e-9 otherwise),
never against values produced by the implementation under test.

End-to-end references come from R SKAT 2.2.5 (``SKAT_Null_Model(y ~ x,
out_type="D", Adjustment=FALSE)``, ``SKAT(g, obj, weights=w,
method="optimal.adj")``) on the matrices regenerated below. Where R's own
``integrate`` misses mass, the reference is R's own integrand
(``SKAT:::SKAT_Optimal_Integrate_Func_Liu``, exact for a rank-one conditional)
integrated with ``integrate`` piecewise at the same breakpoints.
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest
import scipy.integrate
import scipy.stats

from variantcentrifuge.association.backends import davies as dv
from variantcentrifuge.association.backends import skat_math as sm
from variantcentrifuge.association.backends.python_backend import PythonSKATBackend
from variantcentrifuge.association.weights import beta_maf_weights

pytestmark = pytest.mark.skipif(not dv._try_load_davies(), reason="Davies C extension required")

PARAMS = json.loads(
    (pathlib.Path(__file__).parent.parent / "fixtures/skato_integration/params.json").read_text()
)


def _param(case: dict) -> dict:
    return dict(case, tau=np.array(case["tau"]), lambdas=np.array(case["lambdas"]))


def _integrate(case: dict, fn=sm._skato_integrate_davies) -> float:
    return fn(np.array(case["pmin_q"]), _param(case), case["rho"], case["pmin"])


def _reference(case: dict) -> float:
    """Independent adaptive SKAT-O integral (R: SKAT_Optimal_PValue_Davies)."""
    rho = np.array(case["rho"])
    a, b = np.array(case["pmin_q"]) / (1 - rho), np.array(case["tau"]) / (1 - rho)
    lam, mu = np.array(case["lambdas"]), case["mu_q"]
    sd = np.sqrt((case["var_q"] - case["var_remain"]) / case["var_q"])
    cand = [(a[i] - (mu - mu / sd)) / b[i] for i in range(len(a))]
    cand += [(a[i] - a[j]) / (b[i] - b[j]) for i in range(len(a)) for j in range(i)]
    pts = np.sqrt(np.unique([z for z in cand if np.isfinite(z) and 0 < z < 40]))

    def f(t: float) -> float:
        z = (np.min(a - b * t * t) - mu) * sd + mu
        if z <= 0:
            return 0.0
        if len(lam) == 1:
            cdf = scipy.stats.chi2.cdf(z / lam[0], 1)
        else:
            sf, fault = dv.davies_pvalue(float(z), lam, acc=1e-9, lim=10_000_000)
            assert fault == 0
            cdf = 1 - sf
        return cdf * np.sqrt(2 / np.pi) * np.exp(-t * t / 2)

    v, err = scipy.integrate.quad(
        f, 0, np.sqrt(40), points=pts, epsabs=1e-13, epsrel=1e-12, limit=4000
    )
    assert err < 1e-10
    return min(1 - v, case["pmin"] * len(rho))


@pytest.mark.unit
@pytest.mark.parametrize(
    ("name", "tol"),
    [
        # Davies acc=1e-6 in the integrand bounds the multi-eigenvalue error.
        ("support_edge_k11", 1e-7),
        ("multi_k8", 1e-7),
        ("rank_one_moderate", 1e-10),
        ("rank_one_narrow", 1e-10),
        ("severe_imbalance", 1e-10),
    ],
)
def test_integral_matches_independent_adaptive(name: str, tol: float) -> None:
    """Support-edge splits: baseline split only at rho kinks (off by 8e-5 .. 1.5e-2)."""
    case = PARAMS[name]
    assert _integrate(case) == pytest.approx(_reference(case), abs=tol)


@pytest.mark.unit
def test_support_edge_case_needs_the_edge_split() -> None:
    """Captured binary k=11 case: rho-kink-only rule gave 0.3541456 vs 0.3540632."""
    case = PARAMS["support_edge_k11"]
    rho = np.array(case["rho"])
    a, b = np.array(case["pmin_q"]) / (1 - rho), np.array(case["tau"]) / (1 - rho)
    mu, sd = case["mu_q"], np.sqrt((case["var_q"] - case["var_remain"]) / case["var_q"])
    edge = np.min((a - (mu - mu / sd)) / b)
    assert 0 < edge < 40
    assert np.isclose(sm._skato_breakpoints(a, b, np.array([mu - mu / sd])), edge).any()
    assert _integrate(case) == pytest.approx(0.3540631778, abs=1e-7)


@pytest.mark.unit
def test_rank_one_does_not_fall_back_to_liu(monkeypatch: pytest.MonkeyPatch) -> None:
    """Davies cannot evaluate one eigenvalue near the support edge; use the exact chi2(1)."""
    called = []
    monkeypatch.setattr(sm, "_skato_integrate_liu", lambda *a: called.append(a) or 0.0)
    p = _integrate(PARAMS["rank_one_moderate"])
    assert not called
    assert p == pytest.approx(0.220218733, abs=1e-9)


@pytest.mark.unit
@pytest.mark.parametrize("name", ["rank_one_moderate", "rank_one_narrow", "severe_imbalance"])
def test_liu_fallback_keeps_capped_burden_rho(name: str) -> None:
    """For a rank-one conditional the Liu integrand is exact, so it must match the reference.

    The baseline dropped rho=0.999 (``rho < 0.999``) and returned 0.2057 instead of 0.2202.
    """
    case = PARAMS[name]
    assert _integrate(case, sm._skato_integrate_liu) == pytest.approx(_reference(case), abs=1e-9)


@pytest.mark.unit
def test_severe_imbalance_resolution() -> None:
    """Two cases in 60 samples: baseline 0.8579 (< pmin, impossible); R SKAT prints 1.

    R's 1 is numerical: its per-rho p-values come from Liu after Davies ifault (Liu
    returns 1 here; exact Davies gives 0.934..0.980) and ``integrate`` never samples the
    width-1e-3 interval at x ~ 0 that carries the mass. The SKAT-O p-value can never be
    below min_rho p_rho, and the integral matches the independent reference.
    """
    case = PARAMS["severe_imbalance"]
    p = _integrate(case)
    assert case["pmin"] <= p <= 1.0
    assert p == pytest.approx(0.95593010347, abs=1e-9)


@pytest.mark.unit
def test_compute_pvalue_retries_davies_before_approximation() -> None:
    """Davies lim=10_000 fails (ifault=1) for small q with few eigenvalues; R then uses Liu."""
    lam, q = np.array([2.0, 1.0]), 1.0
    assert dv.davies_pvalue(q, lam, acc=1e-6, lim=10_000)[1] == 1
    # Exact P[2 X + Y > 1] for X, Y ~ chi2(1), by one-dimensional convolution.
    inner = scipy.integrate.quad(
        lambda x: scipy.stats.chi2.pdf(x, 1) * scipy.stats.chi2.sf(q - 2 * x, 1),
        0,
        q / 2,
        epsabs=1e-13,
    )[0]
    exact = inner + scipy.stats.chi2.sf(q / 2, 1)
    p, method, converged = dv.compute_pvalue(q, lam)
    assert (method, converged) == ("davies", True)
    assert p == pytest.approx(exact, abs=1e-6)


def _lowrank(n: int, cases: int, signal: int, seed: int):
    """Seeded two-variant binary data (same generator as the R comparison)."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 1))
    y = np.array([1.0] * cases + [0.0] * (n - cases))
    rng.shuffle(y)
    g = rng.binomial(2, [0.05, 0.12], size=(n, 2)).astype(float)
    if signal:
        g[y == 1, :] = rng.binomial(2, [0.3, 0.4], size=(cases, 2))
    return g, x, y, beta_maf_weights(g.mean(0) / 2, a=1, b=25)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("spec", "r_ref"),
    [
        # R SKAT p-value (Davies integrand, no fallback)
        ((300, 10, 0, 51), 0.329416489706),
        # R falls back to Liu (exact for rank one); R integrand split at breakpoints
        ((300, 150, 0, 51), 0.220218733),
        ((60, 30, 1, 51), 0.001997295676),
    ],
)
def test_skato_low_rank_matches_r(spec: tuple[int, int, int, int], r_ref: float) -> None:
    g, x, y, w = _lowrank(*spec)
    backend = PythonSKATBackend()
    null = backend.fit_null_model(y, x, "binary")
    p = backend.test_gene("G", g, null, "SKATO", weights=w)["p_value"]
    assert p == pytest.approx(r_ref, rel=1e-6)


@pytest.mark.unit
def test_davies_failure_in_integrand_uses_liu_like_r(monkeypatch: pytest.MonkeyPatch) -> None:
    """R SKAT_Optimal_PValue_Davies: any Davies ifault -> SKAT_Optimal_PValue_Liu."""
    case = PARAMS["multi_k8"]
    monkeypatch.setattr(sm, "davies_pvalue", lambda *a, **k: (2.0, 1))
    assert _integrate(case) == _integrate(case, sm._skato_integrate_liu)


@pytest.mark.unit
def test_no_rho_below_one_returns_bonferroni_bound() -> None:
    case = dict(PARAMS["multi_k8"], rho=[1.0] * 7, pmin=0.01)
    assert _integrate(case) == pytest.approx(0.07)
