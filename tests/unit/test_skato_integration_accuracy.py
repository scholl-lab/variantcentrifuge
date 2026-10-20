"""SKAT-O omnibus integral accuracy (issue #181).

``tests/fixtures/skato_integration/params.json`` holds SKAT-O integration inputs
(pmin_q, tau, mixture moments, lambdas, pmin) captured from seeded synthetic
binary data sets. Each case is checked against an independent integral
p = int_0^edge f_chi2_1(x) P[Q > g(x)] dx + P[chi2_1 > edge] (adaptive quad split
at the rho kinks, no x <= 40 truncation, no pmin cap; exact survival for one or
two eigenvalues, Davies at acc=1e-9 otherwise), never against values produced by
the implementation under test.

End-to-end references come from R SKAT 2.2.5 (``SKAT_Null_Model(y ~ x,
out_type="D", Adjustment=FALSE)``, ``SKAT(g, obj, weights=w,
method="optimal.adj")``) on the matrices regenerated below. Where R's own
``integrate`` misses mass, the reference is R's own integrand
(``SKAT:::SKAT_Optimal_Integrate_Func_Liu``, exact for a rank-one conditional)
integrated with ``integrate`` piecewise at the same breakpoints.
"""

from __future__ import annotations

import itertools
import json
import pathlib

import numpy as np
import pytest
import scipy.integrate
import scipy.stats

from variantcentrifuge.association.backends import davies as dv
from variantcentrifuge.association.backends import skat_math as sm
from variantcentrifuge.association.backends import skato_integrate
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


def _surv(z: float, lam: np.ndarray) -> float:
    """P[sum(lam * chi2_1) > z]: exact for 1-2 eigenvalues, Davies acc=1e-9 otherwise."""
    if z <= 0:
        return 1.0
    if len(lam) == 1:
        return float(scipy.stats.chi2.sf(z / lam[0], 1))
    if len(lam) == 2:
        # Polar coordinates: r^2 ~ chi2_2 (P[r^2 > s] = exp(-s/2)) independent of theta
        def ang(t: float) -> float:
            return float(np.exp(-z / (2 * (lam[0] * np.cos(t) ** 2 + lam[1] * np.sin(t) ** 2))))

        return scipy.integrate.quad(ang, 0, np.pi, epsabs=0, epsrel=1e-13)[0] / np.pi
    sf, fault = dv.davies_pvalue(float(z), lam, acc=1e-9, lim=10_000_000)
    assert fault == 0
    return min(max(sf, 0.0), 1.0)


def _reference(case: dict, surv=None) -> float:
    """Independent SKAT-O integral p = int_0^edge f1(x) S(g(x)) dx + P[chi2_1 > edge].

    Adaptive quad split at all rho kinks; no x <= 40 truncation and no pmin cap, so it
    shares neither with production (R: 1 - int_0^40 f1 (1 - S), capped at n_rho * pmin).
    """
    rho = np.array(case["rho"])
    a, b = np.array(case["pmin_q"]) / (1 - rho), np.array(case["tau"]) / (1 - rho)
    lam, mu = np.array(case["lambdas"]), case["mu_q"]
    sd = np.sqrt((case["var_q"] - case["var_remain"]) / case["var_q"])
    edge = min(float(np.min((a - (mu - mu / sd)) / b)), 1500.0)
    cand = [(a[i] - a[j]) / (b[i] - b[j]) for i in range(len(a)) for j in range(i)]
    cuts = [0.0, *sorted(c for c in cand if np.isfinite(c) and 0 < c < edge), edge]
    total, eps = 0.0, 1e-10
    for lo, hi in itertools.pairwise(cuts):
        v, err = scipy.integrate.quad(
            lambda x: (
                scipy.stats.chi2.pdf(x, 1)
                * (surv or _surv)((np.min(a - b * x) - mu) * sd + mu, lam)
            ),
            lo,
            hi,
            epsabs=0,
            epsrel=eps,
            limit=1000,
        )
        assert err <= 10 * eps * max(v, 1e-300)
        total += v
    return total + float(scipy.stats.chi2.sf(edge, 1))


@pytest.mark.unit
@pytest.mark.parametrize(
    "name",
    ["support_edge_k11", "multi_k8", "rank_one_moderate", "rank_one_narrow", "severe_imbalance"],
)
def test_integral_matches_independent_adaptive(name: str) -> None:
    """Support-edge splits: baseline split only at rho kinks (off by 8e-5 .. 1.5e-2)."""
    case = PARAMS[name]
    ref = _reference(case)
    assert ref < case["pmin"] * len(case["rho"])  # the shared cap cannot mask an error
    # budget 1e-4 relative (adaptive estimate); achieved errors are ~1e-8
    assert _integrate(case) == pytest.approx(ref, rel=1e-6, abs=0.0)


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
    monkeypatch.setattr(skato_integrate, "_skato_integrate_liu", lambda *a: called.append(a) or 0.0)
    p = _integrate(PARAMS["rank_one_moderate"])
    assert not called
    assert p == pytest.approx(0.220218733, rel=1e-6, abs=0.0)


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
    assert p == pytest.approx(0.95593010347, rel=1e-6, abs=0.0)


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
    assert p == pytest.approx(r_ref, rel=1e-6, abs=0.0)


@pytest.mark.unit
def test_davies_failure_in_integrand_uses_liu_like_r(monkeypatch: pytest.MonkeyPatch) -> None:
    """R SKAT_Optimal_PValue_Davies: any Davies ifault -> SKAT_Optimal_PValue_Liu."""
    case = PARAMS["multi_k8"]
    monkeypatch.setattr(skato_integrate, "davies_pvalue", lambda *a, **k: (2.0, 1))
    assert _integrate(case) == _integrate(case, sm._skato_integrate_liu)


@pytest.mark.unit
def test_no_rho_below_one_returns_bonferroni_bound() -> None:
    case = dict(PARAMS["multi_k8"], rho=[1.0] * 7, pmin=0.01)
    assert _integrate(case) == pytest.approx(0.07, abs=0.0)


def _tail_case(q_rho: np.ndarray, z1: np.ndarray, monkeypatch: pytest.MonkeyPatch):
    """Run _skato_get_pvalue and capture the integration inputs it builds."""
    seen: dict = {}
    real = sm._skato_integrate_davies

    def spy(pq, param, rho, pmin):
        seen.update(pmin_q=list(pq), rho=list(rho), pmin=pmin, **param)
        seen.update(tau=list(param["tau"]), lambdas=list(param["lambdas"]))
        return real(pq, param, rho, pmin)

    monkeypatch.setattr(skato_integrate, "_skato_integrate_davies", spy)
    p, _each = sm._skato_get_pvalue(q_rho, z1, sm._SKATO_RHO_GRID)
    return p, seen


_RHO_CAPPED = np.array([min(r, 0.999) for r in sm._SKATO_RHO_GRID])
_B3 = np.array([[1, 1, 1], [1, -1, 0], [1, 1, -2]], dtype=np.float64)
_B3 /= np.linalg.norm(_B3, axis=1)[:, None]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("z1", "k", "qs"),
    [
        # rank-one conditional: exact chi2(1) survival
        (np.eye(2), 1, [8.0, 20.0, 35.0, 49.0, 60.0]),
        # two eigenvalues: Davies' absolute floor (1e-9) -> tilted tail survival
        (np.diag(np.sqrt([0.01, 0.1, 1.0])) @ _B3, 2, [15.0, 33.0, 42.0, 45.5, 46.0, 50.0]),
    ],
)
def test_statistic_sweep_into_extreme_tail(
    z1: np.ndarray, k: int, qs: list[float], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Review #181: relative accuracy, monotonicity and bounds from moderate p to 1e-12.

    Baseline for q=49 (k=1) / q=33 (k=2): 1.682e-10 / 3.229e-8 (1 - int_0^40 truncation
    and cancellation). f03853b (Davies acc floor 1e-9) at q=45.5/46/50 (k=2):
    3.149e-11 / 7.783e-11 / 1.833e-12, non-monotone. References (no truncation, no cap,
    exact survival, per-rho inputs of 574d36e): 8.844e-11, 2.2228e-8, 3.7739e-11, 2.9265e-11,
    3.8338e-12. R SKAT
    2.2.5 (SKAT_Optimal_Get_Pvalue) gives 1.603e-10 (q=49) and 2.801e-8 (q=33): its
    integrate(0, 40) truncates, so in this tail we deliberately differ from R.
    """
    prev = 1.0
    for q in qs:
        p, case = _tail_case(q * (1 - np.array(sm._SKATO_RHO_GRID)), z1, monkeypatch)
        assert len(case["lambdas"]) == k
        n_rho_pmin = case["pmin"] * len(case["rho"])
        assert case["pmin"] <= p <= min(1.0, n_rho_pmin)
        assert p <= prev
        prev = p
        assert p == pytest.approx(_reference(case), rel=1e-6, abs=0.0)


@pytest.mark.unit
@pytest.mark.parametrize("lam", [[1.0, 0.1], [1.0, 0.3, 0.05]])
def test_tilted_survival_matches_exact(lam: list[float]) -> None:
    """Exponentially tilted Davies keeps relative accuracy where plain Davies is 0."""
    lam_arr = np.array(lam)

    def exact(z: float) -> float:  # k = 3: integrate the exact k = 2 survival
        if len(lam) == 2:
            return _surv(z, lam_arr)
        a = scipy.integrate.quad(
            lambda u: scipy.stats.chi2.pdf(u, 1) * _surv(z - lam[2] * u, lam_arr[:2]),
            0,
            z / lam[2],
            epsabs=0,
            epsrel=1e-11,
        )[0]
        return float(scipy.stats.chi2.sf(z / lam[2], 1) + a)

    for z in (5.0, 20.0, 40.0, 60.0):
        ref = exact(z)
        assert skato_integrate._tilted_sf(z, lam_arr) == pytest.approx(ref, rel=1e-6, abs=0.0)
    assert dv.davies_pvalue(60.0, lam_arr, acc=1e-9, lim=1_000_000)[0] < 1e-9


@pytest.mark.unit
def test_pmin_quantile_survives_tiny_pmin() -> None:
    """chi2.ppf(1 - p) is inf for p < 1e-16 (R: qchisq(1 - pmin)); use isf."""
    lam = [np.array([1.0])] * 7
    pmin, _pval, pmin_q = sm._skato_each_q({}, np.full((1, 7), 90.0), list(_RHO_CAPPED), lam)
    assert pmin[0] < 1e-16 and np.all(np.isfinite(pmin_q))
    assert pmin_q[0, 0] == pytest.approx(90.0, rel=1e-6, abs=0.0)


@pytest.mark.unit
def test_tilted_survival_davies_failure_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """A Davies fault inside the tilted evaluation sends the integral to Liu (R policy)."""
    monkeypatch.setattr(dv, "davies_pvalue", lambda *a, **k: (2.0, 1))
    with pytest.raises(skato_integrate._DaviesFailureError):
        skato_integrate._tilted_sf(20.0, np.array([1.0, 0.3, 0.05]))


@pytest.mark.unit
def test_large_gene_outer_integral_converges(monkeypatch: pytest.MonkeyPatch) -> None:
    """501 variants with an exactly chi2(500) conditional: fixed 16-node rule was +0.54%.

    Reference: exact chi2(500) survival, adaptive quad without truncation or cap
    (with the per-rho p-values of 574d36e: 5.8764202692e-8; 64-node rule 5.8764202648e-8).
    """
    n = 501
    v = np.eye(1, n)[0] - np.ones(n) / np.sqrt(n)
    v /= np.linalg.norm(v)
    z1 = np.eye(n) - 2 * np.outer(v, v)
    z1[0] *= 10
    p, case = _tail_case(3500 * (1 - np.array(sm._SKATO_RHO_GRID)), z1, monkeypatch)
    lam = np.array(case["lambdas"])
    assert len(lam) == 500 and np.ptp(lam) < 1e-10 * lam[0]

    def chi2_500(z: float, lam: np.ndarray) -> float:
        return 1.0 if z <= 0 else float(scipy.stats.chi2.sf(z / lam.mean(), len(lam)))

    assert p == pytest.approx(_reference(case, chi2_500), rel=1e-4, abs=0.0)


@pytest.mark.unit
@pytest.mark.parametrize("c", [1e-12, 1.0, 1e12])
def test_scale_invariance(c: float) -> None:
    """Q and the eigenvalues scale together; p-values must not change (brentq xtol was absolute)."""
    lam = np.array([1.0, 0.1])
    assert skato_integrate._tilted_sf(60 * c, lam * c) == pytest.approx(
        1.0007853115258435e-14, rel=1e-9, abs=0.0
    )
    lam3 = np.array([1.0, 0.1, 0.3])
    assert dv._kuonen_pvalue(60 * c, lam3 * c) == pytest.approx(
        1.3419078340383e-14, rel=1e-9, abs=0.0
    )
    z1 = np.diag(np.sqrt([0.01, 0.1, 1.0])) @ _B3
    grid = np.array(sm._SKATO_RHO_GRID)
    for q in (5.0, 46.0):
        p_ref, _ = sm._skato_get_pvalue(q * (1 - grid), z1, sm._SKATO_RHO_GRID)
        p, _each = sm._skato_get_pvalue(c * q * (1 - grid), np.sqrt(c) * z1, sm._SKATO_RHO_GRID)
        assert p == pytest.approx(p_ref, rel=1e-9, abs=0.0)


@pytest.mark.unit
def test_subdivision_cap_falls_back_to_liu_and_reports_it(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """Outer quadrature over budget -> whole-integral Liu (as for Davies faults), reported."""
    g, x, y, w = _lowrank(300, 150, 0, 51)
    backend = PythonSKATBackend()
    null = backend.fit_null_model(y, x, "binary")
    ok = backend.test_gene("G", g, null, "SKATO", weights=w)
    assert ok["p_converged"] is True
    monkeypatch.setattr(skato_integrate, "_SKATO_QUAD_LIMIT", 1)
    monkeypatch.setattr(skato_integrate, "_SKATO_DAVIES_REL_ACC", 1e-14)
    with caplog.at_level("WARNING", logger="variantcentrifuge"):
        res = backend.test_gene("G", g, null, "SKATO", weights=w)
    assert "SKAT-O integral piece" in caplog.text
    assert res["p_converged"] is False
    # rank-one conditional: the Liu integrand is exact, so the value is unchanged
    assert res["p_value"] == pytest.approx(ok["p_value"], rel=1e-6, abs=0.0)
    assert skato_integrate._skato_last_converged() is False


@pytest.mark.unit
def test_davies_threadsafe() -> None:
    """qfc uses C globals and CFFI releases the GIL: concurrent calls must not interfere."""
    from concurrent.futures import ThreadPoolExecutor

    args = [(5.0 + i % 4, np.array([1.0, 0.1, 0.01])) for i in range(64)]
    seq = [dv.davies_pvalue(q, lam) for q, lam in args]
    with ThreadPoolExecutor(4) as pool:
        par = list(pool.map(lambda a: dv.davies_pvalue(*a), args))
    assert par == seq


@pytest.mark.unit
def test_tail_pvalue_uses_tilted_davies_not_saddlepoint() -> None:
    """Real-data SKAT (public 1000G PKHD1, 35 eigenvalues): p ~ 1.8e-6 is below the proactive
    threshold. Saddlepoint gave 1.90661e-6 (+6.4%); R SKAT 2.2.5 Davies acc=1e-6 gives
    1.78656e-6 (absolute accuracy). References: Davies acc=1e-12, lim=1e7 1.7922060078e-6;
    importance sampling under the tilted law 1.7900e-6 +- 0.0050e-6.
    """
    d = json.loads(
        (
            pathlib.Path(__file__).parent.parent
            / "fixtures/skato_integration/pkhd1_skat_eigen.json"
        ).read_text()
    )
    p, method, converged = dv.compute_pvalue(d["q"], np.array(d["lambdas"]))
    assert (method, converged) == ("davies_tilted", True)
    assert p == pytest.approx(1.7922060078e-6, rel=1e-5, abs=0.0)
