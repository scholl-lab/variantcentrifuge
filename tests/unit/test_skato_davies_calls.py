"""Davies call counts for SKAT-O (0.25 perf fixes 1 + 2): fewer qfc calls, same p-values."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import chi2

from variantcentrifuge.association.backends import davies
from variantcentrifuge.association.backends.davies import compute_pvalue, davies_pvalue
from variantcentrifuge.association.backends.skato_integrate import _skato_each_q

pytestmark = [
    pytest.mark.unit,
    pytest.mark.skipif(not davies._try_load_davies(), reason="Davies C extension not built"),
]

_RHOS = [0.0, 0.01, 0.04, 0.09, 0.25, 0.5, 0.999]
# rho < 0.999: spread spectrum; rho = 0.999: one dominant (burden) eigenvalue
_LAM_SPREAD = np.array([3.0, 2.0, 1.0, 0.5, 0.2])
_LAM_DOMINANT = np.array([40.0, 0.03, 0.01])


class _CountingLib:
    def __init__(self, lib):
        self.lib, self.calls = lib, 0

    def qfc(self, *args):
        self.calls += 1
        return self.lib.qfc(*args)


@pytest.fixture
def qfc_counter(monkeypatch):
    counter = _CountingLib(davies._qfc_lib)
    monkeypatch.setattr(davies, "_qfc_lib", counter)
    return counter


def _calls(counter, fn, *args):
    before = counter.calls
    out = fn(*args)
    return out, counter.calls - before


def test_compute_pvalue_one_qfc_call_where_lim_1e4_hits_term_limit(qfc_counter):
    q = 5.0  # small q under a dominant eigenvalue: lim=10_000 returns ifault=1
    assert davies_pvalue(q, _LAM_DOMINANT, acc=1e-6, lim=10_000)[1] == 1
    p_ref, f_ref = davies_pvalue(q, _LAM_DOMINANT, acc=1e-6, lim=1_000_000)
    assert f_ref == 0
    (p, method, conv), n = _calls(qfc_counter, compute_pvalue, q, _LAM_DOMINANT)
    assert (p, method, conv) == (p_ref, "davies", True)
    assert n == 1


def _q_all(q_spread: float, q_dominant: float) -> np.ndarray:
    return np.array([[q_spread] * 6 + [q_dominant]])


def _lambda_all() -> list[np.ndarray]:
    return [_LAM_SPREAD] * 6 + [_LAM_DOMINANT]


def test_rho_0999_davies_skipped_when_lower_bound_exceeds_min(qfc_counter):
    q_all = _q_all(25.0, 30.0)  # p ~ 1e-4 for rho < 0.999; P[40 chi2_1 > 30] ~ 0.39
    lam = _lambda_all()
    assert chi2.sf(30.0 / 40.0, 1) > compute_pvalue(25.0, _LAM_SPREAD)[0] + 1e-6
    ref, ref_calls = _calls(
        qfc_counter, lambda: [compute_pvalue(q_all[0, i], lam[i]) for i in range(6)]
    )
    (pmin, pval, _), n = _calls(qfc_counter, _skato_each_q, {}, q_all, _RHOS, lam)
    assert n == ref_calls  # rho = 0.999 made no qfc call
    assert list(pval[0, :6]) == [r[0] for r in ref]
    assert pmin[0] == min(r[0] for r in ref)
    assert int(np.argmin(pval[0])) < 6


def test_rho_0999_davies_kept_when_bound_does_not_exceed_min(qfc_counter):
    q_all = _q_all(25.0, 400.0)  # rho = 0.999 is the minimum
    lam = _lambda_all()
    ref, ref_calls = _calls(
        qfc_counter, lambda: [compute_pvalue(q_all[0, i], lam[i]) for i in range(7)]
    )
    (pmin, pval, _), n = _calls(qfc_counter, _skato_each_q, {}, q_all, _RHOS, lam)
    assert n == ref_calls
    assert list(pval[0]) == [r[0] for r in ref]
    assert pmin[0] == pval[0, 6]


def test_lambda_max_chi2_bound_is_a_lower_bound():
    """P[lmax chi2_1 > q] <= P[sum(l chi2_1) > q]: closed form (ties) and Davies (random)."""
    for k in (2, 3, 10):  # equal eigenvalues: P[l chi2_k > q] in closed form
        for q in (0.1, 1.0, 5.0, 30.0):
            assert chi2.sf(q, 1) <= chi2.sf(q, k)
    rng = np.random.default_rng(20251021)
    for _ in range(300):
        lam = rng.exponential(1.0, rng.integers(2, 30))
        if rng.random() < 0.3:  # near-ties at the top
            lam[: rng.integers(2, len(lam) + 1)] = lam.max() * (1 + 1e-9 * rng.random())
        if rng.random() < 0.5:  # dominant eigenvalue, as at rho = 0.999
            lam[0] *= rng.uniform(10, 1000)
        q = float(rng.uniform(0.01, 4.0) * lam.sum())
        p, f = davies_pvalue(q, lam, acc=1e-9, lim=1_000_000)
        if f == 0:
            assert chi2.sf(q / lam.max(), 1) <= p + 1e-9


def test_meta_skato_underflowed_rhos_keep_v024_p():
    """Zero-recovery returns the smallest positive per-rho p: rho=0.999 must stay exact.

    Rhos 0..0.5 underflow to p = 0 here; the expected value is v0.24.0's output.
    """
    from variantcentrifuge.association.meta.combine import meta_skato

    u = [np.sqrt(5000.0) * np.array([1.0, -1.0])]
    assert meta_skato(u, [2.0 * np.eye(2)], np.ones(2)) == (0.11379218388721413, 0.0)
