"""qfc's central chi2_1 fast path (complex product, E2) against the original per-eigenvalue loop.

Reference: the same qfc with one extra term (lb=0, nc=1, n=1). That term adds exact zeros everywhere
(mean, sd, errbd, truncation, cfe, integrate) but is non-central, so qfc runs its original loop.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import chi2

from variantcentrifuge.association.backends import davies

pytestmark = [
    pytest.mark.unit,
    pytest.mark.skipif(not davies._try_load_davies(), reason="Davies C extension not built"),
]


def _qfc(lam, q, acc=1e-6, lim=1_000_000, nc=None, n=None):
    ffi, lib = davies._qfc_ffi, davies._qfc_lib
    lam = np.ascontiguousarray(lam, dtype=np.float64)
    r = len(lam)
    nc = np.zeros(r) if nc is None else np.ascontiguousarray(nc, dtype=np.float64)
    n = np.ones(r, dtype=np.int32) if n is None else np.ascontiguousarray(n, dtype=np.int32)
    trace, ifault, res = ffi.new("double[7]"), ffi.new("int[1]"), ffi.new("double[1]")
    lib.qfc(
        ffi.from_buffer("double[]", lam), ffi.from_buffer("double[]", nc),
        ffi.from_buffer("int[]", n), ffi.new("int[1]", [r]), ffi.new("double[1]", [0.0]),
        ffi.new("double[1]", [q]), ffi.new("int[1]", [lim]), ffi.new("double[1]", [acc]),
        trace, ifault, res,
    )  # fmt: skip
    return 1.0 - res[0], ifault[0], list(trace)


def _qfc_original_loop(lam, q, acc=1e-6, lim=1_000_000):
    nc = np.r_[np.zeros(len(lam)), 1.0]
    return _qfc(np.r_[lam, 0.0], q, acc, lim, nc=nc, n=np.ones(len(lam) + 1, np.int32))


def _exact_sum(lam, q, tr):
    """qfc's single main integration (nt + 1 terms at interval intv) in 80-bit arithmetic."""
    if tr[2] != 1:
        return None  # an auxiliary integration ran too
    lam, intv = np.asarray(lam, np.longdouble), np.longdouble(tr[3])
    u = (np.arange(int(tr[1]), dtype=np.longdouble) + 0.5) * intv
    x = 2 * lam[None, :] * u[:, None]
    s1 = -2 * u * np.longdouble(q) + np.arctan(x).sum(axis=1)
    s3 = -0.25 * np.log1p(x * x).sum(axis=1) - 0.5 * np.longdouble(tr[5]) ** 2 * u * u
    term = intv / np.longdouble(3.14159265358979) * np.exp(s3) / u * np.sin(0.5 * s1)
    return float(0.5 + term[s3 >= -50].sum())  # 1 - (0.5 - intl)


def _assert_same(lam, q, acc=1e-6, lim=1_000_000):
    p, f, tr = _qfc(lam, q, acc, lim)
    p_ref, f_ref, tr_ref = _qfc_original_loop(lam, q, acc, lim)
    assert f == f_ref
    assert tr[1:] == tr_ref[1:]  # same terms, integrations, interval, truncation point, cycles
    assert abs(p - p_ref) <= 2e-12, (p, p_ref)  # rounding only (Davies' own error is >= acc)
    exact = _exact_sum(lam, q, tr) if f == 0 or f == 2 else None
    if exact is not None:  # rounding no worse than the original loop's (whose angle sum grows)
        assert abs(p - exact) <= max(abs(p_ref - exact), 2e-15), (p, p_ref, exact)
    return p, f, tr


def test_reference_term_is_exact_on_the_original_loop():
    lam, nc = np.array([3.0, 1.0, 0.5]), np.array([0.5, 0.0, 0.2])  # non-central: original loop
    p, f, tr = _qfc(lam, 6.0, nc=nc)
    assert (p, f, tr) == _qfc(np.r_[lam, 0.0], 6.0, nc=np.r_[nc, 1.0], n=np.ones(4, np.int32))


def test_reference_matches_closed_form():
    # equal eigenvalues: P[chi2_k > q] exactly; Davies' error is absolute (acc)
    for k, q in [(1, 3.0), (5, 9.0), (40, 60.0)]:
        lam = np.ones(k)
        assert abs(_qfc_original_loop(lam, q)[0] - chi2.sf(q, k)) < 1e-6
        assert abs(_qfc(lam, q)[0] - chi2.sf(q, k)) < 1e-6


def _max_angle(lam, tr):
    return float(np.sum(np.arctan(2 * lam * tr[4])))  # sum atan(x_j) at the truncation point


@pytest.mark.parametrize(
    ("lam", "q", "min_angle"),
    [
        (np.ones(3), 2.0, np.pi),
        (np.ones(5), 8.0, 2 * np.pi),
        (np.ones(50), 70.0, 10 * np.pi),
        (np.linspace(0.01, 5.0, 200), 600.0, 10 * np.pi),
        (np.r_[400.0, np.full(60, 0.5)], 500.0, 10 * np.pi),
        (np.ones(2000), 2100.0, 100 * np.pi),
    ],
)
def test_unwrapped_angle_past_pi_2pi_10pi(lam, q, min_angle):
    # sin(0.5 * sum atan) needs the unwrapped angle; arg(prod) mod 2*pi would flip its sign
    _p, _f, tr = _assert_same(lam, q)
    assert _max_angle(lam, tr) > min_angle


@pytest.mark.parametrize("r", [1000, 5000])
@pytest.mark.parametrize("q_scale", [0.8, 1.0, 1.3])
def test_many_spread_eigenvalues(r, q_scale):
    lam = np.geomspace(1e-3, 1e3, r)
    _assert_same(lam, q_scale * lam.sum())


@pytest.mark.parametrize("r", [1000, 5000])
@pytest.mark.parametrize("acc", [1e-30, 1e-300])
def test_large_products_without_rescaling(r, acc):
    # findu puts the truncation point where sum log(1 + x^2) ~ 4 log(1/acc): |prod|^2 reaches
    # ~1e90-1e100 here (R uses acc 1e-6: ~1e28); overflow (>= 1e308) needs sum3 <= -177, a 0 term
    lam = np.ones(r)
    _p, _f, tr = _assert_same(lam, 1.05 * r, acc)
    assert np.sum(np.log1p((2 * lam * tr[4]) ** 2)) > np.log(1e90)


@pytest.mark.parametrize("r", [1000, 5000])
@pytest.mark.parametrize("q", [5.0, 20.0, 60.0])
def test_many_tiny_eigenvalues_product_near_one(r, q):
    lam = np.r_[10.0, np.full(r, 1e-9)]  # x ~ 1e-6: factors ~ 1, sum of tiny angles and logs
    _p, _f, tr = _assert_same(lam, q)
    assert 2 * 1e-9 * tr[4] < 1e-5


@pytest.mark.parametrize(
    ("lam", "q", "acc", "lim", "ifault"),
    [
        (np.array([40.0, 0.03, 0.01]), 5.0, 1e-6, 10_000, 1),  # term limit: ifault 1
        (np.array([3.0, 2.0, 1.0]), 6.0, 1e-6, 3, 4),  # cycle limit: ifault 4
        (np.array([3.0, 2.0, 1.0]), 6.0, 1e-6, 1_000_000, 0),
        (np.array([40.0, 0.03, 0.01]), 5.0, 1e-6, 1_000_000, 0),  # dominant eigenvalue, many terms
        (np.array([2.5]), 4.0, 1e-6, 1_000_000, 0),  # one eigenvalue (one variant)
        (np.array([2.0, 2.0, 1e-9]), 3.0, 1e-9, 1_000_000, 0),  # near rank-deficient
        (np.array([3.0, 2.0, 1.0]), 1e-12, 1e-6, 1_000_000, 0),  # q ~ 0: p ~ 1
        (np.array([3.0, 2.0, 1.0]), 80.0, 1e-9, 1_000_000, 0),  # tiny p (absolute-accuracy floor)
        (np.array([3.0, 2.0, 1.0]), 1e4, 1e-6, 1_000_000, 0),  # beyond the range: p = 0
    ],
)
def test_ifault_paths_and_edge_inputs(lam, q, acc, lim, ifault):
    assert _assert_same(lam, q, acc, lim)[1] == ifault


class _OriginalLoopLib:
    """qfc with the zero non-central reference term appended: every call runs the original loop."""

    def __init__(self, lib, ffi):
        self.lib, self.ffi = lib, ffi

    def qfc(self, lb1, nc1, n1, r1, sigma, c1, lim1, acc, trace, ifault, res):
        r = r1[0]
        lam = np.r_[np.frombuffer(self.ffi.buffer(lb1, 8 * r), np.float64), 0.0]
        nc = np.r_[np.frombuffer(self.ffi.buffer(nc1, 8 * r), np.float64), 1.0]
        n = np.r_[np.frombuffer(self.ffi.buffer(n1, 4 * r), np.int32), 1].astype(np.int32)
        f = self.ffi.from_buffer
        r_ref = self.ffi.new("int[1]", [r + 1])
        lb, ncs, ns = f("double[]", lam), f("double[]", nc), f("int[]", n)
        self.lib.qfc(lb, ncs, ns, r_ref, sigma, c1, lim1, acc, trace, ifault, res)


def _skato_case(kind, scale, seed=7):
    rng = np.random.default_rng(seed)
    n, p = 400, {"two_variants": 2}.get(kind, 12)  # rank < 2 never reaches SKAT-O (python_backend)
    g = rng.binomial(2, 0.02, size=(n, p)).astype(float)
    if kind == "one_group":  # all carriers among the first 40 samples
        g[40:] = 0.0
        g[:40] = rng.binomial(2, 0.2, size=(40, p))
    if kind == "rank_deficient":
        g[:, 6:] = g[:, :6]
    z1 = (g - g.mean(axis=0)) / np.sqrt(2.0)
    e = rng.standard_normal(n) + scale * g.sum(axis=1)  # scale > 0: an associated burden
    s = z1.T @ e
    return z1, float(s @ s), float(s.sum() ** 2)


_GRIDS = {
    "default": [0.0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0],
    "skat_burden": [0.0, 1.0],
    "skat_only": [0.0],
    "capped": [0.0, 0.5, 0.999],
}


# (kind, burden effect): p ~ 0.3-0.5 (null), 1e-3-1e-6, 1e-11-1e-38 (tilted tail), 0 (underflow)
_CASES = [
    ("random", 0.0), ("random", 0.3), ("random", 0.6),
    ("two_variants", 0.0), ("two_variants", 1.0), ("two_variants", 3.0),
    ("one_group", 0.0), ("one_group", 0.3), ("one_group", 3.0),
    ("rank_deficient", 0.3), ("rank_deficient", 3.0),
]  # fmt: skip


@pytest.mark.parametrize("grid", sorted(_GRIDS))
@pytest.mark.parametrize(("kind", "scale"), _CASES)
def test_skato_pvalue_matches_original_loop(monkeypatch, grid, kind, scale):
    from variantcentrifuge.association.backends import skato_integrate as si

    rho = _GRIDS[grid]
    z1, q_skat, q_burden = _skato_case(kind, scale)
    q_rho = np.array([(1 - r) * q_skat + r * q_burden for r in rho])
    p, each = si._skato_get_pvalue(q_rho, z1, rho)
    conv = si._skato_last_converged()
    monkeypatch.setattr(davies, "_qfc_lib", _OriginalLoopLib(davies._qfc_lib, davies._qfc_ffi))
    p_ref, each_ref = si._skato_get_pvalue(q_rho, z1, rho)
    assert si._skato_last_converged() == conv
    np.testing.assert_allclose(each, each_ref, rtol=1e-10, atol=0)
    assert p == pytest.approx(p_ref, rel=1e-10, abs=0)
