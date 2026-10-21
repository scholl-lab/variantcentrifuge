"""
Imbalance-robust binary score tests: saddlepoint (SPA) and exact (ER) kernel.

Port of R ``SKAT::SKATBinary_Robust`` (SKAT 2.2.5, Zhao et al. 2020 AJHG) with
``SPAtest::Saddle_Prob`` (3.1.2), vectorised over the variants of a gene:

- per variant, the score of the covariate-projected genotype ``G1`` is tested with
  the normal approximation; where ``|z| > 2`` the p-value is replaced by a
  saddlepoint p-value (MAC >= 10) or an exact carrier-enumeration p-value (MAC < 10);
- the calibrated variance ``V~_j = S_j^2 / qchisq(1 - p~_j, 1)`` rescales the score
  covariance, the burden direction is recalibrated by an SPA burden p-value, and
  SKAT / SKAT-O are evaluated on the calibrated (score, covariance) pair.

``variant_pvalues`` is the per-variant kernel (reused for score artifacts);
``robust_skat`` is the gene-level SKAT / SKAT-O p-value.
"""

from __future__ import annotations

import numpy as np
import scipy.special
import scipy.stats

from variantcentrifuge.association.backends.base import NullModelResult
from variantcentrifuge.association.meta.combine import meta_skat, meta_skato

_CUTOFF = 2.0  # |z| above which SPA / ER replace the normal approximation
_ER_MAC = 10  # R: exact p-value when colSums(G) < 10
_TOL = np.finfo(np.float64).eps ** 0.25  # SPAtest getroot_K1 tolerance
_MAXITER = 1000
_ER_EPSILON = 1e-6  # SKATBinary tie tolerance on the exact test statistic


def _k1(t: np.ndarray, mu: np.ndarray, g: np.ndarray, q: np.ndarray) -> np.ndarray:
    """First CGF derivative minus q (SPAtest ``K1_adj``), one column of g per t."""
    mu_c = mu[:, np.newaxis]
    return np.asarray(np.sum(mu_c * g / ((1.0 - mu_c) * np.exp(-g * t) + mu_c), axis=0) - q)


def _k2(t: np.ndarray, mu: np.ndarray, g: np.ndarray) -> np.ndarray:
    """Second CGF derivative (SPAtest ``K2``, ``na.rm = TRUE``)."""
    mu_c = mu[:, np.newaxis]
    e = np.exp(-g * t)
    return np.asarray(
        np.nansum((1.0 - mu_c) * mu_c * g**2 * e / ((1.0 - mu_c) * e + mu_c) ** 2, axis=0)
    )


def _korg(t: np.ndarray, mu: np.ndarray, g: np.ndarray) -> np.ndarray:
    """Cumulant generating function (SPAtest ``Korg``)."""
    mu_c = mu[:, np.newaxis]
    return np.asarray(np.sum(np.log(1.0 - mu_c + mu_c * np.exp(g * t)), axis=0))


def _root(mu: np.ndarray, g: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Saddlepoint roots for all columns at once: SPAtest ``getroot_K1`` iterate-for-iterate."""
    gpos = np.where(g > 0, g, 0.0).sum(axis=0)
    gneg = np.where(g < 0, g, 0.0).sum(axis=0)
    outside = (q >= gpos) | (q <= gneg)
    t = np.where(outside, np.inf, 0.0)
    conv = np.ones(q.size, dtype=bool)
    jump = np.full(q.size, np.inf)
    act = np.flatnonzero(~outside)
    k1 = _k1(t[act], mu, g[:, act], q[act])
    rep = 1
    while act.size:
        t_act = t[act]
        tnew = t_act - k1 / _k2(t_act, mu, g[:, act])
        nan = np.isnan(tnew)
        small = ~nan & (np.abs(tnew - t_act) < _TOL)
        stop = nan | small | (rep == _MAXITER)
        conv[act[stop & ~small]] = False
        keep = ~stop
        act, t_act, tnew, k1 = act[keep], t_act[keep], tnew[keep], k1[keep]
        newk1 = _k1(tnew, mu, g[:, act], q[act])
        flip = np.sign(k1) != np.sign(newk1)
        step = np.abs(tnew - t_act)
        big = flip & (step > jump[act] - _TOL)
        jump[act[flip & ~big]] = step[flip & ~big]
        if big.any():
            b = act[big]
            tnew[big] = t_act[big] + np.sign(newk1[big] - k1[big]) * jump[b] / 2.0
            newk1[big] = _k1(tnew[big], mu, g[:, b], q[b])
            jump[b] /= 2.0
        t[act] = tnew
        k1 = newk1
        rep += 1
    return t, conv


def _tail(
    zeta: np.ndarray, mu: np.ndarray, g: np.ndarray, q: np.ndarray, noadj: np.ndarray
) -> np.ndarray:
    """One-sided saddlepoint tail (SPAtest ``Get_Saddle_Prob``; its error path -> noadj/2)."""
    k1 = _korg(zeta, mu, g)
    k2 = _k2(zeta, mu, g)
    w = np.sign(zeta) * np.sqrt(2.0 * (zeta * q - k1))
    z = w + np.log(zeta * np.sqrt(k2) / w) / w
    p = np.where(z > 0, scipy.stats.norm.sf(z), -scipy.stats.norm.cdf(z))
    p = np.where(np.isnan(z), noadj / 2.0, p)
    return np.where(np.isfinite(k1) & np.isfinite(k2), p, 0.0)


def saddle_pvalues(q: np.ndarray, mu: np.ndarray, g: np.ndarray, cutoff: float) -> np.ndarray:
    """Two-sided SPA p-values of the scores ``q_j = g_j'y`` (SPAtest ``Saddle_Prob``)."""
    m1 = mu @ g
    var1 = (mu * (1.0 - mu)) @ g**2
    noadj = np.asarray(scipy.stats.chi2.sf((q - m1) ** 2 / var1, 1), dtype=np.float64)
    p = noadj.copy()
    spa = np.abs(q - m1) / np.sqrt(var1) >= cutoff
    if not spa.any():
        return p
    gs, qs, na = g[:, spa], q[spa], noadj[spa]
    qinv = 2.0 * m1[spa] - qs
    with np.errstate(all="ignore"):
        r1, c1 = _root(mu, gs, qs)
        r2, c2 = _root(mu, gs, qinv)
        ps = np.abs(_tail(r1, mu, gs, qs, na)) + np.abs(_tail(r2, mu, gs, qinv, na))
        ps = np.where(c1 & c2, ps, na)
        # SPAtest re-runs with a doubled cutoff until it reaches |z|, which ends at the
        # normal p-value whenever the saddlepoint p is over 1000x smaller than it
        ps = np.where((ps != 0) & (na / ps > 1e3), na, ps)
    p[spa] = ps
    return p


def _lchoose(n: int, k: np.ndarray) -> np.ndarray:
    """log C(n, k) with the SKAT C++ convention ``k > n -> 0``."""
    with np.errstate(invalid="ignore"):
        lc = scipy.special.gammaln(n + 1) - scipy.special.gammaln(k + 1)
        lc = lc - scipy.special.gammaln(n - k + 1)
    return np.where(k > n, 0.0, np.where(k < 0, -np.inf, lc))


def exact_pvalue(geno_col: np.ndarray, y: np.ndarray, mu: np.ndarray) -> float:
    """Exact single-variant mid-p (SKATBinary ``method.bin = "ER"``) by carrier enumeration."""
    idx = np.flatnonzero(geno_col > 0)
    k, n = idx.size, y.size
    res = y - mu
    # The beta(1, 25) weight only scales the statistic against the 1e-6 tie tolerance
    z = geno_col[idx] * scipy.stats.beta.pdf(geno_col.mean() / 2.0, 1, 25)
    p1 = mu[idx]
    # P(j cases among carriers | total cases): carriers binned into 10 fitted-probability
    # groups (SKATExactBin.ComputeProb_Group, type.group = 2, k <= 50), Fisher noncentral
    pc = np.minimum(p1, 0.999)
    bins = np.searchsorted(np.arange(1, 10) / 10.0, pc, side="right")
    p2 = np.delete(mu, idx).mean()
    poly = np.ones(1)
    for b in np.unique(bins):
        size = int(np.sum(bins == b))
        pm = pc[bins == b].mean()
        c = np.arange(size + 1)
        poly = np.convolve(
            poly, scipy.special.comb(size, c) * ((pm / (1 - pm)) / (p2 / (1 - p2))) ** c
        )
    ncase = int(np.round(y).sum())
    j = np.arange(k + 1)
    with np.errstate(divide="ignore"):
        lp = np.log(poly) + _lchoose(n - k, ncase - j)
    lp[j > ncase] = -np.inf
    prob_k = np.exp(lp - lp.max())
    prob_k /= prob_k.sum()
    masks = ((np.arange(2**k)[:, np.newaxis] >> np.arange(k)) & 1).astype(bool)
    size = masks.sum(axis=1)
    stat = (masks @ z - z @ p1) ** 2
    fprob = np.prod(np.where(masks, (p1 / (1 - p1))[np.newaxis, :], 1.0), axis=1)
    denom = np.bincount(size, weights=fprob, minlength=k + 1)
    fprob = fprob / denom[size] * prob_k[size]
    fprob /= fprob.sum()
    q_obs = ((res[idx] > 0) @ z - z @ p1) ** 2
    d = q_obs - stat
    d[np.abs(d) <= _ER_EPSILON] = 0.0
    return float(fprob[d <= 0].sum() - fprob[d == 0].sum() / 2.0)


def variant_pvalues(
    geno: np.ndarray, y: np.ndarray, mu: np.ndarray, design: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-variant robust p-values and calibrated variances (R ``SPA_ER_kernel``).

    Columns with MAF > 0.5 are flipped to the minor allele first, as in R. Returns
    ``(score, p, g1, v_cal)``: the projected scores ``G1'(y - mu)``, the SPA/ER-calibrated
    p-values, the projected genotypes ``G1 = G - X (X'WX)^-1 X'W G`` (all in the flipped
    orientation; two-sided p-values are orientation-invariant) and the calibrated
    variances ``V~`` (:func:`calibrated_variance`).

    Per gene the score artifact and the robust test ask for the same columns under the
    same null (``--skat-robust`` with ``--score-artifact``): an input bit-identical to the
    previous call returns that call's (read-only) result.
    """
    geno = np.where(geno.mean(axis=0) / 2.0 > 0.5, 2.0 - geno, geno)
    key = (geno, y, mu, design)
    if _last and all(_identical(a, b) for a, b in zip(_last[0], key, strict=True)):
        hit: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] = _last[1]
        return hit
    result = _variant_pvalues(geno, y, mu, design)
    for arr in result:
        arr.flags.writeable = False  # shared by the artifact hook and the robust test
    _last[:] = [tuple(np.array(a) for a in key), result]  # copies: null arrays may change
    return result


# One-slot memo of variant_pvalues: [(geno, y, mu, design), result] (per process).
_last: list = []


def _identical(a: np.ndarray, b: np.ndarray) -> bool:
    """Same dtype, shape, layout and bytes (so the same arithmetic gives the same bits)."""
    return (
        a.dtype == b.dtype
        and a.shape == b.shape
        and a.strides == b.strides
        and a.tobytes() == b.tobytes()
    )


def _variant_pvalues(
    geno: np.ndarray, y: np.ndarray, mu: np.ndarray, design: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """:func:`variant_pvalues` on minor-allele (already flipped) columns, uncached."""
    xv = design * (mu * (1.0 - mu))[:, np.newaxis]
    g1 = geno - design @ np.linalg.solve(design.T @ xv, xv.T @ geno)
    mac = geno.sum(axis=0)
    g = g1 / np.sqrt(mac)
    q = g.T @ y
    m1 = mu @ g
    var1 = (mu * (1.0 - mu)) @ g**2
    stat = (q - m1) ** 2 / var1
    p = scipy.stats.chi2.sf(stat, 1)
    spa = (stat > _CUTOFF**2) & (mac >= _ER_MAC)
    if spa.any():
        ps = saddle_pvalues(q[spa], mu, g[:, spa], _CUTOFF)
        p[spa] = np.where(ps != 0, ps, p[spa])
    for col in np.flatnonzero((stat > _CUTOFF**2) & (mac < _ER_MAC)):
        p[col] = exact_pvalue(geno[:, col], y, mu)
    score = (q - m1) * np.sqrt(mac)
    # |z| <= 2 keeps the normal p, where V~ is exactly the score variance: return it
    # directly instead of S^2 / qchisq(chi2.sf(z^2)) (round-off at near-zero scores).
    v_cal = np.where(stat > _CUTOFF**2, calibrated_variance(score, p), var1 * mac)
    return score, p, g1, v_cal


def calibrated_variance(score: np.ndarray, p: np.ndarray) -> np.ndarray:
    """``V~ = S^2 / qchisq(1 - p~, 1)``; ``S^2 / 500`` at p~ = 0, 0 where undefined (p~ = 1)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(p > 0, score**2 / scipy.stats.chi2.isf(p, 1), score**2 / 500.0)
    return np.where(np.isfinite(v), v, 0.0)


def null_arrays(null_model: NullModelResult) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(y, mu, design)`` of a fitted binary null model (y recovered as round(r + mu))."""
    mu = np.asarray(null_model.extra["mu_hat"], dtype=np.float64)
    y = np.round(np.asarray(null_model.extra["residuals"], dtype=np.float64) + mu)
    return y, mu, np.asarray(null_model.extra["design_matrix"], dtype=np.float64)


def robust_moments(
    geno: np.ndarray, weights: np.ndarray, y: np.ndarray, mu: np.ndarray, design: np.ndarray
) -> tuple[np.ndarray, np.ndarray | None, float, np.ndarray]:
    """Calibrated unweighted score and covariance of a gene (R ``SKATBinary_spa``).

    Returns ``(score, cov, p_burden, p)``; ``diag(w) score`` and ``diag(w) cov diag(w)``
    are R's ``Score`` / ``Phi``. ``cov`` is None where R falls back to the robust
    burden p-value ``p_burden`` (degenerate burden recalibration); ``p`` are the
    per-variant SPA/ER p-values.
    """
    score, p, g1, v_cal = variant_pvalues(geno, y, mu, design)
    score = np.where(v_cal > 0, score, 0.0)
    cov = g1.T @ (g1 * (mu * (1.0 - mu))[:, np.newaxis])
    d = np.sqrt(v_cal / np.diag(cov))
    cov = cov * np.outer(d, d)
    # Every variant enters the burden; R drops p~ = 1 (exactly zero score), which made p
    # jump between S = 0 and S ~ 1e-9 (a deliberate departure from R).
    g_sum = g1 @ weights
    p_burden = float(saddle_pvalues(np.array([g_sum @ y]), mu, g_sum[:, np.newaxis], _CUTOFF)[0])
    # Burden-direction recalibration: shrink the covariance so w'Σw matches the SPA burden
    var_q = float(weights @ cov @ weights)
    var_b = float((mu * (1.0 - mu)) @ g_sum**2)
    if (g_sum @ (y - mu)) ** 2 < _CUTOFF**2 * var_b:
        # Normal branch: S^2 / qchisq(1 - p, 1) is exactly var_b; the ratio is round-off
        # when opposing variants cancel (S ~ 0, p ~ 1; R's 0/0 at S = 0), as for v_cal in
        # variant_pvalues. Departs from R where R's own round-off decides.
        var_q2 = var_b
    else:
        with np.errstate(divide="ignore"):
            var_q2 = float(score @ weights) ** 2 / scipy.stats.chi2.isf(p_burden, 1)
    r = 1.0 if var_q2 == 0 else min(var_q / var_q2, 1.0)
    if not (r > 0 and np.isfinite(r)):
        return score, None, p_burden, p
    return score, cov / r, p_burden, p


def robust_skat(
    geno: np.ndarray,
    weights: np.ndarray,
    y: np.ndarray,
    mu: np.ndarray,
    design: np.ndarray,
    method: str,
) -> tuple[float | None, float | None, np.ndarray]:
    """Robust binary SKAT ("SKAT") or SKAT-O ("SKATO") p-value (R ``SKATBinary_Robust``).

    Returns ``(p_value, rho, variant_p)``; ``rho`` is None for SKAT; ``variant_p`` are
    the per-variant SPA/ER p-values (for ACAT-V).
    """
    score, cov, p_burden, p = robust_moments(geno, weights, y, mu, design)
    if cov is None:
        return p_burden, None, p  # R: Met_SKAT_Get_Pvalue fails -> robust burden only
    if method.upper() == "SKATO" and geno.shape[1] > 1:
        p_o, rho = meta_skato([score], [cov], weights)
        if p_o is not None:
            return p_o, rho, p
    p_skat = meta_skat([score], [cov], weights)[0]
    return p_skat, (0.0 if method.upper() == "SKATO" else None), p


__all__ = [
    "calibrated_variance",
    "exact_pvalue",
    "null_arrays",
    "robust_moments",
    "robust_skat",
    "saddle_pvalues",
    "variant_pvalues",
]
