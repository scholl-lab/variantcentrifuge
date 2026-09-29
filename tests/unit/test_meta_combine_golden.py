"""Golden-gate tests: single-stratum meta must reproduce the trusted SKAT backend.

If ``meta_skat``/``meta_skato`` on one stratum equal ``PythonSKATBackend._test_skat``/
``_test_skato`` for the same data and weights, the meta score/covariance algebra is
correct by construction. This ties all new statistics to already-validated code.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.stats

from variantcentrifuge.association.backends.python_backend import PythonSKATBackend
from variantcentrifuge.association.meta.artifact import compute_gene_score
from variantcentrifuge.association.meta.combine import meta_burden, meta_skat, meta_skato


def _synthetic_binary(seed: int, n: int = 500, p: int = 18, with_signal: bool = False):
    """Return (G_imputed, y binary, X covariates, variant_keys) — no missingness."""
    rng = np.random.default_rng(seed)
    mafs = rng.uniform(0.01, 0.06, size=p)
    geno = rng.binomial(2, mafs[np.newaxis, :], size=(n, p)).astype(np.float64)
    # Two covariates (e.g. sex, one PC)
    cov = np.column_stack([rng.normal(size=n), rng.integers(0, 2, size=n).astype(float)])
    logit = -1.5 + 0.3 * cov[:, 0]
    if with_signal:
        logit = logit + 0.8 * geno.sum(axis=1) / max(1.0, geno.sum(axis=1).std())
    prob = 1.0 / (1.0 + np.exp(-logit))
    y = (rng.uniform(size=n) < prob).astype(np.float64)
    keys = [f"1:{1000 + j}:A:G" for j in range(p)]
    return geno, y, cov, keys


def _fit(geno, y, cov):
    backend = PythonSKATBackend()
    backend.detect_environment()
    null_model = backend.fit_null_model(phenotype=y, covariates=cov, trait_type="binary")
    return backend, null_model


@pytest.mark.parametrize("seed", [1, 7, 42])
def test_meta_skat_matches_backend_single_stratum(seed):
    geno, y, cov, keys = _synthetic_binary(seed)
    backend, null_model = _fit(geno, y, cov)
    weights = scipy.stats.beta.pdf(geno.mean(axis=0) / 2.0, 1.0, 25.0)

    ref = backend.test_gene("GENE", geno, null_model, method="SKAT", weights=weights)
    gs = compute_gene_score(geno, geno, null_model, keys)
    meta_p, _method = meta_skat([gs.u], [gs.v], weights)

    assert ref["p_value"] is not None and meta_p is not None
    assert meta_p == pytest.approx(ref["p_value"], rel=1e-6, abs=1e-12)


@pytest.mark.parametrize("seed", [2, 11, 43])
def test_meta_skato_matches_backend_single_stratum(seed):
    geno, y, cov, keys = _synthetic_binary(seed, with_signal=True)
    backend, null_model = _fit(geno, y, cov)
    weights = scipy.stats.beta.pdf(geno.mean(axis=0) / 2.0, 1.0, 25.0)

    ref = backend.test_gene("GENE", geno, null_model, method="SKATO", weights=weights)
    gs = compute_gene_score(geno, geno, null_model, keys)
    meta_p, _rho = meta_skato([gs.u], [gs.v], weights)

    assert ref["p_value"] is not None and meta_p is not None
    # SKAT-O involves numerical integration; allow a slightly looser tolerance.
    assert meta_p == pytest.approx(ref["p_value"], rel=1e-4, abs=1e-9)


def test_meta_burden_matches_collapsed_skat():
    """Burden (normal approx) must agree with SKAT on the single collapsed super-variant."""
    geno, y, cov, keys = _synthetic_binary(seed=5, with_signal=True)
    _backend, null_model = _fit(geno, y, cov)
    weights = scipy.stats.beta.pdf(geno.mean(axis=0) / 2.0, 1.0, 25.0)

    gs = compute_gene_score(geno, geno, null_model, keys)
    burden_p, beta, direction = meta_burden([gs.u], [gs.v], weights)

    # Collapse to a single weighted super-variant and run SKAT on it.
    collapsed = (geno @ weights).reshape(-1, 1)
    gs_c = compute_gene_score(collapsed, collapsed, null_model, ["1:1:A:G"])
    skat_p, _method = meta_skat([gs_c.u], [gs_c.v], np.array([1.0]))

    assert burden_p is not None and skat_p is not None
    assert burden_p == pytest.approx(skat_p, rel=5e-3, abs=1e-6)
    assert direction in (-1, 1)
    assert np.isfinite(beta)
