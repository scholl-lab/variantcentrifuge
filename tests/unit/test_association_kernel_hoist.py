"""Hoisted per-gene kernel quantities equal the former per-gene computation exactly (#P4)."""

from __future__ import annotations

import numpy as np
import pytest

from variantcentrifuge.association.backends import null_model as nm
from variantcentrifuge.association.backends.null_model import fit_glm_null_model


def _project_reference(z, null_model):
    """The pre-0.23 per-gene ``_project``: QR and conditioning recomputed on every call."""
    z = z if z.ndim == 2 else z[:, np.newaxis]
    design_x = np.asarray(null_model.extra["design_matrix"], dtype=np.float64)
    if null_model.trait_type == "binary":
        mu_hat = np.asarray(null_model.extra["mu_hat"], dtype=np.float64)
        phi = np.sqrt(mu_hat * (1.0 - mu_hat))[:, np.newaxis]
        z = phi * z
        design_x = phi * design_x
    basis, r_factor = np.linalg.qr(design_x)
    adjusted = z - basis @ (basis.T @ z)
    r_cov = np.asarray(null_model.extra.get("covariate_r", np.empty((0, 0))), dtype=np.float64)
    original = r_factor.copy()
    if r_cov.size:
        original[:, 1:] = r_factor[:, 1:] @ r_cov
    kappa = float(np.linalg.cond(original / np.linalg.norm(original, axis=0)))
    n, k = design_x.shape
    lam = np.sqrt(2.0 * np.log(2.0 * n / 1e-6))
    bound = float(np.finfo(np.float64).eps * lam * np.sqrt(n) * (k + np.sqrt(k) * kappa))
    return adjusted, z, bound


def _null(trait, seed, n=400, k=5):
    rng = np.random.default_rng(seed)
    cov = rng.normal(size=(n, k))
    if trait == "binary":
        y = (rng.random(n) < 1 / (1 + np.exp(-(cov[:, 0] * 0.5 - 1.0)))).astype(float)
    else:
        y = cov[:, 0] + rng.normal(size=n)
    return fit_glm_null_model(y, cov, trait), rng


@pytest.mark.parametrize("trait", ["binary", "quantitative"])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_project_matches_per_gene_reference(trait, seed):
    null, rng = _null(trait, seed)
    for p in (1, 3, 17):
        z = rng.binomial(2, 0.05, size=(null.n_samples, p)).astype(float) * rng.random(p)
        for _ in range(2):  # first call fills the cache, second reuses it
            got, ref = nm._project(z, null), _project_reference(z, null)
            assert np.array_equal(got[0], ref[0]) and np.array_equal(got[1], ref[1])
            assert got[2] == ref[2]
    vec = rng.random(null.n_samples)
    assert np.array_equal(nm._project(vec, null)[0], _project_reference(vec, null)[0])


def test_null_model_qr_computed_once(monkeypatch):
    null, rng = _null("binary", 3)
    calls = []
    real_qr = np.linalg.qr
    monkeypatch.setattr(np.linalg, "qr", lambda *a, **kw: calls.append(1) or real_qr(*a, **kw))
    z = rng.random((null.n_samples, 4))
    for _ in range(5):
        nm.project_out_nuisance(z, null)
        nm.explained_columns(z, null)
    assert len(calls) == 1


def test_cache_follows_reassigned_null_quantities():
    null, rng = _null("quantitative", 4)
    z = rng.random((null.n_samples, 3))
    before = nm._project(z, null)[2]
    near = np.eye(5)
    near[0, 1], near[1, 1] = 1.0, 1e-6  # near-parallel covariates: larger kappa and bound
    null.extra["covariate_r"] = near
    after = nm._project(z, null)
    assert after[2] == _project_reference(z, null)[2] != before


def test_liu_fallback_only_computed_when_davies_fails(monkeypatch):
    from variantcentrifuge.association.backends import davies

    if not davies._try_load_davies():
        pytest.skip("Davies C extension not built")
    lam = np.array([3.0, 1.5, 0.7, 0.2])
    expected = davies.compute_pvalue(4.0, lam)
    calls = []
    real = davies._liu_pvalue
    monkeypatch.setattr(davies, "_liu_pvalue", lambda *a: calls.append(1) or real(*a))
    assert davies.compute_pvalue(4.0, lam) == expected and expected[1] == "davies"
    assert calls == []
    monkeypatch.setattr(davies, "_try_load_davies", lambda: False)  # fallback path unchanged
    p, method, _ = davies.compute_pvalue(4.0, lam)
    assert method in ("saddlepoint", "liu") and 0.0 < p < 1.0
