"""Artifact computation, round-trip, and safe-aggregate-output guarantee (#6)."""

from __future__ import annotations

import json

import numpy as np

from tests.unit.meta_helpers import make_stratum
from variantcentrifuge.association.backends.python_backend import PythonSKATBackend
from variantcentrifuge.association.meta.artifact import (
    compute_gene_score,
    load_artifact,
    save_artifact,
)

KEYS = {"GENEA": [f"1:{100 + j}:A:G" for j in range(8)]}


def test_gene_score_shapes_symmetry_and_psd():
    rng = np.random.default_rng(0)
    n, p = 400, 8
    geno = rng.binomial(2, 0.03, size=(n, p)).astype(float)
    y = (rng.uniform(size=n) < 0.2).astype(float)
    backend = PythonSKATBackend()
    backend.detect_environment()
    nm = backend.fit_null_model(phenotype=y, covariates=None, trait_type="binary")

    gs = compute_gene_score(geno, geno, nm, list(KEYS["GENEA"]))
    assert gs.u.shape == (p,)
    assert gs.v.shape == (p, p)
    np.testing.assert_allclose(gs.v, gs.v.T, atol=1e-10)
    eigvals = np.linalg.eigvalsh(gs.v)
    assert eigvals.min() > -1e-8  # positive semidefinite


def test_ac_an_use_pre_imputation_dosage():
    # Raw dosage has a missing call (NaN) that imputation would have filled.
    raw = np.array([[0.0], [1.0], [np.nan], [2.0]])
    imputed = np.array([[0.0], [1.0], [0.5], [2.0]])
    backend = PythonSKATBackend()
    backend.detect_environment()
    y = np.array([0.0, 1.0, 0.0, 1.0])
    nm = backend.fit_null_model(phenotype=y, covariates=None, trait_type="binary")
    gs = compute_gene_score(imputed, raw, nm, ["1:1:A:G"])
    # AN counts 3 observed samples (2 alleles each); AC excludes the imputed 0.5.
    assert gs.an[0] == 6.0
    assert gs.ac[0] == 3.0


def test_save_load_round_trip(tmp_path):
    art = make_stratum("A", seed=7, gene_keys=KEYS)
    path = str(tmp_path / "stratumA")
    save_artifact(art, path)
    loaded = load_artifact(path)

    assert loaded.stratum_id == art.stratum_id
    assert loaded.mask_hash == art.mask_hash
    gs0, gs1 = art.genes["GENEA"], loaded.genes["GENEA"]
    np.testing.assert_allclose(gs0.u, gs1.u)
    np.testing.assert_allclose(gs0.v, gs1.v)
    assert gs1.variant_keys == gs0.variant_keys


def test_artifact_sidecar_contains_only_safe_aggregates(tmp_path):
    """Acceptance #6: no phenotype rows / raw sample IDs in the artifact metadata."""
    art = make_stratum("A", seed=9, gene_keys=KEYS)
    path = str(tmp_path / "stratumA")
    save_artifact(art, path)
    meta = json.loads((tmp_path / "stratumA.meta.json").read_text())

    # Only hashed sample identifiers, counts, and aggregate keys — no phenotype values.
    assert set(meta) >= {"n_cases", "n_controls", "sample_id_hashes", "genes"}
    assert "phenotype" not in json.dumps(meta).lower()
    # Sample hashes are opaque tokens, never real IDs with PII structure.
    assert all(h.startswith("A_s") for h in meta["sample_id_hashes"])
