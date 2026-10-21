# File: variantcentrifuge/association/meta/artifact_schema.py
"""Score/covariance artifact schema and its (de)serialisation (``.npz`` + ``.meta.json``)."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any

import numpy as np

SCHEMA_VERSION = 3  # 3 adds the optional per-gene v_robust array (#223)
# Artifacts without v_robust are written as schema 2: identical to, and readable by, v0.23.0.
# v1 artifacts predate the recorded policy: they always used bestguess imputation (#200).
V1_MISSING_POLICY = {"imputation": "bestguess", "site_threshold": None, "sample_threshold": None}


@dataclass
class GeneScore:
    """Aggregate score/covariance for one gene in one stratum.

    Fields
    ------
    variant_keys : list[str]
        Canonical ``chrom:pos:ref:alt`` keys, one per column of ``u``/``v``.
    u : np.ndarray, shape (p,)
        Unweighted score vector ``Gᵀ r`` (r = null residuals).
    v : np.ndarray, shape (p, p)
        Unweighted projection-adjusted covariance ``Gᵀ P G`` (== Var(u) for binary).
    ac : np.ndarray, shape (p,)
        Observed ALT-allele counts (pre-imputation).
    an : np.ndarray, shape (p,)
        Observed allele numbers = 2 * non-missing samples (pre-imputation).
    null_converged : bool
        Whether the stratum null model converged.
    n_carriers : int | None
        Samples with >= 1 observed ALT allele over the gene's variants (#194);
        None for artifacts written before it was recorded.
    v_robust : np.ndarray | None, shape (p,)
        SPA/exact-calibrated score variances ``V~ = U^2 / qchisq(1 - p~, 1)``
        (``--skat-robust`` runs, schema 3; see ``association/spa.py``); else None.
    """

    variant_keys: list[str]
    u: np.ndarray
    v: np.ndarray
    ac: np.ndarray
    an: np.ndarray
    null_converged: bool
    n_carriers: int | None = None
    v_robust: np.ndarray | None = None


@dataclass
class ScoreCovarianceArtifact:
    """A stratum's meta-analysis artifact: metadata + per-gene scores."""

    stratum_id: str
    genome_build: str
    annotation_version: str
    mask_id: str
    mask_hash: str
    trait_type: str
    n_cases: int
    n_controls: int
    sample_id_hashes: list[str]
    genes: dict[str, GeneScore]
    schema_version: int = SCHEMA_VERSION
    # imputation ("bestguess" | "missing_as_ref") + site/sample missingness thresholds
    # Required, keyword-only: how missing genotypes were handled must be stated, never
    # defaulted (V1_MISSING_POLICY applies only when loading a v1 artifact).
    missing_policy: dict[str, Any] = field(kw_only=True)
    # Optional v2 fields (#197): readable mask name; annotation provenance (meta/provenance.py)
    mask_label: str | None = None
    provenance: dict[str, Any] = field(default_factory=dict)


def save_artifact(artifact: ScoreCovarianceArtifact, path: str) -> None:
    """Serialise to ``<path>`` (arrays, .npz) + ``<path>.meta.json`` (metadata).

    Written as schema 3 only when some gene has ``v_robust``; otherwise as schema 2.
    """
    arrays: dict[str, np.ndarray] = {}
    gene_meta: dict[str, dict] = {}
    for gene, gs in artifact.genes.items():
        arrays[f"{gene}::u"] = gs.u
        arrays[f"{gene}::v"] = gs.v
        arrays[f"{gene}::ac"] = gs.ac
        arrays[f"{gene}::an"] = gs.an
        if gs.v_robust is not None:
            arrays[f"{gene}::v_robust"] = gs.v_robust
        gene_meta[gene] = {
            "variant_keys": gs.variant_keys,
            "null_converged": gs.null_converged,
            "n_carriers": gs.n_carriers,
        }
    np.savez_compressed(path, **arrays)  # type: ignore[arg-type]
    robust = any(gs.v_robust is not None for gs in artifact.genes.values())
    meta = {
        "schema_version": SCHEMA_VERSION if robust else min(artifact.schema_version, 2),
        "stratum_id": artifact.stratum_id,
        "genome_build": artifact.genome_build,
        "annotation_version": artifact.annotation_version,
        "mask_id": artifact.mask_id,
        "mask_hash": artifact.mask_hash,
        "trait_type": artifact.trait_type,
        "n_cases": artifact.n_cases,
        "n_controls": artifact.n_controls,
        "sample_id_hashes": artifact.sample_id_hashes,
        "missing_policy": artifact.missing_policy,
        "mask_label": artifact.mask_label,
        "provenance": artifact.provenance,
        "genes": gene_meta,
    }
    with open(f"{path}.meta.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)


def load_artifact(path: str) -> ScoreCovarianceArtifact:
    """Load an artifact written by :func:`save_artifact`."""
    with open(f"{path}.meta.json", encoding="utf-8") as fh:
        meta = json.load(fh)
    if meta.get("schema_version") not in (1, 2, SCHEMA_VERSION):
        raise ValueError(
            f"Artifact '{path}' has schema_version {meta.get('schema_version')!r}; "
            f"this version reads 1 to {SCHEMA_VERSION}."
        )
    policy = meta.get("missing_policy") if meta["schema_version"] > 1 else V1_MISSING_POLICY
    if not isinstance(policy, dict) or policy.get("imputation") not in (
        "bestguess",
        "missing_as_ref",
    ):
        raise ValueError(
            f"Artifact '{path}' (schema_version {meta['schema_version']}) has no valid "
            "missing_policy "
            f"(got {policy!r}); refusing to guess how its missing genotypes were handled."
        )
    npz_path = path if path.endswith(".npz") else f"{path}.npz"
    with np.load(npz_path) as data:
        genes: dict[str, GeneScore] = {}
        for gene, gmeta in meta["genes"].items():
            genes[gene] = GeneScore(
                variant_keys=list(gmeta["variant_keys"]),
                u=data[f"{gene}::u"],
                v=data[f"{gene}::v"],
                ac=data[f"{gene}::ac"],
                an=data[f"{gene}::an"],
                null_converged=bool(gmeta["null_converged"]),
                n_carriers=gmeta.get("n_carriers"),
                v_robust=data[f"{gene}::v_robust"] if f"{gene}::v_robust" in data.files else None,
            )
    return ScoreCovarianceArtifact(
        stratum_id=meta["stratum_id"],
        genome_build=meta["genome_build"],
        annotation_version=meta["annotation_version"],
        mask_id=meta["mask_id"],
        mask_hash=meta["mask_hash"],
        trait_type=meta["trait_type"],
        n_cases=int(meta["n_cases"]),
        n_controls=int(meta["n_controls"]),
        sample_id_hashes=list(meta["sample_id_hashes"]),
        genes=genes,
        schema_version=int(meta["schema_version"]),
        missing_policy=dict(policy),
        mask_label=meta.get("mask_label"),
        provenance=dict(meta.get("provenance") or {}),
    )
