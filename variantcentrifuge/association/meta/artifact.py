# File: variantcentrifuge/association/meta/artifact.py
"""
Score/covariance artifact schema and gene-level computation for meta-analysis.

Per stratum we persist, for every gene, the unweighted score vector ``U = Gᵀr``
and the projection-adjusted covariance ``V = Gᵀ P G`` (the same quantity SKAT
computes internally, before weighting), plus per-variant observed ALT-allele
counts ``AC`` and allele numbers ``AN``. These are gene-level aggregates: no
per-sample rows and no phenotype values are ever stored, so artifacts are safe to
share (issue #112 privacy non-goal).

Binary traits only. For Gaussian residuals ``Var(Gᵀr) = σ²·GᵀPG`` so the
``U``/``V`` pairing here is coherent only when ``σ²=1`` (binary). Quantitative
support is deferred (spec revision 1).
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger("variantcentrifuge")

SCHEMA_VERSION = 1


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
    """

    variant_keys: list[str]
    u: np.ndarray
    v: np.ndarray
    ac: np.ndarray
    an: np.ndarray
    null_converged: bool


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


def fit_stratum_null(
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    trait_type: str,
):
    """Fit the stratum null model once (binary only), returning a NullModelResult.

    Raises
    ------
    NotImplementedError
        For non-binary traits (quantitative meta is deferred; spec revision 1).
    """
    if trait_type != "binary":
        raise NotImplementedError(
            "Meta-analysis currently supports binary traits only "
            "(quantitative deferred pending a coherent likelihood-score convention)."
        )
    from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

    backend = PythonSKATBackend()
    backend.detect_environment()
    return backend.fit_null_model(phenotype=phenotype, covariates=covariates, trait_type=trait_type)


def _null_converged(null_model) -> bool:
    """Best-effort statsmodels convergence flag (NOT bool(extra), which is always true)."""
    model = null_model.model
    converged = getattr(model, "converged", None)
    if converged is None:
        retvals = getattr(model, "mle_retvals", None)
        if isinstance(retvals, dict):
            converged = retvals.get("converged", True)
    return bool(True if converged is None else converged)


def compute_gene_score(
    genotype_matrix: np.ndarray,
    raw_dosage: np.ndarray,
    null_model,
    variant_keys: list[str],
) -> GeneScore:
    """Compute the unweighted score vector and projection-adjusted covariance.

    Parameters
    ----------
    genotype_matrix : np.ndarray, shape (n_samples, p)
        Imputed dosage matrix (no NaN) — used for the score/covariance.
    raw_dosage : np.ndarray, shape (n_samples, p)
        Pre-imputation dosage (NaN = missing) — used only for AC/AN.
    null_model : NullModelResult
        From :func:`fit_stratum_null` (binary).
    variant_keys : list[str]
        Canonical keys aligned to the columns of ``genotype_matrix``.

    Returns
    -------
    GeneScore
    """
    geno = np.asarray(genotype_matrix, dtype=np.float64)
    if geno.ndim != 2:
        raise ValueError("genotype_matrix must be 2-D (n_samples, p)")
    _n_samples, p = geno.shape
    if len(variant_keys) != p:
        raise ValueError(f"variant_keys length {len(variant_keys)} != n_variants {p}")

    residuals = np.asarray(null_model.extra["residuals"], dtype=np.float64)
    design_x = np.asarray(null_model.extra["design_matrix"], dtype=np.float64)

    # Unweighted score U = Gᵀ r
    u = geno.T @ residuals

    # Projection-adjusted covariance V = Gᵀ P G (binary form, matches _test_skato)
    mu_hat = np.asarray(null_model.extra["mu_hat"], dtype=np.float64)
    pi_1 = mu_hat * (1.0 - mu_hat)
    phi = np.sqrt(pi_1)
    z_phi = phi[:, np.newaxis] * geno
    x_phi = phi[:, np.newaxis] * design_x
    x_pi = pi_1[:, np.newaxis] * design_x
    # z_adj = phi*G - phi*X (X' diag(pi) X)^-1 (X' diag(pi) G)  [projection-adjusted]
    z_adj = z_phi - x_phi @ np.linalg.solve(x_phi.T @ x_phi, x_pi.T @ geno)
    v = z_adj.T @ z_adj
    v = 0.5 * (v + v.T)  # enforce exact symmetry against round-off

    # AC/AN from the PRE-IMPUTATION dosage (never the imputed matrix)
    raw = np.asarray(raw_dosage, dtype=np.float64)
    if raw.shape != geno.shape:
        raise ValueError("raw_dosage must have the same shape as genotype_matrix")
    observed = ~np.isnan(raw)
    an = 2.0 * observed.sum(axis=0).astype(np.float64)
    ac = np.where(observed, raw, 0.0).sum(axis=0)

    return GeneScore(
        variant_keys=list(variant_keys),
        u=u,
        v=v,
        ac=ac,
        an=an,
        null_converged=_null_converged(null_model),
    )


def build_stratum_artifact(
    stratum_id: str,
    genes: dict[str, tuple[np.ndarray, np.ndarray, list[str]]],
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    trait_type: str,
    *,
    genome_build: str,
    annotation_version: str,
    mask_id: str,
    mask_hash: str,
    n_cases: int,
    n_controls: int,
    sample_id_hashes: list[str] | None = None,
) -> ScoreCovarianceArtifact:
    """Fit the null once and compute a GeneScore for every gene.

    ``genes`` maps gene -> (imputed_matrix, raw_dosage, variant_keys).
    """
    null_model = fit_stratum_null(phenotype, covariates, trait_type)
    gene_scores: dict[str, GeneScore] = {}
    for gene, (imputed, raw, keys) in genes.items():
        gene_scores[gene] = compute_gene_score(imputed, raw, null_model, keys)
    return ScoreCovarianceArtifact(
        stratum_id=stratum_id,
        genome_build=genome_build,
        annotation_version=annotation_version,
        mask_id=mask_id,
        mask_hash=mask_hash,
        trait_type=trait_type,
        n_cases=int(n_cases),
        n_controls=int(n_controls),
        sample_id_hashes=list(sample_id_hashes or []),
        genes=gene_scores,
    )


def save_artifact(artifact: ScoreCovarianceArtifact, path: str) -> None:
    """Serialise to ``<path>`` (arrays, .npz) + ``<path>.meta.json`` (metadata)."""
    arrays: dict[str, np.ndarray] = {}
    gene_meta: dict[str, dict] = {}
    for gene, gs in artifact.genes.items():
        arrays[f"{gene}::u"] = gs.u
        arrays[f"{gene}::v"] = gs.v
        arrays[f"{gene}::ac"] = gs.ac
        arrays[f"{gene}::an"] = gs.an
        gene_meta[gene] = {
            "variant_keys": gs.variant_keys,
            "null_converged": gs.null_converged,
        }
    np.savez_compressed(path, **arrays)  # type: ignore[arg-type]
    meta = {
        "schema_version": artifact.schema_version,
        "stratum_id": artifact.stratum_id,
        "genome_build": artifact.genome_build,
        "annotation_version": artifact.annotation_version,
        "mask_id": artifact.mask_id,
        "mask_hash": artifact.mask_hash,
        "trait_type": artifact.trait_type,
        "n_cases": artifact.n_cases,
        "n_controls": artifact.n_controls,
        "sample_id_hashes": artifact.sample_id_hashes,
        "genes": gene_meta,
    }
    with open(f"{path}.meta.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)


def load_artifact(path: str) -> ScoreCovarianceArtifact:
    """Load an artifact written by :func:`save_artifact`."""
    with open(f"{path}.meta.json", encoding="utf-8") as fh:
        meta = json.load(fh)
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
    )


__all__ = [
    "SCHEMA_VERSION",
    "GeneScore",
    "ScoreCovarianceArtifact",
    "build_stratum_artifact",
    "compute_gene_score",
    "fit_stratum_null",
    "load_artifact",
    "save_artifact",
]
