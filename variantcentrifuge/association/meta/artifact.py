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

import hashlib
import json
import logging
from dataclasses import dataclass
from typing import Any

import numpy as np

from variantcentrifuge.association.backends.base import NullModelResult
from variantcentrifuge.association.backends.null_model import (
    UntestableGeneError,
    explained_columns,
    project_out_nuisance,
)

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
) -> NullModelResult:
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


def _null_converged(null_model: NullModelResult) -> bool:
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
    null_model: NullModelResult,
    variant_keys: list[str],
    ploidy: np.ndarray | None = None,
) -> GeneScore:
    """Compute the unweighted score vector and projection-adjusted covariance.

    Parameters
    ----------
    genotype_matrix : np.ndarray, shape (n_samples, p)
        Imputed ALT-allele dosage matrix (no NaN, NOT minor-allele recoded) —
        used for the score/covariance; its sign must match the ALT keys and AC.
    raw_dosage : np.ndarray, shape (n_samples, p)
        Pre-imputation dosage (NaN = missing) — used only for AC/AN.
    null_model : NullModelResult
        From :func:`fit_stratum_null` (binary).
    variant_keys : list[str]
        Canonical keys aligned to the columns of ``genotype_matrix``.
    ploidy : np.ndarray | None, shape (n_samples, p)
        Per-genotype ploidy for AN (haploid chrX/Y calls count 1); default diploid.

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
    if len(set(variant_keys)) != p:
        raise ValueError("variant_keys must be unique (one column per variant)")

    residuals = np.asarray(null_model.extra["residuals"], dtype=np.float64)

    # Unweighted score U = Gᵀ r; variants in the covariate space contribute exactly
    # zero score and covariance (same rule as the single-cohort backend)
    explained = explained_columns(geno, null_model)
    if explained is None:
        # Never export (or zero) scores whose nuisance projection is not trustworthy
        raise UntestableGeneError("ILL_CONDITIONED_COVARIATES")
    u = np.where(explained, 0.0, geno.T @ residuals)

    # Projection-adjusted covariance V = Gᵀ P G, P = V - VX(X'VX)^-1 X'V (shared QR helper)
    z_adj = project_out_nuisance(geno, null_model)
    z_adj[:, explained] = 0.0
    v = z_adj.T @ z_adj
    v = 0.5 * (v + v.T)  # enforce exact symmetry against round-off

    # AC/AN from the PRE-IMPUTATION dosage (never the imputed matrix)
    raw = np.asarray(raw_dosage, dtype=np.float64)
    if raw.shape != geno.shape:
        raise ValueError("raw_dosage must have the same shape as genotype_matrix")
    observed = ~np.isnan(raw)
    plo = np.full(raw.shape, 2.0) if ploidy is None else np.asarray(ploidy, dtype=np.float64)
    an = np.where(observed, plo, 0.0).sum(axis=0)
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
    if meta.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"Artifact '{path}' has schema_version {meta.get('schema_version')!r}; "
            f"this version reads only {SCHEMA_VERSION}."
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


def _has_cases_and_controls(phenotype: np.ndarray) -> bool:
    return bool((phenotype == 1).any() and (phenotype == 0).any())


def hash_sample_id(sample_id: str) -> str:
    """SHA-256 of a sample ID: lets strata detect shared samples without sharing IDs."""
    return hashlib.sha256(sample_id.encode("utf-8")).hexdigest()


def emit_stratum_artifact(
    path: str,
    builders: dict[str, Any],
    *,
    stratum_id: str,
    genome_build: str,
    annotation_version: str,
    mask_id: str,
    mask_hash: str,
) -> ScoreCovarianceArtifact:
    """Write a stratum artifact from the association stage's per-gene matrix builders.

    Uses the same site/sample missingness QC as the single-cohort tests but skips the
    per-gene MAC<5 filter: a gene too rare to test in one stratum still adds
    information to the meta. Dosages are kept ALT-coded so ``U`` matches the keys.

    A gene whose retained samples (after its missingness QC) are all cases or all
    controls has no estimable null model; it is omitted with a warning. A phenotype
    that is not a usable 0/1 case/control vector for the whole cohort raises.
    """
    from variantcentrifuge.association.genotype_matrix import (
        _infer_missing_ploidy,
        _vectorized_parse_gt_columns,
        build_genotype_matrix,
    )
    from variantcentrifuge.association.meta.align import canonical_key

    first = next(iter(builders.values()), None)
    if first is None or first.phenotype_vector is None:
        raise ValueError("Score artifact needs at least one gene and a case/control phenotype")
    pheno_all = np.asarray(first.phenotype_vector, dtype=np.float64)
    if not np.isin(pheno_all, (0.0, 1.0)).all() or not _has_cases_and_controls(pheno_all):
        raise ValueError(
            "Score artifact needs a 0/1 case/control phenotype with at least one case and "
            f"one control (got {int((pheno_all == 1).sum())} cases, "
            f"{int((pheno_all == 0).sum())} controls, {pheno_all.size} samples)"
        )
    omitted: list[str] = []
    nulls: dict[bytes, NullModelResult] = {}
    gene_scores: dict[str, GeneScore] = {}
    for gene, b in builders.items():
        if b.gene_df.empty:
            continue
        missing_cols = {"CHROM", "POS", "REF", "ALT"} - set(b.gene_df.columns)
        if missing_cols:
            raise ValueError(f"Score artifact needs columns {sorted(missing_cols)}")
        geno, _mafs, sample_mask, _warn, keep = build_genotype_matrix(  # type: ignore[misc]
            b.gene_df,
            b.vcf_samples,
            b.gt_columns,
            is_binary=True,
            missing_site_threshold=b.missing_site_threshold,
            missing_sample_threshold=b.missing_sample_threshold,
            return_keep_mask=True,
        )
        # ponytail: parses GT twice per gene (build_genotype_matrix hides raw/ploidy);
        # fine for an opt-in export, expose them from the builder if it shows in profiles.
        raw, plo, kind = _vectorized_parse_gt_columns(b.gene_df, list(b.gt_columns))
        # Same retained population as build_genotype_matrix's frequencies/recoding/ploidy (#174).
        mask = np.asarray(sample_mask, dtype=bool)
        raw, plo = raw[keep], plo[keep]
        _infer_missing_ploidy(plo, kind[keep], mask)
        geno, raw, plo = geno[mask], raw.T[mask], plo.T[mask]
        observed = ~np.isnan(raw)
        ac = np.where(observed, raw, 0.0).sum(axis=0)
        an = np.where(observed, plo, 0.0).sum(axis=0)
        flipped = (ac > an / 2.0) & (ac < an)  # build_genotype_matrix minor-allele recoding
        geno[:, flipped] = plo[:, flipped] - geno[:, flipped]

        rows = b.gene_df.loc[keep]
        keys = [
            canonical_key(c, p, r, a)
            for c, p, r, a in zip(rows["CHROM"], rows["POS"], rows["REF"], rows["ALT"], strict=True)
        ]
        if not keys:
            continue
        if not _has_cases_and_controls(pheno_all[mask]):
            logger.warning(
                "Score artifact: omitting gene %s: its %d retained samples (after "
                "missingness QC) are not both cases and controls, so no null model fits",
                gene,
                int(mask.sum()),
            )
            omitted.append(gene)
            continue
        null_key = mask.tobytes()
        if null_key not in nulls:
            cov = b.covariate_matrix[mask] if b.covariate_matrix is not None else None
            nulls[null_key] = fit_stratum_null(pheno_all[mask], cov, "binary")
        try:
            gene_scores[gene] = compute_gene_score(geno, raw, nulls[null_key], keys, ploidy=plo)
        except UntestableGeneError as exc:
            logger.warning("Score artifact: omitting gene %s: %s", gene, exc.reason)
            omitted.append(gene)

    if not gene_scores:
        raise ValueError(
            f"Score artifact: no genes survived QC ({len(builders)} genes: empty, removed by "
            f"site missingness QC, or {len(omitted)} untestable)"
        )
    artifact = ScoreCovarianceArtifact(
        stratum_id=stratum_id,
        genome_build=genome_build,
        annotation_version=annotation_version,
        mask_id=mask_id,
        mask_hash=mask_hash,
        trait_type="binary",
        n_cases=int(pheno_all.sum()),
        n_controls=int((1.0 - pheno_all).sum()),
        sample_id_hashes=[hash_sample_id(s) for s in first.vcf_samples],
        genes=gene_scores,
    )
    save_artifact(artifact, path)
    logger.info("Wrote meta-analysis score artifact (%d genes) to %s", len(gene_scores), path)
    return artifact


__all__ = [
    "SCHEMA_VERSION",
    "GeneScore",
    "ScoreCovarianceArtifact",
    "build_stratum_artifact",
    "compute_gene_score",
    "emit_stratum_artifact",
    "fit_stratum_null",
    "hash_sample_id",
    "load_artifact",
    "save_artifact",
]
