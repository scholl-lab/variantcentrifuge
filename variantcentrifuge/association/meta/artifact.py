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
import logging
from typing import Any

import numpy as np

from variantcentrifuge.association.backends.base import NullModelResult
from variantcentrifuge.association.backends.null_model import (
    UntestableGeneError,
    _project,
    classify_projected,
)
from variantcentrifuge.association.meta.artifact_schema import (
    SCHEMA_VERSION,
    V1_MISSING_POLICY,
    GeneScore,
    ScoreCovarianceArtifact,
    load_artifact,
    save_artifact,
)

logger = logging.getLogger("variantcentrifuge")


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


def _allele_counts(raw: np.ndarray, ploidy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-variant AC and AN over the observed (non-NaN) calls of samples x variants arrays.

    Without missing calls the masked ``where`` copy equals the input, so the sums are too.
    """
    observed = ~np.isnan(raw)
    if observed.all() and raw.dtype.kind == ploidy.dtype.kind == "f":  # where() made ints float
        return raw.sum(axis=0), ploidy.sum(axis=0)
    return np.where(observed, raw, 0.0).sum(axis=0), np.where(observed, ploidy, 0.0).sum(axis=0)


def compute_gene_score(
    genotype_matrix: np.ndarray,
    raw_dosage: np.ndarray,
    null_model: NullModelResult,
    variant_keys: list[str],
    ploidy: np.ndarray | None = None,
    robust: bool = False,
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
    robust : bool
        Also store the SPA/exact-calibrated variances ``v_robust`` (``--skat-robust``).

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
    z_adj, z, bound = _project(geno, null_model)  # one projection for both uses
    explained = classify_projected(z_adj, z, bound)
    if explained is None:
        # Never export (or zero) scores whose nuisance projection is not trustworthy
        raise UntestableGeneError("ILL_CONDITIONED_COVARIATES")
    u = np.where(explained, 0.0, geno.T @ residuals)

    # Projection-adjusted covariance V = Gᵀ P G, P = V - VX(X'VX)^-1 X'V (shared QR helper)
    z_adj[:, explained] = 0.0
    v = z_adj.T @ z_adj
    v = 0.5 * (v + v.T)  # enforce exact symmetry against round-off

    v_robust = None
    if robust:
        from variantcentrifuge.association.spa import null_arrays, variant_pvalues

        v_robust = np.zeros(p)
        if not explained.all():
            v_robust[~explained] = variant_pvalues(geno[:, ~explained], *null_arrays(null_model))[3]

    # AC/AN from the PRE-IMPUTATION dosage (never the imputed matrix)
    raw = np.asarray(raw_dosage, dtype=np.float64)
    if raw.shape != geno.shape:
        raise ValueError("raw_dosage must have the same shape as genotype_matrix")
    plo = np.full(raw.shape, 2.0) if ploidy is None else np.asarray(ploidy, dtype=np.float64)
    ac, an = _allele_counts(raw, plo)

    return GeneScore(
        variant_keys=list(variant_keys),
        u=u,
        v=v,
        ac=ac,
        an=an,
        null_converged=_null_converged(null_model),
        n_carriers=int((raw > 0).any(axis=1).sum()),  # NaN (missing) is not a carrier
        v_robust=v_robust,
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
    missing_policy: dict[str, Any],
) -> ScoreCovarianceArtifact:
    """Fit the null once and compute a GeneScore for every gene.

    ``genes`` maps gene -> (imputed_matrix, raw_dosage, variant_keys).
    ``missing_policy`` states how the imputed matrices were filled (see
    :class:`ScoreCovarianceArtifact`); it is recorded, not inferred.
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
        missing_policy=dict(missing_policy),
    )


def _has_cases_and_controls(phenotype: np.ndarray) -> bool:
    return bool((phenotype == 1).any() and (phenotype == 0).any())


def hash_sample_id(sample_id: str) -> str:
    """SHA-256 of a sample ID: lets strata detect shared samples without sharing IDs."""
    return hashlib.sha256(sample_id.encode("utf-8")).hexdigest()


# Distinct per-gene sample masks whose null model is kept; a refit is deterministic,
# so evicting only costs time and memory no longer grows with the number of masks.
_MAX_CACHED_NULLS = 8


def artifact_phenotype(builders: dict[str, Any]) -> np.ndarray:
    """The cohort's 0/1 phenotype from the first builder; raises if unusable."""
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
    return pheno_all


def null_for_mask(
    nulls: dict[bytes, NullModelResult | UntestableGeneError],
    mask: np.ndarray,
    pheno_all: np.ndarray,
    covariates: np.ndarray | None,
) -> NullModelResult:
    """Fitted null model for this sample mask (bounded cache; untestable sets cached too)."""
    key = mask.tobytes()
    if key not in nulls:
        if len(nulls) >= _MAX_CACHED_NULLS:
            del nulls[next(iter(nulls))]
        cov = covariates[mask] if covariates is not None else None
        try:
            nulls[key] = fit_stratum_null(pheno_all[mask], cov, "binary")
        except UntestableGeneError as exc:
            nulls[key] = exc
    null_model = nulls[key]
    if isinstance(null_model, UntestableGeneError):
        raise null_model
    return null_model


def gene_score_record(
    gene: str,
    b: Any,
    built: dict[str, Any],
    pheno_all: np.ndarray,
    nulls: dict[bytes, NullModelResult | UntestableGeneError],
    robust: bool = False,
    sites: list[np.ndarray] | None = None,
) -> tuple[int, int, GeneScore | str | None]:
    """One gene's artifact entry from its builder and that builder's result.

    ``built`` holds ``qc_genotype_matrix`` (retained samples, before the MAC check),
    ``sample_mask``, ``variant_keep_mask`` and ``parsed`` (the dosage/ploidy/kind arrays
    the matrix was built from). ``sites`` are the run's CHROM, POS, REF and ALT arrays
    that ``b.rows`` index (``b.gene_df`` is used without them). Returns (sites, kept
    sites, entry): a GeneScore, the omission reason, or None for a gene with no kept
    variant.
    """
    from variantcentrifuge.association.genotype_matrix import _infer_missing_ploidy
    from variantcentrifuge.association.meta.align import canonical_key

    missing_cols = {"CHROM", "POS", "REF", "ALT"} - set(b.gene_df.columns)
    if missing_cols:
        raise ValueError(f"Score artifact needs columns {sorted(missing_cols)}")
    keep, mask = built["variant_keep_mask"], built["sample_mask"]
    raw, plo, kind = built["parsed"]
    # Same retained population as build_genotype_matrix's frequencies/recoding/ploidy (#174).
    raw, plo = raw[keep], plo[keep]
    _infer_missing_ploidy(plo, kind[keep], mask)
    if b.missing_as_ref:
        raw[np.isnan(raw)] = 0.0  # same REF fill as build_genotype_matrix: full AN
    # C-order copy, as geno[mask] gave: the tests keep their matrix, and BLAS sums alike
    geno = np.array(built["qc_genotype_matrix"], dtype=np.float64, order="C")
    raw, plo = raw.T[mask], plo.T[mask]  # C-order: fractional dosages keep this sum order
    ac, an = _allele_counts(raw, plo)
    flipped = (ac > an / 2.0) & (ac < an)  # build_genotype_matrix minor-allele recoding
    geno[:, flipped] = plo[:, flipped] - geno[:, flipped]

    if sites is not None and b.rows is not None:
        site_cols = [col[b.rows[keep]] for col in sites]  # no pandas per gene
    else:
        rows = b.gene_df.loc[keep]
        site_cols = [rows[c] for c in ("CHROM", "POS", "REF", "ALT")]
    keys = [canonical_key(c, p, r, a) for c, p, r, a in zip(*site_cols, strict=True)]
    n_sites, n_kept = len(keep), int(keep.sum())
    if not keys:
        return n_sites, n_kept, None
    if not _has_cases_and_controls(pheno_all[mask]):
        logger.warning(
            "Score artifact: omitting gene %s: its %d retained samples (after "
            "missingness QC) are not both cases and controls, so no null model fits",
            gene,
            int(mask.sum()),
        )
        return n_sites, n_kept, "ONE_CLASS"
    try:
        null_model = null_for_mask(nulls, mask, pheno_all, b.covariate_matrix)
        score = compute_gene_score(geno, raw, null_model, keys, ploidy=plo, robust=robust)
        return n_sites, n_kept, score
    except UntestableGeneError as exc:
        logger.warning("Score artifact: omitting gene %s: %s", gene, exc.reason)
        return n_sites, n_kept, exc.reason


def _build_for_artifact(b: Any) -> dict[str, Any]:
    """Build one gene's QC'd matrix the way the artifact needs it (no MAC check)."""
    from variantcentrifuge.association.genotype_matrix import build_genotype_matrix

    parsed = b.gene_arrays()
    geno, _mafs, sample_mask, _warn, keep = build_genotype_matrix(  # type: ignore[misc]
        b.gene_df,
        b.vcf_samples,
        b.gt_columns,
        is_binary=True,
        missing_site_threshold=b.missing_site_threshold,
        missing_sample_threshold=b.missing_sample_threshold,
        return_keep_mask=True,
        missing_as_ref=b.missing_as_ref,
        parsed=parsed,
    )
    mask = np.asarray(sample_mask, dtype=bool)
    return {
        "qc_genotype_matrix": geno[mask],
        "sample_mask": mask,
        "variant_keep_mask": keep,
        "parsed": parsed,
    }


def emit_stratum_artifact(
    path: str,
    builders: dict[str, Any],
    *,
    records: dict[str, tuple[int, int, GeneScore | str | None]] | None = None,
    stratum_id: str,
    genome_build: str,
    annotation_version: str,
    mask_id: str,
    mask_hash: str,
    mask_label: str | None = None,
    provenance: dict[str, Any] | None = None,
    robust: bool = False,
) -> ScoreCovarianceArtifact:
    """Write a stratum artifact from the association stage's per-gene matrix builders.

    Uses the same site/sample missingness QC as the single-cohort tests but skips the
    per-gene MAC<5 filter: a gene too rare to test in one stratum still adds
    information to the meta. Dosages are kept ALT-coded so ``U`` matches the keys.
    The builders' missing-data policy is recorded in the artifact; with
    ``missing_as_ref`` missing calls count as REF in AC and AN (#200).

    ``records`` (gene -> ``gene_score_record`` result, computed by the engine's workers
    from the matrices they built for the tests, #212) skips building here for those
    genes; genes are written in ``builders`` order either way. ``robust``
    (``--skat-robust``) also stores each gene's calibrated variances ``v_robust``.

    A gene whose retained samples (after its missingness QC) are all cases or all
    controls has no estimable null model; it is omitted with a warning. A phenotype
    that is not a usable 0/1 case/control vector for the whole cohort raises.
    """
    pheno_all = artifact_phenotype(builders)
    first = next(iter(builders.values()))
    nulls: dict[bytes, NullModelResult | UntestableGeneError] = {}
    omitted: list[str] = []
    gene_scores: dict[str, GeneScore] = {}
    n_sites = n_kept = 0
    for gene, b in builders.items():
        if b.gene_df.empty:
            continue
        record = None if records is None else records.get(gene)
        if record is None:  # no engine result for this gene: build its matrix here
            built = _build_for_artifact(b)
            record = gene_score_record(gene, b, built, pheno_all, nulls, robust)
        sites, kept, entry = record
        n_sites, n_kept = n_sites + sites, n_kept + kept
        if isinstance(entry, GeneScore):
            gene_scores[gene] = entry
        elif entry is not None:
            omitted.append(gene)

    if n_kept < n_sites:
        logger.warning(
            "Score artifact: site missingness QC removed %d of %d variants (missing fraction "
            "> %s). For merged VCFs where ./. means 'no variant called', consider "
            "--missing-as-ref.",
            n_sites - n_kept,
            n_sites,
            first.missing_site_threshold,
        )
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
        missing_policy={
            "imputation": "missing_as_ref" if first.missing_as_ref else "bestguess",
            "site_threshold": first.missing_site_threshold,
            "sample_threshold": first.missing_sample_threshold,
        },
        mask_label=mask_label,
        provenance=dict(provenance or {}),
    )
    save_artifact(artifact, path)
    logger.info("Wrote meta-analysis score artifact (%d genes) to %s", len(gene_scores), path)
    return artifact


__all__ = [
    "SCHEMA_VERSION",
    "V1_MISSING_POLICY",
    "GeneScore",
    "ScoreCovarianceArtifact",
    "artifact_phenotype",
    "build_stratum_artifact",
    "compute_gene_score",
    "emit_stratum_artifact",
    "fit_stratum_null",
    "gene_score_record",
    "hash_sample_id",
    "load_artifact",
    "null_for_mask",
    "save_artifact",
]
