"""Synthetic-stratum builders for meta-analysis tests (not collected by pytest)."""

from __future__ import annotations

import numpy as np

from variantcentrifuge.association.meta.artifact import (
    ScoreCovarianceArtifact,
    build_stratum_artifact,
)


def make_stratum(
    stratum_id: str,
    seed: int,
    *,
    gene_keys: dict[str, list[str]],
    effect_gene: str | None = None,
    effect_sign: float = 0.0,
    effect_size: float = 1.2,
    n: int = 600,
    sample_prefix: str | None = None,
    genome_build: str = "GRCh38",
    annotation_version: str = "snpeff-5.1",
    mask_id: str = "pLoF",
    mask_hash: str = "hash-plof-v1",
) -> ScoreCovarianceArtifact:
    """Build a synthetic binary-trait stratum artifact.

    ``gene_keys`` maps gene -> list of canonical variant keys. An optional
    single-gene burden effect can be injected with a given sign.
    """
    rng = np.random.default_rng(seed)
    all_keys = sorted({k for keys in gene_keys.values() for k in keys})
    key_pos = {k: i for i, k in enumerate(all_keys)}
    p_total = len(all_keys)
    mafs = rng.uniform(0.01, 0.05, size=p_total)
    geno_all = rng.binomial(2, mafs[np.newaxis, :], size=(n, p_total)).astype(np.float64)

    cov = np.column_stack([rng.normal(size=n), rng.integers(0, 2, size=n).astype(float)])
    logit = -1.6 + 0.25 * cov[:, 0]
    if effect_gene is not None and effect_sign != 0.0:
        cols = [key_pos[k] for k in gene_keys[effect_gene]]
        burden = geno_all[:, cols].sum(axis=1)
        scale = burden.std() or 1.0
        logit = logit + effect_sign * effect_size * burden / scale
    prob = 1.0 / (1.0 + np.exp(-logit))
    y = (rng.uniform(size=n) < prob).astype(np.float64)

    genes: dict[str, tuple[np.ndarray, np.ndarray, list[str]]] = {}
    for gene, keys in gene_keys.items():
        cols = [key_pos[k] for k in keys]
        g = geno_all[:, cols]
        genes[gene] = (g, g.copy(), list(keys))

    prefix = sample_prefix if sample_prefix is not None else stratum_id
    sample_hashes = [f"{prefix}_s{i}" for i in range(n)]

    return build_stratum_artifact(
        stratum_id=stratum_id,
        genes=genes,
        phenotype=y,
        covariates=cov,
        trait_type="binary",
        genome_build=genome_build,
        annotation_version=annotation_version,
        mask_id=mask_id,
        mask_hash=mask_hash,
        n_cases=int(y.sum()),
        n_controls=int((1 - y).sum()),
        sample_id_hashes=sample_hashes,
    )
