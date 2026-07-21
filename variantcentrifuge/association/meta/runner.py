# File: variantcentrifuge/association/meta/runner.py
"""
End-to-end meta-analysis orchestration: manifest → per-gene combination → TSV.

For each gene present in the study, align the contributing strata, harmonise
weights from pooled frequency, run the common-effect burden/SKAT/SKAT-O meta,
compute heterogeneity and leave-one-out, then apply exome-wide multiplicity
correction to the primary (SKAT-O) p-value.

This runner handles a single declared mask (all artifacts share ``mask_hash``,
enforced by comparability validation). Multi-mask designs run this once per mask
and ACAT-combine the per-gene primaries at a higher level.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from variantcentrifuge.association.correction import apply_correction
from variantcentrifuge.association.meta.align import align_strata, pooled_frequency_weights
from variantcentrifuge.association.meta.artifact import GeneScore, ScoreCovarianceArtifact
from variantcentrifuge.association.meta.combine import meta_burden, meta_skat, meta_skato
from variantcentrifuge.association.meta.heterogeneity import (
    cochran_q,
    leave_one_out,
    per_stratum_burden,
)
from variantcentrifuge.association.meta.manifest import (
    load_manifest_artifacts,
    read_study_manifest,
    validate_comparability,
)

logger = logging.getLogger("variantcentrifuge")

OUTPUT_COLUMNS = [
    "gene",
    "mask",
    "n_strata",
    "n_cases_total",
    "n_controls_total",
    "n_variants_union",
    "meta_burden_p",
    "meta_burden_beta",
    "meta_burden_direction",
    "meta_skat_p",
    "meta_skato_p",
    "meta_skato_rho",
    "het_q",
    "het_p",
    "het_i2",
    "per_stratum_directions",
    "loo_min_p",
    "primary_p",
    "primary_q",
    "status",
    "warnings",
]


def _gene_row(
    gene: str,
    mask_id: str,
    contributors: list[tuple[ScoreCovarianceArtifact, GeneScore]],
    weight_scheme: str,
) -> dict:
    """Combine one gene's contributing strata into a single result row."""
    warnings: list[str] = []
    kept = [(art, gs) for art, gs in contributors if gs.null_converged]
    dropped = [art.stratum_id for art, gs in contributors if not gs.null_converged]
    if dropped:
        warnings.append(f"excluded_non_converged:{','.join(sorted(dropped))}")

    n_strata = len(kept)
    if n_strata == 0:
        return {
            "gene": gene,
            "mask": mask_id,
            "n_strata": 0,
            "status": "no_converged_strata",
            "warnings": ";".join(warnings) or None,
        }

    gene_scores = [gs for _art, gs in kept]
    u_list, v_list, ac_pooled, an_pooled, union_keys = align_strata(gene_scores)
    weights = pooled_frequency_weights(ac_pooled, an_pooled, weight_scheme)

    burden_p, beta, direction = meta_burden(u_list, v_list, weights)
    skat_p, _method = meta_skat(u_list, v_list, weights)
    skato_p, rho = meta_skato(u_list, v_list, weights)

    betas, ses = per_stratum_burden(u_list, v_list, weights)
    het_q, het_p, het_i2 = cochran_q(betas, ses)
    directions = [int(np.sign(b)) if np.isfinite(b) else 0 for b in betas]

    loo = leave_one_out(
        u_list,
        v_list,
        weights,
        [art.stratum_id for art, _gs in kept],
        lambda u, v, w: meta_burden(u, v, w)[0],
    )
    loo_ps = [p for _sid, p in loo if p is not None]
    loo_min_p = min(loo_ps) if loo_ps else None

    primary_p = skato_p if skato_p is not None else burden_p
    status = "ok" if n_strata >= 2 else "single_stratum"
    if primary_p is None:
        status = "degenerate"

    return {
        "gene": gene,
        "mask": mask_id,
        "n_strata": n_strata,
        "n_cases_total": int(sum(art.n_cases for art, _gs in kept)),
        "n_controls_total": int(sum(art.n_controls for art, _gs in kept)),
        "n_variants_union": len(union_keys),
        "meta_burden_p": burden_p,
        "meta_burden_beta": beta,
        "meta_burden_direction": direction,
        "meta_skat_p": skat_p,
        "meta_skato_p": skato_p,
        "meta_skato_rho": rho,
        "het_q": het_q,
        "het_p": het_p,
        "het_i2": het_i2,
        "per_stratum_directions": ",".join(str(d) for d in directions),
        "loo_min_p": loo_min_p,
        "primary_p": primary_p,
        "primary_q": None,
        "status": status,
        "warnings": ";".join(warnings) or None,
    }


def meta_analyze_artifacts(
    artifacts: list[ScoreCovarianceArtifact],
    weight_scheme: str = "beta:1,25",
    correction: str = "fdr",
) -> pd.DataFrame:
    """Combine already-loaded, already-validated artifacts into a results frame."""
    mask_id = artifacts[0].mask_id
    all_genes = sorted({g for art in artifacts for g in art.genes})

    rows: list[dict] = []
    for gene in all_genes:
        contributors = [(art, art.genes[gene]) for art in artifacts if gene in art.genes]
        rows.append(_gene_row(gene, mask_id, contributors, weight_scheme))

    df = pd.DataFrame(rows).reindex(columns=OUTPUT_COLUMNS)

    testable = df["primary_p"].notna()
    if testable.any():
        corrected = apply_correction(
            df.loc[testable, "primary_p"].to_numpy(dtype=float), correction
        )
        df.loc[testable, "primary_q"] = corrected

    n_sig = int((df["primary_q"] < 0.05).sum())
    logger.info(
        "Meta-analysis complete: %d genes, %d with valid primary, %d significant (q<0.05)",
        len(df),
        int(testable.sum()),
        n_sig,
    )
    return df


def run_meta(
    manifest_path: str,
    output_path: str,
    weight_scheme: str = "beta:1,25",
    correction: str = "fdr",
) -> pd.DataFrame:
    """Full pipeline: manifest → validated artifacts → meta results → TSV."""
    refs = read_study_manifest(manifest_path)
    artifacts = load_manifest_artifacts(refs)
    validate_comparability(artifacts)
    logger.info("Meta-analysis: %d strata, mask '%s'", len(artifacts), artifacts[0].mask_id)

    df = meta_analyze_artifacts(artifacts, weight_scheme=weight_scheme, correction=correction)
    df.to_csv(output_path, sep="\t", index=False)
    logger.info("Wrote meta results to %s", output_path)
    return df


__all__ = ["OUTPUT_COLUMNS", "meta_analyze_artifacts", "run_meta"]
