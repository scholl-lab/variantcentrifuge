# File: variantcentrifuge/association/meta/runner.py
"""
End-to-end meta-analysis orchestration: manifest → per-gene combination → TSV.

For each gene present in the study, align the contributing strata, harmonise
weights from pooled frequency, run the common-effect burden/SKAT/SKAT-O meta,
compute heterogeneity and leave-one-out, then apply exome-wide multiplicity
correction to the primary (SKAT-O) p-value.

Each mask (all its artifacts share ``mask_hash``, enforced by comparability
validation) is combined on its own; with several masks in the manifest
(``mask_id`` column) the per-gene primaries are also ACAT-combined (#196).
"""

from __future__ import annotations

import logging
import os
import re

import numpy as np
import pandas as pd

from variantcentrifuge.association.correction import apply_correction
from variantcentrifuge.association.diagnostics import write_diagnostics
from variantcentrifuge.association.meta.align import align_strata, pooled_frequency_weights
from variantcentrifuge.association.meta.artifact import GeneScore, ScoreCovarianceArtifact
from variantcentrifuge.association.meta.combine import meta_burden, meta_skat, meta_skato
from variantcentrifuge.association.meta.heterogeneity import cochran_q, per_stratum_burden
from variantcentrifuge.association.meta.manifest import (
    group_by_mask,
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
    "loo_burden_by_stratum",
    "loo_skato_by_stratum",
    "loo_max_burden_p",
    "loo_max_burden_stratum",
    "loo_max_skato_p",
    "loo_max_skato_stratum",
    "cmac_total",
    "cmac_by_stratum",
    "n_carriers_total",
    "carriers_by_stratum",
    "primary_p",
    "primary_q",
    "status",
    "warnings",
]


# Columns that do not depend on the weight scheme; the rest get a per-scheme suffix (#202).
_SCHEME_FREE_COLUMNS = [
    "gene",
    "mask",
    "n_strata",
    "n_cases_total",
    "n_controls_total",
    "n_variants_union",
    "cmac_total",
    "cmac_by_stratum",
    "n_carriers_total",
    "carriers_by_stratum",
    "warnings",
]


# Meta p-value columns -> test ids for write_diagnostics (expects ``{test}_pvalue``).
_DIAG_TESTS = {
    "meta_burden_p": "meta_burden",
    "meta_skat_p": "meta_skat",
    "meta_skato_p": "meta_skato",
    "primary_p": "primary",
}


def _fmt(x: float | None) -> str:
    return "NA" if x is None else f"{x:.10g}"


def _mac(ac: np.ndarray, an: np.ndarray, v: np.ndarray) -> float:
    """Cumulative minor-allele count: per-variant min(AC, AN - AC), summed over the
    variants that carry score information (diag(V) > 0); others add nothing to the test.
    """
    info = np.diag(v) > 0
    return float(np.minimum(ac[info], an[info] - ac[info]).sum())


def _leave_one_out(
    gene_scores: list[GeneScore], stratum_ids: list[str], weight_scheme: str
) -> dict:
    """Named leave-one-stratum-out burden and SKAT-O columns.

    Each drop is a full meta-analysis of the remaining strata (re-aligned, weights
    re-derived from their pooled frequency), i.e. identical to running ``meta`` on
    the manifest without that stratum.
    """
    if len(gene_scores) < 2:
        return {}
    burden: list[tuple[str, float | None, int]] = []
    skato: list[tuple[str, float | None]] = []
    for i, sid in enumerate(stratum_ids):
        u, v, ac, an, _keys = align_strata(gene_scores[:i] + gene_scores[i + 1 :])
        w = pooled_frequency_weights(ac, an, weight_scheme)
        p, _beta, direction = meta_burden(u, v, w)
        burden.append((sid, p, direction))
        skato.append((sid, meta_skato(u, v, w)[0]))

    row: dict = {
        "loo_burden_by_stratum": ";".join(f"{s}:{_fmt(p)}:{d:+d}" for s, p, d in burden),
        "loo_skato_by_stratum": ";".join(f"{s}:{_fmt(p)}" for s, p in skato),
    }
    for test, pairs in (("burden", [(s, p) for s, p, _d in burden]), ("skato", skato)):
        valid = [(s, p) for s, p in pairs if p is not None]
        if valid:
            sid, p = max(valid, key=lambda sp: sp[1])
            row[f"loo_max_{test}_p"], row[f"loo_max_{test}_stratum"] = p, sid
    burden_ps = [p for _s, p, _d in burden if p is not None]
    row["loo_min_p"] = min(burden_ps) if burden_ps else None
    return row


def _gene_row(
    gene: str,
    mask_id: str,
    contributors: list[tuple[ScoreCovarianceArtifact, GeneScore]],
    weight_scheme: str,
    min_mac: float = 0.0,
    n_strata_required: int = 0,
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
    directions = [
        f"{art.stratum_id}:{int(np.sign(b)) if np.isfinite(b) else 0}"
        for (art, _gs), b in zip(kept, betas, strict=True)
    ]

    stratum_macs = [_mac(gs.ac, gs.an, gs.v) for gs in gene_scores]
    # Sum of stratum MACs: a pooled min(AC, AN-AC) overcounts sites fixed ALT in a stratum
    cmac_total = float(sum(stratum_macs))
    carriers = [gs.n_carriers for gs in gene_scores]
    primary_p = skato_p if skato_p is not None else burden_p
    status = "ok" if n_strata >= 2 else "single_stratum"
    if primary_p is None:
        status = "degenerate"
    elif sum(m > 0 for m in stratum_macs) < n_strata_required:
        status = "not_in_all_strata"
    elif cmac_total < min_mac:
        status = "low_mac"

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
        "per_stratum_directions": ",".join(directions),
        **_leave_one_out(gene_scores, [art.stratum_id for art, _gs in kept], weight_scheme),
        "cmac_total": cmac_total,
        "cmac_by_stratum": ";".join(
            f"{art.stratum_id}:{m:.0f}" for (art, _gs), m in zip(kept, stratum_macs, strict=True)
        ),
        "n_carriers_total": None if None in carriers else sum(carriers),  # type: ignore[arg-type]
        "carriers_by_stratum": ";".join(
            f"{art.stratum_id}:{'NA' if c is None else c}"
            for (art, _gs), c in zip(kept, carriers, strict=True)
        ),
        "primary_p": primary_p,
        "primary_q": None,
        "status": status,
        "warnings": ";".join(warnings) or None,
    }


def meta_analyze_artifacts(
    artifacts: list[ScoreCovarianceArtifact],
    weight_scheme: str = "beta:1,25",
    correction: str = "fdr",
    min_mac: float = 0.0,
    require_variants_in_all_strata: bool = False,
    mask_name: str | None = None,
) -> pd.DataFrame:
    """Combine already-loaded, already-validated artifacts into a results frame.

    The ``mask`` column is ``mask_name`` (manifest ``mask_id``), else the artifacts'
    ``--mask-label``, else their mask hash prefix.

    Genes with cumulative MAC below ``min_mac`` (``status=low_mac``) or, with
    ``require_variants_in_all_strata``, without an informative minor allele in every
    stratum (``status=not_in_all_strata``) keep their p-values but stay outside the
    multiple-testing family, like the main engine's untestable genes.
    """
    if not (np.isfinite(min_mac) and min_mac >= 0):
        raise ValueError(f"min_mac must be a finite number >= 0, got {min_mac}")
    bad = [a.stratum_id for a in artifacts if set(a.stratum_id) & set(":;,")]
    if bad:
        raise ValueError(
            f"stratum_id must not contain the output separators ':', ';' or ',': {bad}"
        )
    mask_id = mask_name or artifacts[0].mask_label or artifacts[0].mask_id
    all_genes = sorted({g for art in artifacts for g in art.genes})

    rows: list[dict] = []
    for gene in all_genes:
        contributors = [(art, art.genes[gene]) for art in artifacts if gene in art.genes]
        n_required = len(artifacts) if require_variants_in_all_strata else 0
        rows.append(_gene_row(gene, mask_id, contributors, weight_scheme, min_mac, n_required))

    df = pd.DataFrame(rows).reindex(columns=OUTPUT_COLUMNS)

    testable = df["primary_p"].notna() & ~df["status"].isin(["low_mac", "not_in_all_strata"])
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


def _acat_r(p: np.ndarray) -> float:
    """Equal-weight ACAT exactly as the reference R/ACAT.R (yaowuliu/ACAT @ e737d7d).

    Like ``cauchy_combination`` (single-cohort, weighted; #216), p == 1 is kept and
    evaluated as float tan(-pi/2) ~ -1.6e16, and p < 1e-15 uses 1/(p*pi), as in R.
    """
    small = p < 1e-15
    with np.errstate(divide="ignore"):
        terms = np.where(small, 1.0 / (p * np.pi), np.tan((0.5 - p) * np.pi))
    # pcauchy(t, lower.tail=FALSE) without the 0.5 - atan(t)/pi cancellation at large t
    return float(np.arctan2(1.0, terms.mean()) / np.pi)


def add_acat_masks(df: pd.DataFrame, correction: str, suffix: str = "") -> pd.DataFrame:
    """Per-gene ACAT (Cauchy) combination of the masks' ``primary_p`` (#196).

    Only mask rows inside their correction family (``primary_q`` set) contribute.
    ``acat_masks_q`` corrects the ACAT p-values across genes (one test per gene).
    """
    p_col, q_col = f"primary_p{suffix}", f"primary_q{suffix}"
    in_family = df[df[q_col].notna()]
    acat: dict[str, float] = {}
    undefined: list[str] = []
    for gene, rows in in_family.groupby("gene", sort=False):
        ps = rows[p_col].to_numpy(dtype=float)
        if (ps == 0).any() and (ps == 1).any():
            undefined.append(gene)  # R ACAT stops here: the combination is undefined
            continue
        acat[gene] = _acat_r(ps)
    if undefined:
        logger.warning(
            "acat_masks_p%s left empty for %d gene(s) whose masks have p-values of exactly "
            "0 and 1 (ACAT undefined): %s",
            suffix,
            len(undefined),
            ", ".join(undefined),
        )
        flag = f"acat_masks_undefined{suffix}"
        hit = df["gene"].isin(undefined)
        df.loc[hit, "warnings"] = [
            w + ";" + flag if isinstance(w, str) and w else flag for w in df.loc[hit, "warnings"]
        ]
    q = apply_correction(np.array(list(acat.values()), dtype=float), correction) if acat else ()
    df[f"acat_masks_p{suffix}"] = df["gene"].map(acat)
    df[f"acat_masks_q{suffix}"] = df["gene"].map(dict(zip(acat, q, strict=True)))
    return df


def _write_meta_diagnostics(
    df: pd.DataFrame, out_dir: str, artifacts: list[ScoreCovarianceArtifact]
) -> None:
    """lambda/QQ over exactly the correction family; all rows kept for warning counts."""
    diag = df.rename(columns={c: f"{t}_pvalue" for c, t in _DIAG_TESTS.items()})
    diag.loc[df["primary_q"].isna(), [f"{t}_pvalue" for t in _DIAG_TESTS.values()]] = np.nan
    diag["warnings"] = diag["warnings"].fillna("")
    write_diagnostics(
        diag,
        out_dir,
        test_names=list(_DIAG_TESTS.values()),
        n_cases=sum(art.n_cases for art in artifacts),
        n_controls=sum(art.n_controls for art in artifacts),
        cohort_warnings=[],
    )


def _scheme_suffixes(schemes: list[str]) -> list[str]:
    """Column suffix per scheme: none for a single scheme, else e.g. ``_beta_1_25``."""
    if len(schemes) == 1:
        return [""]
    tags = [re.sub(r"[^0-9A-Za-z]+", "_", s).strip("_") for s in schemes]
    if len(set(tags)) != len(tags):
        raise ValueError(f"Weight schemes must be distinct, got {schemes}")
    return [f"_{t}" for t in tags]


def run_meta(
    manifest_path: str,
    output_path: str,
    weight_scheme: str | list[str] = "beta:1,25",
    correction: str = "fdr",
    min_mac: float = 0.0,
    diagnostics_dir: str | None = None,
    require_variants_in_all_strata: bool = False,
    allow_annotation_mismatch: bool = False,
) -> pd.DataFrame:
    """Full pipeline: manifest → validated artifacts → meta results → TSV.

    Each mask in the manifest is validated and combined separately; with several
    masks the rows are stacked and ``acat_masks_p``/``acat_masks_q`` added. With
    several weight schemes, every scheme-dependent column (and the ACAT pair) appears
    once per scheme, suffixed by the scheme name; each scheme is corrected on its own.
    With ``diagnostics_dir``, also writes λ_GC/QQ diagnostics over the genes inside the
    correction family (``primary_q`` set), in a subdirectory per mask and per scheme
    when there are several.
    """
    schemes = [weight_scheme] if isinstance(weight_scheme, str) else list(weight_scheme)
    suffixes = _scheme_suffixes(schemes)
    groups = group_by_mask(read_study_manifest(manifest_path))
    frames: list[pd.DataFrame] = []
    for mask, refs in groups.items():
        artifacts = load_manifest_artifacts(refs)
        validate_comparability(artifacts, allow_annotation_mismatch)
        logger.info(
            "Meta-analysis: %d strata, mask '%s'", len(artifacts), mask or artifacts[0].mask_id
        )
        per_scheme = [
            meta_analyze_artifacts(
                artifacts,
                weight_scheme=scheme,
                correction=correction,
                min_mac=min_mac,
                require_variants_in_all_strata=require_variants_in_all_strata,
                mask_name=mask,
            )
            for scheme in schemes
        ]
        if diagnostics_dir is not None:
            for sfx, scheme_df in zip(suffixes, per_scheme, strict=True):
                parts = [str(mask)] if len(groups) > 1 else []
                sub = os.path.join(diagnostics_dir, *parts, sfx[1:])
                _write_meta_diagnostics(scheme_df, sub, artifacts)
        df = per_scheme[0]
        if len(schemes) > 1:
            df = pd.concat(
                [df[_SCHEME_FREE_COLUMNS]]
                + [
                    f.drop(columns=_SCHEME_FREE_COLUMNS).add_suffix(sfx)
                    for f, sfx in zip(per_scheme, suffixes, strict=True)
                ],
                axis=1,
            )
        frames.append(df)
    df = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    if len(groups) > 1:
        for sfx in suffixes:
            df = add_acat_masks(df, correction, sfx)
    df.to_csv(output_path, sep="\t", index=False)
    logger.info("Wrote meta results to %s", output_path)
    return df


__all__ = ["OUTPUT_COLUMNS", "add_acat_masks", "meta_analyze_artifacts", "run_meta"]
