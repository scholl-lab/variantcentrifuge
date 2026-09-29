# File: variantcentrifuge/association/meta/align.py
"""
Variant alignment across strata and pooled-frequency weight harmonisation.

Two strata sequenced/called separately observe different variant sets within a
gene. Meta-analysis takes the union of variants (by canonical key) and embeds each
stratum's score ``U`` and covariance ``V`` into that union space, with zeros for
variants a stratum did not observe (a variant monomorphic in one stratum
correctly contributes no information there — it is kept, never dropped).

Weights are harmonised once from the POOLED ALT-allele frequency ``ΣAC/ΣAN`` so
every stratum uses the same per-variant weight (cohort-specific MAF weights would
put ``U``/``V`` on incompatible scales). The Beta weight uses ``MAF=min(af,1-af)``.
"""

from __future__ import annotations

import logging

import numpy as np
import scipy.stats

from variantcentrifuge.association.meta.artifact import GeneScore

logger = logging.getLogger("variantcentrifuge")


def canonical_key(chrom: str, pos: int | str, ref: str, alt: str) -> str:
    """Canonical variant key ``chrom:pos:ref:alt`` (chrom normalised, alleles upper).

    Strips a leading ``chr`` prefix so ``chr1`` and ``1`` align. Alleles are
    upper-cased. Callers are responsible for prior normalisation/left-alignment.
    """
    c = str(chrom)
    if c.lower().startswith("chr"):
        c = c[3:]
    return f"{c}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def align_strata(
    gene_scores: list[GeneScore],
) -> tuple[list[np.ndarray], list[np.ndarray], np.ndarray, np.ndarray, list[str]]:
    """Embed each stratum's ``U``/``V`` into the shared union of variant keys.

    Returns
    -------
    (u_list, v_list, ac_pooled, an_pooled, union_keys)
        ``u_list[c]`` has shape (m,), ``v_list[c]`` shape (m, m) where m = |union|;
        ``ac_pooled``/``an_pooled`` are summed across strata (shape (m,)).
    """
    if not gene_scores:
        raise ValueError("align_strata requires at least one GeneScore")

    union_keys = sorted({k for gs in gene_scores for k in gs.variant_keys})
    index = {k: i for i, k in enumerate(union_keys)}
    m = len(union_keys)

    u_list: list[np.ndarray] = []
    v_list: list[np.ndarray] = []
    ac_pooled = np.zeros(m, dtype=np.float64)
    an_pooled = np.zeros(m, dtype=np.float64)

    for gs in gene_scores:
        cols = np.array([index[k] for k in gs.variant_keys], dtype=int)
        u_full = np.zeros(m, dtype=np.float64)
        v_full = np.zeros((m, m), dtype=np.float64)
        u_full[cols] = gs.u
        v_full[np.ix_(cols, cols)] = gs.v
        u_list.append(u_full)
        v_list.append(v_full)
        ac_pooled[cols] += gs.ac
        an_pooled[cols] += gs.an

    return u_list, v_list, ac_pooled, an_pooled, union_keys


def pooled_frequency_weights(
    ac_pooled: np.ndarray,
    an_pooled: np.ndarray,
    scheme: str = "beta:1,25",
) -> np.ndarray:
    """Harmonised per-variant weights from pooled frequency.

    Only ``beta:a,b`` and ``uniform`` are supported — the only schemes reproducible
    from the persisted artifacts. CADD/REVEL/score-column schemes are rejected
    (they need per-variant annotation provenance not stored in the artifact).
    """
    ac = np.asarray(ac_pooled, dtype=np.float64)
    an = np.asarray(an_pooled, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        af = np.where(an > 0, ac / an, 0.0)
    maf = np.minimum(af, 1.0 - af)

    scheme = scheme.strip().lower()
    if scheme == "uniform":
        return np.ones_like(maf)
    if scheme.startswith("beta:"):
        try:
            a_str, b_str = scheme.split(":", 1)[1].split(",")
            a, b = float(a_str), float(b_str)
        except ValueError as exc:
            raise ValueError(
                f"Malformed beta weight scheme '{scheme}' (expected beta:a,b)"
            ) from exc
        return np.asarray(scipy.stats.beta.pdf(maf, a, b), dtype=np.float64)
    raise ValueError(
        f"Unsupported meta weight scheme '{scheme}'. Meta-analysis supports only "
        "'beta:a,b' and 'uniform' (functional-score weights are not reproducible "
        "from score artifacts)."
    )


__all__ = ["align_strata", "canonical_key", "pooled_frequency_weights"]
