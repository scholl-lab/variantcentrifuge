"""Gene sets as extra association testing units (``--gene-set-file``, #201).

Each set is relabelled into its own ``GENE`` value (the set id) holding the union of its
member genes' qualifying variants, so every downstream step (aggregation, genotype
matrices, tests, score artifact) treats it as an ordinary testing unit.
"""

from __future__ import annotations

import dataclasses
import logging
from typing import Any

import pandas as pd

from ..gene_burden import _VARIANT_KEY

logger = logging.getLogger("variantcentrifuge")


def load_gene_sets(path: str) -> dict[str, list[str]]:
    """Read a ``set_id``/``gene`` TSV into ``{set_id: [genes]}`` (file order kept)."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = sorted({"set_id", "gene"} - set(df.columns))
    if missing:
        raise ValueError(f"--gene-set-file {path}: missing column(s) {missing}")
    df = df[["set_id", "gene"]].apply(lambda col: col.str.strip())
    if (df["set_id"] == "").any():
        raise ValueError(f"--gene-set-file {path}: empty set_id")
    if df["set_id"].str.contains("::", regex=False).any():
        raise ValueError(f"--gene-set-file {path}: set_id must not contain '::'")
    dup = df[df.duplicated()]
    if not dup.empty:
        raise ValueError(
            f"--gene-set-file {path}: duplicate set_id/gene rows {dup.values.tolist()[:5]}"
        )
    sets: dict[str, list[str]] = {}
    for set_id, gene in zip(df["set_id"], df["gene"], strict=True):
        sets.setdefault(set_id, [])
        if gene:
            sets[set_id].append(gene)
    return sets


def add_gene_set_units(
    df: pd.DataFrame, gene_sets: dict[str, list[str]]
) -> tuple[pd.DataFrame, list[str]]:
    """Append one relabelled copy of each set's member-gene rows; return (df, set ids added)."""
    genes = df["GENE"].astype(str)
    data_genes = set(genes)
    clash = sorted(set(gene_sets) & data_genes)
    if clash:
        raise ValueError(f"Gene set id(s) {clash} equal gene symbols in the data; rename the sets")
    absent = {g for members in gene_sets.values() for g in members} - data_genes
    if absent:
        logger.warning(
            "%d gene(s) from --gene-set-file are unknown or have no qualifying variants: %s",
            len(absent),
            sorted(absent)[:10],
        )
    key = _VARIANT_KEY if set(_VARIANT_KEY) <= set(df.columns) else None
    parts, added = [], []
    for set_id, members in gene_sets.items():
        rows = df[genes.isin(members)]
        if rows.empty:
            logger.warning(
                "Gene set '%s' has no member gene with qualifying variants; skipped", set_id
            )
            continue
        # A variant annotated to several member genes enters the set once.
        parts.append(rows.assign(GENE=set_id).drop_duplicates(subset=key))
        added.append(set_id)
    if not parts:
        return df, added
    logger.info("Testing %d gene set(s) as extra units: %s", len(added), added)
    return pd.concat([df, *parts], ignore_index=True), added


def run_gene_and_set_units(
    engine: Any,
    test_names: list[str],
    config: Any,
    units: list[dict],
    set_ids: list[str],
    gene_hook: Any = None,
) -> pd.DataFrame:
    """Run genes and sets as separate multiple-testing families; tag rows with ``unit_type``.

    Genes keep exactly the correction they get without sets. Sets are corrected among
    themselves, unweighted (``--gene-prior-weights`` holds per-gene weights).
    ``gene_hook`` runs for both; all its results end up in ``engine.gene_hook_results``.
    """
    from .engine import AssociationEngine

    is_set = set(set_ids)
    set_engine = AssociationEngine.from_names(
        test_names, dataclasses.replace(config, gene_prior_weights=None)
    )
    set_engine.unit_label = "gene sets"
    frames = [
        engine.run_all([d for d in units if d.get("GENE") not in is_set], gene_hook).assign(
            unit_type="gene"
        ),
        set_engine.run_all([d for d in units if d.get("GENE") in is_set], gene_hook).assign(
            unit_type="set"
        ),
    ]
    engine.gene_hook_results.update(set_engine.gene_hook_results)
    engine.run_info += set_engine.run_info
    frames = [f for f in frames if "gene" in f.columns]
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    out.insert(1, "unit_type", out.pop("unit_type"))
    return out
