"""ClinVar PM5 annotation — runtime module.

Loads a pre-built PM5 lookup table (from ``variantcentrifuge-build-pm5``)
and annotates variants with ACMG PM5 criterion information.

PM5: "Novel missense change at an amino acid residue where a different
missense change determined to be pathogenic has been seen before."

Three strength levels:
- PM5_Strong  — 2+ independent Pathogenic variants at same residue
- PM5         — 1 Pathogenic variant at same residue
- PM5_Supporting — only Likely_pathogenic evidence
"""

from __future__ import annotations

import gzip
import json
import logging
import re
from typing import Any

import pandas as pd

logger = logging.getLogger("variantcentrifuge")

# ---------------------------------------------------------------------------
# Amino acid code tables
# ---------------------------------------------------------------------------
AA_3TO1: dict[str, str] = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
}
AA_1TO3: dict[str, str] = {v: k for k, v in AA_3TO1.items()}

# ---------------------------------------------------------------------------
# HGVS protein notation parsers
# ---------------------------------------------------------------------------
# 3-letter: p.Arg123His, p.(Arg123His), p.Arg123*
# The trailing (?![a-z]) prevents matching frameshifts like p.Arg123Profs*45
HGVS_P_3LETTER = re.compile(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)(?![a-z])\)?")
# 1-letter: p.R123H, p.(R123H), p.R123*
HGVS_P_1LETTER = re.compile(r"p\.\(?([A-Z])(\d+)([A-Z]|\*)(?![a-z])\)?")

NOT_TRIGGERED: dict[str, Any] = {
    "triggered": False,
    "strength": None,
    "evidence_count": 0,
    "known_variants": "",
}


def normalize_aa(aa: str) -> str:
    """Normalize amino acid to 1-letter code for comparison."""
    return AA_3TO1.get(aa, aa)


def parse_hgvs_p(hgvs_p: str) -> tuple[str, int, str] | None:
    """Extract (ref_aa, position, alt_aa) from HGVS protein notation.

    Returns None for unparseable notation (p.?, p.=, frameshifts, etc.).
    """
    if not hgvs_p or hgvs_p in ("NA", ".", "p.?", "p.="):
        return None
    for pattern in (HGVS_P_3LETTER, HGVS_P_1LETTER):
        m = pattern.search(str(hgvs_p))
        if m:
            return m.group(1), int(m.group(2)), m.group(3)
    return None


def check_pm5(
    gene: str,
    ref_aa: str,
    aa_pos: int,
    alt_aa: str,
    lookup: dict[str, list[dict[str, str]]],
) -> dict[str, Any]:
    """Check ACMG PM5 criterion for a variant.

    Parameters
    ----------
    gene : str
        Gene symbol.
    ref_aa : str
        Reference amino acid (1- or 3-letter).
    aa_pos : int
        Amino acid position.
    alt_aa : str
        Alternate amino acid (1- or 3-letter).
    lookup : dict
        PM5 lookup table: ``{"GENE:pos": [{"ref_aa", "alt_aa", "clnsig", "review"}]}``.

    Returns
    -------
    dict
        Keys: triggered, strength, evidence_count, known_variants.
    """
    key = f"{gene}:{aa_pos}"
    entries = lookup.get(key, [])
    if not entries:
        return NOT_TRIGGERED.copy()

    query_alt = normalize_aa(alt_aa)
    query_ref = normalize_aa(ref_aa)

    evidence: list[dict[str, str]] = []
    for entry in entries:
        entry_ref = normalize_aa(entry["ref_aa"])
        entry_alt = normalize_aa(entry["alt_aa"])

        # Reference AA must match
        if entry_ref != query_ref:
            continue

        # Same AA change = PS1, NOT PM5 — skip
        if entry_alt == query_alt:
            continue

        evidence.append(entry)

    if not evidence:
        return NOT_TRIGGERED.copy()

    p_count = sum(1 for e in evidence if e["clnsig"] == "P")

    if p_count >= 2:
        strength = "PM5_Strong"
    elif p_count >= 1:
        strength = "PM5"
    else:
        strength = "PM5_Supporting"

    known = ",".join(f"{e['ref_aa']}{aa_pos}{e['alt_aa']}({e['clnsig']})" for e in evidence)

    return {
        "triggered": True,
        "strength": strength,
        "evidence_count": len(evidence),
        "known_variants": known,
    }


def load_pm5_lookup(path: str) -> dict[str, list[dict[str, str]]]:
    """Load a PM5 lookup table from JSON or gzipped JSON.

    Parameters
    ----------
    path : str
        Path to the lookup file (.json or .json.gz).

    Returns
    -------
    dict
        The ``"lookup"`` section of the file.
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        data = json.load(f)

    meta = data.get("_meta", {})
    logger.info(
        "Loaded PM5 lookup: %d positions, %d entries, assembly=%s, built=%s",
        meta.get("position_count", "?"),
        meta.get("entry_count", "?"),
        meta.get("assembly", "?"),
        meta.get("build_date", "?"),
    )
    result: dict[str, list[dict[str, str]]] = data.get("lookup", {})
    return result


def annotate_pm5(
    df: pd.DataFrame,
    lookup: dict[str, list[dict[str, str]]],
    gene_col: str = "GENE",
    hgvs_col: str = "HGVS_P",
) -> pd.DataFrame:
    """Add PM5 annotation columns to a variant DataFrame.

    Adds three columns: ``PM5``, ``PM5_evidence_count``, ``PM5_known_variants``.

    Parameters
    ----------
    df : pd.DataFrame
        Variant table.
    lookup : dict
        PM5 lookup from :func:`load_pm5_lookup`.
    gene_col : str
        Column name for gene symbol.
    hgvs_col : str
        Column name for HGVS protein notation.

    Returns
    -------
    pd.DataFrame
        The input DataFrame with PM5 columns added.
    """
    results: list[dict[str, Any]] = []

    for _, row in df.iterrows():
        gene = str(row.get(gene_col, ""))
        hgvs_p = str(row.get(hgvs_col, ""))
        parsed = parse_hgvs_p(hgvs_p)

        if not parsed or not gene:
            results.append(NOT_TRIGGERED.copy())
            continue

        ref_aa, aa_pos, alt_aa = parsed
        results.append(check_pm5(gene, ref_aa, aa_pos, alt_aa, lookup))

    df = df.copy()
    if not results:
        df["PM5"] = pd.Series(dtype=str)
        df["PM5_evidence_count"] = pd.Series(dtype=int)
        df["PM5_known_variants"] = pd.Series(dtype=str)
        return df

    pm5_df = pd.DataFrame(results)
    df["PM5"] = pm5_df["strength"].fillna("").values
    df["PM5_evidence_count"] = pm5_df["evidence_count"].fillna(0).astype(int).values
    df["PM5_known_variants"] = pm5_df["known_variants"].fillna("").values
    return df
