"""Build a PM5 lookup table from ClinVar variant_summary.txt.gz.

This module provides both a library API and a CLI entry point
(``variantcentrifuge-build-pm5``) for creating the pre-built lookup
table consumed by :mod:`variantcentrifuge.clinvar_pm5` at runtime.

Usage::

    # Download and build in one step
    variantcentrifuge-build-pm5 --download --assembly GRCh38 --output pm5.json.gz

    # From a local file
    variantcentrifuge-build-pm5 --input variant_summary.txt.gz --output pm5.json.gz
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import json
import logging
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger("variantcentrifuge")

CLINVAR_FTP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

# ClinVar review status → star count
REVIEW_STARS: dict[str, int] = {
    "no assertion provided": 0,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no classification for the single variant": 0,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}

# Regex to extract protein change from ClinVar Name field
# Example: "NM_000059.4(BRCA2):c.8167G>C (p.Asp2723His)"
PROTEIN_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)\)")


def _review_stars(status: str) -> int:
    """Map a ReviewStatus string to star count."""
    return REVIEW_STARS.get(status.strip(), 0)


def _simplify_clnsig(clnsig: str) -> str:
    """Simplify clinical significance to P or LP."""
    s = clnsig.strip()
    if "Pathogenic/Likely_pathogenic" in s:
        return "P"
    if "Pathogenic" in s:
        return "P"
    if "Likely_pathogenic" in s or "Likely pathogenic" in s:
        return "LP"
    return s


def download_variant_summary(output_path: str) -> str:
    """Download variant_summary.txt.gz from ClinVar FTP.

    Parameters
    ----------
    output_path : str
        Local path to save the downloaded file.

    Returns
    -------
    str
        Path to the downloaded file.
    """
    logger.info("Downloading ClinVar variant_summary.txt.gz ...")
    urllib.request.urlretrieve(CLINVAR_FTP_URL, output_path)
    logger.info("Downloaded to %s", output_path)
    return output_path


def build_pm5_lookup(
    input_path: str,
    assembly: str = "GRCh38",
    min_stars: int = 1,
    include_lp: bool = True,
) -> dict[str, Any]:
    """Build PM5 lookup table from variant_summary.txt.gz.

    Parameters
    ----------
    input_path : str
        Path to ClinVar variant_summary.txt.gz.
    assembly : str
        Genome assembly to filter on.
    min_stars : int
        Minimum review star count (0-4).
    include_lp : bool
        Include Likely_pathogenic entries for PM5_Supporting.

    Returns
    -------
    dict
        Complete lookup structure with ``_meta`` and ``lookup`` keys.
    """
    logger.info("Reading %s ...", input_path)
    use_cols = [
        "#AlleleID",
        "Type",
        "Name",
        "GeneSymbol",
        "ClinicalSignificance",
        "Origin",
        "Assembly",
        "ReviewStatus",
    ]
    df = pd.read_csv(
        input_path,
        sep="\t",
        compression="gzip",
        usecols=use_cols,
        dtype=str,
        na_filter=False,
    )

    rows_read = len(df)
    logger.info("Read %d rows", rows_read)

    # Step 1: Filter
    mask = (
        (df["Type"] == "single nucleotide variant")
        & (df["Assembly"] == assembly)
        & (df["ClinicalSignificance"].str.contains("athogenic", case=False))
        & (~df["ClinicalSignificance"].str.contains("Conflicting", case=False))
        & (df["Origin"].str.contains("germline", case=False))
        & (df["ReviewStatus"].map(_review_stars) >= min_stars)
    )
    if not include_lp:
        mask = mask & df["ClinicalSignificance"].str.contains("Pathogenic")
        mask = mask & ~df["ClinicalSignificance"].str.contains("Likely_pathogenic")

    filtered = df[mask].copy()
    logger.info("After filtering: %d rows", len(filtered))

    # Step 2: Parse protein change from Name column
    parsed = filtered["Name"].str.extract(PROTEIN_RE.pattern)
    parsed.columns = ["ref_aa", "aa_pos", "alt_aa"]  # type: ignore[assignment]
    parsed_mask = parsed["ref_aa"].notna()
    result = pd.concat(
        [filtered[parsed_mask].reset_index(drop=True), parsed[parsed_mask].reset_index(drop=True)],
        axis=1,
    )

    parse_success = len(result)
    parse_fail = len(filtered) - parse_success
    logger.info("Parse succeeded: %d, failed: %d", parse_success, parse_fail)

    # Step 3: Build lookup dict
    lookup: dict[str, list[dict[str, str]]] = defaultdict(list)
    for _, row in result.iterrows():
        key = f"{row['GeneSymbol']}:{row['aa_pos']}"
        entry = {
            "ref_aa": row["ref_aa"],
            "alt_aa": row["alt_aa"],
            "clnsig": _simplify_clnsig(row["ClinicalSignificance"]),
            "review": row["ReviewStatus"],
        }
        if entry not in lookup[key]:
            lookup[key].append(entry)

    # Step 4: Build output
    output: dict[str, Any] = {
        "_meta": {
            "source": "ClinVar variant_summary.txt.gz",
            "assembly": assembly,
            "min_review_stars": min_stars,
            "include_lp": include_lp,
            "build_date": datetime.date.today().isoformat(),
            "entry_count": sum(len(v) for v in lookup.values()),
            "position_count": len(lookup),
            "rows_read": rows_read,
            "rows_after_filter": len(filtered),
            "parse_success": parse_success,
            "parse_fail": parse_fail,
        },
        "lookup": dict(lookup),
    }
    return output


def write_lookup(data: dict[str, Any], output_path: str) -> None:
    """Write lookup table to JSON or gzipped JSON."""
    if output_path.endswith(".gz"):
        with gzip.open(output_path, "wt", encoding="utf-8") as f:
            json.dump(data, f, separators=(",", ":"))
    else:
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
    logger.info("Wrote lookup to %s", output_path)


def print_stats(data: dict[str, Any]) -> None:
    """Print summary statistics for the lookup table."""
    meta = data["_meta"]
    lookup = data["lookup"]

    counts = [len(v) for v in lookup.values()]
    single = sum(1 for c in counts if c == 1)
    double = sum(1 for c in counts if c == 2)
    triple_plus = sum(1 for c in counts if c >= 3)

    genes: dict[str, int] = defaultdict(int)
    for key in lookup:
        gene = key.rsplit(":", 1)[0]
        genes[gene] += 1

    top_genes = sorted(genes.items(), key=lambda x: -x[1])[:5]

    denom = max(meta["rows_after_filter"], 1)
    pct = meta["parse_success"] / denom * 100

    print(f"""
PM5 Lookup Table Statistics
===========================
Source:           {meta["source"]}
Assembly:         {meta["assembly"]}
Min review stars: {meta["min_review_stars"]}
Include LP:       {"yes" if meta["include_lp"] else "no"}
Build date:       {meta["build_date"]}

Rows read:        {meta["rows_read"]:,}
After filtering:  {meta["rows_after_filter"]:,}
Parse succeeded:  {meta["parse_success"]:,} ({pct:.1f}%)
Parse failed:     {meta["parse_fail"]:,}

Unique positions: {meta["position_count"]:,}
Unique genes:     {len(genes):,}

Positions with 1 variant:  {single:,} ({single / max(len(lookup), 1) * 100:.1f}%)
Positions with 2 variants: {double:,} ({double / max(len(lookup), 1) * 100:.1f}%)
Positions with 3+ variants: {triple_plus:,} ({triple_plus / max(len(lookup), 1) * 100:.1f}%)

Top genes by position count:""")
    for gene, count in top_genes:
        print(f"  {gene:<8s} {count}")


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser for the build tool."""
    parser = argparse.ArgumentParser(
        prog="variantcentrifuge-build-pm5",
        description="Build a PM5 lookup table from ClinVar variant_summary.txt.gz",
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Path to local variant_summary.txt.gz",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download latest variant_summary.txt.gz from ClinVar FTP",
    )
    parser.add_argument(
        "--assembly",
        type=str,
        default="GRCh38",
        choices=["GRCh38", "GRCh37"],
        help="Genome assembly to filter on (default: GRCh38)",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output path for lookup JSON (.json or .json.gz)",
    )
    parser.add_argument(
        "--min-review-stars",
        type=int,
        default=1,
        choices=[0, 1, 2, 3, 4],
        help="Minimum ClinVar review stars (default: 1)",
    )
    parser.add_argument(
        "--include-lp",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include Likely_pathogenic entries (default: True).",
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Print summary statistics after build",
    )
    return parser


def main(args: list[str] | None = None) -> int:
    """CLI entry point for variantcentrifuge-build-pm5."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )

    parser = create_parser()
    parsed = parser.parse_args(args)

    if not parsed.input and not parsed.download:
        parser.error("One of --input or --download is required.")

    input_path = parsed.input
    if parsed.download:
        import tempfile

        if not input_path:
            input_path = str(Path(tempfile.gettempdir()) / "variant_summary.txt.gz")
        download_variant_summary(input_path)

    data = build_pm5_lookup(
        input_path,
        assembly=parsed.assembly,
        min_stars=parsed.min_review_stars,
        include_lp=parsed.include_lp,
    )

    write_lookup(data, parsed.output)

    if parsed.stats:
        print_stats(data)

    return 0


if __name__ == "__main__":
    sys.exit(main())
