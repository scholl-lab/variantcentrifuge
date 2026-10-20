"""Feature loading and file parsing utilities for custom annotations."""

from __future__ import annotations

import json
import logging
import os
import re
from collections.abc import Iterator
from typing import Any

logger = logging.getLogger("variantcentrifuge")

__all__ = [
    "_load_bed_files",
    "_load_gene_lists",
    "_load_json_gene_data",
    "_load_named_bed_files",
    "_parse_bed_line",
    "_sanitize_name",
    "split_bed_spec",
]


def _sanitize_name(name: str) -> str:
    """Sanitize a name to be safe for use in annotations.

    Parameters
    ----------
    name : str
        Raw name string

    Returns
    -------
    str
        Sanitized name with only alphanumeric characters and underscores
    """
    base_name = os.path.splitext(os.path.basename(name))[0]
    sanitized = re.sub(r"[^a-zA-Z0-9_]", "_", base_name)
    sanitized = re.sub(r"_+", "_", sanitized).strip("_")
    return sanitized if sanitized else "unnamed"


def _parse_bed_line(line: str) -> tuple[str, int, int, str] | None:
    """Parse a single BED file line.

    Parameters
    ----------
    line : str
        BED format line

    Returns
    -------
    tuple or None
        Tuple of (chrom, start, end, name) or None if invalid
    """
    line = line.strip()
    if not line or line.startswith(("#", "track", "browser")):
        return None

    parts = line.split("\t")
    if len(parts) < 3:
        return None

    try:
        chrom = parts[0].removeprefix("chr")
        start = int(parts[1])
        end = int(parts[2])
        name = parts[3] if len(parts) > 3 else f"region_{start}_{end}"
        return (chrom, start, end, name)
    except (ValueError, IndexError):
        return None


_BED_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def split_bed_spec(spec: str) -> tuple[str | None, str]:
    """Split ``NAME=path`` into (NAME, path); plain paths (even with '=') give (None, path)."""
    name, sep, path = spec.partition("=")
    if sep and _BED_NAME_RE.match(name) and not os.path.exists(spec):
        return name, path
    return None, spec


def _iter_bed_intervals(file_path: str) -> Iterator[tuple[str, int, int, str]]:
    """Yield valid (chrom, start, end, name) from a plain or gzip BED; raise if unreadable."""
    import gzip

    opener = gzip.open if file_path.endswith(".gz") else open
    with opener(file_path, "rt", encoding="utf-8") as f:  # type: ignore[operator]
        for line_num, line in enumerate(f, 1):
            parsed = _parse_bed_line(line)
            if parsed is None:
                continue
            if parsed[1] >= parsed[2]:
                logger.debug(f"Invalid interval in {file_path}:{line_num}: {parsed[1]}-{parsed[2]}")
                continue
            yield parsed


def _interval_tree_cls() -> Any:
    import variantcentrifuge.annotator as ann_mod

    if not getattr(ann_mod, "INTERVALTREE_AVAILABLE", False):
        logger.warning("intervaltree not available. BED annotation disabled.")
        return None
    cls = getattr(ann_mod, "IntervalTree", None)
    if cls is None:
        logger.warning("IntervalTree class not available.")
    return cls


def _load_bed_files(bed_files: list[str]) -> dict[str, Any]:
    """Load plain (unnamed) BED files into per-chromosome interval trees.

    ``NAME=path`` entries are handled by :func:`_load_named_bed_files` and skipped here.
    Unreadable files raise.
    """
    tree_cls = _interval_tree_cls()
    if tree_cls is None:
        return {}

    regions_by_chrom: dict[str, Any] = {}
    for spec in bed_files:
        name, file_path = split_bed_spec(spec)
        if name is not None:
            continue
        file_name = _sanitize_name(file_path)
        n = 0
        for chrom, start, end, label in _iter_bed_intervals(file_path):
            if chrom not in regions_by_chrom:
                regions_by_chrom[chrom] = tree_cls()
            regions_by_chrom[chrom][start:end] = f"{label}_{file_name}"
            n += 1
        logger.info(f"Loaded {n} regions from BED file: {file_path}")
    return regions_by_chrom


def _load_named_bed_files(bed_files: list[str]) -> dict[str, dict[str, Any]]:
    """Load ``NAME=path`` BED entries into {NAME: {chrom: IntervalTree}}. Unreadable files raise."""
    tree_cls = _interval_tree_cls()
    named: dict[str, dict[str, Any]] = {}
    if tree_cls is None:
        return named
    for spec in bed_files:
        name, file_path = split_bed_spec(spec)
        if name is None:
            continue
        trees = named.setdefault(name, {})
        for chrom, start, end, _ in _iter_bed_intervals(file_path):
            if chrom not in trees:
                trees[chrom] = tree_cls()
            trees[chrom][start:end] = True
        logger.info(f"Loaded named BED {name} from {file_path}")
    return named


def _load_gene_lists(gene_list_files: list[str]) -> dict[str, set[str]]:
    """Load gene lists from text files.

    Parameters
    ----------
    gene_list_files : list[str]
        List of gene list file paths

    Returns
    -------
    dict[str, set[str]]
        Dictionary mapping list names to sets of gene symbols (uppercase)
    """
    gene_lists: dict[str, set[str]] = {}

    for file_path in gene_list_files:
        if not os.path.exists(file_path):
            logger.warning(f"Gene list file not found: {file_path}")
            continue

        list_name = _sanitize_name(file_path)
        genes_loaded = 0

        try:
            with open(file_path, encoding="utf-8") as f:
                genes = set()
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        for gene in line.split(","):
                            gene = gene.strip()
                            if gene:
                                genes.add(gene.upper())
                                genes_loaded += 1

                gene_lists[list_name] = genes

            logger.info(f"Loaded {genes_loaded} genes from list: {file_path}")

        except Exception as e:
            logger.error(f"Failed to load gene list {file_path}: {e}")

    return gene_lists


def _load_json_gene_data(json_files: list[str], mapping_config: str) -> dict[str, dict[str, Any]]:
    """Load structured gene data from JSON files.

    Parameters
    ----------
    json_files : list[str]
        List of JSON file paths
    mapping_config : str
        JSON string specifying field mapping

    Returns
    -------
    dict[str, dict[str, Any]]
        Dictionary mapping gene symbols to their data dictionaries
    """
    if not mapping_config:
        logger.warning("No JSON gene mapping provided. Skipping JSON gene data.")
        return {}

    try:
        mapping = json.loads(mapping_config)
        identifier_key = mapping.get("identifier", "gene_symbol")
        data_fields = mapping.get("dataFields", [])

        if not data_fields:
            logger.warning("No data fields specified in JSON mapping. Skipping JSON gene data.")
            return {}

    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON mapping configuration: {e}")
        return {}

    json_gene_data: dict[str, dict[str, Any]] = {}

    for file_path in json_files:
        if not os.path.exists(file_path):
            logger.warning(f"JSON gene file not found: {file_path}")
            continue

        genes_loaded = 0

        try:
            with open(file_path, encoding="utf-8") as f:
                data = json.load(f)

            items = (
                data if isinstance(data, list) else data.values() if isinstance(data, dict) else []
            )

            for item in items:
                if not isinstance(item, dict):
                    continue

                identifier = item.get(identifier_key)
                if not identifier:
                    continue

                gene_symbol = str(identifier).upper()
                gene_annotations = {}

                for field in data_fields:
                    value = item.get(field)
                    if value is not None:
                        str_value = str(value).strip()
                        if str_value:
                            clean_field = re.sub(r"[^a-zA-Z0-9_]", "_", field)
                            gene_annotations[clean_field] = str_value

                if gene_annotations:
                    json_gene_data[gene_symbol] = gene_annotations
                    genes_loaded += 1

            logger.info(f"Loaded data for {genes_loaded} genes from JSON: {file_path}")

        except Exception as e:
            logger.error(f"Failed to load JSON gene file {file_path}: {e}")

    return json_gene_data
