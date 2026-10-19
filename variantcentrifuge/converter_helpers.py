"""Helper functions for Excel conversion and report enrichment."""

from __future__ import annotations

import json
import logging
import os
import re
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger("variantcentrifuge")

GT_PATTERN = re.compile(r"([^()]+)\(([^)]+)\)")
URL_PREFIXES = (
    "http://",
    "https://",
    "ftp://",
    "ftps://",
    "mailto:",
    "internal:",
    "external:",
    "file://",
)


def fmt_coord(val: Any) -> str:
    """Normalize coordinates, converting float integral values (e.g. 100.0) to '100'."""
    if val is None or pd.isna(val):
        return ""
    if isinstance(val, float) and val.is_integer():
        return str(int(val))
    s = str(val).strip()
    return s[:-2] if s.endswith(".0") and s[:-2].isdigit() else s


def load_igv_lookup(
    base_path: str, igv_enabled: bool = True, search_parents: bool = True
) -> dict[tuple[Any, ...], str]:
    """Find and load IGV reports map, returning lookup dictionary for variants."""
    output_dir = os.path.dirname(base_path)
    igv_map_path = os.path.join(output_dir, "report", "igv", "igv_reports_map.json")
    if not os.path.exists(igv_map_path) and search_parents:
        parent_dir = os.path.dirname(os.path.dirname(output_dir))
        alt = os.path.join(parent_dir, "report", "igv", "igv_reports_map.json")
        igv_map_path = alt if os.path.exists(alt) else igv_map_path

    if not os.path.exists(igv_map_path):
        if igv_enabled:
            logger.warning(f"IGV reports map file not found at: {igv_map_path}")
        return {}

    logger.info(f"Found IGV reports mapping file: {igv_map_path}")
    try:
        with open(igv_map_path, encoding="utf-8") as f:
            igv_map_data = json.load(f)

        igv_map = (
            igv_map_data.get("variants", [])
            if isinstance(igv_map_data, dict)
            else (igv_map_data if isinstance(igv_map_data, list) else [])
        )

        igv_lookup: dict[tuple[Any, ...], str] = {}
        for entry in igv_map:
            c = tuple(fmt_coord(entry.get(k, "")) for k in ("chrom", "pos", "ref", "alt"))
            if reports := (entry.get("sample_reports") or entry.get("sample_links")):
                for sid, rpath in reports.items():
                    igv_lookup[(c[0], c[1], c[2], c[3], sid)] = rpath
            elif "sample_id" in entry and "report_path" in entry:
                igv_lookup[(c[0], c[1], c[2], c[3], entry["sample_id"])] = entry["report_path"]

        logger.info(f"Created lookup for {len(igv_lookup)} IGV report entries")
        return igv_lookup
    except Exception as e:
        logger.error(f"Error loading IGV reports map: {e}")
        return {}


def get_igv_reports_for_variant(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    gt_value: Any,
    igv_lookup: dict[tuple[Any, ...], str],
) -> list[tuple[str, str]]:
    """Return list of (sample_id, report_path) for non-reference samples in gt_value."""
    gt_str = str(gt_value).strip() if gt_value is not None and pd.notna(gt_value) else ""
    if not gt_str:
        return []
    chrom, pos, ref, alt = fmt_coord(chrom), fmt_coord(pos), fmt_coord(ref), fmt_coord(alt)
    igv_reports = []
    for entry in gt_str.split(";"):
        if (m := GT_PATTERN.match(entry.strip())) and m.group(2).strip() not in {
            "0/0",
            "./.",
            "0",
            ".",
            "0|0",
        }:
            sid = m.group(1).strip()
            if (key := (chrom, pos, ref, alt, sid)) in igv_lookup:
                igv_reports.append((sid, igv_lookup[key]))
    return igv_reports


def excel_cell_value(val: Any) -> Any:
    """Normalize values for writing to Excel cells."""
    if val is None:
        return None
    if isinstance(val, (list, tuple, dict, set, np.ndarray)):
        return str(val)
    return None if pd.isna(val) else val


def detect_url_columns(df: pd.DataFrame) -> set[int]:
    """Identify column indices where >= 70% of sample non-empty values start with http:// or https://."""
    if len(df) == 0:
        return set()
    url_columns: set[int] = set()
    for col_idx, col_name in enumerate(df.columns):
        if not col_name:
            continue
        sample_vals = [
            s
            for raw in df.iloc[:4, col_idx]
            if (c := excel_cell_value(raw)) and (s := str(c).strip())
        ]
        if (
            sample_vals
            and sum(v.startswith(("http://", "https://")) for v in sample_vals)
            >= len(sample_vals) * 0.7
        ):
            url_columns.add(col_idx)
    return url_columns
