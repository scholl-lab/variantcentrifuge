# File: variantcentrifuge/converter.py
# Location: variantcentrifuge/variantcentrifuge/converter.py

"""
File conversion module.

This module provides functionality to convert TSV files to XLSX format
and append additional sheets in a single writer session or standalone.
It also supports producing JSON files for the HTML report.
"""

import json
import logging
import os
import re
from typing import Any, cast

import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

logger = logging.getLogger("variantcentrifuge")

# Compile GT parsing pattern once at module level for performance
GT_PATTERN = re.compile(r"([^()]+)\(([^)]+)\)")


def _load_igv_lookup(base_path: str) -> dict[tuple[Any, ...], str]:
    """Find and load IGV reports map, returning lookup dictionary for variants."""
    output_dir = os.path.dirname(base_path)
    igv_map_path = os.path.join(output_dir, "report", "igv", "igv_reports_map.json")
    if not os.path.exists(igv_map_path):
        parent_dir = os.path.dirname(os.path.dirname(output_dir))
        alt_path = os.path.join(parent_dir, "report", "igv", "igv_reports_map.json")
        if os.path.exists(alt_path):
            igv_map_path = alt_path

    if not os.path.exists(igv_map_path):
        logger.warning(f"IGV reports map file not found at: {igv_map_path}")
        return {}

    logger.info(f"Found IGV reports mapping file: {igv_map_path}")
    try:
        with open(igv_map_path, encoding="utf-8") as f:
            igv_map_data = json.load(f)

        igv_map = igv_map_data.get("variants", igv_map_data)
        igv_lookup: dict[tuple[Any, ...], str] = {}
        for entry in igv_map:
            chrom = str(entry.get("chrom", ""))
            pos = str(entry.get("pos", ""))
            ref = str(entry.get("ref", ""))
            alt = str(entry.get("alt", ""))

            if "sample_reports" in entry:
                for sample_id, report_path in entry["sample_reports"].items():
                    key = (chrom, pos, ref, alt, sample_id)
                    igv_lookup[key] = report_path
            elif "sample_id" in entry and "report_path" in entry:
                key = (chrom, pos, ref, alt, entry["sample_id"])
                igv_lookup[key] = entry["report_path"]

        logger.info(f"Created lookup for {len(igv_lookup)} IGV report entries")
        return igv_lookup
    except Exception as e:
        logger.error(f"Error loading IGV reports map: {e}")
        return {}


def _get_igv_reports_for_variant(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    gt_value: str,
    igv_lookup: dict[tuple[Any, ...], str],
) -> list[tuple[str, str]]:
    """Return list of (sample_id, report_path) for non-reference samples in gt_value."""
    if not gt_value:
        return []
    igv_reports = []
    for entry in gt_value.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        m = GT_PATTERN.match(entry)
        if m:
            sample_id = m.group(1).strip()
            genotype = m.group(2).strip()
            if genotype in ["0/0", "./."]:
                continue
            lookup_key = (chrom, pos, ref, alt, sample_id)
            if lookup_key in igv_lookup:
                igv_reports.append((sample_id, igv_lookup[lookup_key]))
    return igv_reports


def _detect_url_columns(df: pd.DataFrame) -> list[str]:
    """Identify columns where >= 70% of sample non-empty values start with http:// or https://."""
    url_columns = []
    for col in df.columns:
        if len(df) == 0:
            continue
        sample_values = [str(val).strip() for val in df[col].dropna().iloc[:4] if str(val).strip()]
        if sample_values:
            url_count = sum(1 for val in sample_values if val.startswith(("http://", "https://")))
            if url_count >= len(sample_values) * 0.7:
                url_columns.append(col)
    return url_columns


def convert_to_excel(
    tsv_file: str,
    cfg: dict[str, Any],
    df: pd.DataFrame | None = None,
    extra_sheets: dict[str, pd.DataFrame] | None = None,
) -> str:
    """
    Convert a TSV file or DataFrame to XLSX format with Results and optional extra sheets.

    Writes all sheets and formatting (freeze panes, autofilters, hyperlinks, IGV links)
    in a single XlsxWriter session without repeated workbook reloads.

    Parameters
    ----------
    tsv_file : str
        Path to the TSV file (or target .xlsx file).
    cfg : dict
        Configuration dictionary containing various settings.
    df : pd.DataFrame, optional
        Optional DataFrame to use instead of reading from disk.
    extra_sheets : dict of str to pd.DataFrame, optional
        Additional sheets to write in the same session (e.g. Metadata, Statistics, Gene Burden).

    Returns
    -------
    str
        The path to the generated XLSX file.
    """
    if df is not None:
        logger.debug("Using provided DataFrame for Excel conversion (skipping disk read)")
    else:
        genomic_dtypes = {
            "CHROM": "str",
            "POS": "str",
            "REF": "str",
            "ALT": "str",
            "QUAL": "str",
            "FILTER": "str",
        }
        df = pd.read_csv(
            tsv_file,
            sep="\t",
            na_values="NA",
            dtype=cast(Any, genomic_dtypes),
            low_memory=False,
        )

    if len(df) == 0:
        logger.warning("TSV file contains only headers or is empty. Creating header-only Excel.")

    if "igv_links" in df.columns:
        logger.debug("Removing raw igv_links column from Excel conversion")
        df = df.drop(columns=["igv_links"])

    xlsx_file = tsv_file if tsv_file.endswith(".xlsx") else os.path.splitext(tsv_file)[0] + ".xlsx"

    igv_enabled = cfg.get("igv_enabled", False)
    igv_lookup = _load_igv_lookup(xlsx_file) if igv_enabled else {}

    with pd.ExcelWriter(xlsx_file, engine="xlsxwriter") as writer:
        df.to_excel(writer, sheet_name="Results", index=False)
        ws = writer.sheets["Results"]
        book = writer.book
        hyperlink_format = book.add_format({"color": "blue", "underline": 1})

        ws.freeze_panes(1, 0)

        # Detect and format URL columns
        url_columns = _detect_url_columns(df)
        if url_columns:
            logger.info(f"Converting {len(url_columns)} URL columns to clickable hyperlinks")
            for col_name in url_columns:
                col_idx = df.columns.get_loc(col_name)
                for row_idx, val in enumerate(df[col_name]):
                    s_val = str(val).strip() if pd.notna(val) else ""
                    if s_val.startswith(("http://", "https://")):
                        ws.write_url(row_idx + 1, col_idx, s_val, hyperlink_format, string=col_name)

        # Handle IGV Report Links column
        num_cols = len(df.columns)
        if igv_enabled and igv_lookup:
            igv_col_idx = num_cols
            ws.write_string(0, igv_col_idx, "IGV Report Links")
            for row_idx in range(len(df)):
                chrom = (
                    str(df["CHROM"].iloc[row_idx])
                    if "CHROM" in df.columns and pd.notna(df["CHROM"].iloc[row_idx])
                    else ""
                )
                pos = (
                    str(df["POS"].iloc[row_idx])
                    if "POS" in df.columns and pd.notna(df["POS"].iloc[row_idx])
                    else ""
                )
                ref = (
                    str(df["REF"].iloc[row_idx])
                    if "REF" in df.columns and pd.notna(df["REF"].iloc[row_idx])
                    else ""
                )
                alt = (
                    str(df["ALT"].iloc[row_idx])
                    if "ALT" in df.columns and pd.notna(df["ALT"].iloc[row_idx])
                    else ""
                )
                gt_val = (
                    str(df["GT"].iloc[row_idx])
                    if "GT" in df.columns and pd.notna(df["GT"].iloc[row_idx])
                    else ""
                )

                igv_reports = _get_igv_reports_for_variant(chrom, pos, ref, alt, gt_val, igv_lookup)
                if len(igv_reports) == 1:
                    sample_id, rpath = igv_reports[0]
                    rel_link = os.path.join("report", rpath)
                    ws.write_url(
                        row_idx + 1,
                        igv_col_idx,
                        f"external:{rel_link}",
                        hyperlink_format,
                        string=sample_id,
                    )
                elif len(igv_reports) > 1:
                    first_sample, _ = igv_reports[0]
                    num_others = len(igv_reports) - 1
                    ws.write_string(
                        row_idx + 1, igv_col_idx, f"{first_sample} (+{num_others} others)"
                    )
                    link_texts = [
                        f"{sid} ({os.path.join('report', rp)})" for sid, rp in igv_reports
                    ]
                    ws.write_comment(
                        row_idx + 1, igv_col_idx, "; ".join(link_texts), {"author": "IGV Reports"}
                    )
                else:
                    ws.write_string(row_idx + 1, igv_col_idx, "N/A")

            # Autofilter across all columns including IGV column
            ws.autofilter(0, 0, len(df), igv_col_idx)
        else:
            if num_cols > 0:
                ws.autofilter(0, 0, len(df), num_cols - 1)

        # Write extra sheets
        if extra_sheets:
            for sheet_name, sheet_df in extra_sheets.items():
                sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)
                extra_ws = writer.sheets[sheet_name]
                extra_ws.freeze_panes(1, 0)
                if len(sheet_df.columns) > 0:
                    extra_ws.autofilter(0, 0, len(sheet_df), len(sheet_df.columns) - 1)

    return xlsx_file


def append_tsv_as_sheet(xlsx_file: str, tsv_file: str, sheet_name: str = "Metadata") -> None:
    """
    Append a TSV file as a new sheet to an existing XLSX file (standalone API).

    Parameters
    ----------
    xlsx_file : str
        Path to the existing XLSX file.
    tsv_file : str
        Path to the TSV file to append as a sheet.
    sheet_name : str, optional
        Name of the new sheet to be appended. Defaults to "Metadata".
    """
    df = pd.read_csv(tsv_file, sep="\t", header=0, low_memory=False)
    with pd.ExcelWriter(xlsx_file, engine="openpyxl", mode="a", if_sheet_exists="new") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)


def finalize_excel_file(xlsx_file: str, cfg: dict[str, Any]) -> None:
    """
    Apply final formatting to all sheets in xlsx_file using openpyxl (standalone API).

    - Freeze top row
    - Enable auto-filter
    - Hyperlinks on URL columns
    - Add IGV Report Links column if enabled and not already present
    """
    wb = load_workbook(xlsx_file)
    igv_enabled = cfg.get("igv_enabled", False)
    igv_lookup = _load_igv_lookup(xlsx_file) if igv_enabled else {}

    for ws in wb.worksheets:
        ws.freeze_panes = "A2"
        max_col_letter = get_column_letter(ws.max_column)
        ws.auto_filter.ref = f"A1:{max_col_letter}1"

        if ws.title != "Results":
            continue

        header_row = [cell.value for cell in next(ws.iter_rows(min_row=1, max_row=1))]
        header_indices = {
            col_name: idx
            for idx, col_name in enumerate(header_row, 1)
            if col_name in ["CHROM", "POS", "REF", "ALT", "GT"]
        }

        # URL columns
        url_columns = []
        if ws.max_row > 1:
            for idx, _col_name in enumerate(header_row, 1):
                sample_cells = [
                    ws.cell(row=r, column=idx) for r in range(2, min(6, ws.max_row + 1))
                ]
                has_hyperlink = any(c.hyperlink is not None for c in sample_cells)
                sample_values = [str(c.value).strip() for c in sample_cells if c.value]
                url_count = sum(
                    1 for val in sample_values if val.startswith(("http://", "https://"))
                )
                if has_hyperlink or (sample_values and url_count >= len(sample_values) * 0.7):
                    url_columns.append(idx)

        if url_columns:
            for row_idx in range(2, ws.max_row + 1):
                for col_idx in url_columns:
                    cell = ws.cell(row=row_idx, column=col_idx)
                    val = str(cell.value or "").strip()
                    if val.startswith(("http://", "https://")):
                        col_name = str(header_row[col_idx - 1])
                        cell.hyperlink = val
                        cell.value = col_name
                    if cell.hyperlink:
                        cell.style = "Hyperlink"

        # IGV links
        if igv_enabled and igv_lookup and "IGV Report Links" not in header_row:
            igv_col = ws.max_column + 1
            ws.cell(row=1, column=igv_col, value="IGV Report Links")
            from openpyxl.comments import Comment

            for row_idx in range(2, ws.max_row + 1):
                chrom = str(ws.cell(row=row_idx, column=header_indices.get("CHROM", 0)).value or "")
                pos = str(ws.cell(row=row_idx, column=header_indices.get("POS", 0)).value or "")
                ref = str(ws.cell(row=row_idx, column=header_indices.get("REF", 0)).value or "")
                alt = str(ws.cell(row=row_idx, column=header_indices.get("ALT", 0)).value or "")
                gt_val = (
                    str(ws.cell(row=row_idx, column=header_indices.get("GT", 0)).value or "")
                    if "GT" in header_indices
                    else ""
                )

                igv_reports = _get_igv_reports_for_variant(chrom, pos, ref, alt, gt_val, igv_lookup)
                igv_cell = ws.cell(row=row_idx, column=igv_col)
                if len(igv_reports) == 1:
                    sample_id, rpath = igv_reports[0]
                    igv_cell.value = sample_id
                    igv_cell.hyperlink = os.path.join("report", rpath)
                    igv_cell.style = "Hyperlink"
                elif len(igv_reports) > 1:
                    first_sample, _ = igv_reports[0]
                    num_others = len(igv_reports) - 1
                    igv_cell.value = f"{first_sample} (+{num_others} others)"
                    link_texts = [
                        f"{sid} ({os.path.join('report', rp)})" for sid, rp in igv_reports
                    ]
                    igv_cell.comment = Comment("; ".join(link_texts), "IGV Reports")
                else:
                    igv_cell.value = "N/A"

    wb.save(xlsx_file)


def produce_report_json(variant_tsv: str, output_dir: str) -> None:
    """
    Produce JSON files (variants.json and summary.json) from a TSV of variants.

    Parameters
    ----------
    variant_tsv : str
        Path to the TSV file containing annotated variants.
    output_dir : str
        Directory to write the JSON outputs (variants.json, summary.json).
    """
    logger.info(f"Producing JSON files for HTML report from {variant_tsv}")

    try:
        df = pd.read_csv(variant_tsv, sep="\t", low_memory=False)
    except pd.errors.EmptyDataError:
        logger.warning("Empty TSV file provided. Generating empty JSON files.")
        df = pd.DataFrame()

    variants = df.to_dict(orient="records") if not df.empty else []

    for variant_dict in variants:
        for key, val in variant_dict.items():
            if pd.isna(val):
                variant_dict[key] = None

    igv_lookup = _load_igv_lookup(os.path.join(output_dir, "variants.tsv"))
    if igv_lookup:
        for variant_dict in variants:
            variant_dict["igv_links"] = []
            chrom = str(variant_dict.get("CHROM", ""))
            pos = str(variant_dict.get("POS", ""))
            ref = str(variant_dict.get("REF", ""))
            alt = str(variant_dict.get("ALT", ""))
            gt_value = variant_dict.get("GT", "")

            igv_reports = _get_igv_reports_for_variant(chrom, pos, ref, alt, gt_value, igv_lookup)
            for sample_id, rpath in igv_reports:
                variant_dict["igv_links"].append({"sample_id": sample_id, "report_path": rpath})
        logger.info("Enriched variants with IGV report links")

    num_variants = len(df) if not df.empty else 0
    unique_genes: Any = df["GENE"].unique() if not df.empty and "GENE" in df.columns else []
    num_genes = len(unique_genes)

    impact_counts: dict[str, int] = {}
    if not df.empty and "IMPACT" in df.columns:
        impact_counts = df["IMPACT"].value_counts().to_dict()

    num_samples = 0
    if not df.empty and "GT" in df.columns:
        sample_ids = set()
        for gt_value in df["GT"]:
            if pd.notna(gt_value):
                sample_entries = str(gt_value).split(";")
                for entry in sample_entries:
                    entry = entry.strip()
                    if entry:
                        m = GT_PATTERN.match(entry)
                        if m:
                            sample_ids.add(m.group(1).strip())
        num_samples = len(sample_ids)

    inheritance_distribution: dict[str, int] = {}
    if not df.empty and "Inheritance_Pattern" in df.columns:
        valid_patterns = df[~df["Inheritance_Pattern"].isin(["none", "reference"])][
            "Inheritance_Pattern"
        ]
        if len(valid_patterns) > 0:
            inheritance_distribution = valid_patterns.value_counts().to_dict()

    top_genes: list[dict[str, Any]] = []
    if not df.empty and "GENE" in df.columns:
        gene_counts = df["GENE"].dropna().value_counts().head(10)
        top_genes = [{"gene": gene, "count": int(count)} for gene, count in gene_counts.items()]

    summary_data = {
        "num_variants": num_variants,
        "num_genes": num_genes,
        "num_samples": num_samples,
        "impact_distribution": impact_counts,
        "inheritance_distribution": inheritance_distribution,
        "top_genes": top_genes,
    }

    report_dir = os.path.join(output_dir, "report")
    os.makedirs(report_dir, exist_ok=True)
    variants_json_path = os.path.join(report_dir, "variants.json")
    with open(variants_json_path, "w", encoding="utf-8") as vjf:
        json.dump(variants, vjf, indent=2)

    summary_json_path = os.path.join(report_dir, "summary.json")
    with open(summary_json_path, "w", encoding="utf-8") as sjf:
        json.dump(summary_data, sjf, indent=2)

    logger.info(f"JSON files produced: {variants_json_path}, {summary_json_path}")
