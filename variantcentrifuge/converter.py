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
import shutil
from typing import Any, cast

import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

from .checkpoint_models import AtomicFileOperation
from .converter_helpers import (
    GT_PATTERN,
    URL_PREFIXES,
)
from .converter_helpers import (
    detect_url_columns as _detect_url_columns,
)
from .converter_helpers import (
    fmt_coord as _fmt_coord,
)
from .converter_helpers import (
    get_igv_reports_for_variant as _get_igv_reports_for_variant,
)
from .converter_helpers import (
    load_igv_lookup as _load_igv_lookup,
)

logger = logging.getLogger("variantcentrifuge")


def convert_to_excel(
    tsv_file: str,
    cfg: dict[str, Any],
    df: pd.DataFrame | None = None,
    extra_sheets: dict[str, pd.DataFrame] | None = None,
    output_file: str | None = None,
) -> str:
    """Convert a TSV file or DataFrame to XLSX format with Results and optional extra sheets."""
    if extra_sheets:
        if conflicts := [s for s in extra_sheets if s.lower() in {"results", "igv_overflow"}]:
            raise ValueError(f"extra_sheets cannot overwrite reserved sheet(s): {conflicts}")
        if len({s.lower() for s in extra_sheets}) != len(extra_sheets):
            raise ValueError("extra_sheets contains duplicate sheet names")

    if df is not None:
        logger.debug("Using provided DataFrame for Excel conversion (skipping disk read)")
    else:
        genomic_dtypes = dict.fromkeys(("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER"), "str")
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

    tgt = output_file if output_file is not None else tsv_file
    xlsx_file = tgt if tgt.endswith(".xlsx") else f"{os.path.splitext(tgt)[0]}.xlsx"

    igv_enabled = cfg.get("igv_enabled", False)
    igv_lookup = _load_igv_lookup(xlsx_file, igv_enabled=igv_enabled) if igv_enabled else {}

    with (
        AtomicFileOperation(xlsx_file) as temp_xlsx,
        pd.ExcelWriter(
            str(temp_xlsx),
            engine="xlsxwriter",
            engine_kwargs={"options": {"strings_to_urls": False}},
        ) as writer,
    ):
        df.to_excel(writer, sheet_name="Results", index=False)
        ws = writer.sheets["Results"]
        book = writer.book
        hyperlink_format = book.add_format({"color": "blue", "underline": 1})

        ws.freeze_panes(1, 0)

        # Format detected URL columns and individual hyperlinks (single registration)
        url_col_indices = _detect_url_columns(df)
        for col_idx, col_name in enumerate(df.columns):
            is_url_col = col_idx in url_col_indices
            col_series = df.iloc[:, col_idx]
            if not is_url_col and (
                pd.api.types.is_numeric_dtype(col_series) or pd.api.types.is_bool_dtype(col_series)
            ):
                continue
            str_col_name = str(col_name)
            for row_idx, val in enumerate(col_series):
                s_val = val.strip() if isinstance(val, str) else ""
                if not s_val:
                    continue
                if is_url_col and s_val.startswith(("http://", "https://")):
                    if (
                        ws.write_url(
                            row_idx + 1, col_idx, s_val, hyperlink_format, string=str_col_name
                        )
                        != 0
                    ):
                        ws.write_string(row_idx + 1, col_idx, s_val)
                elif (
                    s_val.startswith(URL_PREFIXES)
                    and ws.write_url(row_idx + 1, col_idx, s_val) != 0
                ):
                    ws.write_string(row_idx + 1, col_idx, s_val)

        # Handle IGV Report Links column
        num_cols = len(df.columns)
        has_coords = all(c in df.columns for c in ("CHROM", "POS", "REF", "ALT"))
        overflow_rows: list[list[Any]] = []
        if igv_enabled and igv_lookup and has_coords:
            igv_col_idx = num_cols
            ws.write_string(0, igv_col_idx, "IGV Report Links")
            for row_idx in range(len(df)):
                raw_pos, chrom = (
                    df["POS"].iloc[row_idx],
                    _fmt_coord(df["CHROM"].iloc[row_idx]),
                )
                pos, ref, alt = (
                    _fmt_coord(raw_pos),
                    _fmt_coord(df["REF"].iloc[row_idx]),
                    _fmt_coord(df["ALT"].iloc[row_idx]),
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
                    ws.write_comment(
                        row_idx + 1, igv_col_idx, f"Full path: {rel_link}", {"visible": False}
                    )
                    if (
                        ws.write_url(
                            row_idx + 1,
                            igv_col_idx,
                            f"external:{rel_link}",
                            hyperlink_format,
                            string=sample_id,
                        )
                        != 0
                    ):
                        ws.write_string(row_idx + 1, igv_col_idx, rel_link)
                        overflow_rows.append(
                            [row_idx + 2, chrom, raw_pos, ref, alt, sample_id, rel_link]
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
                    comment_text = "; ".join(link_texts)
                    if len(comment_text) > 32000:
                        suffix = "\n... (see IGV_Overflow sheet for all sample links)"
                        truncated = comment_text[: 32000 - len(suffix)] + suffix
                        ws.write_comment(
                            row_idx + 1, igv_col_idx, truncated, {"author": "IGV Reports"}
                        )
                        for sid, rp in igv_reports:
                            overflow_rows.append(
                                [
                                    row_idx + 2,
                                    chrom,
                                    raw_pos,
                                    ref,
                                    alt,
                                    sid,
                                    os.path.join("report", rp),
                                ]
                            )
                    else:
                        ws.write_comment(
                            row_idx + 1, igv_col_idx, comment_text, {"author": "IGV Reports"}
                        )
                else:
                    ws.write_string(row_idx + 1, igv_col_idx, "N/A")

            # Autofilter across all columns including IGV column
            ws.autofilter(0, 0, len(df), igv_col_idx)
        else:
            if num_cols > 0:
                ws.autofilter(0, 0, len(df), num_cols - 1)

        # Write IGV overflow sheet if any comments exceeded limit
        if overflow_rows:
            overflow_df = pd.DataFrame(
                overflow_rows,
                columns=["Row", "CHROM", "POS", "REF", "ALT", "Sample_ID", "IGV_Report_Path"],
            )
            overflow_df.to_excel(writer, sheet_name="IGV_Overflow", index=False)
            overflow_ws = writer.sheets["IGV_Overflow"]
            for ov_idx, r in enumerate(overflow_rows):
                overflow_ws.write_string(ov_idx + 1, 1, str(r[1]))
                overflow_ws.write_string(ov_idx + 1, 5, str(r[5]))
            overflow_ws.freeze_panes(1, 0)
            overflow_ws.autofilter(0, 0, len(overflow_df), len(overflow_df.columns) - 1)

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
    with AtomicFileOperation(xlsx_file) as temp_xlsx:
        shutil.copy2(xlsx_file, temp_xlsx)
        with pd.ExcelWriter(
            str(temp_xlsx), engine="openpyxl", mode="a", if_sheet_exists="new"
        ) as writer:
            df.to_excel(writer, index=False, sheet_name=sheet_name)


def finalize_excel_file(xlsx_file: str, cfg: dict[str, Any]) -> None:
    """
    Apply final formatting to all sheets in xlsx_file using openpyxl (standalone API).

    - Freeze top row
    - Enable auto-filter
    - Hyperlinks on URL columns
    - Add IGV Report Links column if enabled and not already present
    """
    with AtomicFileOperation(xlsx_file) as temp_xlsx:
        shutil.copy2(xlsx_file, temp_xlsx)
        wb = load_workbook(str(temp_xlsx))
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
                    sample_values = [
                        (
                            c.hyperlink.target
                            if c.hyperlink and c.hyperlink.target
                            else str(c.value or "")
                        ).strip()
                        for c in sample_cells
                        if c.value or (c.hyperlink and c.hyperlink.target)
                    ]
                    url_count = sum(
                        1 for val in sample_values if val.startswith(("http://", "https://"))
                    )
                    if sample_values and url_count >= len(sample_values) * 0.7:
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
            has_coords = all(c in header_indices for c in ("CHROM", "POS", "REF", "ALT"))
            overflow_rows: list[list[Any]] = []
            if igv_enabled and igv_lookup and has_coords and "IGV Report Links" not in header_row:
                igv_col = ws.max_column + 1
                ws.cell(row=1, column=igv_col, value="IGV Report Links")
                from openpyxl.comments import Comment

                for row_idx in range(2, ws.max_row + 1):
                    raw_chrom = ws.cell(row=row_idx, column=header_indices["CHROM"]).value
                    raw_pos = ws.cell(row=row_idx, column=header_indices["POS"]).value
                    raw_ref = ws.cell(row=row_idx, column=header_indices["REF"]).value
                    raw_alt = ws.cell(row=row_idx, column=header_indices["ALT"]).value
                    chrom = _fmt_coord(raw_chrom)
                    pos = _fmt_coord(raw_pos)
                    ref = _fmt_coord(raw_ref)
                    alt = _fmt_coord(raw_alt)
                    gt_val = (
                        str(ws.cell(row=row_idx, column=header_indices["GT"]).value or "")
                        if "GT" in header_indices
                        else ""
                    )

                    igv_reports = _get_igv_reports_for_variant(
                        chrom, pos, ref, alt, gt_val, igv_lookup
                    )
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
                        comment_text = "; ".join(link_texts)
                        if len(comment_text) > 32000:
                            suffix = "\n... (see IGV_Overflow sheet for all sample links)"
                            comment_text = comment_text[: 32000 - len(suffix)] + suffix
                            for sid, rp in igv_reports:
                                overflow_rows.append(
                                    [
                                        row_idx,
                                        raw_chrom,
                                        raw_pos,
                                        raw_ref,
                                        raw_alt,
                                        sid,
                                        os.path.join("report", rp),
                                    ]
                                )
                        igv_cell.comment = Comment(comment_text, "IGV Reports")
                    else:
                        igv_cell.value = "N/A"

                if overflow_rows:
                    overflow_headers = [
                        "Row",
                        "CHROM",
                        "POS",
                        "REF",
                        "ALT",
                        "Sample_ID",
                        "IGV_Report_Path",
                    ]
                    if "IGV_Overflow" in wb.sheetnames:
                        overflow_ws = wb["IGV_Overflow"]
                    else:
                        overflow_ws = wb.create_sheet(title="IGV_Overflow")
                        overflow_ws.append(overflow_headers)
                    for r in overflow_rows:
                        overflow_ws.append(r)
                        overflow_ws.cell(row=overflow_ws.max_row, column=2).data_type = "s"
                        overflow_ws.cell(row=overflow_ws.max_row, column=6).data_type = "s"
                    overflow_ws.freeze_panes = "A2"
                    max_letter = get_column_letter(len(overflow_headers))
                    overflow_ws.auto_filter.ref = f"A1:{max_letter}1"

        wb.save(str(temp_xlsx))


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

    igv_lookup = _load_igv_lookup(
        os.path.join(output_dir, "variants.tsv"), igv_enabled=False, search_parents=False
    )
    if igv_lookup:
        for variant_dict in variants:
            variant_dict["igv_links"] = []
            try:
                chrom, pos = (
                    _fmt_coord(variant_dict.get("CHROM")),
                    _fmt_coord(variant_dict.get("POS")),
                )
                ref, alt = _fmt_coord(variant_dict.get("REF")), _fmt_coord(variant_dict.get("ALT"))
                raw_gt = variant_dict.get("GT")
                gt_value = str(raw_gt) if raw_gt is not None and pd.notna(raw_gt) else ""

                igv_reports = _get_igv_reports_for_variant(
                    chrom, pos, ref, alt, gt_value, igv_lookup
                )
                for sample_id, rpath in igv_reports:
                    variant_dict["igv_links"].append({"sample_id": sample_id, "report_path": rpath})
            except Exception as e:
                logger.warning(f"Error enriching variant with IGV links: {e}")
        logger.info("Enriched variants with IGV report links")

    num_variants = len(df) if not df.empty else 0
    unique_genes: Any = df["GENE"].unique() if not df.empty and "GENE" in df.columns else []
    num_genes = len(unique_genes)

    impact_counts: dict[str, int] = {}
    if not df.empty and "IMPACT" in df.columns:
        impact_counts = df["IMPACT"].value_counts().to_dict()

    num_samples = 0
    if not df.empty and "GT" in df.columns:
        sample_ids = {
            m.group(1).strip()
            for gt in df["GT"].dropna()
            for entry in str(gt).split(";")
            if (m := GT_PATTERN.match(entry.strip()))
        }
        num_samples = len(sample_ids)

    inheritance_distribution: dict[str, int] = {}
    if not df.empty and "Inheritance_Pattern" in df.columns:
        valid_patterns = df.loc[
            ~df["Inheritance_Pattern"].isin(["none", "reference"]), "Inheritance_Pattern"
        ]
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
    with (
        AtomicFileOperation(variants_json_path) as temp_vjf,
        open(temp_vjf, "w", encoding="utf-8") as vjf,
    ):
        json.dump(variants, vjf, indent=2)

    summary_json_path = os.path.join(report_dir, "summary.json")
    with (
        AtomicFileOperation(summary_json_path) as temp_sjf,
        open(temp_sjf, "w", encoding="utf-8") as sjf,
    ):
        json.dump(summary_data, sjf, indent=2)

    logger.info(f"JSON files produced: {variants_json_path}, {summary_json_path}")
