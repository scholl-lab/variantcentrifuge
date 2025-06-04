# File: variantcentrifuge/converter.py
# Location: variantcentrifuge/variantcentrifuge/converter.py

"""
File conversion module.

This module provides functionality to convert TSV files to XLSX format
and append additional sheets. It now also supports producing JSON files
for the HTML report.
"""

import hashlib
import logging
import os
from typing import Any, Dict, List, Optional

import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

logger = logging.getLogger("variantcentrifuge")


def convert_to_excel(tsv_file: str, cfg: Dict[str, Any]) -> str:
    """
    Convert a TSV file to XLSX format with a single "Results" sheet.

    Parameters
    ----------
    tsv_file : str
        Path to the TSV file.
    cfg : dict
        Configuration dictionary containing various settings.

    Returns
    -------
    str
        The path to the generated XLSX file.
    """
    df = pd.read_csv(tsv_file, sep="\t", na_values="NA")
    xlsx_file = os.path.splitext(tsv_file)[0] + ".xlsx"
    df.to_excel(xlsx_file, index=False, sheet_name="Results")
    return xlsx_file


def append_tsv_as_sheet(xlsx_file: str, tsv_file: str, sheet_name: str = "Metadata") -> None:
    """
    Append a TSV file as a new sheet to an existing XLSX file.

    This function reads the TSV file and appends it as a new sheet in the existing
    XLSX file. Assumes the TSV file has a header row.

    Parameters
    ----------
    xlsx_file : str
        Path to the existing XLSX file.
    tsv_file : str
        Path to the TSV file to append as a sheet.
    sheet_name : str, optional
        Name of the new sheet to be appended. Defaults to "Metadata".

    Returns
    -------
    None
    """
    df = pd.read_csv(tsv_file, sep="\t", header=0)
    with pd.ExcelWriter(xlsx_file, engine="openpyxl", mode="a", if_sheet_exists="new") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)


def finalize_excel_file(xlsx_file: str, cfg: Dict[str, Any]) -> None:
    """
    Apply final formatting to all sheets in xlsx_file:
      - Freeze the top row
      - Enable auto-filter
      - Only for the 'Results' sheet, generate hyperlinks from CHROM, POS, REF, ALT
        using cfg["links"] (the link template dictionary).
    """
    wb = load_workbook(xlsx_file)
    link_templates = cfg.get("links", {})

    for ws in wb.worksheets:
        # 1. Freeze top row
        ws.freeze_panes = "A2"

        # 2. Enable auto-filter
        max_col_letter = get_column_letter(ws.max_column)
        ws.auto_filter.ref = f"A1:{max_col_letter}1"

        # Only generate links in the 'Results' sheet
        if ws.title != "Results":
            continue

        # 3. Build a map from header -> column index
        headers = [cell.value for cell in ws[1]]
        header_to_index = {str(h).strip(): i + 1 for i, h in enumerate(headers) if h}

        # Identify positions of CHROM, POS, REF, ALT
        chrom_idx = header_to_index.get("CHROM")
        pos_idx = header_to_index.get("POS")
        ref_idx = header_to_index.get("REF")
        alt_idx = header_to_index.get("ALT")

        # If missing any of these columns, skip hyperlink generation
        missing_cols = [c for c in ["CHROM", "POS", "REF", "ALT"] if c not in header_to_index]
        if missing_cols:
            logger.warning(
                f"Cannot create hyperlinks in sheet '{ws.title}': "
                f"missing columns {missing_cols}."
            )
            continue

        # For each link key in cfg['links'], if we have that column, turn it into a hyperlink
        for link_name, link_template in link_templates.items():
            link_col_idx = header_to_index.get(link_name)
            if not link_col_idx:
                # No such column in 'Results', skip it
                continue

            # For each row, build the final hyperlink from CHROM, POS, REF, ALT
            for row_num in range(2, ws.max_row + 1):
                chrom_val = ws.cell(row=row_num, column=chrom_idx).value
                pos_val = ws.cell(row=row_num, column=pos_idx).value
                ref_val = ws.cell(row=row_num, column=ref_idx).value
                alt_val = ws.cell(row=row_num, column=alt_idx).value

                # Safely convert to strings
                chrom_val = str(chrom_val) if chrom_val is not None else ""
                pos_val = str(pos_val) if pos_val is not None else ""
                ref_val = str(ref_val) if ref_val is not None else ""
                alt_val = str(alt_val) if alt_val is not None else ""

                hyperlink_url = link_template.format(
                    CHROM=chrom_val, POS=pos_val, REF=ref_val, ALT=alt_val
                )

                cell = ws.cell(row=row_num, column=link_col_idx)
                cell.hyperlink = hyperlink_url
                cell.style = "Hyperlink"
                # Optionally change display text:
                # cell.value = link_name  # e.g., 'ClinVar'

    wb.save(xlsx_file)


def produce_report_json(variant_tsv: str, output_dir: str) -> None:
    """
    Produce JSON files (variants.json and summary.json) from a TSV of variants.

    The TSV file is expected to have columns including GENE, CHROM, POS, REF, ALT, IMPACT.
    Additional columns can also be included. This function directly converts the variants
    into a JSON array and computes summary statistics for summary.json.

    Parameters
    ----------
    variant_tsv : str
        Path to the TSV file containing variant data.
    output_dir : str
        Directory to write the JSON outputs (variants.json, summary.json).

    Returns
    -------
    None
    """
    logger.info(f"Producing JSON files for HTML report from {variant_tsv}")

    df = pd.read_csv(variant_tsv, sep="\t", header=0)
    variants = df.to_dict(orient="records")

    num_variants = len(df)
    unique_genes = df["GENE"].unique() if "GENE" in df.columns else []
    num_genes = len(unique_genes)

    impact_counts = {}
    if "IMPACT" in df.columns:
        impact_counts = df["IMPACT"].value_counts().to_dict()

    summary_data = {
        "num_variants": num_variants,
        "num_genes": num_genes,
        "impact_distribution": impact_counts,
    }

    variants_json_path = os.path.join(output_dir, "report", "variants.json")
    os.makedirs(os.path.dirname(variants_json_path), exist_ok=True)
    df.to_json(variants_json_path, orient="records", indent=2)

    summary_json_path = os.path.join(output_dir, "report", "summary.json")
    with open(summary_json_path, "w", encoding="utf-8") as sjf:
        import json

        json.dump(summary_data, sjf, indent=2)

    logger.info(f"JSON files produced: {variants_json_path}, {summary_json_path}")
