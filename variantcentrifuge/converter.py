# File: variantcentrifuge/converter.py
# Location: variantcentrifuge/variantcentrifuge/converter.py

"""
File conversion module.

This module provides functionality to convert TSV files to XLSX format
and append additional sheets. It now also supports producing JSON files
for the HTML report.
"""

import os
import logging
from typing import Dict, Any, Optional, List
import pandas as pd
from openpyxl import load_workbook

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
    with pd.ExcelWriter(xlsx_file, engine='openpyxl', mode='a', if_sheet_exists='new') as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)


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
    # Convert all variants to a list of dicts
    variants = df.to_dict(orient="records")

    # Compute summary metrics
    num_variants = len(df)
    unique_genes = df["GENE"].unique() if "GENE" in df.columns else []
    num_genes = len(unique_genes)

    impact_counts = {}
    if "IMPACT" in df.columns:
        impact_counts = df["IMPACT"].value_counts().to_dict()

    summary_data = {
        "num_variants": num_variants,
        "num_genes": num_genes,
        "impact_distribution": impact_counts
        # Add any other summary metrics needed from df here
    }

    # Write variants.json
    variants_json_path = os.path.join(output_dir, "report", "variants.json")
    os.makedirs(os.path.dirname(variants_json_path), exist_ok=True)
    df.to_json(variants_json_path, orient="records", indent=2)

    # Write summary.json
    summary_json_path = os.path.join(output_dir, "report", "summary.json")
    with open(summary_json_path, "w") as sjf:
        import json
        json.dump(summary_data, sjf, indent=2)

    logger.info(f"JSON files produced: {variants_json_path}, {summary_json_path}")
