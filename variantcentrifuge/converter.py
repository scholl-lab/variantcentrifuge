# File: variantcentrifuge/converter.py
# Location: variantcentrifuge/variantcentrifuge/converter.py

"""
File conversion module.

This module provides functionality to convert TSV files to XLSX format
and append additional sheets.
"""

import os
import logging
from typing import Dict, Any, Optional

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
