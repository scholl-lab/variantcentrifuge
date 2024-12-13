# File: variantcentrifuge/converter.py
# Location: variantcentrifuge/variantcentrifuge/converter.py

"""
File conversion module.

This module provides functionality to convert TSV files to XLSX format,
and append additional sheets.
"""

import pandas as pd
import os
from openpyxl import load_workbook

def convert_to_excel(tsv_file, cfg):
    """
    Convert a TSV file to XLSX format with a single "Results" sheet.

    Parameters
    ----------
    tsv_file : str
        Path to the TSV file.
    cfg : dict
        Configuration dictionary.

    Returns
    -------
    str
        Path to the generated XLSX file.
    """
    df = pd.read_csv(tsv_file, sep="\t", na_values="NA")
    xlsx_file = os.path.splitext(tsv_file)[0] + ".xlsx"
    df.to_excel(xlsx_file, index=False, sheet_name="Results")
    return xlsx_file

def append_tsv_as_sheet(xlsx_file, tsv_file, sheet_name="Metadata"):
    """
    Append a TSV file as a new sheet to an existing XLSX file.

    Assumes the TSV has a header row.

    Parameters
    ----------
    xlsx_file : str
        Path to the existing XLSX file.
    tsv_file : str
        Path to the TSV file to append.
    sheet_name : str
        Name of the new sheet.
    """
    df = pd.read_csv(tsv_file, sep="\t", header=0)  # first line is header
    with pd.ExcelWriter(xlsx_file, engine='openpyxl', mode='a', if_sheet_exists='new') as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)
