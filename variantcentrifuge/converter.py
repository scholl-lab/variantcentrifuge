# File: variantcentrifuge/converter.py
# Location: variantcentrifuge/variantcentrifuge/converter.py

"""
File conversion module.

This module provides functionality to convert TSV files to XLSX format.
"""

import pandas as pd
import os


def convert_to_excel(tsv_file, cfg):
    """
    Convert a TSV file to XLSX format.

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
