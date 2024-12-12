# File: variantcentrifuge/extractor.py
# Location: variantcentrifuge/variantcentrifuge/extractor.py

"""
Field extraction module.

This module uses SnpSift extractFields to extract specific fields
from VCF records and save the result to a temporary file.
"""

import tempfile
from .utils import run_command, log_message

def extract_fields(variant_file, fields, cfg):
    """
    Extract specified fields from variant records and save output to a temp file.

    Parameters
    ----------
    variant_file : str
        Path to the VCF file to extract fields from.
    fields : str
        Space-separated list of fields to extract.
    cfg : dict

    Returns
    -------
    str
        Path to the temporary file containing extracted fields (TSV).
    """
    field_list = fields.strip().split()
    output_file = tempfile.mktemp(suffix=".tsv")
    cmd = ["SnpSift", "extractFields", "-s", ",", "-e", "NA", variant_file] + field_list
    log_message("DEBUG", f"Extracting fields to {output_file}")
    run_command(cmd, output_file=output_file)

    # Clean up header line (sed replacement), can do in Python if needed
    # For simplicity, let's just replicate the sed logic here:
    # sed -e '1s/ANN\[0\]\.//g; s/GEN\[\*\]\.//g'
    # We'll do this in Python:
    with open(output_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if lines:
        header = lines[0].replace("ANN[0].", "").replace("GEN[*].", "")
        lines[0] = header

    with open(output_file, "w", encoding="utf-8") as f:
        f.writelines(lines)

    return output_file
