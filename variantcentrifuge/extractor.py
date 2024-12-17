# File: variantcentrifuge/extractor.py
# Location: variantcentrifuge/variantcentrifuge/extractor.py

"""
Field extraction module.

This module uses SnpSift extractFields to extract specific fields
from VCF records and save the result to a temporary file.
"""

import tempfile
import logging
from typing import Dict, Any
from .utils import run_command

logger = logging.getLogger("variantcentrifuge")


def extract_fields(variant_file: str, fields: str, cfg: Dict[str, Any]) -> str:
    """
    Extract specified fields from variant records and save output to a temporary file.

    Parameters
    ----------
    variant_file : str
        Path to the VCF file from which fields should be extracted.
    fields : str
        A space-separated list of fields to extract.
    cfg : dict
        Configuration dictionary that may include tool paths, parameters, etc.

    Returns
    -------
    str
        Path to the temporary file containing extracted fields (TSV).

    Raises
    ------
    RuntimeError
        If the field extraction command fails.
    """
    field_list = fields.strip().split()
    output_file = tempfile.mktemp(suffix=".tsv")

    cmd = ["SnpSift", "extractFields", "-s", ",", "-e", "NA", variant_file] + field_list
    logger.debug("Extracting fields with command: %s", " ".join(cmd))
    logger.debug("Output will be written to: %s", output_file)

    run_command(cmd, output_file=output_file)

    with open(output_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if lines:
        # Clean up the header line by removing certain prefixes
        header = lines[0].replace("ANN[0].", "").replace("GEN[*].", "")
        lines[0] = header

    with open(output_file, "w", encoding="utf-8") as f:
        f.writelines(lines)

    return output_file
