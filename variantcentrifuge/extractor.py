# File: variantcentrifuge/extractor.py
# Location: variantcentrifuge/variantcentrifuge/extractor.py

"""
Field extraction module.

This module uses SnpSift extractFields to extract specific fields
from VCF records and save the result to the specified output file.
"""

import logging
from typing import Dict, Any
from .utils import run_command

logger = logging.getLogger("variantcentrifuge")


def extract_fields(
    variant_file: str,
    fields: str,
    cfg: Dict[str, Any],
    output_file: str
) -> str:
    """
    Extract specified fields from variant records and write them directly to `output_file`.

    Parameters
    ----------
    variant_file : str
        Path to the VCF file from which fields should be extracted.
    fields : str
        A space-separated list of fields to extract.
    cfg : dict
        Configuration dictionary that may include tool paths, parameters, etc.
    output_file : str
        Path to the final TSV file where extracted fields will be written.

    Returns
    -------
    str
        The same `output_file` path that now contains the extracted fields (TSV).

    Raises
    ------
    RuntimeError
        If the field extraction command fails.
    """
    field_list = fields.strip().split()

    cmd = ["SnpSift", "extractFields", "-s", ",", "-e", "NA", variant_file] + field_list
    logger.debug("Extracting fields with command: %s", " ".join(cmd))
    logger.debug("Output will be written to: %s", output_file)

    # Run SnpSift extractFields, writing directly to output_file
    run_command(cmd, output_file=output_file)

    # Now we fix up the header line in place
    with open(output_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if lines:
        # Clean up the header line by removing certain prefixes
        lines[0] = lines[0].replace("ANN[0].", "").replace("GEN[*].", "")

    with open(output_file, "w", encoding="utf-8") as f:
        f.writelines(lines)

    return output_file
