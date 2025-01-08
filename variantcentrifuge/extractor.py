# File: variantcentrifuge/extractor.py
# Location: variantcentrifuge/variantcentrifuge/extractor.py

"""
Field extraction module.

This module uses SnpSift extractFields to extract specific fields
from VCF records and save the result to the specified output file.
"""

import logging
from typing import Dict, Any
from .utils import run_command, normalize_snpeff_headers  # <-- Updated import

logger = logging.getLogger("variantcentrifuge")


def extract_fields(
    variant_file: str,
    fields: str,
    cfg: Dict[str, Any],
    output_file: str
) -> str:
    """
    Extract specified fields from variant records and write them directly to
    `output_file`, controlling the SnpSift field separator if needed.

    Parameters
    ----------
    variant_file : str
        Path to the VCF file from which fields should be extracted.
    fields : str
        A space-separated list of fields to extract (e.g. "CHROM POS REF ALT DP AD").
    cfg : dict
        Configuration dictionary that may include tool paths, parameters, etc.
        - "extract_fields_separator": str
            The separator for multi-sample fields when using SnpSift `-s ...`.
            Often a comma ",". Defaults to "," if not present.
        - "debug_level": str or None
            Optional debug level to control how much we log.
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
    # Pull the user-defined or default multi-sample separator from config
    # E.g. if we want DP,AD per sample, each sample's subfields will be separated by this
    snpsift_sep = cfg.get("extract_fields_separator", ":")
    logger.debug(f"Using SnpSift multi-sample separator: '{snpsift_sep}'")

    field_list = fields.strip().split()
    logger.debug(f"Field list to extract: {field_list}")

    cmd = [
        "SnpSift", "extractFields",
        "-s", snpsift_sep,     # SnpSift subfield separator
        "-e", "NA",            # Replace missing values with "NA"
        variant_file
    ] + field_list

    logger.debug("Running SnpSift with command: %s", " ".join(cmd))
    logger.debug("Extracted fields will be written to: %s", output_file)

    # Run SnpSift extractFields, writing directly to output_file
    run_command(cmd, output_file=output_file)

    # Now fix up the header line in place
    with open(output_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if not lines:
        logger.warning("No lines were written to the output after SnpSift extract. Check input.")
        return output_file

    # Use the new utility function to remove any SnpEff prefixes from the header
    lines = normalize_snpeff_headers(lines)

    # Rewrite the file with updated lines
    with open(output_file, "w", encoding="utf-8") as f:
        f.writelines(lines)

    logger.debug("SnpSift extractFields completed successfully.")
    return output_file
