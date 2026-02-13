# File: variantcentrifuge/extractor.py
# Location: variantcentrifuge/variantcentrifuge/extractor.py

"""
Field extraction module.

This module uses SnpSift extractFields to extract specific fields
from VCF records and save the result to the specified output file.
"""

import gzip
import logging
import os
from typing import Any

from .utils import normalize_vcf_headers, run_command

logger = logging.getLogger("variantcentrifuge")


def extract_fields(variant_file: str, fields: str, cfg: dict[str, Any], output_file: str) -> str:
    """
    Extract specified fields from variant records.

    Write them directly to `output_file`, controlling the SnpSift field separator if needed.

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
        "SnpSift",
        "extractFields",
        "-s",
        snpsift_sep,  # SnpSift subfield separator
        "-e",
        "NA",  # Replace missing values with "NA"
        variant_file,
    ] + field_list

    logger.debug("Running SnpSift with command: %s", " ".join(cmd))

    # Write to a temporary uncompressed file first
    temp_output = output_file.replace(".gz", "")
    logger.debug("Extracted fields will be written to temporary file: %s", temp_output)

    # Run SnpSift extractFields, writing to temp file
    run_command(cmd, output_file=temp_output)

    # Now fix up the header line and compress
    # Use fast compression (level 1) for intermediate files to optimize I/O performance
    def get_open_func():
        if output_file.endswith(".gz"):

            def compressed_open(f, m, **kwargs):
                return gzip.open(f, m, compresslevel=1, **kwargs)

            return compressed_open, "wt"
        else:
            return open, "w"

    open_func, mode = get_open_func()

    with open(temp_output, encoding="utf-8") as f:
        lines = f.readlines()

    if not lines:
        logger.warning("No lines were written to the output after SnpSift extract. Check input.")
        if output_file.endswith(".gz") and os.path.exists(temp_output):
            os.remove(temp_output)
        return output_file

    # Use the utility function to normalize VCF headers, including indexed field renaming
    lines = normalize_vcf_headers(lines)

    # Write the compressed output with optimized compression
    with open_func(output_file, mode, encoding="utf-8") as f:
        f.writelines(lines)

    # Clean up temporary file
    if output_file.endswith(".gz") and os.path.exists(temp_output):
        os.remove(temp_output)

    logger.debug("SnpSift extractFields completed successfully. Output written to: %s", output_file)
    return output_file
