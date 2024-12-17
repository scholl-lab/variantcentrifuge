# File: variantcentrifuge/filters.py
# Location: variantcentrifuge/variantcentrifuge/filters.py

"""
Filtering module.

This module defines functions to extract variants from a VCF file using bcftools
and apply filters via SnpSift. Each function returns a filename containing the output.
"""

import tempfile
import logging
from typing import Dict, Any
from .utils import run_command

logger = logging.getLogger("variantcentrifuge")


def extract_variants(vcf_file: str, bed_file: str, cfg: Dict[str, Any]) -> str:
    """
    Extract variants from a VCF using bcftools and a BED file, saving output to a temp file.

    Parameters
    ----------
    vcf_file : str
        Path to the input VCF file.
    bed_file : str
        Path to the BED file containing genomic regions of interest.
    cfg : dict
        Configuration dictionary that may include paths and parameters for tools.

    Returns
    -------
    str
        Path to the temporary VCF file containing extracted variants.

    Raises
    ------
    RuntimeError
        If the extraction command fails.
    """
    output_file = tempfile.mktemp(suffix=".vcf")
    cmd = ["bcftools", "view", vcf_file, "-R", bed_file]
    logger.debug("Extracting variants with command: %s", " ".join(cmd))
    logger.debug("Output will be written to: %s", output_file)

    run_command(cmd, output_file=output_file)
    return output_file


def apply_snpsift_filter(variant_file: str, filter_string: str, cfg: Dict[str, Any]) -> str:
    """
    Apply a SnpSift filter to a variant file and save output to a temp file.

    Parameters
    ----------
    variant_file : str
        Path to the VCF file with extracted variants.
    filter_string : str
        SnpSift filter expression to apply.
    cfg : dict
        Configuration dictionary that may include paths and parameters for tools.

    Returns
    -------
    str
        Path to the temporary VCF file containing filtered variants.

    Raises
    ------
    RuntimeError
        If the filter command fails.
    """
    output_file = tempfile.mktemp(suffix=".vcf")
    cmd = ["SnpSift", "filter", filter_string, variant_file]
    logger.debug("Applying SnpSift filter with command: %s", " ".join(cmd))
    logger.debug("Filtered output will be written to: %s", output_file)

    run_command(cmd, output_file=output_file)
    return output_file
