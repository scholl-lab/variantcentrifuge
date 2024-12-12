# File: variantcentrifuge/filters.py
# Location: variantcentrifuge/variantcentrifuge/filters.py

"""
Filtering module.

This module defines functions to extract variants from a VCF file using bcftools
and apply filters via SnpSift. Each function returns a filename containing the output.
"""

import tempfile
from .utils import run_command, log_message

def extract_variants(vcf_file, bed_file, cfg):
    """
    Extract variants from a VCF using bcftools and a BED file, saving output to a temp file.

    Parameters
    ----------
    vcf_file : str
        Path to the input VCF.
    bed_file : str
        Path to the BED file of regions.
    cfg : dict

    Returns
    -------
    str
        Path to the temporary file containing extracted variants.
    """
    output_file = tempfile.mktemp(suffix=".vcf")
    cmd = ["bcftools", "view", vcf_file, "-R", bed_file]
    log_message("DEBUG", f"Extracting variants to {output_file}")
    run_command(cmd, output_file=output_file)
    return output_file

def apply_snpsift_filter(variant_file, filter_string, cfg):
    """
    Apply a SnpSift filter to a variant file and save output to a temp file.

    Parameters
    ----------
    variant_file : str
        Path to the VCF file with extracted variants.
    filter_string : str
        SnpSift filter expression.
    cfg : dict

    Returns
    -------
    str
        Path to the temporary file containing filtered variants.
    """
    output_file = tempfile.mktemp(suffix=".vcf")
    cmd = ["SnpSift", "filter", filter_string, variant_file]
    log_message("DEBUG", f"Applying SnpSift filter to produce {output_file}")
    run_command(cmd, output_file=output_file)
    return output_file
