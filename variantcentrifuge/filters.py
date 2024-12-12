# File: variantcentrifuge/filters.py
# Location: variantcentrifuge/variantcentrifuge/filters.py

"""
Filtering module.

This module defines functions to extract variants from a VCF file using
bcftools and apply filters via SnpSift.
"""

from .utils import run_command_stream


def extract_variants(vcf_file, bed_file, cfg):
    """
    Extract variants from a VCF using bcftools and a BED file.

    Parameters
    ----------
    vcf_file : str
        Path to the input VCF.
    bed_file : str
        Path to the BED file of regions.
    cfg : dict
        Configuration dictionary.

    Returns
    -------
    iterator
        An iterator over the lines of the extracted VCF stream.
    """
    cmd = ["bcftools", "view", vcf_file, "-R", bed_file]
    return run_command_stream(cmd)


def apply_snpsift_filter(variant_stream, filter_string, cfg):
    """
    Apply a SnpSift filter to a variant stream.

    Parameters
    ----------
    variant_stream : iterator
        An iterator of VCF lines.
    filter_string : str
        SnpSift filter expression.
    cfg : dict
        Configuration dictionary.

    Returns
    -------
    iterator
        An iterator over filtered variant lines.
    """
    cmd = ["SnpSift", "filter", filter_string]
    return run_command_stream(cmd, input_stream=variant_stream)
