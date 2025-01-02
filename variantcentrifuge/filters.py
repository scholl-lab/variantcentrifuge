# File: variantcentrifuge/filters.py
# Location: variantcentrifuge/variantcentrifuge/filters.py

"""
Filtering module.

This module defines functions to extract variants from a VCF file using bcftools
and apply filters via SnpSift. Each function returns a filename containing the output.
"""

import logging
import tempfile
import os
from typing import Dict, Any
from .utils import run_command

logger = logging.getLogger("variantcentrifuge")


def extract_variants(
    vcf_file: str,
    bed_file: str,
    cfg: Dict[str, Any],
    output_file: str
) -> str:
    """
    Extract variants from a VCF using bcftools and a BED file, writing output
    to the specified compressed VCF ('.vcf.gz'). bcftools is invoked with the
    '-W' option, which writes the index file automatically.

    Parameters
    ----------
    vcf_file : str
        Path to the input VCF file.
    bed_file : str
        Path to the BED file containing genomic regions of interest.
    cfg : dict
        Configuration dictionary that may include paths and parameters for tools.
        Expected keys include:
            - "threads": Number of threads to use with bcftools (default = 1).
    output_file : str
        Path to the final compressed output VCF file ('.vcf.gz').

    Returns
    -------
    str
        Path to the compressed VCF file (.vcf.gz) containing extracted variants.

    Raises
    ------
    RuntimeError
        If the extraction command fails.
    """
    threads = str(cfg.get("threads", 1))

    cmd = [
        "bcftools", "view",
        "--threads", threads,
        "-W",  # writes the index file automatically
        "-R", bed_file,
        "-Oz",  # compressed output
        "-o", output_file,
        vcf_file,
    ]
    logger.debug("Extracting variants with command: %s", " ".join(cmd))
    logger.debug("Output will be written to: %s", output_file)

    run_command(cmd)
    return output_file


def apply_snpsift_filter(
    variant_file: str,
    filter_string: str,
    cfg: Dict[str, Any],
    output_file: str
) -> str:
    """
    Apply a SnpSift filter to a variant file, then compress and index the output.

    Because our run_command function does not support shell pipelines,
    we split it into two steps:
      1) Write SnpSift filter output to a temporary uncompressed file (.vcf).
      2) Compress it with bgzip -@ <threads> to produce the final .vcf.gz.
      3) Index the resulting .vcf.gz with bcftools index.

    Parameters
    ----------
    variant_file : str
        Path to the compressed VCF file with extracted variants.
    filter_string : str
        SnpSift filter expression to apply.
    cfg : dict
        Configuration dictionary that may include paths and parameters for tools.
        Expected keys include:
            - "threads": Number of threads to use with bgzip and bcftools index (default = 1).
    output_file : str
        Path to the compressed VCF file (.vcf.gz) containing filtered variants.

    Returns
    -------
    str
        Path to the compressed VCF file (.vcf.gz) containing filtered variants.

    Raises
    ------
    RuntimeError
        If the filter command fails.
    """
    threads = str(cfg.get("threads", 1))

    # 1) Run SnpSift filter -> write uncompressed .vcf to a temporary file
    tmp_vcf = tempfile.mktemp(suffix=".vcf")
    snpsift_cmd = ["SnpSift", "filter", filter_string, variant_file]
    logger.debug("Applying SnpSift filter to produce uncompressed VCF: %s", tmp_vcf)
    run_command(snpsift_cmd, output_file=tmp_vcf)

    # 2) bgzip compress the temporary VCF
    bgzip_cmd = ["bgzip", "-@", threads, "-c", tmp_vcf]
    logger.debug("bgzip compressing to: %s", output_file)
    run_command(bgzip_cmd, output_file=output_file)

    # 3) bcftools index the resulting .vcf.gz
    index_cmd = ["bcftools", "index", "--threads", threads, output_file]
    logger.debug("Indexing filtered output with command: %s", " ".join(index_cmd))
    run_command(index_cmd)

    # Remove the uncompressed temp file
    if os.path.exists(tmp_vcf):
        os.remove(tmp_vcf)

    logger.debug("Filtered output is now available at: %s", output_file)
    return output_file
