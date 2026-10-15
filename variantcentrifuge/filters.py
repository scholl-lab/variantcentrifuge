# File: variantcentrifuge/filters.py
# Location: variantcentrifuge/variantcentrifuge/filters.py

"""
Filtering module.

This module defines functions to extract variants from a VCF file using bcftools
and apply filters via SnpSift. Each function returns a filename containing the output.

Added:

- A function to filter final TSV rows by genotype (het, hom, comp_het), optionally
  applying different rules per gene based on a separate mapping file.
- Enhanced to append a reason label (e.g., "(comphet)") to each sample genotype if it
  passed the filter because of 'comp_het' or 'het' or 'hom'.

New changes for preserving extra fields (e.g., DP, AD) in parentheses:

- If a sample's entry is "325879(0/1:53,55:108)", we parse "0/1" as the "main" genotype
  but leave the rest of the substring (":53,55:108") intact. If this sample passes
  'het', we produce "325879(0/1:53,55:108)(het)" in the final output.

TSV-based filtering:

- Added filter_tsv_with_expression() to apply filters on TSV data after scoring
  and annotation steps are complete, allowing filtering on computed columns.
"""

import logging
import os
import tempfile
from typing import Any

from .transcript_filter import apply_filter_exceptions_to_vcf, load_filter_exceptions
from .tsv_filters import (
    filter_dataframe_with_query,
    filter_final_tsv_by_genotype,
    filter_tsv_with_expression,
)
from .utils import run_command

logger = logging.getLogger("variantcentrifuge")

__all__ = [
    "apply_bcftools_prefilter",
    "apply_snpsift_filter",
    "extract_variants",
    "filter_dataframe_with_query",
    "filter_final_tsv_by_genotype",
    "filter_tsv_with_expression",
]


def _snpsift_memory_flag(cfg: dict[str, Any]) -> str:
    """Calculate JVM -Xmx flag for SnpSift based on available resources.

    Uses ResourceManager to detect system memory and allocates up to 25%
    for SnpSift, bounded between 2GB and 16GB. Falls back to 4GB on error.

    The conda SnpSift wrapper defaults to only -Xmx1g which is insufficient
    for large multi-sample VCFs (e.g., 5000+ samples).
    """
    try:
        from .memory.resource_manager import ResourceManager

        rm = ResourceManager(config=cfg)
        # Give SnpSift up to 25% of available memory, min 2GB, max 16GB
        memory_gb = max(2, min(16, int(rm.memory_gb * 0.25)))
    except Exception:
        memory_gb = 4  # Safe default
    return f"-Xmx{memory_gb}g"


SNPSIFT_FATAL_STDERR_PATTERNS = (
    "mismatched input",
    "missing ')'",
    "token recognition error",
    "Error parsing",
    "Cannot parse",
    "LexerNoViableAltException",
    "Unknown parameter",
    "INFO field",
    'Exception in thread "main"',
)


def _raise_for_snpsift_filter_stderr(stderr: str | None, filter_string: str) -> None:
    """Raise when SnpSift reports parser/evaluation diagnostics despite exit code 0."""
    if not stderr:
        return

    diagnostic_lines = [
        line.strip()
        for line in stderr.splitlines()
        if any(pattern in line for pattern in SNPSIFT_FATAL_STDERR_PATTERNS)
    ]
    if not diagnostic_lines:
        return

    preview = "\n".join(diagnostic_lines[:5])
    logger.error("SnpSift filter diagnostics detected:\n%s", preview)
    raise RuntimeError(
        "SnpSift filter reported parser diagnostics while exiting successfully. "
        f"Filter expression: {filter_string}\n{preview}"
    )


def apply_bcftools_prefilter(
    input_vcf: str, output_vcf: str, filter_expression: str, cfg: dict[str, Any]
) -> str:
    """
    Apply a bcftools filter expression to a VCF file.

    This function is used for pre-filtering VCF files before more resource-intensive
    operations like SnpSift filtering. It applies the filter expression using
    bcftools view -i option.

    Parameters
    ----------
    input_vcf : str
        Path to the input VCF file.
    output_vcf : str
        Path to the output filtered VCF file.
    filter_expression : str
        bcftools filter expression. Note that bcftools uses different syntax than SnpSift.
        Examples include 'FILTER="PASS"' for only PASS variants, 'INFO/AC<10' for allele
        count less than 10, 'FILTER="PASS" && INFO/AC<10' for combined filters, or
        'INFO/AF<0.01 || INFO/AC<5' for AF less than 1% OR AC less than 5.
    cfg : dict
        Configuration dictionary that may include "threads" (number of threads to use
        with bcftools, default is 1).

    Returns
    -------
    str
        Path to the output VCF file.

    Raises
    ------
    RuntimeError
        If the bcftools command fails.
    """
    threads = str(cfg.get("threads", 1))

    cmd = [
        "bcftools",
        "view",
        "--threads",
        threads,
        "-i",
        filter_expression,
        "-Oz",  # Output compressed VCF
        "-o",
        output_vcf,
        input_vcf,
    ]

    in_count = _count_vcf_records(input_vcf, cached=True)
    logger.info(f"Applying bcftools pre-filter: {filter_expression}")
    logger.debug(f"bcftools prefilter command: {' '.join(cmd)}")
    run_command(cmd)

    # bcftools view does not automatically index, so we must do it
    index_cmd = ["bcftools", "index", "--threads", threads, output_vcf]
    logger.debug(f"Indexing pre-filtered VCF: {' '.join(index_cmd)}")
    run_command(index_cmd)

    out_count = _count_vcf_records(output_vcf)
    logger.info(
        "bcftools pre-filter: %d input -> %d output records (%s)",
        in_count,
        out_count,
        output_vcf,
    )
    return output_vcf


_VCF_RECORD_COUNT_CACHE: dict[tuple[str, int, int], int] = {}


def _count_vcf_records(path: str, *, cached: bool = False) -> int:
    """Count non-header records in a plain or gzipped VCF file (0 if unreadable)."""
    if not path or not os.path.exists(path):
        return 0
    cache_key: tuple[str, int, int] | None = None
    if cached:
        try:
            st = os.stat(path)
            cache_key = (os.path.abspath(path), st.st_mtime_ns, st.st_size)
            if cache_key in _VCF_RECORD_COUNT_CACHE:
                return _VCF_RECORD_COUNT_CACHE[cache_key]
        except OSError:
            cache_key = None
    import gzip

    count = 0
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            count = sum(1 for line in fh if line and not line.startswith("#"))
    except Exception:
        try:
            with open(path, errors="replace") as fh:
                count = sum(1 for line in fh if line and not line.startswith("#"))
        except Exception:
            count = 0
    if cache_key is not None:
        _VCF_RECORD_COUNT_CACHE[cache_key] = count
    return count


def extract_variants(vcf_file: str, bed_file: str, cfg: dict[str, Any], output_file: str) -> str:
    """
    Extract variants from a VCF using bcftools and a BED file.

    Write output to the specified compressed VCF ('.vcf.gz') after decomposing
    multiallelic records with ``bcftools norm -m -any``. bcftools is invoked with the
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
            - "bcftools_prefilter": Optional bcftools filter expression to apply during extraction.
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
    from .vcf_eff_one_per_line import restrict_vcf_ann_to_alt

    threads = str(cfg.get("threads", 1))

    tmp_fd, tmp_vcf = tempfile.mkstemp(suffix=".vcf.gz")
    os.close(tmp_fd)
    try:
        cmd = [
            "bcftools",
            "view",
            "--threads",
            threads,
            "-R",
            bed_file,
        ]

        # Add optional pre-filter expression if provided
        if cfg.get("bcftools_prefilter"):
            cmd.extend(["-i", cfg["bcftools_prefilter"]])
            logger.info(
                f"Applying bcftools pre-filter during extraction: {cfg['bcftools_prefilter']}"
            )

        cmd.extend(
            [
                "-Oz",  # compressed output
                "-o",
                tmp_vcf,
                vcf_file,
            ]
        )

        logger.debug("Extracting variants with command: %s", " ".join(cmd))
        run_command(cmd)

        norm_cmd = [
            "bcftools",
            "norm",
            "--threads",
            threads,
            "-W",  # writes the index file automatically
            "-m",
            "-any",
            "-Oz",
            "-o",
            output_file,
            tmp_vcf,
        ]
        logger.debug("Decomposing multiallelic variants with command: %s", " ".join(norm_cmd))
        logger.debug("Output will be written to: %s", output_file)
        run_command(norm_cmd)
        restrict_vcf_ann_to_alt(output_file)

        out_count = _count_vcf_records(output_file)
        in_count = (
            _count_vcf_records(tmp_vcf) or _count_vcf_records(vcf_file, cached=True) or out_count
        )
        logger.info(
            "Variant extraction (bcftools view + norm): %d input -> %d output records (%s)",
            in_count,
            out_count,
            output_file,
        )
        return output_file
    finally:
        for p in (tmp_vcf, f"{tmp_vcf}.csi", f"{tmp_vcf}.tbi"):
            if os.path.exists(p):
                os.remove(p)


def apply_snpsift_filter(
    variant_file: str, filter_string: str, cfg: dict[str, Any], output_file: str
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
    from .vcf_eff_one_per_line import restrict_vcf_ann_to_alt

    threads = str(cfg.get("threads", 1))
    in_count = _count_vcf_records(variant_file)
    exceptions = cfg.get("filter_exceptions") or load_filter_exceptions(
        cfg.get("filter_exception_file")
    )

    # 1) Run SnpSift filter -> write uncompressed .vcf to a temporary file
    tmp_fd, tmp_vcf = tempfile.mkstemp(suffix=".vcf")
    os.close(tmp_fd)
    rescued_vcf: str | None = None
    ann_vcf: str | None = None
    try:
        snpsift_input = variant_file
        if exceptions and os.path.exists(snpsift_input):
            rescued_vcf = f"{tmp_vcf}.rescued.vcf"
            apply_filter_exceptions_to_vcf(snpsift_input, rescued_vcf, exceptions)
            snpsift_input = rescued_vcf
        if os.path.exists(snpsift_input):
            candidate_ann_vcf = f"{tmp_vcf}.ann.vcf"
            if restrict_vcf_ann_to_alt(snpsift_input, candidate_ann_vcf):
                ann_vcf = candidate_ann_vcf
                snpsift_input = ann_vcf
            elif os.path.exists(candidate_ann_vcf):
                os.remove(candidate_ann_vcf)

        xmx = _snpsift_memory_flag(cfg)
        snpsift_cmd = ["SnpSift", xmx, "filter", filter_string, snpsift_input]
        logger.debug("Applying SnpSift filter (%s) to produce uncompressed VCF: %s", xmx, tmp_vcf)
        result = run_command(snpsift_cmd, output_file=tmp_vcf, return_result=True)
        _raise_for_snpsift_filter_stderr(result.stderr, filter_string)

        # 2) bgzip compress the temporary VCF
        bgzip_cmd = ["bgzip", "-@", threads, "-c", tmp_vcf]
        logger.debug("bgzip compressing to: %s", output_file)
        run_command(bgzip_cmd, output_file=output_file)

        # 3) bcftools index the resulting .vcf.gz
        index_cmd = ["bcftools", "index", "--threads", threads, output_file]
        logger.debug("Indexing filtered output with command: %s", " ".join(index_cmd))
        run_command(index_cmd)

        out_count = _count_vcf_records(tmp_vcf) or _count_vcf_records(output_file)
        logger.info(
            "SnpSift filter: %d input -> %d output records (%s)",
            in_count,
            out_count,
            output_file,
        )
        logger.debug("Filtered output is now available at: %s", output_file)
        return output_file
    finally:
        # Remove the uncompressed temp files
        for p in (tmp_vcf, rescued_vcf, ann_vcf):
            if p and os.path.exists(p):
                os.remove(p)
