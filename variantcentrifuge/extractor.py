# File: variantcentrifuge/extractor.py
# Location: variantcentrifuge/variantcentrifuge/extractor.py

"""
Field extraction module.

This module provides two field extraction backends:
1. extract_fields_bcftools() - Fast C-based extraction using bcftools query (default)
2. extract_fields_snpsift() - Legacy Java-based extraction using SnpSift (fallback)
"""

import gzip
import logging
import os
import re
import subprocess
import tempfile
from collections.abc import Callable
from typing import Any

import pandas as pd

from .extractor_parsing import (  # noqa: F401 - re-exported for existing callers
    ANN_FIELD_POSITIONS,
    ANN_SUBFIELD_MAP,
    NMD_SUBFIELD_MAP,
    parse_ann_subfields,
    parse_nmd_subfields,
)
from .utils import normalize_vcf_headers, run_command

logger = logging.getLogger("variantcentrifuge")

_WILDCARD_ANN_FIELD = re.compile(r"^(ANN|EFF)\[\*\]\.")


def validate_extraction_fields(fields: list[str]) -> None:
    """Reject ANN[*]/EFF[*] subfields: only ANN[0].<FIELD> is parsed; a wildcard would pass
    through as raw annotation text and be grouped as a bogus gene name."""
    unsupported = [f for f in fields if _WILDCARD_ANN_FIELD.match(f)]
    if unsupported:
        raise ValueError(
            f"Unsupported extraction field(s): {', '.join(unsupported)}. Wildcard ANN[*]/EFF[*] "
            "subfields are not parsed by the bcftools extractor; use the supported "
            "ANN[0].<FIELD> form instead (e.g. ANN[0].GENE, ANN[0].EFFECT)."
        )


def build_bcftools_format_string(
    fields: list[str], vcf_samples: list[str] | None = None
) -> tuple[str, list[str]]:
    """
    Build bcftools query format string from field list.

    Parameters
    ----------
    fields : list[str]
        List of field names to extract (e.g., ["CHROM", "POS", "ANN[0].GENE", "GEN[*].GT"])
    vcf_samples : list[str], optional
        List of sample names from VCF for per-sample column naming

    Returns
    -------
    tuple[str, list[str]]
        Format string for bcftools query and ordered column name list
    """
    validate_extraction_fields(fields)
    format_parts = []
    column_names = []
    seen_ann = False
    seen_nmd = False
    per_sample_fields = []

    # Fixed VCF fields mapping
    fixed_fields = {
        "CHROM": "%CHROM",
        "POS": "%POS",
        "REF": "%REF",
        "ALT": "%ALT",
        "ID": "%ID",
        "FILTER": "%FILTER",
        "QUAL": "%QUAL",
    }

    for field in fields:
        if field in fixed_fields:
            # Fixed VCF field
            format_parts.append(fixed_fields[field])
            column_names.append(field)

        elif field.startswith("ANN[0]."):
            # ANN subfield - extract raw ANN once
            if not seen_ann:
                format_parts.append("%INFO/ANN")
                column_names.append("ANN")
                seen_ann = True

        elif field.startswith("NMD[0]."):
            # NMD subfield - extract raw NMD once
            if not seen_nmd:
                format_parts.append("%INFO/NMD")
                column_names.append("NMD")
                seen_nmd = True

        elif field.startswith("GEN[*]."):
            # Per-sample FORMAT field
            format_field = field.replace("GEN[*].", "")
            per_sample_fields.append(format_field)

        elif (
            field.startswith("dbNSFP_")
            or field.startswith("splice_")
            or field
            in [
                "AC",
                "RAW_FILTER",
                "ClinVar_CLNSIG",
                "hgmd_CLASS",
            ]
        ):
            # INFO field
            info_field = field
            format_parts.append(f"%INFO/{info_field}")
            column_names.append(field)

        else:
            # Unknown field - try as INFO field
            logger.warning(f"Unknown field type: {field}, treating as INFO field")
            format_parts.append(f"%INFO/{field}")
            column_names.append(field)

    # Join fixed/INFO fields with tabs
    format_string = "\\t".join(format_parts)

    # Append per-sample fields directly (NOT via format_parts join, to avoid double tab)
    if per_sample_fields:
        if vcf_samples is not None:
            if vcf_samples:
                sample_format = "\\t".join(f"%{f}" for f in per_sample_fields)
                format_string += f"[\\t{sample_format}]"
                for sample_idx in range(len(vcf_samples)):
                    for format_field in per_sample_fields:
                        column_names.append(f"GEN[{sample_idx}].{format_field}")
        else:
            sample_format = "\\t".join(f"%{f}" for f in per_sample_fields)
            format_string += f"[\\t{sample_format}]"
            logger.warning("No vcf_samples provided - per-sample columns will have generic names")
            for format_field in per_sample_fields:
                column_names.append(f"GEN[*].{format_field}")

    format_string += "\\n"

    return format_string, column_names


def get_expected_extraction_columns(
    fields: str | list[str], vcf_samples: list[str] | None = None
) -> list[str]:
    """Return the finalized column names for an extraction schema."""
    field_list = fields.strip().split() if isinstance(fields, str) else fields
    _, column_names = build_bcftools_format_string(field_list, vcf_samples)
    df = pd.DataFrame(columns=column_names)
    df = parse_ann_subfields(df, field_list)
    df = parse_nmd_subfields(df, field_list)
    return list(df.columns)


def extract_fields_bcftools(
    variant_file: str,
    fields: str,
    cfg: dict[str, Any],
    output_file: str,
    vcf_samples: list[str] | None = None,
) -> str:
    """
    Extract specified fields from variant records using bcftools query.

    This is 19x faster than SnpSift extractFields for large cohorts.

    Parameters
    ----------
    variant_file : str
        Path to the VCF file from which fields should be extracted.
    fields : str
        A space-separated list of fields to extract.
    cfg : dict
        Configuration dictionary. ``bcftools_samples`` (raw VCF names, VCF order) limits
        per-sample columns to those samples (``bcftools query -S``); ``vcf_samples`` must
        then name the same samples.
    output_file : str
        Path to the final TSV file where extracted fields will be written.
    vcf_samples : list[str], optional
        List of sample names from VCF for per-sample column naming.

    Returns
    -------
    str
        The output_file path that now contains the extracted fields (TSV).

    Raises
    ------
    RuntimeError
        If the bcftools query command fails.
    """
    field_list = fields.strip().split()
    logger.debug(f"Field list to extract: {field_list}")

    # Build bcftools format string
    format_string, column_names = build_bcftools_format_string(field_list, vcf_samples)
    logger.debug(f"bcftools format string: {format_string}")
    logger.debug(f"Column names: {column_names}")

    # Create temporary file for bcftools output
    temp_fd, temp_output = tempfile.mkstemp(
        suffix=".tsv", prefix="bcftools_extract_", dir=cfg.get("temp_dir")
    )
    os.close(temp_fd)
    samples_file = None

    try:
        # Run bcftools query
        cmd = ["bcftools", "query", "-u", "-f", format_string]
        if cfg.get("bcftools_samples"):  # analysed-sample subset (#212), names in VCF order
            samples_file = f"{temp_output}.samples"
            with open(samples_file, "w") as fh:
                fh.write("".join(f"{s}\n" for s in cfg["bcftools_samples"]))
            cmd.extend(["-S", samples_file])
        cmd.append(variant_file)

        logger.debug(f"Running bcftools query: {' '.join(cmd)}")

        with open(temp_output, "w") as out_f:
            result = subprocess.run(
                cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, check=False
            )

        if result.returncode != 0:
            error_msg = result.stderr if result.stderr else "Unknown error"
            raise RuntimeError(
                f"bcftools query failed with exit code {result.returncode}: {error_msg}"
            )

        # Read raw output with pandas
        logger.debug(f"Reading bcftools output from {temp_output}")

        # Handle empty output (0 variants after filtering)
        if os.path.getsize(temp_output) == 0:
            logger.warning(
                "Field extraction (bcftools query): 0 variants extracted from %s", variant_file
            )
            df = pd.DataFrame(columns=column_names)
            df = parse_ann_subfields(df, field_list)
            df = parse_nmd_subfields(df, field_list)
            # Write header-only output
            if output_file.endswith(".gz"):
                with gzip.open(output_file, "wt", compresslevel=1) as f:
                    df.to_csv(f, sep="\t", index=False, na_rep="NA")
            else:
                df.to_csv(output_file, sep="\t", index=False, na_rep="NA")
            return output_file

        df = pd.read_csv(
            temp_output,
            sep="\t",
            header=None,
            dtype=str,
            keep_default_na=False,
        )

        # Assign column names
        if len(column_names) != len(df.columns):
            logger.warning(
                f"Column count mismatch: expected {len(column_names)}, got {len(df.columns)}"
            )
            # Adjust column names if mismatch
            if len(column_names) < len(df.columns):
                # Add generic names for extra columns
                for i in range(len(column_names), len(df.columns)):
                    column_names.append(f"EXTRA_{i}")
            else:
                # Truncate column names
                column_names = column_names[: len(df.columns)]

        df.columns = column_names
        if (
            "ALT" in df.columns
            and df["ALT"].astype(str).str.contains(",", na=False).any()
            and any(c.endswith(("__AF", ".AF", "__AD", ".AD")) for c in df.columns)
        ):
            raise ValueError("Unsplit multiallelic ALT record encountered during AF/AD extraction")

        # Parse ANN subfields if present
        df = parse_ann_subfields(df, field_list)

        # Parse NMD subfields if present
        df = parse_nmd_subfields(df, field_list)

        # Normalize missing values: replace empty strings and "." with "NA"
        # (preserve "." in GT columns)
        df = df.fillna("NA")
        df = df.replace("", "NA")
        gt_cols = {
            c
            for c in df.columns
            if c == "GT" or c.endswith((".GT", "__GT")) or re.match(r"^GEN\[\d+\]\.GT$", c)
        }
        non_gt_cols = [c for c in df.columns if c not in gt_cols]
        if non_gt_cols:
            df[non_gt_cols] = df[non_gt_cols].replace(".", "NA")

        # Write output TSV
        logger.debug(f"Writing output to {output_file}")

        # Determine if we need compression
        if output_file.endswith(".gz"):
            with gzip.open(output_file, "wt", compresslevel=1) as f:
                df.to_csv(f, sep="\t", index=False, na_rep="NA")
        else:
            df.to_csv(output_file, sep="\t", index=False, na_rep="NA")

        logger.info(
            "Field extraction (bcftools query): %d variants extracted from %s to %s",
            len(df),
            variant_file,
            output_file,
        )

    finally:
        # Clean up temp files
        for p in (temp_output, samples_file):
            if p and os.path.exists(p):
                os.remove(p)

    return output_file


def extract_fields_snpsift(
    variant_file: str, fields: str, cfg: dict[str, Any], output_file: str
) -> str:
    """
    Extract specified fields from variant records using SnpSift extractFields (legacy).

    This is the original Java-based extraction method, kept as a fallback.

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

    from .filters import _java_tmp_flags, _snpsift_memory_flag

    xmx = _snpsift_memory_flag(cfg)
    cmd = [
        "SnpSift",
        xmx,
        *_java_tmp_flags(cfg),
        "extractFields",
        "-s",
        snpsift_sep,  # SnpSift subfield separator
        "-e",
        "NA",  # Replace missing values with "NA"
        variant_file,
        *field_list,
    ]

    logger.debug("Running SnpSift (%s) with command: %s", xmx, " ".join(cmd))

    # Write to a temporary uncompressed file first
    temp_output = output_file.replace(".gz", "")
    logger.debug("Extracted fields will be written to temporary file: %s", temp_output)

    # Run SnpSift extractFields, writing to temp file
    run_command(cmd, output_file=temp_output)

    # Now fix up the header line and compress
    # Use fast compression (level 1) for intermediate files to optimize I/O performance
    def get_open_func() -> tuple[Callable[..., Any], str]:
        if output_file.endswith(".gz"):

            def compressed_open(f: Any, m: str, **kwargs: Any) -> Any:
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


# For backwards compatibility, keep extract_fields as alias to SnpSift version
extract_fields = extract_fields_snpsift
