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
import subprocess
import tempfile
from typing import Any

import pandas as pd

from .utils import normalize_vcf_headers, run_command

logger = logging.getLogger("variantcentrifuge")

# SnpEff ANN field positions (pipe-delimited)
# See: http://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files
ANN_FIELD_POSITIONS = {
    "Allele": 0,
    "Annotation": 1,  # EFFECT
    "Annotation_Impact": 2,  # IMPACT
    "Gene_Name": 3,  # GENE
    "Gene_ID": 4,
    "Feature_Type": 5,
    "Feature_ID": 6,  # FEATUREID
    "Transcript_BioType": 7,
    "Rank": 8,
    "HGVS.c": 9,  # HGVS_C
    "HGVS.p": 10,  # HGVS_P
    "cDNA.pos/cDNA.length": 11,
    "CDS.pos/CDS.length": 12,
    "AA.pos/AA.length": 13,  # Contains both AA_POS and AA_LEN
    "Distance": 14,
    "ERRORS/WARNINGS/INFO": 15,
}

# Map config field names to ANN positions
ANN_SUBFIELD_MAP = {
    "ANN[0].GENE": ("Gene_Name", 3),
    "ANN[0].FEATUREID": ("Feature_ID", 6),
    "ANN[0].EFFECT": ("Annotation", 1),
    "ANN[0].IMPACT": ("Annotation_Impact", 2),
    "ANN[0].HGVS_C": ("HGVS.c", 9),
    "ANN[0].HGVS_P": ("HGVS.p", 10),
    "ANN[0].AA_POS": ("AA.pos/AA.length", 13),  # Will need split
    "ANN[0].AA_LEN": ("AA.pos/AA.length", 13),  # Will need split
}

# NMD field positions (pipe-delimited, similar to ANN)
NMD_SUBFIELD_MAP = {
    "NMD[0].PERC": 0,  # First field is the percentage
}


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

    # Add per-sample fields at the end
    if per_sample_fields:
        # Build per-sample format string
        # Format: [\t%GT\t%DP] outputs GT, DP for each sample in separate columns
        sample_format = "\\t".join(f"%{f}" for f in per_sample_fields)
        format_parts.append(f"[\\t{sample_format}]")

        # Add column names for each sample
        if vcf_samples:
            for sample in vcf_samples:
                for format_field in per_sample_fields:
                    # Column name format: GEN[0].GT, GEN[1].GT, etc.
                    sample_idx = vcf_samples.index(sample)
                    column_names.append(f"GEN[{sample_idx}].{format_field}")
        else:
            logger.warning("No vcf_samples provided - per-sample columns will have generic names")
            # Just use field names without sample info
            for format_field in per_sample_fields:
                column_names.append(f"GEN[*].{format_field}")

    # Join with tabs
    format_string = "\\t".join(format_parts) + "\\n"

    return format_string, column_names


def parse_ann_subfields(df: pd.DataFrame, fields: list[str]) -> pd.DataFrame:
    """
    Parse ANN subfields from raw ANN column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with raw "ANN" column containing pipe-delimited annotations
    fields : list[str]
        Original field list to determine which ANN subfields to extract

    Returns
    -------
    pd.DataFrame
        DataFrame with ANN subfields added and raw ANN column removed
    """
    if "ANN" not in df.columns:
        return df

    # Determine which ANN subfields are needed
    ann_fields_needed = [f for f in fields if f.startswith("ANN[0].")]

    if not ann_fields_needed:
        # No ANN subfields requested, just drop the raw column
        return df.drop(columns=["ANN"])

    # Take first annotation (index [0])
    # Multiple annotations are comma-separated
    first_ann = df["ANN"].fillna("").str.split(",", n=1).str[0]

    # Split by pipe
    ann_split = first_ann.str.split("|", expand=True)

    # Extract each requested subfield
    for field in ann_fields_needed:
        if field in ANN_SUBFIELD_MAP:
            _, pos = ANN_SUBFIELD_MAP[field]

            # Special handling for AA_POS and AA_LEN (combined in one field)
            if field == "ANN[0].AA_POS":
                # Extract position part (before "/")
                if pos < ann_split.shape[1]:
                    df[field] = ann_split[pos].fillna("").str.split("/").str[0]
                else:
                    df[field] = "NA"
            elif field == "ANN[0].AA_LEN":
                # Extract length part (after "/")
                if pos < ann_split.shape[1]:
                    df[field] = ann_split[pos].fillna("").str.split("/").str[1]
                else:
                    df[field] = "NA"
            else:
                # Normal extraction
                if pos < ann_split.shape[1]:
                    df[field] = ann_split[pos].fillna("NA")
                else:
                    df[field] = "NA"

    # Drop raw ANN column
    df = df.drop(columns=["ANN"])

    return df


def parse_nmd_subfields(df: pd.DataFrame, fields: list[str]) -> pd.DataFrame:
    """
    Parse NMD subfields from raw NMD column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with raw "NMD" column containing pipe-delimited data
    fields : list[str]
        Original field list to determine which NMD subfields to extract

    Returns
    -------
    pd.DataFrame
        DataFrame with NMD subfields added and raw NMD column removed
    """
    if "NMD" not in df.columns:
        return df

    # Determine which NMD subfields are needed
    nmd_fields_needed = [f for f in fields if f.startswith("NMD[0].")]

    if not nmd_fields_needed:
        return df.drop(columns=["NMD"])

    # Take first NMD annotation
    first_nmd = df["NMD"].fillna("").str.split(",", n=1).str[0]

    # Split by pipe
    nmd_split = first_nmd.str.split("|", expand=True)

    # Extract each requested subfield
    for field in nmd_fields_needed:
        if field in NMD_SUBFIELD_MAP:
            pos = NMD_SUBFIELD_MAP[field]
            if pos < nmd_split.shape[1]:
                df[field] = nmd_split[pos].fillna("NA")
            else:
                df[field] = "NA"

    # Drop raw NMD column
    df = df.drop(columns=["NMD"])

    return df


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
        Configuration dictionary.
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
    temp_fd, temp_output = tempfile.mkstemp(suffix=".tsv", prefix="bcftools_extract_")
    os.close(temp_fd)

    try:
        # Run bcftools query
        cmd = ["bcftools", "query", "-u", "-f", format_string, variant_file]

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
        df = pd.read_csv(
            temp_output,
            sep="\t",
            header=None,
            dtype=str,
            na_values=["."],
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

        # Parse ANN subfields if present
        df = parse_ann_subfields(df, field_list)

        # Parse NMD subfields if present
        df = parse_nmd_subfields(df, field_list)

        # Normalize missing values: replace empty strings and "." with "NA"
        df = df.fillna("NA")
        df = df.replace("", "NA")
        df = df.replace(".", "NA")

        # Write output TSV
        logger.debug(f"Writing output to {output_file}")

        # Determine if we need compression
        if output_file.endswith(".gz"):
            with gzip.open(output_file, "wt", compresslevel=1) as f:
                df.to_csv(f, sep="\t", index=False, na_rep="NA")
        else:
            df.to_csv(output_file, sep="\t", index=False, na_rep="NA")

        logger.debug(
            f"bcftools field extraction completed successfully. Output written to: {output_file}"
        )

    finally:
        # Clean up temp file
        if os.path.exists(temp_output):
            os.remove(temp_output)

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

    cmd = [
        "SnpSift",
        "extractFields",
        "-s",
        snpsift_sep,  # SnpSift subfield separator
        "-e",
        "NA",  # Replace missing values with "NA"
        variant_file,
        *field_list,
    ]

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


# For backwards compatibility, keep extract_fields as alias to SnpSift version
extract_fields = extract_fields_snpsift
