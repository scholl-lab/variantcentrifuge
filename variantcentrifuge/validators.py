# File: variantcentrifuge/validators.py
# Location: variantcentrifuge/variantcentrifuge/validators.py

"""
Validation module for variantcentrifuge.

This module provides functions to validate:
- VCF files (existence, non-empty)
- Phenotype files (existence, non-empty, required columns)
- Mandatory parameters (reference, filters, fields)
- IGV local FASTA and related files (existence, readability)

These validations ensure that all critical inputs and parameters
are provided correctly before proceeding with the analysis.
"""

import logging
import os
import sys
from typing import Optional

logger = logging.getLogger("variantcentrifuge")


def validate_vcf_file(vcf_path: Optional[str], logger: logging.Logger) -> None:
    """
    Validate that the input VCF file exists, is non-empty, and is readable.

    Parameters
    ----------
    vcf_path : str or None
        Path to the VCF file to validate.
    logger : logging.Logger
        Logger instance for logging errors and debug information.

    Raises
    ------
    SystemExit
        If the VCF file is missing or empty.
    """
    if not vcf_path or not os.path.exists(vcf_path):
        logger.error("VCF file not found: %s", vcf_path)
        sys.exit(1)
    if os.path.getsize(vcf_path) == 0:
        logger.error("VCF file %s is empty.", vcf_path)
        sys.exit(1)


def validate_phenotype_file(
    phenotype_file: Optional[str],
    sample_col: str,
    value_col: str,
    logger: logging.Logger,
) -> None:
    """
    Validate phenotype file presence, non-empty status, and required columns.

    Parameters
    ----------
    phenotype_file : str or None
        Path to the phenotype file.
    sample_col : str
        The expected column name for samples in the phenotype file.
    value_col : str
        The expected column name for phenotype values in the phenotype file.
    logger : logging.Logger
        Logger instance for logging errors and debug information.

    Raises
    ------
    SystemExit
        If the phenotype file is missing, empty, lacks the required columns,
        or contains no data rows.
    """
    if not phenotype_file:
        return

    if not os.path.exists(phenotype_file):
        logger.error("Phenotype file not found: %s", phenotype_file)
        sys.exit(1)
    if os.path.getsize(phenotype_file) == 0:
        logger.error("Phenotype file %s is empty.", phenotype_file)
        sys.exit(1)

    with open(phenotype_file, "r", encoding="utf-8") as pf:
        header = pf.readline().strip()
        if not header:
            logger.error("Phenotype file %s has no header.", phenotype_file)
            sys.exit(1)

        # Determine delimiter by checking for tabs or commas
        if "\t" in header:
            columns = header.split("\t")
        else:
            columns = header.split(",")

        if sample_col not in columns:
            logger.error(
                "Phenotype sample column '%s' not found in %s.",
                sample_col,
                phenotype_file,
            )
            sys.exit(1)
        if value_col not in columns:
            logger.error(
                "Phenotype value column '%s' not found in %s.",
                value_col,
                phenotype_file,
            )
            sys.exit(1)

        data_line = pf.readline()
        if not data_line.strip():
            logger.error("Phenotype file %s contains only a header and no data.", phenotype_file)
            sys.exit(1)


def validate_mandatory_parameters(
    reference: Optional[str], filters: Optional[str], fields: Optional[str]
) -> None:
    """
    Validate that mandatory parameters (reference, filters, fields) are provided.

    Parameters
    ----------
    reference : str or None
        Reference database parameter.
    filters : str or None
        Filters for variant extraction.
    fields : str or None
        Fields to extract from the VCF.

    Raises
    ------
    SystemExit
        If any of the mandatory parameters are missing.
    """
    if not reference:
        sys.stderr.write(
            "A reference database must be specified either via --reference or "
            "in the configuration file.\n"
        )
        sys.exit(1)

    if not filters:
        sys.stderr.write("No filters provided. Provide via --filters or in config.\n")
        sys.exit(1)

    if not fields:
        sys.stderr.write("No fields to extract provided. Provide via --fields or in config.\n")
        sys.exit(1)


# MODIFIED: Start of local IGV FASTA feature
def validate_igv_files(local_fasta: Optional[str], ideogram: Optional[str]) -> None:
    """Validate the existence and readability of files for IGV integration.

    Parameters
    ----------
    local_fasta : str or None
        Path to the local FASTA file. If provided, its index file is expected
        to be co-located with the same name plus '.fai' extension (e.g., genome.fa.fai).
    ideogram : str or None
        Path to the ideogram file.

    Raises
    ------
    SystemExit
        If any of the specified files do not exist or are not readable.
    """
    if local_fasta:
        if not os.path.exists(local_fasta):
            sys.stderr.write(f"Local FASTA file not found: {local_fasta}\n")
            sys.exit(1)
        if os.path.getsize(local_fasta) == 0:
            sys.stderr.write(f"Local FASTA file is empty: {local_fasta}\n")
            sys.exit(1)

        # Check for FASTA index with standard naming convention (e.g., genome.fa.fai)
        index_path = f"{local_fasta}.fai"
        index_description = "FASTA index"
        if not os.path.exists(index_path):
            sys.stderr.write(f"{index_description} file not found: {index_path}\n")
            # Provide helpful suggestion on how to create the index file
            sys.stderr.write(f"Please index your FASTA file with 'samtools faidx {local_fasta}'\n")
            sys.exit(1)

        if os.path.getsize(index_path) == 0:
            sys.stderr.write(f"{index_description} file is empty: {index_path}\n")
            sys.exit(1)

    if ideogram:
        if not os.path.exists(ideogram):
            sys.stderr.write(f"Ideogram file not found: {ideogram}\n")
            sys.exit(1)
        if os.path.getsize(ideogram) == 0:
            sys.stderr.write(f"Ideogram file is empty: {ideogram}\n")
            sys.exit(1)


# MODIFIED: End of local IGV FASTA feature
