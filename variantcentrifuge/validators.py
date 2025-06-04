# File: variantcentrifuge/validators.py
# Location: variantcentrifuge/variantcentrifuge/validators.py

"""
Validation module for variantcentrifuge.

This module provides functions to validate:
- VCF files (existence, non-empty)
- Phenotype files (existence, non-empty, required columns)
- Mandatory parameters (reference, filters, fields)

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
