# variantcentrifuge/validators.py

"""
Validation module for variantcentrifuge.

This module provides functions to validate:
- VCF files (existence, non-empty)
- Phenotype files (existence, non-empty, required columns)
- Mandatory parameters (reference, filters, fields)

These validations ensure that all critical inputs and parameters
are provided correctly before proceeding with the analysis.
"""

import sys
import os

def validate_vcf_file(vcf_path, logger):
    """
    Validate that the input VCF file exists, is non-empty, and is readable.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file to validate.
    logger : logging.Logger
        Logger instance for logging errors and debug information.

    Raises
    ------
    SystemExit
        If the VCF file is missing or empty.
    """
    if not vcf_path or not os.path.exists(vcf_path):
        logger.error(f"VCF file not found: {vcf_path}")
        sys.exit(1)
    if os.path.getsize(vcf_path) == 0:
        logger.error(f"VCF file {vcf_path} is empty.")
        sys.exit(1)


def validate_phenotype_file(phenotype_file, sample_col, value_col, logger):
    """
    Validate phenotype file presence, non-empty status, and required columns.

    Parameters
    ----------
    phenotype_file : str
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
        logger.error(f"Phenotype file not found: {phenotype_file}")
        sys.exit(1)
    if os.path.getsize(phenotype_file) == 0:
        logger.error(f"Phenotype file {phenotype_file} is empty.")
        sys.exit(1)

    with open(phenotype_file, "r", encoding="utf-8") as pf:
        header = pf.readline().strip()
        if not header:
            logger.error(f"Phenotype file {phenotype_file} has no header.")
            sys.exit(1)
        columns = header.split("\t") if "\t" in header else header.split(",")
        if sample_col not in columns:
            logger.error(
                f"Phenotype sample column '{sample_col}' not found in {phenotype_file}."
            )
            sys.exit(1)
        if value_col not in columns:
            logger.error(
                f"Phenotype value column '{value_col}' not found in {phenotype_file}."
            )
            sys.exit(1)
        data_line = pf.readline()
        if not data_line.strip():
            logger.error(
                f"Phenotype file {phenotype_file} contains only a header and no data."
            )
            sys.exit(1)


def validate_mandatory_parameters(reference, filters, fields):
    """
    Validate that mandatory parameters (reference, filters, fields) are provided.

    Parameters
    ----------
    reference : str
        Reference database parameter.
    filters : str
        Filters for variant extraction.
    fields : str
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
        sys.stderr.write(
            "No filters provided. Provide via --filters or in config.\n"
        )
        sys.exit(1)

    if not fields:
        sys.stderr.write(
            "No fields to extract provided. Provide via --fields or in config.\n"
        )
        sys.exit(1)
