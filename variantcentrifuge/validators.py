"""Validation module for variantcentrifuge.

This module provides functions to validate:
- VCF files (existence, non-empty)
- Phenotype files (existence, non-empty, required columns)
- Mandatory parameters (reference, filters, fields)
- IGV local FASTA and related files (existence, readability)
- CLI cross-argument dependencies

These validations ensure that all critical inputs and parameters
are provided correctly before proceeding with the analysis.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Any

logger = logging.getLogger("variantcentrifuge")


def validate_vcf_file(
    vcf_path: str | None,
    logger: logging.Logger | None = None,
    exit_on_error: bool = False,
) -> None:
    """Validate that the input VCF file exists, is non-empty, and is readable."""
    if not vcf_path or not os.path.exists(vcf_path):
        msg = f"VCF file not found: {vcf_path}"
        if logger:
            logger.error("VCF file not found: %s", vcf_path)
        if exit_on_error:
            sys.exit(1)
        raise FileNotFoundError(msg) if vcf_path else ValueError(msg)

    if os.path.getsize(vcf_path) == 0:
        msg = f"VCF file {vcf_path} is empty."
        if logger:
            logger.error("VCF file %s is empty.", vcf_path)
        if exit_on_error:
            sys.exit(1)
        raise ValueError(msg)


def validate_phenotype_file(
    phenotype_file: str | None,
    sample_col: str,
    value_col: str,
    logger: logging.Logger | None = None,
    exit_on_error: bool = False,
) -> None:
    """Validate phenotype file presence, non-empty status, and required columns."""
    if not phenotype_file:
        return

    if not os.path.exists(phenotype_file):
        msg = f"Phenotype file not found: {phenotype_file}"
        if logger:
            logger.error("Phenotype file not found: %s", phenotype_file)
        if exit_on_error:
            sys.exit(1)
        raise FileNotFoundError(msg)

    if os.path.getsize(phenotype_file) == 0:
        msg = f"Phenotype file {phenotype_file} is empty."
        if logger:
            logger.error("Phenotype file %s is empty.", phenotype_file)
        if exit_on_error:
            sys.exit(1)
        raise ValueError(msg)

    with open(phenotype_file, encoding="utf-8") as pf:
        header = pf.readline().strip()
        if not header:
            msg = f"Phenotype file {phenotype_file} has no header."
            if logger:
                logger.error("Phenotype file %s has no header.", phenotype_file)
            if exit_on_error:
                sys.exit(1)
            raise ValueError(msg)

        columns = header.split("\t") if "\t" in header else header.split(",")

        if sample_col not in columns:
            msg = f"Phenotype sample column '{sample_col}' not found in {phenotype_file}."
            if logger:
                logger.error(
                    "Phenotype sample column '%s' not found in %s.",
                    sample_col,
                    phenotype_file,
                )
            if exit_on_error:
                sys.exit(1)
            raise ValueError(msg)

        if value_col not in columns:
            msg = f"Phenotype value column '{value_col}' not found in {phenotype_file}."
            if logger:
                logger.error(
                    "Phenotype value column '%s' not found in %s.",
                    value_col,
                    phenotype_file,
                )
            if exit_on_error:
                sys.exit(1)
            raise ValueError(msg)

        data_line = pf.readline()
        if not data_line.strip():
            msg = f"Phenotype file {phenotype_file} contains only a header and no data."
            if logger:
                logger.error(
                    "Phenotype file %s contains only a header and no data.", phenotype_file
                )
            if exit_on_error:
                sys.exit(1)
            raise ValueError(msg)


def validate_mandatory_parameters(
    reference: str | None,
    filters: str | None,
    fields: str | None,
    exit_on_error: bool = False,
) -> None:
    """Validate that mandatory parameters (reference, filters, fields) are provided."""
    if not reference:
        msg = (
            "A reference database must be specified either via --reference or "
            "in the configuration file."
        )
        if exit_on_error:
            sys.stderr.write(f"{msg}\n")
            sys.exit(1)
        raise ValueError(msg)

    if not filters:
        msg = "No filters provided. Provide via --filters or in config."
        if exit_on_error:
            sys.stderr.write(f"{msg}\n")
            sys.exit(1)
        raise ValueError(msg)

    if not fields:
        msg = "No fields to extract provided. Provide via --fields or in config."
        if exit_on_error:
            sys.stderr.write(f"{msg}\n")
            sys.exit(1)
        raise ValueError(msg)


def validate_igv_files(
    local_fasta: str | None,
    ideogram: str | None,
    exit_on_error: bool = False,
) -> None:
    """Validate the existence and readability of files for IGV integration."""
    if local_fasta:
        if not os.path.exists(local_fasta):
            msg = f"Local FASTA file not found: {local_fasta}"
            if exit_on_error:
                sys.stderr.write(f"{msg}\n")
                sys.exit(1)
            raise FileNotFoundError(msg)

        if os.path.getsize(local_fasta) == 0:
            msg = f"Local FASTA file is empty: {local_fasta}"
            if exit_on_error:
                sys.stderr.write(f"{msg}\n")
                sys.exit(1)
            raise ValueError(msg)

        index_path = f"{local_fasta}.fai"
        if not os.path.exists(index_path):
            msg = (
                f"FASTA index file not found: {index_path}\n"
                f"Please index your FASTA file with 'samtools faidx {local_fasta}'"
            )
            if exit_on_error:
                sys.stderr.write(f"{msg}\n")
                sys.exit(1)
            raise FileNotFoundError(msg)

        if os.path.getsize(index_path) == 0:
            msg = f"FASTA index file is empty: {index_path}"
            if exit_on_error:
                sys.stderr.write(f"{msg}\n")
                sys.exit(1)
            raise ValueError(msg)

    if ideogram:
        if not os.path.exists(ideogram):
            msg = f"Ideogram file not found: {ideogram}"
            if exit_on_error:
                sys.stderr.write(f"{msg}\n")
                sys.exit(1)
            raise FileNotFoundError(msg)

        if os.path.getsize(ideogram) == 0:
            msg = f"Ideogram file is empty: {ideogram}"
            if exit_on_error:
                sys.stderr.write(f"{msg}\n")
                sys.exit(1)
            raise ValueError(msg)


def validate_cli_args(
    args: argparse.Namespace,
    cfg: dict[str, Any],
    parser: argparse.ArgumentParser | None = None,
) -> None:
    """Validate cross-argument dependencies and file existence."""
    fields_file = getattr(args, "fields_file", None) or cfg.get("fields_file")
    if fields_file and not Path(fields_file).is_file():
        msg = f"Fields file not found: {fields_file}"
        if parser is not None:
            parser.error(msg)
        raise FileNotFoundError(msg)

    tx_file = getattr(args, "transcript_file", None) or cfg.get("transcript_file")
    if tx_file and not Path(tx_file).is_file():
        msg = f"Transcript file not found: {tx_file}"
        if parser is not None:
            parser.error(msg)
        raise FileNotFoundError(msg)

    phenotype_file = getattr(args, "phenotype_file", None) or cfg.get("phenotype_file")
    if phenotype_file and (
        getattr(args, "phenotype_sample_column", None) is None
        or getattr(args, "phenotype_value_column", None) is None
    ):
        msg = (
            "If a phenotype file is provided, both --phenotype-sample-column "
            "and --phenotype-value-column must be specified."
        )
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    if getattr(args, "resume", False) and not getattr(args, "enable_checkpoint", False):
        msg = "--resume requires --enable-checkpoint to be set"
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    if getattr(args, "resume_from", None) and not getattr(args, "enable_checkpoint", False):
        msg = "--resume-from requires --enable-checkpoint to be set"
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    if getattr(args, "interactive_resume", False) and not getattr(args, "enable_checkpoint", False):
        msg = "--interactive-resume requires --enable-checkpoint to be set"
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    resume_options = [
        getattr(args, "resume", False),
        getattr(args, "resume_from", None),
        getattr(args, "interactive_resume", False),
    ]
    if sum(bool(x) for x in resume_options) > 1:
        msg = (
            "Cannot use multiple resume options simultaneously: "
            "--resume, --resume-from, --interactive-resume"
        )
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    if (getattr(args, "pseudonymize", False) or cfg.get("pseudonymize")) and not getattr(
        args, "pseudonymize_table", None
    ):
        msg = "--pseudonymize-table is required when using --pseudonymize"
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)
