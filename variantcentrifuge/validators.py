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
import re
import sys
from pathlib import Path
from typing import Any

logger = logging.getLogger("variantcentrifuge")


def _is_bgzf(header: bytes) -> bool:
    """Check if header bytes correspond to a BGZF block."""
    if len(header) < 18:
        return False
    if header[:3] != b"\x1f\x8b\x08":
        return False
    flg = header[3]
    if not (flg & 0x04):  # FEXTRA flag
        return False
    xlen = int.from_bytes(header[10:12], byteorder="little")
    if len(header) < 12 + xlen:
        return False
    pos = 12
    end = 12 + xlen
    while pos + 4 <= end:
        si1 = header[pos]
        si2 = header[pos + 1]
        slen = int.from_bytes(header[pos + 2 : pos + 4], byteorder="little")
        if si1 == 66 and si2 == 67 and slen == 2:  # 'B', 'C'
            return True
        pos += 4 + slen
    return False


def validate_vcf_file(
    vcf_path: str | None,
    logger: logging.Logger | None = None,
    exit_on_error: bool = False,
) -> None:
    """Validate that the input VCF file exists, is BGZF/BCF compressed, and has a valid index."""
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

    # Check compression: must be BGZF or BCF
    with open(vcf_path, "rb") as f:
        header = f.read(64)

    is_bcf = header.startswith(b"BCF")
    if not is_bcf and not _is_bgzf(header):
        if header.startswith(b"\x1f\x8b"):
            msg = (
                f"VCF file '{vcf_path}' is standard gzip compressed, not BGZF compressed. "
                "bcftools requires BGZF compression for random access indexing. "
                f"Please re-compress with 'bgzip -c {vcf_path} > {vcf_path}.bgz' "
                "and index with tabix."
            )
        else:
            msg = (
                f"VCF file '{vcf_path}' is uncompressed plain text. "
                "bcftools requires BGZF compressed VCF (.vcf.gz) or BCF (.bcf) with an index. "
                f"Please compress with 'bgzip {vcf_path}' "
                f"and index with 'tabix -p vcf {vcf_path}.gz'."
            )
        if logger:
            logger.error(msg)
        if exit_on_error:
            sys.exit(1)
        raise ValueError(msg)

    # Check index: must have non-empty, readable .tbi or .csi
    candidates = [
        f"{vcf_path}.tbi",
        f"{vcf_path}.csi",
    ]
    if vcf_path.endswith(".gz"):
        candidates.append(f"{vcf_path[:-3]}.tbi")
        candidates.append(f"{vcf_path[:-3]}.csi")
    elif vcf_path.endswith(".bcf"):
        candidates.append(f"{vcf_path[:-4]}.csi")

    found_idx = None
    for cand in candidates:
        if os.path.exists(cand):
            found_idx = cand
            break

    if not found_idx:
        msg = (
            f"VCF index file not found for '{vcf_path}'. "
            f"Expected {vcf_path}.tbi or {vcf_path}.csi. "
            f"Please create an index using 'tabix -p vcf {vcf_path}' "
            f"or 'bcftools index {vcf_path}'."
        )
        if logger:
            logger.error(msg)
        if exit_on_error:
            sys.exit(1)
        raise FileNotFoundError(msg)

    if os.path.getsize(found_idx) == 0:
        msg = f"VCF index file '{found_idx}' is empty (0 bytes)."
        if logger:
            logger.error(msg)
        if exit_on_error:
            sys.exit(1)
        raise ValueError(msg)

    if not os.access(found_idx, os.R_OK):
        msg = f"VCF index file '{found_idx}' is not readable."
        if logger:
            logger.error(msg)
        if exit_on_error:
            sys.exit(1)
        raise ValueError(msg)

    import gzip

    try:
        with gzip.open(found_idx, "rb") as gz_in:
            magic = gz_in.read(4)
    except Exception as e:
        msg = f"VCF index file '{found_idx}' is corrupted or not a valid gzip/BGZF index: {e}"
        if logger:
            logger.error(msg)
        if exit_on_error:
            sys.exit(1)
        raise ValueError(msg) from e

    valid_magic = False
    if found_idx.endswith(".tbi"):
        valid_magic = magic == b"TBI\x01"
    elif found_idx.endswith(".csi"):
        valid_magic = magic in (b"CSI\x01", b"CSI\x02")
    else:
        valid_magic = magic in (b"TBI\x01", b"CSI\x01", b"CSI\x02")

    if not valid_magic:
        msg = f"VCF index file '{found_idx}' has invalid format or magic header: {magic!r}"
        if logger:
            logger.error(msg)
        if exit_on_error:
            sys.exit(1)
        raise ValueError(msg)


def validate_required_tools(
    args: argparse.Namespace,
    cfg: dict[str, Any] | None = None,
    logger: logging.Logger | None = None,
    exit_on_error: bool = False,
) -> None:
    """Validate that required external tools are available in PATH before processing."""
    import shutil

    cfg = cfg or {}
    gene_name = getattr(args, "gene_name", None) or cfg.get("gene_name") or ""
    regions_bed = getattr(args, "regions_bed", None) or cfg.get("regions_bed")
    gene_file = getattr(args, "gene_file", None) or cfg.get("gene_file")
    annotate_bed = (
        getattr(args, "annotate_bed", None)
        or cfg.get("annotate_bed")
        or cfg.get("annotate_bed_files")
    )
    uses_regions_bed_directly = (
        (not gene_name or str(gene_name).strip().lower() == "all")
        and not gene_file
        and bool(regions_bed)
    )

    required_tools = ["bcftools"]

    late_filtering = bool(
        getattr(args, "late_filtering", False) or cfg.get("late_filtering", False)
    )
    if not late_filtering:
        required_tools.append("SnpSift")

    if not uses_regions_bed_directly:
        required_tools.append("snpEff")
        for tool in ("bedtools", "sortBed"):
            if tool not in required_tools:
                required_tools.append(tool)
    elif annotate_bed:
        for tool in ("bedtools", "sortBed"):
            if tool not in required_tools:
                required_tools.append(tool)

    pca = getattr(args, "pca", None) or cfg.get("pca") or cfg.get("pca_tool")
    if pca == "akt" and "akt" not in required_tools:
        required_tools.append("akt")

    missing = [tool for tool in required_tools if not shutil.which(tool)]
    if missing:
        msg = f"Required external tool(s) not found in PATH: {', '.join(missing)}"
        if logger:
            logger.error(msg)
        if exit_on_error:
            sys.stderr.write(f"{msg}\n")
            sys.exit(1)
        raise RuntimeError(msg)


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


_RESERVED_COLUMNS = {"CHROM", "POS", "ID", "REF", "ALT", "GT", "GENE", "Custom_Annotation"}


def validate_annotate_bed_specs(specs: list[str], fields: str | None = None) -> None:
    """Validate ``--annotate-bed`` entries: NAME unique, not reserved, BED readable."""
    from .feature_loaders import split_bed_spec

    # extracted columns keep the last dotted segment, e.g. ANN[0].IMPACT -> IMPACT
    extracted = {f.split(".")[-1] for f in re.split(r"[\s,]+", fields or "") if f}
    seen: set[str] = set()
    for spec in specs or []:
        name, path = split_bed_spec(spec)
        if name is not None:
            if name in seen or name in _RESERVED_COLUMNS or name in extracted:
                raise ValueError(
                    f"--annotate-bed name '{name}' is duplicated, reserved or an extracted field"
                )
            seen.add(name)
        if not Path(path).is_file():
            raise FileNotFoundError(f"BED file not found: {path}")


def validate_cli_args(
    args: argparse.Namespace,
    cfg: dict[str, Any],
    parser: argparse.ArgumentParser | None = None,
) -> None:
    """Validate cross-argument dependencies and file existence."""
    fields_file = getattr(args, "fields_file", None) or cfg.get("fields_file")
    if fields_file and not getattr(args, "fields", None) and not Path(fields_file).is_file():
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

    bed_specs = getattr(args, "annotate_bed", None) or cfg.get("annotate_bed_files") or []
    try:
        validate_annotate_bed_specs(bed_specs, getattr(args, "fields", None) or cfg.get("fields"))
    except (ValueError, FileNotFoundError) as exc:
        if parser is not None:
            parser.error(str(exc))
        raise

    phenotype_file = getattr(args, "phenotype_file", None) or cfg.get("phenotype_file")
    sample_col = getattr(args, "phenotype_sample_column", None) or cfg.get(
        "phenotype_sample_column"
    )
    val_col = getattr(args, "phenotype_value_column", None) or cfg.get("phenotype_value_column")
    if phenotype_file and (sample_col is None or val_col is None):
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

    if (getattr(args, "pseudonymize", False) or cfg.get("pseudonymize")) and not (
        getattr(args, "pseudonymize_table", None) or cfg.get("pseudonymize_table")
    ):
        msg = "--pseudonymize-table is required when using --pseudonymize"
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)
