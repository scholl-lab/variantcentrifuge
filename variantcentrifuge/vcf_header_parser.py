"""VCF header parsing for annotation field discovery.

Parses ##INFO and ##FORMAT lines from VCF headers to let users discover
available fields before configuring extraction.  Uses ``bcftools view
--header-only`` so that compressed / indexed VCFs are handled transparently.
"""

from __future__ import annotations

import json
import logging
import re
import subprocess
from dataclasses import dataclass

logger = logging.getLogger("variantcentrifuge")

# Regex for standard VCF meta-information lines (INFO / FORMAT).
# Handles optional extra attributes (Source=, Version=) after Description.
HEADER_RE = re.compile(
    r"##(?P<source>INFO|FORMAT)=<"
    r"ID=(?P<id>[^,]+),"
    r"Number=(?P<number>[^,]+),"
    r"Type=(?P<type>[^,]+),"
    r'Description="(?P<desc>[^"]*)"'
)


@dataclass
class AnnotationField:
    """A single annotation field parsed from a VCF header."""

    id: str
    number: str
    type: str
    description: str
    source: str  # "INFO" or "FORMAT"


def parse_vcf_header(vcf_path: str) -> list[AnnotationField]:
    """Parse INFO and FORMAT fields from a VCF header.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file (plain or compressed).

    Returns
    -------
    list[AnnotationField]
        Parsed fields sorted by (source, id).
    """
    header_text = _read_header(vcf_path)
    if not header_text:
        logger.warning("Empty or missing VCF header for %s", vcf_path)
        return []

    fields: list[AnnotationField] = []
    for line in header_text.splitlines():
        m = HEADER_RE.match(line)
        if m:
            fields.append(
                AnnotationField(
                    id=m.group("id"),
                    number=m.group("number"),
                    type=m.group("type"),
                    description=m.group("desc"),
                    source=m.group("source"),
                )
            )

    fields.sort(key=lambda f: (f.source, f.id.lower()))
    return fields


def _read_header(vcf_path: str) -> str:
    """Run ``bcftools view --header-only`` and return stdout."""
    cmd = ["bcftools", "view", "--header-only", vcf_path]
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout
    except FileNotFoundError:
        logger.error("bcftools not found. Ensure it is installed and on PATH.")
        raise
    except subprocess.CalledProcessError as exc:
        logger.error("bcftools failed (exit %d): %s", exc.returncode, exc.stderr.strip())
        raise


def format_fields_table(
    fields: list[AnnotationField],
    *,
    pattern: str | None = None,
    info_only: bool = False,
    format_only: bool = False,
) -> str:
    """Format parsed fields as a human-readable table.

    Parameters
    ----------
    fields : list[AnnotationField]
        Fields to display.
    pattern : str, optional
        Case-insensitive substring filter on field ID.
    info_only : bool
        Show only INFO fields.
    format_only : bool
        Show only FORMAT fields.

    Returns
    -------
    str
        Formatted table text.
    """
    filtered = _filter_fields(fields, pattern=pattern, info_only=info_only, format_only=format_only)

    sections: list[str] = []
    for source in ("INFO", "FORMAT"):
        group = [f for f in filtered if f.source == source]
        if not group:
            continue
        lines = [f"{source} fields ({len(group)} total):"]
        lines.append(f"  {'ID':<40s} {'Type':<10s} {'Number':<8s} Description")
        for f in group:
            desc = f.description[:60] + "..." if len(f.description) > 60 else f.description
            lines.append(f"  {f.id:<40s} {f.type:<10s} {f.number:<8s} {desc}")
        sections.append("\n".join(lines))

    return "\n\n".join(sections) if sections else "No matching fields found."


def format_fields_json(
    fields: list[AnnotationField],
    *,
    pattern: str | None = None,
    info_only: bool = False,
    format_only: bool = False,
) -> str:
    """Format parsed fields as JSON.

    Parameters
    ----------
    fields : list[AnnotationField]
        Fields to display.
    pattern : str, optional
        Case-insensitive substring filter.
    info_only : bool
        Show only INFO fields.
    format_only : bool
        Show only FORMAT fields.

    Returns
    -------
    str
        JSON string.
    """
    filtered = _filter_fields(fields, pattern=pattern, info_only=info_only, format_only=format_only)
    data = [
        {
            "id": f.id,
            "source": f.source,
            "type": f.type,
            "number": f.number,
            "description": f.description,
        }
        for f in filtered
    ]
    return json.dumps(data, indent=2)


def _filter_fields(
    fields: list[AnnotationField],
    *,
    pattern: str | None = None,
    info_only: bool = False,
    format_only: bool = False,
) -> list[AnnotationField]:
    """Apply source and pattern filters."""
    result = fields
    if info_only:
        result = [f for f in result if f.source == "INFO"]
    elif format_only:
        result = [f for f in result if f.source == "FORMAT"]
    if pattern:
        pat = pattern.lower()
        result = [f for f in result if pat in f.id.lower()]
    return result


def validate_requested_fields(
    requested_fields: str,
    available: list[AnnotationField],
) -> list[str]:
    """Check which requested extraction fields are missing from the VCF header.

    Handles SnpEff compound fields (``ANN[*].GENE`` → checks ``ANN`` exists)
    and FORMAT genotype fields (``GEN[*].DP`` → checks ``DP`` in FORMAT).

    Parameters
    ----------
    requested_fields : str
        Space-separated field names as passed to SnpSift extractFields.
    available : list[AnnotationField]
        Parsed fields from the VCF header.

    Returns
    -------
    list[str]
        Field names that could not be validated against the header.
        Empty list means all fields validated successfully.
    """
    info_ids = {f.id for f in available if f.source == "INFO"}
    format_ids = {f.id for f in available if f.source == "FORMAT"}

    # Standard VCF fixed columns are always valid
    fixed_cols = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"}

    missing: list[str] = []
    for raw_field in requested_fields.strip().split():
        field_name = raw_field.strip()
        if not field_name:
            continue

        # Fixed VCF columns
        if field_name in fixed_cols:
            continue

        # ANN[*].X / EFF[*].X / LOF[*].X / NMD[*].X — check parent INFO field
        ann_match = re.match(r"(ANN|EFF|LOF|NMD)\[.*\]\.\w+", field_name)
        if ann_match:
            parent = ann_match.group(1)
            if parent in info_ids:
                continue
            missing.append(field_name)
            continue

        # GEN[*].X — check FORMAT field after the dot
        gen_match = re.match(r"GEN\[.*\]\.(\w+)", field_name)
        if gen_match:
            fmt_field = gen_match.group(1)
            if fmt_field in format_ids or fmt_field == "GT":
                continue
            missing.append(field_name)
            continue

        # Direct INFO field
        if field_name in info_ids:
            continue

        # Direct FORMAT field
        if field_name in format_ids:
            continue

        missing.append(field_name)

    return missing
