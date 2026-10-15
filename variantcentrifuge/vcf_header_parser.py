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
from dataclasses import dataclass, field

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
TUMOR_SAMPLE_RE = re.compile(r"^##tumor_sample=(?P<sample>\S+)\s*$")
NORMAL_SAMPLE_RE = re.compile(r"^##normal_sample=(?P<sample>\S+)\s*$")


@dataclass
class AnnotationField:
    """A single annotation field parsed from a VCF header."""

    id: str
    number: str
    type: str
    description: str
    source: str  # "INFO" or "FORMAT"


@dataclass
class VcfSampleHeader:
    """Sample names and somatic sample roles parsed from a VCF header."""

    samples: list[str]
    tumor_sample: str | None = None
    normal_sample: str | None = None
    tumor_samples: list[str] = field(default_factory=list)


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


def parse_vcf_sample_header(vcf_path: str, header_text: str | None = None) -> VcfSampleHeader:
    """Parse ``#CHROM`` samples and ``##tumor_sample``/``##normal_sample`` from a VCF header."""
    text = _read_header(vcf_path) if header_text is None else header_text
    samples: list[str] = []
    tumor_samples: list[str] = []
    normal_sample: str | None = None
    for line in text.splitlines():
        if m := TUMOR_SAMPLE_RE.match(line):
            tumor_samples.append(m.group("sample"))
        elif m := NORMAL_SAMPLE_RE.match(line):
            normal_sample = m.group("sample")
        elif line.startswith("#CHROM"):
            cols = line.rstrip("\r\n").split("\t")
            samples = cols[9:] if len(cols) > 9 else []
    return VcfSampleHeader(
        samples=samples,
        tumor_sample=tumor_samples[-1] if tumor_samples else None,
        normal_sample=normal_sample,
        tumor_samples=tumor_samples,
    )


def resolve_tumor_normal_indices(
    vcf_path: str | None = None,
    *,
    cli_tumor_idx: int | str | None = None,
    cli_normal_idx: int | str | None = None,
    cfg_tumor_idx: int | str | None = None,
    cfg_normal_idx: int | str | None = None,
    header_text: str | None = None,
) -> tuple[str, str]:
    """Resolve ``(tumor_idx, normal_idx)`` from VCF header, CLI, or config."""
    hdr = VcfSampleHeader(samples=[])
    if header_text is not None:
        hdr = parse_vcf_sample_header(vcf_path or "", header_text=header_text)
    elif vcf_path:
        hdr = parse_vcf_sample_header(vcf_path)

    explicit_t = (
        int(cli_tumor_idx)
        if cli_tumor_idx is not None
        else (int(cfg_tumor_idx) if cfg_tumor_idx is not None else None)
    )
    explicit_n = (
        int(cli_normal_idx)
        if cli_normal_idx is not None
        else (int(cfg_normal_idx) if cfg_normal_idx is not None else None)
    )

    header_t: int | None = None
    tumor_list = hdr.tumor_samples or ([hdr.tumor_sample] if hdr.tumor_sample is not None else [])
    if tumor_list:
        header_ts: list[int] = []
        for ts in tumor_list:
            if ts not in hdr.samples:
                raise ValueError(f"##tumor_sample={ts!r} not found in #CHROM samples {hdr.samples}")
            header_ts.append(hdr.samples.index(ts))
        if len(header_ts) > 1:
            if explicit_t is None:
                raise ValueError(
                    f"Multiple ##tumor_sample headers found ({tumor_list}); "
                    "specify --tumor-sample-index explicitly"
                )
            if explicit_t not in header_ts:
                raise ValueError(
                    f"tumor_sample_index={explicit_t} conflicts with "
                    f"##tumor_sample headers {tumor_list} at indices {header_ts}"
                )
            header_t = explicit_t
        else:
            header_t = header_ts[0]
            if explicit_t is not None and explicit_t != header_t:
                raise ValueError(
                    f"tumor_sample_index={explicit_t} conflicts with "
                    f"##tumor_sample={tumor_list[0]!r} at index {header_t}"
                )

    header_n: int | None = None
    if hdr.normal_sample is not None:
        if hdr.normal_sample not in hdr.samples:
            raise ValueError(
                f"##normal_sample={hdr.normal_sample!r} not found in #CHROM samples {hdr.samples}"
            )
        header_n = hdr.samples.index(hdr.normal_sample)
        if explicit_n is not None and explicit_n != header_n:
            raise ValueError(
                f"normal_sample_index={explicit_n} conflicts with "
                f"##normal_sample={hdr.normal_sample!r} at index {header_n}"
            )

    if explicit_t is not None and (
        explicit_t < 0 or (hdr.samples and explicit_t >= len(hdr.samples))
    ):
        raise ValueError(
            f"tumor_sample_index={explicit_t} out of range for VCF samples {hdr.samples}"
        )
    if explicit_n is not None and (
        explicit_n < 0 or (hdr.samples and explicit_n >= len(hdr.samples))
    ):
        raise ValueError(
            f"normal_sample_index={explicit_n} out of range for VCF samples {hdr.samples}"
        )

    if header_t is not None:
        t_idx = header_t
    elif explicit_t is not None:
        t_idx = explicit_t
    elif header_n in (0, 1) and len(hdr.samples) == 2:
        t_idx = 1 - header_n
    elif len(hdr.samples) == 1:
        t_idx = 0
    else:
        t_idx = 1

    if header_n is not None:
        n_idx = header_n
    elif explicit_n is not None:
        n_idx = explicit_n
    elif t_idx in (0, 1) and len(hdr.samples) == 2:
        n_idx = 1 - t_idx
    else:
        n_idx = 0

    if (
        (header_t is not None or explicit_t is not None)
        and (header_n is not None or explicit_n is not None)
        and t_idx == n_idx
    ):
        raise ValueError(f"Tumor and normal sample indices cannot be identical ({t_idx})")

    return str(t_idx), str(n_idx)


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
