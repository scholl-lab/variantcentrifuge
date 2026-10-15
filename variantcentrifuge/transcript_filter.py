"""Transcript-level filtering for SnpEff ANN annotations."""

from __future__ import annotations

import gzip
import logging
import os
import subprocess
import tempfile
from contextlib import ExitStack
from pathlib import Path
from typing import TextIO

logger = logging.getLogger(__name__)

ANN_FEATURE_ID_INDEX = 6


def _base_transcript_id(transcript_id: str) -> str:
    """Return transcript ID stripped of any RefSeq/Ensembl version suffix."""
    return transcript_id.strip().split(".")[0]


def _format_id_preview(items: list[str], limit: int = 10) -> str:
    if not items:
        return "none"
    if len(items) <= limit:
        return ", ".join(items)
    return f"{', '.join(items[:limit])}, ... (+{len(items) - limit} more)"


def load_transcript_ids(
    transcript_list: str | None = None,
    transcript_file: str | None = None,
) -> set[str]:
    """Load transcript IDs from a comma-separated list and/or one-per-line file."""
    transcript_ids: set[str] = set()

    if transcript_list:
        transcript_ids.update(item.strip() for item in transcript_list.split(",") if item.strip())

    if transcript_file:
        path = Path(transcript_file)
        if not path.exists():
            raise FileNotFoundError(f"Transcript file not found: {transcript_file}")

        with path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if line and not line.startswith("#"):
                    transcript_ids.add(line.split()[0])

    if (transcript_list or transcript_file) and not transcript_ids:
        raise ValueError("No transcript IDs were provided")

    return transcript_ids


def filter_ann_value_by_transcripts(
    ann_value: str,
    transcript_ids: set[str],
    base_ids: set[str] | None = None,
) -> str | None:
    """Return ANN entries whose Feature_ID matches the requested transcripts."""
    if not ann_value or ann_value == ".":
        return None

    if base_ids is None:
        base_ids = {_base_transcript_id(tid) for tid in transcript_ids if tid.strip()}

    kept_entries: list[str] = []
    for entry in ann_value.split(","):
        parts = entry.split("|")
        if len(parts) > ANN_FEATURE_ID_INDEX:
            feature_id = parts[ANN_FEATURE_ID_INDEX].strip()
            if feature_id and (
                feature_id in transcript_ids or _base_transcript_id(feature_id) in base_ids
            ):
                kept_entries.append(entry)

    if not kept_entries:
        return None

    return ",".join(kept_entries)


def filter_vcf_line_by_transcripts(
    line: str,
    transcript_ids: set[str],
    *,
    base_ids: set[str] | None = None,
    matched_features: set[str] | None = None,
) -> str | None:
    """Filter one VCF record line by ANN Feature_ID, emitting one VCF line per kept ANN entry."""
    if line.startswith("#"):
        return line

    columns = line.rstrip("\n").split("\t")
    if len(columns) < 8:
        return line

    info_parts = columns[7].split(";") if columns[7] else []
    ann_index: int | None = None
    kept_entries: list[str] = []

    for idx, info_part in enumerate(info_parts):
        if info_part.startswith("ANN="):
            ann_index = idx
            filtered_ann = filter_ann_value_by_transcripts(
                info_part.removeprefix("ANN="),
                transcript_ids,
                base_ids=base_ids,
            )
            if filtered_ann is None:
                return None
            kept_entries = filtered_ann.split(",")
            break

    if ann_index is None or not kept_entries:
        return None

    emitted_lines: list[str] = []
    for entry in kept_entries:
        if matched_features is not None:
            parts = entry.split("|")
            if len(parts) > ANN_FEATURE_ID_INDEX and parts[ANN_FEATURE_ID_INDEX].strip():
                matched_features.add(parts[ANN_FEATURE_ID_INDEX].strip())
        rewritten_info = list(info_parts)
        rewritten_info[ann_index] = f"ANN={entry}"
        row_cols = list(columns)
        row_cols[7] = ";".join(rewritten_info)
        emitted_lines.append("\t".join(row_cols) + "\n")

    return "".join(emitted_lines)


def _open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def filter_vcf_to_transcripts(
    input_vcf: str,
    output_vcf: str,
    transcript_ids: set[str],
    *,
    index_output: bool = False,
) -> int:
    """Write a VCF containing one record per matching ANN transcript ID."""
    input_records = 0
    kept_records = 0
    base_ids = {_base_transcript_id(tid) for tid in transcript_ids if tid.strip()}
    requested_by_base: dict[str, set[str]] = {}
    for tid in transcript_ids:
        if tid.strip():
            requested_by_base.setdefault(_base_transcript_id(tid), set()).add(tid.strip())
    matched_features: set[str] = set()

    output_path = Path(output_vcf)
    write_path = output_path
    temp_uncompressed: Path | None = None

    with ExitStack() as stack:
        input_handle = stack.enter_context(_open_text(input_vcf))
        if index_output and output_vcf.endswith(".gz"):
            fd, temp_name = tempfile.mkstemp(
                dir=output_path.parent,
                suffix=".vcf",
            )
            os.close(fd)
            temp_uncompressed = Path(temp_name)
            output_handle = stack.enter_context(temp_uncompressed.open("w", encoding="utf-8"))
            write_path = temp_uncompressed
        elif output_vcf.endswith(".gz"):
            output_handle = stack.enter_context(
                gzip.open(output_vcf, "wt", encoding="utf-8", compresslevel=1)
            )
        else:
            output_handle = stack.enter_context(open(output_vcf, "w", encoding="utf-8"))

        for line in input_handle:
            if not line.startswith("#") and line.strip():
                input_records += 1
            filtered = filter_vcf_line_by_transcripts(
                line,
                transcript_ids,
                base_ids=base_ids,
                matched_features=matched_features,
            )
            if filtered is None:
                continue
            if not filtered.startswith("#"):
                kept_records += filtered.count("\n")
            output_handle.write(filtered)

    if index_output and output_vcf.endswith(".gz"):
        try:
            subprocess.run(
                ["bcftools", "view", "-Oz", "-o", output_vcf, str(write_path)],
                check=True,
            )
            subprocess.run(["bcftools", "index", "-f", output_vcf], check=True)
        finally:
            if temp_uncompressed is not None:
                temp_uncompressed.unlink(missing_ok=True)

    matched_requested: set[str] = set()
    skewed_matches: set[str] = set()
    for feature_id in matched_features:
        for req_id in requested_by_base.get(_base_transcript_id(feature_id), set()):
            matched_requested.add(req_id)
            if req_id != feature_id:
                skewed_matches.add(f"{req_id}->{feature_id}")

    unmatched = sorted(transcript_ids - matched_requested)
    skew_note = (
        f"; version-skewed/base matches: {_format_id_preview(sorted(skewed_matches))}"
        if skewed_matches
        else ""
    )
    logger.info(
        "Transcript filtering retained %d records from %d input records for %d transcript IDs "
        "(matched %d/%d: %s; unmatched %d: %s%s)",
        kept_records,
        input_records,
        len(transcript_ids),
        len(matched_requested),
        len(transcript_ids),
        _format_id_preview(sorted(matched_requested)),
        len(unmatched),
        _format_id_preview(unmatched),
        skew_note,
    )
    return kept_records
