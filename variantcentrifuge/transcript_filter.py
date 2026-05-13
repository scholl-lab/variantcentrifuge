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
                    transcript_ids.add(line)

    if (transcript_list or transcript_file) and not transcript_ids:
        raise ValueError("No transcript IDs were provided")

    return transcript_ids


def filter_ann_value_by_transcripts(
    ann_value: str,
    transcript_ids: set[str],
) -> str | None:
    """Return ANN entries whose Feature_ID matches the requested transcripts."""
    if not ann_value or ann_value == ".":
        return None

    kept_entries: list[str] = []
    for entry in ann_value.split(","):
        parts = entry.split("|")
        if len(parts) > ANN_FEATURE_ID_INDEX and parts[ANN_FEATURE_ID_INDEX] in transcript_ids:
            kept_entries.append(entry)

    if not kept_entries:
        return None

    return ",".join(kept_entries)


def filter_vcf_line_by_transcripts(line: str, transcript_ids: set[str]) -> str | None:
    """Filter one VCF record line by ANN Feature_ID."""
    if line.startswith("#"):
        return line

    columns = line.rstrip("\n").split("\t")
    if len(columns) < 8:
        return line

    info_parts = columns[7].split(";") if columns[7] else []
    rewritten_parts: list[str] = []
    found_ann = False

    for info_part in info_parts:
        if info_part.startswith("ANN="):
            found_ann = True
            filtered_ann = filter_ann_value_by_transcripts(
                info_part.removeprefix("ANN="),
                transcript_ids,
            )
            if filtered_ann is None:
                return None
            rewritten_parts.append(f"ANN={filtered_ann}")
        else:
            rewritten_parts.append(info_part)

    if not found_ann:
        return None

    columns[7] = ";".join(rewritten_parts)
    return "\t".join(columns) + "\n"


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
    """Write a VCF containing only records with matching ANN transcript IDs."""
    kept_records = 0
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
            filtered = filter_vcf_line_by_transcripts(line, transcript_ids)
            if filtered is None:
                continue
            if not filtered.startswith("#"):
                kept_records += 1
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

    logger.info(
        "Transcript filtering retained %d records for %d transcript IDs",
        kept_records,
        len(transcript_ids),
    )
    return kept_records
