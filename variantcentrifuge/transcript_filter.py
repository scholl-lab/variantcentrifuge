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

ANN_GENE_NAME_INDEX = 3
ANN_FEATURE_ID_INDEX = 6
DROPPED_GENES_TSV = "transcript_filter_dropped_genes.tsv"


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
    alt: str | None = None,
    keep_genes: set[str] | None = None,
) -> str | None:
    """Return ANN entries whose Feature_ID matches the requested transcripts (or gene is kept)."""
    if not ann_value or ann_value == ".":
        return None

    if base_ids is None:
        base_ids = {_base_transcript_id(tid) for tid in transcript_ids if tid.strip()}

    entries = [e for e in ann_value.split(",") if e]
    if alt and "," not in alt:
        entries = [e for e in entries if e.split("|", 1)[0] == alt]
        if not entries:
            return None

    parsed: list[tuple[str, str, str, str]] = []
    for entry in entries:
        parts = entry.split("|")
        if len(parts) > ANN_FEATURE_ID_INDEX:
            feature_id = parts[ANN_FEATURE_ID_INDEX].strip()
            if feature_id:
                gene = parts[ANN_GENE_NAME_INDEX]
                parsed.append((entry, feature_id, _base_transcript_id(feature_id), gene))

    exact_bases = {base for _, fid, base, _ in parsed if fid in transcript_ids}
    kept_entries = [
        entry
        for entry, fid, base, gene in parsed
        if fid in transcript_ids
        or (base in base_ids and base not in exact_bases)
        or (keep_genes is not None and gene in keep_genes)
    ]

    if not kept_entries:
        return None

    return ",".join(kept_entries)


def filter_vcf_line_by_transcripts(
    line: str,
    transcript_ids: set[str],
    *,
    base_ids: set[str] | None = None,
    matched_features: set[str] | None = None,
    genes_before: set[str] | None = None,
    genes_after: set[str] | None = None,
    keep_genes: set[str] | None = None,
) -> str | None:
    """Filter one VCF record line by ANN Feature_ID, emitting one VCF line per kept ANN entry."""
    if line.startswith("#"):
        return line

    columns = line.rstrip("\n").split("\t")
    if len(columns) < 8:
        return line

    info_parts = columns[7].split(";") if columns[7] else []
    ann_index: int | None = None
    all_ann_entries: list[str] = []
    kept_entries: list[str] = []

    for idx, info_part in enumerate(info_parts):
        if info_part.startswith("ANN="):
            ann_index = idx
            raw_ann = info_part.removeprefix("ANN=")
            all_ann_entries = [e for e in raw_ann.split(",") if e]
            if genes_before is not None:
                genes_before.update(_entry_genes(all_ann_entries))
            filtered_ann = filter_ann_value_by_transcripts(
                raw_ann,
                transcript_ids,
                base_ids=base_ids,
                alt=columns[4],
                keep_genes=keep_genes,
            )
            if filtered_ann is None:
                return None
            kept_entries = filtered_ann.split(",")
            break

    if ann_index is None or not kept_entries:
        return None
    if genes_after is not None:
        genes_after.update(_entry_genes(kept_entries))

    from .vcf_eff_one_per_line import _eligible_lof_nmd_gene_keys, _filter_gene_scoped_info_value

    emitted_lines: list[str] = []
    for entry in kept_entries:
        if matched_features is not None:
            parts = entry.split("|")
            if len(parts) > ANN_FEATURE_ID_INDEX and parts[ANN_FEATURE_ID_INDEX].strip():
                matched_features.add(parts[ANN_FEATURE_ID_INDEX].strip())
        gene_keys = _eligible_lof_nmd_gene_keys([entry], all_ann_entries)
        rewritten_info: list[str] = []
        for idx, part in enumerate(info_parts):
            if idx == ann_index:
                rewritten_info.append(f"ANN={entry}")
            elif part.startswith(("LOF=", "NMD=")):
                key, val = part.split("=", 1)
                filtered_val = _filter_gene_scoped_info_value(val, gene_keys, key)
                if filtered_val is not None:
                    rewritten_info.append(f"{key}={filtered_val}")
            else:
                rewritten_info.append(part)
        row_cols = list(columns)
        row_cols[7] = ";".join(rewritten_info)
        emitted_lines.append("\t".join(row_cols) + "\n")

    return "".join(emitted_lines)


def _entry_genes(entries: list[str]) -> set[str]:
    genes = set()
    for entry in entries:
        parts = entry.split("|")
        if len(parts) > ANN_GENE_NAME_INDEX and parts[ANN_GENE_NAME_INDEX]:
            genes.add(parts[ANN_GENE_NAME_INDEX])
    return genes


def _genes_without_matched_transcript(
    input_vcf: str, transcript_ids: set[str], base_ids: set[str]
) -> set[str]:
    """First pass for --transcript-fallback keep-gene: genes none of whose ANN entries match."""
    genes_before: set[str] = set()
    genes_after: set[str] = set()
    with _open_text(input_vcf) as handle:
        for line in handle:
            if not line.startswith("#") and line.strip():
                filter_vcf_line_by_transcripts(
                    line,
                    transcript_ids,
                    base_ids=base_ids,
                    genes_before=genes_before,
                    genes_after=genes_after,
                )
    return genes_before - genes_after


def merge_filter_stats(into: dict, other: dict) -> dict:
    """Add one filter run's stats (record counts, ID/gene sets) to an aggregate (chunks)."""
    for key, value in other.items():
        if key not in into:
            into[key] = value
        elif isinstance(value, set):
            into[key] = into[key] | value
        else:
            into[key] += value
    return into


def log_transcript_filter_summary(
    stats: dict, transcript_ids: set[str], diagnostics_dir: str | None = None
) -> None:
    """Log the transcript-filter summary and dropped genes; write the TSV to diagnostics_dir."""
    requested_by_base: dict[str, set[str]] = {}
    for tid in transcript_ids:
        if tid.strip():
            requested_by_base.setdefault(_base_transcript_id(tid), set()).add(tid.strip())
    matched_requested: set[str] = set()
    skewed_matches: set[str] = set()
    for feature_id in stats["matched_features"]:
        for req_id in requested_by_base.get(_base_transcript_id(feature_id), set()):
            matched_requested.add(req_id)
            if req_id != feature_id:
                skewed_matches.add(f"{req_id}->{feature_id}")
    n_exact = len(matched_requested & stats["matched_features"])

    unmatched = sorted(transcript_ids - matched_requested)
    skew_note = (
        f"; version-skewed/base matches: {_format_id_preview(sorted(skewed_matches))}"
        if skewed_matches
        else ""
    )
    logger.info(
        "Transcript filtering retained %d records from %d input records for %d transcript IDs "
        "(matched %d/%d (exact %d, version-skewed %d): %s; unmatched %d: %s%s); genes %d -> %d",
        stats["kept_records"],
        stats["input_records"],
        len(transcript_ids),
        len(matched_requested),
        len(transcript_ids),
        n_exact,
        len(matched_requested) - n_exact,
        _format_id_preview(sorted(matched_requested)),
        len(unmatched),
        _format_id_preview(unmatched),
        skew_note,
        len(stats["genes_before"]),
        len(stats["genes_after"]),
    )

    dropped = sorted(stats["genes_before"] - stats["genes_after"])
    if diagnostics_dir:
        out = Path(diagnostics_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / DROPPED_GENES_TSV).write_text("".join(f"{g}\n" for g in ["gene", *dropped]))
    if dropped:
        full_list = (
            f"full list: {Path(diagnostics_dir) / DROPPED_GENES_TSV}"
            if diagnostics_dir
            else "full list with --diagnostics-output"
        )
        logger.warning(
            "Transcript filtering dropped %d genes that had annotated variants but no matched "
            "transcript (%s; keep them with --transcript-fallback keep-gene): %s",
            len(dropped),
            full_list,
            _format_id_preview(dropped),
        )


def _open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def load_filter_exceptions(
    exception_file: str | None,
) -> dict[tuple[str, str, str, str], set[str]]:
    """Load 5-column hotspot FILTER exception TSV (CHROM, POS, REF, ALT, ALLOWED_FILTERS)."""
    import re

    exceptions: dict[tuple[str, str, str, str], set[str]] = {}
    if not exception_file:
        return exceptions
    with open(exception_file, encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                cols = line.split()
            if len(cols) < 5 or (
                cols[0].upper() in {"CHROM", "#CHROM"} and cols[1].upper() == "POS"
            ):
                continue
            chrom, pos, ref, alt, raw_allowed = (c.strip() for c in cols[:5])
            allowed = {tok.strip() for tok in re.split(r"[,;]", raw_allowed) if tok.strip()}
            alt_chrom = chrom.removeprefix("chr") if chrom.startswith("chr") else f"chr{chrom}"
            for ch in (chrom, alt_chrom):
                exceptions.setdefault((ch, pos, ref, alt), set()).update(allowed)
    return exceptions


def apply_filter_exceptions_to_line(
    line: str,
    exceptions: dict[tuple[str, str, str, str], set[str]],
) -> str:
    """Rescue matching hotspot VCF line to FILTER=PASS and record INFO/RAW_FILTER."""
    if not exceptions or line.startswith("#"):
        return line
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 8:
        return line
    allowed = exceptions.get((cols[0], cols[1], cols[3], cols[4]))
    if not allowed:
        return line
    raw_filter = cols[6].strip()
    if raw_filter in {"", ".", "PASS"}:
        return line
    tokens = [t.strip() for t in raw_filter.split(";") if t.strip()]
    permitted = allowed | {"PASS", "on_target", "clustered_events"}
    if not tokens or not all(t in permitted for t in tokens):
        return line
    raw_filter_info = ",".join(tokens)
    if not cols[7] or cols[7] == ".":
        cols[7] = f"RAW_FILTER={raw_filter_info}"
    elif "RAW_FILTER=" not in cols[7]:
        cols[7] = f"{cols[7]};RAW_FILTER={raw_filter_info}"
    cols[6] = "PASS"
    logging.getLogger("variantcentrifuge").info(
        "Rescued hotspot variant %s:%s %s>%s (RAW_FILTER=%s -> PASS)",
        cols[0],
        cols[1],
        cols[3],
        cols[4],
        raw_filter_info,
    )
    return "\t".join(cols) + ("\n" if line.endswith("\n") else "")


def apply_filter_exceptions_to_vcf(
    input_vcf: str,
    output_vcf: str,
    exceptions: dict[tuple[str, str, str, str], set[str]],
) -> int:
    """Rewrite VCF applying hotspot FILTER exceptions and return number of rescued records."""
    from .utils import smart_open

    rescued = 0
    has_raw_filter_hdr = False
    with smart_open(input_vcf, "r") as fin, smart_open(output_vcf, "w") as fout:
        for line in fin:
            if line.startswith("##INFO=<ID=RAW_FILTER,"):
                has_raw_filter_hdr = True
            elif line.startswith("#CHROM") and not has_raw_filter_hdr:
                fout.write(
                    "##INFO=<ID=RAW_FILTER,Number=.,Type=String,"
                    'Description="Raw FILTER before hotspot rescue">\n'
                )
            new_line = apply_filter_exceptions_to_line(line, exceptions)
            if new_line != line:
                rescued += 1
            fout.write(new_line)
    return rescued


def filter_vcf_to_transcripts(
    input_vcf: str,
    output_vcf: str,
    transcript_ids: set[str],
    *,
    index_output: bool = False,
    filter_exception_file: str | None = None,
    filter_exceptions: dict[tuple[str, str, str, str], set[str]] | None = None,
    fallback: str = "none",
    stats: dict | None = None,
) -> int:
    """Write a VCF containing one record per matching ANN transcript ID.

    With ``stats`` the run's counts are merged into it and the caller logs the summary
    (``log_transcript_filter_summary``); otherwise it is logged here.
    ``fallback="keep-gene"`` keeps every ANN entry of a gene with no matched transcript.
    """
    exceptions = filter_exceptions or load_filter_exceptions(filter_exception_file)
    input_records = 0
    kept_records = 0
    base_ids = {_base_transcript_id(tid) for tid in transcript_ids if tid.strip()}
    matched_features: set[str] = set()
    genes_before: set[str] = set()
    genes_after: set[str] = set()
    # ponytail: decided per input file, so per chunk in the parallel path; a gene whose
    # variants span two chunks may be kept in one and dropped in the other.
    keep_genes = (
        _genes_without_matched_transcript(input_vcf, transcript_ids, base_ids)
        if fallback == "keep-gene"
        else None
    )

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

        has_raw_filter_hdr = False
        for line in input_handle:
            if line.startswith("##INFO=<ID=RAW_FILTER,"):
                has_raw_filter_hdr = True
            elif line.startswith("#CHROM") and exceptions and not has_raw_filter_hdr:
                output_handle.write(
                    "##INFO=<ID=RAW_FILTER,Number=.,Type=String,"
                    'Description="Raw FILTER before hotspot rescue">\n'
                )
            if not line.startswith("#") and line.strip():
                input_records += 1
                if exceptions:
                    line = apply_filter_exceptions_to_line(line, exceptions)
            filtered = filter_vcf_line_by_transcripts(
                line,
                transcript_ids,
                base_ids=base_ids,
                matched_features=matched_features,
                genes_before=genes_before,
                genes_after=genes_after,
                keep_genes=keep_genes,
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

    run_stats = {
        "input_records": input_records,
        "kept_records": kept_records,
        "matched_features": matched_features,
        "genes_before": genes_before,
        "genes_after": genes_after,
    }
    if stats is None:
        log_transcript_filter_summary(run_stats, transcript_ids)
    else:
        merge_filter_stats(stats, run_stats)
    return kept_records
