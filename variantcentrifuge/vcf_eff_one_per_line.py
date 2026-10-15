#!/usr/bin/env python3

# File: variantcentrifuge/vcf_eff_one_per_line.py
# Location: variantcentrifuge/variantcentrifuge/vcf_eff_one_per_line.py

"""
Module: vcf_eff_one_per_line.

Replicates the functionality of the original vcfEffOnePerLine.pl script in Python.
It reads VCF lines (uncompressed or bgzipped) that contain SNPeff or SnpSift
annotations (EFF or ANN in the INFO field), and splits multiple annotations
into separate lines. All other information is repeated to preserve the original
records.

This can be used either as a standalone CLI or imported as a module in a Python
pipeline.

Example usage (standalone):
---------------------------
    # From stdin to stdout
    cat input.vcf | python vcf_eff_one_per_line.py > output.vcf

    # Specifying input and output files
    python vcf_eff_one_per_line.py -i input.vcf -o output.vcf

    # Using bgzipped files
    python vcf_eff_one_per_line.py -i input.vcf.gz -o output.vcf

Functions:
----------
    split_vcf_effects(line: str) -> list[str]:
        Splits a single VCF line if multiple EFF or ANN annotations are present.

    process_vcf_file(input_file: str, output_file: str = None):
        Reads an uncompressed or bgzipped VCF file line by line and writes the
        split lines to stdout or a specified output file.

Command-line arguments:
-----------------------
    -i, --input    Path to input VCF file (or "-" for stdin).
    -o, --output   Path to output file (or "-" for stdout).
"""

import argparse
import contextlib
import gzip
import io
import logging
import os
import re
import shutil
import subprocess
import sys

logger = logging.getLogger(__name__)

ANN_GENE_NAME_INDEX = 3
ANN_GENE_ID_INDEX = 4
_NON_ALLELE_INFO_KEYS = {"ANN", "EFF", "LOF", "NMD", "CSQ"}


def _ann_gene_keys(ann_entry: str) -> set[str]:
    """Return non-empty gene name and gene ID keys from a SnpEff ANN entry."""
    parts = ann_entry.split("|")
    gene_keys = set()
    if len(parts) > ANN_GENE_NAME_INDEX and parts[ANN_GENE_NAME_INDEX]:
        gene_keys.add(parts[ANN_GENE_NAME_INDEX])
    if len(parts) > ANN_GENE_ID_INDEX and parts[ANN_GENE_ID_INDEX]:
        gene_keys.add(parts[ANN_GENE_ID_INDEX])
    return gene_keys


def _filter_gene_scoped_info_value(value: str, gene_keys: set[str], info_key: str) -> str | None:
    """Keep LOF/NMD entries whose gene symbol or ID matches the retained ANN gene."""
    kept_entries = []
    for entry in value.split(","):
        stripped = entry.strip()
        parts = stripped.strip("()").split("|")
        if len(parts) < 2:
            logger.warning("Dropping malformed %s entry from ANN split row: %s", info_key, entry)
            continue
        if parts[0] in gene_keys or parts[1] in gene_keys:
            kept_entries.append(entry)

    if not kept_entries:
        return None
    return ",".join(kept_entries)


def _filter_ann_for_alt(ann_value: str, alt: str, all_alts: list[str]) -> str | None:
    """Filter comma-separated ANN entries to those matching the given ALT allele."""
    entries = [e for e in ann_value.split(",") if e]
    if len(entries) <= 1:
        return ann_value
    matching = [e for e in entries if e.split("|", 1)[0] == alt]
    if matching:
        return ",".join(matching)
    if len(all_alts) > 1 and len(entries) == len(all_alts) and alt in all_alts:
        return entries[all_alts.index(alt)]
    if any(e.split("|", 1)[0] in set(all_alts) for e in entries):
        return None
    return ann_value


def _decompose_sample_column(
    fmt_keys: list[str], sample_str: str, alt_idx: int, num_alts: int
) -> str:
    """Decompose Number=A/R FORMAT fields and normalize GT for one ALT allele."""
    vals = sample_str.split(":")
    if len(vals) == 1 and vals[0] in {".", "./.", ".|."}:
        return sample_str

    has_alt_support = False
    for k, v in zip(fmt_keys, vals, strict=False):
        if k == "AF":
            parts = v.split(",")
            if len(parts) == num_alts and parts[alt_idx] not in {".", ""}:
                with contextlib.suppress(ValueError):
                    has_alt_support = has_alt_support or float(parts[alt_idx]) > 0
        elif k == "AD":
            parts = v.split(",")
            if len(parts) == num_alts + 1 and parts[alt_idx + 1] not in {".", ""}:
                with contextlib.suppress(ValueError):
                    has_alt_support = has_alt_support or int(parts[alt_idx + 1]) > 0

    out_vals: list[str] = []
    target = str(alt_idx + 1)
    for k, v in zip(fmt_keys, vals, strict=False):
        if k == "GT":
            sep = "|" if "|" in v else "/"
            alleles = re.split(r"[/|]", v)
            if all(a == "." for a in alleles):
                out_vals.append(v)
                continue
            mapped = ["." if a == "." else ("1" if a == target else "0") for a in alleles]
            if (
                num_alts > 1
                and alt_idx > 0
                and "1" not in mapped
                and "1" in alleles
                and all(a in {"0", "1", "."} for a in alleles)
                and has_alt_support
            ):
                mapped = ["0", "1"]
            elif len(mapped) > 2:
                if "1" in mapped and "0" in mapped:
                    mapped = ["0", "1"]
                elif all(x == "1" for x in mapped):
                    mapped = ["1", "1"]
                elif all(x == "0" for x in mapped):
                    mapped = ["0", "0"]
                else:
                    mapped = mapped[:2]
            out_vals.append(sep.join(mapped))
        elif "," in v and num_alts > 1:
            parts = v.split(",")
            if k == "AD" or (len(parts) == num_alts + 1 and k not in {"FT", "PGT", "PID"}):
                out_vals.append(f"{parts[0]},{parts[alt_idx + 1]}")
            elif k in {"AF", "VAF", "FA"} or (
                len(parts) == num_alts and k not in {"FT", "PGT", "PID", "AD"}
            ):
                out_vals.append(parts[alt_idx])
            else:
                out_vals.append(v)
        else:
            out_vals.append(v)
    if len(vals) > len(fmt_keys):
        out_vals.extend(vals[len(fmt_keys) :])
    return ":".join(out_vals)


def decompose_multiallelic_line(line: str) -> list[str]:
    """Decompose a multiallelic VCF line into biallelic lines with allele-matched ANN and AF/AD."""
    stripped = line.rstrip("\n")
    if stripped.startswith("#"):
        return [stripped]
    columns = stripped.split("\t")
    if len(columns) < 8:
        return [stripped]

    alts = columns[4].split(",")
    num_alts = len(alts)
    infos = columns[7].split(";")
    fmt_keys = columns[8].split(":") if len(columns) > 8 else []

    out_lines: list[str] = []
    for alt_idx, alt in enumerate(alts):
        row_cols = list(columns)
        row_cols[4] = alt
        new_infos: list[str] = []
        for inf in infos:
            if "=" not in inf:
                new_infos.append(inf)
                continue
            key, val = inf.split("=", 1)
            if key == "ANN":
                filtered_ann = _filter_ann_for_alt(val, alt, alts)
                if filtered_ann is not None:
                    new_infos.append(f"ANN={filtered_ann}")
            elif num_alts > 1 and key not in _NON_ALLELE_INFO_KEYS and "," in val:
                parts = val.split(",")
                if len(parts) == num_alts:
                    new_infos.append(f"{key}={parts[alt_idx]}")
                elif len(parts) == num_alts + 1:
                    new_infos.append(f"{key}={parts[0]},{parts[alt_idx + 1]}")
                else:
                    new_infos.append(inf)
            else:
                new_infos.append(inf)
        row_cols[7] = ";".join(new_infos) if new_infos else "."
        if fmt_keys and len(row_cols) > 9:
            for s_idx in range(9, len(row_cols)):
                row_cols[s_idx] = _decompose_sample_column(
                    fmt_keys, row_cols[s_idx], alt_idx, num_alts
                )
        out_lines.append("\t".join(row_cols))
    return out_lines


def _split_biallelic_line_effects(line: str) -> list[str]:
    """Split multiple EFF or ANN annotations on a single biallelic VCF line."""
    if line.startswith("#"):
        return [line]

    columns = line.strip().split("\t")
    info_field_num = 7
    if len(columns) <= info_field_num:
        return [line]

    info_str = columns[info_field_num]
    infos = info_str.split(";")

    effs: list[str] = []
    field_name = None
    other_info_parts = []
    gene_scoped_info_parts: list[tuple[str, str]] = []

    for inf in infos:
        match_eff = re.match(r"^(EFF)=(.*)", inf)
        match_ann = re.match(r"^(ANN)=(.*)", inf)
        match_gene_scoped = re.match(r"^(LOF|NMD)=(.*)", inf)
        if match_eff:
            field_name = match_eff.group(1)
            effs = match_eff.group(2).split(",")
        elif match_ann:
            field_name = match_ann.group(1)
            effs = match_ann.group(2).split(",")
        elif match_gene_scoped:
            gene_scoped_info_parts.append((match_gene_scoped.group(1), match_gene_scoped.group(2)))
        else:
            other_info_parts.append(inf)

    if not effs or len(effs) <= 1:
        return [line]

    pre_cols = columns[:info_field_num]
    post_cols = columns[info_field_num + 1 :]

    split_lines = []
    pre_string = "\t".join(pre_cols)
    post_string = "\t" + "\t".join(post_cols) if post_cols else ""

    for eff in effs:
        new_info_parts = list(filter(None, other_info_parts))
        if field_name == "ANN":
            new_info_parts.append(f"{field_name}={eff}")
            gene_keys = _ann_gene_keys(eff)
            for info_key, value in gene_scoped_info_parts:
                filtered_value = _filter_gene_scoped_info_value(value, gene_keys, info_key)
                if filtered_value is not None:
                    new_info_parts.append(f"{info_key}={filtered_value}")
        else:
            new_info_parts.extend(
                f"{info_key}={value}" for info_key, value in gene_scoped_info_parts
            )
            new_info_parts.append(f"{field_name}={eff}")

        new_info = ";".join(new_info_parts)
        new_line = f"{pre_string}\t{new_info}{post_string}"
        split_lines.append(new_line)

    return split_lines


def split_vcf_effects(line: str) -> list[str]:
    """Decompose multiallelic ALT records and split multiple EFF/ANN annotations one per line."""
    if line.startswith("#"):
        return [line]
    split_lines: list[str] = []
    for biallelic_line in decompose_multiallelic_line(line):
        split_lines.extend(_split_biallelic_line_effects(biallelic_line))
    return split_lines


def vcf_has_multiallelic(path: str) -> bool:
    """Return True if the VCF at `path` contains any multiallelic record (`ALT` with `,`)."""
    if not path or not os.path.exists(path):
        return False
    try:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[call-overload]
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t", 5)
                if len(cols) >= 5 and "," in cols[4]:
                    return True
    except Exception:
        return False
    return False


def process_vcf_file(
    input_file: str,
    output_file: str | None = None,
    *,
    split_transcripts: bool = True,
    index_output: bool = True,
) -> None:
    """Process a VCF file, decomposing multiallelic records and optionally splitting ANN/EFF."""
    use_bgzip = bool(
        output_file is not None
        and output_file != "-"
        and output_file.endswith(".gz")
        and shutil.which("bgzip")
    )
    use_gzip_fallback = bool(
        output_file is not None
        and output_file != "-"
        and output_file.endswith(".gz")
        and not use_bgzip
    )

    bgzip_proc = None
    out_fh = None
    if output_file is None or output_file == "-":
        out_handle = sys.stdout
    elif use_bgzip:
        out_fh = open(output_file, "wb")  # noqa: SIM115
        bgzip_proc = subprocess.Popen(
            ["bgzip", "-c"],
            stdin=subprocess.PIPE,
            stdout=out_fh,
            stderr=subprocess.PIPE,
        )
        assert bgzip_proc.stdin is not None
        out_handle = io.TextIOWrapper(bgzip_proc.stdin, encoding="utf-8")
    elif use_gzip_fallback:
        out_handle = gzip.open(output_file, mode="wt", encoding="utf-8")  # noqa: SIM115
    else:
        out_handle = open(output_file, "w", encoding="utf-8")  # noqa: SIM115

    if input_file == "-":
        in_handle = sys.stdin
    elif input_file.endswith(".gz"):
        in_handle = gzip.open(input_file, mode="rt", encoding="utf-8")  # noqa: SIM115
    else:
        in_handle = open(input_file, encoding="utf-8")  # noqa: SIM115

    line_fn = split_vcf_effects if split_transcripts else decompose_multiallelic_line
    try:
        with in_handle:
            for line in in_handle:
                for out_line in line_fn(line.rstrip("\n")):
                    out_handle.write(out_line + "\n")
    finally:
        if out_handle not in (sys.stdout, sys.stdin):
            out_handle.close()
        if use_bgzip and bgzip_proc is not None:
            bgzip_proc.wait()
            if out_fh is not None:
                out_fh.close()
            if bgzip_proc.returncode == 0 and index_output and shutil.which("tabix"):
                assert output_file is not None
                subprocess.run(["tabix", "-f", "-p", "vcf", output_file], check=True)


def normalize_vcf_multiallelics(
    input_file: str, output_file: str, *, index_output: bool = True
) -> None:
    """Decompose multiallelic VCF records into biallelic rows without splitting transcripts."""
    process_vcf_file(input_file, output_file, split_transcripts=False, index_output=index_output)


def main():
    """Serve as CLI entry point for splitting multiple EFF/ANN annotations in a VCF file."""
    parser = argparse.ArgumentParser(
        description=(
            "Split multiple EFF or ANN annotations in a VCF file into separate lines. "
            "Replicates the functionality of the original vcfEffOnePerLine.pl script."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        default="-",
        help="Path to input VCF file (uncompressed or bgzipped). Use '-' for stdin.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Path to output file. Use '-' for stdout. If '.gz' extension is used, output is bgzipped + tabixed.",
    )

    args = parser.parse_args()
    process_vcf_file(args.input, args.output)


if __name__ == "__main__":
    main()
