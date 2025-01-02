#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File: variantcentrifuge/vcf_eff_one_per_line.py
# Location: variantcentrifuge/variantcentrifuge/vcf_eff_one_per_line.py


"""
Module: vcf_eff_one_per_line
============================
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

import sys
import re
import gzip
import argparse


def split_vcf_effects(line: str) -> list[str]:
    """
    Splits multiple EFF or ANN annotations in the INFO field of a VCF line
    into separate lines. If no EFF/ANN is found or only one is found,
    returns a list containing the original line.

    :param line: A single VCF record (line).
    :return: A list of one or more lines, each containing one EFF/ANN annotation.
    """
    # If this is a header line, return as is
    if line.startswith('#'):
        return [line]

    columns = line.strip().split('\t')
    # VCF columns: CHROM (0), POS (1), ID (2), REF (3), ALT (4),
    #              QUAL (5), FILTER (6), INFO (7), FORMAT (8), ...
    INFO_FIELD_NUM = 7

    # Guard: If the line doesn't have enough columns, return it as is
    if len(columns) <= INFO_FIELD_NUM:
        return [line]

    info_str = columns[INFO_FIELD_NUM]
    infos = info_str.split(';')

    # We look for EFF= or ANN= in the INFO field
    effs = []
    field_name = None
    other_info_parts = []

    for inf in infos:
        match_eff = re.match(r'^(EFF)=(.*)', inf)
        match_ann = re.match(r'^(ANN)=(.*)', inf)
        if match_eff:
            field_name = match_eff.group(1)
            effs = match_eff.group(2).split(',')
        elif match_ann:
            field_name = match_ann.group(1)
            effs = match_ann.group(2).split(',')
        else:
            other_info_parts.append(inf)

    # If no EFF/ANN found or only one annotation, return original line
    if not effs or len(effs) <= 1:
        return [line]

    # Otherwise, replicate line for each EFF/ANN
    pre_cols = columns[:INFO_FIELD_NUM]  # columns before INFO
    post_cols = columns[INFO_FIELD_NUM + 1:]  # columns after INFO

    split_lines = []
    pre_string = '\t'.join(pre_cols)
    post_string = '\t' + '\t'.join(post_cols) if post_cols else ''

    # Reconstruct lines with each effect/annotation on separate lines
    for eff in effs:
        # Build the final INFO field: combine leftover info with the single EFF/ANN
        new_info = ';'.join(filter(None, other_info_parts))
        if new_info:
            new_info += f';{field_name}={eff}'
        else:
            new_info = f'{field_name}={eff}'

        new_line = f'{pre_string}\t{new_info}{post_string}'
        split_lines.append(new_line)

    return split_lines


def process_vcf_file(input_file: str, output_file: str = None) -> None:
    """
    Processes a VCF file (bgzipped or uncompressed), splits multiple EFF/ANN
    annotations into separate lines, and writes the output to a file or stdout.

    :param input_file: Path to an uncompressed or bgzipped VCF file. Use '-'
                       to read from stdin.
    :param output_file: Path to output file. Use '-' for stdout. If None, stdout
                        is used by default.
    """
    if output_file is None or output_file == '-':
        out_handle = sys.stdout
    else:
        out_handle = open(output_file, 'w', encoding='utf-8')

    if input_file == '-':
        in_handle = sys.stdin
    elif input_file.endswith('.gz'):
        in_handle = gzip.open(input_file, mode='rt', encoding='utf-8')
    else:
        in_handle = open(input_file, 'r', encoding='utf-8')

    with in_handle, out_handle:
        for line in in_handle:
            # Strip only trailing newlines (avoid removing other whitespace)
            line = line.rstrip('\n')
            split_lines = split_vcf_effects(line)
            for out_line in split_lines:
                out_handle.write(out_line + '\n')


def main():
    """
    CLI entry point for splitting multiple EFF/ANN annotations in a VCF file.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Split multiple EFF or ANN annotations in a VCF file into separate lines. "
            "Replicates the functionality of the original vcfEffOnePerLine.pl script."
        )
    )
    parser.add_argument(
        '-i',
        '--input',
        default='-',
        help="Path to input VCF file (uncompressed or bgzipped). Use '-' for stdin."
    )
    parser.add_argument(
        '-o',
        '--output',
        default='-',
        help="Path to output file. Use '-' for stdout."
    )

    args = parser.parse_args()
    process_vcf_file(args.input, args.output)


if __name__ == '__main__':
    main()
