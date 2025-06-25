# File: variantcentrifuge/phenotype_filter.py
# Location: variantcentrifuge/variantcentrifuge/phenotype_filter.py

"""
Phenotype filtering module.

This module applies filters based on a list of sample IDs and a specified
column name. It replicates the logic from the filter_phenotypes.sh script,
filtering CSV/TSV lines so that only lines with sample IDs in the chosen column
are retained.
"""

import logging
from typing import Any, Dict, Iterator

logger = logging.getLogger("variantcentrifuge")


def filter_phenotypes(lines: Iterator[str], cfg: Dict[str, Any]) -> Iterator[str]:
    """
    Filter variants based on phenotype information (sample IDs).

    This function reads lines from an input iterator (representing a CSV/TSV file),
    identifies the column specified by `column_name` in the configuration, and then
    filters out all lines whose sample ID (in that column) is not in the provided
    sample list.

    Parameters
    ----------
    lines : Iterator[str]
        Iterator of lines (strings) representing CSV/TSV data.
    cfg : dict
        Configuration dictionary. Must contain:

        - "column_name": str
            The name of the column containing sample IDs.
        - "sample_list": str, optional
            Comma-separated list of sample IDs.
        - "sample_file": str, optional
            Path to a file containing sample IDs, one per line.
        - "input_delimiter": str, optional
            Delimiter used in the input file (default "\t").
        - "output_delimiter": str, optional
            Delimiter used in the output file (defaults to input_delimiter).

    Yields
    ------
    str
        Filtered lines matching the given sample criteria.

    Raises
    ------
    ValueError
        If "column_name" is not provided, or if neither "sample_list" nor "sample_file"
        is provided, or if "column_name" is not found in the header.
    """
    column_name = cfg.get("column_name")
    if not column_name:
        raise ValueError("Phenotype filtering requires 'column_name' in cfg.")

    input_delim = cfg.get("input_delimiter", "\t")
    output_delim = cfg.get("output_delimiter", input_delim)

    # Load samples
    sample_list_str = cfg.get("sample_list", "")
    sample_file = cfg.get("sample_file", "")
    samples = []

    if sample_list_str and sample_file:
        # Both provided; sample_list_str takes precedence as per instructions
        logger.debug("Both sample_list and sample_file provided. Using sample_list.")

    if sample_list_str:
        samples = [x.strip() for x in sample_list_str.split(",") if x.strip()]
    elif sample_file:
        with open(sample_file, "r", encoding="utf-8") as sf:
            for line in sf:
                s = line.strip()
                if s:
                    samples.append(s)
    else:
        raise ValueError("No sample_list or sample_file provided for phenotype filtering.")

    samples_set = set(samples)

    lines_iter = iter(lines)
    header = next(lines_iter, None)
    if header is None:
        # No data
        return

    header_fields = header.rstrip("\n").split(input_delim)
    if column_name not in header_fields:
        raise ValueError(f"Column '{column_name}' not found in header.")

    column_number = header_fields.index(column_name)

    yield output_delim.join(header_fields)

    for line in lines_iter:
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split(input_delim)
        if len(fields) <= column_number:
            continue
        sample_id = fields[column_number]
        if sample_id in samples_set:
            yield output_delim.join(fields)
