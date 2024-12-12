# File: variantcentrifuge/phenotype_filter.py
# Location: variantcentrifuge/variantcentrifuge/phenotype_filter.py

"""
Phenotype filtering module.

This module applies filters based on a list of sample IDs and a specified
column name. It replicates the logic from the filter_phenotypes.sh script,
filtering CSV/TSV lines so that only lines with sample IDs in the chosen column
are retained.
"""

def filter_phenotypes(lines, cfg):
    """
    Filter variants based on phenotype information (sample IDs).

    Parameters
    ----------
    lines : iterator
        Iterator of variant lines (strings).
    cfg : dict
        Configuration dictionary containing:
        - column_name: str, required
        - sample_file: str, optional
        - sample_list: str, optional (comma-separated)
        - output_file: str, optional (not directly needed here)
        - Delimiter info (optional):
          - input_delimiter: str (e.g., "," or "\t")
          - output_delimiter: str (e.g., "," or "\t")

    Returns
    -------
    iterator
        Iterator of filtered lines.
    """
    column_name = cfg.get("column_name")
    if not column_name:
        raise ValueError("Phenotype filtering requires 'column_name' in cfg.")

    # Determine input and output delimiters
    input_delim = cfg.get("input_delimiter", "\t")
    output_delim = cfg.get("output_delimiter", input_delim)

    # Load samples
    sample_list_str = cfg.get("sample_list", "")
    sample_file = cfg.get("sample_file", "")
    samples = []

    if sample_list_str and sample_file:
        # Warning: both sample_list and sample_file given, use sample_list_str
        pass

    if sample_list_str:
        samples = [x.strip() for x in sample_list_str.split(",") if x.strip()]
    elif sample_file:
        # Read samples from file
        with open(sample_file, "r", encoding="utf-8") as sf:
            for line in sf:
                s = line.strip()
                if s:
                    samples.append(s)
    else:
        raise ValueError("No sample_list or sample_file provided for phenotype filtering.")

    samples_set = set(samples)

    # Convert lines to a list for column detection (or we can peek first line)
    lines_iter = iter(lines)
    header = next(lines_iter, None)
    if header is None:
        # No data
        return

    # Parse header and find column_number
    header_fields = header.rstrip("\n").split(input_delim)
    if column_name not in header_fields:
        raise ValueError(f"Column '{column_name}' not found in header.")

    column_number = header_fields.index(column_name)

    # Yield the header line, converting to output_delim
    yield output_delim.join(header_fields)

    # For each subsequent line, check if sample ID matches
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
