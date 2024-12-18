# File: variantcentrifuge/links.py
# Location: variantcentrifuge/variantcentrifuge/links.py

"""
Module for adding configurable links to output tables.

This module provides functionality to add new columns containing URL links
to each variant line in the final output. The URLs and their substitution
logic are defined in the configuration file (config.json) under "links".

For example, if you have a link defined as:
"links": {
    "SpliceAI": "https://spliceailookup.broadinstitute.org/#variant={CHROM}-{POS}-{REF}-{ALT}&hg=38&bc=basic&distance=500&mask=0&ra=0"
}

This will create a new column named "SpliceAI" with a URL where {CHROM}, {POS}, {REF}, {ALT}
are replaced by the corresponding values from each row.
"""

from typing import List, Dict


def add_links_to_table(lines: List[str], link_configs: Dict[str, str]) -> List[str]:
    """
    Add link columns to the table lines based on provided link configurations.

    Parameters
    ----------
    lines : list of str
        The lines of the final table (including header).
    link_configs : dict of str to str
        A dictionary where keys are link column names and values are URL templates.
        URL templates contain placeholders like {CHROM}, {POS}, {REF}, {ALT} that
        are replaced with values from the corresponding columns.

    Returns
    -------
    list of str
        The modified lines with new link columns appended.

    Notes
    -----
    - Assumes the first line is a header line.
    - Each link column is appended to the end of the table.
    - Placeholders in the link template must match header column names exactly (case-sensitive).
    - If a placeholder column doesn't exist in the header, the link will be empty.
    """
    if not link_configs:
        return lines

    if not lines:
        return lines

    header = lines[0].rstrip("\n").split("\t")
    # Map column names to their indices for substitution
    col_index = {name: i for i, name in enumerate(header)}

    # Extend the header with new link columns
    new_header = header[:]
    link_names = sorted(link_configs.keys())
    for link_name in link_names:
        new_header.append(link_name)

    modified_lines = ["\t".join(new_header)]

    for line in lines[1:]:
        line = line.rstrip("\n")
        if not line.strip():
            # Empty or blank line
            modified_lines.append(line)
            continue
        fields = line.split("\t")

        # For each link, perform substitution
        for link_name in link_names:
            template = link_configs[link_name]
            # Replace placeholders with corresponding column values
            # If a column is missing, produce empty string for that link
            link_url = template
            for placeholder in col_index:
                placeholder_token = "{" + placeholder + "}"
                if placeholder_token in template:
                    if col_index[placeholder] < len(fields):
                        link_url = link_url.replace(placeholder_token, fields[col_index[placeholder]])
                    else:
                        # Column index out of range, set link_url empty
                        link_url = ""
                        break
            fields.append(link_url)

        modified_lines.append("\t".join(fields))

    return modified_lines
