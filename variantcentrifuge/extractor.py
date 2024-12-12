# File: variantcentrifuge/extractor.py
# Location: variantcentrifuge/variantcentrifuge/extractor.py

"""
Field extraction module.

This module uses SnpSift extractFields to extract specific fields
from VCF records.
"""

from .utils import run_command_stream


def extract_fields(variant_stream, fields, cfg):
    """
    Extract specified fields from variant records.

    Parameters
    ----------
    variant_stream : iterator
        Iterator over variant lines.
    fields : str
        Space-separated list of fields to extract.
    cfg : dict
        Configuration dictionary.

    Returns
    -------
    iterator
        Iterator over lines with extracted fields.
    """
    field_list = fields.strip().split()
    cmd = ["SnpSift", "extractFields", "-s", ",", "-e", "NA", "-"] + field_list
    return run_command_stream(cmd, input_stream=variant_stream)
