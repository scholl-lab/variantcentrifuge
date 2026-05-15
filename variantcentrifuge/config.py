# File: variantcentrifuge/config.py
# Location: variantcentrifuge/variantcentrifuge/config.py

"""
Configuration management module.

This module handles loading configuration from a JSON file.
All default values reside in config.json, which is included in
the installed package directory.

If no config_file is provided, this module attempts to load the default
config.json from the package installation directory.
"""

import json
import os
from collections.abc import Iterable
from typing import Any


def load_config(config_file: str | None = None) -> dict[str, Any]:
    """
    Load configuration from a JSON file.

    If no config_file is provided, the function attempts to load the
    'config.json' from the installed package directory. If it fails
    to find or parse the file, it raises an error.

    Parameters
    ----------
    config_file : str, optional
        Path to a configuration file in JSON format. If None, defaults to
        the package-installed 'config.json'.

    Returns
    -------
    dict
        Configuration dictionary loaded from the JSON file.

    Raises
    ------
    FileNotFoundError
        If the specified configuration file does not exist.
    ValueError
        If there is an error parsing the JSON configuration file.
    """
    if not config_file:
        # Use the package's installed config.json
        config_file = os.path.join(os.path.dirname(__file__), "config.json")

    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file '{config_file}' not found.")

    with open(config_file, encoding="utf-8") as f:
        try:
            config = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Error parsing JSON configuration: {e}") from e

    if not isinstance(config, dict):
        raise ValueError(
            f"Configuration file must contain a JSON object, got {type(config).__name__}"
        )
    return config


def _annotation_value_as_list(value: Any) -> list[str]:
    """Normalize annotation config values to a list of non-empty strings."""
    if value is None:
        return []
    if isinstance(value, str):
        stripped = value.strip()
        return [stripped] if stripped else []
    if isinstance(value, Iterable):
        normalized = []
        for item in value:
            if item is None:
                continue
            item_str = str(item).strip()
            if item_str:
                normalized.append(item_str)
        return normalized
    item_str = str(value).strip()
    return [item_str] if item_str else []


def _first_non_empty_annotation_value(config: dict[str, Any], keys: list[str]) -> list[str]:
    """Return the first non-empty normalized annotation value for key aliases."""
    for key in keys:
        value = _annotation_value_as_list(config.get(key))
        if value:
            return value
    return []


def normalize_annotation_config(config: dict[str, Any]) -> dict[str, Any]:
    """Normalize custom annotation config aliases in-place.

    Canonical internal keys:
    - annotate_bed_files
    - annotate_gene_lists

    Compatibility aliases:
    - annotate_bed
    - annotate_gene_list
    - annotate_gene_list_files
    """
    bed_files = _first_non_empty_annotation_value(
        config,
        ["annotate_bed_files", "annotate_bed"],
    )
    gene_lists = _first_non_empty_annotation_value(
        config,
        ["annotate_gene_lists", "annotate_gene_list_files", "annotate_gene_list"],
    )

    config["annotate_bed_files"] = bed_files
    config["annotate_bed"] = bed_files
    config["annotate_gene_lists"] = gene_lists
    config["annotate_gene_list_files"] = gene_lists
    config["annotate_gene_list"] = gene_lists
    return config
