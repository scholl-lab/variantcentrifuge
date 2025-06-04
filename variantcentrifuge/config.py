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
from typing import Any, Dict, Optional


def load_config(config_file: Optional[str] = None) -> Dict[str, Any]:
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

    with open(config_file, "r", encoding="utf-8") as f:
        try:
            config = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Error parsing JSON configuration: {e}")

    return config
