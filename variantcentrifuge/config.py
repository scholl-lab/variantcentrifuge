# File: variantcentrifuge/config.py
# Location: variantcentrifuge/variantcentrifuge/config.py

"""
Configuration management module.

This module handles loading configuration from a JSON file.
All default values should now reside in config.json.
"""

import os
import json


def load_config(config_file=None):
    """
    Load configuration from a JSON file.

    If no config_file is provided, it attempts to load 'config.json'
    from the current directory. If it fails to find or parse the file,
    it raises an error.

    Parameters
    ----------
    config_file : str, optional
        Path to a configuration file in JSON format.

    Returns
    -------
    dict
        Configuration dictionary loaded from JSON.
    """
    if not config_file:
        config_file = "config.json"

    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file '{config_file}' not found.")

    with open(config_file, "r", encoding="utf-8") as f:
        try:
            config = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Error parsing JSON configuration: {e}")

    # config should now contain all keys loaded from the JSON file
    return config
