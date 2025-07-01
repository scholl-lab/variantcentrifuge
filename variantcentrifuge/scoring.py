"""
Scoring module for applying custom formulas to variant data.

This module reads scoring configurations and applies them to a pandas DataFrame
of annotated variants to generate custom scores.
"""

import json
import logging
import os
from typing import Any, Dict

import pandas as pd

logger = logging.getLogger("variantcentrifuge")


def read_scoring_config(config_path: str) -> Dict[str, Any]:
    """
    Read and parse the scoring configuration files from a directory.

    Expects two files:
    - variable_assignment_config.json: Maps DataFrame columns to shorter variable names.
    - formula_config.json: Contains the scoring formulas.

    Parameters
    ----------
    config_path : str
        Path to the scoring configuration directory.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing the parsed variables and formulas.
    """
    try:
        var_assign_path = os.path.join(config_path, "variable_assignment_config.json")
        formula_path = os.path.join(config_path, "formula_config.json")

        logger.debug(f"Reading scoring config from: {config_path}")

        with open(var_assign_path, "r", encoding="utf-8") as f:
            variable_assignment = json.load(f)

        with open(formula_path, "r", encoding="utf-8") as f:
            formula_config = json.load(f)

        return {
            "variables": variable_assignment.get("variables", {}),
            "formulas": formula_config.get("formulas", []),
        }
    except FileNotFoundError as e:
        logger.error(f"Scoring configuration file not found: {e.filename}")
        raise
    except json.JSONDecodeError as e:
        logger.error(f"Error parsing scoring configuration JSON: {e}")
        raise


def convert_to_numeric(series: pd.Series, default: float = 0.0) -> pd.Series:
    """
    Convert a series to numeric, handling empty strings and other non-numeric values.

    Parameters
    ----------
    series : pd.Series
        The series to convert.
    default : float
        The default value to use for non-numeric entries.

    Returns
    -------
    pd.Series
        The numeric series.
    """
    # Replace empty strings with NaN, then convert to numeric
    series = series.replace("", default).replace(".", default)
    return pd.to_numeric(series, errors="coerce").fillna(default)


def apply_scoring(df: pd.DataFrame, scoring_config: Dict[str, Any]) -> pd.DataFrame:
    """
    Apply scoring formulas to the DataFrame of variants.

    This function iterates through the scoring formulas, renames columns to match
    the variables used in the formulas, evaluates the formulas using pandas.eval(),
    and appends the resulting scores as new columns to the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing annotated variant data.
    scoring_config : Dict[str, Any]
        The parsed scoring configuration.

    Returns
    -------
    pd.DataFrame
        The DataFrame with new score columns added.
    """
    variables = scoring_config.get("variables", {})
    formulas = scoring_config.get("formulas", [])

    if not formulas:
        logger.debug("No scoring formulas provided. Skipping scoring.")
        return df

    logger.info(f"Applying {len(formulas)} scoring formula(s)...")

    # Create a copy to avoid modifying the original DataFrame in place
    scored_df = df.copy()

    # Create a mapping from original column names to formula variable names
    # and handle default values for missing columns.
    rename_map = {}
    created_vars = []  # Track variables created for missing columns
    for original_col, var_config in variables.items():
        if isinstance(var_config, str):  # Legacy format "target|default:value"
            parts = var_config.split("|")
            target_name = parts[0]
            default_val_str = (
                parts[1].split(":")[1] if len(parts) > 1 and "default" in parts[1] else "0"
            )
        elif isinstance(var_config, dict):  # New object format
            target_name = var_config.get("target")
            default_val_str = str(var_config.get("default", "0"))
        else:
            continue

        if original_col in scored_df.columns:
            # Convert numeric columns to proper numeric type
            if "AF" in original_col or "CADD" in original_col:
                scored_df[original_col] = convert_to_numeric(
                    scored_df[original_col], float(default_val_str)
                )
            rename_map[original_col] = target_name
        else:
            # If the column is missing, create it with the default value
            logger.warning(
                f"Scoring variable column '{original_col}' not found. "
                f"Using default value '{default_val_str}'."
            )
            try:
                default_value = float(default_val_str)
            except ValueError:
                default_value = 0
            scored_df[target_name] = default_value
            created_vars.append(target_name)

    # Rename columns for formula evaluation
    if rename_map:
        scored_df.rename(columns=rename_map, inplace=True)

    # Apply each formula
    for formula_item in formulas:
        for score_name, formula_str in formula_item.items():
            logger.debug(f"Evaluating formula for '{score_name}': {formula_str}")
            try:
                # Use pandas.eval for safe and efficient evaluation
                scored_df[score_name] = scored_df.eval(formula_str, engine="python")
                logger.info(f"Successfully calculated score column: '{score_name}'")
            except Exception as e:
                logger.error(f"Error evaluating formula for '{score_name}': {e}")
                # Add a column with NaN to indicate failure without crashing
                scored_df[score_name] = float("nan")

    # After scoring, we need to clean up:
    # 1. Rename columns back to their original names
    # 2. Remove any variables that were created for missing columns
    # 3. Keep the new score columns

    # First, rename columns back to original names
    inverse_rename_map = {v: k for k, v in rename_map.items()}
    scored_df.rename(columns=inverse_rename_map, inplace=True)

    # Now remove only the created variables (for missing columns)
    if created_vars:
        cols_to_keep = [col for col in scored_df.columns if col not in created_vars]
        scored_df = scored_df[cols_to_keep]

    return scored_df
