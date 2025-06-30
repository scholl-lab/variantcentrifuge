"""
PED file reader for pedigree information.

This module provides functionality to parse standard 6-column PED files
used in genetic analysis to define family relationships.
"""

import logging
from typing import Dict, Any
import pandas as pd

logger = logging.getLogger(__name__)


def read_pedigree(file_path: str) -> Dict[str, Dict[str, Any]]:
    """
    Parses a PED file into a dictionary keyed by sample ID.

    Args:
        file_path: Path to the PED file

    Returns:
        Dictionary mapping sample IDs to their pedigree information

    Raises:
        ValueError: If the PED file is invalid or cannot be parsed
    """
    try:
        # First, read without column names to check format
        try:
            raw_df = pd.read_csv(file_path, sep=r"\s+", header=None, dtype=str, comment="#")
        except pd.errors.EmptyDataError:
            raise ValueError("PED file is empty")

        if raw_df.empty:
            raise ValueError("PED file is empty")

        # Check if all rows have exactly 6 columns
        if raw_df.shape[1] != 6:
            raise ValueError(f"PED file must have exactly 6 columns, found {raw_df.shape[1]}")

        # Now read with proper column names
        ped_df = pd.read_csv(
            file_path,
            sep=r"\s+",
            header=None,
            names=["family_id", "sample_id", "father_id", "mother_id", "sex", "affected_status"],
            dtype=str,
            comment="#",
        )

        pedigree_data = {}
        for _, row in ped_df.iterrows():
            sample_id = row["sample_id"]
            if pd.isna(sample_id) or sample_id == "":
                logger.warning("Skipping row with empty sample ID")
                continue

            pedigree_data[sample_id] = {
                "family_id": (
                    row["family_id"]
                    if not pd.isna(row["family_id"]) and row["family_id"] != "."
                    else ""
                ),
                "sample_id": sample_id,
                "father_id": (
                    row["father_id"]
                    if not pd.isna(row["father_id"]) and row["father_id"] != "."
                    else "0"
                ),
                "mother_id": (
                    row["mother_id"]
                    if not pd.isna(row["mother_id"]) and row["mother_id"] != "."
                    else "0"
                ),
                "sex": row["sex"] if not pd.isna(row["sex"]) and row["sex"] != "." else "0",
                "affected_status": (
                    row["affected_status"]
                    if not pd.isna(row["affected_status"]) and row["affected_status"] != "."
                    else "0"
                ),
            }

        logger.info(f"Successfully parsed PED file with {len(pedigree_data)} individuals")
        return pedigree_data

    except Exception as e:
        raise ValueError(f"Failed to parse PED file: {e}")


def get_parents(sample_id: str, pedigree_data: Dict[str, Dict[str, Any]]) -> tuple:
    """
    Get the parent IDs for a given sample.

    Args:
        sample_id: The sample ID to get parents for
        pedigree_data: The pedigree data dictionary

    Returns:
        Tuple of (father_id, mother_id), where '0' indicates no parent
    """
    if sample_id not in pedigree_data:
        return (None, None)

    sample_info = pedigree_data[sample_id]
    father_id = sample_info.get("father_id", "0")
    mother_id = sample_info.get("mother_id", "0")

    return (father_id if father_id != "0" else None, mother_id if mother_id != "0" else None)


def is_affected(sample_id: str, pedigree_data: Dict[str, Dict[str, Any]]) -> bool:
    """
    Check if a sample is affected according to the pedigree.

    Args:
        sample_id: The sample ID to check
        pedigree_data: The pedigree data dictionary

    Returns:
        True if affected (status = 2), False otherwise
    """
    if sample_id not in pedigree_data:
        return False

    status = pedigree_data[sample_id].get("affected_status", "0")
    return status == "2"


def get_family_members(sample_id: str, pedigree_data: Dict[str, Dict[str, Any]]) -> list:
    """
    Get all family members for a given sample.

    Args:
        sample_id: The sample ID to get family for
        pedigree_data: The pedigree data dictionary

    Returns:
        List of sample IDs in the same family
    """
    if sample_id not in pedigree_data:
        return []

    family_id = pedigree_data[sample_id]["family_id"]
    return [sid for sid, info in pedigree_data.items() if info["family_id"] == family_id]
