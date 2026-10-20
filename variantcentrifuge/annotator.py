"""
Unified Annotation Module for VariantCentrifuge.

This module provides comprehensive custom annotation capabilities that unify
BED file annotations, gene list annotations, and structured JSON gene data
into a single, powerful annotation system.

Key Features:
- Efficient interval-based annotation using IntervalTrees for BED files
- Gene list annotation with flexible naming
- Structured JSON gene data integration
- Unified output in a single 'Custom_Annotation' column
- Full integration with filtering and scoring systems

The module generates semicolon-separated key-value pairs in the format:
Region=promoter_X;InGeneList=cancer_panel;Panel=HereditaryCancer
"""

import json
import logging
import os
from typing import Any

import pandas as pd

# Handle optional dependency gracefully
try:
    from intervaltree import Interval, IntervalTree

    INTERVALTREE_AVAILABLE = True
except ImportError:
    INTERVALTREE_AVAILABLE = False
    IntervalTree = None
    Interval = None

logger = logging.getLogger("variantcentrifuge")


from .feature_loaders import (
    _load_bed_files,
    _load_gene_lists,
    _load_json_gene_data,
    _load_named_bed_files,
    _parse_bed_line,
    _sanitize_name,
    split_bed_spec,
)

__all__ = [
    "INTERVALTREE_AVAILABLE",
    "_add_json_annotations_as_columns",
    "_extract_genes_from_row",
    "_find_gene_list_matches",
    "_find_json_gene_matches",
    "_find_region_overlaps",
    "_load_bed_files",
    "_load_gene_lists",
    "_load_json_gene_data",
    "_parse_bed_line",
    "_sanitize_name",
    "annotate_dataframe_with_features",
    "get_annotation_summary",
    "load_custom_features",
    "validate_annotation_config",
]


def load_custom_features(cfg: dict[str, Any]) -> dict[str, Any]:
    """
    Load all custom annotation features from configuration.

    Parameters
    ----------
    cfg : Dict[str, Any]
        Configuration dictionary containing annotation file paths

    Returns
    -------
    Dict[str, Any]
        Dictionary containing all loaded annotation features
    """
    logger.info("Loading custom annotation features...")

    features: dict[str, Any] = {"regions_by_chrom": {}, "gene_lists": {}, "json_gene_data": {}}
    features["json_genes_as_columns"] = cfg.get("json_genes_as_columns", False)

    # Load BED files
    bed_files = cfg.get("annotate_bed_files", [])
    if bed_files:
        logger.info(f"Loading {len(bed_files)} BED file(s) for region annotation...")
        features["regions_by_chrom"] = _load_bed_files(bed_files)
        features["named_regions"] = _load_named_bed_files(bed_files)

    # Load gene lists
    gene_list_files = cfg.get("annotate_gene_lists", [])
    if gene_list_files:
        logger.info(f"Loading {len(gene_list_files)} gene list file(s)...")
        features["gene_lists"] = _load_gene_lists(gene_list_files)

    # Load JSON gene data
    json_files = cfg.get("annotate_json_genes", [])
    if json_files:
        logger.info(f"Loading {len(json_files)} JSON gene file(s)...")
        features["json_gene_data"] = _load_json_gene_data(
            json_files, cfg.get("json_gene_mapping", "")
        )

    # Log summary
    total_regions = sum(len(tree) for tree in features["regions_by_chrom"].values())
    total_gene_lists = len(features["gene_lists"])
    total_json_genes = len(features["json_gene_data"])

    logger.info(
        f"Custom annotation features loaded: {total_regions} regions, "
        f"{total_gene_lists} gene lists, {total_json_genes} JSON gene entries"
    )

    return features


def _extract_genes_from_row(row: pd.Series) -> set[str]:
    """
    Extract gene symbols from a variant row.

    Parameters
    ----------
    row : pd.Series
        Pandas Series representing a variant row

    Returns
    -------
    Set[str]
        Set of uppercase gene symbols
    """
    genes: set[str] = set()

    # Check common gene annotation columns
    gene_columns = ["GENE", "Gene", "gene_symbol", "symbol"]

    for col in gene_columns:
        if col in row:
            gene_value = str(row[col])
            if gene_value and gene_value.lower() not in ("nan", "none", ""):
                # Handle multiple genes separated by common delimiters
                for separator in [",", ";", "|", "&"]:
                    if separator in gene_value:
                        gene_parts = gene_value.split(separator)
                        genes.update(g.strip().upper() for g in gene_parts if g.strip())
                        break
                else:
                    # Single gene
                    genes.add(gene_value.strip().upper())

    return genes


def _overlapping(tree: Any, pos: int, ref: Any) -> Any:
    """Intervals hit by the REF span; BED is 0-based half-open, VCF POS is 1-based."""
    return tree.overlap(pos - 1, pos - 1 + max(len(str(ref)) if ref is not None else 1, 1))


def _find_region_overlaps(row: pd.Series, regions_by_chrom: dict[str, Any]) -> list[str]:
    """
    Find overlapping regions for a variant.

    Parameters
    ----------
    row : pd.Series
        Variant row
    regions_by_chrom : Dict[str, Any]
        Dictionary of interval trees by chromosome

    Returns
    -------
    List[str]
        List of region annotation strings
    """
    annotations: list[str] = []

    if not regions_by_chrom or not INTERVALTREE_AVAILABLE:
        return annotations

    # Extract position information
    chrom = str(row.get("CHROM", "")).removeprefix("chr")

    try:
        pos = int(row.get("POS", 0))
        if pos <= 0:
            return annotations
    except (ValueError, TypeError):
        return annotations

    if chrom in regions_by_chrom:
        for interval in _overlapping(regions_by_chrom[chrom], pos, row.get("REF")):
            annotations.append(f"Region={interval.data}")

    return annotations


def _find_gene_list_matches(genes: set[str], gene_lists: dict[str, set[str]]) -> list[str]:
    """
    Find gene list matches for a set of genes.

    Parameters
    ----------
    genes : Set[str]
        Set of gene symbols
    gene_lists : Dict[str, Set[str]]
        Dictionary of gene lists

    Returns
    -------
    List[str]
        List of gene list annotation strings
    """
    annotations = []

    for list_name, gene_set in gene_lists.items():
        if genes & gene_set:  # Intersection check
            annotations.append(f"InGeneList={list_name}")

    return annotations


def _find_json_gene_matches(
    genes: set[str], json_gene_data: dict[str, dict[str, Any]]
) -> list[str]:
    """
    Find JSON gene data matches for a set of genes.

    Parameters
    ----------
    genes : Set[str]
        Set of gene symbols
    json_gene_data : Dict[str, Dict[str, Any]]
        Dictionary of gene data

    Returns
    -------
    List[str]
        List of gene data annotation strings
    """
    annotations = []

    for gene in genes:
        if gene in json_gene_data:
            for key, value in json_gene_data[gene].items():
                annotations.append(f"{key}={value}")

    return annotations


def _add_json_annotations_as_columns(
    df: pd.DataFrame, json_gene_data: dict[str, dict[str, Any]]
) -> pd.DataFrame:
    """Add annotations from JSON data as separate columns to the DataFrame."""
    if not json_gene_data:
        return df

    # Identify all possible new column names from the JSON data
    all_new_columns: set[str] = set()
    for data in json_gene_data.values():
        all_new_columns.update(data.keys())

    if not all_new_columns:
        return df

    # Prepare a list of dictionaries, one for each row in the DataFrame
    annotations_for_df = []
    for row in df.itertuples(index=True):
        # Convert row to Series for _extract_genes_from_row
        row_series = df.loc[row.Index]
        assert isinstance(row_series, pd.Series)
        variant_genes = _extract_genes_from_row(row_series)
        row_annotations = {}
        # Find the first gene in the list that has an entry in our JSON data
        for gene in variant_genes:
            if gene in json_gene_data:
                row_annotations = json_gene_data[gene]
                break  # Use the first match
        annotations_for_df.append(row_annotations)

    # Create a temporary DataFrame from the annotations
    if annotations_for_df:
        # Ensure the index matches the original DataFrame for proper alignment
        annotations_df = pd.DataFrame(annotations_for_df, index=df.index)

        # Add any missing columns from all_new_columns to ensure all are present
        for col in all_new_columns:
            if col not in annotations_df.columns:
                annotations_df[col] = pd.NA

        # Concatenate the new columns to the original DataFrame
        # Sort column names for deterministic ordering
        df = pd.concat([df, annotations_df[sorted(all_new_columns)]], axis=1)

    return df


def annotate_dataframe_with_features(df: pd.DataFrame, features: dict[str, Any]) -> pd.DataFrame:
    """
    Add Custom_Annotation column to DataFrame with unified annotation data.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with variant data
    features : Dict[str, Any]
        Dictionary of loaded annotation features

    Returns
    -------
    pd.DataFrame
        DataFrame with added Custom_Annotation column
    """
    if df.empty:
        logger.warning("Empty DataFrame provided for annotation")
        df["Custom_Annotation"] = ""
        for name in features.get("named_regions", {}):
            df[name] = ""
        return df

    # Make a copy to work on
    df = df.copy()

    for name, trees in features.get("named_regions", {}).items():
        if name in df.columns:
            raise ValueError(f"--annotate-bed name '{name}' clashes with an existing column")
        refs = df["REF"] if "REF" in df.columns else ["N"] * len(df)
        chroms = df["CHROM"].astype(str).str.removeprefix("chr")
        df[name] = [
            "yes" if c in trees and _overlapping(trees[c], int(p), r) else "no"
            for c, p, r in zip(chroms, df["POS"], refs, strict=True)
        ]

    # Handle JSON as columns if flag is set
    if features.get("json_genes_as_columns") and features.get("json_gene_data"):
        logger.info("Annotating with JSON data as new columns...")
        df = _add_json_annotations_as_columns(df, features["json_gene_data"])

    # Handle all other annotations (and JSON if not done above)
    def annotate_variant(row: Any) -> str:
        annotations = []
        variant_genes = _extract_genes_from_row(row)

        annotations.extend(_find_region_overlaps(row, features["regions_by_chrom"]))
        annotations.extend(_find_gene_list_matches(variant_genes, features["gene_lists"]))

        # Conditionally add JSON data to the Custom_Annotation column
        if not features.get("json_genes_as_columns"):
            annotations.extend(_find_json_gene_matches(variant_genes, features["json_gene_data"]))

        return ";".join(sorted(set(annotations))) if annotations else ""

    df["Custom_Annotation"] = df.apply(annotate_variant, axis=1)

    annotated_variants = (df["Custom_Annotation"] != "").sum()
    logger.info(
        f"Custom annotation complete: {annotated_variants}/{len(df)} variants annotated "
        f"in 'Custom_Annotation' column."
    )

    return df


def get_annotation_summary(df: pd.DataFrame) -> dict[str, Any]:
    """
    Generate a summary of annotation results.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with Custom_Annotation column

    Returns
    -------
    Dict[str, Any]
        Dictionary with annotation statistics
    """
    if "Custom_Annotation" not in df.columns:
        return {"error": "Custom_Annotation column not found"}

    total_variants = len(df)
    annotated_variants = (df["Custom_Annotation"] != "").sum()

    # Count annotation types
    annotation_types: dict[str, int] = {}
    for annotations in df["Custom_Annotation"]:
        if annotations:
            for annotation in annotations.split(";"):
                if "=" in annotation:
                    key = annotation.split("=")[0]
                    annotation_types[key] = annotation_types.get(key, 0) + 1

    return {
        "total_variants": total_variants,
        "annotated_variants": annotated_variants,
        "annotation_rate": annotated_variants / total_variants if total_variants > 0 else 0,
        "annotation_types": annotation_types,
    }


def validate_annotation_config(cfg: dict[str, Any]) -> list[str]:
    """
    Validate annotation configuration and return any errors.

    Parameters
    ----------
    cfg : Dict[str, Any]
        Configuration dictionary

    Returns
    -------
    List[str]
        List of validation error messages
    """
    errors = []

    # Check BED files
    for bed_spec in cfg.get("annotate_bed_files", []):
        bed_file = split_bed_spec(bed_spec)[1]
        if not os.path.exists(bed_file):
            errors.append(f"BED file not found: {bed_file}")
        elif not os.path.isfile(bed_file):
            errors.append(f"BED path is not a file: {bed_file}")

    # Check gene list files
    for gene_file in cfg.get("annotate_gene_lists", []):
        if not os.path.exists(gene_file):
            errors.append(f"Gene list file not found: {gene_file}")
        elif not os.path.isfile(gene_file):
            errors.append(f"Gene list path is not a file: {gene_file}")

    # Check JSON files and mapping
    json_files = cfg.get("annotate_json_genes", [])
    json_mapping = cfg.get("json_gene_mapping", "")

    if json_files and not json_mapping:
        errors.append("JSON gene files specified but no json_gene_mapping provided")
    elif json_mapping and not json_files:
        errors.append("JSON gene mapping provided but no JSON files specified")

    if json_mapping:
        try:
            mapping = json.loads(json_mapping)
            if not isinstance(mapping, dict):
                errors.append("JSON gene mapping must be a JSON object")
            elif "identifier" not in mapping:
                errors.append("JSON gene mapping must include 'identifier' field")
        except json.JSONDecodeError:
            errors.append("Invalid JSON syntax in gene mapping")

    for json_file in json_files:
        if not os.path.exists(json_file):
            errors.append(f"JSON gene file not found: {json_file}")
        elif not os.path.isfile(json_file):
            errors.append(f"JSON gene path is not a file: {json_file}")

    # Check intervaltree availability if BED files are specified
    if cfg.get("annotate_bed_files") and not INTERVALTREE_AVAILABLE:
        errors.append("intervaltree package required for BED file annotation but not installed")

    return errors
