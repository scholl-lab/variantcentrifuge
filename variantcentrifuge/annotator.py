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
import re
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


def _sanitize_name(name: str) -> str:
    """
    Sanitize a name to be safe for use in annotations.

    Parameters
    ----------
    name : str
        Raw name string

    Returns
    -------
    str
        Sanitized name with only alphanumeric characters and underscores
    """
    # Remove extension and path
    base_name = os.path.splitext(os.path.basename(name))[0]
    # Replace non-alphanumeric characters with underscores
    sanitized = re.sub(r"[^a-zA-Z0-9_]", "_", base_name)
    # Remove leading/trailing underscores and collapse multiple underscores
    sanitized = re.sub(r"_+", "_", sanitized).strip("_")
    # Ensure not empty
    return sanitized if sanitized else "unnamed"


def _parse_bed_line(line: str) -> tuple | None:
    """
    Parse a single BED file line.

    Parameters
    ----------
    line : str
        BED format line

    Returns
    -------
    tuple or None
        Tuple of (chrom, start, end, name) or None if invalid
    """
    line = line.strip()
    if not line or line.startswith(("#", "track", "browser")):
        return None

    parts = line.split("\t")
    if len(parts) < 3:
        return None

    try:
        chrom = parts[0].replace("chr", "")
        start = int(parts[1])
        end = int(parts[2])
        name = parts[3] if len(parts) > 3 else f"region_{start}_{end}"
        return (chrom, start, end, name)
    except (ValueError, IndexError):
        return None


def _load_bed_files(bed_files: list[str]) -> dict[str, Any]:
    """
    Load BED files into interval trees for efficient overlap detection.

    Parameters
    ----------
    bed_files : List[str]
        List of BED file paths

    Returns
    -------
    Dict[str, Any]
        Dictionary with chromosome-keyed interval trees
    """
    if not INTERVALTREE_AVAILABLE:
        logger.warning("intervaltree not available. BED annotation disabled.")
        return {}

    regions_by_chrom = {}

    for file_path in bed_files:
        if not os.path.exists(file_path):
            logger.warning(f"BED file not found: {file_path}")
            continue

        file_name = _sanitize_name(file_path)
        regions_loaded = 0

        try:
            with open(file_path, encoding="utf-8") as f:
                for line_num, line in enumerate(f, 1):
                    parsed = _parse_bed_line(line)
                    if parsed is None:
                        continue

                    chrom, start, end, name = parsed

                    # Ensure interval is valid
                    if start >= end:
                        logger.debug(f"Invalid interval in {file_path}:{line_num}: {start}-{end}")
                        continue

                    if chrom not in regions_by_chrom:
                        regions_by_chrom[chrom] = IntervalTree()

                    # Store both region name and file source
                    annotation_name = f"{name}_{file_name}"
                    regions_by_chrom[chrom][start:end] = annotation_name
                    regions_loaded += 1

            logger.info(f"Loaded {regions_loaded} regions from BED file: {file_path}")

        except Exception as e:
            logger.error(f"Failed to load BED file {file_path}: {e}")

    return regions_by_chrom


def _load_gene_lists(gene_list_files: list[str]) -> dict[str, set[str]]:
    """
    Load gene lists from text files.

    Parameters
    ----------
    gene_list_files : List[str]
        List of gene list file paths

    Returns
    -------
    Dict[str, Set[str]]
        Dictionary mapping list names to sets of gene symbols (uppercase)
    """
    gene_lists = {}

    for file_path in gene_list_files:
        if not os.path.exists(file_path):
            logger.warning(f"Gene list file not found: {file_path}")
            continue

        list_name = _sanitize_name(file_path)
        genes_loaded = 0

        try:
            with open(file_path, encoding="utf-8") as f:
                genes = set()
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        # Handle comma-separated genes on a single line
                        for gene in line.split(","):
                            gene = gene.strip()
                            if gene:
                                genes.add(gene.upper())
                                genes_loaded += 1

                gene_lists[list_name] = genes

            logger.info(f"Loaded {genes_loaded} genes from list: {file_path}")

        except Exception as e:
            logger.error(f"Failed to load gene list {file_path}: {e}")

    return gene_lists


def _load_json_gene_data(json_files: list[str], mapping_config: str) -> dict[str, dict[str, Any]]:
    """
    Load structured gene data from JSON files.

    Parameters
    ----------
    json_files : List[str]
        List of JSON file paths
    mapping_config : str
        JSON string specifying field mapping

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Dictionary mapping gene symbols to their data dictionaries
    """
    if not mapping_config:
        logger.warning("No JSON gene mapping provided. Skipping JSON gene data.")
        return {}

    try:
        mapping = json.loads(mapping_config)
        identifier_key = mapping.get("identifier", "gene_symbol")
        data_fields = mapping.get("dataFields", [])

        if not data_fields:
            logger.warning("No data fields specified in JSON mapping. Skipping JSON gene data.")
            return {}

    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON mapping configuration: {e}")
        return {}

    json_gene_data = {}

    for file_path in json_files:
        if not os.path.exists(file_path):
            logger.warning(f"JSON gene file not found: {file_path}")
            continue

        genes_loaded = 0

        try:
            with open(file_path, encoding="utf-8") as f:
                data = json.load(f)

            # Handle both list and dict formats
            items = (
                data if isinstance(data, list) else data.values() if isinstance(data, dict) else []
            )

            for item in items:
                if not isinstance(item, dict):
                    continue

                identifier = item.get(identifier_key)
                if not identifier:
                    continue

                gene_symbol = str(identifier).upper()
                gene_annotations = {}

                for field in data_fields:
                    value = item.get(field)
                    if value is not None:
                        # Convert value to string and sanitize
                        str_value = str(value).strip()
                        if str_value:
                            # Sanitize field name
                            clean_field = re.sub(r"[^a-zA-Z0-9_]", "_", field)
                            gene_annotations[clean_field] = str_value

                if gene_annotations:
                    json_gene_data[gene_symbol] = gene_annotations
                    genes_loaded += 1

            logger.info(f"Loaded data for {genes_loaded} genes from JSON: {file_path}")

        except Exception as e:
            logger.error(f"Failed to load JSON gene file {file_path}: {e}")

    return json_gene_data


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
    chrom = str(row.get("CHROM", "")).replace("chr", "")

    try:
        pos = int(row.get("POS", 0))
        if pos <= 0:
            return annotations
    except (ValueError, TypeError):
        return annotations

    # Find overlapping intervals
    if chrom in regions_by_chrom:
        overlaps = regions_by_chrom[chrom][pos]
        for interval in overlaps:
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
        return df

    # Make a copy to work on
    df = df.copy()

    # Handle JSON as columns if flag is set
    if features.get("json_genes_as_columns") and features.get("json_gene_data"):
        logger.info("Annotating with JSON data as new columns...")
        df = _add_json_annotations_as_columns(df, features["json_gene_data"])

    # Handle all other annotations (and JSON if not done above)
    def annotate_variant(row):
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
    for bed_file in cfg.get("annotate_bed_files", []):
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
