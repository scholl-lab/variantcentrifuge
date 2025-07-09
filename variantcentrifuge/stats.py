# File: variantcentrifuge/stats.py
# Location: variantcentrifuge/variantcentrifuge/stats.py

"""
Statistics module for variantcentrifuge.

Provides functions to compute:
- Basic variant-level statistics.
- Comprehensive gene-level statistics.
- Impact and variant type summaries.
- Merging and formatting of these stats.

All functions return DataFrames suitable for further processing.
"""

import logging
import os
from typing import Dict, Optional, Set

import pandas as pd

from .helpers import extract_sample_and_genotype
from .stats_engine import StatsEngine

logger = logging.getLogger("variantcentrifuge")


def compute_basic_stats(df: pd.DataFrame, all_samples: Set[str]) -> pd.DataFrame:
    """
    Compute basic statistics about the dataset.

    Including:
    - Number of variants
    - Number of samples
    - Number of genes
    - Het/Hom genotype counts
    - Variant type and impact counts (if columns are present)

    Parameters
    ----------
    df : pd.DataFrame
        Input variants DataFrame. Expected to have columns "GENE", "GT" and optionally "EFFECT", "IMPACT".
    all_samples : set of str
        Set of all sample names.

    Returns
    -------
    pd.DataFrame
        A DataFrame with 'metric' and 'value' columns listing basic stats.
    """
    logger.debug("Computing basic stats...")
    num_variants = len(df)
    num_samples = len(all_samples)
    num_genes = df["GENE"].nunique()

    het_counts = 0
    hom_counts = 0
    for val in df["GT"]:
        if isinstance(val, str) and val.strip():
            for g in val.split(";"):
                g = g.strip()
                _, genotype = extract_sample_and_genotype(g)
                if genotype in ["0/1", "1/0"]:
                    het_counts += 1
                elif genotype == "1/1":
                    hom_counts += 1

    metric_rows = [
        ["Number of variants", str(num_variants)],
        ["Number of samples", str(num_samples)],
        ["Number of genes", str(num_genes)],
        ["Het counts", str(het_counts)],
        ["Hom counts", str(hom_counts)],
    ]

    # Count variant types if "EFFECT" present
    if "EFFECT" in df.columns:
        variant_types = df["EFFECT"].value_counts()
        for vt, count in variant_types.items():
            metric_rows.append([f"Variant_type_{vt}", str(count)])

    # Count impact types if "IMPACT" present
    if "IMPACT" in df.columns:
        impact_types = df["IMPACT"].value_counts()
        for it, count in impact_types.items():
            metric_rows.append([f"Impact_type_{it}", str(count)])

    basic_stats_df = pd.DataFrame(metric_rows, columns=["metric", "value"])
    logger.debug("Basic stats computed.")
    return basic_stats_df


def compute_gene_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute gene-level aggregated stats (sum of proband/control counts and alleles).

    Parameters
    ----------
    df : pd.DataFrame
        Input variants DataFrame with assigned case/control counts.
        Expected columns: "GENE", "proband_count", "control_count", "proband_allele_count", "control_allele_count".

    Returns
    -------
    pd.DataFrame
        Gene-level summary DataFrame with columns: GENE, proband_count, control_count, proband_allele_count, control_allele_count.
    """
    logger.debug("Computing gene-level stats...")
    for col in [
        "proband_count",
        "control_count",
        "proband_allele_count",
        "control_allele_count",
    ]:
        if col not in df.columns:
            df[col] = 0
    grouped = (
        df.groupby("GENE")
        .agg(
            {
                "proband_count": "sum",
                "control_count": "sum",
                "proband_allele_count": "sum",
                "control_allele_count": "sum",
            }
        )
        .reset_index()
    )
    logger.debug("Gene-level stats computed.")
    return grouped


def compute_impact_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute a per-gene impact summary if the "IMPACT" column exists.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with "GENE" and "IMPACT" columns.

    Returns
    -------
    pd.DataFrame
        A pivoted table of gene vs. impact counts. Columns for each impact type.
        If columns are missing, returns an empty DataFrame.
    """
    logger.debug("Computing impact summary...")
    if "GENE" not in df.columns or "IMPACT" not in df.columns:
        logger.debug("No IMPACT or GENE column, returning empty impact summary.")
        return pd.DataFrame()
    impact_counts = df.groupby(["GENE", "IMPACT"]).size().reset_index(name="count")
    pivot_impact = (
        impact_counts.pivot(index="GENE", columns="IMPACT", values="count").fillna(0).reset_index()
    )
    logger.debug("Impact summary computed.")
    return pivot_impact


def compute_variant_type_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute a per-gene variant type summary if the "EFFECT" column exists.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with "GENE" and "EFFECT" columns.

    Returns
    -------
    pd.DataFrame
        A pivoted table of gene vs. variant types. Columns for each variant type.
        If columns are missing, returns an empty DataFrame.
    """
    logger.debug("Computing variant type summary...")
    if "GENE" not in df.columns or "EFFECT" not in df.columns:
        logger.debug("No EFFECT or GENE column, returning empty variant type summary.")
        return pd.DataFrame()
    type_counts = df.groupby(["GENE", "EFFECT"]).size().reset_index(name="count")
    pivot_types = (
        type_counts.pivot(index="GENE", columns="EFFECT", values="count").fillna(0).reset_index()
    )
    logger.debug("Variant type summary computed.")
    return pivot_types


def compute_custom_stats(
    df: pd.DataFrame, stats_config: Optional[str] = None
) -> Dict[str, pd.DataFrame]:
    """
    Compute custom statistics based on configuration.

    Parameters
    ----------
    df : pd.DataFrame
        Input variants DataFrame
    stats_config : Optional[str]
        Path to statistics configuration JSON file.
        If None, uses default configuration.

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary with statistics results by category (dataset, genes, groups)
    """
    if stats_config is None:
        # Use default configuration
        default_config_path = os.path.join(os.path.dirname(__file__), "default_stats_config.json")
        stats_config = default_config_path

    try:
        engine = StatsEngine(stats_config)
        results = engine.compute(df)
        logger.info(f"Computed custom statistics using config: {stats_config}")
        return results
    except Exception as e:
        logger.error(f"Failed to compute custom statistics: {e}")
        return {}


def merge_and_format_stats(
    gene_stats: pd.DataFrame,
    impact_summary: pd.DataFrame,
    variant_type_summary: pd.DataFrame,
    custom_stats: Optional[Dict[str, pd.DataFrame]] = None,
) -> pd.DataFrame:
    """
    Merge gene_stats with impact_summary and variant_type_summary into a single DataFrame.

    Parameters
    ----------
    gene_stats : pd.DataFrame
        DataFrame of gene-level aggregated stats.
    impact_summary : pd.DataFrame
        DataFrame of gene vs. impact counts.
    variant_type_summary : pd.DataFrame
        DataFrame of gene vs. variant type counts.
    custom_stats : Optional[Dict[str, pd.DataFrame]]
        Dictionary of custom statistics from StatsEngine

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with all gene-level stats, filling missing values with 0.
    """
    logger.debug("Merging gene stats with impact and variant type summaries...")
    merged = gene_stats.copy()
    if not impact_summary.empty:
        merged = pd.merge(merged, impact_summary, on="GENE", how="left")
    if not variant_type_summary.empty:
        merged = pd.merge(merged, variant_type_summary, on="GENE", how="left")

    # Merge custom gene-level stats if available
    if custom_stats and "genes" in custom_stats and not custom_stats["genes"].empty:
        custom_gene_stats = custom_stats["genes"]
        if "GENE" in custom_gene_stats.columns:
            merged = pd.merge(merged, custom_gene_stats, on="GENE", how="left")
            logger.debug("Merged custom gene-level statistics")

    merged = merged.fillna(0)
    logger.debug("All gene-level stats merged.")
    return merged
