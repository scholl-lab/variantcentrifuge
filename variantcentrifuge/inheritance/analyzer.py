"""
Main inheritance analyzer orchestrator.

This module coordinates the inheritance pattern analysis by combining
pattern deduction, compound heterozygous analysis, and pattern prioritization.
"""

import logging
import json
from typing import Dict, List, Any, Optional
import pandas as pd
from .deducer import deduce_patterns_for_variant
from .comp_het import analyze_gene_for_compound_het, create_variant_key
from .prioritizer import prioritize_patterns, get_pattern_description
from .segregation_checker import calculate_segregation_score

logger = logging.getLogger(__name__)


def analyze_inheritance(
    df: pd.DataFrame, pedigree_data: Dict[str, Dict[str, Any]], sample_list: List[str]
) -> pd.DataFrame:
    """
    Orchestrates inheritance analysis and adds results to the DataFrame.

    This function performs three main passes:
    1. Per-variant pattern deduction
    2. Compound heterozygous analysis (per gene)
    3. Pattern prioritization and finalization

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing variant data with sample genotypes
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information dictionary
    sample_list : List[str]
        List of sample IDs to analyze

    Returns
    -------
    pd.DataFrame
        DataFrame with added inheritance columns:
        - Inheritance_Pattern: The highest priority pattern
        - Inheritance_Details: JSON string with detailed information
    """
    if df.empty:
        logger.warning("Empty DataFrame provided for inheritance analysis")
        df["Inheritance_Pattern"] = "none"
        df["Inheritance_Details"] = "{}"
        return df

    # Handle empty pedigree data differently - allow analysis to proceed
    if not pedigree_data:
        pedigree_data = {}
        logger.info("No pedigree data provided - treating samples as unrelated individuals")

    logger.info(
        f"Starting inheritance analysis for {len(df)} variants across {len(sample_list)} samples"
    )

    # Initialize new columns
    df["_inheritance_patterns"] = None
    df["_comp_het_info"] = None
    df["Inheritance_Pattern"] = "none"
    df["Inheritance_Details"] = "{}"

    # Pass 1: Per-Variant Pattern Deduction
    logger.info("Pass 1: Deducing inheritance patterns per variant")
    df["_inheritance_patterns"] = df.apply(
        lambda row: deduce_patterns_for_variant(row.to_dict(), pedigree_data, sample_list), axis=1
    )

    # Pass 2: Compound Heterozygous Analysis
    logger.info("Pass 2: Analyzing compound heterozygous patterns")

    # Group by gene and analyze
    comp_het_results_by_gene = {}
    if "GENE" in df.columns:
        gene_counts = df.groupby("GENE").size()
        logger.debug(f"Genes with multiple variants: {gene_counts[gene_counts > 1].to_dict()}")

        for gene, gene_df in df.groupby("GENE"):
            if pd.isna(gene) or gene == "":
                continue

            if len(gene_df) > 1:
                logger.debug(f"Analyzing gene {gene} with {len(gene_df)} variants for compound het")

            comp_het_results = analyze_gene_for_compound_het(gene_df, pedigree_data, sample_list)

            if comp_het_results:
                comp_het_results_by_gene[gene] = comp_het_results
                logger.debug(
                    f"Found compound het patterns in gene {gene}: {len(comp_het_results)} variants"
                )

    # Apply compound het results to DataFrame
    for idx, row in df.iterrows():
        variant_key = create_variant_key(row)
        gene = row.get("GENE", "")

        if gene in comp_het_results_by_gene and variant_key in comp_het_results_by_gene[gene]:
            df.at[idx, "_comp_het_info"] = comp_het_results_by_gene[gene][variant_key]

    # Pass 3: Prioritize and Finalize
    logger.info("Pass 3: Prioritizing patterns and creating final output")

    for idx, row in df.iterrows():
        # Get all patterns including compound het
        all_patterns = list(row["_inheritance_patterns"] or [])

        # Add compound het pattern if applicable
        comp_het_info = row["_comp_het_info"]
        comp_het_patterns = []
        if comp_het_info:
            for sample_id, info in comp_het_info.items():
                if info.get("is_compound_het"):
                    comp_het_type = info.get("comp_het_type", "compound_heterozygous")
                    if comp_het_type != "not_compound_heterozygous":
                        comp_het_patterns.append(comp_het_type)

        # Add unique compound het patterns
        for pattern in comp_het_patterns:
            if pattern not in all_patterns:
                all_patterns.append(pattern)

        # Prepare variant info for prioritization
        variant_info = prepare_variant_info(row)

        # Calculate segregation scores for all patterns if we have family data
        segregation_results = None
        if pedigree_data and len(pedigree_data) > 1:
            segregation_results = calculate_segregation_score(
                all_patterns, row.to_dict(), pedigree_data, sample_list
            )

        # Get the best pattern with segregation consideration
        best_pattern, confidence = prioritize_patterns(
            all_patterns, variant_info, None, segregation_results
        )

        # Create detailed inheritance information
        details = create_inheritance_details(
            row, best_pattern, all_patterns, confidence, comp_het_info, pedigree_data, sample_list
        )

        # Set final values
        df.at[idx, "Inheritance_Pattern"] = best_pattern
        df.at[idx, "Inheritance_Details"] = json.dumps(details)

    # Clean up temporary columns
    df = df.drop(columns=["_inheritance_patterns", "_comp_het_info"])

    # Log summary statistics
    pattern_counts = df["Inheritance_Pattern"].value_counts()
    logger.info(f"Inheritance analysis complete. Pattern distribution: {pattern_counts.to_dict()}")

    return df


def prepare_variant_info(row: pd.Series) -> Dict[str, Any]:
    """
    Prepare variant information for pattern prioritization.

    Parameters
    ----------
    row : pd.Series
        Variant row from DataFrame

    Returns
    -------
    Dict[str, Any]
        Dictionary with variant properties
    """
    variant_info = {}

    # Check rarity based on population frequencies
    af_fields = ["AF", "gnomAD_AF", "MAX_AF", "gnomADg_AF_POPMAX"]
    min_af = 1.0

    for field in af_fields:
        if field in row and pd.notna(row[field]):
            try:
                af = float(row[field])
                min_af = min(min_af, af)
            except (ValueError, TypeError):
                pass

    variant_info["allele_frequency"] = min_af
    variant_info["is_rare"] = min_af < 0.01
    variant_info["is_very_rare"] = min_af < 0.001

    # Check deleteriousness
    impact = row.get("IMPACT", "").upper()
    variant_info["is_deleterious"] = impact in ["HIGH", "MODERATE"]
    variant_info["is_lof"] = impact == "HIGH"

    # Add clinical significance if available
    clin_sig = row.get("CLIN_SIG", "")
    variant_info["is_pathogenic"] = "pathogenic" in str(clin_sig).lower()

    return variant_info


def create_inheritance_details(
    row: pd.Series,
    best_pattern: str,
    all_patterns: List[str],
    confidence: float,
    comp_het_info: Optional[Dict[str, Any]],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Dict[str, Any]:
    """
    Create detailed inheritance information dictionary.

    Parameters
    ----------
    row : pd.Series
        Variant row
    best_pattern : str
        The prioritized pattern
    all_patterns : List[str]
        All possible patterns
    confidence : float
        Confidence score
    comp_het_info : Optional[Dict[str, Any]]
        Compound het information if applicable
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree data
    sample_list : List[str]
        Sample list

    Returns
    -------
    Dict[str, Any]
        Dictionary with detailed inheritance information
    """
    details = {
        "primary_pattern": best_pattern,
        "all_patterns": list(set(all_patterns)),
        "confidence": confidence,
        "pattern_description": get_pattern_description(best_pattern),
        "samples_with_pattern": [],
    }

    # Add sample-specific information
    for sample_id in sample_list:
        if sample_id not in row:
            continue

        gt = str(row[sample_id])
        if gt == "./." or pd.isna(gt):
            continue

        sample_info = {
            "sample_id": sample_id,
            "genotype": gt,
            "affected": pedigree_data.get(sample_id, {}).get("affected_status") == "2",
        }

        # Add compound het partner if applicable
        if comp_het_info and sample_id in comp_het_info:
            ch_info = comp_het_info[sample_id]
            if ch_info.get("is_compound_het"):
                sample_info["compound_het_partner"] = ch_info.get("partner_variant")
                sample_info["compound_het_gene"] = ch_info.get("gene")
                sample_info["compound_het_configuration"] = ch_info.get("inheritance_type")

        details["samples_with_pattern"].append(sample_info)

    # Add summary counts
    details["affected_count"] = sum(1 for s in details["samples_with_pattern"] if s["affected"])
    details["carrier_count"] = sum(
        1
        for s in details["samples_with_pattern"]
        if s["genotype"] in ["0/1", "1/0"] and not s["affected"]
    )

    return details


def get_inheritance_summary(df: pd.DataFrame) -> Dict[str, Any]:
    """
    Generate a summary of inheritance analysis results.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with inheritance analysis results

    Returns
    -------
    Dict[str, Any]
        Summary dictionary
    """
    summary = {
        "total_variants": len(df),
        "pattern_counts": df["Inheritance_Pattern"].value_counts().to_dict(),
        "high_confidence_patterns": 0,
        "compound_het_genes": set(),
        "de_novo_variants": 0,
    }

    # Analyze details
    for _, row in df.iterrows():
        try:
            details = json.loads(row["Inheritance_Details"])

            # Count high confidence patterns
            if details.get("confidence", 0) > 0.8:
                summary["high_confidence_patterns"] += 1

            # Track compound het genes
            for sample in details.get("samples_with_pattern", []):
                if "compound_het_gene" in sample:
                    summary["compound_het_genes"].add(sample["compound_het_gene"])

            # Count de novo
            if row["Inheritance_Pattern"] == "de_novo":
                summary["de_novo_variants"] += 1

        except (json.JSONDecodeError, KeyError):
            continue

    summary["compound_het_genes"] = list(summary["compound_het_genes"])
    summary["compound_het_gene_count"] = len(summary["compound_het_genes"])

    return summary


def filter_by_inheritance_pattern(
    df: pd.DataFrame, patterns: List[str], min_confidence: Optional[float] = None
) -> pd.DataFrame:
    """
    Filter DataFrame by inheritance patterns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with inheritance analysis
    patterns : List[str]
        List of patterns to include
    min_confidence : Optional[float]
        Minimum confidence threshold

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame
    """
    # Filter by pattern
    mask = df["Inheritance_Pattern"].isin(patterns)

    # Optionally filter by confidence
    if min_confidence is not None:
        confidence_mask = df["Inheritance_Details"].apply(
            lambda x: json.loads(x).get("confidence", 0) >= min_confidence
        )
        mask = mask & confidence_mask

    return df[mask].copy()


def export_inheritance_report(
    df: pd.DataFrame, output_path: str, sample_list: Optional[List[str]] = None
) -> None:
    """
    Export a detailed inheritance report.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with inheritance analysis
    output_path : str
        Path for output file
    sample_list : Optional[List[str]]
        Optional subset of samples to include
    """
    report_data = []

    for _, row in df.iterrows():
        variant_info = {
            "chromosome": row.get("CHROM", ""),
            "position": row.get("POS", ""),
            "reference": row.get("REF", ""),
            "alternate": row.get("ALT", ""),
            "gene": row.get("GENE", ""),
            "impact": row.get("IMPACT", ""),
            "inheritance_pattern": row["Inheritance_Pattern"],
        }

        try:
            details = json.loads(row["Inheritance_Details"])
            variant_info.update(details)
        except (json.JSONDecodeError, KeyError):
            pass

        report_data.append(variant_info)

    # Save as JSON
    with open(output_path, "w") as f:
        json.dump(report_data, f, indent=2)

    logger.info(f"Inheritance report saved to {output_path}")
