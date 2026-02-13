"""
Main inheritance analyzer orchestrator.

This module coordinates the inheritance pattern analysis by combining
pattern deduction, compound heterozygous analysis, and pattern prioritization.
"""

import json
import logging
from typing import Any

import pandas as pd

from .comp_het import analyze_gene_for_compound_het, create_variant_key
from .deducer import deduce_patterns_for_variant

try:
    from .comp_het_vectorized import analyze_gene_for_compound_het_vectorized

    VECTORIZED_AVAILABLE = True
except ImportError:
    VECTORIZED_AVAILABLE = False
from .prioritizer import get_pattern_description, prioritize_patterns
from .segregation_checker import calculate_segregation_score

logger = logging.getLogger(__name__)


def analyze_inheritance(
    df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
    use_vectorized_comp_het: bool = True,
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
    use_vectorized_comp_het : bool, optional
        Whether to use the vectorized compound het implementation for better performance.
        Defaults to True if available, falls back to original if not.

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

    # Create a copy to avoid fragmentation and initialize new columns
    df = df.copy()
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
    implementation = (
        "vectorized" if (use_vectorized_comp_het and VECTORIZED_AVAILABLE) else "original"
    )
    logger.info(
        f"Pass 2: Analyzing compound heterozygous patterns (using {implementation} implementation)"
    )

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

            # Use vectorized implementation if available and requested
            if use_vectorized_comp_het and VECTORIZED_AVAILABLE:
                comp_het_results = analyze_gene_for_compound_het_vectorized(
                    gene_df, pedigree_data, sample_list
                )
            else:
                comp_het_results = analyze_gene_for_compound_het(
                    gene_df, pedigree_data, sample_list
                )

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
            for _sample_id, info in comp_het_info.items():
                if info.get("is_compound_het"):
                    comp_het_type = info.get("comp_het_type", "compound_heterozygous")
                    if comp_het_type != "not_compound_heterozygous":
                        comp_het_patterns.append(comp_het_type)

        # Add unique compound het patterns
        for pattern in comp_het_patterns:
            if pattern not in all_patterns:
                all_patterns.append(pattern)

        # Calculate segregation scores for all patterns if we have family data
        segregation_results = None
        if pedigree_data and len(pedigree_data) > 1:
            segregation_results = calculate_segregation_score(
                all_patterns, row.to_dict(), pedigree_data, sample_list
            )

        # Get the best pattern with segregation consideration
        best_pattern, confidence = prioritize_patterns(all_patterns, segregation_results)

        # Create detailed inheritance information
        details = create_inheritance_details(
            row,
            best_pattern,
            all_patterns,
            confidence,
            comp_het_info,
            pedigree_data,
            sample_list,
            segregation_results,
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


def create_inheritance_details(
    row: pd.Series,
    best_pattern: str,
    all_patterns: list[str],
    confidence: float,
    comp_het_info: dict[str, Any] | None,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
    segregation_results: dict[str, tuple] | None = None,
) -> dict[str, Any]:
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
    details: dict[str, Any] = {
        "primary_pattern": best_pattern,
        "all_patterns": list(set(all_patterns)),
        "confidence": confidence,
        "pattern_description": get_pattern_description(best_pattern),
        "samples_with_pattern": [],
    }

    # Add segregation p-value if available
    if segregation_results and best_pattern in segregation_results:
        # The 'confidence' from the segregation check is our p-value equivalent.
        details["segregation_p_value"] = segregation_results[best_pattern][1]

    # Add sample-specific information
    for sample_id in sample_list:
        if sample_id not in row:
            continue

        gt = str(row[sample_id])
        if gt in ["./.", "0/0"] or pd.isna(gt):
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
                # Handle both old single partner and new multiple partners format
                if "partner_variants" in ch_info:
                    sample_info["compound_het_partners"] = ch_info["partner_variants"]
                elif "partner_variant" in ch_info:
                    # Backward compatibility - single partner
                    sample_info["compound_het_partner"] = ch_info["partner_variant"]
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


def get_inheritance_summary(df: pd.DataFrame) -> dict[str, Any]:
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
    df: pd.DataFrame, patterns: list[str], min_confidence: float | None = None
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
    df: pd.DataFrame, output_path: str, sample_list: list[str] | None = None
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


def process_inheritance_output(
    df: pd.DataFrame, mode: str, preserve_details_for_scoring: bool = False
) -> pd.DataFrame:
    """
    Process inheritance analysis output based on the selected mode.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with inheritance analysis results
    mode : str
        Output mode: 'simple', 'columns', or 'full'
    preserve_details_for_scoring : bool, default False
        If True, preserve Inheritance_Details even in simple mode for scoring

    Returns
    -------
    pd.DataFrame
        Processed DataFrame based on the selected mode
    """
    # If no inheritance columns exist, return as-is
    if "Inheritance_Pattern" not in df.columns or "Inheritance_Details" not in df.columns:
        logger.debug("No inheritance columns found in DataFrame")
        return df

    if mode == "full":
        # Return as-is with full JSON details
        return df

    elif mode == "simple":
        # Drop the details column unless needed for scoring
        if "Inheritance_Details" in df.columns and not preserve_details_for_scoring:
            logger.debug("Dropping Inheritance_Details column in simple mode")
            df = df.drop(columns=["Inheritance_Details"])
        elif preserve_details_for_scoring:
            logger.debug("Preserving Inheritance_Details column for scoring configuration")
        return df

    elif mode == "columns":
        # Create new columns from the JSON details
        for idx, row in df.iterrows():
            try:
                details = json.loads(row["Inheritance_Details"])

                # Extract key information
                confidence = details.get("confidence", 0)
                pattern_desc = details.get("pattern_description", "")
                samples_info = details.get("samples_with_pattern", [])

                # Format samples with pattern
                sample_strings = []
                for sample in samples_info:
                    sample_id = sample.get("sample_id", "")
                    genotype = sample.get("genotype", "")

                    # Basic format: sample_id(genotype)
                    sample_str = f"{sample_id}({genotype})"

                    # Add compound het partner(s) if present
                    if "compound_het_partners" in sample:
                        # New format - list of partners
                        partners = sample["compound_het_partners"]
                        if partners:
                            if len(partners) == 1:
                                sample_str += f", partner:{partners[0]}"
                            else:
                                partners_str = ", ".join(partners)
                                sample_str += f", partners:[{partners_str}]"
                    elif "compound_het_partner" in sample:
                        # Old format - single partner (backward compatibility)
                        partner = sample["compound_het_partner"]
                        sample_str += f", partner:{partner}"

                    sample_strings.append(sample_str)

                # Set the values
                df.at[idx, "Inheritance_Confidence"] = f"{confidence:.2f}"
                df.at[idx, "Inheritance_Description"] = pattern_desc
                df.at[idx, "Inheritance_Samples"] = "; ".join(sample_strings)

            except (json.JSONDecodeError, KeyError, TypeError) as e:
                logger.debug(f"Error parsing inheritance details for row {idx}: {e}")
                # Set empty values for this row
                df.at[idx, "Inheritance_Confidence"] = ""
                df.at[idx, "Inheritance_Description"] = ""
                df.at[idx, "Inheritance_Samples"] = ""

        # Drop the original details column
        if "Inheritance_Details" in df.columns:
            df = df.drop(columns=["Inheritance_Details"])

        return df

    else:
        logger.warning(f"Unknown inheritance mode: {mode}. Returning data as-is.")
        return df
