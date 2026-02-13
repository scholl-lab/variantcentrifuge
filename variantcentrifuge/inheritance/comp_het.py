"""
Compound heterozygous analyzer for inheritance analysis.

This module identifies compound heterozygous variants where an individual
has two different heterozygous variants in the same gene, typically one
from each parent.
"""

import logging
from typing import Any

import pandas as pd

from ..genotype_utils import is_het, is_missing, is_ref, is_variant
from ..ped_reader import get_parents, is_affected

logger = logging.getLogger(__name__)


def analyze_gene_for_compound_het(
    gene_df: pd.DataFrame, pedigree_data: dict[str, dict[str, Any]], sample_list: list[str]
) -> dict[str, dict[str, Any]]:
    """
    Analyze a gene for compound heterozygous patterns.

    Parameters
    ----------
    gene_df : pd.DataFrame
        DataFrame containing all variants in a single gene
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information
    sample_list : List[str]
        List of sample IDs to analyze

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Dictionary mapping variant keys to compound het details
    """
    comp_het_results: dict[str, dict[str, Any]] = {}

    # Skip if too few variants
    if len(gene_df) < 2:
        return comp_het_results

    # Deduplicate by genomic position to prevent false positives from
    # split-snpeff-lines creating multiple rows per variant (one per transcript)
    gene_df_unique = gene_df.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"])
    if len(gene_df_unique) < 2:
        return comp_het_results

    # Get gene name
    gene_name = gene_df_unique.iloc[0].get("GENE", "Unknown")

    # Analyze each sample using deduplicated variants
    for sample_id in sample_list:
        # For samples without pedigree data, still check if affected
        if sample_id not in pedigree_data:
            # Create minimal pedigree entry for single samples
            pedigree_data[sample_id] = {
                "sample_id": sample_id,
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "2",  # Assume affected if no pedigree
            }

        # Find compound het pairs for this sample (using deduplicated DataFrame)
        comp_het_pairs = find_compound_het_pairs(
            sample_id, gene_df_unique, pedigree_data, sample_list
        )

        # Store results for each variant in a pair
        for var1_idx, var2_idx in comp_het_pairs:
            var1_key = create_variant_key(gene_df_unique.iloc[var1_idx])
            var2_key = create_variant_key(gene_df_unique.iloc[var2_idx])

            # Determine compound het type
            trans_config, comp_het_type = determine_compound_het_type(
                sample_id, var1_idx, var2_idx, gene_df_unique, pedigree_data, sample_list
            )

            # Store compound het info for both variants
            comp_het_info = {
                "sample_id": sample_id,
                "gene": gene_name,
                "partner_variant": var2_key,
                "is_compound_het": True,
                "comp_het_type": comp_het_type,
                "inheritance_type": trans_config,
            }

            if var1_key not in comp_het_results:
                comp_het_results[var1_key] = {}
            comp_het_results[var1_key][sample_id] = comp_het_info

            # Store reciprocal information
            comp_het_info_reciprocal = comp_het_info.copy()
            comp_het_info_reciprocal["partner_variant"] = var1_key

            if var2_key not in comp_het_results:
                comp_het_results[var2_key] = {}
            comp_het_results[var2_key][sample_id] = comp_het_info_reciprocal

    return comp_het_results


def find_compound_het_pairs(
    sample_id: str,
    gene_df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
) -> list[tuple[int, int]]:
    """
    Find all compound heterozygous variant pairs for a sample in a gene.

    Parameters
    ----------
    sample_id : str
        The sample to analyze
    gene_df : pd.DataFrame
        DataFrame with variants in the gene
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information
    sample_list : List[str]
        Available samples

    Returns
    -------
    List[Tuple[int, int]]
        List of tuples (variant1_index, variant2_index) representing comp het pairs
    """
    pairs: list[tuple[int, int]] = []

    # Get indices of heterozygous variants for this sample
    het_indices = []
    for pos_idx, (idx, row) in enumerate(gene_df.iterrows()):
        if sample_id in row and is_het(str(row[sample_id])):
            het_indices.append(pos_idx)

    # Need at least 2 heterozygous variants
    if len(het_indices) < 2:
        return pairs

    # Check all pairs
    for i in range(len(het_indices)):
        for j in range(i + 1, len(het_indices)):
            idx1, idx2 = het_indices[i], het_indices[j]

            # Check if this could be a compound het pair
            if is_potential_compound_het(
                sample_id, idx1, idx2, gene_df, pedigree_data, sample_list
            ):
                pairs.append((idx1, idx2))

    return pairs


def is_potential_compound_het(
    sample_id: str,
    var1_idx: int,
    var2_idx: int,
    gene_df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
) -> bool:
    """
    Check if two variants could form a compound heterozygous pair.

    Parameters
    ----------
    sample_id : str
        The sample to check
    var1_idx : int
        Index of first variant
    var2_idx : int
        Index of second variant
    gene_df : pd.DataFrame
        DataFrame with gene variants
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information
    sample_list : List[str]
        Available samples

    Returns
    -------
    bool
        True if variants could be compound heterozygous
    """
    var1 = gene_df.iloc[var1_idx]
    var2 = gene_df.iloc[var2_idx]

    # Both must be heterozygous in the sample
    if not (is_het(str(var1.get(sample_id, "./."))) and is_het(str(var2.get(sample_id, "./.")))):
        return False

    # Get parents
    father_id, mother_id = get_parents(sample_id, pedigree_data)

    # If we have parent data, check for trans configuration
    if father_id and mother_id and father_id in sample_list and mother_id in sample_list:
        # Check if variants come from different parents (trans)
        father_var1_gt = str(var1.get(father_id, "./."))
        mother_var1_gt = str(var1.get(mother_id, "./."))
        father_var2_gt = str(var2.get(father_id, "./."))
        mother_var2_gt = str(var2.get(mother_id, "./."))

        # Classic compound het: one variant from each parent
        pattern1 = (
            is_variant(father_var1_gt)
            and is_ref(mother_var1_gt)
            and is_ref(father_var2_gt)
            and is_variant(mother_var2_gt)
        )

        pattern2 = (
            is_ref(father_var1_gt)
            and is_variant(mother_var1_gt)
            and is_variant(father_var2_gt)
            and is_ref(mother_var2_gt)
        )

        if pattern1 or pattern2:
            return True

    # If no parent data or unclear inheritance, still consider as potential
    # if the sample is affected
    return is_affected(sample_id, pedigree_data)


def determine_compound_het_type(
    sample_id: str,
    var1_idx: int,
    var2_idx: int,
    gene_df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
) -> tuple[str, str]:
    """
    Determine the type of compound heterozygous pattern.

    Parameters
    ----------
    sample_id : str
        The sample to analyze
    var1_idx : int
        Index of first variant
    var2_idx : int
        Index of second variant
    gene_df : pd.DataFrame
        DataFrame with gene variants
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information
    sample_list : List[str]
        Available samples

    Returns
    -------
    Tuple[str, str]
        Tuple of (inheritance_type, comp_het_type)
        inheritance_type: "trans", "cis", "unknown", "ambiguous"
        comp_het_type: "compound_heterozygous", "compound_heterozygous_possible", etc.
    """
    var1 = gene_df.iloc[var1_idx]
    var2 = gene_df.iloc[var2_idx]

    father_id, mother_id = get_parents(sample_id, pedigree_data)

    # No pedigree data
    if not pedigree_data or (not father_id and not mother_id):
        return ("unknown", "compound_heterozygous_possible_no_pedigree")

    # Missing parents
    if not (father_id and mother_id):
        return ("unknown", "compound_heterozygous_possible_missing_parents")

    # Parents not in sample list
    if father_id not in sample_list or mother_id not in sample_list:
        return ("unknown", "compound_heterozygous_possible_missing_parents")

    # Get genotypes
    father_var1_gt = str(var1.get(father_id, "./."))
    mother_var1_gt = str(var1.get(mother_id, "./."))
    father_var2_gt = str(var2.get(father_id, "./."))
    mother_var2_gt = str(var2.get(mother_id, "./."))

    # Check for missing genotypes
    missing_genotypes = (
        is_missing(father_var1_gt)
        or is_missing(mother_var1_gt)
        or is_missing(father_var2_gt)
        or is_missing(mother_var2_gt)
    )

    if missing_genotypes:
        return ("unknown", "compound_heterozygous_possible_missing_parent_genotypes")

    # Determine parent of origin for each variant
    var1_from_father = is_variant(father_var1_gt)
    var1_from_mother = is_variant(mother_var1_gt)
    var2_from_father = is_variant(father_var2_gt)
    var2_from_mother = is_variant(mother_var2_gt)

    # Clear trans pattern
    if (
        var1_from_father and not var1_from_mother and not var2_from_father and var2_from_mother
    ) or (not var1_from_father and var1_from_mother and var2_from_father and not var2_from_mother):
        return ("trans", "compound_heterozygous")

    # Clear cis pattern (both from same parent)
    if (
        var1_from_father and var2_from_father and not var1_from_mother and not var2_from_mother
    ) or (var1_from_mother and var2_from_mother and not var1_from_father and not var2_from_father):
        return ("cis", "not_compound_heterozygous")

    # Ambiguous case - both parents have both variants
    if var1_from_father and var1_from_mother and var2_from_father and var2_from_mother:
        return ("ambiguous", "compound_heterozygous_possible")

    # One variant has ambiguous origin
    if (var1_from_father and var1_from_mother) or (var2_from_father and var2_from_mother):
        return ("ambiguous", "compound_heterozygous_possible")

    # Default to possible
    return ("unknown", "compound_heterozygous_possible")


def create_variant_key(variant_row: pd.Series) -> str:
    """
    Create a unique key for a variant.

    Parameters
    ----------
    variant_row : pd.Series
        Series containing variant information

    Returns
    -------
    str
        Unique variant identifier
    """
    chrom = variant_row.get("CHROM", "chr?")
    pos = variant_row.get("POS", "0")
    ref = variant_row.get("REF", "N")
    alt = variant_row.get("ALT", "N")

    return f"{chrom}:{pos}:{ref}>{alt}"


def get_compound_het_summary(
    comp_het_results: dict[str, dict[str, Any]], sample_id: str
) -> dict[str, Any]:
    """
    Get a summary of compound heterozygous findings for a sample.

    Parameters
    ----------
    comp_het_results : Dict[str, Dict[str, Any]]
        Complete compound het analysis results
    sample_id : str
        Sample to summarize

    Returns
    -------
    Dict[str, Any]
        Summary dictionary
    """
    summary: dict[str, Any] = {"total_comp_het_variants": 0, "comp_het_pairs": [], "genes_with_comp_het": set()}

    processed_pairs = set()

    for variant_key, sample_data in comp_het_results.items():
        if sample_id in sample_data:
            info = sample_data[sample_id]
            if info["is_compound_het"]:
                summary["total_comp_het_variants"] += 1
                summary["genes_with_comp_het"].add(info["gene"])

                # Track unique pairs
                pair_key = tuple(sorted([variant_key, info["partner_variant"]]))
                if pair_key not in processed_pairs:
                    processed_pairs.add(pair_key)
                    summary["comp_het_pairs"].append(
                        {
                            "variant1": pair_key[0],
                            "variant2": pair_key[1],
                            "gene": info["gene"],
                            "configuration": info["inheritance_type"],
                        }
                    )

    summary["genes_with_comp_het"] = list(summary["genes_with_comp_het"])
    return summary


def filter_high_confidence_compound_hets(
    comp_het_results: dict[str, dict[str, Any]], min_quality: float = 30.0
) -> dict[str, dict[str, Any]]:
    """
    Filter compound het results to keep only high-confidence calls.

    Parameters
    ----------
    comp_het_results : Dict[str, Dict[str, Any]]
        All compound het results
    min_quality : float
        Minimum quality score

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Filtered results
    """
    # This is a placeholder for quality filtering
    # In practice, you would check variant quality scores
    return comp_het_results
