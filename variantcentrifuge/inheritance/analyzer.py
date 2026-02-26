"""
Main inheritance analyzer orchestrator.

This module coordinates the inheritance pattern analysis by combining
pattern deduction, compound heterozygous analysis, and pattern prioritization.
"""

import json
import logging
import time
from typing import Any

import numpy as np
import pandas as pd

from .comp_het_vectorized import (
    analyze_gene_for_compound_het_vectorized,
    build_variant_keys_array,
)
from .prioritizer import get_pattern_description, prioritize_patterns
from .segregation_checker import calculate_segregation_score
from .vectorized_deducer import vectorized_deduce_patterns

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
        Deprecated parameter (kept for backward compatibility).
        Vectorized implementation is now always used.

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

    # Pre-populate pedigree_data for unknown samples (consistent with parallel_analyzer)
    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            pedigree_data[sample_id] = {
                "sample_id": sample_id,
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "2",
            }

    # Create a copy to avoid fragmentation and initialize new columns
    df = df.copy()
    df["_inheritance_patterns"] = None
    df["_comp_het_info"] = None
    df["Inheritance_Pattern"] = "none"
    df["Inheritance_Details"] = "{}"

    # Pass 1: Per-Variant Pattern Deduction (vectorized)
    # Request genotype matrix back so Pass 2 can reuse it (avoids redundant encoding)
    logger.info("Pass 1: Deducing inheritance patterns per variant (vectorized)")
    t_pass1_start = time.monotonic()
    pass1_result = vectorized_deduce_patterns(
        df, pedigree_data, sample_list, return_genotype_matrix=True
    )
    if isinstance(pass1_result, tuple):
        patterns_list_all, gt_matrix, sample_to_idx = pass1_result
    else:
        patterns_list_all = pass1_result
        gt_matrix, sample_to_idx = None, None
    df["_inheritance_patterns"] = patterns_list_all
    t_pass1 = time.monotonic() - t_pass1_start
    logger.info(
        f"Pass 1 completed in {t_pass1:.2f}s ({len(df)} variants, {len(sample_list)} samples)"
    )

    # Pass 2: Compound Heterozygous Analysis
    logger.info("Pass 2: Analyzing compound heterozygous patterns")
    t_pass2_start = time.monotonic()

    # Build positional row index array for slicing the genotype matrix per-gene
    row_positions = np.arange(len(df))

    # Pre-compute variant keys for all rows (reused per-gene and for result application)
    all_variant_keys = build_variant_keys_array(df)

    # Group by gene and analyze
    comp_het_results_by_gene = {}
    if "GENE" in df.columns:
        gene_counts = df.groupby("GENE", observed=True).size()
        logger.debug(f"Genes with multiple variants: {gene_counts[gene_counts > 1].to_dict()}")

        for gene, gene_df in df.groupby("GENE", observed=True):
            if pd.isna(gene) or gene == "":
                continue

            if len(gene_df) > 1:
                logger.debug(f"Analyzing gene {gene} with {len(gene_df)} variants for compound het")

            # Get positional indices for this gene's rows in the full DataFrame
            gene_iloc = df.index.get_indexer(gene_df.index)
            gene_row_idx = row_positions[gene_iloc]

            # Slice pre-computed variant keys for this gene
            gene_vkeys = all_variant_keys[gene_iloc]

            # Use vectorized implementation with pre-built matrix
            comp_het_results = analyze_gene_for_compound_het_vectorized(
                gene_df,
                pedigree_data,
                sample_list,
                gt_matrix=gt_matrix,
                sample_to_idx=sample_to_idx,
                gene_row_indices=gene_row_idx,
                variant_keys=gene_vkeys,
            )

            if comp_het_results:
                comp_het_results_by_gene[gene] = comp_het_results
                logger.debug(
                    f"Found compound het patterns in gene {gene}: {len(comp_het_results)} variants"
                )

    n_genes_with_comphet = len(comp_het_results_by_gene)
    t_pass2_analysis = time.monotonic() - t_pass2_start
    logger.info(
        f"Pass 2 gene analysis completed in {t_pass2_analysis:.2f}s "
        f"({n_genes_with_comphet} genes with compound het patterns)"
    )

    # Apply compound het results to DataFrame — batch assignment
    if comp_het_results_by_gene:
        # Reuse pre-computed variant keys (already a numpy string array)
        variant_keys_series = pd.Series(all_variant_keys, index=df.index)

        # Get genes array - handle categorical dtype
        if isinstance(df["GENE"].dtype, pd.CategoricalDtype):
            genes = df["GENE"].astype(str).values
        else:
            genes = df["GENE"].values.astype(str)

        # Collect all (index, info) pairs, then assign in one pass
        batch_indices: list[int] = []
        batch_values: list[dict] = []

        for gene, gene_results in comp_het_results_by_gene.items():
            gene_mask = genes == gene
            if not np.any(gene_mask):
                continue

            # Get variant keys for this gene
            gene_variant_keys = variant_keys_series[gene_mask]

            # For each variant in gene_results, find matching rows
            for vk, info in gene_results.items():
                vk_mask = gene_variant_keys == vk
                if vk_mask.any():
                    matching_indices = gene_variant_keys.index[vk_mask]
                    for idx in matching_indices:
                        batch_indices.append(idx)
                        batch_values.append(info)

        # Single batch assignment
        if batch_indices:
            for i, idx in enumerate(batch_indices):
                df.at[idx, "_comp_het_info"] = batch_values[i]

    t_pass2_total = time.monotonic() - t_pass2_start
    logger.info(f"Pass 2 total (analysis + apply) completed in {t_pass2_total:.2f}s")

    # Pass 3: Prioritize and Finalize
    logger.info("Pass 3: Prioritizing patterns and creating final output")
    t_pass3_start = time.monotonic()

    _finalize_inheritance_patterns(df, pedigree_data, sample_list)

    t_pass3 = time.monotonic() - t_pass3_start
    logger.info(f"Pass 3 completed in {t_pass3:.2f}s ({len(df)} variants prioritized)")

    # Clean up temporary columns
    df = df.drop(columns=["_inheritance_patterns", "_comp_het_info"])

    # Log summary statistics and timing breakdown
    t_total = t_pass1 + t_pass2_total + t_pass3
    pattern_counts = df["Inheritance_Pattern"].value_counts()
    logger.info(
        f"Inheritance analysis complete in {t_total:.2f}s. "
        f"Timing: Pass1={t_pass1:.2f}s ({t_pass1 / t_total * 100:.0f}%), "
        f"Pass2={t_pass2_total:.2f}s ({t_pass2_total / t_total * 100:.0f}%), "
        f"Pass3={t_pass3:.2f}s ({t_pass3 / t_total * 100:.0f}%). "
        f"Pattern distribution: {pattern_counts.to_dict()}"
    )

    return df


def _finalize_inheritance_patterns(
    df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
) -> None:
    """
    Pass 3: prioritize patterns and create final inheritance output (in-place).

    Replaces iterrows() with column-based iteration to avoid creating
    a full pd.Series per row for wide DataFrames (1000+ columns).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with ``_inheritance_patterns`` and ``_comp_het_info`` columns.
        Modified in-place: ``Inheritance_Pattern`` and ``Inheritance_Details``
        are overwritten.
    pedigree_data : dict
        Pedigree information dictionary.
    sample_list : list[str]
        List of sample IDs to analyze.
    """
    all_inheritance_patterns = df["_inheritance_patterns"].tolist()
    all_comp_het_info = df["_comp_het_info"].tolist()
    needs_segregation = bool(pedigree_data and len(pedigree_data) > 1)

    # Pre-extract only the columns needed for create_inheritance_details:
    # each sample column (for genotype lookup).  Metadata columns are not
    # accessed by the function — it only iterates sample_list.
    sample_cols_present = [s for s in sample_list if s in df.columns]
    sample_col_data: dict[str, list] = {col: df[col].tolist() for col in sample_cols_present}

    n_rows = len(df)
    inheritance_patterns_result: list[str] = []
    inheritance_details_result: list[str] = []

    for idx in range(n_rows):
        patterns_list = list(all_inheritance_patterns[idx] or [])
        comp_info = all_comp_het_info[idx]

        # Add compound het patterns if applicable
        if comp_info:
            for _sample_id, info in comp_info.items():
                if info.get("is_compound_het"):
                    comp_het_type = info.get("comp_het_type", "compound_heterozygous")
                    if (
                        comp_het_type != "not_compound_heterozygous"
                        and comp_het_type not in patterns_list
                    ):
                        patterns_list.append(comp_het_type)

        # Calculate segregation scores if needed
        segregation_results = None
        if needs_segregation:
            row_dict = {col: sample_col_data[col][idx] for col in sample_cols_present}
            segregation_results = calculate_segregation_score(
                patterns_list, row_dict, pedigree_data, sample_list
            )

        best_pattern, confidence = prioritize_patterns(patterns_list, segregation_results)

        # Pass a dict instead of pd.Series — create_inheritance_details
        # only uses `sample_id in row` and `row[sample_id]`, both of
        # which work on dicts.  Eliminates 294k Series allocations.
        row_dict = {col: sample_col_data[col][idx] for col in sample_cols_present}

        details = _create_inheritance_details(
            row_dict,
            best_pattern,
            patterns_list,
            confidence,
            comp_info,
            pedigree_data,
            sample_list,
            segregation_results,
        )

        inheritance_patterns_result.append(best_pattern)
        inheritance_details_result.append(json.dumps(details))

    df["Inheritance_Pattern"] = inheritance_patterns_result
    df["Inheritance_Details"] = inheritance_details_result


def _create_inheritance_details(
    row: dict[str, Any] | pd.Series,
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
    row : dict or pd.Series
        Variant row (only sample genotype columns are accessed via
        ``sample_id in row`` and ``row[sample_id]``)
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

        raw_gt = row[sample_id]
        if raw_gt is None or pd.isna(raw_gt):
            continue
        gt = str(raw_gt)
        if gt in ["./.", "0/0"]:
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
        if s["genotype"] in ("0/1", "1/0", "0|1", "1|0") and not s["affected"]
    )

    return details


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
        for row in df.itertuples(index=True):
            try:
                details = json.loads(getattr(row, "Inheritance_Details", "{}"))

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
                df.at[row.Index, "Inheritance_Confidence"] = f"{confidence:.2f}"
                df.at[row.Index, "Inheritance_Description"] = pattern_desc
                df.at[row.Index, "Inheritance_Samples"] = "; ".join(sample_strings)

            except (json.JSONDecodeError, KeyError, TypeError, AttributeError) as e:
                logger.debug(f"Error parsing inheritance details for row {row.Index}: {e}")
                # Set empty values for this row
                df.at[row.Index, "Inheritance_Confidence"] = ""
                df.at[row.Index, "Inheritance_Description"] = ""
                df.at[row.Index, "Inheritance_Samples"] = ""

        # Drop the original details column
        if "Inheritance_Details" in df.columns:
            df = df.drop(columns=["Inheritance_Details"])

        return df

    else:
        logger.warning(f"Unknown inheritance mode: {mode}. Returning data as-is.")
        return df
