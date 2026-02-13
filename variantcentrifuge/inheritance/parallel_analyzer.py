"""
Parallel inheritance analyzer for improved performance.

This module provides a parallel implementation of inheritance analysis
that processes genes concurrently for compound heterozygous detection.
"""

import json
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any

import pandas as pd

from .analyzer import create_inheritance_details
from .comp_het import analyze_gene_for_compound_het, create_variant_key
from .deducer import deduce_patterns_for_variant
from .prioritizer import prioritize_patterns
from .segregation_checker import calculate_segregation_score

try:
    from .comp_het_vectorized import analyze_gene_for_compound_het_vectorized

    VECTORIZED_AVAILABLE = True
except ImportError:
    VECTORIZED_AVAILABLE = False

logger = logging.getLogger(__name__)


def _process_gene_group(
    gene: str,
    gene_df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
    use_vectorized: bool = True,
) -> tuple[str, dict[str, Any]]:
    """
    Process a single gene's variants for compound heterozygous patterns.

    This function is designed to be run in a worker process.

    Parameters
    ----------
    gene : str
        Gene name
    gene_df : pd.DataFrame
        DataFrame containing variants for this gene
    pedigree_data : Dict
        Pedigree information
    sample_list : List[str]
        List of sample IDs
    use_vectorized : bool
        Whether to use vectorized implementation

    Returns
    -------
    Tuple[str, Dict[str, Any]]
        Gene name and compound het results
    """
    if pd.isna(gene) or gene == "":
        return gene, {}

    if len(gene_df) <= 1:
        return gene, {}

    # Use vectorized implementation if available and requested
    if use_vectorized and VECTORIZED_AVAILABLE:
        comp_het_results = analyze_gene_for_compound_het_vectorized(
            gene_df, pedigree_data, sample_list
        )
    else:
        comp_het_results = analyze_gene_for_compound_het(gene_df, pedigree_data, sample_list)

    return gene, comp_het_results


def analyze_inheritance_parallel(
    df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
    use_vectorized_comp_het: bool = True,
    n_workers: int | None = None,
    min_variants_for_parallel: int = 100,
) -> pd.DataFrame:
    """
    Parallel version of inheritance analysis with concurrent gene processing.

    This function performs the same analysis as analyze_inheritance but
    processes genes in parallel for improved performance on multi-core systems.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing variant data with sample genotypes
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information dictionary
    sample_list : List[str]
        List of sample IDs to analyze
    use_vectorized_comp_het : bool, optional
        Whether to use the vectorized compound het implementation
    n_workers : int, optional
        Number of parallel workers. If None, uses CPU count
    min_variants_for_parallel : int, optional
        Minimum number of variants to use parallel processing

    Returns
    -------
    pd.DataFrame
        DataFrame with added inheritance columns
    """
    import time

    start_time = time.time()

    if df.empty:
        logger.warning("Empty DataFrame provided for inheritance analysis")
        df["Inheritance_Pattern"] = "none"
        df["Inheritance_Details"] = "{}"
        return df

    # Handle empty pedigree data
    if not pedigree_data:
        pedigree_data = {}
        logger.info("No pedigree data provided - treating samples as unrelated individuals")

    logger.info(
        f"Starting parallel inheritance analysis for {len(df)} variants "
        f"across {len(sample_list)} samples"
    )

    # Create a copy to avoid fragmentation
    df = df.copy()
    df["_inheritance_patterns"] = None
    df["_comp_het_info"] = None
    df["Inheritance_Pattern"] = "none"
    df["Inheritance_Details"] = "{}"

    # Track timing for each pass
    pass_times = {}

    # Pass 1: Per-Variant Pattern Deduction
    pass1_start = time.time()
    logger.info("Pass 1: Deducing inheritance patterns per variant")
    df["_inheritance_patterns"] = df.apply(
        lambda row: deduce_patterns_for_variant(row.to_dict(), pedigree_data, sample_list), axis=1
    )
    pass_times["pattern_deduction"] = time.time() - pass1_start

    # Pass 2: Compound Heterozygous Analysis (PARALLEL)
    pass2_start = time.time()
    implementation = (
        "vectorized" if (use_vectorized_comp_het and VECTORIZED_AVAILABLE) else "original"
    )

    comp_het_results_by_gene = {}

    # Check if we should use parallel processing
    use_parallel = (
        (n_workers is None or n_workers > 1)
        and len(df) >= min_variants_for_parallel
        and "GENE" in df.columns
    )

    if use_parallel and "GENE" in df.columns:
        # Group by gene
        gene_groups = df.groupby("GENE")
        genes_with_multiple_variants = [
            (gene, group_df)
            for gene, group_df in gene_groups
            if len(group_df) > 1 and not pd.isna(gene) and gene != ""
        ]

        num_genes = len(genes_with_multiple_variants)
        if num_genes > 0:
            logger.info(
                f"Pass 2: Analyzing {num_genes} genes for compound heterozygous patterns "
                f"in parallel (using {implementation} implementation)"
            )

            # Process genes in parallel
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                # Submit all gene processing tasks
                future_to_gene = {
                    executor.submit(
                        _process_gene_group,
                        gene,
                        gene_df,
                        pedigree_data,
                        sample_list,
                        use_vectorized_comp_het,
                    ): gene
                    for gene, gene_df in genes_with_multiple_variants
                }

                # Collect results as they complete
                completed = 0
                for future in as_completed(future_to_gene):
                    gene = future_to_gene[future]
                    try:
                        gene_name, comp_het_results = future.result()
                        if comp_het_results:
                            comp_het_results_by_gene[gene_name] = comp_het_results
                            logger.debug(
                                f"Found compound het patterns in gene {gene_name}: "
                                f"{len(comp_het_results)} variants"
                            )
                        completed += 1
                        if completed % 100 == 0:
                            logger.info(f"Processed {completed}/{num_genes} genes")
                    except Exception as e:
                        logger.error(f"Error processing gene {gene}: {e}")
    else:
        # Fall back to sequential processing
        logger.info(
            f"Pass 2: Analyzing compound heterozygous patterns sequentially "
            f"(using {implementation} implementation)"
        )

        if "GENE" in df.columns:
            for gene, gene_df in df.groupby("GENE"):
                if pd.isna(gene) or gene == "" or len(gene_df) <= 1:
                    continue

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

    pass_times["compound_het_analysis"] = time.time() - pass2_start

    # Apply compound het results to DataFrame
    apply_start = time.time()
    for idx, row in df.iterrows():
        variant_key = create_variant_key(row)
        gene = row.get("GENE", "")

        if gene in comp_het_results_by_gene and variant_key in comp_het_results_by_gene[gene]:
            df.at[idx, "_comp_het_info"] = comp_het_results_by_gene[gene][variant_key]
    pass_times["apply_comp_het_results"] = time.time() - apply_start

    # Pass 3: Prioritize and Finalize
    pass3_start = time.time()
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

    pass_times["prioritization"] = time.time() - pass3_start

    # Clean up temporary columns
    df = df.drop(columns=["_inheritance_patterns", "_comp_het_info"])

    # Log summary statistics
    total_time = time.time() - start_time
    pattern_counts = df["Inheritance_Pattern"].value_counts()
    logger.info(
        f"Parallel inheritance analysis complete in {total_time:.1f}s. "
        f"Pattern distribution: {pattern_counts.to_dict()}"
    )

    # Log timing breakdown
    logger.debug("Inheritance analysis timing breakdown:")
    for pass_name, duration in pass_times.items():
        pct = (duration / total_time) * 100 if total_time > 0 else 0
        logger.debug(f"  {pass_name}: {duration:.1f}s ({pct:.1f}%)")

    return df
