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
from .comp_het_vectorized import analyze_gene_for_compound_het_vectorized
from .prioritizer import prioritize_patterns
from .segregation_checker import calculate_segregation_score
from .vectorized_deducer import vectorized_deduce_patterns

logger = logging.getLogger(__name__)


def _process_gene_group(
    gene: str,
    gene_df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
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

    Returns
    -------
    Tuple[str, Dict[str, Any]]
        Gene name and compound het results
    """
    if pd.isna(gene) or gene == "":
        return gene, {}

    if len(gene_df) <= 1:
        return gene, {}

    # Use vectorized implementation
    comp_het_results = analyze_gene_for_compound_het_vectorized(gene_df, pedigree_data, sample_list)

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
        Deprecated parameter (kept for backward compatibility).
        Vectorized implementation is now always used.
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

    # Pass 1: Per-Variant Pattern Deduction (vectorized)
    pass1_start = time.time()
    logger.info("Pass 1: Deducing inheritance patterns per variant (vectorized)")
    df["_inheritance_patterns"] = vectorized_deduce_patterns(df, pedigree_data, sample_list)
    pass_times["pattern_deduction"] = time.time() - pass1_start

    # Pass 2: Compound Heterozygous Analysis (PARALLEL)
    pass2_start = time.time()

    comp_het_results_by_gene = {}

    # Check if we should use parallel processing
    from ..memory import ResourceManager

    rm = ResourceManager()

    # Use min_variants_for_parallel if explicitly provided (for backward compatibility),
    # otherwise use ResourceManager's threshold
    if min_variants_for_parallel != 100:  # 100 is the default
        should_parallelize = len(df) >= min_variants_for_parallel
    else:
        should_parallelize = rm.should_parallelize(len(df))

    use_parallel = (
        (n_workers is None or n_workers > 1) and should_parallelize and "GENE" in df.columns
    )

    if use_parallel and "GENE" in df.columns:
        # Group by gene
        gene_groups = df.groupby("GENE", observed=True)
        genes_with_multiple_variants = [
            (gene, group_df)
            for gene, group_df in gene_groups
            if len(group_df) > 1 and not pd.isna(gene) and gene != ""
        ]

        num_genes = len(genes_with_multiple_variants)
        if num_genes > 0:
            # Sort genes by variant count descending (largest first for load balancing)
            genes_with_multiple_variants.sort(key=lambda x: len(x[1]), reverse=True)
            max_gene_size = (
                len(genes_with_multiple_variants[0][1]) if genes_with_multiple_variants else 0
            )

            # Auto-detect worker count if not specified
            if n_workers is None:
                memory_per_gene_gb = rm.estimate_memory(max_gene_size, len(sample_list))
                n_workers = rm.auto_workers(
                    task_count=num_genes, memory_per_task_gb=memory_per_gene_gb
                )

            logger.info(
                f"Pass 2: Analyzing {num_genes} genes for compound heterozygous patterns "
                f"in parallel ({n_workers} workers, largest gene: {max_gene_size} variants)"
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
        logger.info("Pass 2: Analyzing compound heterozygous patterns sequentially")

        if "GENE" in df.columns:
            for gene, gene_df in df.groupby("GENE", observed=True):
                if pd.isna(gene) or gene == "" or len(gene_df) <= 1:
                    continue

                try:
                    # Use vectorized implementation
                    comp_het_results = analyze_gene_for_compound_het_vectorized(
                        gene_df, pedigree_data, sample_list
                    )

                    if comp_het_results:
                        comp_het_results_by_gene[gene] = comp_het_results
                except Exception as e:
                    logger.error(f"Error processing gene {gene}: {e}")

    pass_times["compound_het_analysis"] = time.time() - pass2_start

    # Apply compound het results to DataFrame (vectorized)
    apply_start = time.time()
    if comp_het_results_by_gene:
        # Build variant keys for all rows at once
        variant_keys = (
            df["CHROM"].astype(str)
            + ":"
            + df["POS"].astype(str)
            + ":"
            + df["REF"].astype(str)
            + ">"
            + df["ALT"].astype(str)
        )

        # Get genes array - handle categorical dtype
        if isinstance(df["GENE"].dtype, pd.CategoricalDtype):
            genes = df["GENE"].astype(str).values
        else:
            genes = df["GENE"].values.astype(str)

        # Apply comp het results per gene
        for gene, gene_results in comp_het_results_by_gene.items():
            gene_mask = genes == gene
            if not gene_mask.any():
                continue

            # Get variant keys for this gene
            gene_variant_keys = variant_keys[gene_mask]

            # For each variant in gene_results, find matching rows
            for vk, info in gene_results.items():
                vk_mask = gene_variant_keys == vk
                if vk_mask.any():
                    # Get matching indices in the original DataFrame
                    matching_indices = df.index[gene_mask][vk_mask]
                    for idx in matching_indices:
                        df.at[idx, "_comp_het_info"] = info

    pass_times["apply_comp_het_results"] = time.time() - apply_start

    # Pass 3: Prioritize and Finalize (optimized)
    pass3_start = time.time()
    logger.info("Pass 3: Prioritizing patterns and creating final output")

    # Pre-extract column data to avoid repeated df.at[] lookups
    all_inheritance_patterns = df["_inheritance_patterns"].tolist()
    all_comp_het_info = df["_comp_het_info"].tolist()

    # Pre-compute whether segregation is needed
    needs_segregation = bool(pedigree_data and len(pedigree_data) > 1)

    # Accumulate results in lists for bulk assignment
    inheritance_patterns_result = []
    inheritance_details_result = []

    # Use iterrows since we need Series for create_inheritance_details
    for idx_pos, (_df_idx, row_series) in enumerate(df.iterrows()):
        # Get patterns from pre-extracted list
        patterns_list = list(all_inheritance_patterns[idx_pos] or [])
        comp_info = all_comp_het_info[idx_pos]

        # Add compound het patterns if applicable
        if comp_info:
            for _sample_id, info in comp_info.items():
                if info.get("is_compound_het"):
                    comp_het_type = info.get("comp_het_type", "compound_heterozygous")
                    is_not_negative = comp_het_type != "not_compound_heterozygous"
                    if is_not_negative and comp_het_type not in patterns_list:
                        patterns_list.append(comp_het_type)

        # Calculate segregation scores if needed
        segregation_results = None
        if needs_segregation:
            row_dict = row_series.to_dict()
            segregation_results = calculate_segregation_score(
                patterns_list, row_dict, pedigree_data, sample_list
            )

        # Prioritize patterns
        best_pattern, confidence = prioritize_patterns(patterns_list, segregation_results)

        # Create detailed inheritance information
        details = create_inheritance_details(
            row_series,
            best_pattern,
            patterns_list,
            confidence,
            comp_info,
            pedigree_data,
            sample_list,
            segregation_results,
        )

        # Accumulate results
        inheritance_patterns_result.append(best_pattern)
        inheritance_details_result.append(json.dumps(details))

    # Assign results in bulk (avoids per-row df.at[] overhead)
    df["Inheritance_Pattern"] = inheritance_patterns_result
    df["Inheritance_Details"] = inheritance_details_result

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
