"""
Parallel inheritance analyzer for improved performance.

This module provides a parallel implementation of inheritance analysis
that processes genes concurrently for compound heterozygous detection.
"""

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

import numpy as np
import pandas as pd

from .analyzer import _finalize_inheritance_patterns
from .comp_het_vectorized import (
    _process_gene_group_arrays,
    analyze_gene_for_compound_het_vectorized,
    build_variant_keys_array,
)
from .vectorized_deducer import vectorized_deduce_patterns

logger = logging.getLogger(__name__)


def build_pedigree_arrays(
    sample_list: list[str],
    pedigree_data: dict[str, dict[str, Any]],
    sample_to_idx: dict[str, int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Pre-compute pedigree arrays for fast worker-side lookups.

    Builds integer arrays indexed by sample position (same order as sample_list)
    so workers can do O(1) array lookups instead of dict.get() calls.

    Parameters
    ----------
    sample_list : list[str]
        Ordered list of sample IDs (defines array positions).
    pedigree_data : dict[str, dict[str, Any]]
        Pedigree information keyed by sample ID.
    sample_to_idx : dict[str, int]
        Mapping from sample ID to column index in the genotype matrix.

    Returns
    -------
    father_idx_arr : np.ndarray (int32)
        Index of each sample's father in sample_list (-1 if none/missing).
    mother_idx_arr : np.ndarray (int32)
        Index of each sample's mother in sample_list (-1 if none/missing).
    affected_arr : np.ndarray (int8)
        1 if sample is affected (status "2"), else 0.
    """
    n = len(sample_list)
    father_idx_arr = np.full(n, -1, dtype=np.int32)
    mother_idx_arr = np.full(n, -1, dtype=np.int32)
    affected_arr = np.zeros(n, dtype=np.int8)

    for i, sample_id in enumerate(sample_list):
        ped = pedigree_data.get(sample_id, {})

        # Affected status
        if ped.get("affected_status") == "2":
            affected_arr[i] = 1

        # Father
        father_id = ped.get("father_id", "0")
        if father_id and father_id != "0" and father_id in sample_to_idx:
            father_idx_arr[i] = sample_to_idx[father_id]

        # Mother
        mother_id = ped.get("mother_id", "0")
        if mother_id and mother_id != "0" and mother_id in sample_to_idx:
            mother_idx_arr[i] = sample_to_idx[mother_id]

    return father_idx_arr, mother_idx_arr, affected_arr


def _get_batch_size(n_workers: int) -> int:
    """
    Determine the batch size for parallel gene processing.

    Checks ``VARIANTCENTRIFUGE_BATCH_SIZE`` environment variable first.
    Falls back to ``min(4 * n_workers, 500)`` if not set or invalid.

    Parameters
    ----------
    n_workers : int
        Number of parallel workers.

    Returns
    -------
    int
        Batch size (number of genes per submission batch).
    """
    env_val = os.environ.get("VARIANTCENTRIFUGE_BATCH_SIZE")
    if env_val is not None:
        try:
            batch_size = int(env_val)
            if batch_size >= 1:
                return batch_size
            logger.warning(
                f"VARIANTCENTRIFUGE_BATCH_SIZE={env_val!r} must be >= 1; "
                "falling back to auto batch size"
            )
        except ValueError:
            logger.warning(
                f"VARIANTCENTRIFUGE_BATCH_SIZE={env_val!r} is not a valid integer; "
                "falling back to auto batch size"
            )
    return min(4 * n_workers, 500)


def _process_gene_group(
    gene: str,
    gene_df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
    gt_matrix: np.ndarray | None = None,
    sample_to_idx: dict[str, int] | None = None,
    gene_row_indices: np.ndarray | None = None,
    variant_keys: np.ndarray | None = None,
) -> tuple[str, dict[str, Any]]:
    """
    Process a single gene's variants for compound heterozygous patterns.

    This function is designed to be run in a worker thread.

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
    gt_matrix : np.ndarray, optional
        Pre-built genotype matrix from Pass 1 (read-only, thread-safe).
    sample_to_idx : Dict[str, int], optional
        Sample-to-column mapping for gt_matrix.
    gene_row_indices : np.ndarray, optional
        Row indices into gt_matrix for this gene.
    variant_keys : np.ndarray, optional
        Pre-computed variant key strings for this gene's rows.

    Returns
    -------
    Tuple[str, Dict[str, Any]]
        Gene name and compound het results
    """
    if pd.isna(gene) or gene == "":
        return gene, {}

    if len(gene_df) <= 1:
        return gene, {}

    # Use vectorized implementation with optional pre-built matrix
    comp_het_results = analyze_gene_for_compound_het_vectorized(
        gene_df,
        pedigree_data,
        sample_list,
        gt_matrix=gt_matrix,
        sample_to_idx=sample_to_idx,
        gene_row_indices=gene_row_indices,
        variant_keys=variant_keys,
    )

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

    # Pre-populate pedigree_data for unknown samples to prevent thread race
    # condition in comp_het_vectorized.py (which mutates the dict for missing samples)
    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            pedigree_data[sample_id] = {
                "sample_id": sample_id,
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "2",
            }

    # Create a copy to avoid fragmentation
    df = df.copy()
    df["_inheritance_patterns"] = None
    df["_comp_het_info"] = None
    df["Inheritance_Pattern"] = "none"
    df["Inheritance_Details"] = "{}"

    # Track timing for each pass
    pass_times = {}

    # Pass 1: Per-Variant Pattern Deduction (vectorized)
    # Request genotype matrix back so Pass 2 can reuse it (avoids redundant encoding)
    pass1_start = time.time()
    logger.info("Pass 1: Deducing inheritance patterns per variant (vectorized)")
    pass1_result = vectorized_deduce_patterns(
        df, pedigree_data, sample_list, return_genotype_matrix=True
    )
    if isinstance(pass1_result, tuple):
        patterns_list_all, gt_matrix, sample_to_idx = pass1_result
    else:
        patterns_list_all = pass1_result
        gt_matrix, sample_to_idx = None, None
    df["_inheritance_patterns"] = patterns_list_all
    pass_times["pattern_deduction"] = time.time() - pass1_start

    # Build positional row index array for slicing the genotype matrix per-gene
    row_positions = np.arange(len(df))

    # Pre-compute variant keys for all rows (reused per-gene and for result application)
    all_variant_keys = build_variant_keys_array(df)

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
        # Group by gene — pre-dispatch deduplication in the main thread so
        # workers receive only NumPy arrays (no DataFrame references).
        gene_groups = df.groupby("GENE", observed=True)
        # Each entry: (gene_name_str, n_unique, gene_row_idx_deduped, gene_vkeys_deduped)
        genes_to_dispatch: list[tuple[str, int, np.ndarray, np.ndarray]] = []
        for gene, group_df in gene_groups:
            if pd.isna(gene) or gene == "":
                continue
            gene_iloc = df.index.get_indexer(group_df.index)
            # Dedup in main thread — reuse all_variant_keys for positional mask
            dedup_mask = ~group_df.duplicated(subset=["CHROM", "POS", "REF", "ALT"])
            n_unique = int(dedup_mask.sum())
            if n_unique < 2:
                continue
            dedup_iloc = gene_iloc[dedup_mask.values]
            gene_row_idx_deduped = row_positions[dedup_iloc]
            gene_vkeys_deduped = all_variant_keys[dedup_iloc]
            genes_to_dispatch.append(
                (str(gene), n_unique, gene_row_idx_deduped, gene_vkeys_deduped)
            )

        num_genes = len(genes_to_dispatch)
        if num_genes > 0:
            # Sort by n_unique descending (largest genes first for load balancing)
            genes_to_dispatch.sort(key=lambda x: x[1], reverse=True)
            max_gene_size = genes_to_dispatch[0][1]

            # Auto-detect worker count if not specified
            if n_workers is None:
                memory_per_gene_gb = rm.estimate_memory(max_gene_size, len(sample_list))
                n_workers = rm.auto_workers(
                    task_count=num_genes, memory_per_task_gb=memory_per_gene_gb
                )

            # Build pedigree arrays once in main thread — workers use array lookups
            # instead of dict.get() calls (O(1), no lock contention).
            # sample_to_idx is guaranteed non-None when gt_matrix is returned from Pass 1.
            eff_sample_to_idx: dict[str, int] = sample_to_idx if sample_to_idx is not None else {}
            father_idx_arr, mother_idx_arr, affected_arr = build_pedigree_arrays(
                sample_list, pedigree_data, eff_sample_to_idx
            )

            logger.info(
                f"Pass 2: Analyzing {num_genes} genes for compound heterozygous patterns "
                f"in parallel ({n_workers} workers, largest gene: {max_gene_size} variants)"
            )

            # Process genes in parallel using threads (NumPy releases GIL)
            batch_size = _get_batch_size(n_workers)
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                completed = 0
                # Submit in batches to bound concurrent work
                for batch_start in range(0, num_genes, batch_size):
                    batch = genes_to_dispatch[batch_start : batch_start + batch_size]
                    future_to_gene = {
                        executor.submit(
                            _process_gene_group_arrays,
                            gene_name,
                            n_unique,
                            gene_row_idx,
                            gene_vkeys,
                            gt_matrix,
                            eff_sample_to_idx,
                            father_idx_arr,
                            mother_idx_arr,
                            affected_arr,
                            sample_list,
                        ): gene_name
                        for gene_name, n_unique, gene_row_idx, gene_vkeys in batch
                    }

                    # Collect results as they complete
                    for future in as_completed(future_to_gene):
                        gene_name_key = future_to_gene[future]
                        try:
                            result_gene_name, comp_het_results = future.result()
                            if comp_het_results:
                                comp_het_results_by_gene[result_gene_name] = comp_het_results
                                logger.debug(
                                    f"Found compound het patterns in gene {result_gene_name}: "
                                    f"{len(comp_het_results)} variants"
                                )
                            completed += 1
                            if completed % 100 == 0:
                                logger.info(f"Processed {completed}/{num_genes} genes")
                        except Exception as e:
                            logger.error(f"Error processing gene {gene_name_key!r}: {e}")
    else:
        # Fall back to sequential processing
        logger.info("Pass 2: Analyzing compound heterozygous patterns sequentially")

        if "GENE" in df.columns:
            for gene, gene_df in df.groupby("GENE", observed=True):
                if pd.isna(gene) or gene == "" or len(gene_df) <= 1:
                    continue

                try:
                    gene_iloc = df.index.get_indexer(gene_df.index)
                    gene_row_idx = row_positions[gene_iloc]
                    gene_vkeys = all_variant_keys[gene_iloc]

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
                        comp_het_results_by_gene[str(gene)] = comp_het_results
                except Exception as e:
                    logger.error(f"Error processing gene {gene!r}: {e}")

    pass_times["compound_het_analysis"] = time.time() - pass2_start

    # Apply compound het results to DataFrame — batch assignment
    apply_start = time.time()
    if comp_het_results_by_gene:
        # Reuse pre-computed variant keys
        variant_keys_series = pd.Series(all_variant_keys, index=df.index)

        # Get genes array - handle categorical dtype
        if isinstance(df["GENE"].dtype, pd.CategoricalDtype):
            genes = df["GENE"].astype(str).values
        else:
            genes = df["GENE"].values.astype(str)

        # Collect all (index, info) pairs, then assign
        batch_indices: list[int] = []
        batch_values: list[dict] = []

        for gene, gene_results in comp_het_results_by_gene.items():
            gene_mask = genes == gene
            if not np.any(gene_mask):
                continue

            gene_variant_keys = variant_keys_series[gene_mask]

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

    pass_times["apply_comp_het_results"] = time.time() - apply_start

    # Pass 3: Prioritize and Finalize
    pass3_start = time.time()
    logger.info("Pass 3: Prioritizing patterns and creating final output")

    _finalize_inheritance_patterns(df, pedigree_data, sample_list)

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
