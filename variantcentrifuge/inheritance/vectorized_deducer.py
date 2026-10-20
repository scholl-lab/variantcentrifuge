"""
Vectorized inheritance pattern deduction using NumPy boolean masks.

This module provides a vectorized implementation of Pass 1 (pattern deduction)
for inheritance analysis. Instead of using df.apply() which calls Python
per-row, it operates on entire columns using NumPy boolean masks for 10-100x
speedup.

Performance benefits:
- NumPy vectorized operations eliminate Python interpreter overhead
- Boolean masks operate on entire arrays at once
- Memory-efficient int8 encoding for genotypes
"""

import logging
import time
from typing import Any

import numpy as np
import pandas as pd

from ..ped_reader import get_parents, is_affected
from .comp_het_vectorized import encode_genotypes
from .non_autosomal_deducer import (
    _check_mitochondrial_vectorized,
    _check_x_fallback_vectorized,
    _check_x_linked_vectorized,
    _check_xld_vectorized,
    _check_xlr_female_vectorized,
    _check_xlr_male_vectorized,
    autosomal_mask,
)

__all__ = [
    "_check_mitochondrial_vectorized",
    "_check_x_fallback_vectorized",
    "_check_x_linked_vectorized",
    "_check_xld_vectorized",
    "_check_xlr_female_vectorized",
    "_check_xlr_male_vectorized",
    "vectorized_deduce_patterns",
]

logger = logging.getLogger(__name__)


def vectorized_deduce_patterns(
    df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
    return_genotype_matrix: bool = False,
) -> list[list[str]] | tuple[list[list[str]], np.ndarray, dict[str, int]]:
    """
    Vectorized inheritance pattern deduction for ALL variants at once.

    This function replaces the df.apply(deduce_patterns_for_variant) approach
    with NumPy boolean mask operations for significant performance improvement.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing variant data with sample genotypes
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information dictionary
    sample_list : List[str]
        List of sample IDs to analyze
    return_genotype_matrix : bool, optional
        If True, also return the genotype matrix and sample-to-index mapping
        built during Pass 1 so they can be reused by Pass 2 (compound het).

    Returns
    -------
    List[List[str]] or Tuple[List[List[str]], np.ndarray, Dict[str, int]]
        List of pattern lists (one per variant), matching the output
        format of the original deduce_patterns_for_variant.
        When return_genotype_matrix is True, returns a tuple of
        (patterns, gt_matrix, sample_to_idx).
    """
    t_start = time.monotonic()
    n_variants = len(df)

    # Initialize result structure - list of lists matching original format
    patterns_per_variant: list[list[str]] = [[] for _ in range(n_variants)]

    # Handle single sample / no pedigree case
    if not pedigree_data or len(sample_list) == 1:
        patterns = _deduce_single_sample_patterns_vectorized(df, sample_list, n_variants)
        if return_genotype_matrix:
            gt_matrix, sample_to_idx = _build_genotype_matrix(df, sample_list)
            return patterns, gt_matrix, sample_to_idx
        return patterns

    # Step 1: Build genotype matrix (n_variants x n_samples)
    # This encodes ALL sample genotypes ONCE before pattern checking
    t_matrix_start = time.monotonic()
    gt_matrix, sample_to_idx = _build_genotype_matrix(df, sample_list)
    t_matrix = time.monotonic() - t_matrix_start
    logger.info(
        f"Vectorized deducer: genotype matrix built in {t_matrix:.3f}s "
        f"({n_variants} variants x {len(sample_list)} samples, "
        f"matrix shape={gt_matrix.shape}, dtype={gt_matrix.dtype})"
    )

    # Compute autosome mask
    chrom_clean_series = df["CHROM"].astype(str).str.strip().str.upper()
    chrom_array = chrom_clean_series.to_numpy()
    pos = df["POS"] if "POS" in df else [np.nan] * n_variants
    autosome_mask = autosomal_mask(df["CHROM"], pos)
    # PAR rows are autosomal: hide them from the X-linked checks
    x_chrom_array = np.where(autosome_mask, "", chrom_array)

    # Step 2: Per-sample vectorized pattern deduction
    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            continue

        # Get sample index in matrix
        if sample_id not in sample_to_idx:
            continue

        sample_idx = sample_to_idx[sample_id]
        sample_gts = gt_matrix[:, sample_idx]

        # Create has_variant mask (gt > 0 means het or hom_alt)
        has_variant_mask = sample_gts > 0

        # Get parent information
        father_id, mother_id = get_parents(sample_id, pedigree_data)
        affected = is_affected(sample_id, pedigree_data)

        # Get parent genotypes if available
        father_gts = None
        mother_gts = None
        if father_id and father_id in sample_to_idx:
            father_gts = gt_matrix[:, sample_to_idx[father_id]]
        if mother_id and mother_id in sample_to_idx:
            mother_gts = gt_matrix[:, sample_to_idx[mother_id]]

        # Track which patterns were added by this sample (for fallback logic)
        patterns_before_sample = [len(p) for p in patterns_per_variant]

        # a. De novo check (vectorized)
        _check_de_novo_vectorized(
            patterns_per_variant,
            sample_gts,
            father_gts,
            mother_gts,
            father_id,
            mother_id,
            sample_list,
            has_variant_mask,
            affected,
        )

        # b. Dominant check (vectorized)
        _check_dominant_vectorized(
            patterns_per_variant,
            sample_gts,
            father_gts,
            mother_gts,
            father_id,
            mother_id,
            pedigree_data,
            has_variant_mask,
            affected,
            autosome_mask=autosome_mask,
        )

        # c. Recessive check (vectorized)
        _check_recessive_vectorized(
            patterns_per_variant,
            sample_gts,
            father_gts,
            mother_gts,
            father_id,
            mother_id,
            affected,
            autosome_mask=autosome_mask,
        )

        # d. X-linked check (vectorized)
        _check_x_linked_vectorized(
            patterns_per_variant,
            sample_id,
            sample_gts,
            father_gts,
            mother_gts,
            father_id,
            mother_id,
            pedigree_data,
            sample_list,
            x_chrom_array,
            has_variant_mask,
            affected,
        )

        # e. Mitochondrial check (vectorized)
        _check_mitochondrial_vectorized(
            patterns_per_variant,
            sample_gts,
            mother_gts,
            mother_id,
            sample_list,
            chrom_array,
            has_variant_mask,
        )

        # f. Fallback for variants with no patterns FROM THIS SAMPLE
        # This matches original logic: if this sample added no patterns but has variant,
        # add "unknown" (affected) or "carrier" (unaffected)
        _apply_fallback_vectorized(
            patterns_per_variant, sample_gts, has_variant_mask, affected, patterns_before_sample
        )

    t_patterns = time.monotonic() - t_start - t_matrix
    logger.info(
        f"Vectorized deducer: pattern checks completed in {t_patterns:.3f}s "
        f"({len(sample_list)} samples x {n_variants} variants)"
    )

    # Step 3: Deduplication per variant (preserve order)
    for i in range(n_variants):
        if patterns_per_variant[i]:
            # Remove duplicates while preserving order
            seen = set()
            unique_patterns = []
            for pattern in patterns_per_variant[i]:
                if pattern not in seen:
                    seen.add(pattern)
                    unique_patterns.append(pattern)
            patterns_per_variant[i] = unique_patterns

    t_total = time.monotonic() - t_start
    n_with_patterns = sum(1 for p in patterns_per_variant if p)
    logger.info(
        f"Vectorized deducer complete in {t_total:.3f}s: "
        f"{n_with_patterns}/{n_variants} variants have patterns "
        f"(matrix={t_matrix:.3f}s, patterns={t_patterns:.3f}s)"
    )

    if return_genotype_matrix:
        return patterns_per_variant, gt_matrix, sample_to_idx
    return patterns_per_variant


def _build_genotype_matrix(
    df: pd.DataFrame, sample_list: list[str]
) -> tuple[np.ndarray, dict[str, int]]:
    """
    Build a genotype matrix (n_variants x n_samples) with int8 encoding.

    Parameters
    ----------
    df : pd.DataFrame
        Variant DataFrame
    sample_list : List[str]
        List of sample IDs

    Returns
    -------
    Tuple[np.ndarray, Dict[str, int]]
        - Genotype matrix (n_variants x n_samples)
        - Mapping from sample_id to column index
    """
    n_variants = len(df)
    n_samples = len(sample_list)

    # Initialize matrix with -1 (missing)
    gt_matrix = np.full((n_variants, n_samples), -1, dtype=np.int8)

    # Build sample to index mapping
    sample_to_idx = {}

    for i, sample_id in enumerate(sample_list):
        sample_to_idx[sample_id] = i

        # Encode genotypes for this sample if column exists
        if sample_id in df.columns:
            gt_matrix[:, i] = encode_genotypes(df[sample_id])

    return gt_matrix, sample_to_idx


def _deduce_single_sample_patterns_vectorized(
    df: pd.DataFrame, sample_list: list[str], n_variants: int
) -> list[list[str]]:
    """
    Vectorized single sample pattern deduction (no pedigree data).

    Parameters
    ----------
    df : pd.DataFrame
        Variant DataFrame
    sample_list : List[str]
        List of sample IDs
    n_variants : int
        Number of variants

    Returns
    -------
    List[List[str]]
        Pattern list per variant
    """
    patterns_per_variant: list[list[str]] = [[] for _ in range(n_variants)]

    for sample_id in sample_list:
        if sample_id not in df.columns:
            continue

        # Encode genotypes
        sample_gts = encode_genotypes(df[sample_id])

        # Create masks for each genotype type
        ref_mask = sample_gts == 0
        hom_alt_mask = sample_gts == 2
        het_mask = sample_gts == 1

        # Apply patterns based on genotype
        # Reference: append "reference"
        ref_indices = np.where(ref_mask)[0]
        for idx in ref_indices:
            patterns_per_variant[idx].append("reference")

        # Homozygous alt: append "homozygous"
        hom_alt_indices = np.where(hom_alt_mask)[0]
        for idx in hom_alt_indices:
            patterns_per_variant[idx].append("homozygous")

        # Heterozygous: append "unknown"
        het_indices = np.where(het_mask)[0]
        for idx in het_indices:
            patterns_per_variant[idx].append("unknown")

    # For variants with no pattern from any sample, append "unknown"
    for i in range(n_variants):
        if not patterns_per_variant[i]:
            patterns_per_variant[i].append("unknown")

    return patterns_per_variant


def _check_de_novo_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    father_gts: np.ndarray | None,
    mother_gts: np.ndarray | None,
    father_id: str | None,
    mother_id: str | None,
    sample_list: list[str],
    has_variant_mask: np.ndarray,
    affected: bool,
) -> None:
    """
    Vectorized de novo pattern detection.

    De novo: child has variant, both parents are ref (0/0).
    De novo candidate: child has variant, at least one parent has missing GT.
    """
    # Both parents must be in sample list
    if not (father_id and mother_id and father_id in sample_list and mother_id in sample_list):
        return

    if father_gts is None or mother_gts is None:
        return

    # Classic de novo: child has variant, both parents ref
    de_novo_mask = has_variant_mask & (father_gts == 0) & (mother_gts == 0)
    de_novo_indices = np.where(de_novo_mask)[0]
    for idx in de_novo_indices:
        patterns_per_variant[idx].append("de_novo")

    # De novo candidate: child has variant, at least one parent missing GT, child affected,
    # and neither parent carries a variant
    if affected:
        de_novo_candidate_mask = (
            has_variant_mask
            & ((father_gts == -1) | (mother_gts == -1))
            & (father_gts <= 0)
            & (mother_gts <= 0)
        )
        de_novo_candidate_indices = np.where(de_novo_candidate_mask)[0]
        for idx in de_novo_candidate_indices:
            patterns_per_variant[idx].append("de_novo_candidate")


def _check_dominant_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    father_gts: np.ndarray | None,
    mother_gts: np.ndarray | None,
    father_id: str | None,
    mother_id: str | None,
    pedigree_data: dict[str, dict[str, Any]],
    has_variant_mask: np.ndarray,
    affected: bool,
    autosome_mask: np.ndarray | None = None,
) -> None:
    """
    Vectorized autosomal dominant pattern detection.

    Classic AD: affected child has variant, at least one affected parent has variant.
    AD possible: affected child has variant, parent has variant but not affected,
                 or parent data incomplete.
    """
    # Must be affected
    if not affected:
        return

    if autosome_mask is not None:
        has_variant_mask = has_variant_mask & autosome_mask

    # Must have variant
    if not has_variant_mask.any():
        return

    # Check father contribution
    father_has_variant_mask = np.zeros(len(sample_gts), dtype=bool)
    father_affected: np.bool_ = np.bool_(False)
    father_data_complete = np.ones(len(sample_gts), dtype=bool)

    if father_id and father_gts is not None:
        father_affected = np.bool_(is_affected(father_id, pedigree_data))
        father_has_variant_mask = father_gts > 0  # type: ignore[assignment]
        father_data_complete = father_gts != -1  # type: ignore[assignment]

    # Check mother contribution
    mother_has_variant_mask = np.zeros(len(sample_gts), dtype=bool)
    mother_affected: np.bool_ = np.bool_(False)
    mother_data_complete = np.ones(len(sample_gts), dtype=bool)

    if mother_id and mother_gts is not None:
        mother_affected = np.bool_(is_affected(mother_id, pedigree_data))
        mother_has_variant_mask = mother_gts > 0  # type: ignore[assignment]
        mother_data_complete = mother_gts != -1  # type: ignore[assignment]

    # Classic dominant: at least one parent affected and has variant
    classic_ad_mask = has_variant_mask & (
        (father_has_variant_mask & father_affected) | (mother_has_variant_mask & mother_affected)
    )
    classic_ad_indices = np.where(classic_ad_mask)[0]
    for idx in classic_ad_indices:
        patterns_per_variant[idx].append("autosomal_dominant")

    # Parent has variant but not affected - incomplete penetrance
    parent_variant_not_affected_mask = has_variant_mask & (
        (father_has_variant_mask & ~father_affected) | (mother_has_variant_mask & ~mother_affected)
    )
    # Exclude indices already classified as classic AD
    parent_variant_not_affected_mask[classic_ad_indices] = False
    parent_variant_not_affected_indices = np.where(parent_variant_not_affected_mask)[0]
    for idx in parent_variant_not_affected_indices:
        patterns_per_variant[idx].append("autosomal_dominant_possible")

    # Incomplete parent data
    incomplete_data_mask = has_variant_mask & (~father_data_complete | ~mother_data_complete)
    # Exclude indices already classified
    incomplete_data_mask[classic_ad_indices] = False
    incomplete_data_mask[parent_variant_not_affected_indices] = False
    incomplete_data_indices = np.where(incomplete_data_mask)[0]
    for idx in incomplete_data_indices:
        patterns_per_variant[idx].append("autosomal_dominant_possible")


def _check_recessive_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    father_gts: np.ndarray | None,
    mother_gts: np.ndarray | None,
    father_id: str | None,
    mother_id: str | None,
    affected: bool,
    autosome_mask: np.ndarray | None = None,
) -> None:
    """
    Vectorized autosomal recessive pattern detection.

    Classic AR: affected child is hom_alt (1/1), both parents are carriers (have variant).
    AR possible: affected child is hom_alt, incomplete parent data consistent with AR.
    """
    # Must be affected
    if not affected:
        return

    # Must be homozygous alt
    hom_alt_mask = sample_gts == 2
    if autosome_mask is not None:
        hom_alt_mask = hom_alt_mask & autosome_mask

    if not hom_alt_mask.any():
        return

    # Check parent carrier status
    father_carrier_mask = np.zeros(len(sample_gts), dtype=bool)
    father_data_complete = np.ones(len(sample_gts), dtype=bool)

    if father_id and father_gts is not None:
        father_carrier_mask = father_gts > 0  # type: ignore[assignment]
        father_data_complete = father_gts != -1  # type: ignore[assignment]

    mother_carrier_mask = np.zeros(len(sample_gts), dtype=bool)
    mother_data_complete = np.ones(len(sample_gts), dtype=bool)

    if mother_id and mother_gts is not None:
        mother_carrier_mask = mother_gts > 0  # type: ignore[assignment]
        mother_data_complete = mother_gts != -1  # type: ignore[assignment]

    # Classic recessive: both parents are carriers
    classic_ar_mask = hom_alt_mask & father_carrier_mask & mother_carrier_mask
    classic_ar_indices = np.where(classic_ar_mask)[0]
    for idx in classic_ar_indices:
        patterns_per_variant[idx].append("autosomal_recessive")

    # Incomplete parent data consistent with AR
    # (parent data incomplete but what we have is consistent)
    incomplete_ar_mask = (
        hom_alt_mask
        & ~classic_ar_mask
        & (~father_data_complete | ~mother_data_complete)
        & (father_carrier_mask | ~father_data_complete)
        & (mother_carrier_mask | ~mother_data_complete)
    )
    incomplete_ar_indices = np.where(incomplete_ar_mask)[0]
    for idx in incomplete_ar_indices:
        patterns_per_variant[idx].append("autosomal_recessive_possible")


def _apply_fallback_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    has_variant_mask: np.ndarray,
    affected: bool,
    patterns_before_sample: list[int],
) -> None:
    """
    Apply fallback patterns for variants where THIS SAMPLE has variant but added no patterns.

    This matches original per-sample logic: if a sample has a variant but contributed
    no specific patterns, add "unknown" (affected) or "carrier" (unaffected).

    Parameters
    ----------
    patterns_per_variant : List[List[str]]
        Pattern lists per variant
    sample_gts : np.ndarray
        Sample genotypes
    has_variant_mask : np.ndarray
        Boolean mask where sample has variant
    affected : bool
        Whether sample is affected
    patterns_before_sample : List[int]
        Number of patterns per variant before this sample's analysis
    """
    # Find variants where sample has variant
    variant_indices = np.where(has_variant_mask)[0]

    for idx in variant_indices:
        # Check if THIS SAMPLE added any patterns
        # (compare current count to count before this sample)
        patterns_added_by_sample = len(patterns_per_variant[idx]) - patterns_before_sample[idx]

        if patterns_added_by_sample == 0:
            # This sample has variant but added no patterns - apply fallback
            if affected:
                patterns_per_variant[idx].append("unknown")
            else:
                patterns_per_variant[idx].append("carrier")
