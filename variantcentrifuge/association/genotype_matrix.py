# File: variantcentrifuge/association/genotype_matrix.py
# Location: variantcentrifuge/variantcentrifuge/association/genotype_matrix.py
"""
Genotype matrix construction for rare variant association burden tests.

Provides two public functions:

- ``parse_gt_to_dosage``: Parse a VCF GT string to an additive dosage (0/1/2)
  with a multi-allelic flag. Returns ``(dosage, is_multi_allelic)``.
- ``build_genotype_matrix``: Build an imputed ``(n_samples, n_variants)``
  float64 matrix from a per-gene DataFrame, applying a two-layer
  missing-data strategy that matches the SKAT R package convention.

Design notes
------------
- ``parse_gt_to_dosage`` is a NEW implementation. It must NOT delegate to
  ``gene_burden._gt_to_dosage()`` which returns 0 for '1/2' (wrong for
  regression — het-equivalent dosage 1 is correct).
- MAF for imputation is computed from ALL samples combined (never stratified
  by phenotype) to match SKAT R convention and avoid label leakage.
- Genotype matrices are never stored across genes; callers should build per
  gene and allow the matrix to be garbage-collected after the test result
  is captured.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger("variantcentrifuge")


def _vectorized_parse_gt_columns(
    gene_df: pd.DataFrame,
    gt_columns_list: list[str],
) -> tuple[np.ndarray, np.ndarray, int]:
    """Parse all GT columns to a dosage matrix and ploidy matrix using vectorized pandas operations.

    Replaces the nested iterrows() loop with column-wise vectorized string matching.
    Produces identical results to per-cell parse_gt_to_dosage() calls.

    Parameters
    ----------
    gene_df : pd.DataFrame
        Per-gene DataFrame with GT string columns.
    gt_columns_list : list[str]
        Ordered GT column names (one per sample).

    Returns
    -------
    raw : np.ndarray, shape (n_variants, n_samples), float64
        Dosage matrix. NaN for missing genotypes.
    ploidy : np.ndarray, shape (n_variants, n_samples), float64
        Ploidy matrix (1.0 for haploid, 2.0 for diploid).
    multi_allelic_count : int
        Total number of multi-allelic genotype entries detected.
    """
    n_variants = len(gene_df)
    n_samples = len(gt_columns_list)
    raw = np.full((n_variants, n_samples), np.nan, dtype=np.float64)
    ploidy = np.full((n_variants, n_samples), 2.0, dtype=np.float64)
    multi_allelic_count = 0
    row_has_haploid = np.zeros(n_variants, dtype=bool)
    row_has_diploid = np.zeros(n_variants, dtype=bool)
    generic_missing = np.zeros((n_variants, n_samples), dtype=bool)

    for s_idx, col in enumerate(gt_columns_list):
        s_obj = (
            gene_df[col].astype(object)
            if isinstance(gene_df[col].dtype, pd.CategoricalDtype)
            else gene_df[col]
        )
        s_str = s_obj.fillna("").astype(str).str.strip()
        is_gen_miss = s_obj.isna() | s_str.str.lower().isin({"", "na", "nan", "none", "null"})
        generic_missing[:, s_idx] = is_gen_miss.to_numpy()

        gt_series = s_str.str.replace("|", "/", regex=False)
        m_diploid = ~is_gen_miss & gt_series.str.contains("/", regex=False)
        row_has_diploid |= m_diploid.to_numpy()

        # Standard diploid genotypes — exact match for speed
        m_00 = gt_series == "0/0"
        m_01 = (gt_series == "0/1") | (gt_series == "1/0")
        m_11 = gt_series == "1/1"

        # Haploid genotypes
        m_hap_0 = gt_series == "0"
        m_hap_1 = gt_series == "1"
        m_hap_miss = gt_series == "."
        row_has_haploid |= (m_hap_0 | m_hap_1 | m_hap_miss).to_numpy()

        # Multi-allelic: any allele > 1 (e.g. 1/2, 0/2, 2/2, 0/3)
        multi_mask = gt_series.str.match(r"^[0-9]+/[0-9]+$", na=False) & ~(m_00 | m_01 | m_11)
        hap_multi_mask = (
            ~is_gen_miss
            & ~m_diploid
            & gt_series.str.match(r"^[0-9]+$", na=False)
            & ~(m_hap_0 | m_hap_1)
        )

        # Filter: only flag as multi if actually contains allele > 1
        if multi_mask.any():
            candidates = gt_series[multi_mask]
            left = candidates.str.split("/").str[0]
            right = candidates.str.split("/").str[1]
            truly_multi = (left.astype(int) > 1) | (right.astype(int) > 1)
            multi_idx = multi_mask[multi_mask].index
            truly_multi_values = truly_multi.to_numpy(dtype=bool)
            not_truly = multi_idx[~truly_multi_values]
            multi_mask.loc[not_truly] = False
            multi_allelic_count += int(truly_multi.sum())

        if hap_multi_mask.any():
            candidates_hap = gt_series[hap_multi_mask]
            truly_hap_multi = candidates_hap.astype(int) > 1
            hap_multi_idx = hap_multi_mask[hap_multi_mask].index
            truly_hap_multi_values = truly_hap_multi.to_numpy(dtype=bool)
            not_truly_hap = hap_multi_idx[~truly_hap_multi_values]
            hap_multi_mask.loc[not_truly_hap] = False
            multi_allelic_count += int(truly_hap_multi.sum())
            row_has_haploid |= hap_multi_mask.to_numpy()

        # Assign dosages and ploidy
        col_dosage = np.full(n_variants, np.nan, dtype=np.float64)
        col_ploidy = np.full(n_variants, 2.0, dtype=np.float64)

        col_dosage[m_00.values] = 0.0
        col_dosage[m_01.values] = 1.0
        col_dosage[m_11.values] = 2.0
        col_dosage[multi_mask.values] = 1.0  # het-equivalent for multi-allelic

        col_dosage[m_hap_0.values] = 0.0
        col_dosage[m_hap_1.values] = 1.0
        col_dosage[hap_multi_mask.values] = 1.0

        col_ploidy[m_hap_0.values] = 1.0
        col_ploidy[m_hap_1.values] = 1.0
        col_ploidy[hap_multi_mask.values] = 1.0
        col_ploidy[m_hap_miss.values] = 1.0

        raw[:, s_idx] = col_dosage
        ploidy[:, s_idx] = col_ploidy

    haploid_rows = row_has_haploid & ~row_has_diploid
    if haploid_rows.any():
        ploidy[haploid_rows[:, None] & generic_missing] = 1.0

    return raw, ploidy, multi_allelic_count


def parse_gt_to_dosage(gt: str | None) -> tuple[int | None, bool]:
    """
    Parse a VCF GT string to additive dosage and a multi-allelic flag.

    Parameters
    ----------
    gt : str | None
        VCF genotype string. May be phased (``|``) or unphased (``/``).

    Returns
    -------
    dosage : int | None
        Additive dosage: 0 (hom-ref/hap-ref), 1 (het/hap-alt), 2 (hom-alt),
        or ``None`` (missing — any allele is ``.``).
    is_multi_allelic : bool
        ``True`` when any allele integer > 1 (e.g. ``1/2``, ``0/2``,
        ``2/2``, ``2``). Multi-allelic genotypes are assigned het-equivalent
        dosage 1.

    Rules
    -----
    - ``'0/0'``, ``'0|0'``, ``'0'``  -> ``(0, False)``
    - ``'0/1'``, ``'1/0'``, phased, ``'1'`` -> ``(1, False)``
    - ``'1/1'``, ``'1|1'``           -> ``(2, False)``
    - ``'./.'``, ``'.|.'``, ``'.'``, None, '' -> ``(None, False)``
    - ``'./1'``, ``'1/.'``           -> ``(None, False)``  (partial missing)
    - ``'1/2'``, ``'0/2'``, ``'2/2'``, ``'2'`` -> ``(1, True)``  (multi-allelic)

    Notes
    -----
    Phased genotypes (``|`` separator) are treated identically to unphased
    (``/``): the separator is normalised before parsing.
    """
    if not gt:
        return None, False

    # Normalise separator so we only need to handle '/'
    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)

    # Haploid call (single allele)
    if len(parts) == 1:
        if parts[0] in (".", ""):
            return None, False
        try:
            val = int(parts[0])
            if val == 0:
                return 0, False
            elif val == 1:
                return 1, False
            else:
                return 1, True  # multi-allelic
        except ValueError:
            return None, False

    if len(parts) != 2:
        return None, False

    # Parse individual allele strings
    try:
        a1 = None if parts[0] == "." else int(parts[0])
        a2 = None if parts[1] == "." else int(parts[1])
    except ValueError:
        return None, False

    # Partial or full missing
    if a1 is None or a2 is None:
        return None, False

    # Multi-allelic: any allele index > 1
    if a1 > 1 or a2 > 1:
        return 1, True

    # Standard diploid dosage
    return a1 + a2, False


def build_genotype_matrix(
    gene_df: pd.DataFrame,
    vcf_samples: Sequence[str],
    gt_columns: Sequence[str],
    is_binary: bool = True,
    missing_site_threshold: float = 0.10,
    missing_sample_threshold: float = 0.80,
    phenotype_vector: np.ndarray | None = None,
    return_keep_mask: bool = False,
) -> (
    tuple[np.ndarray, np.ndarray, list[bool], list[str]]
    | tuple[np.ndarray, np.ndarray, list[bool], list[str], np.ndarray]
):
    """
    Build a fully imputed genotype matrix from a per-gene variant DataFrame.

    Parameters
    ----------
    gene_df : pd.DataFrame
        Rows = variants for one gene. Must contain columns listed in
        ``gt_columns``.
    vcf_samples : sequence of str
        Ordered sample identifiers matching the column order in ``gt_columns``.
    gt_columns : sequence of str
        Column names that hold per-sample GT strings (e.g. ``GEN_0__GT``).
        The index ``i`` of a column corresponds to ``vcf_samples[i]``.
    is_binary : bool
        Determines imputation mode using ploidy-scaled mean imputation.
        - ``True`` (binary trait): impute missing as ``round(ploidy * alt_allele_freq)``
          (bestguess, SKAT R convention; dosage units {0, 1, 2} for diploid, {0, 1} for haploid).
        - ``False`` (quantitative): impute missing as ``ploidy * alt_allele_freq``
          (mean; ``2 * alt_freq`` for diploid, ``1 * alt_freq`` for haploid).
    missing_site_threshold : float
        Variants with a fraction of missing genotypes > this threshold
        across all samples are removed before imputation. Default: 0.10.
    missing_sample_threshold : float
        Samples with a fraction of missing genotypes > this threshold
        across the *kept* variants are flagged in ``sample_mask``.
        Default: 0.80.
    phenotype_vector : np.ndarray | None, shape (n_samples,)
        Optional binary phenotype (0/1). When provided, differential
        missingness warnings are computed per variant. Ignored if None.
    return_keep_mask : bool
        When True, append the site-level variant keep mask as a fifth return value.
        Default False preserves the historical four-tuple return.

    Returns
    -------
    G : np.ndarray, shape (n_kept_samples, n_kept_variants), float64
        Fully imputed genotype matrix. ``np.isnan(G).sum() == 0`` is
        guaranteed.  Note: currently all samples are included in ``G``
        (``sample_mask`` reports the problematic ones — callers may choose
        to exclude them).
    mafs : np.ndarray, shape (n_kept_variants,), float64
        Per-variant minor allele frequency computed from **all** samples
        combined (never stratified by phenotype — SKAT R convention).
    sample_mask : list[bool], length len(vcf_samples)
        ``True`` for each sample that survived the ``missing_sample_threshold``
        filter. ``False`` samples have high missingness and should be
        considered for exclusion by the caller.
    warnings_list : list[str]
        Diagnostic messages accumulated during construction (e.g. multi-
        allelic detection, differential missingness).
    keep_variants_mask : np.ndarray, optional
        Boolean mask over input rows, returned only when ``return_keep_mask=True``.

    Notes
    -----
    MAF is computed before imputation (from observed values only) to avoid
    circular dependency.  Variants where every sample is missing receive
    MAF=0 and are imputed to 0 (ref/ref).
    """
    vcf_samples_list = list(vcf_samples)
    gt_columns_list = list(gt_columns)
    n_samples = len(vcf_samples_list)
    warnings_list: list[str] = []

    # ------------------------------------------------------------------ #
    # Step 1: Parse raw dosages -> (n_variants, n_samples), NaN=missing   #
    # Vectorized: processes one GT column at a time using pandas string    #
    # ops instead of per-cell Python function calls.                       #
    # ------------------------------------------------------------------ #
    raw, ploidy, multi_allelic_count = _vectorized_parse_gt_columns(gene_df, gt_columns_list)

    if multi_allelic_count > 0:
        warnings_list.append(
            f"Multi-allelic genotypes detected at {multi_allelic_count} variant(s). "
            "Recommend running bcftools norm -m- to split multi-allelic sites."
        )

    # ------------------------------------------------------------------ #
    # Step 2 — Layer 1a: Remove variants with >missing_site_threshold     #
    # fraction missing site-wide                                          #
    # ------------------------------------------------------------------ #
    missing_per_variant = np.isnan(raw).mean(axis=1)  # shape (n_variants,)
    keep_variants_mask = missing_per_variant <= missing_site_threshold
    raw = raw[keep_variants_mask]  # shape (n_kept_variants, n_samples)
    ploidy = ploidy[keep_variants_mask]

    # ------------------------------------------------------------------ #
    # Step 3 — Layer 1b: Flag samples with >missing_sample_threshold      #
    # fraction missing across kept variants                               #
    # ------------------------------------------------------------------ #
    if raw.shape[0] > 0:
        missing_per_sample = np.isnan(raw).mean(axis=0)  # shape (n_samples,)
    else:
        missing_per_sample = np.zeros(n_samples, dtype=np.float64)

    sample_mask: list[bool] = [
        bool(missing_per_sample[i] <= missing_sample_threshold) for i in range(n_samples)
    ]

    # ------------------------------------------------------------------ #
    # Step 4 (optional): Differential missingness warning per variant     #
    # ------------------------------------------------------------------ #
    if phenotype_vector is not None and raw.shape[0] > 0:
        pheno = np.asarray(phenotype_vector, dtype=np.float64)
        case_mask = pheno == 1
        ctrl_mask = pheno == 0
        n_cases = int(case_mask.sum())
        n_ctrls = int(ctrl_mask.sum())
        if n_cases > 0 and n_ctrls > 0:
            for v_idx in range(raw.shape[0]):
                miss_v = np.isnan(raw[v_idx])
                case_missing_rate = miss_v[case_mask].mean()
                ctrl_missing_rate = miss_v[ctrl_mask].mean()
                if abs(case_missing_rate - ctrl_missing_rate) > 0.05:
                    warnings_list.append(
                        f"Variant index {v_idx}: differential missingness "
                        f"(cases={case_missing_rate:.2f}, "
                        f"controls={ctrl_missing_rate:.2f})."
                    )

    # ------------------------------------------------------------------ #
    # Step 5 — Layer 2: Compute MAFs and impute missing values            #
    # MAF from ALL samples combined (never stratified) — SKAT R / Pitfall5#
    # ------------------------------------------------------------------ #
    geno = raw.copy()
    n_kept_variants = geno.shape[0]
    mafs = np.zeros(n_kept_variants, dtype=np.float64)

    for v_idx in range(n_kept_variants):
        obs_mask = ~np.isnan(geno[v_idx])
        observed_dosages = geno[v_idx, obs_mask]
        observed_ploidies = ploidy[v_idx, obs_mask]
        if len(observed_dosages) == 0:
            # Every sample missing: MAF=0, impute 0
            mafs[v_idx] = 0.0
            geno[v_idx] = 0.0
            continue
        total_ac = float(observed_dosages.sum())
        total_an = float(observed_ploidies.sum())
        if total_an <= 0.0 or total_ac >= total_an:
            # Monomorphic for ALT (or no observed chromosomes): MAF=0, zero out dosages
            mafs[v_idx] = 0.0
            geno[v_idx] = 0.0
            continue
        alt_freq = total_ac / total_an if total_an > 0 else 0.0
        mafs[v_idx] = min(alt_freq, 1.0 - alt_freq)
        missing_mask = np.isnan(geno[v_idx])
        if missing_mask.any():
            for s in np.where(missing_mask)[0]:
                s_ploidy = ploidy[v_idx, s]
                imputed_val = round(s_ploidy * maf) if is_binary else s_ploidy * maf
                geno[v_idx, s] = imputed_val

    # ------------------------------------------------------------------ #
    # Step 6: Transpose to (n_samples, n_variants)                        #
    # ------------------------------------------------------------------ #
    geno = geno.T  # shape (n_samples, n_kept_variants)

    if return_keep_mask:
        return geno, mafs, sample_mask, warnings_list, keep_variants_mask
    return geno, mafs, sample_mask, warnings_list
