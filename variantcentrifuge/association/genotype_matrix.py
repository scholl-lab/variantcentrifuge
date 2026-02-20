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
        Additive dosage: 0 (hom-ref), 1 (het), 2 (hom-alt), or ``None``
        (missing — any allele is ``.``).
    is_multi_allelic : bool
        ``True`` when any allele integer > 1 (e.g. ``1/2``, ``0/2``,
        ``2/2``). Multi-allelic genotypes are assigned het-equivalent
        dosage 1.

    Rules
    -----
    - ``'0/0'``, ``'0|0'``           -> ``(0, False)``
    - ``'0/1'``, ``'1/0'``, phased   -> ``(1, False)``
    - ``'1/1'``, ``'1|1'``           -> ``(2, False)``
    - ``'./.'``, ``'.|.'``, None, '' -> ``(None, False)``
    - ``'./1'``, ``'1/.'``           -> ``(None, False)``  (partial missing)
    - ``'1/2'``, ``'0/2'``, ``'2/2'``-> ``(1, True)``     (multi-allelic)

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
) -> tuple[np.ndarray, np.ndarray, list[bool], list[str]]:
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
        Determines imputation mode.
        - ``True`` (binary trait): impute missing as ``round(2 * MAF)``
          (bestguess, SKAT R convention).
        - ``False`` (quantitative): impute missing as ``2 * MAF`` (mean).
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

    Notes
    -----
    MAF is computed before imputation (from observed values only) to avoid
    circular dependency.  Variants where every sample is missing receive
    MAF=0 and are imputed to 0 (ref/ref).
    """
    vcf_samples_list = list(vcf_samples)
    gt_columns_list = list(gt_columns)
    n_samples = len(vcf_samples_list)
    n_variants = len(gene_df)
    warnings_list: list[str] = []

    # ------------------------------------------------------------------ #
    # Step 1: Parse raw dosages -> (n_variants, n_samples), NaN=missing   #
    # ------------------------------------------------------------------ #
    raw = np.full((n_variants, n_samples), np.nan, dtype=np.float64)
    multi_allelic_count = 0

    for v_idx, (_, row) in enumerate(gene_df.iterrows()):
        for s_idx, col in enumerate(gt_columns_list):
            gt_val = row.get(col, "./.")
            dosage, is_multi = parse_gt_to_dosage(str(gt_val) if gt_val is not None else "")
            if dosage is not None:
                raw[v_idx, s_idx] = float(dosage)
            if is_multi:
                multi_allelic_count += 1

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
        observed = geno[v_idx, ~np.isnan(geno[v_idx])]
        if len(observed) == 0:
            # Every sample missing: MAF=0, impute 0
            mafs[v_idx] = 0.0
            geno[v_idx] = 0.0
            continue
        maf = float(observed.mean()) / 2.0
        mafs[v_idx] = maf
        missing_mask = np.isnan(geno[v_idx])
        if missing_mask.any():
            imputed_val = round(2.0 * maf) if is_binary else 2.0 * maf
            geno[v_idx, missing_mask] = imputed_val

    # ------------------------------------------------------------------ #
    # Step 6: Transpose to (n_samples, n_variants)                        #
    # ------------------------------------------------------------------ #
    geno = geno.T  # shape (n_samples, n_kept_variants)

    return geno, mafs, sample_mask, warnings_list
