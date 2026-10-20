# File: variantcentrifuge/association/genotype_matrix.py
# Location: variantcentrifuge/variantcentrifuge/association/genotype_matrix.py
"""
Genotype matrix construction for rare variant association burden tests.

Provides two public functions:

- ``parse_gt_to_dosage``: Parse a VCF GT string to an additive dosage (0/1/2).
  Returns ``(dosage, is_multi_allelic)``; allele indices >1 raise ``ValueError``.
- ``build_genotype_matrix``: Build an imputed ``(n_samples, n_variants)``
  float64 matrix from a per-gene DataFrame, applying a two-layer
  missing-data strategy that matches the SKAT R package convention.

Design notes
------------
- ``parse_gt_to_dosage`` is a NEW implementation. It must NOT delegate to
  ``gene_burden._gt_to_dosage()`` which returns 0 for '1/2' (wrong for
  regression — het-equivalent dosage 1 is correct).
- MAF for imputation is computed from all *retained* samples combined (never
  stratified by phenotype) to match SKAT R convention and avoid label leakage.
  Samples flagged by the missingness filter never influence frequencies,
  imputation or minor-allele recoding (#174).
- Genotype matrices are never stored across genes; callers should build per
  gene and allow the matrix to be garbage-collected after the test result
  is captured.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger("variantcentrifuge")

_SPLIT_HINT = (
    "allele indices >1 cannot be encoded as an additive ALT dosage. Split multi-allelic "
    "sites before association testing, e.g. `bcftools norm -m- in.vcf.gz -Oz -o split.vcf.gz`."
)


def _multiallelic_error(
    gene_df: pd.DataFrame, row: int, column: str, value: object = None
) -> ValueError:
    """Build the actionable #175 error naming the gene, variant and sample column."""
    record = gene_df.iloc[row]
    value = record[column] if value is None else value
    site = ":".join(str(record[c]) for c in ("CHROM", "POS", "REF", "ALT") if c in gene_df)
    gene = f"gene {record['GENE']}, " if "GENE" in gene_df else ""
    return ValueError(
        f"Unsplit multi-allelic genotype {value!r} ({gene}variant "
        f"{site or f'row {gene_df.index[row]!r}'}, column {column!r}): {_SPLIT_HINT}"
    )


_MISSING_TOKENS = {"", "na", "nan", "none", "null"}
_BLOCK_CELLS = 1_000_000
# Per-cell call kinds returned by the parser (row ploidy inference for missing calls).
_HAPLOID, _DIPLOID, _GENERIC_MISSING = 1, 2, 3
_GT_DOSAGE = {"0/0": 0.0, "0/1": 1.0, "1/0": 1.0, "1/1": 2.0, "0": 0.0, "1": 1.0}


def _classify_gt_vocabulary(vocab: pd.Series) -> tuple[np.ndarray, ...]:
    """Classify each distinct GT string once (dosage, ploidy, flags; see caller)."""
    s_str = vocab.str.strip()
    gen_miss = s_str.str.lower().isin(_MISSING_TOKENS).to_numpy()
    gt = s_str.str.replace("|", "/", regex=False)
    diploid = ~gen_miss & gt.str.contains("/", regex=False).to_numpy()
    dosage = gt.map(_GT_DOSAGE).to_numpy(dtype=np.float64)
    haploid = gt.isin({"0", "1", "."}).to_numpy()
    kind = np.select([haploid, diploid, gen_miss], [_HAPLOID, _DIPLOID, _GENERIC_MISSING], 0)
    # #175: any observed allele index >1 (0/2, 2, ./2, 0/1/2, ...) cannot be encoded as
    # additive ALT dosage; flag it so the caller rejects instead of recoding/missing it.
    multi = np.array([_has_allele_above_one(g.split("/")) for g in gt], dtype=bool)
    return dosage, np.where(haploid, 1.0, 2.0), kind.astype(np.int8), multi


def _has_allele_above_one(tokens: list[str]) -> bool:
    return any(t.isascii() and t.isdigit() and int(t) > 1 for t in tokens)


def _infer_missing_ploidy(ploidy: np.ndarray, kind: np.ndarray, retained: np.ndarray) -> None:
    """Give generic missing calls ploidy 1 on rows whose *retained* calls are all haploid.

    In place. Explicit calls keep their encoded ploidy; only retained samples decide the
    row ploidy, so excluded samples cannot change retained imputation (#174).
    """
    kept = kind[:, retained]
    haploid_rows = (kept == _HAPLOID).any(axis=1) & ~(kept == _DIPLOID).any(axis=1)
    ploidy[haploid_rows[:, None] & (kind == _GENERIC_MISSING)] = 1.0


def _column_strings(column: pd.Series) -> np.ndarray:
    """Normalise one GT column to strings exactly as the per-column parser did."""
    if isinstance(column.dtype, pd.CategoricalDtype):
        column = column.astype(object)
    return np.asarray(column.fillna("").astype(str), dtype=object)


def _vectorized_parse_gt_columns(
    gene_df: pd.DataFrame,
    gt_columns_list: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse all GT columns to dosage, explicit ploidy and per-call kind matrices.

    Every cell is normalised like ``Series.fillna("").astype(str)`` per column
    (categoricals via their object values); per sample-column block of about
    ``_BLOCK_CELLS`` cells (bounded memory), the distinct strings are factorised
    and classified once and the per-string results are mapped back by code
    instead of running ~15 pandas string operations per sample column (#176).

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
        Explicit ploidy (1.0 for haploid calls, else 2.0). Generic missing calls
        (None/NA/empty) stay 2.0 until :func:`_infer_missing_ploidy` resolves them
        from the retained samples.
    kind : np.ndarray, shape (n_variants, n_samples), int8
        Call kind for that inference (haploid / diploid / generic missing / other).

    Raises
    ------
    ValueError
        If any genotype carries an allele index >1 (unsplit multi-allelic site).
    """
    n_variants = len(gene_df)
    n_samples = len(gt_columns_list)
    raw = np.full((n_variants, n_samples), np.nan, dtype=np.float64)
    ploidy = np.full((n_variants, n_samples), 2.0, dtype=np.float64)
    kind = np.zeros((n_variants, n_samples), dtype=np.int8)
    # Sample-column blocks of ~_BLOCK_CELLS cells bound the transient object/code arrays.
    cols_per_block = max(1, _BLOCK_CELLS // max(1, n_variants))
    for c0 in range(0, n_samples, cols_per_block):
        c1 = min(c0 + cols_per_block, n_samples)
        codes, vocab = _factorize_gt_block(gene_df[gt_columns_list[c0:c1]])
        dosage, plo, cell_kind, multi = _classify_gt_vocabulary(pd.Series(vocab, dtype=object))
        if multi.any():
            # Column-major within ordered blocks: first offending sample, then row.
            first = int(np.flatnonzero(multi[codes])[0])
            raise _multiallelic_error(
                gene_df, first % n_variants, gt_columns_list[c0 + first // n_variants]
            )
        shape = (n_variants, c1 - c0)
        raw[:, c0:c1] = dosage[codes].reshape(shape, order="F")
        ploidy[:, c0:c1] = plo[codes].reshape(shape, order="F")
        kind[:, c0:c1] = cell_kind[codes].reshape(shape, order="F")
    return raw, ploidy, kind


def _factorize_gt_block(block: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """Codes (column-major) and distinct normalised GT strings for one column block."""
    try:
        codes, vocab = pd.factorize(block.to_numpy(dtype=object).ravel(order="F"))
        all_strings = all(isinstance(v, str) for v in vocab)
    except TypeError:  # unhashable cells (list, dict, ndarray): stringify per column
        all_strings = False
    if all_strings:
        codes[codes < 0] = len(vocab)  # NA cells normalise to ""
        return codes, [*vocab, ""]
    # Non-string values (numbers, bytes): factorize equates 1, 1.0 and True, and the
    # column-wise fillna may downcast numeric objects, so normalise per column.
    codes, uniques = pd.factorize(
        np.concatenate([_column_strings(block.iloc[:, j]) for j in range(block.shape[1])])
    )
    return codes, list(uniques)


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
        Always ``False``: allele indices >1 raise instead (kept for the
        historical tuple shape).

    Raises
    ------
    ValueError
        For allele indices >1 (``1/2``, ``0/2``, ``2/2``, haploid ``2``): unsplit
        multi-allelic genotypes cannot be represented as additive ALT dosage.

    Rules
    -----
    - ``'0/0'``, ``'0|0'``, ``'0'``  -> ``(0, False)``
    - ``'0/1'``, ``'1/0'``, phased, ``'1'`` -> ``(1, False)``
    - ``'1/1'``, ``'1|1'``           -> ``(2, False)``
    - ``'./.'``, ``'.|.'``, ``'.'``, None, '' -> ``(None, False)``
    - ``'./1'``, ``'1/.'``           -> ``(None, False)``  (partial missing)
    - ``'1/2'``, ``'2'``, ``'./2'``, ``'0/1/2'`` (any allele >1) -> ``ValueError``

    Notes
    -----
    Phased genotypes (``|`` separator) are treated identically to unphased
    (``/``): the separator is normalised before parsing.
    """
    if not gt:
        return None, False
    if _has_allele_above_one(re.split(r"[/|]", gt)):
        raise ValueError(f"Unsplit multi-allelic genotype {gt!r}: {_SPLIT_HINT}")

    # Normalise separator so we only need to handle '/'
    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)

    # Haploid call (single allele)
    if len(parts) == 1:
        if parts[0] in (".", ""):
            return None, False
        try:
            val = int(parts[0])
        except ValueError:
            return None, False
        if val > 1:
            raise ValueError(f"Unsplit multi-allelic genotype {gt!r}: {_SPLIT_HINT}")
        return (val, False) if val >= 0 else (None, False)

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
        raise ValueError(f"Unsplit multi-allelic genotype {gt!r}: {_SPLIT_HINT}")

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
    missing_as_ref: bool = False,
    parsed: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
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
    missing_as_ref : bool
        When True, every missing call (including partial calls such as ``./1``) is
        set to reference dosage 0 *before* site/sample missingness QC, so nothing is
        imputed or dropped for missingness (``bcftools +setGT -t . -n 0``). Meant for
        merged single-sample VCFs where ``./.`` means "not called as variant".
        Default False keeps the SKAT conventions (bestguess / mean imputation).
    parsed : tuple of np.ndarray, optional
        Pre-parsed ``(dosage, ploidy, kind)`` arrays (``GenotypeCodes.gene_arrays``); when
        given, ``gt_columns`` is not read (#219).

    Returns
    -------
    G : np.ndarray, shape (n_kept_samples, n_kept_variants), float64
        Fully imputed genotype matrix. ``np.isnan(G).sum() == 0`` is
        guaranteed.  Note: currently all samples are included in ``G``
        (``sample_mask`` reports the problematic ones — callers may choose
        to exclude them).
    mafs : np.ndarray, shape (n_kept_variants,), float64
        Per-variant minor allele frequency computed from the retained samples
        (``sample_mask``) combined, never stratified by phenotype (SKAT R
        convention). Imputation and minor-allele recoding use the same frequency.
    sample_mask : list[bool], length len(vcf_samples)
        ``True`` for each sample that survived the ``missing_sample_threshold``
        filter. ``False`` samples have high missingness and should be
        considered for exclusion by the caller.
    warnings_list : list[str]
        Diagnostic messages accumulated during construction (e.g. differential
        missingness).
    keep_variants_mask : np.ndarray, optional
        Boolean mask over input rows, returned only when ``return_keep_mask=True``.

    Notes
    -----
    MAF is computed before imputation (from observed values only) to avoid
    circular dependency.  Variants where every retained sample is missing receive
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
    if parsed is None:
        parsed = _vectorized_parse_gt_columns(gene_df, gt_columns_list)
    raw, ploidy, kind = parsed
    if missing_as_ref:
        # Missing ploidy is still inferred below (generic missing calls on haploid rows -> 1).
        raw[np.isnan(raw)] = 0.0

    # ------------------------------------------------------------------ #
    # Step 2 — Layer 1a: Remove variants with >missing_site_threshold     #
    # fraction missing site-wide                                          #
    # ------------------------------------------------------------------ #
    missing_per_variant = np.isnan(raw).mean(axis=1)  # shape (n_variants,)
    keep_variants_mask = missing_per_variant <= missing_site_threshold
    raw = raw[keep_variants_mask]  # shape (n_kept_variants, n_samples)
    ploidy = ploidy[keep_variants_mask]
    kind = kind[keep_variants_mask]

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
    # Case/control diagnostic: meaningless for quantitative outcomes that happen to be 0/1.
    if is_binary and phenotype_vector is not None and raw.shape[0] > 0:
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
    # MAF from all RETAINED samples combined (never stratified; #174)     #
    # ------------------------------------------------------------------ #
    retained = np.asarray(sample_mask, dtype=bool)
    _infer_missing_ploidy(ploidy, kind, retained)
    geno = raw.copy()
    n_kept_variants = geno.shape[0]
    mafs = np.zeros(n_kept_variants, dtype=np.float64)

    for v_idx in range(n_kept_variants):
        obs_mask = ~np.isnan(geno[v_idx]) & retained
        observed_dosages = geno[v_idx, obs_mask]
        observed_ploidies = ploidy[v_idx, obs_mask]
        if len(observed_dosages) == 0:
            # Every retained sample missing: MAF=0, impute 0
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
        if np.any(missing_mask):
            v = ploidy[v_idx, missing_mask] * alt_freq
            geno[v_idx, missing_mask] = np.round(v) if is_binary else v
        if alt_freq > 0.5:
            geno[v_idx] = ploidy[v_idx] - geno[v_idx]

    # ------------------------------------------------------------------ #
    # Step 6: Transpose to (n_samples, n_variants)                        #
    # ------------------------------------------------------------------ #
    geno = geno.T  # shape (n_samples, n_kept_variants)

    if return_keep_mask:
        return geno, mafs, sample_mask, warnings_list, keep_variants_mask
    return geno, mafs, sample_mask, warnings_list
