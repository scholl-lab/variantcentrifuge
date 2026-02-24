# File: variantcentrifuge/association/correction.py
# Location: variantcentrifuge/variantcentrifuge/association/correction.py
"""
Multiple testing correction for association analyses.

Provides:
- apply_correction() — standalone wrapper around statsmodels multipletests
  for standard (unweighted) FDR/Bonferroni correction.
- load_gene_weights() — read a TSV of per-gene prior weights for weighted BH.
- apply_weighted_correction() — weighted Benjamini-Hochberg (Genovese 2006)
  that renormalizes weights to mean=1.0 and applies weight-adjusted p-values.

This module is intentionally leaf-level: it imports only stdlib, numpy,
statsmodels, and pandas. No imports from gene_burden.py or other
variantcentrifuge modules.
"""

from __future__ import annotations

import logging

import numpy as np

try:
    import statsmodels.stats.multitest as smm
except ImportError:
    smm = None  # type: ignore[assignment]

logger = logging.getLogger("variantcentrifuge")


def apply_correction(
    pvals: list[float] | np.ndarray,
    method: str = "fdr",
) -> np.ndarray:
    """
    Apply multiple testing correction to a sequence of p-values.

    Produces output identical to the inline smm.multipletests calls in
    gene_burden.py (lines 506-510) for the same inputs and method.

    Parameters
    ----------
    pvals : list of float or np.ndarray
        Raw p-values to correct. Must be in [0, 1].
    method : str
        Correction method: "fdr" (Benjamini-Hochberg, default) or
        "bonferroni". Any other value is treated as "fdr".

    Returns
    -------
    np.ndarray
        Corrected p-values in the same order as input. If statsmodels is
        unavailable, returns the raw p-values unchanged (with a warning).
    """
    pvals_array = np.asarray(pvals, dtype=float)

    if len(pvals_array) == 0:
        return pvals_array

    if smm is None:
        logger.warning("statsmodels not available. Skipping multiple testing correction.")
        return pvals_array

    if method == "bonferroni":
        corrected: np.ndarray = smm.multipletests(pvals_array, method="bonferroni")[1]
    else:
        corrected = smm.multipletests(pvals_array, method="fdr_bh")[1]

    return corrected


def load_gene_weights(
    filepath: str,
    weight_column: str = "weight",
) -> dict[str, float]:
    """
    Read a TSV file of per-gene prior weights for weighted BH FDR correction.

    The file must have a header row. The first column is the gene identifier;
    the column named ``weight_column`` contains the weight values. All weights
    must be strictly positive (> 0).

    Parameters
    ----------
    filepath : str
        Path to the TSV file.
    weight_column : str
        Name of the column containing weight values. Default: "weight".

    Returns
    -------
    dict[str, float]
        Mapping from gene name (str) to raw weight (float). Weights are NOT
        normalized here — normalization happens in apply_weighted_correction()
        on the subset of actually testable genes.

    Raises
    ------
    ValueError
        If the weight column is missing, or if any weight is <= 0.
    """
    import pandas as pd

    df = pd.read_csv(filepath, sep="\t", header=0)

    # First column = gene identifier
    gene_col = df.columns[0]

    if weight_column not in df.columns:
        raise ValueError(
            f"Weight column '{weight_column}' not found in '{filepath}'. "
            f"Available columns: {list(df.columns)}"
        )

    gene_weights: dict[str, float] = {}
    bad_genes: list[str] = []

    for _, row in df.iterrows():
        gene = str(row[gene_col])
        w = float(row[weight_column])
        if w <= 0:
            bad_genes.append(f"{gene}={w}")
        gene_weights[gene] = w

    if bad_genes:
        shown = bad_genes[:10]
        n_more = len(bad_genes) - len(shown)
        detail = ", ".join(shown)
        if n_more > 0:
            detail += f" and {n_more} more"
        raise ValueError(
            f"Gene weights must be strictly positive (> 0). "
            f"Found zero or negative weights: {detail}"
        )

    logger.info(f"Loaded {len(gene_weights)} gene weights from {filepath}")
    return gene_weights


def apply_weighted_correction(
    pvals: list[float] | np.ndarray,
    genes: list[str],
    weight_map: dict[str, float],
    method: str = "fdr",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply weighted Benjamini-Hochberg FDR correction (Genovese et al. 2006).

    Weights are renormalized to mean=1.0 over the tested genes, then the
    adjusted p-values (p_i / w_i) are passed to the standard BH procedure.
    This preserves the FDR guarantee while increasing power for genes with
    high prior weights and decreasing it for genes with low prior weights.

    Parameters
    ----------
    pvals : list of float or np.ndarray
        Raw p-values aligned 1:1 with ``genes``.
    genes : list of str
        Gene names, one per p-value. Order must match ``pvals``.
    weight_map : dict[str, float]
        Mapping from gene name to raw weight (from load_gene_weights()).
        Genes not in the map receive weight=1.0 (neutral/default).
    method : str
        Correction method: "fdr" (default, Benjamini-Hochberg) or
        "bonferroni" (weighted Bonferroni: min(p_i * m / w_i, 1.0)).

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (corrected_pvals, normalized_weights)
        - corrected_pvals : array of corrected p-values, aligned with input.
        - normalized_weights : array of per-gene weights after mean=1.0
          renormalization, aligned with input genes.

    Notes
    -----
    - If statsmodels is unavailable, returns raw p-values unchanged (with warning).
    - Single gene input is returned unchanged (weighting has no effect for m=1).
    - Coverage warnings are emitted when >50% or >80% of tested genes are
      absent from weight_map.
    """
    pvals_array = np.asarray(pvals, dtype=float)
    m = len(pvals_array)

    # Uniform default normalized weights (will be overwritten below)
    default_norm = np.ones(m, dtype=float)

    if m == 0:
        return pvals_array, default_norm

    # Single gene edge case — weighting has no statistical effect
    if m == 1:
        logger.info("Single gene tested — FDR weighting has no effect, skipping")
        return pvals_array.copy(), default_norm

    # Look up raw weights; track coverage
    raw_weights = np.ones(m, dtype=float)
    n_missing = 0
    for i, gene in enumerate(genes):
        if gene in weight_map:
            raw_weights[i] = weight_map[gene]
        else:
            n_missing += 1

    # Coverage warnings
    missing_frac = n_missing / m
    if missing_frac > 0.80:
        logger.warning(
            f"STRONG WARNING: >80% of tested genes ({n_missing}/{m}) are absent "
            "from the weight file and will receive weight=1.0 (neutral). "
            "Weighted FDR may not reflect your biological priors."
        )
    elif missing_frac > 0.50:
        logger.warning(
            f"More than 50% of tested genes ({n_missing}/{m}) are absent from "
            "the weight file and will receive weight=1.0 (neutral)."
        )

    # Renormalize weights to mean=1.0: w_norm = raw * m / sum(raw)
    weight_sum = raw_weights.sum()
    if weight_sum == 0:
        # Pathological: all weights zero (should have been caught by load_gene_weights)
        logger.warning("All weights are zero after lookup; falling back to uniform weights.")
        w_norm = np.ones(m, dtype=float)
    else:
        w_norm = raw_weights * m / weight_sum

    # Warn about extreme weight ratios
    w_min = float(w_norm.min())
    w_max = float(w_norm.max())
    if w_min > 0 and w_max / w_min > 100:
        logger.warning(
            f"Extreme weight ratio detected: max/min = {w_max / w_min:.1f} "
            "(threshold: 100). Highly variable weights may distort FDR."
        )
    if w_max > 10:
        logger.warning(
            f"Largest normalized weight = {w_max:.2f} > 10. "
            "A single gene with very high weight can dominate the correction."
        )

    # Effective number of tests: n_eff = (sum w)^2 / sum(w^2)
    n_eff = float(w_norm.sum() ** 2 / (w_norm**2).sum())
    logger.info(
        f"FDR weighting: {m} genes tested, {n_missing} at default weight=1.0, "
        f"effective number of tests = {n_eff:.1f}"
    )

    # Graceful degradation when statsmodels is absent
    if smm is None:
        logger.warning("statsmodels not available. Skipping weighted multiple testing correction.")
        return pvals_array.copy(), w_norm

    # Apply weighted correction
    if method == "bonferroni":
        # Weighted Bonferroni: corrected_i = min(p_i * m / w_i, 1.0)
        corrected: np.ndarray = np.clip(pvals_array * m / w_norm, 0.0, 1.0)
    else:
        # Weighted BH: pass p_i / w_i to standard BH, then clip to [0, 1]
        p_adj = np.clip(pvals_array / w_norm, 0.0, 1.0)
        corrected = smm.multipletests(p_adj, method="fdr_bh")[1]

    return corrected, w_norm
