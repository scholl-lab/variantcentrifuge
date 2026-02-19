# File: variantcentrifuge/association/correction.py
# Location: variantcentrifuge/variantcentrifuge/association/correction.py
"""
Multiple testing correction for association analyses.

Provides apply_correction() — a standalone wrapper around
statsmodels.stats.multitest that supports FDR (Benjamini-Hochberg) and
Bonferroni correction, with graceful degradation when statsmodels is absent.

This module is intentionally leaf-level: it imports only stdlib, numpy, and
statsmodels. No imports from gene_burden.py or other variantcentrifuge modules.
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

    if smm is None:
        logger.warning("statsmodels not available. Skipping multiple testing correction.")
        return pvals_array

    if method == "bonferroni":
        corrected: np.ndarray = smm.multipletests(pvals_array, method="bonferroni")[1]
    else:
        corrected = smm.multipletests(pvals_array, method="fdr_bh")[1]

    return corrected
