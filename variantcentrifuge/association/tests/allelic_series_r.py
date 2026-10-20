# File: variantcentrifuge/association/tests/allelic_series_r.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/allelic_series_r.py
"""
R AllelicSeries::COAST() invocation bridge via rpy2.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from variantcentrifuge.association.base import AssociationConfig

logger = logging.getLogger("variantcentrifuge")

# COAST-specific R GC interval (same as RSKATTest)
_GC_INTERVAL = 100


def check_r_dependencies() -> None:
    """Verify that rpy2 is importable and the R AllelicSeries package is installed.

    Raises
    ------
    ImportError
        If rpy2 is not installed, R is not found, or the AllelicSeries
        R package is not installed.
    """
    try:
        import rpy2.robjects  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "COASTTest requires rpy2. Install with: pip install rpy2\n"
            "Also ensure R is installed and R_HOME is set correctly."
        ) from exc

    try:
        from rpy2.robjects.packages import importr

        importr("AllelicSeries")
    except Exception as exc:
        raise ImportError(
            "COASTTest requires the R AllelicSeries package. "
            "Install in R with: install.packages('AllelicSeries')\n"
            "Or from GitHub: remotes::install_github('zhangz19/AllelicSeries')\n"
            f"Original error: {exc}"
        ) from exc


def _extract_pval(pvals_dict: dict[str, float | None], keys: list[str]) -> float | None:
    """Extract p-value from dict by trying multiple key names."""
    for key in keys:
        if key in pvals_dict:
            return pvals_dict[key]
    return None


def run_r_coast(
    geno_filtered: np.ndarray,
    anno_codes_filtered: np.ndarray,
    phenotype: np.ndarray,
    covariate_matrix: np.ndarray | None,
    config: AssociationConfig,
) -> tuple[float | None, float | None, float | None]:
    """Call R AllelicSeries::COAST() via rpy2 and extract p-values.

    Parameters
    ----------
    geno_filtered : np.ndarray (n_samples, n_variants_filtered)
        Genotype matrix for COAST-eligible variants only.
    anno_codes_filtered : np.ndarray (n_variants_filtered,)
        Integer annotation codes (1=BMV, 2=DMV, 3=PTV).
    phenotype : np.ndarray (n_samples,)
        Phenotype vector.
    covariate_matrix : np.ndarray or None
        Covariate matrix (n_samples, k) or None.
    config : AssociationConfig
        Runtime configuration.

    Returns
    -------
    (burden_p, skat_p, omnibus_p)
        Three p-values from COAST Pvals slot.
    """
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    allelic_series = importr("AllelicSeries")

    n_samples, n_variants_filtered = geno_filtered.shape

    # Build R matrix: AllelicSeries expects genotype as (n_samples, n_variants)
    r_geno = ro.r.matrix(
        ro.FloatVector(geno_filtered.T.ravel()),
        nrow=n_samples,
        ncol=n_variants_filtered,
    )
    r_anno = ro.IntVector(anno_codes_filtered.tolist())
    r_pheno = ro.FloatVector(phenotype.tolist())

    coast_weights = getattr(config, "coast_weights", None) or [1.0, 2.0, 3.0]
    r_weights = ro.FloatVector(coast_weights)

    is_binary = getattr(config, "trait_type", "binary") == "binary"
    r_is_binary = ro.BoolVector([is_binary])

    # Build covariate argument
    if covariate_matrix is not None and covariate_matrix.ndim == 2:
        n_cov = covariate_matrix.shape[1]
        r_covar = ro.r.matrix(
            ro.FloatVector(covariate_matrix.T.ravel()),
            nrow=n_samples,
            ncol=n_cov,
        )
    else:
        r_covar = ro.NULL

    result = allelic_series.COAST(
        anno=r_anno,
        geno=r_geno,
        pheno=r_pheno,
        covar=r_covar,
        weights=r_weights,
        is_pheno_binary=r_is_binary,
    )

    # Extract p-values from Pvals slot.
    pvals_slot = result.slots["Pvals"]

    # Convert R data.frame to Python dict keyed by test name
    pvals_dict: dict[str, float | None] = {}
    try:
        test_names = list(pvals_slot.rx2("test"))
        pval_values = list(pvals_slot.rx2("pval"))
        import rpy2.rinterface as ri

        for name, val in zip(test_names, pval_values, strict=False):
            if hasattr(ri, "NA_Real") and val is ri.NA_Real:
                pvals_dict[str(name)] = None
            elif hasattr(val, "__float__"):
                pvals_dict[str(name)] = float(val)
            else:
                pvals_dict[str(name)] = float(val) if val is not None else None
    except Exception:
        pvals_names = list(pvals_slot.names) if hasattr(pvals_slot, "names") else []
        pvals_values = list(pvals_slot)
        for name, val in zip(pvals_names, pvals_values, strict=False):
            try:
                pvals_dict[str(name)] = float(val) if val is not None else None
            except (TypeError, ValueError):
                pvals_dict[str(name)] = None

    logger.debug(f"COAST R Pvals: {pvals_dict}")

    burden_p = _extract_pval(
        pvals_dict,
        ["Burden", "Burden_p", "burden", "burden_p", "baseline"],
    )
    skat_p = _extract_pval(
        pvals_dict,
        ["SKAT", "SKAT_p", "skat", "skat_p", "allelic_skat"],
    )
    omnibus_p = _extract_pval(
        pvals_dict,
        ["Omni", "Omni_p", "omni", "omni_p", "Omnibus", "omnibus"],
    )

    if omnibus_p is None:
        omnibus_p = _extract_pval(pvals_dict, ["O", "o"])

    if omnibus_p is None and pvals_dict:
        last_key = list(pvals_dict.keys())[-1]
        omnibus_p = pvals_dict[last_key]
        logger.debug(
            f"COAST: omnibus p-value not found by standard keys; "
            f"using last key '{last_key}' = {omnibus_p}. "
            f"Available keys: {list(pvals_dict.keys())}"
        )

    return burden_p, skat_p, omnibus_p
