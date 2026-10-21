"""
Covariate pre-flight checks, run once before the association gene loop (#222).

Ill-conditioned covariates used to turn every logistic burden fit into NA without
stopping the run.  :func:`preflight_covariates` drops near-constant covariates and
fails up front, naming the columns, on collinearity and on covariates that separate
a binary phenotype.  :func:`check_convergence_failures` fails the stage when one
logistic convergence-failure reason hits most tested genes.
"""

from __future__ import annotations

import logging
import warnings
from collections import Counter

import numpy as np
import pandas as pd

logger = logging.getLogger("variantcentrifuge")

# A covariate with fewer non-modal values carries (almost) no information and can
# quasi-separate the phenotype by itself, forcing every gene into the Firth fallback.
MIN_NON_MODAL = 5
# Condition-number limit of the intercept + centred, unit-norm covariate design.
MAX_CONDITION = 1e8
# Prefix of the per-gene convergence-failure codes in the ``warnings`` column.
CONVERGENCE_PREFIX = "CONVERGENCE_"


def preflight_covariates(
    covariates: np.ndarray | None,
    names: list[str],
    phenotype: np.ndarray | None,
    binary: bool,
) -> tuple[np.ndarray | None, list[str]]:
    """
    Check the covariate matrix once; return it (unchanged object when nothing is dropped).

    Drops (with a warning) covariates with fewer than ``MIN_NON_MODAL`` non-modal
    values.  Raises ``ValueError`` naming the columns when the remaining design is
    rank-deficient or its condition number after centring and scaling exceeds
    ``MAX_CONDITION``, or (binary traits) a covariate separates cases from controls.
    """
    if covariates is None or covariates.shape[1] == 0:
        return covariates, names
    n_non_modal = [
        len(col) - int(np.unique(col, return_counts=True)[1].max()) for col in covariates.T
    ]
    keep = [j for j, k in enumerate(n_non_modal) if k >= MIN_NON_MODAL]
    if len(keep) < covariates.shape[1]:
        dropped = [f"{names[j]} ({n_non_modal[j]})" for j in range(len(names)) if j not in keep]
        logger.warning(
            "Covariate pre-flight: dropped %d covariate(s) with fewer than %d values "
            "differing from the most common value (non-modal count): %s",
            len(dropped),
            MIN_NON_MODAL,
            ", ".join(dropped),
        )
        covariates, names = covariates[:, keep], [names[j] for j in keep]
        if not keep:
            return None, []

    centred = covariates - covariates.mean(axis=0)
    scaled = centred / np.linalg.norm(centred, axis=0)
    n = len(scaled)
    _u, sv, vt = np.linalg.svd(np.column_stack([np.full(n, n**-0.5), scaled]))
    condition = sv[0] / sv[-1] if sv[-1] > 0 else np.inf
    if condition > MAX_CONDITION:
        loadings = np.abs(vt[-1, 1:])
        involved = [nm for nm, v in zip(names, loadings, strict=True) if v >= 0.1 * loadings.max()]
        raise ValueError(
            f"Covariates are rank-deficient or collinear (condition number {condition:.3g} "
            f"after centring and scaling, limit {MAX_CONDITION:.0e}): {', '.join(involved)}. "
            "Remove or combine these columns."
        )

    if binary and phenotype is not None and np.ptp(phenotype) > 0:
        cases, controls = covariates[phenotype == 1], covariates[phenotype == 0]
        separating = [
            nm
            for nm, lo_ca, hi_ca, lo_co, hi_co in zip(
                names,
                cases.min(axis=0),
                cases.max(axis=0),
                controls.min(axis=0),
                controls.max(axis=0),
                strict=True,
            )
            if hi_ca <= lo_co or hi_co <= lo_ca
        ]
        if separating:
            raise ValueError(
                "Covariate(s) separate cases from controls, so no logistic model "
                f"converges: {', '.join(separating)}. Remove these columns."
            )
    return covariates, names


def check_convergence_failures(
    results_df: pd.DataFrame,
    names: list[str],
    phenotype: np.ndarray | None,
    covariates: np.ndarray | None,
) -> None:
    """
    Log logistic burden convergence-failure reasons; raise when one hits > 50 % of genes.

    Genes tested = genes with a logistic burden p-value plus genes that failed to
    converge.  The error lists every reason count and the covariates whose null-model
    standard error exceeds the burden test's separation threshold.
    """
    if "logistic_burden_pvalue" not in results_df.columns:
        return
    per_gene = [
        [c for c in str(w).split(";") if c.startswith(CONVERGENCE_PREFIX)]
        for w in results_df["warnings"]
    ]
    reasons = Counter(c for codes in per_gene for c in codes)
    if not reasons:
        return
    n_failed = sum(1 for codes in per_gene if codes)
    n_tested = int(results_df["logistic_burden_pvalue"].notna().sum()) + n_failed
    summary = ", ".join(f"{code}={count}" for code, count in reasons.most_common())
    logger.warning(
        "logistic_burden: %d of %d tested genes failed to converge (%s)",
        n_failed,
        n_tested,
        summary,
    )
    if reasons.most_common(1)[0][1] <= 0.5 * n_tested:
        return
    raise RuntimeError(
        f"logistic_burden failed to converge for {n_failed} of {n_tested} tested genes "
        f"({summary}); suspect covariates: {_suspect_covariates(names, phenotype, covariates)}. "
        "Check the covariates (scale, near-constant columns, separation) before rerunning."
    )


def _suspect_covariates(
    names: list[str], phenotype: np.ndarray | None, covariates: np.ndarray | None
) -> str:
    """Name covariates with a non-finite or separation-sized SE in the null logistic fit."""
    if covariates is None or phenotype is None:
        return "none (no covariates)"
    import statsmodels.api as sm

    from variantcentrifuge.association.tests.logistic_burden import _SEPARATION_BSE_THRESHOLD

    design = np.column_stack([np.ones(len(covariates)), covariates])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            bse = np.asarray(sm.Logit(phenotype, design).fit(disp=False, maxiter=100).bse)[1:]
        except np.linalg.LinAlgError:
            bse = np.full(len(names), np.inf)
    suspects = [
        nm for nm, se in zip(names, bse, strict=True) if not se <= _SEPARATION_BSE_THRESHOLD
    ]
    return ", ".join(suspects) or f"none identified (covariates: {', '.join(names)})"
