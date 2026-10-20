"""
Inverse-variance meta of per-stratum Firth burden log-ORs (#198, step 1).

At emission each stratum fits, per gene, a Firth-penalised logistic regression of
case status on the weighted ALT burden ``G @ w`` plus the stratum's covariates
(``w`` from the stratum's own MAF with the run's ``beta:a,b``/``uniform`` scheme) and
stores the burden log-OR and its Wald SE. The meta pools those with a fixed-effect
inverse-variance estimate: an alternative burden p-value that stays well behaved
under strong case:control imbalance and separation, where the score test's normal
approximation does not. It never feeds ``primary_p`` or the multiple-testing family.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import scipy.stats

from variantcentrifuge.association.backends.null_model import full_rank_design
from variantcentrifuge.association.meta.align import pooled_frequency_weights
from variantcentrifuge.association.tests.logistic_burden import _firth_logistic

logger = logging.getLogger("variantcentrifuge")

FIRTH_COLUMNS = [
    "burden_firth_ivw_beta",
    "burden_firth_ivw_se",
    "burden_firth_ivw_p",
    "burden_firth_ivw_n_strata",
]


def stratum_firth_burden(
    geno: np.ndarray,
    ac: np.ndarray,
    an: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    weights: str,
) -> dict[str, Any] | None:
    """Firth burden fit for one gene in one stratum: ``{"beta", "se", "weights"}``.

    ``geno`` is the ALT-coded dosage matrix behind the stratum's score (so the sign
    matches ``meta_burden``); ``ac``/``an`` give the stratum MAF for the weights.
    Returns None when the scheme is not reproducible at meta time, the burden has no
    variance, or the fit does not converge.
    """
    try:
        w = pooled_frequency_weights(ac, an, weights)
    except ValueError:
        return None
    burden = np.asarray(geno, dtype=np.float64) @ w
    if np.ptp(burden) == 0:
        return None
    # Full-rank nuisance basis [1, Q] (drops constant/collinear covariates, as the null
    # model does); the burden log-OR and SE are invariant to this reparametrisation.
    nuisance, _r = full_rank_design(covariates, len(burden))
    design = np.column_stack([nuisance[:, :1], burden, nuisance[:, 1:]])
    # Tighter than the single-cohort 1e-4 (1e-6 misses R parity at rel 1e-6); 1e-8 can stall
    # step-halving on log-likelihood rounding noise in large strata.
    fit = _firth_logistic(np.asarray(phenotype, dtype=np.float64), design, max_iter=50, tol=1e-7)
    if fit is None or not np.isfinite(fit.bse[1]):  # type: ignore[attr-defined]
        return None
    return {"beta": float(fit.params[1]), "se": float(fit.bse[1]), "weights": weights}  # type: ignore[attr-defined]


def ivw_meta(
    betas: list[float], ses: list[float]
) -> tuple[float | None, float | None, float | None]:
    """Fixed-effect inverse-variance (beta, se, two-sided p) over finite strata."""
    b, s = np.asarray(betas, dtype=np.float64), np.asarray(ses, dtype=np.float64)
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    if not ok.any():
        return None, None, None
    w = 1.0 / s[ok] ** 2
    beta = float(np.sum(w * b[ok]) / np.sum(w))
    se = float(1.0 / np.sqrt(np.sum(w)))
    return beta, se, float(2.0 * scipy.stats.norm.sf(abs(beta / se)))


def _scheme_matches(firth: dict[str, Any], scheme: str) -> bool:
    return str(firth.get("weights", "")).strip().lower() == scheme.strip().lower()


def warn_firth_exclusions(artifacts: list[Any], weight_scheme: str) -> None:
    """Log one warning per stratum whose genes leave the Firth meta (no fit / other weights)."""
    for art in artifacts:
        fits = [gs.firth for gs in art.genes.values()]
        missing = sum(f is None for f in fits)
        other = sorted({f["weights"] for f in fits if f and not _scheme_matches(f, weight_scheme)})
        n_other = sum(1 for f in fits if f and not _scheme_matches(f, weight_scheme))
        if missing or n_other:
            logger.warning(
                "Firth burden meta: stratum %s contributes no fit for %d of %d genes "
                "(%d with no Firth fit: failed or not stored; %d fitted with weights %s, "
                "not --weights %s)",
                art.stratum_id,
                missing + n_other,
                len(fits),
                missing,
                n_other,
                ",".join(other) or "-",
                weight_scheme,
            )


def firth_ivw_columns(gene_scores: list[Any], weight_scheme: str) -> dict[str, Any]:
    """Output columns from the strata whose Firth fit used the meta's weight scheme.

    Strata without a fit (v1 artifacts, failed fits) or fitted with another scheme
    are left out; ``burden_firth_ivw_n_strata`` counts the strata used.
    """
    fits = [
        gs.firth
        for gs in gene_scores
        if gs.firth is not None and _scheme_matches(gs.firth, weight_scheme)
    ]
    beta, se, p = ivw_meta([f["beta"] for f in fits], [f["se"] for f in fits])
    return dict(zip(FIRTH_COLUMNS, (beta, se, p, len(fits)), strict=True))


__all__ = [
    "FIRTH_COLUMNS",
    "firth_ivw_columns",
    "ivw_meta",
    "stratum_firth_burden",
    "warn_firth_exclusions",
]
