# File: variantcentrifuge/association/backends/coast_python.py
# Location: variantcentrifuge/variantcentrifuge/association/backends/coast_python.py
"""
Pure Python COAST backend implementing the allelic series association test.

PythonCOASTBackend implements the Combined Omnibus Association Test (COAST,
McCaw et al. AJHG 2023) using only Python-native scientific libraries — no R
or rpy2 required. The implementation decomposes into three layers:

Layer 1 — Burden component (6 tests):
  Three models (baseline, sum, max) x two encodings (count, indicator)
  - Baseline: 3-df joint chi-squared over all three category coefficients
  - Sum: 1-df regression on weighted sum of category scores
  - Max: 1-df regression on maximum weighted category score
  Continuous traits: OLS (Wald). Binary traits: LRT (deviance difference).

Layer 2 — Allelic SKAT-O (1 test):
  SKAT-Optimal (Lee et al. 2012) with annotation-aware weights:
    w_j = sqrt(coast_weight_{A_j} / (aaf_j * (1 - aaf_j)))
  where A_j is the annotation category of variant j. SKAT-O optimises over
  a rho grid [0, 0.01, ..., 1] blending SKAT (rho=0) and burden (rho=1).
  This matches R AllelicSeries::ASKAT() which calls SKAT::SKAT(method="SKATO").

Layer 3 — Cauchy omnibus (1 combined p-value):
  7 component p-values combined via cauchy_combination() with weights
  [1, 1, 1, 1, 1, 1, 6] (equal burden weights, 6x SKAT weight to achieve
  50/50 burden vs SKAT evidence split).

Weight conventions:
  - Baseline burden: UNIFORM weights [1, 1, 1] (NOT coast_weights)
  - Sum/Max burden: coast_weights (default [1.0, 2.0, 3.0])
  - Allelic SKAT: per-variant sqrt(coast_weight / aaf*(1-aaf))

Thread safety:
  PythonCOASTBackend is thread-safe — pure Python/numpy/scipy, no rpy2.

References:
  McCaw et al. (2023) AJHG. "An allelic-series rare-variant association test
  for candidate-gene discovery."
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import numpy as np
import scipy.stats

if TYPE_CHECKING:
    from variantcentrifuge.association.backends.base import NullModelResult
    from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

logger = logging.getLogger("variantcentrifuge")

# Default COAST Cauchy combination weights: 6 burden tests x weight 1, 1 SKAT x weight 6.
# Achieves 50/50 burden vs SKAT evidence split (6 x 1 = 6 x 1 = 6).
_COAST_CAUCHY_WEIGHTS = [1, 1, 1, 1, 1, 1, 6]

# Category codes used by classify_variants()
_BMV_CODE = 1
_DMV_CODE = 2
_PTV_CODE = 3
_CATEGORY_ORDER = [_BMV_CODE, _DMV_CODE, _PTV_CODE]

# Minimum variance threshold to skip monomorphic variants in SKAT weight computation
_MIN_VARIANT_VARIANCE = 1e-8


def _aggregate_by_category(
    geno: np.ndarray,
    anno_codes: np.ndarray,
    weights: list[float],
    indicator: bool = False,
    method: str = "none",
) -> np.ndarray:
    """
    Aggregate genotype dosages by annotation category.

    For each COAST category (BMV=1, DMV=2, PTV=3), sums genotype dosages across
    all variants in that category per sample. Optionally converts to carrier
    indicators and applies further aggregation.

    Parameters
    ----------
    geno : np.ndarray, shape (n_samples, n_variants)
        Genotype dosage matrix (0/1/2 encoded).
    anno_codes : np.ndarray, shape (n_variants,)
        Integer annotation codes: 1=BMV, 2=DMV, 3=PTV.
    weights : list of float
        Per-category weights in order [BMV_weight, DMV_weight, PTV_weight].
        Used for ``method="sum"`` and ``method="max"``; ignored for ``method="none"``.
    indicator : bool
        If True, convert per-category allele counts to 0/1 carrier indicators.
    method : str
        - ``"none"``: return 3-column matrix (n_samples, 3) — baseline 3-df test
        - ``"sum"``: return weighted sum as (n_samples, 1)
        - ``"max"``: return element-wise max over weighted categories as (n_samples, 1)

    Returns
    -------
    np.ndarray
        Shape (n_samples, 3) for method="none", or (n_samples, 1) for sum/max.
    """
    n_samples = geno.shape[0]
    # Build per-category columns (n_samples, 3)
    cat_matrix = np.zeros((n_samples, 3), dtype=np.float64)
    for col_idx, cat_code in enumerate(_CATEGORY_ORDER):
        variant_mask = anno_codes == cat_code
        if np.any(variant_mask):
            cat_matrix[:, col_idx] = geno[:, variant_mask].sum(axis=1)

    if indicator:
        cat_matrix = (cat_matrix > 0).astype(np.float64)

    if method == "none":
        return cat_matrix

    # Apply per-category weights for sum/max aggregation
    w_arr = np.array(weights, dtype=np.float64)  # shape (3,)

    if method == "sum":
        # Weighted sum: scalar per sample
        weighted_sum = cat_matrix @ w_arr  # (n_samples,)
        return weighted_sum[:, np.newaxis]  # (n_samples, 1)

    if method == "max":
        # Element-wise max of weighted categories: max(w1*N1, w2*N2, w3*N3)
        weighted_cats = cat_matrix * w_arr[np.newaxis, :]  # (n_samples, 3)
        max_vals: np.ndarray = weighted_cats.max(axis=1)  # (n_samples,)
        return max_vals[:, np.newaxis]  # (n_samples, 1)

    raise ValueError(f"Unknown aggregation method: {method!r}. Expected 'none', 'sum', or 'max'.")


def _run_burden_test(
    predictor: np.ndarray,
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    trait_type: str,
) -> float | None:
    """
    Run a single burden regression test.

    Fits a regression of phenotype on predictor (+ covariates + intercept).
    For baseline (3-df) tests, uses an F-test or LRT comparing all three
    burden coefficients to zero jointly.

    Parameters
    ----------
    predictor : np.ndarray, shape (n_samples, k)
        Predictor matrix. k=3 for baseline (3-df test); k=1 for sum/max (1-df).
    phenotype : np.ndarray, shape (n_samples,)
        Phenotype vector (binary 0/1 or continuous).
    covariates : np.ndarray or None, shape (n_samples, p)
        Covariate matrix. None = intercept only.
    trait_type : str
        "quantitative" -> OLS with Wald/F-test; "binary" -> Logit with LRT.

    Returns
    -------
    float or None
        P-value for the joint burden effect, or None on numerical failure.
    """
    import statsmodels.api as sm

    n_samples = len(phenotype)
    n_burden_cols = predictor.shape[1]

    # Build full design matrix: intercept + predictor + covariates
    design_with_burden = sm.add_constant(predictor, prepend=True, has_constant="add")
    if covariates is not None and covariates.ndim == 2 and covariates.shape[1] > 0:
        design_with_burden = np.hstack([design_with_burden, covariates])

    # Build null design matrix: intercept + covariates only (no burden predictor)
    if covariates is not None and covariates.ndim == 2 and covariates.shape[1] > 0:
        design_null = np.hstack([np.ones((n_samples, 1)), covariates])
    else:
        design_null = np.ones((n_samples, 1), dtype=np.float64)

    if trait_type == "quantitative":
        try:
            full_fit = sm.OLS(phenotype, design_with_burden).fit()
        except Exception as exc:
            logger.debug(f"OLS burden test failed: {exc}")
            return None

        if n_burden_cols == 1:
            # 1-df test: direct p-value from coefficient
            return float(full_fit.pvalues[1])
        else:
            # 3-df test: F-test on burden coefficients (indices 1..n_burden_cols)
            # Restriction matrix: one row per burden coefficient, selecting that coefficient
            n_total_params = design_with_burden.shape[1]
            r_matrix = np.zeros((n_burden_cols, n_total_params), dtype=np.float64)
            for i in range(n_burden_cols):
                r_matrix[i, i + 1] = 1.0  # +1 to skip the intercept at index 0
            try:
                f_test = full_fit.f_test(r_matrix)
                return float(f_test.pvalue)
            except Exception as exc:
                logger.debug(f"OLS F-test failed: {exc}")
                return None

    else:  # binary
        try:
            null_fit = sm.Logit(phenotype, design_null).fit(disp=False, method="bfgs", maxiter=100)
            full_fit = sm.Logit(phenotype, design_with_burden).fit(
                disp=False, method="bfgs", maxiter=100
            )
        except Exception as exc:
            logger.debug(f"Logit burden test failed (convergence): {exc}")
            return None

        # Likelihood ratio test
        lrt = 2.0 * (full_fit.llf - null_fit.llf)
        lrt = max(lrt, 0.0)  # Clamp to non-negative (numerical noise)
        return float(scipy.stats.chi2.sf(lrt, df=n_burden_cols))


def _compute_allelic_skat_weights(
    geno_filtered: np.ndarray,
    anno_codes_filtered: np.ndarray,
    coast_weights: list[float],
) -> np.ndarray:
    """
    Compute per-variant SKAT weights for the allelic series test.

    Implements: w_j = sqrt(coast_weight_{A_j} / (aaf_j * (1 - aaf_j)))

    Note: variance denominator uses aaf*(1-aaf), NOT 2*aaf*(1-aaf).
    This matches the AllelicSeries R source (GenVar computes var(G)/2 implicitly
    by the factor of 2 in the SKAT Q statistic formula Q = score'score/2).

    Parameters
    ----------
    geno_filtered : np.ndarray, shape (n_samples, n_variants)
        Genotype matrix for COAST-eligible variants (BMV/DMV/PTV only).
    anno_codes_filtered : np.ndarray, shape (n_variants,)
        Annotation codes for filtered variants (1/2/3 only, no 0s).
    coast_weights : list of float
        Category weights in order [BMV_weight, DMV_weight, PTV_weight].

    Returns
    -------
    np.ndarray, shape (n_variants,)
        Per-variant SKAT weights (square root of weight-to-variance ratio).
    """
    n_variants = geno_filtered.shape[1]
    aaf = geno_filtered.mean(axis=0) / 2.0
    aaf = np.clip(aaf, 1e-8, 1.0 - 1e-8)

    # Per-variant annotation-based weight
    cat_weights = np.zeros(n_variants, dtype=np.float64)
    for col_idx, cat_code in enumerate(_CATEGORY_ORDER):
        mask = anno_codes_filtered == cat_code
        if np.any(mask):
            cat_weights[mask] = coast_weights[col_idx]

    # Variance: aaf * (1 - aaf)  — NOT 2*aaf*(1-aaf)
    variance = aaf * (1.0 - aaf)
    variance = np.maximum(variance, _MIN_VARIANT_VARIANCE)

    result: np.ndarray = np.sqrt(cat_weights / variance)
    return result


class PythonCOASTBackend:
    """
    Pure Python COAST backend (no R/rpy2 dependency).

    Implements the full COAST algorithm (McCaw et al. AJHG 2023):
    - 6 burden tests (3 models x 2 encodings)
    - 1 allelic SKAT-O test
    - Cauchy combination of all 7 p-values

    The null model (from PythonSKATBackend.fit_null_model()) is passed in
    by PurePythonCOASTTest, which fits it lazily and caches it.

    Thread safety
    -------------
    Thread-safe — pure Python/numpy/scipy, no rpy2 or shared mutable state.
    """

    def __init__(self) -> None:
        self._skat_backend: PythonSKATBackend | None = None

    def _ensure_skat_backend(self) -> None:
        """Lazily initialize the PythonSKATBackend (used for eigenvalue computation)."""
        if self._skat_backend is None:
            from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

            self._skat_backend = PythonSKATBackend()
            self._skat_backend.detect_environment()

    def _run_allelic_skat(
        self,
        geno_filtered: np.ndarray,
        skat_weights: np.ndarray,
        null_model: NullModelResult,
    ) -> float:
        """
        Run the allelic series SKAT-O component.

        R AllelicSeries::ASKAT() uses SKAT::SKAT(method="SKATO"), which runs
        SKAT-Optimal (Lee et al. 2012) — an optimized blend of SKAT and burden
        over a rho grid.  This method mirrors that behaviour using the existing
        ``_skato_get_pvalue`` infrastructure from PythonSKATBackend.

        Parameters
        ----------
        geno_filtered : np.ndarray, shape (n_samples, n_variants)
            Genotype matrix for COAST-eligible variants.
        skat_weights : np.ndarray, shape (n_variants,)
            Per-variant SKAT weights from _compute_allelic_skat_weights().
        null_model : NullModelResult
            Fitted null model from PythonSKATBackend.fit_null_model().

        Returns
        -------
        float
            P-value for the allelic SKAT-O component (1.0 if degenerate).
        """
        from variantcentrifuge.association.backends.python_backend import (
            _SKATO_RHO_GRID,
            _skato_get_pvalue,
        )

        n_variants = geno_filtered.shape[1]
        if n_variants == 0:
            return 1.0

        residuals: np.ndarray = null_model.extra["residuals"]
        trait_type: str = null_model.trait_type

        # Apply per-variant allelic-series weights
        geno_weighted = geno_filtered * skat_weights[np.newaxis, :]

        # Rank check — SKAT-O needs rank >= 2
        rank = int(np.linalg.matrix_rank(geno_weighted))
        if rank < 2:
            logger.debug("PythonCOASTBackend allelic SKAT-O: rank < 2; p=1.0")
            return 1.0

        # ── Step 1: Compute Q_rho for each rho ──
        score_vec = geno_weighted.T @ residuals  # (n_variants,)
        score_sum = float(np.sum(score_vec))

        q_rho = np.zeros(len(_SKATO_RHO_GRID), dtype=np.float64)
        for i, rho in enumerate(_SKATO_RHO_GRID):
            q1 = (1.0 - rho) * float(score_vec @ score_vec)
            q2 = rho * score_sum**2
            q_rho[i] = (q1 + q2) / 2.0

        # ── Step 2: Build projection-adjusted Z1 ──
        if trait_type == "binary":
            mu_hat: np.ndarray = null_model.extra["mu_hat"]
            design_x: np.ndarray = null_model.extra["design_matrix"]
            pi_1 = mu_hat * (1.0 - mu_hat)
            sqrt_pi = np.sqrt(pi_1)

            z_phi = geno_weighted * sqrt_pi[:, np.newaxis]
            x_phi = design_x * sqrt_pi[:, np.newaxis]
            x_pi = design_x * pi_1[:, np.newaxis]
            xtx_pi = x_phi.T @ x_phi
            xtx_pi_inv = np.linalg.inv(xtx_pi)
            z1_adj = z_phi - x_phi @ (xtx_pi_inv @ (x_pi.T @ geno_weighted))
        else:
            design_x = null_model.extra["design_matrix"]
            sigma2: float = null_model.extra["sigma2"]
            xtx = design_x.T @ design_x
            xtx_inv = np.linalg.inv(xtx)
            z1_adj = geno_weighted - design_x @ (xtx_inv @ (design_x.T @ geno_weighted))
            q_rho /= sigma2

        z1_half = z1_adj / np.sqrt(2.0)

        # ── Step 3-5: Full SKAT-O omnibus p-value ──
        pvalue, pval_each = _skato_get_pvalue(q_rho, z1_half, _SKATO_RHO_GRID)

        best_idx = int(np.argmin(pval_each))
        best_rho = _SKATO_RHO_GRID[best_idx]
        logger.debug(f"PythonCOASTBackend allelic SKAT-O: p={pvalue:.6e}, rho={best_rho:.3f}")

        if np.isnan(pvalue) or np.isinf(pvalue):
            return 1.0

        return float(np.clip(pvalue, 0.0, 1.0))

    def test_gene(
        self,
        gene: str,
        geno_filtered: np.ndarray,
        anno_codes_filtered: np.ndarray,
        phenotype: np.ndarray,
        covariates: np.ndarray | None,
        coast_weights: list[float],
        trait_type: str,
        null_model: NullModelResult,
    ) -> dict[str, Any]:
        """
        Run the full COAST test for a single gene.

        Computes 6 burden p-values + 1 allelic SKAT p-value, then combines all
        7 via Cauchy combination with weights [1, 1, 1, 1, 1, 1, 6].

        Parameters
        ----------
        gene : str
            Gene symbol (for logging).
        geno_filtered : np.ndarray, shape (n_samples, n_variants)
            Genotype matrix for COAST-eligible variants (BMV/DMV/PTV only).
        anno_codes_filtered : np.ndarray, shape (n_variants,)
            Annotation codes for filtered variants (1=BMV, 2=DMV, 3=PTV).
        phenotype : np.ndarray, shape (n_samples,)
            Phenotype vector (binary 0/1 or continuous).
        covariates : np.ndarray or None, shape (n_samples, k)
            Covariate matrix. None = intercept only.
        coast_weights : list of float
            Category weights [BMV_weight, DMV_weight, PTV_weight].
        trait_type : str
            "binary" or "quantitative".
        null_model : NullModelResult
            Fitted null model (from PythonSKATBackend.fit_null_model()).

        Returns
        -------
        dict with keys:
            - p_value : float | None — omnibus Cauchy p-value
            - burden_p_values : list[float | None] — 6 burden component p-values
            - skat_p_value : float | None — allelic SKAT p-value
            - burden_labels : list[str] — labels for the 6 burden components
        """
        from variantcentrifuge.association.tests.acat import cauchy_combination

        # Uniform weights for baseline (3-df) tests: equal per-category weight
        uniform_weights = [1.0, 1.0, 1.0]

        # ── Layer 1: 6 Burden tests ──────────────────────────────────────────────

        # Baseline models: 3-df tests (uniform weights, aggregate but don't weight)
        base_count = _aggregate_by_category(
            geno_filtered, anno_codes_filtered, uniform_weights, indicator=False, method="none"
        )
        base_ind = _aggregate_by_category(
            geno_filtered, anno_codes_filtered, uniform_weights, indicator=True, method="none"
        )

        # Sum models: 1-df tests (coast_weights applied)
        sum_count = _aggregate_by_category(
            geno_filtered, anno_codes_filtered, coast_weights, indicator=False, method="sum"
        )
        sum_ind = _aggregate_by_category(
            geno_filtered, anno_codes_filtered, coast_weights, indicator=True, method="sum"
        )

        # Max models: 1-df tests (coast_weights applied)
        max_count = _aggregate_by_category(
            geno_filtered, anno_codes_filtered, coast_weights, indicator=False, method="max"
        )
        max_ind = _aggregate_by_category(
            geno_filtered, anno_codes_filtered, coast_weights, indicator=True, method="max"
        )

        burden_p_values: list[float | None] = []
        burden_labels = [
            "base_count",
            "base_ind",
            "sum_count",
            "sum_ind",
            "max_count",
            "max_ind",
        ]

        for predictor, label in zip(
            [base_count, base_ind, sum_count, sum_ind, max_count, max_ind],
            burden_labels,
            strict=True,
        ):
            p = _run_burden_test(predictor, phenotype, covariates, trait_type)
            burden_p_values.append(p)
            logger.debug(f"COAST [{gene}] burden {label}: p={p}")

        # ── Layer 2: Allelic SKAT ────────────────────────────────────────────────

        skat_weights = _compute_allelic_skat_weights(
            geno_filtered, anno_codes_filtered, coast_weights
        )

        skat_p_value: float | None
        try:
            skat_p_value = self._run_allelic_skat(geno_filtered, skat_weights, null_model)
        except Exception as exc:
            logger.debug(f"COAST [{gene}] allelic SKAT failed: {exc}")
            skat_p_value = None

        logger.debug(f"COAST [{gene}] allelic SKAT: p={skat_p_value}")

        # ── Layer 3: Cauchy omnibus combination ──────────────────────────────────

        all_7_pvals: list[float | None] = [*burden_p_values, skat_p_value]
        omnibus_p = cauchy_combination(all_7_pvals, weights=_COAST_CAUCHY_WEIGHTS)

        logger.debug(
            f"COAST [{gene}] omnibus: p={omnibus_p} "
            f"(burden: {burden_p_values}, skat: {skat_p_value})"
        )

        return {
            "p_value": omnibus_p,
            "burden_p_values": burden_p_values,
            "skat_p_value": skat_p_value,
            "burden_labels": burden_labels,
        }
