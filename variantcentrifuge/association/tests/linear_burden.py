# File: variantcentrifuge/association/tests/linear_burden.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/linear_burden.py
"""
Linear burden test for quantitative phenotypes.

Implements a weighted burden collapsing test for quantitative traits using
ordinary least squares (statsmodels.OLS) with optional covariate adjustment.
Reports the beta coefficient (burden effect size), standard error, and 95% CI.

conf_int() return type
----------------------
In statsmodels 0.14.4, ``OLSResults.conf_int()`` returns a numpy ndarray of
shape ``(n_params, 2)``, NOT a DataFrame. Use ``ci[1, 0]`` and ``ci[1, 1]`` for
the burden coefficient (index 1, after intercept at index 0).

Warning codes
-------------
``NO_GENOTYPE_MATRIX``
    ``genotype_matrix`` key not found in ``contingency_data``.
``NO_PHENOTYPE_OR_EMPTY_MATRIX``
    ``phenotype_vector`` is None or genotype matrix is empty.
``NONFINITE_PHENOTYPE`` / ``CONSTANT_PHENOTYPE``
    Outcome contains NaN/inf or has zero variance; no regression is fitted.
``NONFINITE_COVARIATES``
    Covariate matrix contains NaN/inf.
``BURDEN_NOT_ESTIMABLE``
    Burden is constant or lies in the span of the intercept and covariates.
``OUTCOME_EXPLAINED_BY_COVARIATES``
    Covariates reproduce the outcome exactly (no residual variation to test).
``ILL_CONDITIONED_COVARIATES``
    Covariates too ill-conditioned to decide whether the burden is estimable.
``OLS_FIT_FAILED`` / ``NONFINITE_STATISTICS``
    OLS raised, or returned a non-finite beta, SE, p-value or CI.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

from variantcentrifuge.association.backends.base import NullModelResult
from variantcentrifuge.association.backends.null_model import (
    explained_columns,
    full_rank_design,
    outcome_explained,
    outcome_problem,
)
from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult

logger = logging.getLogger("variantcentrifuge")


class LinearBurdenTest(AssociationTest):
    """
    Weighted burden collapsing test for quantitative traits using OLS.

    Collapses per-sample variant genotypes into a single burden score using
    Beta(MAF; 1, 25) weights (or user-specified scheme), then fits ordinary
    least squares regression with optional covariate adjustment. Reports the
    beta coefficient (effect size) + SE + 95% CI.

    Input requirements (via ``contingency_data``)
    ---------------------------------------------
    ``genotype_matrix`` : np.ndarray, shape (n_samples, n_variants)
        Imputed genotype dosage matrix (0/1/2, no NaN).
    ``variant_mafs`` : np.ndarray, shape (n_variants,)
        Per-variant minor allele frequencies (for weight computation).
    ``phenotype_vector`` : np.ndarray, shape (n_samples,)
        Continuous phenotype values.
    ``covariate_matrix`` : np.ndarray, shape (n_samples, k), or None
        Covariate matrix (already encoded, no intercept column).

    If ``genotype_matrix`` is absent from ``contingency_data``, returns
    ``p_value=None`` with warning code ``NO_GENOTYPE_MATRIX``.
    """

    parallel_safe: bool = True  # Thread-safe -- pure statsmodels, no R

    @property
    def name(self) -> str:
        """Short identifier for test registry and output column prefixes."""
        return "linear_burden"

    def effect_column_names(self) -> dict[str, str | None]:
        """Beta/SE column naming for linear burden (OLS coefficient scale)."""
        return {
            "effect": "beta",
            "se": "se",
            "ci_lower": "beta_ci_lower",
            "ci_upper": "beta_ci_upper",
        }

    def check_dependencies(self) -> None:
        """Raise ImportError if statsmodels is not available."""
        try:
            import statsmodels.api  # noqa: F401
        except ImportError as e:
            raise ImportError(
                "LinearBurdenTest requires statsmodels. Install with: pip install statsmodels"
            ) from e

    def _skip(
        self,
        gene: str,
        codes: list[str],
        n_cases: int,
        n_controls: int,
        n_variants: int,
        n_carriers: int | None = None,
    ) -> TestResult:
        """Return a skipped result (no p-value) with warning ``codes``."""
        extra: dict[str, Any] = {"warnings": codes}
        if n_carriers is not None:
            extra["n_carriers"] = n_carriers
        return TestResult(
            gene=gene,
            test_name=self.name,
            p_value=None,
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            se=None,
            n_cases=n_cases,
            n_controls=n_controls,
            n_variants=n_variants,
            extra=extra,
        )

    def run(
        self,
        gene: str,
        contingency_data: dict[str, Any],
        config: AssociationConfig,
    ) -> TestResult:
        """
        Run linear burden test for one gene.

        Parameters
        ----------
        gene : str
        contingency_data : dict
            Must contain ``genotype_matrix``, ``variant_mafs``,
            ``phenotype_vector``. Optional: ``covariate_matrix``.
        config : AssociationConfig

        Returns
        -------
        TestResult
            ``p_value=None`` when test is skipped or fails. ``effect_size``
            is the OLS beta coefficient for the burden score.
        """
        import statsmodels.api as sm

        from variantcentrifuge.association.weights import get_weights

        n_cases = int(contingency_data.get("proband_count", 0))
        n_controls = int(contingency_data.get("control_count", 0))
        n_variants = int(contingency_data.get("n_qualifying_variants", 0))

        # Require genotype matrix — not present for Fisher-only pipelines
        if "genotype_matrix" not in contingency_data:
            return self._skip(gene, ["NO_GENOTYPE_MATRIX"], n_cases, n_controls, n_variants)

        geno: np.ndarray = contingency_data["genotype_matrix"]
        mafs: np.ndarray = contingency_data["variant_mafs"]
        phenotype: np.ndarray | None = contingency_data.get("phenotype_vector")
        covariate_matrix: np.ndarray | None = contingency_data.get("covariate_matrix")

        if phenotype is None or geno.shape[0] == 0 or geno.shape[1] == 0:
            return self._skip(
                gene, ["NO_PHENOTYPE_OR_EMPTY_MATRIX"], n_cases, n_controls, n_variants
            )
        # Non-finite or zero-variance outcomes have no regression to test (#186)
        problem = outcome_problem(phenotype)
        if problem is not None:
            return self._skip(gene, [problem], n_cases, n_controls, n_variants)

        # Compute weighted burden score per sample: (n_samples,)
        weights = get_weights(
            mafs,
            config.variant_weights,
            cadd_scores=contingency_data.get("cadd_scores"),
            revel_scores=contingency_data.get("revel_scores"),
            score_values=contingency_data.get("score_values"),
            variant_effects=contingency_data.get("variant_effects"),
            weight_params=config.variant_weight_params,
        )
        burden = geno @ weights

        # Carrier statistics (useful for downstream interpretation)
        carriers = burden > 0
        n_carriers = int(carriers.sum())
        warning_codes: list[str] = []

        if covariate_matrix is not None and (
            covariate_matrix.ndim != 2 or covariate_matrix.shape[1] == 0
        ):
            covariate_matrix = None
        if covariate_matrix is not None and not np.isfinite(covariate_matrix).all():
            return self._skip(
                gene, ["NONFINITE_COVARIATES"], n_cases, n_controls, n_variants, n_carriers
            )

        # One nuisance basis and rank policy for estimability AND fitting: the
        # SKAT null-model basis [1, Q] (intercept explicit, redundant covariates
        # dropped).  Burden keeps its units, so beta/SE equal lm(y ~ b + x).
        n_obs = len(burden)
        nuisance, covariate_r = full_rank_design(covariate_matrix, n_obs)
        nuisance_model = NullModelResult(
            model=None,
            trait_type="quantitative",
            n_samples=n_obs,
            extra={"design_matrix": nuisance, "covariate_r": covariate_r},
        )
        if outcome_explained(phenotype, nuisance_model, phenotype):
            return self._skip(
                gene,
                ["OUTCOME_EXPLAINED_BY_COVARIATES"],
                n_cases,
                n_controls,
                n_variants,
                n_carriers,
            )
        explained = explained_columns(burden[:, np.newaxis], nuisance_model)
        if explained is None:
            return self._skip(
                gene, ["ILL_CONDITIONED_COVARIATES"], n_cases, n_controls, n_variants, n_carriers
            )
        if explained[0]:  # constant burden or burden in span(1, covariates)
            return self._skip(
                gene, ["BURDEN_NOT_ESTIMABLE"], n_cases, n_controls, n_variants, n_carriers
            )
        design = np.column_stack([nuisance[:, :1], burden, nuisance[:, 1:]])

        # Fit OLS
        try:
            fit_result = sm.OLS(phenotype, design).fit()
        except Exception as exc:
            logger.warning(f"Gene {gene}: OLS.fit() raised {type(exc).__name__}: {exc}")
            return self._skip(gene, ["OLS_FIT_FAILED"], n_cases, n_controls, n_variants, n_carriers)

        # Extract burden coefficient (index 1: intercept at 0, burden at 1)
        beta = float(fit_result.params[1])
        se = float(fit_result.bse[1])
        p_value = float(fit_result.pvalues[1])

        # conf_int() returns ndarray shape (n_params, 2) — NOT DataFrame (Pitfall 1)
        ci = fit_result.conf_int()
        ci_lower = float(ci[1, 0])
        ci_upper = float(ci[1, 1])
        if not np.isfinite([beta, se, p_value, ci_lower, ci_upper]).all():
            return self._skip(
                gene, ["NONFINITE_STATISTICS"], n_cases, n_controls, n_variants, n_carriers
            )

        return TestResult(
            gene=gene,
            test_name=self.name,
            p_value=p_value,
            corrected_p_value=None,
            effect_size=beta,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            se=se,
            n_cases=n_cases,
            n_controls=n_controls,
            n_variants=n_variants,
            extra={
                "se": se,
                "n_carriers": n_carriers,
                "warnings": warning_codes,
            },
        )
