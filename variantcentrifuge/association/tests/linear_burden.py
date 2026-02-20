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

Note on ``effect_size`` column naming
--------------------------------------
The engine names the effect size column ``{test_name}_or`` for all tests.
For ``linear_burden``, this column contains the OLS beta (not an odds ratio).
Column renaming to ``linear_burden_beta`` is deferred to Phase 22 when
output column standardization is addressed.

Warning codes
-------------
``NO_GENOTYPE_MATRIX``
    ``genotype_matrix`` key not found in ``contingency_data``.
``NO_PHENOTYPE_OR_EMPTY_MATRIX``
    ``phenotype_vector`` is None or genotype matrix is empty.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

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

    @property
    def name(self) -> str:
        """Short identifier for test registry and output column prefixes."""
        return "linear_burden"

    def check_dependencies(self) -> None:
        """Raise ImportError if statsmodels is not available."""
        try:
            import statsmodels.api  # noqa: F401
        except ImportError as e:
            raise ImportError(
                "LinearBurdenTest requires statsmodels. Install with: pip install statsmodels"
            ) from e

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
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"warnings": ["NO_GENOTYPE_MATRIX"]},
            )

        geno: np.ndarray = contingency_data["genotype_matrix"]
        mafs: np.ndarray = contingency_data["variant_mafs"]
        phenotype: np.ndarray | None = contingency_data.get("phenotype_vector")
        covariate_matrix: np.ndarray | None = contingency_data.get("covariate_matrix")

        if phenotype is None or geno.shape[0] == 0 or geno.shape[1] == 0:
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"warnings": ["NO_PHENOTYPE_OR_EMPTY_MATRIX"]},
            )

        # Compute weighted burden score per sample: (n_samples,)
        weights = get_weights(mafs, config.variant_weights)
        burden = geno @ weights

        # Build design matrix: [intercept, burden, covariates...]
        design = sm.add_constant(burden.reshape(-1, 1))
        has_covariates = (
            covariate_matrix is not None
            and covariate_matrix.ndim == 2
            and covariate_matrix.shape[1] > 0
        )
        if has_covariates:
            design = np.column_stack([design, covariate_matrix])

        # Carrier statistics (useful for downstream interpretation)
        carriers = burden > 0
        n_carriers = int(carriers.sum())
        warning_codes: list[str] = []

        # Fit OLS
        try:
            fit_result = sm.OLS(phenotype, design).fit()
        except Exception as exc:
            logger.warning(f"Gene {gene}: OLS.fit() raised {type(exc).__name__}: {exc}")
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=None,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                n_cases=n_cases,
                n_controls=n_controls,
                n_variants=n_variants,
                extra={"warnings": ["OLS_FIT_FAILED"], "n_carriers": n_carriers},
            )

        # Extract burden coefficient (index 1: intercept at 0, burden at 1)
        beta = float(fit_result.params[1])
        se = float(fit_result.bse[1])
        p_value = float(fit_result.pvalues[1])

        # conf_int() returns ndarray shape (n_params, 2) — NOT DataFrame (Pitfall 1)
        ci = fit_result.conf_int()
        ci_lower = float(ci[1, 0])
        ci_upper = float(ci[1, 1])

        return TestResult(
            gene=gene,
            test_name=self.name,
            p_value=p_value,
            corrected_p_value=None,
            effect_size=beta,  # beta, not OR (linear test)
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            n_cases=n_cases,
            n_controls=n_controls,
            n_variants=n_variants,
            extra={
                "se": se,
                "n_carriers": n_carriers,
                "warnings": warning_codes,
            },
        )
