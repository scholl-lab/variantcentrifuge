# File: variantcentrifuge/association/tests/logistic_burden.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/logistic_burden.py
"""
Logistic burden test with Firth penalized likelihood fallback.

Implements a weighted burden collapsing test for binary phenotypes using
statsmodels.Logit with covariate adjustment. When perfect or quasi-separation
is detected (or convergence fails), automatically retries with Firth's
penalized logistic regression.

Separation detection
--------------------
- statsmodels Logit may not raise an exception on separation. Detection uses:
  - ``result.mle_retvals.get("converged", True)`` (False = not converged)
  - ``result.bse.max() > 100.0`` (BSE inflation indicates separation)
  Both checks together catch the cases where statsmodels emits a
  PerfectSeparationWarning but returns a result object anyway.

conf_int() return type
----------------------
In statsmodels 0.14.4, ``LogitResults.conf_int()`` returns a numpy ndarray of
shape ``(n_params, 2)``, NOT a DataFrame. Use ``ci[1, 0]`` and ``ci[1, 1]`` for
the burden coefficient (index 1, after intercept at index 0).

Warning codes
-------------
Structured codes added to ``TestResult.extra["warnings"]``:

``PERFECT_SEPARATION``
    All carriers are cases or all are controls.
``QUASI_SEPARATION``
    >=90% of carriers are in one group.
``ZERO_CARRIERS_ONE_GROUP``
    No carriers in cases or no carriers in controls.
``LOW_CARRIER_COUNT``
    Fewer than 3 carriers total.
``FIRTH_CONVERGE_FAIL``
    Firth fallback attempted and failed to converge.
``NO_GENOTYPE_MATRIX``
    ``genotype_matrix`` key not found in ``contingency_data``.
"""

from __future__ import annotations

import logging
import warnings
from typing import Any

import numpy as np

from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult

logger = logging.getLogger("variantcentrifuge")

# BSE threshold above which separation is assumed (Pitfall 2 from RESEARCH.md)
_SEPARATION_BSE_THRESHOLD = 100.0


def _firth_loglik(beta: np.ndarray, y: np.ndarray, x: np.ndarray) -> float:
    """
    Firth penalized log-likelihood: log-lik + 0.5 * log(det(I(beta))).

    Parameters
    ----------
    beta : np.ndarray, shape (k,)
    y : np.ndarray, shape (n,), binary 0/1
    x : np.ndarray, shape (n, k)

    Returns
    -------
    float
        Penalized log-likelihood. Returns plain log-likelihood if Fisher
        information matrix is singular (det <= 0).
    """
    pi = 1.0 / (1.0 + np.exp(-x @ beta))
    ll = float(np.sum(y * np.log(pi + 1e-15) + (1 - y) * np.log(1 - pi + 1e-15)))
    w = pi * (1 - pi)
    xtwx = x.T @ (w[:, None] * x)
    sign, logdet = np.linalg.slogdet(xtwx)
    if sign <= 0:
        return ll
    return ll + 0.5 * logdet


def _firth_logistic(
    y: np.ndarray,
    x: np.ndarray,
    max_iter: int = 25,
    tol: float = 1e-4,
) -> object | None:
    """
    Firth penalized logistic regression via Newton-Raphson with step-halving.

    Implements the Jeffreys prior penalty to obtain finite estimates under
    complete or quasi-complete separation. Self-contained — no external
    Firth package required.

    Algorithm adapted from:
    John Lees, https://gist.github.com/johnlees/3e06380965f367e4894ea20fbae2b90d

    Parameters
    ----------
    y : np.ndarray, shape (n,)
        Binary outcome vector (0/1).
    x : np.ndarray, shape (n, k)
        Design matrix (including intercept column).
    max_iter : int
        Maximum Newton-Raphson iterations. Default: 25.
    tol : float
        Convergence tolerance on max(|delta_beta|). Default: 1e-4.

    Returns
    -------
    object | None
        Namespace with ``.params``, ``.bse``, ``.pvalues``, and
        ``.conf_int(alpha)`` method — mimics statsmodels result interface.
        Returns ``None`` if Fisher information matrix is singular at any step
        or convergence fails.

    Notes
    -----
    The returned object's ``.conf_int()`` returns a numpy ndarray of shape
    ``(k, 2)`` (Wald confidence intervals), consistent with how statsmodels
    returns CIs. Callers extract ``ci[1, 0]`` and ``ci[1, 1]`` for the
    burden coefficient.
    """
    from scipy.stats import norm as sp_norm

    _n, k = x.shape
    beta = np.zeros(k)

    for _iteration in range(max_iter):
        pi = 1.0 / (1.0 + np.exp(-x @ beta))
        w = pi * (1 - pi)

        # Fisher information matrix
        xtwx = x.T @ (w[:, None] * x)
        try:
            var_cov = np.linalg.inv(xtwx)
        except np.linalg.LinAlgError:
            return None

        # Hat matrix diagonal (memory-efficient: diag(W^0.5 X V X^T W^0.5))
        wsqrt_x = np.sqrt(w)[:, None] * x
        h_diag = np.sum(wsqrt_x @ var_cov * wsqrt_x, axis=1)

        # Penalized score (Jeffreys prior gradient)
        u_pen = x.T @ (y - pi + h_diag * (0.5 - pi))
        step = var_cov @ u_pen
        beta_new = beta + step

        # Step-halving: ensure penalized log-likelihood increases
        old_ll = _firth_loglik(beta, y, x)
        for _halving in range(10):
            new_ll = _firth_loglik(beta_new, y, x)
            if new_ll > old_ll:
                break
            step *= 0.5
            beta_new = beta + step

        converged = bool(np.max(np.abs(beta_new - beta)) < tol)
        beta = beta_new
        if converged:
            break

    # Final Fisher information at converged beta
    pi = 1.0 / (1.0 + np.exp(-x @ beta))
    w = pi * (1 - pi)
    xtwx = x.T @ (w[:, None] * x)
    try:
        var_cov = np.linalg.inv(xtwx)
    except np.linalg.LinAlgError:
        return None

    bse = np.sqrt(np.diag(var_cov))
    z = beta / bse
    pvals = 2.0 * (1.0 - sp_norm.cdf(np.abs(z)))

    # Capture in closure for the result object
    _beta = beta.copy()
    _bse = bse.copy()
    _pvals = pvals.copy()
    _var_cov = var_cov.copy()

    class FirthResult:
        """Minimal namespace mimicking statsmodels LogitResults interface."""

        params = _beta
        bse = _bse
        pvalues = _pvals

        def conf_int(self, alpha: float = 0.05) -> np.ndarray:
            """
            Wald confidence intervals.

            Returns
            -------
            np.ndarray, shape (k, 2)
                ``ci[i, 0]`` = lower bound for parameter ``i``;
                ``ci[i, 1]`` = upper bound.
            """
            from scipy.stats import norm as _n

            z_val = _n.ppf(1.0 - alpha / 2.0)
            return np.column_stack([_beta - z_val * _bse, _beta + z_val * _bse])

    return FirthResult()


class LogisticBurdenTest(AssociationTest):
    """
    Weighted burden collapsing test for binary traits using logistic regression.

    Collapses per-sample variant genotypes into a single burden score using
    Beta(MAF; 1, 25) weights (or user-specified scheme), then fits logistic
    regression with optional covariate adjustment. Reports beta coefficient
    (log-odds) + SE, matching SKAT/SAIGE-GENE convention.

    Firth fallback is triggered automatically when:
    - Standard Logit fails to converge (``mle_retvals["converged"] == False``)
    - OR standard errors exceed 100 (separation indicator)

    Input requirements (via ``contingency_data``)
    ---------------------------------------------
    ``genotype_matrix`` : np.ndarray, shape (n_samples, n_variants)
        Imputed genotype dosage matrix (0/1/2, no NaN).
    ``variant_mafs`` : np.ndarray, shape (n_variants,)
        Per-variant minor allele frequencies (for weight computation).
    ``phenotype_vector`` : np.ndarray, shape (n_samples,)
        Binary phenotype (0 = control, 1 = case).
    ``covariate_matrix`` : np.ndarray, shape (n_samples, k), or None
        Covariate matrix (already encoded, no intercept column).

    If ``genotype_matrix`` is absent from ``contingency_data``, returns
    ``p_value=None`` with warning code ``NO_GENOTYPE_MATRIX``.
    """

    @property
    def name(self) -> str:
        """Short identifier for test registry and output column prefixes."""
        return "logistic_burden"

    def effect_column_names(self) -> dict[str, str]:
        """Beta/SE column naming for logistic burden (log-odds scale)."""
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
                "LogisticBurdenTest requires statsmodels. Install with: pip install statsmodels"
            ) from e

    def run(
        self,
        gene: str,
        contingency_data: dict[str, Any],
        config: AssociationConfig,
    ) -> TestResult:
        """
        Run logistic burden test for one gene.

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
            is the beta coefficient (log-odds) for the burden score.
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
                se=None,
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
                se=None,
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

        # Carrier statistics
        carriers = burden > 0
        n_carriers = int(carriers.sum())
        n_carriers_cases = int((carriers & (phenotype == 1)).sum())
        n_carriers_controls = int((carriers & (phenotype == 0)).sum())

        # Pre-flight warning codes
        warning_codes: list[str] = []
        if n_carriers < 3:
            warning_codes.append("LOW_CARRIER_COUNT")
        if n_carriers > 0:
            case_rate = n_carriers_cases / n_carriers
            if case_rate >= 0.999 or case_rate <= 0.001:
                warning_codes.append("PERFECT_SEPARATION")
            elif case_rate >= 0.9 or case_rate <= 0.1:
                warning_codes.append("QUASI_SEPARATION")
        if (n_carriers_cases == 0 or n_carriers_controls == 0) and (
            "PERFECT_SEPARATION" not in warning_codes
        ):
            warning_codes.append("ZERO_CARRIERS_ONE_GROUP")

        # Fit standard logistic regression
        try:
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                fit_result = sm.Logit(phenotype, design).fit(disp=False, maxiter=100)

            converged = bool(fit_result.mle_retvals.get("converged", True))
            bse_max = float(fit_result.bse.max())
            needs_firth = (not converged) or (bse_max > _SEPARATION_BSE_THRESHOLD)

        except Exception as exc:
            logger.debug(f"Gene {gene}: Logit.fit() raised {type(exc).__name__}: {exc}")
            needs_firth = True
            fit_result = None

        if needs_firth:
            no_sep_warning = (
                "PERFECT_SEPARATION" not in warning_codes
                and "QUASI_SEPARATION" not in warning_codes
            )
            if no_sep_warning:
                warning_codes.append("PERFECT_SEPARATION")
            logger.debug(f"Gene {gene}: Logit non-convergence/separation — trying Firth fallback")
            firth_result = _firth_logistic(phenotype, design, max_iter=config.firth_max_iter)
            if firth_result is None:
                warning_codes.append("FIRTH_CONVERGE_FAIL")
                logger.warning(
                    f"Gene {gene}: Both Logit and Firth failed to converge — reporting NA"
                )
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
                    extra={
                        "warnings": warning_codes,
                        "n_carriers": n_carriers,
                        "n_carriers_cases": n_carriers_cases,
                        "n_carriers_controls": n_carriers_controls,
                    },
                )
            fit_result = firth_result

        # Extract burden coefficient (index 1: intercept at 0, burden at 1)
        beta = float(fit_result.params[1])
        se = float(fit_result.bse[1])
        p_value = float(fit_result.pvalues[1])

        # conf_int() returns ndarray shape (n_params, 2) — NOT DataFrame (Pitfall 1)
        ci = fit_result.conf_int()
        ci_lower = float(ci[1, 0])
        ci_upper = float(ci[1, 1])

        # Report beta (log-odds) + SE directly — matching SKAT/SAIGE-GENE convention.
        # OR = exp(beta) is NOT reported because per-unit burden OR is non-intuitive
        # when Beta(MAF;1,25) weights make each variant contribute ~24 units.
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
                "beta": beta,
                "se": se,
                "n_carriers": n_carriers,
                "n_carriers_cases": n_carriers_cases,
                "n_carriers_controls": n_carriers_controls,
                "warnings": warning_codes,
            },
        )
