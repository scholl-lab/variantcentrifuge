"""
SKAT null-model fitting and score-test projection for the Python backend.

The null model and every score statistic built on it share one full-rank
design basis: an intercept plus the covariate columns that are not linearly
dependent on it or on each other.  Using that basis for both the GLM fit and
the nuisance projection keeps redundant covariates from changing results.
"""

from __future__ import annotations

import logging

import numpy as np
import scipy.linalg

from variantcentrifuge.association.backends.base import NullModelResult

logger = logging.getLogger("variantcentrifuge")

# Pivot tolerance for collinearity among centred, unit-norm covariates (R lm()/glm() qr tol).
_RANK_TOL = 1e-7


class UntestableGeneError(ValueError):
    """A gene cannot be tested reliably; ``reason`` is its skip code."""

    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


def outcome_problem(phenotype: np.ndarray) -> str | None:
    """Return a skip-reason code when an outcome cannot be tested, else None."""
    y = np.asarray(phenotype, dtype=np.float64)
    if y.size == 0 or not np.isfinite(y).all():
        return "NONFINITE_PHENOTYPE"
    if np.ptp(y) == 0.0:
        return "CONSTANT_PHENOTYPE"
    return None


def full_rank_design(covariates: np.ndarray | None, n_obs: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Return an orthonormal-covariate basis ``[1, Q]`` for ``[1, covariates]`` and ``R``.

    Each covariate is first referenced to its first observation, so an exactly
    representable translation leaves the data bit-identical; a column is
    constant only if that leaves it exactly zero.  Varying columns are centred,
    scaled to unit norm and rank-selected by column-pivoted QR (R ``lm`` tol
    1e-7); ``Q`` is the orthonormal factor of the kept columns.  ``[1, Q]`` is
    used for the GLM fit, the residual df and every projection.  ``R`` (with
    ``Q R`` = the kept normalized covariates) lets callers recover the
    conditioning of the ORIGINAL covariates, which limits how precisely ``Q``
    spans them.  Dropped columns are logged once per fit.
    """
    ones = np.ones((n_obs, 1), dtype=np.float64)
    if covariates is None or covariates.ndim != 2 or covariates.shape[1] == 0:
        return ones, np.empty((0, 0))
    cov = np.asarray(covariates, dtype=np.float64)
    if not np.isfinite(cov).all():
        raise ValueError("covariate matrix must be finite")

    shifted = cov - cov[0]
    varies = np.flatnonzero(np.any(shifted != 0.0, axis=0))
    basis = np.empty((n_obs, 0))
    kept: list[int] = []
    r_sel = np.empty((0, 0))
    if len(varies) > 0:
        centred = shifted[:, varies] - shifted[:, varies].mean(axis=0)
        scaled = centred / np.linalg.norm(centred, axis=0)
        q, r, piv = scipy.linalg.qr(scaled, mode="economic", pivoting=True)
        diag = np.abs(np.diag(r))
        rank = int(np.sum(diag > _RANK_TOL * diag[0]))
        basis, r_sel = q[:, :rank], r[:rank, :rank]
        kept = sorted(varies[piv[:rank]].tolist())
    if len(kept) < cov.shape[1]:
        dropped = sorted(set(range(cov.shape[1])) - set(kept))
        logger.warning(
            "SKAT null model: dropped %d covariate column(s) that are constant or "
            "linearly dependent on other covariates (0-based): %s",
            len(dropped),
            dropped,
        )
    return np.column_stack([ones, basis]), r_sel


def outcome_explained(y: np.ndarray, nuisance: NullModelResult, residuals: np.ndarray) -> bool:
    """
    True when the covariates leave no testable outcome variation.

    Quantitative: ``||P y|| / ||y||`` is within the projection round-off bound
    (same reasoning as :func:`explained_columns`; unit invariant).  Binary: every
    response residual is below 0.5, i.e. the fitted linear predictor separates
    cases from controls, so no finite null MLE exists.
    """
    if nuisance.trait_type == "binary":
        return bool(np.all(np.abs(residuals) < 0.5))
    y = np.asarray(y, dtype=np.float64)
    # Exact translation (the design has an intercept): keeps offsets out of ||y||.
    adjusted, z, bound = _project(y - y[0], nuisance)
    return bool(np.linalg.norm(adjusted) <= bound * np.linalg.norm(z))


def fit_glm_null_model(
    phenotype: np.ndarray,
    covariates: np.ndarray | None,
    trait_type: str,
) -> NullModelResult:
    """
    Fit the SKAT null model (Binomial or Gaussian GLM) on a full-rank design.

    Matches R ``SKAT_Null_Model`` (``Adjustment=FALSE``): response residuals
    ``y - mu_hat``; for quantitative traits ``sigma2 = RSS / (n - rank(X))``
    (``summary(lm)$sigma^2``); binary traits use ``sigma2 = 1``.

    Raises
    ------
    ValueError
        Non-finite or constant outcome, or non-finite covariates.
    UntestableGeneError
        ``INSUFFICIENT_RESIDUAL_DF`` (samples <= independent design columns) or
        ``OUTCOME_EXPLAINED_BY_COVARIATES``; callers skip the affected gene/sample set.
    """
    import statsmodels.api as sm

    problem = outcome_problem(phenotype)
    if problem is not None:
        raise ValueError(
            f"Cannot fit SKAT null model: {problem} "
            "(outcome must be finite and vary across samples)"
        )
    y = np.asarray(phenotype, dtype=np.float64)
    n_obs = len(y)
    design_matrix, covariate_r = full_rank_design(covariates, n_obs)
    df_resid = n_obs - design_matrix.shape[1]
    if df_resid <= 0:
        logger.warning(
            f"Cannot fit SKAT null model: {n_obs} samples for "
            f"{design_matrix.shape[1]} independent design columns"
        )
        raise UntestableGeneError("INSUFFICIENT_RESIDUAL_DF")

    family = sm.families.Binomial() if trait_type == "binary" else sm.families.Gaussian()
    # The design has an intercept, so shifting a quantitative outcome by its first value is
    # an exact, model-equivalent translation that keeps large offsets out of the round-off.
    shift = y[0] if trait_type == "quantitative" else 0.0
    fit_result = sm.GLM(y - shift, design_matrix, family=family).fit()

    residuals = np.asarray(fit_result.resid_response, dtype=np.float64)
    nuisance = NullModelResult(
        model=None,
        trait_type=trait_type,
        n_samples=n_obs,
        extra={"design_matrix": design_matrix, "covariate_r": covariate_r},
    )
    if outcome_explained(y - shift, nuisance, residuals):
        raise UntestableGeneError("OUTCOME_EXPLAINED_BY_COVARIATES")
    sigma2 = float(residuals @ residuals) / df_resid if trait_type == "quantitative" else 1.0
    mu_hat = np.asarray(fit_result.fittedvalues, dtype=np.float64) + shift
    logger.debug(f"Null model fit: trait_type={trait_type}, n_samples={n_obs}, sigma2={sigma2:.6f}")

    return NullModelResult(
        model=fit_result,
        trait_type=trait_type,
        n_samples=n_obs,
        adjustment=False,
        extra={
            "residuals": residuals,
            "sigma2": sigma2,
            "mu_hat": mu_hat,
            "design_matrix": design_matrix,
            "covariate_r": covariate_r,
        },
    )


# Davies p-values are computed to acc = 1e-6; a projected column known only to a
# worse relative accuracy cannot support a test, and one whose projected part is
# below a round-off bound of at most 1e-6 is negligible at that accuracy.
_DECISION_ACC = 1e-6


def _project(
    geno_weighted: np.ndarray, null_model: NullModelResult
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return (projected Z, variance-weighted Z, relative round-off bound of P z_j)."""
    z: np.ndarray = geno_weighted if geno_weighted.ndim == 2 else geno_weighted[:, np.newaxis]
    design_x = np.asarray(null_model.extra["design_matrix"], dtype=np.float64)
    if null_model.trait_type == "binary":
        mu_hat = np.asarray(null_model.extra["mu_hat"], dtype=np.float64)
        phi = np.sqrt(mu_hat * (1.0 - mu_hat))[:, np.newaxis]
        z = phi * z
        design_x = phi * design_x
    basis, r_factor = np.linalg.qr(design_x)
    adjusted: np.ndarray = z - basis @ (basis.T @ z)
    # Conditioning of the normalized, (variance-)weighted ORIGINAL nuisance design
    # V^1/2 [1, C] = basis @ r_factor @ blockdiag(1, R_cov), from exact singular values.
    r_cov = np.asarray(null_model.extra.get("covariate_r", np.empty((0, 0))), dtype=np.float64)
    original = r_factor.copy()
    if r_cov.size:
        original[:, 1:] = r_factor[:, 1:] @ r_cov
    kappa = float(np.linalg.cond(original / np.linalg.norm(original, axis=0)))
    # Round-off bound for ||P z_j|| / ||z_j|| (n samples, k nuisance columns):
    #  Householder QR / projection of length-n columns has backward error
    #  ~ lambda sqrt(n) eps per column (probabilistic analysis, Higham & Mary 2019;
    #  lambda = sqrt(2 ln(2n / 1e-6)), failure probability ~1e-6).  It enters
    #  (i) directly via the k basis columns: lambda k sqrt(n) eps, and (ii) through
    #  the tilt of span(Q) from the original covariates, which Wedin's theorem
    #  amplifies by kappa: lambda sqrt(k n) eps kappa.  Not a certified bound:
    #  a scan of >1000 exactly-in-span constructions (n 300-10000, k 2-8,
    #  separations 2^-8..2^-30, both traits; 1,794 constructions) peaked at
    #  0.36 sqrt(k n) eps kappa (max ratio/bound 0.054); consumers treat
    #  "indeterminate" fail-safe (skip / UntestableGeneError, never zeroing).
    # ponytail: conservative, so strongly near-collinear covariates at large n
    # (e.g. kappa ~1e7, n = 2e5) make genes skip as ILL_CONDITIONED_COVARIATES even
    # when variants are plainly informative; a certified tighter bound would recover them.
    n, k = design_x.shape
    lam = np.sqrt(2.0 * np.log(2.0 * n / _DECISION_ACC))
    bound = float(np.finfo(np.float64).eps * lam * np.sqrt(n) * (k + np.sqrt(k) * kappa))
    return adjusted, z, bound


def project_out_nuisance(geno_weighted: np.ndarray, null_model: NullModelResult) -> np.ndarray:
    """
    Return ``A`` with ``A.T @ A = Z.T P Z`` for the null-model score covariance.

    Binary: ``P = V - V X (X'VX)^-1 X'V`` with ``V = diag(mu (1 - mu))``, computed
    as the residual of ``V^1/2 Z`` against an orthonormal QR basis of ``V^1/2 X``.
    Quantitative: ``P = I - X (X'X)^-1 X'`` (score covariance is ``sigma2 * P``).
    """
    adjusted = _project(geno_weighted, null_model)[0]
    return adjusted if geno_weighted.ndim == 2 else adjusted[:, 0]


def projection_noise_floor(geno_weighted: np.ndarray, null_model: NullModelResult) -> float:
    """
    Singular-value floor below which ``project_out_nuisance`` output is round-off.

    The per-column round-off bound of :func:`_project` times ``sigma_max`` of the
    UNprojected (variance-weighted) matrix.  Use after :func:`explained_columns`
    has zeroed individually explained columns, so they cannot inflate ``sigma_max``.
    """
    _adjusted, z, bound = _project(geno_weighted, null_model)
    return float(bound * np.linalg.norm(z, 2))


def explained_columns(geno: np.ndarray, null_model: NullModelResult) -> np.ndarray | None:
    """
    Classify unweighted genotype columns against the nuisance (covariate) space.

    With ``ratio_j = ||P z_j|| / ||z_j||`` and ``bound`` the projection round-off
    bound (see :func:`_project`):

    - explained: ``ratio_j <= bound <= 1e-6`` -- callers zero the column;
    - informative: ``bound <= 1e-6 * ratio_j`` -- its projection is accurate to 1e-6;
    - otherwise numerically indeterminate: returns ``None`` and callers skip the
      gene (``ILL_CONDITIONED_COVARIATES``) instead of zeroing anything.

    Zero columns count as explained.  Returns the boolean explained mask.
    """
    adjusted, z, bound = _project(geno, null_model)
    z_norm = np.linalg.norm(z, axis=0)
    ratio = np.linalg.norm(adjusted, axis=0) / np.maximum(z_norm, np.finfo(np.float64).tiny)
    explained = ((ratio <= bound) & (bound <= _DECISION_ACC)) | (z_norm == 0.0)
    informative = bound <= _DECISION_ACC * ratio
    if not np.all(explained | informative):
        return None
    return np.asarray(explained, dtype=bool)
