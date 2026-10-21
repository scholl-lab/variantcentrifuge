"""Empirical p-values for the top genes of a single-cohort binary run (#195).

Parametric bootstrap from a fixed null fit: y* ~ Bernoulli(mu_hat) of the covariate null
model, never refitted.  With ``A`` the nuisance-projected weighted genotypes
(``A.T @ A = Z' P Z``, see ``null_model.project_out_nuisance``) and ``E`` the standardised
bootstrap residuals ``(y* - mu) / sqrt(mu (1 - mu))``, a batch of B score vectors is one
GEMM ``A.T @ E``; they have the null covariance ``Z' P Z`` of the observed score.  Each
gene stops adaptively (Besag & Clifford 1991) after ``H_STOP`` exceedances.  Meta-analysis
permutation needs genotypes from every cohort and is not supported (#195).
"""

from __future__ import annotations

import hashlib
import logging
import os
import time
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np
import pandas as pd

from variantcentrifuge.association.backends.davies import compute_pvalue
from variantcentrifuge.association.backends.null_model import (
    UntestableGeneError,
    fit_glm_null_model,
    project_out_nuisance,
)
from variantcentrifuge.association.backends.skat_math import compute_eigenvalues_filtered
from variantcentrifuge.association.engine_workers import (
    fork_available,
    fork_map,
    single_blas_thread,
)
from variantcentrifuge.association.meta.artifact import _MAX_CACHED_NULLS
from variantcentrifuge.association.weights import get_weights

logger = logging.getLogger("variantcentrifuge")

H_STOP = 10  # Besag-Clifford: stop at the 10th replicate statistic >= the observed one
FIRST_BATCH = 128  # batches double up to MAX_BATCH: early stops stay cheap
MAX_BATCH = 1024  # ponytail: 2 x n x 1024 float64 per batch (~100 MB at 6,000 samples)
BURDEN_TESTS = ("logistic_burden", "linear_burden")

_STATISTICS: dict[str, Callable[[np.ndarray], np.ndarray]] = {
    "skat": lambda s: (s * s).sum(axis=0),  # SKAT Q (up to the constant 1/2)
    "burden": lambda s: s.sum(axis=0) ** 2,  # squared weighted burden score
}
METHOD_LABEL = {"skat": "bootstrap_skat_q", "burden": "bootstrap_burden_score"}


def permutation_statistic(
    test_names: Sequence[str], primary_id: str | None, skat_method: str, trait_type: str
) -> str:
    """Statistic behind the empirical p: "skat" (Q) or "burden" (score).

    The primary's own statistic when it is SKAT or burden; for ACAT-O / SKAT-O primaries
    (not rebuildable from the score matrix cheaply) SKAT Q if SKAT runs, else the burden
    score.  Raises ValueError when permutation cannot run.
    """
    if trait_type != "binary":
        raise ValueError("--permutations supports binary traits only")
    skat_stat = "burden" if skat_method == "Burden" else "skat"
    if primary_id in BURDEN_TESTS:
        return "burden"
    if {"skat", "skat_python"} & set(test_names):
        return skat_stat
    if "logistic_burden" in test_names:
        return "burden"
    raise ValueError("--permutations requires 'skat' or 'logistic_burden' in --association-tests")


def source_pvalue_column(stat: str, test_names: Sequence[str], columns: Sequence[str]) -> str:
    """Output column of the test whose statistic is permuted: genes are selected on it."""
    if stat == "skat":
        return "skat_o_pvalue" if "skat_o_pvalue" in columns else "skat_pvalue"
    burden = [t for t in BURDEN_TESTS if t in test_names]
    return f"{burden[0]}_pvalue" if burden else "skat_pvalue"  # else SKAT --skat-method Burden


def is_primary_statistic(stat: str, primary_id: str, skat_method: str) -> bool:
    """Whether ``stat`` is the statistic behind the primary p-value itself."""
    if primary_id in BURDEN_TESTS:
        return stat == "burden"
    return primary_id == "skat" and skat_method in ("SKAT", "Burden")


def gene_rng(seed: int, gene: str) -> np.random.Generator:
    """Generator for one gene from (seed, gene name): independent of worker order."""
    digest = int.from_bytes(hashlib.sha256(gene.encode()).digest()[:8], "little")
    return np.random.default_rng([seed, digest])


def _standardise(r: np.ndarray, mu: np.ndarray) -> np.ndarray:
    """(r) / sqrt(mu (1 - mu)); 0 where mu is 0 or 1 and r is 0 (zero null weight, as in SKAT)."""
    phi = np.sqrt(mu * (1.0 - mu))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where((phi == 0) & (r == 0), 0.0, r / phi)


def empirical_p(
    stat: Callable[[np.ndarray], np.ndarray],
    a: np.ndarray,
    mu: np.ndarray,
    e_obs: np.ndarray,
    rng: np.random.Generator,
    n_perm: int,
    h: int = H_STOP,
) -> tuple[float, int]:
    """Return (empirical p, replicates used); Besag-Clifford h / L, else (r + 1) / (N + 1).

    (nan, 0) when the observed or a replicate statistic is not finite: never the p floor.
    """
    t_obs = float(stat(a.T @ e_obs))
    if not np.isfinite(t_obs):
        return np.nan, 0
    # Ties count as exceedances; the margin absorbs round-off between GEMV and GEMM.
    thresh = t_obs - 1e-9 * abs(t_obs)
    case, control = _standardise(1.0 - mu, mu)[:, None], _standardise(-mu, mu)[:, None]
    mu_col = mu[:, None]
    done = exceed = 0
    batch = FIRST_BATCH
    while done < n_perm:
        batch = min(batch, n_perm - done)
        e = np.where(rng.random((mu.size, batch)) < mu_col, case, control)
        t = stat(a.T @ e)  # the batch's one GEMM
        if not np.isfinite(t).all():
            return np.nan, 0
        cum = exceed + np.cumsum(t >= thresh)
        if cum[-1] >= h:
            used = done + int(np.argmax(cum >= h)) + 1
            return h / used, used
        exceed, done = int(cum[-1]), done + batch
        batch = min(2 * batch, MAX_BATCH)
    return (exceed + 1) / (n_perm + 1), n_perm


def _gene_empirical(
    gene: str, builder: Any, config: Any, stat: str, nulls: dict, q_below: float | None = None
) -> tuple[float, int] | None:
    """Empirical p for one gene from its rebuilt matrix; None when it cannot be tested.

    ``q_below``: skip the gene unless SKAT Q's own analytic p is below it (SKAT-O and
    ``--skat-robust`` runs).
    """
    built = builder()
    geno, y, cov = (
        built["genotype_matrix"],
        built["phenotype_vector"],
        built.get("covariate_matrix"),
    )
    if y is None or geno.shape[1] == 0:
        return None
    key = (y.tobytes(), None if cov is None else cov.tobytes())
    if key not in nulls:
        if len(nulls) >= _MAX_CACHED_NULLS:  # bounded, as artifact.null_for_mask
            del nulls[next(iter(nulls))]
        try:
            nulls[key] = fit_glm_null_model(y, cov, "binary")
        except UntestableGeneError:
            nulls[key] = None
    null = nulls[key]
    if null is None:
        return None
    weights = get_weights(
        built["variant_mafs"],
        config.variant_weights,
        cadd_scores=built.get("cadd_scores"),
        revel_scores=built.get("revel_scores"),
        score_values=built.get("score_values"),
        variant_effects=built.get("variant_effects"),
        weight_params=config.variant_weight_params,
    )
    a = project_out_nuisance(geno * weights[np.newaxis, :], null)
    mu = null.extra["mu_hat"]
    e_obs = _standardise(null.extra["residuals"], mu)
    if q_below is not None:
        lambdas = compute_eigenvalues_filtered(geno * weights[np.newaxis, :], null)
        q = float(_STATISTICS["skat"](a.T @ e_obs)) / 2.0  # SKAT Q = ||Z'r||^2 / 2
        if len(lambdas) == 0 or not compute_pvalue(q, lambdas)[0] < q_below:
            return None
    return empirical_p(
        _STATISTICS[stat], a, mu, e_obs, gene_rng(config.seed, gene), config.permutations
    )


def _gene_range(state: tuple, lo: int, hi: int) -> list[tuple[float, int] | None]:
    genes, builders, config, stat, q_below = state
    nulls: dict = {}
    return [_gene_empirical(g, builders[g], config, stat, nulls, q_below) for g in genes[lo:hi]]


def add_empirical_p(
    result_df: pd.DataFrame,
    builders: dict[str, Any],
    test_names: Sequence[str],
    config: Any,
    primary_id: str,
) -> pd.DataFrame:
    """Add empirical_p, n_permutations and permutation_method for genes below the threshold.

    Genes are selected on the p of the test whose statistic is permuted; a gene that
    test gave no p is not permuted and ``permutation_method`` says why.
    """
    stat = permutation_statistic(test_names, primary_id, config.skat_method, config.trait_type)
    column = source_pvalue_column(stat, test_names, list(result_df.columns))
    source_p = pd.to_numeric(result_df[column], errors="coerce")
    # Every SKAT mode selects on SKAT Q's standard analytic p: the reported skat_pvalue for
    # default SKAT, else computed in the worker for every gene with a p (a SKAT-O or robust
    # p is a different statistic's p, so it cannot pre-filter).
    robust = getattr(config, "skat_robust", False)
    own_q = stat == "skat" and (column == "skat_o_pvalue" or robust)
    q_below = config.permute_below_p if own_q else None
    below = source_p.notna() if q_below else source_p < config.permute_below_p
    genes = [g for g in result_df.loc[below, "gene"] if builders.get(g) is not None]
    reason = result_df.get("skat_skip_reason") if stat == "skat" else None
    skipped = [
        f"skipped_{r}" if isinstance(r, str) and r else "skipped_no_pvalue"
        for r in (reason if reason is not None else [None] * len(result_df))
    ]
    if not is_primary_statistic(stat, primary_id, config.skat_method):
        logger.warning(
            f"Permutation: primary result '{primary_id}' cannot be rebuilt from score "
            f"matrices; empirical_p is for the {METHOD_LABEL[stat]} statistic instead"
        )
    state = (genes, builders, config, stat, q_below)
    workers = config.association_workers
    workers = min((os.cpu_count() or 1) if workers == -1 else workers, len(genes))
    t0 = time.perf_counter()
    with single_blas_thread():  # identical floating point for any number of workers
        if workers > 1 and fork_available():
            chunks = fork_map(_gene_range, state, len(genes), workers)
            found = [r for _lo, _hi, rs in chunks for r in rs]
        else:
            found = _gene_range(state, 0, len(genes))
    by_gene = {g: r for g, r in zip(genes, found, strict=True) if r is not None}
    elapsed = time.perf_counter() - t0
    logger.info(
        f"Permutation: {len(by_gene)} of {len(genes)} genes with {column} < "
        f"{config.permute_below_p:g} got an empirical p ({METHOD_LABEL[stat]}, N <= "
        f"{config.permutations}, seed {config.seed}) in {elapsed:.1f}s"
    )
    hits = result_df["gene"].map(by_gene)
    return result_df.assign(
        empirical_p=[r[0] if isinstance(r, tuple) else np.nan for r in hits],
        n_permutations=pd.array(
            [r[1] if isinstance(r, tuple) and r[1] else None for r in hits], dtype="Int64"
        ),
        permutation_method=[
            ("skipped_nonfinite" if np.isnan(r[0]) else METHOD_LABEL[stat])
            if isinstance(r, tuple)
            else (skip if np.isnan(p) else None)
            for r, p, skip in zip(hits, source_p, skipped, strict=True)
        ],
    )
