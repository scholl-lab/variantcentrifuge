# File: variantcentrifuge/association/backends/base.py
# Location: variantcentrifuge/variantcentrifuge/association/backends/base.py
"""
Abstract base class and data containers for SKAT backend implementations.

The SKATBackend ABC defines the interface all backend implementations must
satisfy. NullModelResult is a lightweight container for the fitted null model
that is passed between fit_null_model() and test_gene() calls.

Current implementations
-----------------------
RSKATBackend (r_backend.py) — R SKAT package via rpy2 (Phase 20)

Planned
-------
PurePythonSKATBackend — Davies ctypes + saddlepoint + Liu fallback (Phase 21)
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

import numpy as np

logger = logging.getLogger("variantcentrifuge")


@dataclass
class NullModelResult:
    """
    Container for a fitted SKAT null model.

    Holds the backend-specific null model object (an R object for RSKATBackend)
    alongside metadata needed to validate compatibility when the model is reused
    across genes.

    Fields
    ------
    model : Any
        The fitted null model object. For RSKATBackend this is an rpy2
        ListVector wrapping the R SKAT_Null_Model result.
    trait_type : str
        "binary" or "quantitative". Determines whether SKATBinary or the
        continuous-trait SKAT function was used during fitting.
    n_samples : int
        Number of samples used to fit the null model. Used to validate that
        per-gene genotype matrices have the correct number of rows.
    adjustment : bool
        Whether small-sample moment adjustment was applied during fitting
        (only relevant for SKATBinary with n < 2000).
    """

    model: Any
    trait_type: str
    n_samples: int
    adjustment: bool = True
    extra: dict[str, Any] = field(default_factory=dict)


class SKATBackend(ABC):
    """
    Abstract interface for SKAT computation backends.

    All implementations must support the three-step lifecycle:

    1. detect_environment() — verify the backend is available (eager check)
    2. fit_null_model() — fit once before the gene loop
    3. test_gene() — called once per gene, reusing the null model
    4. cleanup() — release R memory after gene loop completes

    Thread safety
    -------------
    Backend implementations that are not thread-safe (e.g. RSKATBackend due
    to rpy2) must raise RuntimeError from a non-main-thread context rather than
    silently producing incorrect results or causing a segfault.
    """

    @abstractmethod
    def detect_environment(self) -> None:
        """
        Verify that this backend's runtime dependencies are available.

        Called eagerly at engine construction (before any data processing)
        so users receive clear diagnostics early.

        Raises
        ------
        ImportError
            If a required library (e.g. rpy2) or R package (e.g. SKAT) is not
            installed. The error message must include actionable install
            instructions.
        RuntimeError
            If the environment is detected but misconfigured (e.g. R_HOME
            points to a non-existent path).
        """
        ...

    @abstractmethod
    def log_environment(self) -> None:
        """
        Log version information for the backend and its dependencies.

        Called after detect_environment() succeeds. Emits INFO-level messages
        with library versions and warns if versions are below recommended
        minimums.
        """
        ...

    @abstractmethod
    def fit_null_model(
        self,
        phenotype: np.ndarray,
        covariates: np.ndarray | None,
        trait_type: str,
    ) -> NullModelResult:
        """
        Fit the SKAT null model on the cohort-level phenotype + covariates.

        The null model is fit once and reused for all genes (per SKAT paper
        design: the null model does not depend on the variant data).

        Parameters
        ----------
        phenotype : np.ndarray, shape (n_samples,)
            Outcome vector. Binary: 0/1 (control/case). Quantitative: float.
        covariates : np.ndarray, shape (n_samples, k), or None
            Covariate matrix (no intercept column; backend adds intercept).
            None = intercept only.
        trait_type : str
            "binary" → SKATBinary (moment-adjusted). "quantitative" → SKAT
            with linear kernel.

        Returns
        -------
        NullModelResult
            Fitted null model container.

        Raises
        ------
        NotImplementedError
            In skeleton implementations (Plan 20-01) until Plan 20-02
            completes the backend logic.
        RuntimeError
            If called from a non-main thread in thread-unsafe backends.
        """
        ...

    @abstractmethod
    def test_gene(
        self,
        gene: str,
        genotype_matrix: np.ndarray,
        null_model: NullModelResult,
        method: str,
        weights_beta: tuple[float, float],
    ) -> dict[str, Any]:
        """
        Run the SKAT test for a single gene.

        Parameters
        ----------
        gene : str
            Gene symbol (used only for logging).
        genotype_matrix : np.ndarray, shape (n_samples, n_variants)
            Imputed dosage matrix (0/1/2, no NaN). Rows = samples in the
            same order as phenotype passed to fit_null_model().
        null_model : NullModelResult
            Fitted null model from fit_null_model().
        method : str
            SKAT variant: "SKAT" (default), "Burden", or "SKATO".
        weights_beta : tuple of (float, float)
            Beta distribution parameters (a1, a2) for variant weights.
            SKAT convention: (1, 25).

        Returns
        -------
        dict with keys:
            p_value : float | None
                Raw p-value. None if test skipped (e.g. no variants).
            rho : float | None
                Optimal rho from SKAT-O (None for SKAT/Burden).
            n_variants : int
                Number of variants in the genotype matrix.
            n_marker_test : int | None
                Number of variants actually tested (after any backend filtering).
            warnings : list[str]
                Backend-generated warning messages.

        Raises
        ------
        NotImplementedError
            In skeleton implementations until Plan 20-02 is complete.
        RuntimeError
            If called from a non-main thread in thread-unsafe backends.
        """
        ...

    @abstractmethod
    def cleanup(self) -> None:
        """
        Release backend resources after the gene loop completes.

        For RSKATBackend: triggers R garbage collection and Python gc.collect().
        Should be idempotent (safe to call even if no test was run).
        """
        ...
