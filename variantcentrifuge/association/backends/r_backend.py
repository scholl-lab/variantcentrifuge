# File: variantcentrifuge/association/backends/r_backend.py
# Location: variantcentrifuge/variantcentrifuge/association/backends/r_backend.py
"""
R SKAT backend via rpy2.

RSKATBackend wraps the R SKAT package through rpy2's robjects interface.
All rpy2 imports are deferred to method bodies — never at module level —
because ``import rpy2.robjects`` initializes the R runtime immediately and
would break environments where R is not installed.

Thread safety
-------------
rpy2 is NOT thread-safe. Calling rpy2 functions from any thread other than
the main thread can cause segfaults with no Python traceback. RSKATBackend
enforces main-thread-only execution via _assert_main_thread(). The stage
that instantiates this backend must set parallel_safe=False.

Version requirements (warnings emitted but not enforced)
---------------------------------------------------------
R >= 4.0        — required for SKAT 2.x compatibility
SKAT >= 2.0     — SKATBinary moment adjustment improvements
rpy2 >= 3.5.0   — stable numpy <-> R vector conversion
"""

from __future__ import annotations

import gc
import logging
import os
import threading
import time
from typing import Any

import numpy as np

from variantcentrifuge.association.backends.base import NullModelResult, SKATBackend

logger = logging.getLogger("variantcentrifuge")

# Minimum recommended versions (warnings only, not hard requirements)
_MIN_R_VERSION = "4.0"
_MIN_SKAT_VERSION = "2.0"
_MIN_RPY2_VERSION = "3.5.0"


def _version_lt(version_str: str, minimum: str) -> bool:
    """
    Return True if version_str is strictly less than minimum.

    Compares dot-separated integer tuples. Non-integer components are ignored
    (treated as 0) so "4.0.3" vs "4.0" works correctly.

    Parameters
    ----------
    version_str : str
        Detected version (e.g. "3.6.5").
    minimum : str
        Minimum version to compare against (e.g. "3.5.0").

    Returns
    -------
    bool
    """
    try:
        # Attempt packaging.version first (available in most venvs)
        from packaging.version import Version

        return Version(version_str) < Version(minimum)
    except Exception:
        pass

    # Fallback: integer tuple comparison (handles "4.0" vs "4.0.3")
    def to_tuple(v: str) -> tuple[int, ...]:
        parts = []
        for p in v.split("."):
            try:
                parts.append(int(p))
            except ValueError:
                parts.append(0)
        return tuple(parts)

    v_tuple = to_tuple(version_str)
    m_tuple = to_tuple(minimum)
    # Pad shorter tuple with zeros
    length = max(len(v_tuple), len(m_tuple))
    v_padded = v_tuple + (0,) * (length - len(v_tuple))
    m_padded = m_tuple + (0,) * (length - len(m_tuple))
    return v_padded < m_padded


class RSKATBackend(SKATBackend):
    """
    SKAT backend using the R SKAT package via rpy2.

    Lifecycle
    ---------
    1. Instantiate RSKATBackend()
    2. Call detect_environment() — raises ImportError with diagnostics if R or
       SKAT package is unavailable
    3. Call log_environment() — emits INFO version line + version warnings
    4. Call fit_null_model() — fit once per cohort
    5. Call test_gene() for each gene
    6. Call cleanup() to release R memory

    All rpy2 calls must occur on the main thread (enforced by _assert_main_thread).

    Constants
    ---------
    GC_INTERVAL : int
        Number of genes between R gc() calls (100 — fixed, not configurable).
    R_HEAP_WARNING_GB : float
        Threshold in GB for R heap size WARNING log (4.0 GB).
    LARGE_PANEL_THRESHOLD : int
        Gene count above which a WARNING is emitted before the gene loop (2000).
    """

    # Memory management constants
    GC_INTERVAL: int = 100
    R_HEAP_WARNING_GB: float = 4.0
    LARGE_PANEL_THRESHOLD: int = 2000

    def __init__(self) -> None:
        # Loaded R package handles (set by detect_environment)
        self._skat_pkg: Any = None
        self._base_pkg: Any = None

        # Cached R gc() function (set by detect_environment, used by _run_r_gc)
        self._r_gc_func: Any = None

        # Version strings (set by detect_environment)
        self._rpy2_version: str = "<unknown>"
        self._r_version: str = "<unknown>"
        self._skat_version: str = "<unknown>"
        self._r_home: str = "<not set>"

        # Tracking for cleanup() summary
        self._genes_processed: int = 0
        self._start_time: float = time.time()

    def _assert_main_thread(self) -> None:
        """
        Raise RuntimeError if called from a non-main thread.

        rpy2 is not thread-safe. Calling R from a worker thread created by
        ThreadPoolExecutor causes segfaults with no Python traceback. This
        guard makes the failure explicit and debuggable.

        Raises
        ------
        RuntimeError
            If the current thread is not the main thread.
        """
        if threading.current_thread() is not threading.main_thread():
            raise RuntimeError(
                "RSKATBackend called from a non-main thread. rpy2/R is not thread-safe "
                "and must only be used from the main thread. Ensure the stage using "
                "RSKATBackend declares parallel_safe=False."
            )

    def detect_environment(self) -> None:
        """
        Verify that rpy2 is importable and the R SKAT package is installed.

        Detection steps
        ---------------
        1. Import rpy2.robjects (deferred — initializes R runtime on first call)
        2. Check that SKAT package is installed in the R library
        3. Load SKAT and base R packages; cache as instance attributes
        4. Extract and store rpy2, R, and SKAT version strings

        Raises
        ------
        ImportError
            If rpy2 cannot be imported (R not installed or R_HOME misconfigured)
            or if the SKAT package is not installed in R.

        Notes
        -----
        All rpy2 calls are deferred to this method body (never at module level)
        to ensure that simply importing this module does not initialize R.
        """
        self._assert_main_thread()

        r_home = os.environ.get("R_HOME", "<not set>")
        self._r_home = r_home

        # Step 1: Try importing rpy2.robjects
        try:
            import rpy2.robjects as ro
            import rpy2.robjects.packages as rpacks

            _ = (ro, rpacks)  # probe imports; used later in this function
        except Exception as exc:
            raise ImportError(
                f"rpy2 import failed — R backend is not available.\n"
                f"  R_HOME: {r_home}\n"
                f"  PATH: {os.environ.get('PATH', '<not set>')}\n"
                f"\n"
                f"To fix:\n"
                f"  1. Install R from https://www.r-project.org/ (>= 4.0 recommended)\n"
                f"  2. Set R_HOME to point to your R installation, e.g.:\n"
                f"       export R_HOME=$(R RHOME)\n"
                f"  3. Install rpy2: pip install 'rpy2>=3.5.0'\n"
                f"\n"
                f"Original error: {type(exc).__name__}: {exc}"
            ) from exc

        # Step 2: Check SKAT package is installed in R
        import rpy2.robjects.packages as rpacks

        if not rpacks.isinstalled("SKAT"):
            raise ImportError(
                "R found but SKAT package is not installed.\n"
                "Install with: install.packages('SKAT') in R\n"
                "\n"
                "Or from the command line:\n"
                "  Rscript -e \"install.packages('SKAT', repos='https://cloud.r-project.org')\""
            )

        # Step 3: Load packages and cache handles
        from rpy2.robjects.packages import importr

        self._skat_pkg = importr("SKAT")
        self._base_pkg = importr("base")

        # Step 4: Extract version strings
        import rpy2

        self._rpy2_version = str(rpy2.__version__)

        import rpy2.robjects as ro

        r_ver_list = ro.r("R.Version()")
        # R.Version() returns a named list; extract "major" and "minor"
        try:
            r_ver_names = list(r_ver_list.names)
            major_idx = r_ver_names.index("major")
            minor_idx = r_ver_names.index("minor")
            r_major = str(r_ver_list[major_idx][0])
            r_minor = str(r_ver_list[minor_idx][0])
            self._r_version = f"{r_major}.{r_minor}"
        except Exception:
            try:
                self._r_version = (
                    str(ro.r("R.version$major")[0]) + "." + str(ro.r("R.version$minor")[0])
                )
            except Exception:
                self._r_version = "<unknown>"

        # Extract SKAT version from packageVersion()
        try:
            skat_ver = ro.r("as.character(packageVersion('SKAT'))")
            self._skat_version = str(skat_ver[0])
        except Exception:
            self._skat_version = "<unknown>"

        # Cache R gc() function for efficient memory management
        self._r_gc_func = ro.r["gc"]

    def log_environment(self) -> None:
        """
        Log R environment info at INFO level; warn if versions below recommended.

        Emits one INFO line summarizing all versions, then version-specific
        WARNING lines if any component is below the recommended minimum.

        Should be called immediately after detect_environment() succeeds.
        """
        logger.info(
            f"R SKAT backend: rpy2={self._rpy2_version}, "
            f"R={self._r_version}, "
            f"SKAT={self._skat_version}, "
            f"R_HOME={self._r_home}"
        )

        if self._r_version != "<unknown>" and _version_lt(self._r_version, _MIN_R_VERSION):
            logger.warning(
                f"R version {self._r_version} may not support all SKAT 2.x features. "
                f"Recommend R >= {_MIN_R_VERSION}."
            )

        if self._skat_version != "<unknown>" and _version_lt(self._skat_version, _MIN_SKAT_VERSION):
            logger.warning(
                f"SKAT version {self._skat_version} predates SKATBinary improvements. "
                f"Recommend SKAT >= {_MIN_SKAT_VERSION}."
            )

        if self._rpy2_version != "<unknown>" and _version_lt(self._rpy2_version, _MIN_RPY2_VERSION):
            logger.warning(
                f"rpy2 version {self._rpy2_version} may have unstable numpy <-> R "
                f"vector conversion. Recommend rpy2 >= {_MIN_RPY2_VERSION}."
            )

    def fit_null_model(
        self,
        phenotype: np.ndarray,
        covariates: np.ndarray | None,
        trait_type: str,
    ) -> NullModelResult:
        """
        Fit the SKAT null model using SKAT_Null_Model with Adjustment=TRUE.

        The null model is fit once per cohort (all samples) and reused for
        every gene. Adjustment=TRUE enables moment adjustment for small samples
        (n < 2000), using the unified SKAT 2.x API.

        Parameters
        ----------
        phenotype : np.ndarray, shape (n_samples,)
            Binary (0/1) or continuous phenotype vector. Order must match
            the genotype matrix rows and covariate matrix rows.
        covariates : np.ndarray or None, shape (n_samples, k)
            Covariate matrix or None for intercept-only model.
        trait_type : str
            "binary" -> out_type="D" (dichotomous/SKATBinary)
            "quantitative" -> out_type="C" (continuous/SKAT)

        Returns
        -------
        NullModelResult
            Wraps the R null model object with metadata.

        Raises
        ------
        RuntimeError
            If called from a non-main thread.
        """
        self._assert_main_thread()

        import rpy2.robjects as ro

        # Map Python trait_type to R SKAT out_type parameter
        out_type = "D" if trait_type == "binary" else "C"

        # Build phenotype vector in R global env
        ro.globalenv["._vc_y"] = ro.FloatVector(phenotype.tolist())

        # Build formula string and optionally set covariate matrix
        if covariates is not None and covariates.shape[1] > 0:
            r_cov = ro.r["matrix"](
                ro.FloatVector(covariates.ravel(order="F").tolist()),
                nrow=len(phenotype),
                ncol=covariates.shape[1],
            )
            ro.globalenv["._vc_x"] = r_cov
            formula_str = "._vc_y ~ ._vc_x"
        else:
            formula_str = "._vc_y ~ 1"

        # Fit null model — Adjustment=TRUE enables moment adjustment for small samples
        # Uses the unified SKAT 2.x API (not deprecated SKAT_Null_Model_MomentAdjust)
        null_obj = self._skat_pkg.SKAT_Null_Model(
            ro.Formula(formula_str),
            out_type=out_type,
            Adjustment=True,
        )

        # Clean up R global env (but keep the null model object — it's returned)
        ro.r("rm(list=ls(pattern='\\._vc_'))")

        logger.debug(
            f"Null model fit: trait_type={trait_type}, out_type={out_type}, "
            f"n_samples={len(phenotype)}, adjustment=True"
        )

        # Reset tracking for this cohort
        self._genes_processed = 0
        self._start_time = time.time()

        return NullModelResult(
            model=null_obj,
            trait_type=trait_type,
            n_samples=len(phenotype),
            adjustment=True,
        )

    def test_gene(
        self,
        gene: str,
        genotype_matrix: np.ndarray,
        null_model: NullModelResult,
        method: str,
        weights_beta: tuple[float, float],
    ) -> dict[str, Any]:
        """
        Run SKAT for a single gene.

        Dispatches to SKATBinary (binary phenotype) or SKAT (continuous).
        SKAT-O is selected via method="SKATO" and returns the optimal rho.

        R warnings are captured via withCallingHandlers and returned in the
        result dict. Statistical NA (R NA_Real) is mapped to p_value=None.

        Parameters
        ----------
        gene : str
            Gene symbol (used only for logging).
        genotype_matrix : np.ndarray, shape (n_samples, n_variants)
            Dosage matrix (0/1/2). Must use Fortran/column-major order when
            passed to R (handled internally).
        null_model : NullModelResult
            Fitted null model from fit_null_model().
        method : str
            SKAT variant: "SKAT" (default), "Burden", or "SKATO" (omnibus).
        weights_beta : tuple of (float, float)
            Beta distribution parameters (a1, a2) for variant weighting.
            Use (1.0, 25.0) for SKAT default, (1.0, 1.0) for uniform weights.

        Returns
        -------
        dict with keys:
            p_value : float | None
                P-value from SKAT. None if R returns statistical NA.
            rho : float | None
                Optimal rho from SKAT-O (method="SKATO"), else None.
            n_variants : int
                Number of variants in the genotype matrix.
            n_marker_test : int | None
                Number of markers actually used in the test (after internal QC).
            warnings : list of str
                R warnings captured via withCallingHandlers.

        Raises
        ------
        RuntimeError
            If called from a non-main thread.
        Exception
            Any R-level exception propagates directly (infrastructure failure).
        """
        self._assert_main_thread()

        import rpy2.robjects as ro
        from rpy2.rinterface import NA_Real

        n_samples, n_variants = genotype_matrix.shape
        a1, a2 = weights_beta

        # Convert numpy matrix to R matrix (column-major order for R)
        r_z = ro.r["matrix"](
            ro.FloatVector(genotype_matrix.ravel(order="F").tolist()),
            nrow=n_samples,
            ncol=n_variants,
        )

        # Assign to R global env
        ro.globalenv["._vc_Z"] = r_z
        ro.globalenv["._vc_obj"] = null_model.model

        # Dispatch: binary -> SKATBinary, continuous -> SKAT
        func_name = "SKATBinary" if null_model.trait_type == "binary" else "SKAT"

        # Build R code with withCallingHandlers for warning capture
        # Weights vector uses the Beta distribution parameters
        r_code = (
            "local({"
            "  warns <- character(0);"
            f"  result <- withCallingHandlers("
            f"    {func_name}(._vc_Z, ._vc_obj, kernel='linear.weighted',"
            f"                method='{method}', weights.beta=c({a1},{a2})),"
            "    warning = function(w) {"
            "      warns <<- c(warns, conditionMessage(w));"
            "      invokeRestart('muffleWarning')"
            "    }"
            "  );"
            "  list(result=result, warnings=warns)"
            "})"
        )

        # Execute — R exceptions propagate (abort on infrastructure failure)
        outer = ro.r(r_code)

        # Extract inner result and warnings
        result = outer.rx2("result")
        r_warnings_vec = outer.rx2("warnings")
        warnings_list = list(r_warnings_vec) if len(r_warnings_vec) > 0 else []

        # Extract p-value with NA_Real guard
        p_val_r = result.rx2("p.value")[0]
        p_value = None if p_val_r is NA_Real else float(p_val_r)

        # Extract SKAT-O rho (only present when method="SKATO")
        rho: float | None = None
        if method == "SKATO":
            try:
                rho_r = result.rx2("param").rx2("rho")[0]
                rho = None if rho_r is NA_Real else float(rho_r)
            except Exception:
                rho = None

        # Extract n_marker_test from param (number of markers actually tested)
        n_marker_test: int | None = None
        try:
            nmt_r = result.rx2("param").rx2("n.marker.test")[0]
            n_marker_test = None if nmt_r is NA_Real else int(nmt_r)
        except Exception:
            n_marker_test = None

        # Clean up R globals
        ro.r("rm(list=ls(pattern='\\._vc_'))")

        # Release rpy2 wrapper for the R matrix (frees Python reference)
        del r_z

        self._genes_processed += 1

        logger.debug(
            f"Gene {gene}: {func_name}(method={method}) -> p={p_value}, "
            f"n_variants={n_variants}, n_marker_test={n_marker_test}, "
            f"warnings={len(warnings_list)}"
        )

        return {
            "p_value": p_value,
            "rho": rho,
            "n_variants": n_variants,
            "n_marker_test": n_marker_test,
            "warnings": warnings_list,
        }

    def _run_r_gc(self) -> None:
        """
        Run R garbage collection and Python gc.

        Calls R's gc() via the cached function handle and then Python's
        gc.collect(). Logs at DEBUG level.

        Called every GC_INTERVAL genes to prevent R heap accumulation.
        """
        try:
            if self._r_gc_func is not None:
                self._r_gc_func()
            else:
                # Fallback if gc func not cached yet
                import rpy2.robjects as ro

                ro.r["gc"]()
        except Exception as exc:
            logger.debug(f"RSKATBackend._run_r_gc: R gc() failed: {exc}")
        gc.collect()
        logger.debug(f"RSKATBackend: R GC triggered after {self._genes_processed} genes")

    def _check_r_heap(self) -> float | None:
        """
        Check current R heap usage in MB.

        Calls ``gc(verbose=FALSE)`` in R and extracts Vcells used (row 2,
        col 1 of the returned matrix). Vcells are 8-byte double-precision
        cells, so heap_mb = vcells_used * 8 / 1024^2.

        Returns None on any exception (non-fatal; monitoring only).
        Logs WARNING if heap exceeds R_HEAP_WARNING_GB threshold.

        Returns
        -------
        float | None
            Heap usage in MB, or None if extraction failed.
        """
        try:
            import rpy2.robjects as ro

            gc_result = ro.r("gc(verbose=FALSE)")
            # gc() returns a matrix with rows: Ncells, Vcells
            # columns: used, gc trigger, max used, limit
            # Vcells (row index 1, 0-based) used (col index 0)
            vcells_used = gc_result.rx(2, 1)[0]  # row 2 (1-based), col 1 (1-based)
            heap_mb = float(vcells_used) * 8.0 / (1024.0 * 1024.0)
            heap_gb = heap_mb / 1024.0
            if heap_gb > self.R_HEAP_WARNING_GB:
                logger.warning(
                    f"R heap usage {heap_gb:.1f} GB exceeds "
                    f"{self.R_HEAP_WARNING_GB} GB threshold — "
                    "consider reducing panel size or increasing available memory"
                )
            return heap_mb
        except Exception as exc:
            logger.debug(f"RSKATBackend._check_r_heap: failed to check heap: {exc}")
            return None

    def cleanup(self) -> None:
        """
        Release R-side memory and trigger Python garbage collection.

        Calls R's gc() to release any rpy2-held R objects and then runs
        Python's gc.collect() to free Python-side rpy2 wrappers.

        Logs a summary of total genes processed and elapsed time.

        Safe to call even if detect_environment() was never called (no-op if
        no rpy2 session is active).
        """
        elapsed = time.time() - self._start_time
        logger.info(
            f"RSKATBackend.cleanup: {self._genes_processed} genes processed "
            f"in {elapsed:.1f}s"
        )
        self._run_r_gc()
