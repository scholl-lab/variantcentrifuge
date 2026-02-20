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
    """

    def __init__(self) -> None:
        # Loaded R package handles (set by detect_environment)
        self._skat_pkg: Any = None
        self._base_pkg: Any = None

        # Version strings (set by detect_environment)
        self._rpy2_version: str = "<unknown>"
        self._r_version: str = "<unknown>"
        self._skat_version: str = "<unknown>"
        self._r_home: str = "<not set>"

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
        Fit the SKAT null model.

        Stub implementation — raises NotImplementedError.
        Full implementation added in Plan 20-02.

        Raises
        ------
        NotImplementedError
            Always (skeleton for Plan 20-01).
        RuntimeError
            If called from a non-main thread.
        """
        self._assert_main_thread()
        raise NotImplementedError(
            "RSKATBackend.fit_null_model() is implemented in Plan 20-02. "
            "This skeleton (Plan 20-01) establishes the interface only."
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

        Stub implementation — raises NotImplementedError.
        Full implementation added in Plan 20-02.

        Raises
        ------
        NotImplementedError
            Always (skeleton for Plan 20-01).
        RuntimeError
            If called from a non-main thread.
        """
        self._assert_main_thread()
        raise NotImplementedError(
            "RSKATBackend.test_gene() is implemented in Plan 20-02. "
            "This skeleton (Plan 20-01) establishes the interface only."
        )

    def cleanup(self) -> None:
        """
        Release R-side memory and trigger Python garbage collection.

        Calls R's gc() to release any rpy2-held R objects and then runs
        Python's gc.collect() to free Python-side rpy2 wrappers.

        Safe to call even if detect_environment() was never called (no-op if
        no rpy2 session is active).
        """
        try:
            import rpy2.robjects as ro

            ro.r["gc"]()
        except Exception as exc:
            logger.debug(f"RSKATBackend.cleanup: R gc() skipped: {exc}")
        gc.collect()
