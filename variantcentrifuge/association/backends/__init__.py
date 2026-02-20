# File: variantcentrifuge/association/backends/__init__.py
# Location: variantcentrifuge/variantcentrifuge/association/backends/__init__.py
"""
SKAT backend factory and public API for the backends subpackage.

Usage
-----
>>> from variantcentrifuge.association.backends import get_skat_backend
>>> backend = get_skat_backend("r")
>>> backend.detect_environment()
>>> backend.log_environment()

Importing this module does NOT initialize R or import rpy2 — all backend
imports are deferred to get_skat_backend() so that environments without R
can still import the association package without errors.
"""

from __future__ import annotations

from variantcentrifuge.association.backends.base import NullModelResult, SKATBackend

__all__ = ["NullModelResult", "SKATBackend", "get_skat_backend"]


def get_skat_backend(backend_name: str) -> SKATBackend:
    """
    Factory function for SKAT backend instances.

    Returns an uninitialised SKATBackend instance. Callers must call
    ``detect_environment()`` before using the backend to verify that all
    required runtime dependencies are available.

    Parameters
    ----------
    backend_name : str
        Backend selector:
        - ``"r"``      — RSKATBackend (R SKAT package via rpy2). Phase 20.
        - ``"python"`` — Pure Python SKAT backend. Phase 21 (not yet implemented).
        - ``"auto"``   — Try ``"r"`` first; fall back to ``"python"`` if R is
                         unavailable. Raises ImportError if neither is available.

    Returns
    -------
    SKATBackend
        Uninitialised backend instance (detect_environment() not yet called).

    Raises
    ------
    ValueError
        If ``backend_name`` is not one of the recognised values.
    NotImplementedError
        If ``backend_name`` is ``"python"`` (Phase 21) or if ``"auto"`` mode
        finds no available backend.
    ImportError
        Propagated from RSKATBackend.detect_environment() when the R backend
        is explicitly requested but R or SKAT is not installed. In ``"auto"``
        mode, ImportError from the R backend triggers fall-through to Python.

    Notes
    -----
    All backend class imports are deferred to this function body so that
    importing the backends subpackage never initialises R.
    """
    if backend_name == "r":
        # Deferred import: importing RSKATBackend does not initialise R,
        # but keep it lazy here for symmetry and to avoid any side effects.
        from variantcentrifuge.association.backends.r_backend import RSKATBackend

        return RSKATBackend()

    if backend_name == "python":
        raise NotImplementedError(
            "Pure Python SKAT backend is not yet implemented. "
            "It is planned for Phase 21 of the v0.15.0 milestone. "
            "Use backend_name='r' for R-based SKAT."
        )

    if backend_name == "auto":
        # Try R first; fall through to Python on ImportError
        try:
            from variantcentrifuge.association.backends.r_backend import RSKATBackend

            backend = RSKATBackend()
            # Probe availability without caching; caller will call detect_environment()
            # again in proper context (main thread, after logging setup).
            # We try a lightweight check here: just import rpy2 to see if it's present.
            try:
                import rpy2.robjects

                _ = rpy2.robjects  # probe only; caller will call detect_environment()
                return backend
            except Exception:
                pass
        except Exception:
            pass

        # Fall back to Python backend (will raise NotImplementedError until Phase 21)
        raise NotImplementedError(
            "auto backend: R backend (rpy2) is not available, and the pure Python "
            "SKAT backend is not yet implemented (Phase 21). "
            "To use SKAT, install R and rpy2:\n"
            "  1. Install R >= 4.0: https://www.r-project.org/\n"
            "  2. Install SKAT in R: install.packages('SKAT')\n"
            "  3. Install rpy2: pip install 'rpy2>=3.5.0'"
        )

    raise ValueError(
        f"Unknown SKAT backend: '{backend_name}'. Valid values: 'r', 'python', 'auto'."
    )
