# File: variantcentrifuge/association/backends/__init__.py
# Location: variantcentrifuge/variantcentrifuge/association/backends/__init__.py
"""
SKAT backend factory and public API for the backends subpackage.

Usage
-----
>>> from variantcentrifuge.association.backends import get_skat_backend
>>> backend = get_skat_backend("python")
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
        - ``"python"`` — PythonSKATBackend (pure Python via numpy/scipy/statsmodels). Phase 21.
        - ``"auto"``   — Try ``"r"`` first; fall back to ``"python"`` if R is
                         unavailable. The Python backend is always available.

    Returns
    -------
    SKATBackend
        Uninitialised backend instance (detect_environment() not yet called).

    Raises
    ------
    ValueError
        If ``backend_name`` is not one of the recognised values.
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
        from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

        return PythonSKATBackend()

    if backend_name == "auto":
        # Try R first; fall through to Python on ImportError / RuntimeError
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

        # R backend unavailable — fall back to Python backend (always available)
        from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

        return PythonSKATBackend()

    raise ValueError(
        f"Unknown SKAT backend: '{backend_name}'. Valid values: 'r', 'python', 'auto'."
    )
