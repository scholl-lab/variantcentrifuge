# File: variantcentrifuge/association/tests/__init__.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/__init__.py
"""
Association test implementations.

Tests are loaded lazily so that optional dependencies (scipy, statsmodels)
are only imported when a specific test is actually instantiated.
"""

from __future__ import annotations


def __getattr__(name: str) -> object:
    if name == "FisherExactTest":
        from variantcentrifuge.association.tests.fisher import FisherExactTest

        return FisherExactTest
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["FisherExactTest"]
