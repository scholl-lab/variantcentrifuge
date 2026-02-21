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
    if name == "LogisticBurdenTest":
        from variantcentrifuge.association.tests.logistic_burden import LogisticBurdenTest

        return LogisticBurdenTest
    if name == "LinearBurdenTest":
        from variantcentrifuge.association.tests.linear_burden import LinearBurdenTest

        return LinearBurdenTest
    if name == "RSKATTest":
        from variantcentrifuge.association.tests.skat_r import RSKATTest

        return RSKATTest
    if name == "PurePythonSKATTest":
        from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest

        return PurePythonSKATTest
    if name == "COASTTest":
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        return COASTTest
    if name in ("cauchy_combination", "compute_acat_o"):
        from variantcentrifuge.association.tests.acat import cauchy_combination, compute_acat_o

        if name == "cauchy_combination":
            return cauchy_combination
        return compute_acat_o
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "COASTTest",
    "FisherExactTest",
    "LinearBurdenTest",
    "LogisticBurdenTest",
    "PurePythonSKATTest",
    "RSKATTest",
    "cauchy_combination",
    "compute_acat_o",
]
