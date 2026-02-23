# File: variantcentrifuge/association/tests/_utils.py
# Location: variantcentrifuge/variantcentrifuge/association/tests/_utils.py
"""
Shared utilities for association test implementations.

This module contains helpers shared across multiple test classes to avoid
code duplication and prevent transitive dependency issues (e.g., importing
_parse_weights_beta from skat_r.py would transitively import rpy2).
"""

from __future__ import annotations

import logging

logger = logging.getLogger("variantcentrifuge")


def parse_weights_beta(variant_weights: str) -> tuple[float, float]:
    """
    Parse weight specification string to Beta distribution parameters.

    Handles ``"beta:a,b"`` format and ``"uniform"`` special case. Falls back
    to (1.0, 25.0) — the SKAT paper default — if the string cannot be parsed.

    Parameters
    ----------
    variant_weights : str
        Weight scheme, e.g. ``"beta:1,25"``, ``"uniform"``, or any string.

    Returns
    -------
    tuple of (float, float)
        Beta distribution parameters (a1, a2).

    Examples
    --------
    >>> parse_weights_beta("beta:1,25")
    (1.0, 25.0)
    >>> parse_weights_beta("uniform")
    (1.0, 1.0)
    >>> parse_weights_beta("beta:0.5,0.5")
    (0.5, 0.5)
    """
    if variant_weights == "uniform":
        return 1.0, 1.0
    if variant_weights.startswith("beta:"):
        try:
            params = variant_weights[5:].split(",")
            return float(params[0]), float(params[1])
        except (IndexError, ValueError):
            logger.warning(
                f"Could not parse variant_weights '{variant_weights}'. "
                "Using SKAT default Beta(1, 25)."
            )
    return 1.0, 25.0
