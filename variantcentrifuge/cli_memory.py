"""Memory calculation helpers for VariantCentrifuge CLI."""

from __future__ import annotations

import logging

logger = logging.getLogger("variantcentrifuge")


def calculate_sort_memory_limit(sort_limit: str, max_memory_gb: float | None = None) -> str:
    """Calculate intelligent sort memory limit based on available memory."""
    if sort_limit != "auto":
        return sort_limit

    try:
        import psutil

        available_gb = (
            max_memory_gb if max_memory_gb else psutil.virtual_memory().available / (1024**3)
        )
        sort_memory_gb = max(1, min(available_gb * 0.10, 50))
        return (
            f"{int(sort_memory_gb)}G" if sort_memory_gb >= 1 else f"{int(sort_memory_gb * 1024)}M"
        )
    except Exception as e:
        logger.debug(f"Could not calculate sort memory limit: {e}. Using default 4G")
        return "4G"
