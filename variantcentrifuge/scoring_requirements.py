"""Scoring requirement analysis helper for VariantCentrifuge pipeline."""

from __future__ import annotations

import argparse
import logging
from typing import Any

logger = logging.getLogger("variantcentrifuge")


def check_scoring_requires_inheritance(args: argparse.Namespace, config: dict[str, Any]) -> bool:
    """Check if scoring configuration requires inheritance analysis.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments
    config : dict
        Configuration dictionary

    Returns
    -------
    bool
        True if inheritance analysis is required for scoring
    """
    scoring_config_path = config.get("scoring_config_path") or getattr(
        args, "scoring_config_path", None
    )
    if not scoring_config_path:
        return False

    try:
        from .scoring import read_scoring_config

        scoring_config = read_scoring_config(scoring_config_path)

        if scoring_config and "formulas" in scoring_config:
            for formula_dict in scoring_config["formulas"]:
                for formula_name, formula_expr in formula_dict.items():
                    if isinstance(formula_expr, str) and any(
                        var in formula_expr
                        for var in [
                            "pattern",
                            "details",
                            "Inheritance_Pattern",
                            "Inheritance_Details",
                        ]
                    ):
                        logger.info(
                            f"Scoring formula '{formula_name}' requires inheritance analysis"
                        )
                        return True

        if scoring_config and "variables" in scoring_config:
            variables = scoring_config["variables"]
            if any(col in variables for col in ["Inheritance_Pattern", "Inheritance_Details"]):
                logger.info("Scoring configuration requires inheritance columns")
                return True

    except Exception as e:
        logger.warning(f"Could not check scoring config for inheritance dependencies: {e}")

    return False
