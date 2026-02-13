"""Configurable statistics computation engine."""

import json
import logging
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger("variantcentrifuge")


class StatsEngine:
    """Computes statistics based on JSON configuration."""

    def __init__(self, config: dict[str, Any] | str):
        """
        Initialize the stats engine with configuration.

        Parameters
        ----------
        config : Union[Dict[str, Any], str]
            Either a configuration dictionary or path to JSON config file
        """
        if isinstance(config, str):
            with open(config) as f:
                self.config = json.load(f)
        else:
            self.config = config

        self.results = {}

    def compute(self, df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        """
        Compute all configured statistics.

        Parameters
        ----------
        df : pd.DataFrame
            The variant DataFrame to compute statistics on

        Returns
        -------
        Dict[str, pd.DataFrame]
            Dictionary with keys 'dataset', 'genes', 'groups' containing results
        """
        # Dataset-level stats
        if "dataset_stats" in self.config:
            self.results["dataset"] = self._compute_dataset_stats(df)

        # Gene-level stats
        if "gene_stats" in self.config:
            self.results["genes"] = self._compute_gene_stats(df)

        # Custom groupings
        if "grouped_stats" in self.config:
            self.results["groups"] = self._compute_grouped_stats(df)

        return self.results

    def _check_required_columns(self, df: pd.DataFrame, required: list[str]) -> bool:
        """Check if required columns exist in DataFrame."""
        missing = [col for col in required if col not in df.columns]
        if missing:
            logger.warning(f"Missing required columns: {missing}")
            return False
        return True

    def _safe_eval(self, df: pd.DataFrame, expression: str, context: str = "dataset") -> Any:
        """
        Safely evaluate an expression in the DataFrame context.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame to evaluate expression on
        expression : str
            Expression to evaluate
        context : str
            Context for error messages

        Returns
        -------
        Any
            Result of the expression evaluation
        """
        try:
            # Create a safe namespace with common functions
            namespace = {
                "df": df,
                "group_df": df,  # Support both df and group_df names
                "pd": pd,
                "np": np,
                "len": len,
                "float": float,
                "str": str,
                "int": int,
                "bool": bool,
                "size": lambda: len(df),  # For grouped operations
            }

            # Use eval with restricted namespace
            result = eval(expression, {"__builtins__": {}}, namespace)
            return result

        except Exception as e:
            logger.error(f"Failed to evaluate {context} expression '{expression}': {e}")
            return None

    def _compute_dataset_stats(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute dataset-level statistics."""
        stats_list = []

        for stat_config in self.config.get("dataset_stats", []):
            name = stat_config["name"]
            expression = stat_config["expression"]
            required_cols = stat_config.get("required_columns", [])

            # Check required columns
            if required_cols and not self._check_required_columns(df, required_cols):
                logger.warning(f"Skipping dataset stat '{name}' due to missing columns")
                continue

            # Evaluate expression
            value = self._safe_eval(df, expression, f"dataset stat '{name}'")

            if value is not None:
                stats_list.append({"metric": name, "value": value})

        return pd.DataFrame(stats_list)

    def _compute_gene_stats(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute gene-level statistics."""
        if "GENE" not in df.columns:
            logger.warning("No GENE column found, skipping gene-level stats")
            return pd.DataFrame()

        gene_stats = {}

        for stat_config in self.config.get("gene_stats", []):
            name = stat_config["name"]
            expression = stat_config["expression"]
            required_cols = stat_config.get("required_columns", [])
            groupby_cols = stat_config.get("groupby", "GENE")

            # Ensure groupby is a list
            if isinstance(groupby_cols, str):
                groupby_cols = [groupby_cols]

            # Check required columns
            if required_cols and not self._check_required_columns(df, required_cols):
                logger.warning(f"Skipping gene stat '{name}' due to missing columns")
                continue

            try:
                # Group by specified columns
                grouped = df.groupby(groupby_cols)

                # Apply the expression to each group
                if "size()" in expression:
                    # Special case for counting
                    result = grouped.size()
                else:
                    # General case - apply expression to each group
                    try:
                        # Try with include_groups for newer pandas versions
                        result = grouped.apply(
                            lambda group_df: self._safe_eval(
                                group_df, expression, f"gene stat '{name}'"
                            ),
                            include_groups=False,
                        )
                    except TypeError:
                        # Fall back for older pandas versions
                        result = grouped.apply(
                            lambda group_df: self._safe_eval(
                                group_df, expression, f"gene stat '{name}'"
                            )
                        )

                gene_stats[name] = result

            except Exception as e:
                logger.error(f"Failed to compute gene stat '{name}': {e}")
                continue

        # Combine all stats into a single DataFrame
        if gene_stats:
            try:
                # Try to create DataFrame directly first
                stats_df = pd.DataFrame(gene_stats)
                stats_df = stats_df.reset_index()
                return stats_df
            except ValueError as e:
                if "If using all scalar values, you must pass an index" in str(e):
                    # Handle case where all values are scalars
                    logger.debug("All gene stats are scalar values, creating single-row DataFrame")
                    stats_df = pd.DataFrame([gene_stats])
                    return stats_df
                elif "All arrays must be of the same length" in str(e):
                    # Handle mixed scalar/Series case by converting scalars to single-element Series
                    logger.debug("Mixed scalar/Series gene stats detected, normalizing to Series")
                    normalized_stats = {}
                    for key, value in gene_stats.items():
                        if hasattr(value, "index"):
                            # Already a Series
                            normalized_stats[key] = value
                        else:
                            # Convert scalar to Series with single unnamed index
                            normalized_stats[key] = pd.Series([value], index=["total"])

                    stats_df = pd.DataFrame(normalized_stats)
                    stats_df = stats_df.reset_index()
                    return stats_df
                else:
                    # Re-raise other ValueError types
                    raise
        else:
            return pd.DataFrame()

    def _compute_grouped_stats(self, df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        """Compute statistics with custom groupings."""
        grouped_results = {}

        for stat_config in self.config.get("grouped_stats", []):
            name = stat_config["name"]
            expression = stat_config["expression"]
            groupby_cols = stat_config["groupby"]
            output_format = stat_config.get("output_format", "long")
            required_cols = stat_config.get("required_columns", [])

            # Check required columns
            all_required = required_cols + groupby_cols
            if not self._check_required_columns(df, all_required):
                logger.warning(f"Skipping grouped stat '{name}' due to missing columns")
                continue

            try:
                # Group by specified columns
                grouped = df.groupby(groupby_cols)

                # Apply the expression
                if "size()" in expression:
                    result = grouped.size()
                else:
                    try:
                        # Try with include_groups for newer pandas versions
                        result = grouped.apply(
                            lambda group_df: self._safe_eval(
                                group_df, expression, f"grouped stat '{name}'"
                            ),
                            include_groups=False,
                        )
                    except TypeError:
                        # Fall back for older pandas versions
                        result = grouped.apply(
                            lambda group_df: self._safe_eval(
                                group_df, expression, f"grouped stat '{name}'"
                            )
                        )

                # Format output
                if output_format == "pivot" and len(groupby_cols) == 2:
                    # Create pivot table for 2D groupings
                    result_df = result.unstack(fill_value=0)
                else:
                    # Keep as series/long format
                    result_df = result.to_frame(name="value")
                    result_df = result_df.reset_index()

                grouped_results[name] = result_df

            except Exception as e:
                logger.error(f"Failed to compute grouped stat '{name}': {e}")
                continue

        return grouped_results

    def get_formatted_output(self) -> str:
        """Get formatted string representation of all results."""
        output_lines = []

        # Dataset stats
        if "dataset" in self.results and not self.results["dataset"].empty:
            output_lines.append("=== Dataset Statistics ===")
            for _, row in self.results["dataset"].iterrows():
                output_lines.append(f"{row['metric']}: {row['value']}")
            output_lines.append("")

        # Gene stats
        if "genes" in self.results and not self.results["genes"].empty:
            output_lines.append("=== Gene Statistics ===")
            output_lines.append(self.results["genes"].to_string(index=False))
            output_lines.append("")

        # Grouped stats
        if "groups" in self.results:
            output_lines.append("=== Grouped Statistics ===")
            for name, df in self.results["groups"].items():
                output_lines.append(f"\n{name}:")
                output_lines.append(df.to_string())
            output_lines.append("")

        return "\n".join(output_lines)
