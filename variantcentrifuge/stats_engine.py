"""Configurable statistics computation engine."""

import ast
import json
import logging
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger("variantcentrifuge")

# AST node types allowed in stats expressions.  Anything outside this set
# (e.g. Import, FunctionDef, Global, Exec, Yield, …) will be rejected
# *before* the expression reaches eval().
_ALLOWED_AST_NODES: frozenset[type] = frozenset(
    {
        # Literals & containers
        ast.Constant,
        ast.List,
        ast.Tuple,
        ast.Dict,
        ast.Set,
        # Variables & attribute access (dunder attrs blocked separately)
        ast.Name,
        ast.Attribute,
        ast.Subscript,
        ast.Starred,
        ast.Index,  # kept for Python 3.8 compat, no-op on 3.9+
        # Expressions
        ast.Expression,
        # Comprehensions
        ast.ListComp,
        ast.SetComp,
        ast.DictComp,
        ast.GeneratorExp,
        ast.comprehension,
        # Operators & comparisons
        ast.BinOp,
        ast.UnaryOp,
        ast.BoolOp,
        ast.Compare,
        ast.IfExp,
        # Operator tokens
        ast.Add,
        ast.Sub,
        ast.Mult,
        ast.Div,
        ast.FloorDiv,
        ast.Mod,
        ast.Pow,
        ast.BitAnd,
        ast.BitOr,
        ast.BitXor,
        ast.Invert,
        ast.Not,
        ast.UAdd,
        ast.USub,
        ast.And,
        ast.Or,
        # Comparison tokens
        ast.Eq,
        ast.NotEq,
        ast.Lt,
        ast.LtE,
        ast.Gt,
        ast.GtE,
        ast.Is,
        ast.IsNot,
        ast.In,
        ast.NotIn,
        # Function calls & keyword args
        ast.Call,
        ast.keyword,
        # Slicing
        ast.Slice,
        # Context (Load/Store/Del markers)
        ast.Load,
        ast.Store,
        ast.Del,
        # f-strings / JoinedStr (harmless in expression context)
        ast.JoinedStr,
        ast.FormattedValue,
    }
)

# Names the restricted eval namespace exposes — only these may appear
# as bare ``ast.Name`` references.
_ALLOWED_NAMES: frozenset[str] = frozenset(
    {
        "df",
        "group_df",
        "pd",
        "np",
        "len",
        "float",
        "str",
        "int",
        "bool",
        "size",
        "True",
        "False",
        "None",
        "sum",
        "min",
        "max",
        "abs",
        "round",
        "sorted",
        "set",
        "list",
        "dict",
        "tuple",
        "range",
        "zip",
        "map",
        "filter",
        "any",
        "all",
        "enumerate",
        "isinstance",
        # Comprehension iteration variables (single-letter by convention)
        *list("abcdefghijklmnopqrstuvwxyz"),
        # Common iteration variable names
        "col",
        "row",
        "val",
        "value",
        "item",
        "idx",
        "key",
        "_",
    }
)


def _validate_expression(expression: str) -> None:
    """Validate a stats expression against the AST allowlist.

    Raises ``ValueError`` if the expression contains disallowed constructs
    such as imports, dunder attribute access, or unsupported AST node types.
    """
    try:
        tree = ast.parse(expression, mode="eval")
    except SyntaxError as exc:
        raise ValueError(f"Syntax error in expression: {exc}") from exc

    for node in ast.walk(tree):
        node_type = type(node)

        # 1. Reject disallowed AST node types
        if node_type not in _ALLOWED_AST_NODES:
            raise ValueError(f"Disallowed construct {node_type.__name__} in expression")

        # 2. Block dunder attribute access (e.g. __class__, __subclasses__)
        if isinstance(node, ast.Attribute) and node.attr.startswith("_"):
            raise ValueError(f"Access to private/dunder attribute '{node.attr}' is not allowed")

        # 3. Restrict bare names to the known namespace
        if isinstance(node, ast.Name) and node.id not in _ALLOWED_NAMES:
            raise ValueError(
                f"Unknown name '{node.id}' in expression; "
                f"allowed names: df, group_df, pd, np, len, float, str, int, bool, size, …"
            )


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

        self.results: dict[str, Any] = {}

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

        The expression is first validated against an AST allowlist that blocks
        imports, dunder attribute access, and other dangerous constructs before
        being passed to ``eval()`` with a restricted namespace.

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
            # Validate the expression AST before executing
            _validate_expression(expression)

            # Create a restricted namespace with only safe objects
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

            # eval with empty __builtins__ AND prior AST validation
            result = eval(expression, {"__builtins__": {}}, namespace)
            return result

        except ValueError as e:
            logger.error(f"Blocked unsafe {context} expression '{expression}': {e}")
            return None
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
                grouped = df.groupby(groupby_cols, observed=True)

                # Apply the expression to each group
                if "size()" in expression:
                    # Special case for counting
                    result = grouped.size()
                else:
                    # General case - apply expression to each group
                    try:
                        # Try with include_groups for newer pandas versions
                        result = grouped.apply(
                            lambda group_df, expression=expression, name=name: self._safe_eval(
                                group_df, expression, f"gene stat '{name}'"
                            ),
                            include_groups=False,
                        )
                    except TypeError:
                        # Fall back for older pandas versions
                        result = grouped.apply(
                            lambda group_df, expression=expression, name=name: self._safe_eval(
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
                grouped = df.groupby(groupby_cols, observed=True)

                # Apply the expression
                if "size()" in expression:
                    result = grouped.size()
                else:
                    try:
                        # Try with include_groups for newer pandas versions
                        result = grouped.apply(
                            lambda group_df, expression=expression, name=name: self._safe_eval(
                                group_df, expression, f"grouped stat '{name}'"
                            ),
                            include_groups=False,
                        )
                    except TypeError:
                        # Fall back for older pandas versions
                        result = grouped.apply(
                            lambda group_df, expression=expression, name=name: self._safe_eval(
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
            for row in self.results["dataset"].itertuples(index=False):
                output_lines.append(f"{getattr(row, 'metric', '')}: {getattr(row, 'value', '')}")
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
