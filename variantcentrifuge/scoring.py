"""
Scoring module for applying custom formulas to variant data.

This module reads scoring configurations and applies them to a pandas DataFrame
of annotated variants to generate custom scores.
"""

import ast
import json
import logging
import math
import os
import re
from collections.abc import Callable
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger("variantcentrifuge")


def read_scoring_config(config_path: str) -> dict[str, Any]:
    """Read and parse scoring configuration files from a directory."""
    try:
        logger.debug(f"Reading scoring config from: {config_path}")
        p1 = os.path.join(config_path, "variable_assignment_config.json")
        p2 = os.path.join(config_path, "formula_config.json")
        with open(p1, encoding="utf-8") as f1, open(p2, encoding="utf-8") as f2:
            v_cfg, f_cfg = json.load(f1), json.load(f2)
        variables = v_cfg.get("variables", {}) if isinstance(v_cfg, dict) else None
        formulas = f_cfg.get("formulas", []) if isinstance(f_cfg, dict) else None
        if not isinstance(variables, (dict, list)):
            raise ValueError("Scoring configuration 'variables' must be a JSON object.")
        if not isinstance(formulas, list):
            raise ValueError("Scoring configuration 'formulas' must be a JSON array.")
        if not all(isinstance(item, dict) for item in formulas):
            raise ValueError("Each entry in scoring 'formulas' must be a non-empty object.")
        # "_"-prefixed keys are comments, as at the top level of formula_config.json
        formulas = [{k: v for k, v in item.items() if not k.startswith("_")} for item in formulas]
        for item in formulas:
            if not item:
                raise ValueError("Each entry in scoring 'formulas' must be a non-empty object.")
            for name, expr in item.items():
                if not isinstance(expr, str) or not expr.strip():
                    raise ValueError(f"Formula '{name}' must be a non-empty string expression.")
                try:
                    ast.parse(expr, mode="eval")
                except SyntaxError as syn_err:
                    raise ValueError(f"Invalid syntax in formula '{name}': {syn_err}") from syn_err
        return {
            "variables": variables,
            "formulas": formulas,
            "output_scores": f_cfg.get("output_scores"),
        }
    except FileNotFoundError as e:
        logger.error(f"Scoring configuration file not found: {e.filename}")
        raise
    except json.JSONDecodeError as e:
        logger.error(f"Error parsing scoring configuration JSON: {e}")
        raise


KNOWN_CATEGORICAL_COLS = set(
    "#CHROM ALT BIOTYPE CHROM CLINVAR_CLNSIG CONSEQUENCE EFFECT FEATURE_ID "  # noqa: SIM905
    "FILTER GENE GT HGVS_C HGVS_P ID IMPACT INHERITANCE_DETAILS INHERITANCE_PATTERN REF".split()
)
CATEGORICAL_SUFFIXES = ("_CLNSIG", "_CONSEQUENCE", "_EFFECT", "_IMPACT", "_PATTERN", "_PRED")
_NUMERIC_TYPES = ("numeric", "number", "float", "int")
_TEXT_TYPES = ("categorical", "string", "str")
_ORDERING = (ast.Lt, ast.LtE, ast.Gt, ast.GtE)
_BITWISE = (ast.BitAnd, ast.BitOr, ast.BitXor)
_CMP_METHODS = ("isin", "eq", "ne", "lt", "le", "gt", "ge")

ANN_CAT_RE = re.compile(
    r"^ANN\[\d+\]\.(IMPACT|EFFECT|GENE|FEATUREID|BIOTYPE|HGVS_C|HGVS_P)$",
    re.IGNORECASE,
)

_MISSING_TOKENS = frozenset({"", ".", "na", "n/a", "nan", "none", "null", "inf", "-inf", "+inf"})
_TRUE_STRS = frozenset({"true", "1", "1.0", "yes", "t", "y"})
_FALSE_STRS = frozenset({"false", "0", "0.0", "no", "f", "n"})


def _parse_default_bool(default_val: Any) -> bool:
    if isinstance(default_val, str):
        return default_val.strip().lower() in _TRUE_STRS
    return False if pd.isna(default_val) else bool(default_val)


def _convert_to_bool(series: pd.Series, default_val: Any) -> pd.Series:
    def_b = _parse_default_bool(default_val)

    def _to_b(x: Any) -> bool:
        if pd.isna(x):
            return def_b
        if isinstance(x, (int, float, np.number, bool, np.bool_)):
            return bool(x)
        s = str(x).strip().lower()
        return True if s in _TRUE_STRS else (False if s in _FALSE_STRS else def_b)

    s_obj = series.astype(object) if isinstance(series.dtype, pd.CategoricalDtype) else series
    return s_obj.apply(_to_b).astype(bool)


def _is_known_categorical_column(col_name: str) -> bool:
    u = col_name.upper().strip()
    return bool(
        u in KNOWN_CATEGORICAL_COLS or u.endswith(CATEGORICAL_SUFFIXES) or ANN_CAT_RE.match(u)
    )


def _formula_usage(name: str, formulas: list[dict[str, Any]], typed: dict[str, str]) -> set[str]:
    """Return how formulas use variable *name* ("numeric", "categorical", ...).

    *typed* maps every variable to its explicit or known-column type ("" when untyped).
    ponytail: transition guard for rule-based variable typing (only detects disagreement,
    never chooses a type); delete in a later release once configs carry explicit types.
    """
    norm = lambda s: re.sub(r"\W", "_", s)  # noqa: E731
    name = norm(name)

    def has_var(e: ast.AST, anywhere: bool = False) -> bool:
        nodes = ast.walk(e) if anywhere else [e]
        return any(isinstance(x, ast.Name) and x.id == name for x in nodes)

    def use_of(other: ast.AST, op: ast.AST | None = None) -> str:
        if isinstance(other, (ast.List, ast.Tuple, ast.Set)):
            uses = {use_of(e) for e in other.elts}
            return "categorical" if "categorical" in uses else max(uses, default="")
        val = other.value if isinstance(other, ast.Constant) else None
        if isinstance(val, str):
            return "" if val.strip().lower() in _MISSING_TOKENS else "categorical"
        if type(val) in (int, float) or isinstance(op, _ORDERING):
            return "numeric"
        if isinstance(other, ast.Name) and other.id in typed:
            return typed[other.id] or "comparison with an untyped variable"
        return ""

    used: set[str] = set()
    for expr in (str(e) for f in formulas for e in f.values()):
        try:
            tree = ast.parse(re.sub(r"`([^`]+)`", lambda m: norm(m.group(1)), expr), mode="eval")
        except SyntaxError:
            continue
        for n in ast.walk(tree):
            if isinstance(n, ast.BinOp) and not isinstance(n.op, _BITWISE):
                used.update(["numeric"] if has_var(n.left) or has_var(n.right) else [])
            elif isinstance(n, ast.Compare):
                sides = [n.left, *n.comparators]
                for lhs, op, rhs in zip(sides, n.ops, sides[1:], strict=False):
                    for me, other in ((lhs, rhs), (rhs, lhs)):
                        if has_var(me) or (has_var(me, True) and use_of(other) == "categorical"):
                            used.add(use_of(other, op))
            elif isinstance(n, ast.Attribute) and has_var(n.value) and n.attr == "str":
                used.add("categorical")
            elif (
                isinstance(n, ast.Call)
                and isinstance(n.func, ast.Attribute)
                and has_var(n.func.value)
                and n.func.attr in _CMP_METHODS
            ):
                args = [*n.args, *(k.value for k in n.keywords)]
                used.update(use_of(c) for a in args for c in ast.walk(a))
    return used - {""}


def _is_numeric_default(val_str: str) -> bool:
    s = val_str.strip()
    if len(s) >= 2 and s[0] in ("'", '"') and s[-1] == s[0]:
        return False
    return math.isfinite(_to_finite_float(s, float("nan"))) or s.lower() in ("nan", "+nan", "-nan")


_REDUCTION_FUNCS: dict[str, Callable[[list[float]], float]] = {
    "max": max,
    "min": min,
    "sum": sum,
    "mean": lambda x: float(sum(x) / len(x)),
    "first": lambda x: x[0],
    "last": lambda x: x[-1],
}


def _to_finite_float(val: Any, default: float) -> float:
    try:
        f = float(val)
        return f if math.isfinite(f) else default
    except (ValueError, TypeError):
        return default


def _parse_default_float(val: Any) -> float:
    """Parse default float, preserving explicit nan strings/floats."""
    is_nan = (isinstance(val, float) and np.isnan(val)) or (
        isinstance(val, str) and val.strip().lower() in ("nan", "+nan", "-nan", ".nan")
    )
    return float("nan") if is_nan else _to_finite_float(val, 0.0)


def _parse_numeric_scalar(
    val: Any,
    default: float,
    reduction: str | None,
    col_name: str,
) -> float:
    """Parse a single value or delimited string into a numeric float."""
    if val is None or pd.isna(val):
        return default
    if isinstance(val, (int, float, np.number)):
        return _to_finite_float(val, default)
    if isinstance(val, str):
        if (val := val.strip()) in _MISSING_TOKENS:
            return default
        if "," in val or ";" in val:
            tokens = [t.strip() for t in val.replace(";", ",").split(",")]
            valid_floats = [
                f
                for t in tokens
                if t not in _MISSING_TOKENS
                if math.isfinite(f := _to_finite_float(t, float("nan")))
            ]
            if not valid_floats:
                return default
            if len(valid_floats) == 1:
                return valid_floats[0]
            if not reduction:
                raise ValueError(
                    f"Ambiguous multi-value numeric string '{val}' in column '{col_name}' "
                    f"requires an explicit reduction policy ('max', 'min', etc.)."
                )
            red = reduction.lower()
            if red in _REDUCTION_FUNCS:
                return _REDUCTION_FUNCS[red](valid_floats)
            raise ValueError(f"Unsupported reduction policy: '{reduction}'")
    return _to_finite_float(val, default)


def convert_to_numeric(
    series: pd.Series,
    default: float = 0.0,
    reduction: str | None = None,
) -> pd.Series:
    """Convert a series to numeric, handling missing tokens and multi-transcript values."""
    if (
        pd.api.types.is_numeric_dtype(series)
        and not pd.api.types.is_bool_dtype(series)
        and reduction is None
    ):
        return series.astype(float).replace([np.inf, -np.inf], default).fillna(default)
    col_name = str(series.name) if series.name is not None else ""
    s_obj = series.astype(object) if isinstance(series.dtype, pd.CategoricalDtype) else series
    return s_obj.apply(lambda val: _parse_numeric_scalar(val, default, reduction, col_name)).astype(
        float
    )


def _clean_default_str(val: str) -> str:
    s = str(val).strip()
    return s[1:-1] if len(s) >= 2 and s[0] in ("'", '"') and s[-1] == s[0] else s


def _parse_var_config(
    original_col: str, var_config: Any
) -> tuple[str, str, str | None, str | None, bool] | None:
    if not isinstance(var_config, (str, dict)):
        return None
    target_name = str(original_col)
    default_val_str, reduction_policy, var_type, has_default = "0", None, None, False
    if isinstance(var_config, str):
        parts = var_config.split("|")
        target_name = parts[0].strip() or target_name
        for part in parts[1:]:
            if ":" in part:
                k, v = (x.strip() for x in part.split(":", 1))
                if (kl := k.lower()) == "default":
                    default_val_str, has_default = v, True
                elif kl in ("reduce", "reduction"):
                    reduction_policy = v
                elif kl == "type":
                    var_type = v.lower()
    else:
        target_name = str(var_config.get("target", original_col))
        if "default" in var_config:
            default_val_str, has_default = str(var_config["default"]), True
        reduction_policy = var_config.get("reduce") or var_config.get("reduction")
        var_type = str(var_config.get("type", "")).lower() or None
    return target_name, default_val_str, reduction_policy, var_type, has_default


def apply_scoring(
    df: pd.DataFrame, scoring_config: dict[str, Any], *, strict: bool = False
) -> pd.DataFrame:
    """Apply scoring formulas to variant DataFrame, evaluating and adding score columns."""
    variables = scoring_config.get("variables", {})
    formulas = scoring_config.get("formulas", [])

    if not formulas:
        logger.debug("No scoring formulas provided. Skipping scoring.")
        return df

    logger.info(f"Applying {len(formulas)} scoring formula(s)...")

    scored_df = df.copy()
    calculated_columns: list[str] = []
    parsed_vars: dict[str, tuple[str, str, str | None, str | None, bool]] = {}
    for col, cfg in variables.items():
        parsed = _parse_var_config(str(col), cfg)
        if parsed:
            parsed_vars[str(col)] = parsed

    # Fixed order: explicit type, bool dtype, known text column, default value type, numeric.
    resolved: dict[str, str] = {}
    for original_col, (target_name, default_str, _r, var_type, has_default) in parsed_vars.items():
        col = scored_df.get(original_col)
        if var_type in ("bool", "boolean") or (
            var_type is None and col is not None and pd.api.types.is_bool_dtype(col)
        ):
            resolved[target_name] = "bool"
        elif var_type in _NUMERIC_TYPES:
            resolved[target_name] = "numeric"
        elif (
            var_type in _TEXT_TYPES
            or _is_known_categorical_column(original_col)
            or (has_default and not _is_numeric_default(default_str))
        ):
            resolved[target_name] = "categorical"
        else:
            resolved[target_name] = "numeric"
    typed = {
        cfg[0]: resolved[cfg[0]] if cfg[3] or _is_known_categorical_column(col) else ""
        for col, cfg in parsed_vars.items()
    }
    # ponytail: transition guard (see _formula_usage) - raise where the old AST guesser's
    # precedence (string use > numeric use/reduce/numeric dtype > variable comparison) differs.
    for original_col, (target_name, _d, red, var_type, _h) in parsed_vars.items():
        res, col = resolved[target_name], scored_df.get(original_col)
        used = _formula_usage(target_name, formulas, typed)
        if red is not None or (col is not None and col.dtype.kind in "iuf"):
            used.add("numeric")
        old = next((u for u in ("categorical", "numeric") if u in used), min(used, default=res))
        if var_type or _is_known_categorical_column(original_col) or old == res:
            continue
        if res != "bool" or old == "categorical":
            raise ValueError(
                f"Scoring variable '{target_name}' is ambiguous (formulas use it as {old}, "
                f"type rules resolve it to {res}): add '|type:numeric' or "
                "'|type:categorical' to its assignment in variable_assignment_config.json"
            )

    rename_map, created_vars = {}, []
    for original_col, cfg_tuple in parsed_vars.items():
        target_name, default_val_str, reduction_policy, var_type, has_default = cfg_tuple
        var_res = resolved[target_name]
        val_to_use = _clean_default_str(default_val_str)
        if original_col in scored_df.columns:
            if var_res == "bool":
                scored_df[original_col] = _convert_to_bool(scored_df[original_col], default_val_str)
            elif var_res == "numeric":
                s_col = scored_df[original_col]
                if not (
                    pd.api.types.is_integer_dtype(s_col)
                    and not s_col.isna().any()
                    and reduction_policy is None
                    and var_type not in ("float",)
                ):
                    scored_df[original_col] = convert_to_numeric(
                        s_col, _parse_default_float(default_val_str), reduction_policy
                    )
            else:
                fill_str = val_to_use if has_default else ""
                s_orig = scored_df[original_col]
                s_cat = (
                    s_orig.astype(object)
                    .map(lambda x, _f=fill_str: _f if pd.isna(x) else str(x))
                    .astype(str)
                )
                scored_df[original_col] = s_cat.replace([".", "NA", ""], fill_str)
            rename_map[original_col] = target_name
        else:
            logger.warning(
                f"Variable column '{original_col}' not found. Using default '{default_val_str}'."
            )
            scored_df[target_name] = (
                pd.Series(_parse_default_bool(default_val_str), index=scored_df.index, dtype=bool)
                if var_res == "bool"
                else (
                    _parse_default_float(default_val_str)
                    if var_res == "numeric"
                    else (val_to_use if has_default else "")
                )
            )
            created_vars.append(target_name)

    if rename_map:
        scored_df.rename(columns=rename_map, inplace=True)

    for formula_item in formulas:
        for score_name, formula_str in formula_item.items():
            logger.debug(f"Evaluating formula for '{score_name}': {formula_str}")
            try:
                scored_df[score_name] = scored_df.eval(
                    formula_str, engine="python", local_dict={"np": np}
                )
                if not pd.api.types.is_numeric_dtype(scored_df[score_name]):
                    scored_df[score_name] = pd.to_numeric(scored_df[score_name], errors="coerce")
                logger.info(f"Successfully calculated score column: '{score_name}'")
            except Exception as e:
                logger.error(f"Error evaluating formula for '{score_name}': {e}")
                if strict:
                    raise
                scored_df[score_name] = float("nan")
            calculated_columns.append(score_name)

    if rename_map:
        scored_df.rename(columns={v: k for k, v in rename_map.items()}, inplace=True)
    if created_vars:
        scored_df.drop(columns=created_vars, inplace=True)

    output_scores = scoring_config.get("output_scores")
    if output_scores is not None:
        columns_to_drop = [col for col in calculated_columns if col not in output_scores]
        if columns_to_drop:
            logger.info(f"Dropping intermediate columns: {columns_to_drop}")
            scored_df.drop(columns=columns_to_drop, inplace=True)

    return scored_df
