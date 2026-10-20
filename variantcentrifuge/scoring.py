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
from collections.abc import Callable, Sequence
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
        for item in formulas:
            if not isinstance(item, dict) or not item:
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
NUMERIC_SUFFIXES = tuple(
    "_AC _AF _AN _COUNT _P _PHRED _PROB _PVAL _RANKSCORE _SCORE".split()  # noqa: SIM905
)
CATEGORICAL_SUFFIXES = ("_CLNSIG", "_CONSEQUENCE", "_EFFECT", "_IMPACT", "_PATTERN", "_PRED")
NUMERIC_ATTRS = {"abs", "clip", "log", "max", "mean", "min", "pow", "round", "sqrt", "sum"}
STR_ATTRS = set(
    (  # noqa: SIM905
        "astype contains endswith eq fillna ge gt isin le lower lt map ne replace startswith"
        " str strip upper"
    ).split()
)


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


def _is_numeric_col_name(col_name: str) -> bool:
    if _is_known_categorical_column(col_name):
        return False
    u = col_name.upper().strip()
    return u.endswith(NUMERIC_SUFFIXES) or any(
        s in u for s in ("AF", "CADD", "NGS", "SCORE", "PHRED")
    )


def _is_non_bool_numeric(df: pd.DataFrame, col: str | None) -> bool:
    if not col or col not in df.columns:
        return False
    return pd.api.types.is_numeric_dtype(df[col]) and not pd.api.types.is_bool_dtype(df[col])


def _norm_ident(name: str) -> str:
    s = re.sub(r"\W+", "_", name.strip())
    return f"_{s}" if s and s[0].isdigit() else s


def _norm_formula(f: str) -> str:
    return re.sub(r"`([^`]+)`", lambda m: _norm_ident(m.group(1)), f)


def _formula_strings(formulas: list[dict[str, Any]]) -> list[str]:
    return [
        _norm_formula(str(v))
        for f in formulas
        for v in (f.values() if isinstance(f, dict) else [str(f)])
    ]


def _call_args(node: ast.Call) -> list[ast.expr]:
    return [*node.args, *(kw.value for kw in node.keywords)]


def _ast_names(node: ast.AST) -> list[str]:
    return [n.id for n in ast.walk(node) if isinstance(n, ast.Name)]


def _container_elts(op: ast.AST) -> Sequence[ast.AST | None]:
    if isinstance(op, ast.Dict):
        return (*op.keys, *op.values)
    return op.elts if isinstance(op, (ast.List, ast.Tuple, ast.Set)) else [op]


def _has_str_constant(operands: Sequence[ast.AST]) -> bool:
    return any(
        isinstance(e, ast.Constant)
        and isinstance(e.value, str)
        and e.value.strip().lower() not in _MISSING_TOKENS
        for op in operands
        for e in _container_elts(op)
        if e is not None
    )


def _has_num_constant(operands: Sequence[ast.AST]) -> bool:
    return any(
        isinstance(e, ast.Constant) and type(e.value) in (int, float)
        for op in operands
        for e in _container_elts(op)
        if e is not None
    )


def _unwrap_bitwise(e: ast.expr, get_right: bool) -> ast.expr:
    while isinstance(e, ast.BinOp) and isinstance(e.op, (ast.BitAnd, ast.BitOr, ast.BitXor)):
        e = e.right if get_right else e.left
    return e


def _compare_pairs(node: ast.Compare) -> list[tuple[ast.expr, ast.cmpop, ast.expr]]:
    return [
        (_unwrap_bitwise(lhs, True), op, _unwrap_bitwise(rhs, False))
        for lhs, op, rhs in zip(
            [node.left, *node.comparators], node.ops, node.comparators, strict=False
        )
    ]


def _is_compared_to_string_in_formulas(target_name: str, formulas: list[dict[str, Any]]) -> bool:
    target_ident = _norm_ident(target_name)
    formula_strings = _formula_strings(formulas)
    for f_str in formula_strings:
        try:
            tree = ast.parse(f_str, mode="eval")
            for node in ast.walk(tree):
                if isinstance(node, ast.Compare):
                    for lhs, _op, rhs in _compare_pairs(node):
                        operands = [lhs, rhs]
                        if any(
                            target_ident in _ast_names(op) for op in operands
                        ) and _has_str_constant(operands):
                            return True
                elif isinstance(node, ast.Attribute):
                    if getattr(node.value, "id", None) == target_ident and node.attr == "str":
                        return True
                elif isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
                    ops = [node.func.value, *_call_args(node)]
                    if (
                        node.func.attr in STR_ATTRS
                        and any(target_ident in _ast_names(op) for op in ops)
                        and _has_str_constant(ops)
                    ):
                        return True
        except SyntaxError:
            pass

    return False


def _numeric_expr_names(node: ast.AST) -> set[str]:
    if isinstance(node, ast.Name):
        return {node.id}
    if isinstance(node, ast.BinOp):
        return _numeric_expr_names(node.left) | _numeric_expr_names(node.right)
    if isinstance(node, ast.UnaryOp):
        return _numeric_expr_names(node.operand)
    if isinstance(node, ast.Call):
        names = {n for a in _call_args(node) for n in _numeric_expr_names(a)}
        if isinstance(node.func, ast.Name):
            return names
        if isinstance(node.func, ast.Attribute) and node.func.attr in NUMERIC_ATTRS:
            return _numeric_expr_names(node.func.value) | names
    return set()


def _variable_requires_numeric(target_name: str, formulas: list[dict[str, Any]]) -> bool:
    """Check if target_name is used in numeric operations (+, -, *, /, <, math calls)."""
    target_ident = _norm_ident(target_name)
    for f_str in _formula_strings(formulas):
        try:
            tree = ast.parse(f_str, mode="eval")
            nodes: list[ast.AST] = [tree.body]
            while nodes:
                node = nodes.pop()
                if isinstance(node, ast.Compare):
                    for lhs, op, rhs in _compare_pairs(node):
                        operands = [lhs, rhs]
                        if _has_str_constant(operands):
                            continue
                        if (
                            isinstance(op, (ast.Lt, ast.LtE, ast.Gt, ast.GtE))
                            or _has_num_constant(operands)
                        ) and any(target_ident in _numeric_expr_names(o) for o in operands):
                            return True
                        nodes.extend(operands)
                    continue
                if isinstance(node, (ast.BinOp, ast.UnaryOp, ast.Call)):
                    if target_ident in _numeric_expr_names(node):
                        return True
                    if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
                        attr = node.func.attr
                        c_args = _call_args(node)
                        ops = [node.func.value, *c_args]
                        if (
                            attr in ("isin", "eq", "ne")
                            and _has_num_constant(c_args)
                            and target_ident in _numeric_expr_names(node.func.value)
                        ) or (
                            attr in ("lt", "le", "gt", "ge")
                            and not _has_str_constant(ops)
                            and any(target_ident in _numeric_expr_names(op) for op in ops)
                        ):
                            return True
                nodes.extend(ast.iter_child_nodes(node))
        except SyntaxError:
            pass
    return False


def _get_compared_variables(target_name: str, formulas: list[dict[str, Any]]) -> set[str]:
    """Find other variables compared via ==, !=, in, not in, or pandas methods."""
    target_ident = _norm_ident(target_name)
    compared: set[str] = set()
    for f_str in _formula_strings(formulas):
        try:
            tree = ast.parse(f_str, mode="eval")
            for node in ast.walk(tree):
                if isinstance(node, ast.Compare):
                    for lhs, op, rhs in _compare_pairs(node):
                        if isinstance(op, (ast.Eq, ast.NotEq, ast.In, ast.NotIn)):
                            names = [*_ast_names(lhs), *_ast_names(rhs)]
                            if target_ident in names:
                                compared.update(n for n in names if n != target_ident)
                elif (
                    isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr in ("isin", "eq", "ne", "lt", "le", "gt", "ge")
                ):
                    names = [
                        n for op in [node.func.value, *_call_args(node)] for n in _ast_names(op)
                    ]
                    if target_ident in names:
                        compared.update(n for n in names if n != target_ident)
        except SyntaxError:
            pass
    return compared


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
    var_name_to_col: dict[str, str] = {}
    for col, cfg in variables.items():
        parsed = _parse_var_config(str(col), cfg)
        if parsed:
            parsed_vars[str(col)] = parsed
            var_name_to_col[parsed[0]] = str(col)

    norm_to_var = {_norm_ident(k): k for k in var_name_to_col}
    var_partners = {
        cfg[0]: {norm_to_var.get(p, p) for p in _get_compared_variables(cfg[0], formulas)}
        for cfg in parsed_vars.values()
    }

    def get_p_type(p: str) -> str | None:
        c = var_name_to_col.get(p, p)
        if (r := resolved.get(p)) or _is_known_categorical_column(c):
            return r or "categorical"
        return "numeric" if _is_non_bool_numeric(scored_df, c) or _is_numeric_col_name(c) else None

    resolved: dict[str, str | None] = {}
    for original_col, cfg_tuple in parsed_vars.items():
        target_name, default_val_str, reduction_policy, var_type, has_default = cfg_tuple
        is_b = var_type in ("bool", "boolean") or (
            var_type is None
            and original_col in scored_df.columns
            and pd.api.types.is_bool_dtype(scored_df[original_col])
            and not _is_compared_to_string_in_formulas(target_name, formulas)
        )
        if is_b:
            resolved[target_name] = "bool"
        elif var_type in ("numeric", "number", "float", "int"):
            resolved[target_name] = "numeric"
        elif (
            var_type in ("categorical", "string", "str")
            or _is_known_categorical_column(original_col)
            or _is_compared_to_string_in_formulas(target_name, formulas)
        ):
            resolved[target_name] = "categorical"
        elif (
            reduction_policy is not None
            or _is_non_bool_numeric(scored_df, original_col)
            or _variable_requires_numeric(target_name, formulas)
            or _is_numeric_col_name(original_col)
        ):
            resolved[target_name] = "numeric"
        elif var_partners.get(target_name):
            resolved[target_name] = None
        elif not has_default or _is_numeric_default(default_val_str):
            resolved[target_name] = "numeric"
        else:
            resolved[target_name] = "categorical"

    for _ in range(len(parsed_vars)):
        changed = False
        for t, r in resolved.items():
            if r is None:
                pts = {get_p_type(p) for p in var_partners.get(t, ())}
                if "categorical" in pts:
                    resolved[t], changed = "categorical", True
                elif "numeric" in pts:
                    resolved[t], changed = "numeric", True
        if not changed:
            break
    resolved = {t: r or "categorical" for t, r in resolved.items()}

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
                    s_orig.astype(object).map(lambda x, _f=fill_str: _f if pd.isna(x) else str(x))
                    if pd.api.types.is_bool_dtype(s_orig)
                    else s_orig.astype(object).fillna(fill_str).astype(str)
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
