# ruff: noqa: N999
"""Regression tests for parenthesized categorical comparison operands (Issue #136)."""

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.scoring import apply_scoring


@pytest.mark.unit
def test_parenthesized_variable_comparison_operand():
    """Verify parenthesized variable operand ((v) == 'HIGH') * 10 treats v as categorical."""
    df = pd.DataFrame({"CUSTOM_CAT": ["HIGH", "LOW"]})
    config = {
        "variables": {"CUSTOM_CAT": "v|default:0"},
        "formulas": [{"score": "((v) == 'HIGH') * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_parenthesized_string_comparison_operand():
    """Verify parenthesized string operand (v == ('HIGH')) * 10 treats v as categorical."""
    df = pd.DataFrame({"CUSTOM_CAT": ["HIGH", "LOW"]})
    config = {
        "variables": {"CUSTOM_CAT": "v|default:0"},
        "formulas": [{"score": "(v == ('HIGH')) * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_numeric_isin_comparison():
    """Verify v.isin([20, 30]) * 10 treats numeric string column as numeric."""
    df = pd.DataFrame({"CADD": ["20", "30"]})
    config = {
        "variables": {"CADD": "v|default:0"},
        "formulas": [{"score": "v.isin([20, 30]) * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 10]


@pytest.mark.unit
def test_variable_to_variable_categorical_comparison():
    """Verify (a == b) * 10 compares custom categorical string columns
    without zeroing non-numerics.
    """
    df = pd.DataFrame(
        {
            "ANNOTATION_A": ["D", "D"],
            "ANNOTATION_B": ["D", "T"],
        }
    )
    config = {
        "variables": {
            "ANNOTATION_A": "a|default:0",
            "ANNOTATION_B": "b|default:0",
        },
        "formulas": [{"score": "(a == b) * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]
    assert scored["ANNOTATION_A"].tolist() == ["D", "D"]
    assert scored["ANNOTATION_B"].tolist() == ["D", "T"]


@pytest.mark.unit
def test_all_sentinel_single_row_chunk_numeric_default():
    """Verify all-sentinel single-row chunk with numeric default maintains is_numeric=True."""
    df_single = pd.DataFrame({"dbNSFP_CADD_phred": ["not_available"]})
    config_direct = {
        "variables": {"dbNSFP_CADD_phred": "v|default:7"},
        "formulas": [{"score": "v"}],
    }
    scored_single = apply_scoring(df_single, config_direct)
    assert scored_single["score"].tolist() == [7.0]

    config_abs = {
        "variables": {"dbNSFP_CADD_phred": "v|default:7"},
        "formulas": [{"score": "abs(v)"}],
    }
    scored_single_abs = apply_scoring(df_single, config_abs)
    assert scored_single_abs["score"].tolist() == [7.0]

    df_mixed = pd.DataFrame({"dbNSFP_CADD_phred": ["25", "not_available"]})
    scored_mixed = apply_scoring(df_mixed, config_direct)
    assert scored_mixed["score"].tolist() == [25.0, 7.0]


@pytest.mark.unit
def test_categorical_variable_comparison_independent_of_chunk_composition():
    """Verify (a == b) * 10 produces identical results for full df and row chunks
    even when values contain digits like '1'.
    """
    df = pd.DataFrame(
        {
            "CUSTOM_A": ["D", "D", "1"],
            "CUSTOM_B": ["D", "T", "1"],
        }
    )
    config = {
        "variables": {
            "CUSTOM_A": "a|default:0",
            "CUSTOM_B": "b|default:0",
        },
        "formulas": [{"score": "(a == b) * 10"}],
    }
    # Full dataframe
    scored_full = apply_scoring(df, config)
    assert scored_full["score"].tolist() == [10, 0, 10]

    # Row-by-row chunks
    chunk_scores = []
    for i in range(len(df)):
        chunk = df.iloc[[i]].copy()
        scored_chunk = apply_scoring(chunk, config)
        chunk_scores.extend(scored_chunk["score"].tolist())
    assert chunk_scores == [10, 0, 10]


@pytest.mark.unit
def test_explicit_nan_default_preserves_nan():
    """Verify v|default:nan preserves NaN for missing values and absent columns."""
    df_present = pd.DataFrame({"dbNSFP_CADD_phred": [".", "20"]})
    config = {
        "variables": {"dbNSFP_CADD_phred": "v|default:nan"},
        "formulas": [{"score": "v"}],
    }
    scored_present = apply_scoring(df_present, config)
    scores_present = scored_present["score"].tolist()
    assert np.isnan(scores_present[0])
    assert scores_present[1] == 20.0

    df_absent = pd.DataFrame({"other_col": [1, 2]})
    scored_absent = apply_scoring(df_absent, config)
    scores_absent = scored_absent["score"].tolist()
    assert len(scores_absent) == 2
    assert all(np.isnan(s) for s in scores_absent)


@pytest.mark.unit
def test_categorical_membership_isin():
    """Verify v.isin(w) * 10 and v.isin(['D']) * 10 treat non-numeric variables as categorical."""
    df = pd.DataFrame({"CUSTOM": ["D", "T"], "EXPECTED": ["D", "D"]})
    config_var_partner = {
        "variables": {"CUSTOM": "v", "EXPECTED": "w"},
        "formulas": [{"score": "v.isin(w) * 10"}],
    }
    scored_partner = apply_scoring(df, config_var_partner)
    assert scored_partner["score"].tolist() == [10, 0]

    config_literal = {
        "variables": {"CUSTOM": "v"},
        "formulas": [{"score": "v.isin(['D']) * 10"}],
    }
    scored_literal = apply_scoring(df, config_literal)
    assert scored_literal["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_hgvs_p_categorical_inference():
    """Verify HGVS_P ending with _P is treated as categorical when compared to partner variable."""
    df = pd.DataFrame(
        {
            "HGVS_P": ["p.Ala1Gly", "p.Val2Ala"],
            "EXPECTED": ["p.Ala1Gly", "p.Val3Ala"],
        }
    )
    config = {
        "variables": {"HGVS_P": "v", "EXPECTED": "e"},
        "formulas": [{"score": "(v == e) * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_categorical_string_transformation_formulas():
    """Verify categorical string method transformation (e.g. .replace) evaluates correctly."""
    df = pd.DataFrame({"CUSTOM": ["HIGH", "LOW"]})
    config = {
        "variables": {"CUSTOM": "v|default:0"},
        "formulas": [{"score": "(v.replace('HIGH', 'H') == 'H') * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_numeric_substring_and_partner_numeric_dtype_inference():
    """Verify numeric substrings (e.g. CADD) and partner numeric dtypes infer numeric type."""
    df1 = pd.DataFrame({"CADD": ["20", "30"], "EXPECTED": [20, 30]})
    config1 = {
        "variables": {"CADD": "v|default:0", "EXPECTED": "e|default:0"},
        "formulas": [{"score": "(v == e) * 10"}],
    }
    scored1 = apply_scoring(df1, config1)
    assert scored1["score"].tolist() == [10, 10]

    df2 = pd.DataFrame({"CADD_1": ["0.10"], "CADD_2": ["0.1"]})
    config2 = {
        "variables": {"CADD_1": "v", "CADD_2": "e"},
        "formulas": [{"score": "(v == e) * 10"}],
    }
    scored2 = apply_scoring(df2, config2)
    assert scored2["score"].tolist() == [10]


@pytest.mark.unit
def test_preserve_boolean_columns():
    """Verify boolean input columns are preserved as boolean unless explicit numeric type."""
    df = pd.DataFrame({"FLAG": [True, False]})
    config_and = {
        "variables": {"FLAG": "v"},
        "formulas": [{"score": "(v & v) * 10"}],
    }
    scored_and = apply_scoring(df, config_and)
    assert scored_and["score"].tolist() == [10, 0]

    config_not = {
        "variables": {"FLAG": "v"},
        "formulas": [{"score": "(~v) * 10"}],
    }
    scored_not = apply_scoring(df, config_not)
    assert scored_not["score"].tolist() == [0, 10]

    config_explicit_numeric = {
        "variables": {"FLAG": "v|type:numeric"},
        "formulas": [{"score": "(v + 1) * 10"}],
    }
    scored_num = apply_scoring(df, config_explicit_numeric)
    assert scored_num["score"].tolist() == [20.0, 10.0]


@pytest.mark.unit
def test_missing_value_guard_in_numeric_formula():
    """Verify missing-value guard (v != '.') does not force numeric variable to categorical."""
    df = pd.DataFrame({"CADD": ["20", "."]})
    config = {
        "variables": {"CADD": "v|default:0"},
        "formulas": [{"score": "(v != '.') * v / 40"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [0.5, 0.0]


@pytest.mark.unit
def test_explicit_boolean_parsing():
    """Verify explicit boolean parsing for string tokens, nullable series, and missing column."""
    df_str = pd.DataFrame({"FLAG": ["False", "True", "."]})
    config_bool = {
        "variables": {"FLAG": "v|type:boolean|default:False"},
        "formulas": [{"score": "v * 10"}],
    }
    scored_str = apply_scoring(df_str, config_bool)
    assert scored_str["score"].tolist() == [0, 10, 0]

    df_nullable = pd.DataFrame({"FLAG": pd.Series([True, False, None], dtype="boolean")})
    scored_nullable = apply_scoring(df_nullable, config_bool)
    assert scored_nullable["score"].tolist() == [10, 0, 0]

    df_empty = pd.DataFrame(index=[0, 1, 2])
    scored_missing = apply_scoring(df_empty, config_bool)
    assert scored_missing["score"].tolist() == [0, 0, 0]


@pytest.mark.unit
def test_categorical_comparison_and_str_method_in_numeric_function():
    """Verify numeric-function inference does not force categorical vars to numeric."""
    df1 = pd.DataFrame({"IMPACT": ["HIGH", "LOW"], "CADD": ["16", "9"]})
    config1 = {
        "variables": {"IMPACT": "impact|default:0", "CADD": "cadd|default:0"},
        "formulas": [{"score": 'sqrt((impact == "HIGH") * cadd)'}],
    }
    scored1 = apply_scoring(df1, config1)
    assert scored1["score"].tolist() == [4.0, 0.0]

    df2 = pd.DataFrame({"IMPACT": ["HIGH", "LOW"]})
    config2 = {
        "variables": {"IMPACT": "impact|default:''"},
        "formulas": [{"score": "abs(impact.str.len() - 4)"}],
    }
    scored2 = apply_scoring(df2, config2)
    assert scored2["score"].tolist() == [0, 1]


@pytest.mark.unit
def test_ordered_comparison_against_string_literal():
    """Verify ordered comparisons against string literals are treated as categorical."""
    df = pd.DataFrame({"CUSTOM_GRADE": ["A", "C"]})
    config = {
        "variables": {"CUSTOM_GRADE": "g|default:''"},
        "formulas": [{"score": "(g < 'C') * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_partner_variable_explicit_categorical():
    """Verify partner variable with explicit type categorical prevents numeric inference."""
    df = pd.DataFrame({"CUSTOM": ["HIGH", "LOW"], "CUSTOM_SCORE": ["HIGH", "HIGH"]})
    config = {
        "variables": {
            "CUSTOM": "a|default:''",
            "CUSTOM_SCORE": {"target": "b", "type": "categorical", "default": ""},
        },
        "formulas": [{"score": "(a == b) * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_categorical_preserves_none_pattern():
    """Verify categorical columns preserve 'none' value (e.g. Inheritance_Pattern)."""
    df = pd.DataFrame({"Inheritance_Pattern": ["none", "carrier"]})
    config = {
        "variables": {"Inheritance_Pattern": "p"},
        "formulas": [{"score": "(p == 'none') * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_pandas_comparison_methods():
    """Verify pandas comparison methods (.eq, .ne, .lt) recognize categorical strings."""
    df = pd.DataFrame({"CUSTOM_CLASS": ["HIGH", "LOW"]})
    config_eq = {
        "variables": {"CUSTOM_CLASS": "v|default:0"},
        "formulas": [{"score": "v.eq('HIGH') * 10"}],
    }
    scored_eq = apply_scoring(df, config_eq)
    assert scored_eq["score"].tolist() == [10, 0]

    config_ne = {
        "variables": {"CUSTOM_CLASS": "v|default:0"},
        "formulas": [{"score": "v.ne('LOW') * 10"}],
    }
    scored_ne = apply_scoring(df, config_ne)
    assert scored_ne["score"].tolist() == [10, 0]

    config_lt = {
        "variables": {"CUSTOM_CLASS": "v|default:0"},
        "formulas": [{"score": "v.lt('LOW') * 10"}],
    }
    scored_lt = apply_scoring(df, config_lt)
    assert scored_lt["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_boolean_column_explicit_categorical():
    """Verify boolean column with explicit type categorical behaves as strings."""
    df = pd.DataFrame({"FLAG": [True, False]})
    config = {
        "variables": {"FLAG": {"target": "v", "type": "categorical"}},
        "formulas": [{"score": "(v == 'True') * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]

    df_nullable = pd.DataFrame({"FLAG": pd.Series([True, False, None], dtype="boolean")})
    config_nullable = {
        "variables": {"FLAG": {"target": "v", "type": "categorical", "default": "unknown"}},
        "formulas": [{"score": "(v == 'unknown') * 10"}],
    }
    scored_nullable = apply_scoring(df_nullable, config_nullable)
    assert scored_nullable["score"].tolist() == [0, 0, 10]


@pytest.mark.unit
def test_backtick_quoted_pandas_identifiers():
    """Verify backtick-quoted variables in formulas infer categorical type correctly."""
    df = pd.DataFrame({"CUSTOM_CAT": ["HIGH", "LOW"]})
    config = {
        "variables": {"CUSTOM_CAT": "v|default:0"},
        "formulas": [{"score": "(`v` == 'HIGH') * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_nullable_integer_fractional_default():
    """Verify Int64 nullable integer column works with fractional float default."""
    df = pd.DataFrame({"QUAL": pd.Series([1, None], dtype="Int64")})
    config = {
        "variables": {"QUAL": {"target": "v", "default": 0.5, "type": "numeric"}},
        "formulas": [{"score": "v * 2"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [2.0, 1.0]


@pytest.mark.unit
def test_integer_column_astype_str():
    """Verify non-null integer column preserves integer representation in astype('str')."""
    df = pd.DataFrame({"CUSTOM": pd.Series([1, 2], dtype="int64")})
    config = {
        "variables": {"CUSTOM": "v"},
        "formulas": [{"score": "(v.astype('str') == '1') * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_keyword_arguments_in_ast_inspection():
    """Verify methods with keyword arguments (.eq(other=...), .isin(values=...)) are recognized."""
    df = pd.DataFrame({"CUSTOM_CAT": ["HIGH", "LOW"]})
    config_eq = {
        "variables": {"CUSTOM_CAT": "v|default:0"},
        "formulas": [{"score": "v.eq(other='HIGH') * 10"}],
    }
    scored_eq = apply_scoring(df, config_eq)
    assert scored_eq["score"].tolist() == [10, 0]

    config_isin = {
        "variables": {"CUSTOM_CAT": "v|default:0"},
        "formulas": [{"score": "v.isin(values=['HIGH']) * 10"}],
    }
    scored_isin = apply_scoring(df, config_isin)
    assert scored_isin["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_backtick_variable_name_with_spaces():
    """Verify backtick-quoted variable name with spaces is recognized as categorical in formula."""
    df = pd.DataFrame({"CUSTOM_CAT": ["HIGH", "LOW"]})
    config = {
        "variables": {"CUSTOM_CAT": "Custom Class|default:0"},
        "formulas": [{"score": "(`Custom Class` == 'HIGH') * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_string_concatenation_in_comparison():
    """Verify string concatenation within a string comparison is not treated as numeric."""
    df = pd.DataFrame({"CUSTOM_CAT": ["HIGH", "LOW"]})
    config = {
        "variables": {"CUSTOM_CAT": "v|default:0"},
        "formulas": [{"score": "(v + v == 'HIGHHIGH') * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_chained_or_bitwise_glued_mixed_comparisons():
    """Verify unparenthesized bitwise-glued and chained comparisons infer numeric
    and string types separately.
    """
    from variantcentrifuge.scoring import (
        _is_compared_to_string_in_formulas,
        _variable_requires_numeric,
    )

    df = pd.DataFrame({"CUSTOM_NUM": ["25.0", "10.0"], "CUSTOM_CAT": ["HIGH", "LOW"]})
    config = {
        "variables": {"CUSTOM_NUM": "cadd_v", "CUSTOM_CAT": "impact_v|default:0"},
        "formulas": [{"score": "(cadd_v > 20 & impact_v == 'HIGH') * 10"}],
    }
    assert _variable_requires_numeric("cadd_v", config["formulas"]) is True
    assert _is_compared_to_string_in_formulas("cadd_v", config["formulas"]) is False
    assert _is_compared_to_string_in_formulas("impact_v", config["formulas"]) is True
    assert _variable_requires_numeric("impact_v", config["formulas"]) is False
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_partner_type_inference_is_order_independent_and_transitive():
    """Verify multi-hop comparison partner type inference does not depend on variable dict order."""
    df1 = pd.DataFrame({"A": ["1.50", "2.0"], "B": ["1.5", "3.0"], "C": ["1.5", "2.0"]})
    formulas = [{"score": "(a == b) * 10 + (b == c) * 5 + (c > 1.0) * 1"}]
    cfg_forward = {
        "variables": {"A": "a", "B": "b", "C": "c"},
        "formulas": formulas,
    }
    cfg_reverse = {
        "variables": {"C": "c", "B": "b", "A": "a"},
        "formulas": formulas,
    }
    res_fwd = apply_scoring(df1.copy(), cfg_forward)
    res_rev = apply_scoring(df1.copy(), cfg_reverse)
    assert res_fwd["score"].tolist() == [16.0, 1.0]
    assert res_rev["score"].tolist() == [16.0, 1.0]

    df_cat = pd.DataFrame({"A": ["HIGH", "LOW"], "B": ["HIGH", "MED"], "C": ["HIGH", "HIGH"]})
    formulas_cat = [{"score": "(a == b) * 10 + (b == c) * 5 + (c == 'HIGH') * 1"}]
    cat_fwd = apply_scoring(
        df_cat.copy(),
        {
            "variables": {"A": "a|default:0", "B": "b|default:0", "C": "c|default:0"},
            "formulas": formulas_cat,
        },
    )
    cat_rev = apply_scoring(
        df_cat.copy(),
        {
            "variables": {"C": "c|default:0", "B": "b|default:0", "A": "a|default:0"},
            "formulas": formulas_cat,
        },
    )
    assert cat_fwd["score"].tolist() == [16, 1]
    assert cat_rev["score"].tolist() == [16, 1]


@pytest.mark.unit
def test_unparenthesized_bitwise_mixed_comparisons():
    """Verify cadd > 20 & impact == 'HIGH' evaluates correctly as boolean series."""
    df = pd.DataFrame({"CADD": ["25", "15"], "IMPACT": ["HIGH", "LOW"]})
    config = {
        "variables": {"CADD": "cadd", "IMPACT": "impact"},
        "formulas": [{"score": "cadd > 20 & impact == 'HIGH'"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [True, False]


@pytest.mark.unit
def test_order_independent_comparison_chain():
    """Verify (a == b) * 10 + (b == c) is order-independent when one variable has numeric type."""
    df = pd.DataFrame({"A": ["0.1"], "B": ["0.10"], "C": ["0.1"]})
    cfg_abc = {
        "variables": {"A": "a", "B": "b", "C": "c|type:numeric"},
        "formulas": [{"score": "(a == b) * 10 + (b == c)"}],
    }
    cfg_cba = {
        "variables": {"C": "c|type:numeric", "B": "b", "A": "a"},
        "formulas": [{"score": "(a == b) * 10 + (b == c)"}],
    }
    assert apply_scoring(df, cfg_abc)["score"].tolist() == [11.0]
    assert apply_scoring(df, cfg_cba)["score"].tolist() == [11.0]
