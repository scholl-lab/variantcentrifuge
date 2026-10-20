"""Tests for scoring module."""

import pandas as pd
import pytest

from variantcentrifuge.scoring import apply_scoring


@pytest.fixture
def sample_df():
    """Provide a sample DataFrame for testing."""
    data = {
        "dbNSFP_CADD_phred": [10.0, 25.0, 5.0],
        "dbNSFP_gnomAD_exomes_AF": [0.1, 0.001, 0.0],
        "IMPACT": ["MODERATE", "HIGH", "LOW"],
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_scoring_config():
    """Provide a sample scoring configuration."""
    return {
        "variables": {
            "dbNSFP_CADD_phred": "cadd",
            "dbNSFP_gnomAD_exomes_AF": "freq",
            "IMPACT": "impact_val",
        },
        "formulas": [
            {"pathogenicity_score": "cadd + (1 - freq) * 10"},
            {
                "impact_score": (
                    "((impact_val == 'HIGH') * 10) + ((impact_val == 'MODERATE') * 5) + "
                    "((impact_val == 'LOW') * 1)"
                )
            },
        ],
    }


def test_apply_scoring_calculates_correct_scores(sample_df, sample_scoring_config):
    """Test that apply_scoring adds correct score columns to the DataFrame."""
    scored_df = apply_scoring(sample_df, sample_scoring_config)

    assert "pathogenicity_score" in scored_df.columns
    assert "impact_score" in scored_df.columns

    # Verify calculations
    # Row 1: 10.0 + (1 - 0.1) * 10 = 10 + 9 = 19.0
    # Row 2: 25.0 + (1 - 0.001) * 10 = 25 + 9.99 = 34.99
    # Row 3: 5.0 + (1 - 0.0) * 10 = 5 + 10 = 15.0
    pd.testing.assert_series_equal(
        scored_df["pathogenicity_score"], pd.Series([19.0, 34.99, 15.0]), check_names=False
    )

    # Verify impact scores
    pd.testing.assert_series_equal(
        scored_df["impact_score"], pd.Series([5, 10, 1]), check_names=False
    )


def test_apply_scoring_handles_missing_columns(sample_df, sample_scoring_config):
    """Test that scoring handles missing columns by using defaults."""
    # Remove a column that the scoring config expects
    df_missing_col = sample_df.drop(columns=["dbNSFP_gnomAD_exomes_AF"])

    # Modify config to provide a default
    config_with_default = sample_scoring_config.copy()
    config_with_default["variables"]["dbNSFP_gnomAD_exomes_AF"] = "freq|default:0.5"

    scored_df = apply_scoring(df_missing_col, config_with_default)

    assert "pathogenicity_score" in scored_df.columns
    # All rows should use the default freq of 0.5
    # Row 1: 10 + (1 - 0.5) * 10 = 15
    # Row 2: 25 + (1 - 0.5) * 10 = 30
    # Row 3: 5 + (1 - 0.5) * 10 = 10
    pd.testing.assert_series_equal(
        scored_df["pathogenicity_score"], pd.Series([15.0, 30.0, 10.0]), check_names=False
    )


def test_scoring_hides_intermediate_columns():
    """Test that scoring correctly hides intermediate columns when output_scores is specified."""
    # Create a simple DataFrame
    df = pd.DataFrame(
        {
            "value1": [1.0, 2.0, 3.0],
            "value2": [4.0, 5.0, 6.0],
        }
    )

    # Multi-step config with output_scores
    config = {
        "variables": {
            "value1": "v1",
            "value2": "v2",
        },
        "output_scores": ["final_score"],  # Only keep final_score
        "formulas": [
            {"intermediate1": "v1 * 2"},
            {"intermediate2": "v2 * 3"},
            {"final_score": "intermediate1 + intermediate2"},
        ],
    }

    scored_df = apply_scoring(df, config)

    # Check that only final_score was kept
    assert "final_score" in scored_df.columns
    assert "intermediate1" not in scored_df.columns
    assert "intermediate2" not in scored_df.columns

    # Verify the final score is correct
    # Row 1: (1 * 2) + (4 * 3) = 2 + 12 = 14
    # Row 2: (2 * 2) + (5 * 3) = 4 + 15 = 19
    # Row 3: (3 * 2) + (6 * 3) = 6 + 18 = 24
    pd.testing.assert_series_equal(
        scored_df["final_score"], pd.Series([14.0, 19.0, 24.0]), check_names=False
    )

    # Original columns should still be present
    assert "value1" in scored_df.columns
    assert "value2" in scored_df.columns


def test_scoring_backward_compatibility_without_output_scores():
    """Test that scoring keeps all columns when output_scores is not specified."""
    # Create a simple DataFrame
    df = pd.DataFrame(
        {
            "value1": [1.0, 2.0, 3.0],
            "value2": [4.0, 5.0, 6.0],
        }
    )

    # Multi-step config WITHOUT output_scores
    config = {
        "variables": {
            "value1": "v1",
            "value2": "v2",
        },
        "formulas": [
            {"intermediate1": "v1 * 2"},
            {"intermediate2": "v2 * 3"},
            {"final_score": "intermediate1 + intermediate2"},
        ],
    }

    scored_df = apply_scoring(df, config)

    # Check that all calculated columns are kept
    assert "final_score" in scored_df.columns
    assert "intermediate1" in scored_df.columns
    assert "intermediate2" in scored_df.columns

    # Verify all calculations
    pd.testing.assert_series_equal(
        scored_df["intermediate1"], pd.Series([2.0, 4.0, 6.0]), check_names=False
    )
    pd.testing.assert_series_equal(
        scored_df["intermediate2"], pd.Series([12.0, 15.0, 18.0]), check_names=False
    )
    pd.testing.assert_series_equal(
        scored_df["final_score"], pd.Series([14.0, 19.0, 24.0]), check_names=False
    )


def test_nephro_score_equal_af_reproducibility():
    """Verify legacy and v5 produce identical scores when exome AF == genome AF == joint AF."""
    from variantcentrifuge.scoring import read_scoring_config

    cfg_legacy = read_scoring_config("scoring/nephro_candidate_score")
    cfg_v5 = read_scoring_config("scoring/nephro_candidate_score_v5")

    df_legacy = pd.DataFrame(
        [
            {
                "dbNSFP_gnomAD_exomes_AF": 0.0005,
                "dbNSFP_gnomAD_genomes_AF": 0.0005,
                "dbNSFP_CADD_phred": 20.0,
                "EFFECT": "missense_variant",
                "IMPACT": "MODERATE",
                "ngs": 0.0,
                "Inheritance_Pattern": "unknown",
                "Inheritance_Details": "",
            }
        ]
    )
    df_v5 = pd.DataFrame(
        [
            {
                "dbNSFP_gnomAD4.1_joint_AF": 0.0005,
                "dbNSFP_CADD_phred": 20.0,
                "EFFECT": "missense_variant",
                "IMPACT": "MODERATE",
                "ngs": 0.0,
                "Inheritance_Pattern": "unknown",
                "Inheritance_Details": "",
            }
        ]
    )

    res_legacy = apply_scoring(df_legacy, cfg_legacy)
    res_v5 = apply_scoring(df_v5, cfg_v5)

    score_legacy = res_legacy["nephro_variant_score"].iloc[0]
    score_v5 = res_v5["nephro_variant_score"].iloc[0]

    assert pytest.approx(score_legacy, rel=1e-6) == 0.056695344
    assert pytest.approx(score_v5, rel=1e-6) == 0.056695344
    assert pytest.approx(score_legacy, rel=1e-9) == score_v5


def test_nephro_score_unequal_af_reproducibility():
    """Verify legacy and v5 behavior under unequal AF inputs (Issue #144)."""
    from variantcentrifuge.scoring import read_scoring_config

    cfg_legacy = read_scoring_config("scoring/nephro_candidate_score")
    cfg_v5 = read_scoring_config("scoring/nephro_candidate_score_v5")

    # Exome AF = 0, Genome AF = 0.001, Joint AF = 0.0005
    df_legacy = pd.DataFrame(
        [
            {
                "dbNSFP_gnomAD_exomes_AF": 0.0,
                "dbNSFP_gnomAD_genomes_AF": 0.001,
                "dbNSFP_CADD_phred": 20.0,
                "EFFECT": "missense_variant",
                "IMPACT": "MODERATE",
                "ngs": 0.0,
                "Inheritance_Pattern": "unknown",
                "Inheritance_Details": "",
            }
        ]
    )
    df_v5 = pd.DataFrame(
        [
            {
                "dbNSFP_gnomAD4.1_joint_AF": 0.0005,
                "dbNSFP_CADD_phred": 20.0,
                "EFFECT": "missense_variant",
                "IMPACT": "MODERATE",
                "ngs": 0.0,
                "Inheritance_Pattern": "unknown",
                "Inheritance_Details": "",
            }
        ]
    )

    res_legacy = apply_scoring(df_legacy, cfg_legacy)
    res_v5 = apply_scoring(df_v5, cfg_v5)

    score_legacy = res_legacy["nephro_variant_score"].iloc[0]
    score_v5 = res_v5["nephro_variant_score"].iloc[0]

    assert pytest.approx(score_legacy, rel=1e-6) == 0.443311077
    assert pytest.approx(score_v5, rel=1e-6) == 0.056695344


def test_scoring_numeric_string_conversion_revel_splice_qual():
    """Test that REVEL, splice-score, and QUAL strings are converted to numeric and evaluated."""
    for col in ["dbNSFP_REVEL_score", "splice_dbscSNV_rf_score", "QUAL"]:
        df = pd.DataFrame({col: ["0.92", "."]})
        config = {
            "variables": {col: "v|default:0.0"},
            "formulas": [{"score": "v * 2"}],
        }
        scored_df = apply_scoring(df, config)
        assert "score" in scored_df.columns
        pd.testing.assert_series_equal(
            scored_df["score"], pd.Series([1.84, 0.0]), check_names=False
        )


def test_scoring_reduction_policy_max_and_min():
    """Test that configured max/min reductions work on comma and semicolon lists."""
    df = pd.DataFrame({"score_list": ["0.15,0.92", "25.0;30.0"]})

    # Test max reduction via object config
    config_max = {
        "variables": {"score_list": {"target": "val", "default": 0.0, "reduce": "max"}},
        "formulas": [{"res": "val"}],
    }
    scored_max = apply_scoring(df, config_max)
    pd.testing.assert_series_equal(scored_max["res"], pd.Series([0.92, 30.0]), check_names=False)

    # Test min reduction via string config
    config_min = {
        "variables": {"score_list": "val|default:0.0|reduce:min"},
        "formulas": [{"res": "val"}],
    }
    scored_min = apply_scoring(df, config_min)
    pd.testing.assert_series_equal(scored_min["res"], pd.Series([0.15, 25.0]), check_names=False)


def test_scoring_mixed_tokens_and_missing_and_ambiguous():
    """Test mixed tokens, all-missing inputs, and ambiguous lists without policy."""
    # Mixed missing/numeric retains valid observations
    df_mixed = pd.DataFrame({"score_val": ["0.15,.", ".,0.92", ".;30.0"]})
    config_mixed = {
        "variables": {"score_val": "val|default:0.0"},
        "formulas": [{"res": "val"}],
    }
    scored_mixed = apply_scoring(df_mixed, config_mixed)
    pd.testing.assert_series_equal(
        scored_mixed["res"], pd.Series([0.15, 0.92, 30.0]), check_names=False
    )

    # All-missing numeric inputs use configured default
    df_missing = pd.DataFrame({"score_val": [".", ".,.", ""]})
    config_missing = {
        "variables": {"score_val": "val|default:5.0"},
        "formulas": [{"res": "val"}],
    }
    scored_missing = apply_scoring(df_missing, config_missing)
    pd.testing.assert_series_equal(
        scored_missing["res"], pd.Series([5.0, 5.0, 5.0]), check_names=False
    )

    # Ambiguous lists without policy produce an explicit error
    df_ambiguous = pd.DataFrame({"score_val": ["0.15,0.92"]})
    config_no_policy = {
        "variables": {"score_val": "val|default:0.0"},
        "formulas": [{"res": "val"}],
    }
    with pytest.raises(ValueError, match="reduction"):
        apply_scoring(df_ambiguous, config_no_policy)


def test_scoring_documented_af_cadd_revel_formula():
    """Test that the documented AF/CADD/REVEL formula produces 0.7697 from string inputs."""
    df = pd.DataFrame(
        {
            "dbNSFP_gnomAD_exomes_AF": ["0.001"],
            "dbNSFP_CADD_phred": ["20"],
            "dbNSFP_REVEL_score": ["0.9"],
        }
    )
    config = {
        "variables": {
            "dbNSFP_gnomAD_exomes_AF": "af_exomes|default:0.0",
            "dbNSFP_CADD_phred": "cadd_score|default:0.0",
            "dbNSFP_REVEL_score": "revel_score|default:0.0",
            "ANN[0].IMPACT": "impact|default:'MODIFIER'",
        },
        "formulas": [
            {
                "my_variant_score": (
                    "(cadd_score / 40) * 0.4 + (revel_score * 0.3) + ((1 - af_exomes) * 0.3)"
                )
            }
        ],
    }
    orig_df = df.copy()
    scored_df = apply_scoring(df, config)

    # Caller DataFrame must remain unmodified
    pd.testing.assert_frame_equal(df, orig_df)

    assert "my_variant_score" in scored_df.columns
    assert pytest.approx(scored_df["my_variant_score"].iloc[0], rel=1e-4) == 0.7697


def test_scoring_implicit_numeric_default_when_absent_or_missing():
    """Verify implicit numeric default 0.0 when column is absent
    or all-missing without |default:.
    """
    # Absent column
    df_absent = pd.DataFrame({"dummy": [1, 2]})
    config_absent = {
        "variables": {"absent_col": "v"},
        "formulas": [{"score": "v + 1"}],
    }
    scored_absent = apply_scoring(df_absent, config_absent)
    pd.testing.assert_series_equal(scored_absent["score"], pd.Series([1.0, 1.0]), check_names=False)

    # Column present with only missing/dot values
    df_missing = pd.DataFrame({"present_col": [".", "."]})
    config_missing = {
        "variables": {"present_col": "v"},
        "formulas": [{"score": "v + 1"}],
    }
    scored_missing = apply_scoring(df_missing, config_missing)
    pd.testing.assert_series_equal(
        scored_missing["score"], pd.Series([1.0, 1.0]), check_names=False
    )


def test_scoring_categorical_dtype_preserves_default():
    """Verify category dtype columns fill nulls properly without 'nan' literal string."""
    df = pd.DataFrame(
        {
            "ANN_IMPACT": pd.Series(["HIGH", None], dtype="category"),
        }
    )
    config = {
        "variables": {"ANN_IMPACT": "impact|default:'MODIFIER'"},
        "formulas": [{"is_modifier": "impact == 'MODIFIER'"}],
    }
    scored = apply_scoring(df, config)
    assert scored["ANN_IMPACT"].tolist() == ["HIGH", "MODIFIER"]
    assert scored["is_modifier"].tolist() == [0.0, 1.0]


def test_scoring_coerces_invalid_cells_with_numeric_default():
    """Verify numeric default coerces non-numeric strings like 'N/A' to default."""
    df = pd.DataFrame({"CADD": ["25", "N/A"]})
    config = {
        "variables": {"CADD": "v|default:0"},
        "formulas": [{"score": "v * 2"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [50.0, 0.0]


def test_scoring_type_inference_row_order_invariant():
    """Verify type inference checks entire column and is invariant to row order > 20 rows."""
    # 25 dots followed by a string
    vals1 = ["."] * 25 + ["SOME_TEXT"]
    df1 = pd.DataFrame({"MY_COL": vals1})
    config1 = {
        "variables": {"MY_COL": "val|default:''"},
        "formulas": [{"res": "val == 'SOME_TEXT'"}],
    }
    scored1 = apply_scoring(df1, config1)
    assert scored1["res"].iloc[-1] == 1.0

    # Reversed order: string followed by 25 dots
    vals2 = ["SOME_TEXT"] + ["."] * 25
    df2 = pd.DataFrame({"MY_COL": vals2})
    scored2 = apply_scoring(df2, config1)
    assert scored2["res"].iloc[0] == 1.0


def test_scoring_categorical_dtype_numeric_values_and_none():
    """Verify categorical dtype with numeric values and None uses configured numeric default."""
    df = pd.DataFrame({"SCORE": pd.Series(["0.92", None], dtype="category")})
    config = {
        "variables": {"SCORE": "v|default:7"},
        "formulas": [{"res": "v * 2"}],
    }
    scored = apply_scoring(df, config)
    assert scored["res"].tolist() == [1.84, 14.0]


def test_scoring_target_name_with_categorical_keyword_preserves_numeric():
    """Verify mapping numeric column to target_name containing 'impact' honors numeric default."""
    df = pd.DataFrame({"CADD": ["25", "."]})
    config = {
        "variables": {"CADD": "impact_score|default:0"},
        "formulas": [{"score": "impact_score * 2"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [50.0, 0.0]


@pytest.mark.parametrize(
    "col_name,var_mapping,formula",
    [
        ("IMPACT", "v|default:0", "(v == 'HIGH') * 10"),
        ("EFFECT", "e|default:0", "(e == 'HIGH') * 10"),
        ("Inheritance_Pattern", "ip|default:0", "(ip == 'HIGH') * 10"),
    ],
)
def test_legacy_categorical_mapping_with_numeric_default(col_name, var_mapping, formula):
    """Verify legacy string mapping with numeric default preserves categorical strings and produces
    [10, 0].
    """
    df = pd.DataFrame({col_name: ["HIGH", "."]})
    config = {
        "variables": {col_name: var_mapping},
        "formulas": [{"score": formula}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.parametrize(
    ("values", "assignment", "expected"),
    [
        (["1", "x"], "v|default:''|type:numeric", [1.0, 0.0]),  # 1. explicit beats default
        (["1", "x"], "v|default:0|type:categorical", ["1", "x"]),
        ([True, False], "v", [True, False]),  # 2. bool dtype
        (["HIGH", None], "v|default:0", ["HIGH", "0"]),  # 3. known text column (IMPACT)
        (["1", "x"], "v|default:0.0", [1.0, 0.0]),  # 4. numeric default
        (["1", None], "v|default:''", ["1", ""]),  # 4. string default
        (["1", "x"], "v", [1.0, 0.0]),  # 5. fallback numeric
    ],
)
def test_variable_type_rules(values, assignment, expected):
    """Verify the fixed type order: explicit, bool dtype, known text column, default, numeric."""
    col = "IMPACT" if "HIGH" in values else "CUSTOM"
    scored = apply_scoring(
        pd.DataFrame({col: values}),
        {"variables": {col: assignment}, "formulas": [{"s": "v.notna()"}]},
    )
    assert scored[col].tolist() == expected


@pytest.mark.parametrize(
    ("variables", "formula", "message"),
    [
        ({"CUSTOM": "v|default:''"}, "v * 2", r"'v' is ambiguous \(formulas use it as numeric"),
        ({"CUSTOM": "v|default:0"}, "v == 'A'", "use it as categorical, type rules .* numeric"),
        ({"CUSTOM": "v"}, "v.isin(['A'])", "'v' is ambiguous"),
        ({"CUSTOM": "v|default:''|reduce:max"}, "v.notna()", "use it as numeric"),
        ({"CUSTOM": "v", "OTHER": "w"}, "v == w", "comparison with an untyped variable"),
    ],
)
def test_ambiguous_variable_type_raises(variables, formula, message):
    """Verify the transition guard raises instead of silently changing a variable's type."""
    df = pd.DataFrame({"CUSTOM": ["A", "1"], "OTHER": ["A", "2"]})
    with pytest.raises(ValueError, match=message):
        apply_scoring(df, {"variables": variables, "formulas": [{"s": formula}]})
    typed = {k: f"{v}|type:categorical" for k, v in variables.items()}
    apply_scoring(df, {"variables": typed, "formulas": [{"s": formula}]})  # explicit type fixes it


def test_ambiguity_guard_ignores_missing_tokens_and_syntax_errors():
    """Verify '.' guards and unparsable formulas do not trigger the ambiguity error."""
    df = pd.DataFrame({"CUSTOM": ["20", "."]})
    config = {"variables": {"CUSTOM": "v|default:0"}, "formulas": [{"s": "(v != '.') * v"}]}
    assert apply_scoring(df, config)["s"].tolist() == [20.0, 0.0]
    bad = {"variables": {"CUSTOM": "v"}, "formulas": [{"s": "v +"}]}
    assert apply_scoring(df, bad)["s"].isna().all()


def test_shipped_custom_template_loads_with_comment_keys():
    """'_comment' keys inside formula entries are ignored, so the shipped template loads."""
    from pathlib import Path

    from variantcentrifuge.scoring import read_scoring_config

    root = Path(__file__).resolve().parents[2] / "scoring/coast_classification/custom_template"
    cfg = read_scoring_config(str(root))
    assert all(not k.startswith("_") for item in cfg["formulas"] for k in item)
    assert "coast_category" in cfg["formulas"][-1]
