"""Tests for scoring module categorical mapping with numeric default."""

import pandas as pd
import pytest

from variantcentrifuge.scoring import apply_scoring


@pytest.mark.unit
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


@pytest.mark.unit
def test_chrom_scoring_chunk_invariance():
    """Verify CHROM is treated as categorical even in single-row chunks with numeric chromosome
    names.
    """
    df = pd.DataFrame({"CHROM": ["1", "X"]})
    config = {
        "variables": {"CHROM": "chrom"},
        "formulas": [{"score": "(chrom == '1') * 10"}],
    }

    # Whole DataFrame
    scored_full = apply_scoring(df.copy(), config)
    assert scored_full["score"].tolist() == [10, 0]
    assert scored_full["CHROM"].tolist() == ["1", "X"]

    # Chunk by chunk (single row)
    chunk0 = apply_scoring(df.iloc[[0]].copy(), config)
    chunk1 = apply_scoring(df.iloc[[1]].copy(), config)
    assert chunk0["score"].tolist() == [10]
    assert chunk0["CHROM"].tolist() == ["1"]
    assert chunk1["score"].tolist() == [0]
    assert chunk1["CHROM"].tolist() == ["X"]


@pytest.mark.unit
def test_dbnsfp_sift_pred_with_numeric_default():
    """Verify dbNSFP_SIFT_pred with numeric default preserves strings and evaluates formulas."""
    df = pd.DataFrame({"dbNSFP_SIFT_pred": ["D", "T"]})
    config = {
        "variables": {"dbNSFP_SIFT_pred": "p|default:0"},
        "formulas": [{"score": "(p == 'D') * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]
    assert scored["dbNSFP_SIFT_pred"].tolist() == ["D", "T"]


@pytest.mark.unit
def test_inheritance_score_numeric_dtype_preserved():
    """Verify inheritance_score numeric column is treated as numeric despite 'inheritance' in
    name.
    """
    df = pd.DataFrame({"inheritance_score": [0.95, 0.8]})
    config = {
        "variables": {"inheritance_score": "p|default:0"},
        "formulas": [{"score": "p * 2"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == pytest.approx([1.9, 1.6])


@pytest.mark.unit
def test_numeric_scoring_chunk_invariance_with_unparseable_token():
    """Verify dbNSFP_CADD_phred with unparseable tokens and explicit default
    scores identically on whole DataFrame and 1-row chunks.
    """
    df = pd.DataFrame({"dbNSFP_CADD_phred": ["25", "not_available", "."]})
    config = {
        "variables": {"dbNSFP_CADD_phred": "v|default:7"},
        "formulas": [{"score": "v * 2"}],
    }
    scored_full = apply_scoring(df.copy(), config)
    assert scored_full["score"].tolist() == pytest.approx([50.0, 14.0, 14.0])

    chunks = [
        apply_scoring(df.iloc[[0]].copy(), config),
        apply_scoring(df.iloc[[1]].copy(), config),
        apply_scoring(df.iloc[[2]].copy(), config),
    ]
    chunk_scores = [c["score"].iloc[0] for c in chunks]
    assert chunk_scores == pytest.approx([50.0, 14.0, 14.0])


@pytest.mark.unit
def test_formula_string_in_comparison_preserves_categorical():
    """Verify variables compared using 'in' with string literals are treated as categorical."""
    df = pd.DataFrame({"CUSTOM_CAT": ["VAL_A", "VAL_B"]})
    config = {
        "variables": {"CUSTOM_CAT": "cat|default:0|type:categorical"},
        "formulas": [{"score": "(cat in ['VAL_A']) * 10"}],
    }
    scored = apply_scoring(df, config)
    assert scored["score"].tolist() == [10, 0]


@pytest.mark.unit
def test_numeric_scoring_chunk_invariance_without_explicit_default():
    """Verify numeric column without explicit default scores identically on whole DataFrame
    and single-row chunks even with unparseable tokens.
    """
    df = pd.DataFrame({"dbNSFP_CADD_phred": ["25", "not_available", "."]})
    config = {
        "variables": {"dbNSFP_CADD_phred": "v"},
        "formulas": [{"score": "v * 2"}],
    }
    scored_full = apply_scoring(df.copy(), config)
    assert scored_full["score"].tolist() == pytest.approx([50.0, 0.0, 0.0])

    chunks = [
        apply_scoring(df.iloc[[0]].copy(), config),
        apply_scoring(df.iloc[[1]].copy(), config),
        apply_scoring(df.iloc[[2]].copy(), config),
    ]
    chunk_scores = [c["score"].iloc[0] for c in chunks]
    assert chunk_scores == pytest.approx([50.0, 0.0, 0.0])
