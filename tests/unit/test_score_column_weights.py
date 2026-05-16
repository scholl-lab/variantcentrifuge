"""Tests for arbitrary score-column variant weights."""

from __future__ import annotations

import logging

import numpy as np
import pytest

from variantcentrifuge.association.weights import (
    beta_maf_weights,
    get_weights,
    resolve_score_weight_column,
    score_column_weights,
)


@pytest.fixture
def mafs() -> np.ndarray:
    return np.array([0.01, 0.02, 0.05], dtype=np.float64)


def test_resolve_score_weight_column_from_inline_column():
    assert resolve_score_weight_column("column:nephro_candidate_score", None) == (
        "nephro_candidate_score"
    )


def test_resolve_score_weight_column_inline_wins_over_separate(caplog):
    caplog.set_level(logging.DEBUG, logger="variantcentrifuge")

    result = resolve_score_weight_column("column:inline_score", "separate_score")

    assert result == "inline_score"
    assert "separate_score" in caplog.text


def test_resolve_score_weight_column_from_friendly_form():
    assert resolve_score_weight_column("score_column", "nephro_candidate_score") == (
        "nephro_candidate_score"
    )


def test_resolve_score_weight_column_requires_column_for_friendly_form():
    with pytest.raises(ValueError, match="score_column requires variant_weight_column"):
        resolve_score_weight_column("score_column", None)


def test_resolve_score_weight_column_rejects_empty_inline_column():
    with pytest.raises(ValueError, match="column:<name> requires a non-empty column name"):
        resolve_score_weight_column("column:", None)


def test_score_column_default_combines_beta_with_unit_scale_scores(mafs):
    scores = np.array([0.2, 0.5, 1.0], dtype=np.float64)

    result = score_column_weights(mafs, scores)

    expected = beta_maf_weights(mafs, a=1.0, b=25.0) * scores
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_score_column_explicit_range_floor_and_beta_params(mafs):
    scores = np.array([0.0, 5.0, 10.0], dtype=np.float64)

    result = score_column_weights(
        mafs,
        scores,
        weight_params={
            "score_min": 0,
            "score_max": 10,
            "floor": 0.1,
            "beta_a": 0.5,
            "beta_b": 0.5,
        },
    )

    functional = np.array([0.1, 0.5, 1.0], dtype=np.float64)
    expected = beta_maf_weights(mafs, a=0.5, b=0.5) * functional
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_score_column_raw_score_only_mode(mafs):
    scores = np.array([0.2, 0.5, 2.0], dtype=np.float64)

    result = score_column_weights(
        mafs,
        scores,
        weight_params={"combine_with_beta": False, "ceiling": 1.0},
    )

    np.testing.assert_allclose(result, np.array([0.2, 0.5, 1.0]), rtol=1e-12)


def test_score_column_missing_default_is_neutral_when_combined_with_beta(mafs):
    scores = np.array([0.2, np.nan, "."], dtype=object)

    result = score_column_weights(mafs, scores)

    expected = beta_maf_weights(mafs) * np.array([0.2, 1.0, 1.0])
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_score_column_missing_default_is_floor_in_raw_mode(mafs):
    scores = np.array([0.2, np.nan, "."], dtype=object)

    result = score_column_weights(
        mafs,
        scores,
        weight_params={"combine_with_beta": False, "floor": 0.1},
    )

    np.testing.assert_allclose(result, np.array([0.2, 0.1, 0.1]), rtol=1e-12)


def test_score_column_missing_neutral_invalid_in_raw_mode(mafs):
    scores = np.array([0.2, np.nan, 0.5], dtype=object)

    with pytest.raises(ValueError, match=r"missing='neutral'.*combine_with_beta=false"):
        score_column_weights(
            mafs,
            scores,
            weight_params={"combine_with_beta": False, "missing": "neutral"},
        )


@pytest.mark.parametrize(
    "params, message",
    [
        ({"score_min": 10, "score_max": 0}, "score_min must be less than score_max"),
        ({"floor": -0.1}, "floor must be >= 0"),
        ({"ceiling": 0}, "ceiling must be > 0"),
        ({"floor": 0.8, "ceiling": 0.5}, "floor must be <= ceiling"),
        ({"beta_a": 0}, "beta_a must be > 0"),
        ({"beta_b": -1}, "beta_b must be > 0"),
        ({"missing": "mean"}, "missing must be one of"),
    ],
)
def test_score_column_invalid_params_raise(mafs, params, message):
    with pytest.raises(ValueError, match=message):
        score_column_weights(mafs, np.array([0.1, 0.2, 0.3]), weight_params=params)


def test_score_column_length_mismatch_raises(mafs):
    with pytest.raises(ValueError, match="same length"):
        score_column_weights(mafs, np.array([0.1, 0.2]))


def test_get_weights_dispatches_inline_column(mafs):
    scores = np.array([0.2, 0.5, 1.0], dtype=np.float64)

    result = get_weights(mafs, "column:nephro_candidate_score", score_values=scores)

    expected = score_column_weights(mafs, scores)
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_get_weights_dispatches_score_column_alias(mafs):
    scores = np.array([0.2, 0.5, 1.0], dtype=np.float64)

    result = get_weights(mafs, "score_column", score_values=scores)

    expected = score_column_weights(mafs, scores)
    np.testing.assert_allclose(result, expected, rtol=1e-12)


def test_get_weights_score_column_requires_score_values(mafs):
    with pytest.raises(ValueError, match="requires score_values"):
        get_weights(mafs, "score_column")


def test_get_weights_column_empty_name_raises(mafs):
    with pytest.raises(ValueError, match="column:<name> requires a non-empty column name"):
        get_weights(mafs, "column:", score_values=np.array([0.1, 0.2, 0.3]))
