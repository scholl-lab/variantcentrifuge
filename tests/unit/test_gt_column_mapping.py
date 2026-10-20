"""Per-sample GT columns map to samples by their GEN index, never by position (#219).

Every case uses a subset GEN_0/GEN_2 for samples s0, s1, s2: positional pairing would
give s1 the genotypes of s2.
"""

from __future__ import annotations

import pandas as pd
import pytest

from variantcentrifuge.association.aggregation import _aggregate_gene_burden_from_columns
from variantcentrifuge.gene_burden import perform_gene_burden_analysis
from variantcentrifuge.gt_reconstruction import (
    _collapse_per_sample_gt_columns_for_output,
    _copy_gt_columns_to_samples,
    gt_columns_match_samples,
    reconstruct_gt_column,
)

SAMPLES = ["s0", "s1", "s2"]


def _subset_df(**extra):
    return pd.DataFrame(
        {"GENE": ["G"], "GEN_0__GT": ["0/1"], "GEN_2__GT": ["1/1"], **extra}, index=[0]
    )


@pytest.mark.parametrize(
    ("columns", "ok"),
    [
        (["GEN_0__GT", "GEN_1__GT", "GEN_2__GT"], True),
        (["GEN[0].GT", "GEN[1].GT", "GEN[2].GT"], True),
        (["GEN_0__GT", "GEN_2__GT"], False),
        (["GEN_1__GT", "GEN_0__GT", "GEN_2__GT"], False),
        (["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"], False),
        (["s0", "s1", "s2"], False),
        ([], False),
    ],
)
def test_gt_columns_match_samples(columns, ok):
    assert gt_columns_match_samples(columns, SAMPLES) is ok
    assert gt_columns_match_samples(columns, []) is False


def test_column_aggregation_rejects_unmatched_columns():
    with pytest.raises(ValueError, match="do not map 1:1"):
        _aggregate_gene_burden_from_columns(
            _subset_df(), {"s1"}, {"s0", "s2"}, SAMPLES, ["GEN_0__GT", "GEN_2__GT"]
        )


def test_gene_burden_falls_back_to_packed_gt_for_unmatched_columns():
    cfg = {"gene_burden_mode": "samples"}
    with_cols = _subset_df(GT=["s0(0/1);s2(1/1)"])
    packed_only = with_cols.drop(columns=["GEN_0__GT", "GEN_2__GT"])
    args = ({"s1"}, {"s0", "s2"}, SAMPLES)
    a = perform_gene_burden_analysis(with_cols, cfg, *args)
    b = perform_gene_burden_analysis(packed_only, cfg, *args)
    pd.testing.assert_frame_equal(a, b)
    assert a[["proband_carrier_count", "control_carrier_count"]].values.tolist() == [[0, 2]]


def test_reconstruct_gt_column_maps_by_index():
    out = reconstruct_gt_column(_subset_df(), SAMPLES)
    assert out["GT"].tolist() == ["s0(0/1);s2(1/1)"]


def test_reconstruct_gt_column_extra_fields_follow_the_index():
    df = _subset_df(GEN_0__DP=["10"], GEN_2__DP=["30"])
    out = reconstruct_gt_column(df, SAMPLES, extra_sample_fields=["GEN[*].DP"])
    assert out["GT"].tolist() == ["s0(0/1:10);s2(1/1:30)"]


def test_reconstruct_gt_column_skips_columns_beyond_the_sample_list(caplog):
    out = reconstruct_gt_column(_subset_df(), ["s0", "s1"])
    assert out["GT"].tolist() == ["s0(0/1)"]
    assert "no VCF sample" in caplog.text


def test_collapse_for_output_maps_by_index():
    out = _collapse_per_sample_gt_columns_for_output(_subset_df(), SAMPLES)
    assert out["GT"].tolist() == ["s0(0/1);s2(1/1)"]
    assert not [c for c in out.columns if c.startswith("GEN_")]


def test_collapse_for_output_learns_pseudonyms_by_index():
    df = _subset_df(GT=["p0(0/1);p2(1/1)"], GEN_2__DP=["30"])
    out = _collapse_per_sample_gt_columns_for_output(df, SAMPLES, extra_sample_fields=["DP"])
    assert out["GT"].tolist() == ["p0(0/1:);p2(1/1:30)"]


def test_copy_gt_columns_to_samples_by_index():
    df = _subset_df()
    assert _copy_gt_columns_to_samples(df, ["GEN_0__GT", "GEN_2__GT"], SAMPLES) == 2
    assert df[["s0", "s2"]].values.tolist() == [["0/1", "1/1"]]
    assert "s1" not in df.columns


def test_from_tsv_rejects_gapped_anonymous_columns():
    from variantcentrifuge.from_tsv import infer_and_validate_from_tsv_samples

    df = _subset_df()
    with pytest.raises(ValueError, match="GEN"):
        infer_and_validate_from_tsv_samples(df, {})
