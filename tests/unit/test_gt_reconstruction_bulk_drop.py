"""#173: packed GT reconstruction drops per-sample columns in one operation."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.gt_reconstruction import reconstruct_gt_column


def _frame(categorical: bool) -> pd.DataFrame:
    gts = {
        "GEN[0].GT": ["0/0", "0/1", None, " 1|1 ", np.str_("0/0")],
        "GEN[1].GT": ["./.", "0|0", "1", "NA", np.nan],
        "GEN[2].GT": ["1/1", "0/0", "0/1", "", "0|1"],
    }
    df = pd.DataFrame({"CHROM": ["1"] * 5, "POS": [1, 2, 3, 4, 5]})
    for col, vals in gts.items():
        df[col] = pd.Categorical(vals) if categorical else pd.Series(vals, dtype=object)
    df["GEN[0].DP"] = ["10", "11", None, "13", "14"]
    df["GEN[2].DP"] = ["7", "8", "9", " 10 ", None]
    df["QUAL"] = [50.0, 60.0, 70.0, 80.0, 90.0]
    return df


EXPECTED_GT = [
    "B(./.);C(1/1)",
    "A(0/1)",
    "A(./.);B(1);C(0/1)",
    "A(1|1);B(./.);C(./.)",
    "B(./.);C(0|1)",
]
EXPECTED_GT_DP = [
    "A(0/0:10);B(./.:);C(1/1:7)",
    "A(0/1:11);B(0|0:);C(0/0:8)",
    "A(./.:);B(1:);C(0/1:9)",
    "A(1|1:13);B(./.:);C(./.:10)",
    "A(0/0:14);B(./.:);C(0|1:)",
]


@pytest.mark.unit
@pytest.mark.parametrize("categorical", [False, True])
def test_reconstruction_output_columns_and_values(categorical: bool) -> None:
    df = _frame(categorical)
    before = df.copy()
    out = reconstruct_gt_column(df, ["S0", "S1", "S2"], sample_mapping={"S0": "A", "S1": "B"})
    assert out["GT"].tolist() == [s.replace("C(", "S2(") for s in EXPECTED_GT]
    assert list(out.columns) == ["CHROM", "POS", "GEN[0].DP", "GEN[2].DP", "QUAL", "GT"]
    assert out["QUAL"].dtype == np.float64 and out["POS"].dtype == before["POS"].dtype
    pd.testing.assert_frame_equal(df, before)  # input frame untouched


@pytest.mark.unit
def test_reconstruction_with_extra_fields_drops_them_too() -> None:
    df = _frame(categorical=False)
    out = reconstruct_gt_column(
        df, ["A", "B", "C"], extra_sample_fields=["GEN[*].DP"], delimiter=":"
    )
    assert out["GT"].tolist() == EXPECTED_GT_DP
    assert list(out.columns) == ["CHROM", "POS", "QUAL", "GT"]


@pytest.mark.unit
def test_reconstruction_deletes_sample_columns_in_one_operation(monkeypatch) -> None:
    """Per-column ``del`` re-splits the remaining blocks every time (quadratic)."""
    calls = []
    original = pd.DataFrame.__delitem__

    def _counting_delitem(self, key):
        calls.append(key)
        original(self, key)

    monkeypatch.setattr(pd.DataFrame, "__delitem__", _counting_delitem)
    df = _frame(categorical=True)
    out = reconstruct_gt_column(df, ["A", "B", "C"])
    assert out["GT"].tolist() == EXPECTED_GT
    assert len(calls) <= 1  # one bulk drop/delete, never one per sample column
    # Kept columns stay shared with the input (shallow-copy memory contract).
    assert np.shares_memory(df["QUAL"].to_numpy(), out["QUAL"].to_numpy())
