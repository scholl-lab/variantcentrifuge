"""#176: vocabulary-factorised GT parsing is bit-identical to the per-column parser."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.genotype_matrix import (
    _multiallelic_error,
    _vectorized_parse_gt_columns,
)


def _reference_parse(gene_df: pd.DataFrame, cols: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Frozen per-column implementation (v0.19.0 semantics plus #175 rejection)."""
    n_variants, n_samples = len(gene_df), len(cols)
    raw = np.full((n_variants, n_samples), np.nan, dtype=np.float64)
    ploidy = np.full((n_variants, n_samples), 2.0, dtype=np.float64)
    row_has_haploid = np.zeros(n_variants, dtype=bool)
    row_has_diploid = np.zeros(n_variants, dtype=bool)
    generic_missing = np.zeros((n_variants, n_samples), dtype=bool)
    for s_idx, col in enumerate(cols):
        s_obj = (
            gene_df[col].astype(object)
            if isinstance(gene_df[col].dtype, pd.CategoricalDtype)
            else gene_df[col]
        )
        s_str = s_obj.fillna("").astype(str).str.strip()
        is_gen_miss = s_obj.isna() | s_str.str.lower().isin({"", "na", "nan", "none", "null"})
        generic_missing[:, s_idx] = is_gen_miss.to_numpy()
        gt_series = s_str.str.replace("|", "/", regex=False)
        m_diploid = ~is_gen_miss & gt_series.str.contains("/", regex=False)
        row_has_diploid |= m_diploid.to_numpy()
        m_00 = gt_series == "0/0"
        m_01 = (gt_series == "0/1") | (gt_series == "1/0")
        m_11 = gt_series == "1/1"
        m_hap_0 = gt_series == "0"
        m_hap_1 = gt_series == "1"
        m_hap_miss = gt_series == "."
        row_has_haploid |= (m_hap_0 | m_hap_1 | m_hap_miss).to_numpy()
        alleles = gt_series.str.extract(r"^([0-9]+)(?:/([0-9]+))?$").astype(float)
        multi = (alleles > 1).any(axis=1).to_numpy()
        if multi.any():
            raise _multiallelic_error(gene_df, int(np.flatnonzero(multi)[0]), col)
        col_dosage = np.full(n_variants, np.nan, dtype=np.float64)
        col_ploidy = np.full(n_variants, 2.0, dtype=np.float64)
        col_dosage[m_00.values] = 0.0
        col_dosage[m_01.values] = 1.0
        col_dosage[m_11.values] = 2.0
        col_dosage[m_hap_0.values] = 0.0
        col_dosage[m_hap_1.values] = 1.0
        col_ploidy[m_hap_0.values] = 1.0
        col_ploidy[m_hap_1.values] = 1.0
        col_ploidy[m_hap_miss.values] = 1.0
        raw[:, s_idx] = col_dosage
        ploidy[:, s_idx] = col_ploidy
    haploid_rows = row_has_haploid & ~row_has_diploid
    if haploid_rows.any():
        ploidy[haploid_rows[:, None] & generic_missing] = 1.0
    return raw, ploidy


EDGE_STRINGS = [
    "0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1", " 0/1 ", "\t1/1", "0/1\n",
    "./.", ".|.", ".", "./1", "1/.", ".|1", "0", "1", " 1", "", " ", "NA", "na", "nan",
    "NaN", "None", "NULL", "null", "0/01", "01", "00", "0/1/1", "0//1", "/", "|", "abc",
    "-1", "+1", "1/-1", "0/0:12", "0/1 0/1",
]  # fmt: skip
NA_OBJECTS = [None, np.nan, pd.NA, pd.NaT, np.str_("0/1"), np.str_(" 1 ")]
NON_STR_OBJECTS = [0, 1, 0.0, 1.0, True, False, b"0/1", np.int64(1), np.float32(0.0)]
MULTI_STRINGS = ["0/2", "1/2", "2/2", "2|1", "2", " 3 ", "0/10", "10"]


def _assert_same(df: pd.DataFrame, cols: list[str]) -> None:
    try:
        expected = _reference_parse(df, cols)
    except ValueError as exc:
        with pytest.raises(ValueError) as got:
            _vectorized_parse_gt_columns(df, cols)
        assert str(got.value) == str(exc)
        return
    raw, ploidy = _vectorized_parse_gt_columns(df, cols)
    for new, ref in zip((raw, ploidy), expected, strict=True):
        assert new.dtype == ref.dtype and new.shape == ref.shape
        np.testing.assert_array_equal(np.isnan(new), np.isnan(ref))
        np.testing.assert_array_equal(new, ref)


def _random_frame(rng: np.random.Generator, n_rows: int, n_cols: int, pool: list) -> pd.DataFrame:
    data = {}
    for j in range(n_cols):
        # Mostly standard calls with a few edge values, so haploid/diploid rows mix.
        vals = [pool[k] for k in rng.integers(0, len(pool), n_rows)]
        std = rng.choice(["0/0", "0/1", "1/1", "./.", "0|1", "0", "1", "."], n_rows)
        keep_std = rng.random(n_rows) < 0.6
        col = [s if k else v for s, k, v in zip(std, keep_std, vals, strict=True)]
        kind = j % 4
        if kind == 1 and all(isinstance(v, str) for v in col):
            data[f"c{j}"] = pd.Categorical(col)
        elif kind == 2 and all(isinstance(v, str) or v is None for v in col):
            data[f"c{j}"] = pd.array(col, dtype="string")
        else:
            data[f"c{j}"] = pd.Series(col, dtype=object)
    return pd.DataFrame(data)


@pytest.mark.unit
@pytest.mark.parametrize("seed", range(40))
def test_randomized_vocabulary_matches_reference(seed: int) -> None:
    rng = np.random.default_rng(seed)
    # Even seeds: strings + NA only (factorised fast path); odd: non-string objects too.
    pool = EDGE_STRINGS + NA_OBJECTS + (NON_STR_OBJECTS if seed % 2 else [])
    pool += MULTI_STRINGS if seed % 5 == 0 else []
    df = _random_frame(rng, int(rng.integers(1, 30)), int(rng.integers(1, 25)), pool)
    if seed % 3 == 0:
        df.index = [0] * len(df)  # duplicate index labels
    _assert_same(df, list(df.columns))


@pytest.mark.unit
@pytest.mark.parametrize(
    "column",
    [
        pd.Series([0, 1, 1, 0], dtype="int64"),
        pd.Series([0.0, 1.0, np.nan, 1.0], dtype="float64"),
        pd.Series([0.0, 1.0, 2.5, 1.0], dtype="float64"),
        pd.Series([True, False, True, True]),
        pd.Series(["0/1", None, "1", "."], dtype="string"),
        pd.Series(pd.Categorical(["0/1", None, "1/1", "0/0"], categories=["0/0", "0/1", "1/1"])),
        pd.Series(pd.Categorical([0, 1, None, 1])),
        pd.Series([1, 0, 1.0, None], dtype=object),
        pd.Series([0, 1, 0, 1], dtype=object),
        # No NA: the per-column fillna downcasts to float64 -> "1.0", never haploid "1".
        pd.Series([1, 1.0, 0, 1], dtype=object),
        pd.Series([0, 2, 0, 0], dtype=object),
        pd.Series([b"0/1", b"1", None, "1/1"], dtype=object),
    ],
)
def test_non_string_dtypes_match_reference(column: pd.Series) -> None:
    df = pd.DataFrame({"a": ["0/1", "0/0", "1", "./."]})
    df["b"] = column.values
    _assert_same(df, ["a", "b"])


@pytest.mark.unit
@pytest.mark.parametrize("gt", MULTI_STRINGS)
def test_multiallelic_first_offender_matches_reference(gt: str) -> None:
    """The reported offender (first sample column, then row) is unchanged."""
    df = pd.DataFrame({"a": ["0/0", "0/1", gt], "b": [gt, "0/0", "0/0"], "c": ["1", "0", gt]})
    _assert_same(df, ["c", "a", "b"])
    _assert_same(df, ["a", "b"])


@pytest.mark.unit
def test_empty_inputs() -> None:
    df = pd.DataFrame({"a": pd.Series([], dtype=object)})
    raw, ploidy = _vectorized_parse_gt_columns(df, ["a"])
    assert raw.shape == ploidy.shape == (0, 1)
    raw, ploidy = _vectorized_parse_gt_columns(pd.DataFrame({"a": ["0/1"]}), [])
    assert raw.shape == ploidy.shape == (1, 0)
