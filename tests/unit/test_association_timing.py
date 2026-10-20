"""Tests for association sub-step timing hooks and the benchmark harness's --compare gate."""

import json
import logging

import numpy as np
import pandas as pd
import pytest

from scripts.bench_association import compare, main, make_df, parse_args, run_stage
from variantcentrifuge.association.timing import substep


def test_substep_accumulates_and_logs(caplog):
    timings: dict[str, float] = {}
    with caplog.at_level(logging.DEBUG, logger="variantcentrifuge.association.timing"):
        with substep(timings, "a"):
            pass
        first = timings["a"]
        with substep(timings, "a"):
            pass
    assert timings["a"] >= first >= 0.0
    assert list(timings) == ["a"]
    assert "Association sub-step a" in caplog.text


def test_substep_records_time_when_body_raises():
    timings: dict[str, float] = {}
    with pytest.raises(ValueError), substep(timings, "boom"):
        raise ValueError("x")
    assert "boom" in timings


@pytest.mark.parametrize("layout", ["packed", "columns"])
def test_stage_exposes_substep_seconds(tmp_path, layout):
    df, samples = make_df(3, 40, layout, seed=1)
    argv = "--seed 1 --case-frac 0.3 --unassigned 0.1 --tests fisher,logistic_burden"
    args = parse_args(argv.split())
    res, wall, timings = run_stage(df, samples, args, tmp_path)
    assert len(res) == 3
    expected = {"aggregation", "matrix_builders", "engine", "write_output"}
    assert expected <= set(timings)
    assert ("gt_unpack" in timings) == (layout == "packed")
    assert sum(timings.values()) <= wall


def test_target_workload_presets_and_covariates(tmp_path):
    args = parse_args(["--size", "B", "--genes", "3", "--samples", "60"])
    assert (args.cases, args.controls, args.covariates) == (790, 100, 14)
    assert args.missing_as_ref and args.artifact and args.tests == "logistic_burden,skat"
    assert parse_args(["--size", "A", "--cases", "7"]).cases == 7  # explicit flag wins
    df, samples = make_df(3, 60, "packed", seed=2)
    args.cases, args.controls, args.covariates, args.artifact = 20, 30, 4, False
    res, _, timings = run_stage(df, samples, args, tmp_path)
    assert len(res) == 3 and res["n_cases"].max() <= 20
    assert (tmp_path / "covariates.tsv").read_text().startswith("sample\tage\tsex\tPC1\tPC2\n")
    assert "covariates_phenotype" in timings


def _write(path, rows):
    path.write_text("gene\tp\tnote\n" + "".join(f"{g}\t{p}\t{n}\n" for g, p, n in rows))
    return str(path)


def test_compare_identical(tmp_path, capsys):
    rows = [("G1", "0.1", ""), ("G2", "nan", "LOW")]
    a, b = _write(tmp_path / "a.tsv", rows), _write(tmp_path / "b.tsv", rows)
    assert compare(a, b) == 0
    assert main(["--compare", a, b]) == 0
    out = capsys.readouterr().out
    assert "IDENTICAL" in out and "NOTE no frame pickles" in out


def test_compare_reports_float_and_string_diffs(tmp_path, capsys):
    a = _write(tmp_path / "a.tsv", [("G1", "0.1", ""), ("G2", "0.5", "LOW")])
    b = _write(tmp_path / "b.tsv", [("G1", "0.10000000000000002", ""), ("G2", "0.5", "")])
    assert compare(a, b) == 1
    out = capsys.readouterr().out
    assert "DIFF p: 1 rows, max_abs=1.39e-17" in out
    assert "DIFF note: 1 rows, 1 non-numeric/NaN mismatches" in out


def test_compare_rejects_shape_or_column_order(tmp_path, capsys):
    a = _write(tmp_path / "a.tsv", [("G1", "0.1", "")])
    b = tmp_path / "b.tsv"
    b.write_text("p\tgene\tnote\n0.1\tG1\t\n")
    assert compare(a, str(b)) == 1
    assert "DIFF shape/columns" in capsys.readouterr().out


def _frame_pair(tmp_path, fa, fb):
    """Write two frames as --out does (TSV + pickle) and return their TSV paths."""
    paths = []
    for name, frame in (("a", fa), ("b", fb)):
        out = str(tmp_path / f"{name}.tsv")
        frame.to_csv(out, sep="\t", index=False)
        frame.to_pickle(out + ".pkl")
        paths.append(out)
    return paths


def test_compare_pickle_catches_dtype_and_none_vs_nan(tmp_path, capsys):
    f64 = pd.DataFrame({"gene": ["G1"], "p": np.array([0.5], dtype="float64")})
    a, b = _frame_pair(tmp_path, f64, f64.astype({"p": "float32"}))
    assert compare(a, b) == 1  # TSV text is identical ("0.5"); only the pickle differs
    assert "DIFF frame" in capsys.readouterr().out
    a, b = _frame_pair(
        tmp_path, pd.DataFrame({"w": [None]}, dtype=object), pd.DataFrame({"w": [np.nan]})
    )
    assert compare(a, b) == 1
    assert "DIFF frame" in capsys.readouterr().out
    a, b = _frame_pair(tmp_path, f64, f64)
    assert compare(a, b) == 0


def test_compare_requires_pickle_on_both_sides(tmp_path, capsys):
    a, b = _frame_pair(tmp_path, pd.DataFrame({"x": [1]}), pd.DataFrame({"x": [1]}))
    (tmp_path / "b.tsv.pkl").unlink()
    assert compare(a, b) == 1
    assert "present only for A" in capsys.readouterr().out


def _artifact(out, arrays, meta=None):
    np.savez_compressed(out + ".scores.npz", **arrays)
    with open(out + ".scores.npz.meta.json", "w") as fh:
        json.dump(meta or {"n_cases": 5, "provenance": {"snpeff_version": None}}, fh)


@pytest.mark.parametrize(
    ("arrays_b", "meta_b", "message"),
    [
        ({"G1::u": np.array([1.0, np.nan]), "G1::ac": np.array([3])}, None, None),
        ({"G1::u": np.array([1.5, np.nan]), "G1::ac": np.array([3])}, None, "artifact array G1::u"),
        (
            {"G1::u": np.array([1.0, np.nan], "f4"), "G1::ac": np.array([3])},
            None,
            "artifact array G1::u: float64",
        ),
        ({"G1::u": np.array([1.0, np.nan])}, None, "artifact keys"),
        (
            {"G1::u": np.array([1.0, np.nan]), "G1::ac": np.array([3]), "G2::u": np.ones(1)},
            None,
            "artifact keys",
        ),
        (
            {"G1::u": np.array([1.0, np.nan]), "G1::ac": np.array([3])},
            {"n_cases": 6, "provenance": {"snpeff_version": None}},
            "artifact metadata: ['n_cases']",
        ),
    ],
)
def test_compare_score_artifacts(tmp_path, capsys, arrays_b, meta_b, message):
    a, b = _frame_pair(tmp_path, pd.DataFrame({"x": [1]}), pd.DataFrame({"x": [1]}))
    _artifact(a, {"G1::u": np.array([1.0, np.nan]), "G1::ac": np.array([3])})
    _artifact(b, arrays_b, meta_b)
    rc = compare(a, b)
    out = capsys.readouterr().out
    if message is None:
        assert rc == 0, out
    else:
        assert rc == 1 and f"DIFF {message}" in out, out


def test_out_writes_pickle_and_artifact_that_compare_equal(tmp_path, capsys):
    argv = [
        "--genes",
        "3",
        "--samples",
        "200",
        "--tests",
        "fisher,logistic_burden",
        "--artifact",
        "--out",
    ]
    a, b = str(tmp_path / "a.tsv"), str(tmp_path / "b.tsv")
    assert main([*argv, a]) == 0 and main([*argv, b]) == 0
    for suffix in (".pkl", ".scores.npz", ".scores.npz.meta.json"):
        assert (tmp_path / f"a.tsv{suffix}").exists()
    capsys.readouterr()
    assert compare(a, b) == 0
