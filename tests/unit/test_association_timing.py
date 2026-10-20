"""Tests for association sub-step timing hooks and the benchmark harness's --compare gate."""

import logging

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
    assert "IDENTICAL 2 rows x 3 columns" in capsys.readouterr().out


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
