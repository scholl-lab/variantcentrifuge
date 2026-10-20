"""Run-level genotype-code matrix (#219): one parse, identical results for every GT layout."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scripts.bench_association import make_df, parse_args, run_stage
from variantcentrifuge.association.aggregation import (
    _aggregate_gene_burden_from_columns,
    _aggregate_gene_burden_from_gt,
)
from variantcentrifuge.association.gt_codes import aggregate_gene_counts, parse_genotype_codes
from variantcentrifuge.gt_reconstruction import _find_per_sample_gt_columns
from variantcentrifuge.stages.genotype_matrix_builder import _GenotypeMatrixBuilder


def _cells(gc):
    return [[gc.vocab[c] for c in row] for row in gc.codes]


@pytest.mark.parametrize("missing_as_ref", [False, True])
def test_stage_identical_with_and_without_per_sample_gt_columns(tmp_path, missing_as_ref):
    df, samples = make_df(6, 80, "columns", seed=3, miss=0.03)
    gen_cols = _find_per_sample_gt_columns(df)
    argv = "--seed 3 --case-frac 0.3 --unassigned 0.1 --tests fisher,logistic_burden,skat"
    args = parse_args(argv.split() + (["--missing-as-ref"] if missing_as_ref else []))
    layouts = {
        "packed_only": df.drop(columns=gen_cols),
        "packed_and_columns": df,
        "columns_only": df.drop(columns=["GT"]),
    }
    out = {}
    for name, frame in layouts.items():
        res, _, timings = run_stage(frame, samples, args, tmp_path)
        assert "gt_parse" in timings
        out[name] = res.to_csv(sep="\t", index=False)
    assert out["packed_only"] == out["packed_and_columns"] == out["columns_only"]


def test_packed_parse_named_positional_and_empty_rows():
    gt = pd.Series(
        [
            " A ( 0/1:12:3,9); B(./.:0) ;Z(1/1)",  # spaces, extra fields, unknown sample
            "0/1,./.,1|1",  # positional: every VCF sample in order
            "NA",
            None,
            "C(1)",
        ]
    )
    gc = parse_genotype_codes(pd.DataFrame({"GT": gt}), ["A", "C", "B"], ["A", "B", "C"], [])
    assert gc.codes.dtype == np.int8
    assert _cells(gc) == [
        ["0/1", "0/0", "./."],
        ["0/1", "1|1", "./."],
        ["0/0", "0/0", "0/0"],
        ["0/0", "0/0", "0/0"],
        ["0/0", "1", "0/0"],
    ]
    assert gc.labels == ["A", "C", "B"]


def test_positional_rows_shorter_than_samples_are_missing():
    # Transcript de-duplication leaves a non-contiguous index; codes follow row order.
    df = pd.DataFrame({"GT": ["0/1", "A(1/1)", "0/1,./.,./1"]}, index=[0, 2, 5])
    gc = parse_genotype_codes(df, ["A", "B", "C"], ["A", "B", "C"], [])
    assert _cells(gc) == [["0/1", "./.", "./."], ["1/1", "0/0", "0/0"], ["0/1", "./.", "./1"]]


def test_column_parse_categorical_object_and_absent_columns(caplog):
    df = pd.DataFrame(
        {
            "GEN_0__GT": pd.Categorical(["0/1", None, "1/1"]),
            "GEN_1__GT": ["0|1", np.nan, 1],
        }
    )
    gc = parse_genotype_codes(df, ["s0", "s1", "s2"], ["s0", "s1", "s2"], list(df))
    assert _cells(gc) == [["0/1", "0|1", ""], ["", "", ""], ["1/1", "1", ""]]
    assert gc.labels == ["GEN_0__GT", "GEN_1__GT", "s2"]
    assert "GT column count (2) != VCF sample count (3)" in caplog.text


def test_codes_widen_to_int16_for_large_vocabularies():
    gt = pd.Series([f"A({i}/0:x)" for i in range(200)])
    gc = parse_genotype_codes(pd.DataFrame({"GT": gt}), ["A"], ["A"], [])
    assert gc.codes.dtype == np.int16
    assert gc.vocab[gc.codes[150, 0]] == "150/0"


def test_column_codes_widen_to_int16():
    df = pd.DataFrame({f"GEN_{j}__GT": [f"{j}/0"] for j in range(130)})
    samples = [f"s{j}" for j in range(130)]
    gc = parse_genotype_codes(df, samples, samples, list(df))
    assert gc.codes.dtype == np.int16
    assert gc.vocab[gc.codes[0, 129]] == "129/0"


# Inputs on which the base column and packed collapsing paths agree.
_AGREEING = [
    ["0/1", "0/0", "1/1", "./.", "0|1", "1/2", "0/0", "1"],
    ["0/0", "0/0", "1/1", "0/0", "1|1", "0/0", "./1", "0"],
    ["0/1", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "."],
]


@pytest.mark.parametrize("missing_as_ref", [False, True])
def test_counts_match_base_column_and_packed_aggregation(missing_as_ref):
    samples = [f"S{i}" for i in range(8)]
    cases, controls = set(samples[:4]), set(samples[4:])
    cols = {f"GEN_{j}__GT": [r[j] for r in _AGREEING] for j in range(8)}
    df = pd.DataFrame({"GENE": ["G1", "G1", "G2"], **cols})
    packed = [
        ";".join(f"{s}({g})" for s, g in zip(samples, r, strict=True) if g != "0/0")
        for r in _AGREEING
    ]
    gc = parse_genotype_codes(df, samples, samples, list(cols))
    rows = df.groupby("GENE").indices
    is_case = np.array([s in cases for s in samples])
    got = aggregate_gene_counts(gc, rows, is_case, ~is_case, missing_as_ref)
    by_col = _aggregate_gene_burden_from_columns(
        df, cases, controls, samples, list(cols), missing_as_ref
    )
    by_gt = _aggregate_gene_burden_from_gt(df.assign(GT=packed), cases, controls, missing_as_ref)
    assert got == by_col == by_gt
    packed_gc = parse_genotype_codes(df.assign(GT=packed), samples, samples, [])
    assert aggregate_gene_counts(packed_gc, rows, is_case, ~is_case, missing_as_ref) == got


def _builder(gc, gene_df, rows, samples):
    return _GenotypeMatrixBuilder(
        gene_df=gene_df,
        vcf_samples=samples,
        gt_columns=[],
        is_binary=True,
        missing_site_threshold=0.5,
        missing_sample_threshold=0.8,
        phenotype_vector=np.array([1.0, 0.0, 0.0]),
        covariate_matrix=None,
        genotypes=gc,
        rows=rows,
    )


def test_builder_from_codes_matches_column_builder():
    samples = ["a", "b", "c"]
    df = pd.DataFrame(
        {
            "GENE": ["G"] * 3,
            "EFFECT": ["x", "y", "z"],
            "GEN_0__GT": ["0/1", "1/1", "./."],
            "GEN_1__GT": ["0/0", "0/1", "0/1"],
            "GEN_2__GT": ["1/1", "0/0", "0/1"],
        }
    )
    gt_cols = _find_per_sample_gt_columns(df)
    gc = parse_genotype_codes(df, samples, samples, gt_cols)
    rows = np.arange(3)
    from_codes = _builder(gc, df[["GENE", "EFFECT"]], rows, samples)
    from_cols = _builder(None, df, None, samples)
    from_cols.gt_columns = gt_cols
    a, b = from_codes(), from_cols()
    assert a.keys() == b.keys()
    for key in ("genotype_matrix", "variant_mafs", "variant_keep_mask"):
        np.testing.assert_array_equal(a[key], b[key])


def test_builder_from_codes_rejects_multiallelic_with_sample_label():
    df = pd.DataFrame({"GENE": ["G"], "CHROM": ["1"], "POS": [5], "GT": ["b(1/2);c(0/1)"]})
    gc = parse_genotype_codes(df, ["a", "b", "c"], ["a", "b", "c"], [])
    with pytest.raises(ValueError, match=r"'1/2' \(gene G, variant 1:5, column 'b'\)"):
        _builder(gc, df, np.arange(1), ["a", "b", "c"])()


_HASHSEED_SCRIPT = r"""
import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock

from scripts.bench_association import make_df
from variantcentrifuge.gt_reconstruction import _find_per_sample_gt_columns
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.association_analysis_stages import AssociationAnalysisStage

tmp = Path(sys.argv[1])
df, samples = make_df(4, 60, "columns", seed=5, miss=0.02)
layouts = {"packed": df.drop(columns=_find_per_sample_gt_columns(df)),
           "columns": df.drop(columns=["GT"]), "both": df}
for name, frame in layouts.items():
    ws = Mock(spec=Workspace)
    ws.output_dir = ws.intermediate_dir = tmp
    cfg = {
        "perform_association": True,
        # "ABSENT" is a listed case that is not in the VCF
        "case_samples": samples[:15] + ["ABSENT"],
        "control_samples": samples[15:],
        "association_tests": ["fisher", "logistic_burden"], "missing_as_ref": True,
        "output_dir": str(tmp), "output_file_base": name, "gzip_intermediates": False,
    }
    ctx = PipelineContext(args=Namespace(), config=cfg, workspace=ws)
    ctx.current_dataframe, ctx.vcf_samples = frame, samples
    res = AssociationAnalysisStage()._process(ctx).association_results
    print(name, res.to_csv(sep="\t", index=False).replace("\n", "|"))
"""


def test_absent_case_sample_counts_are_deterministic_and_layout_independent(tmp_path):
    """Only genotyped samples count, whatever the hash seed or GT layout (#219)."""
    import os
    import subprocess
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[2]
    outputs = []
    for seed in ("0", "1"):
        env = {**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": str(root)}
        proc = subprocess.run(
            [sys.executable, "-c", _HASHSEED_SCRIPT, str(tmp_path)],
            capture_output=True,
            text=True,
            env=env,
            cwd=root,
            check=True,
        )
        tables = dict(line.split(" ", 1) for line in proc.stdout.splitlines())
        assert tables["packed"] == tables["columns"] == tables["both"]
        outputs.append(tables["packed"])
    assert outputs[0] == outputs[1]
    header, *rows = outputs[0].strip("|").split("|")
    n_cases = header.split("\t").index("n_cases")
    assert {r.split("\t")[n_cases] for r in rows} == {"15"}


def test_packed_parse_categorical_gt_and_unknown_samples_only():
    gt = pd.Series(pd.Categorical(["X(0/1)", None, "A(1/1)"]))
    gc = parse_genotype_codes(pd.DataFrame({"GT": gt}), ["A"], ["A", "X"], [])
    assert _cells(gc) == [["0/0"], ["0/0"], ["1/1"]]
    only_unknown = pd.DataFrame({"GT": ["X(0/1)"]})
    assert _cells(parse_genotype_codes(only_unknown, ["A"], ["A", "X"], [])) == [["0/0"]]


def test_stage_without_genotypes_uses_legacy_counts(tmp_path, caplog):
    """No packed GT and no per-sample columns: pre-computed per-variant counts (legacy)."""
    from argparse import Namespace
    from unittest.mock import Mock

    from variantcentrifuge.pipeline_core.context import PipelineContext
    from variantcentrifuge.pipeline_core.workspace import Workspace
    from variantcentrifuge.stages.association_analysis_stages import AssociationAnalysisStage

    df = pd.DataFrame(
        {
            "GENE": ["G1", "G1"],
            "proband_count": [20, 20],
            "control_count": [40, 40],
            "proband_variant_count": [3, 1],
            "control_variant_count": [1, 0],
            "proband_allele_count": [3, 1],
            "control_allele_count": [1, 0],
        }
    )
    ws = Mock(spec=Workspace)
    ws.output_dir = ws.intermediate_dir = tmp_path
    cfg = {
        "perform_association": True,
        "case_samples": [f"c{i}" for i in range(20)],
        "control_samples": [f"n{i}" for i in range(40)],
        "association_tests": ["fisher"],
        "output_dir": str(tmp_path),
        "gzip_intermediates": False,
    }
    ctx = PipelineContext(args=Namespace(), config=cfg, workspace=ws)
    ctx.current_dataframe = df
    with caplog.at_level("INFO"):
        out = AssociationAnalysisStage()._process(ctx)
    assert "pre-computed per-variant counts (legacy mode)" in caplog.text
    assert "gt_parse" not in out.stage_results["association_analysis"]["substep_seconds"]
