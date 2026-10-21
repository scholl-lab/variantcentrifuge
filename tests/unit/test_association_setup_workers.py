"""Association set-up (GT parse, counts) follows the engine's worker decision (#212).

Stage results at workers > 1 equal workers = 1 for quantitative traits and COAST + gene sets.
"""

from __future__ import annotations

from argparse import Namespace
from unittest.mock import Mock

import numpy as np
import pyarrow as pa
import pytest

import variantcentrifuge.association.gt_codes as gt_codes
from scripts.bench_association import make_df
from variantcentrifuge.association.engine_workers import fork_available
from variantcentrifuge.association.tests.fisher import FisherExactTest
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.association_analysis_stages import AssociationAnalysisStage

needs_fork = pytest.mark.skipif(not fork_available(), reason="fork_map needs fork (Linux)")


def _stage(tmp_path, n_genes, workers, **config):
    df, samples = make_df(n_genes, 60, "packed", seed=7, vmin=2, vmax=6, miss=0.02)
    df["IMPACT"] = np.where(df["EFFECT"] == "stop_gained", "HIGH", "MODERATE")
    # COAST (SIFT/PolyPhen model) annotation columns
    df["ANN_0__EFFECT"], df["ANN_0__IMPACT"] = df["EFFECT"], df["IMPACT"]
    df["dbNSFP_SIFT_pred"] = np.where(df.index % 2 == 0, "D", "T")
    df["dbNSFP_Polyphen2_HDIV_pred"] = np.where(df.index % 3 == 0, "D", "B")
    ws = Mock(spec=Workspace)
    ws.output_dir = ws.intermediate_dir = tmp_path
    cfg = {
        "perform_association": True,
        "case_samples": samples[:25],
        "control_samples": samples[25:],
        "association_tests": ["fisher"],
        "association_min_cases": 1,
        "output_dir": str(tmp_path),
        "output_file_base": f"w{workers}",
        "gzip_intermediates": False,
        "association_workers": workers,
        **config,
    }
    ctx = PipelineContext(args=Namespace(), config=cfg, workspace=ws)
    ctx.current_dataframe = df
    ctx.vcf_samples = samples
    rng = np.random.default_rng(3)
    values = rng.normal(size=60)
    ctx.phenotype_data = {s: {f"{v:.6f}"} for s, v in zip(samples, values, strict=True)}
    AssociationAnalysisStage()._process(ctx)
    return ctx.association_results.to_csv(sep="\t", index=False)


def _forbid_fork(monkeypatch):
    def boom(*_a, **_k):
        raise AssertionError("set-up forked although the engine runs single-process")

    monkeypatch.setattr(gt_codes, "fork_map", boom)


@needs_fork
def test_setup_sequential_below_engine_parallel_threshold(tmp_path, monkeypatch):
    _forbid_fork(monkeypatch)
    _stage(tmp_path, 20, 4)  # < MIN_GENES_FOR_PARALLEL: the engine runs sequentially


@needs_fork
def test_setup_sequential_when_a_test_is_not_parallel_safe(tmp_path, monkeypatch):
    monkeypatch.setattr(FisherExactTest, "parallel_safe", False)  # like R SKAT / R COAST
    _forbid_fork(monkeypatch)
    _stage(tmp_path, 60, 4)


@needs_fork
def test_setup_forks_when_engine_does(tmp_path, monkeypatch):
    calls = []
    real = gt_codes.fork_map

    def spy(*a, **k):
        calls.append(a[0].__name__)
        yield from real(*a, **k)

    monkeypatch.setattr(gt_codes, "fork_map", spy)
    assert _stage(tmp_path / "w4", 60, 4) == _stage(tmp_path / "w1", 60, 1)
    assert calls == ["_parse_packed_rows", "_count_genes"]


@needs_fork
def test_quantitative_identical_across_workers(tmp_path):
    cfg = {
        "trait_type": "quantitative",
        "phenotype_value_column": "y",
        "association_tests": ["linear_burden", "skat"],
    }
    assert _stage(tmp_path / "w3", 60, 3, **cfg) == _stage(tmp_path / "w1", 60, 1, **cfg)


@needs_fork
def test_coast_and_gene_sets_identical_across_workers(tmp_path):
    sets = tmp_path / "sets.tsv"
    sets.write_text(
        "set_id\tgene\n" + "".join(f"SET{i // 10}\tG{i:05d}\n" for i in range(0, 60, 3))
    )
    cfg = {"association_tests": ["fisher", "coast"], "gene_set_file": str(sets)}
    w3 = _stage(tmp_path / "w3", 60, 3, **cfg)
    assert w3 == _stage(tmp_path / "w1", 60, 1, **cfg)
    assert "SET1" in w3 and "coast_pvalue" in w3.splitlines()[0]


def test_parallel_task_keeps_batch_char_cap(monkeypatch):
    seen = []
    real = gt_codes.pa.array

    def spy(values, *a, **k):
        if k.get("type") == pa.large_string() and len(values) and isinstance(values[0], str):
            seen.append(sum(map(len, values)))
        return real(values, *a, **k)

    text = np.array([f"A({i % 2}/1);B(1/1)" for i in range(50)], dtype=object)
    state = (text, real(["A", "B"], type=pa.large_string()))
    monkeypatch.setattr(gt_codes, "_BATCH_CHARS", 40)  # ~3 rows per batch
    monkeypatch.setattr(gt_codes.pa, "array", spy)
    named, rows, cols, gts, idx = gt_codes._parse_packed_rows(state, 0, 50)
    assert max(seen) <= 40 + max(map(len, text))
    assert named.all() and len(rows) == len(cols) == len(idx) == 100
    assert [gts[i] for i in idx[:4]] == ["0/1", "1/1", "1/1", "1/1"]


@needs_fork
@pytest.mark.parametrize("workers", [0, 1, -1, 2])
@pytest.mark.parametrize("n_genes", [49, 50])
def test_setup_and_engine_take_the_same_decision(tmp_path, monkeypatch, workers, n_genes):
    import variantcentrifuge.stages.association_analysis_stages as stage_mod
    from variantcentrifuge.association.engine import AssociationEngine

    decisions, setup = [], []
    real_decide, real_parse = AssociationEngine.runs_parallel, stage_mod.parse_genotype_codes

    def decide(self, n_units):
        decisions.append((n_units, real_decide(self, n_units)))
        return decisions[-1][1]

    def parse(*a):
        setup.append(a[4])
        return real_parse(*a)

    monkeypatch.setattr(AssociationEngine, "runs_parallel", decide)
    monkeypatch.setattr(stage_mod, "parse_genotype_codes", parse)
    _stage(tmp_path, n_genes, workers)
    # first call: the stage's gene count (then gene sets if needed); last: the engine's loop
    (n_stage, stage_says), (n_engine, engine_says) = decisions[0], decisions[-1]
    assert n_stage == n_engine == n_genes
    assert stage_says == engine_says == (n_genes >= 50 and workers not in (0, 1))
    assert (setup[0] > 1) == engine_says
