"""Stage table, chunking-decision and parallel_safe=False reason logging (#220)."""

import argparse
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.runner import PipelineRunner
from variantcentrifuge.pipeline_core.stage import Stage
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.data_loading_stages import DataFrameLoadingStage


class _S(Stage):
    def __init__(self, name, parallel_safe=True, deps=()):
        self._n, self._ps, self._d = name, parallel_safe, set(deps)

    name = property(lambda self: self._n)
    description = property(lambda self: self._n)
    parallel_safe = property(lambda self: self._ps)
    dependencies = property(lambda self: self._d)

    def _process(self, context):
        sum(i * i for i in range(20000))
        return context


def _ctx(tmp_path, config=None):
    return PipelineContext(
        args=argparse.Namespace(), config=config or {}, workspace=Workspace(tmp_path, "t")
    )


@pytest.mark.unit
def test_stage_table_logged_and_written(tmp_path, caplog):
    diag = tmp_path / "diag"
    ctx = _ctx(tmp_path, {"diagnostics_output": str(diag)})
    with caplog.at_level(logging.INFO):
        PipelineRunner().run([_S("a"), _S("b", deps=["a"])], ctx)
    text = caplog.text
    assert "Wall(s)" in text and "CPU(s)" in text and "high-water" in text
    lines = (diag / "stage_summary.tsv").read_text().splitlines()
    assert lines[0] == "stage\twall_s\tcpu_s\tpeak_rss_mib"
    assert [ln.split("\t")[0] for ln in lines[1:]] == ["a", "b"]
    assert all(len(ln.split("\t")) == 4 and float(ln.split("\t")[3]) > 0 for ln in lines[1:])


@pytest.mark.unit
def test_stage_table_includes_process_pool_stages(tmp_path):
    diag = tmp_path / "diag"
    ctx = _ctx(tmp_path, {"diagnostics_output": str(diag)})
    runner = PipelineRunner(max_workers=2, executor_type="process")
    runner.run([_S("root"), _S("p1", deps=["root"]), _S("p2", deps=["root"])], ctx)
    rows = [ln.split("\t")[0] for ln in (diag / "stage_summary.tsv").read_text().splitlines()[1:]]
    assert sorted(rows) == ["p1", "p2", "root"]
    assert {"p1", "p2"} <= set(runner._execution_times)


@pytest.mark.unit
def test_stage_table_empty_is_noop(caplog):
    from variantcentrifuge.pipeline_core.runner_planner import log_stage_table

    with caplog.at_level(logging.INFO):
        log_stage_table({})
    assert "Wall(s)" not in caplog.text


@pytest.mark.unit
def test_stage_table_no_file_without_diagnostics(tmp_path, caplog):
    with caplog.at_level(logging.INFO):
        PipelineRunner().run([_S("a")], _ctx(tmp_path))
    assert "Wall(s)" in caplog.text
    assert not list(tmp_path.rglob("stage_summary.tsv"))


@pytest.mark.unit
def test_parallel_safe_false_reason_logged(tmp_path, caplog):
    stages = [_S("seq", parallel_safe=False), _S("par1"), _S("par2")]
    with caplog.at_level(logging.INFO):
        PipelineRunner().run(stages, _ctx(tmp_path))
    assert "Stage 'seq' has parallel_safe=False: running sequentially" in caplog.text


def _chunk_decision(tmp_path, size_mb, samples, caplog):
    f = tmp_path / "x.tsv"
    with open(f, "wb") as fh:
        fh.truncate(int(size_mb * 1024 * 1024))
    ctx = SimpleNamespace(config={}, vcf_samples=[f"s{i}" for i in range(samples)])
    with caplog.at_level(logging.INFO):
        return DataFrameLoadingStage()._should_use_chunks(ctx, Path(f))


@pytest.mark.unit
@pytest.mark.parametrize(("size_mb", "on"), [(1, False), (600, True)])
@pytest.mark.parametrize("samples", [10, 0])
def test_chunk_decision_logged(tmp_path, caplog, size_mb, on, samples):
    assert _chunk_decision(tmp_path, size_mb, samples, caplog) is on
    msg = next(r.message for r in caplog.records if "Chunked processing" in r.message)
    assert ("ON" if on else "OFF") in msg
    assert "file size" in msg and "estimated memory" in msg and "threshold" in msg
