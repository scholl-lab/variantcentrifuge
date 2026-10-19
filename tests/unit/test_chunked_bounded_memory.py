"""Tests verifying bounded chunk memory in ChunkedAnalysisStage."""

from __future__ import annotations

import threading
from argparse import Namespace
from pathlib import Path

from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.chunked_analysis_stages import ChunkedAnalysisStage


def _create_mock_tsv(path: Path, num_genes: int = 25) -> None:
    lines = ["CHROM\tPOS\tREF\tALT\tGENE\tGT\n"]
    for i in range(num_genes):
        gene = f"GENE_{i:02d}"
        lines.append(f"1\t{i * 100 + 1}\tA\tT\t{gene}\tS1(0/1)\n")
        lines.append(f"1\t{i * 100 + 2}\tC\tG\t{gene}\tS1(1/1)\n")
    path.write_text("".join(lines))


def test_chunked_sequential_bounded_memory(tmp_path: Path):
    """Verify that sequential chunk processing maintains at most 1 active chunk in memory."""
    tsv_file = tmp_path / "input.tsv"
    _create_mock_tsv(tsv_file, num_genes=25)

    workspace = Workspace(tmp_path / "work", "test")
    config = {
        "use_chunked_processing": True,
        "threads": 1,
        "parallel_chunks": False,
        "chunk_size": 2,
        "gzip_intermediates": False,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.data = tsv_file
    context.vcf_samples = ["S1"]

    stage = ChunkedAnalysisStage()

    active_chunks = 0
    max_active_observed = 0

    orig_process_single = stage._process_single_chunk

    def instrumented_process(chunk_df, ctx, idx):
        nonlocal active_chunks, max_active_observed
        active_chunks += 1
        if active_chunks > max_active_observed:
            max_active_observed = active_chunks
        try:
            return orig_process_single(chunk_df, ctx, idx)
        finally:
            active_chunks -= 1

    stage._process_single_chunk = instrumented_process  # type: ignore[method-assign]
    stage._process(context)

    # In sequential mode, at most 1 chunk is processed at a time
    assert max_active_observed == 1
    assert context.current_dataframe is not None
    assert len(context.current_dataframe) == 50


def test_chunked_parallel_bounded_memory(tmp_path: Path):
    """Verify that parallel chunk processing caps active/pending chunks across 20+ chunks."""
    tsv_file = tmp_path / "input.tsv"
    _create_mock_tsv(tsv_file, num_genes=25)

    workspace = Workspace(tmp_path / "work", "test")
    max_workers = 2
    config = {
        "use_chunked_processing": True,
        "threads": max_workers,
        "parallel_chunks": True,
        "chunk_size": 2,
        "gzip_intermediates": False,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.data = tsv_file
    context.vcf_samples = ["S1"]

    stage = ChunkedAnalysisStage()

    lock = threading.Lock()
    active_chunks = 0
    max_active_observed = 0

    orig_process_single = stage._process_single_chunk

    def instrumented_process(chunk_df, ctx, idx):
        nonlocal active_chunks, max_active_observed
        with lock:
            active_chunks += 1
            if active_chunks > max_active_observed:
                max_active_observed = active_chunks
        try:
            return orig_process_single(chunk_df, ctx, idx)
        finally:
            with lock:
                active_chunks -= 1

    stage._process_single_chunk = instrumented_process  # type: ignore[method-assign]
    stage._process(context)

    # In parallel mode with max_workers=2, active workers <= max_workers <= max_workers * 2
    assert max_active_observed <= max_workers
    assert max_active_observed <= max_workers * 2
    assert context.current_dataframe is not None
    assert len(context.current_dataframe) == 50


def test_chunked_analysis_quoting_json_roundtrip(tmp_path: Path):
    """Verify JSON strings and special characters survive chunk write and reload."""
    tsv_file = tmp_path / "json_input.tsv"
    tsv_file.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\tGT\tInheritance_Details\n"
        '1\t100\tA\tT\tGENE1\tS1(0/1)\t{"pattern": "AD", "desc": "simple"}\n'
        '1\t200\tC\tG\tGENE1\tS1(1/1)\t{"pattern": "AR", "nested": {"key": "val"}}\n'
    )

    workspace = Workspace(tmp_path / "work", "test")
    config = {
        "use_chunked_processing": True,
        "threads": 1,
        "parallel_chunks": False,
        "chunk_size": 10,
        "gzip_intermediates": False,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.data = tsv_file
    context.vcf_samples = ["S1"]

    stage = ChunkedAnalysisStage()
    stage._process(context)

    assert context.current_dataframe is not None
    df = context.current_dataframe
    assert len(df) == 2
    assert df["Inheritance_Details"].iloc[0] == '{"pattern": "AD", "desc": "simple"}'
    assert df["Inheritance_Details"].iloc[1] == '{"pattern": "AR", "nested": {"key": "val"}}'


def test_chunked_parallel_out_of_order_backlog_bounded(tmp_path: Path):
    """Verify that slow chunk 0 caps in-flight plus buffered results to max_pending."""
    import time
    from unittest.mock import patch

    from variantcentrifuge.memory import ResourceManager

    tsv_file = tmp_path / "input.tsv"
    _create_mock_tsv(tsv_file, num_genes=20)

    workspace = Workspace(tmp_path / "work", "test")
    max_workers = 2
    max_pending = max(max_workers * 2, 4)
    config = {
        "use_chunked_processing": True,
        "threads": max_workers,
        "parallel_chunks": True,
        "chunk_size": 2,
        "gzip_intermediates": False,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.data = tsv_file
    context.vcf_samples = ["S1"]

    stage = ChunkedAnalysisStage()

    submitted_indices: list[int] = []
    lock = threading.Lock()
    orig_process_single = stage._process_single_chunk

    chunk0_started = threading.Event()
    block_chunk0 = threading.Event()

    def delayed_process(chunk_df, ctx, idx):
        with lock:
            submitted_indices.append(idx)
        if idx == 0:
            chunk0_started.set()
            block_chunk0.wait(timeout=5.0)
        return orig_process_single(chunk_df, ctx, idx)

    stage._process_single_chunk = delayed_process  # type: ignore[method-assign]

    with patch.object(ResourceManager, "auto_chunk_size", return_value=2):
        t = threading.Thread(target=stage._process, args=(context,))
        t.start()

        assert chunk0_started.wait(timeout=5.0)
        time.sleep(0.05)

        with lock:
            count = len(submitted_indices)
            assert count <= max_pending + max_workers
            assert count < 20

        block_chunk0.set()
        t.join(timeout=10.0)
        assert not t.is_alive()

    # When chunk 0 was sleeping, chunk submissions must pause once backlog reaches max_pending
    assert context.current_dataframe is not None
    assert len(context.current_dataframe) == 40
