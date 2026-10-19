"""Unit tests for parallel stage resource budgeting and worker heap/threads (#156)."""

from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.parallel_processing_stage import ParallelCompleteProcessingStage


@pytest.fixture
def mock_context(tmp_path: Path) -> PipelineContext:
    ws = Mock(spec=Workspace)
    ws.base_name = "test_sample"
    ws.intermediate_dir = tmp_path / "inter"
    ws.intermediate_dir.mkdir(parents=True, exist_ok=True)
    ws.get_intermediate_path = lambda name: ws.intermediate_dir / name

    config = {
        "threads": 8,
        "max_memory_gb": 64.0,
        "vcf_file": str(tmp_path / "input.vcf.gz"),
        "gzip_intermediates": False,
        "fields_to_extract": "CHROM POS REF ALT",
    }
    ctx = PipelineContext(args=Namespace(), config=config, workspace=ws)
    ctx.vcf_samples = ["S1", "S2"]
    return ctx


def test_resource_allocation_fresh_run(mock_context: PipelineContext) -> None:
    """Verify 8 chunks with 8 threads and 64GB RAM allocate 1 thread and 2GB heap per worker."""
    stage = ParallelCompleteProcessingStage()
    resources = stage._calculate_chunk_worker_resources(
        total_threads=8,
        num_chunks=8,
        shared_memory_gb=64.0,
    )
    assert resources["active_workers"] == 8
    assert resources["threads_per_worker"] == 1
    assert resources["snpsift_memory_gb"] == 2


def test_resource_allocation_memory_constrains_workers(mock_context: PipelineContext) -> None:
    """Verify memory constraint reduces pool concurrency when memory cannot fit 2GB per worker."""
    stage = ParallelCompleteProcessingStage()
    # 16 GB memory -> SnpSift total heap = 16 * 0.25 = 4 GB -> max 2 workers with 2GB each
    resources = stage._calculate_chunk_worker_resources(
        total_threads=8,
        num_chunks=8,
        shared_memory_gb=16.0,
    )
    assert resources["active_workers"] == 2
    assert resources["threads_per_worker"] == 4  # 8 threads // 2 workers
    assert resources["snpsift_memory_gb"] == 2


def test_resource_allocation_low_memory_single_worker(mock_context: PipelineContext) -> None:
    """Verify low memory clamps to 1 worker with minimum 2GB heap."""
    stage = ParallelCompleteProcessingStage()
    resources = stage._calculate_chunk_worker_resources(
        total_threads=4,
        num_chunks=4,
        shared_memory_gb=4.0,
    )
    assert resources["active_workers"] == 1
    assert resources["threads_per_worker"] == 4
    assert resources["snpsift_memory_gb"] == 2


def test_resource_allocation_partial_resume(mock_context: PipelineContext) -> None:
    """Verify only 2 missing chunks allocate 4 threads each with 8 threads total."""
    stage = ParallelCompleteProcessingStage()
    resources = stage._calculate_chunk_worker_resources(
        total_threads=8,
        num_chunks=2,
        shared_memory_gb=64.0,
    )
    assert resources["active_workers"] == 2
    assert resources["threads_per_worker"] == 4
    assert resources["snpsift_memory_gb"] == 8


def test_resource_allocation_single_chunk(mock_context: PipelineContext) -> None:
    """Verify 1 chunk gets all threads and full SnpSift heap budget."""
    stage = ParallelCompleteProcessingStage()
    resources = stage._calculate_chunk_worker_resources(
        total_threads=8,
        num_chunks=1,
        shared_memory_gb=64.0,
    )
    assert resources["active_workers"] == 1
    assert resources["threads_per_worker"] == 8
    assert resources["snpsift_memory_gb"] == 16


def test_no_unconditional_two_thread_floor() -> None:
    """Verify thread calculation allows 1 thread per worker without unconditional 2-thread floor."""
    stage = ParallelCompleteProcessingStage()
    resources = stage._calculate_chunk_worker_resources(
        total_threads=4,
        num_chunks=4,
        shared_memory_gb=64.0,
    )
    assert resources["active_workers"] == 4
    assert resources["threads_per_worker"] == 1  # Not 2!


def test_zero_missing_chunks_handled_cleanly(mock_context: PipelineContext, tmp_path: Path) -> None:
    """Verify 0 missing chunks (all exist) skips pool execution cleanly."""
    stage = ParallelCompleteProcessingStage()
    mock_context.config["vcf_file"] = str(tmp_path / "dummy.vcf")
    Path(mock_context.config["vcf_file"]).write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n")

    bed_chunks = [tmp_path / "c0.bed", tmp_path / "c1.bed"]
    for b in bed_chunks:
        b.write_text("chr1\t1\t10\n")

    # Mock chunk validation returning True for all chunks
    with (
        patch.object(stage, "_validate_chunk_tsv", return_value=True),
        patch.object(stage, "_process_single_chunk") as mock_chunk,
    ):
        result = stage._process_chunks_parallel(mock_context, bed_chunks)
        assert mock_chunk.call_count == 0
        assert len(result) == 2
