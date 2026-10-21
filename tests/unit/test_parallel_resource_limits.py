"""Unit tests for parallel stage resource budgeting and worker heap/threads (#156)."""

from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from variantcentrifuge.filters import _snpsift_memory_flag
from variantcentrifuge.memory.resource_manager import ResourceManager
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
    """Verify 8 chunks with 8 threads and 64GB RAM allocate 1 thread and 4GB heap per worker."""
    stage = ParallelCompleteProcessingStage()
    resources = stage._calculate_chunk_worker_resources(
        total_threads=8,
        num_chunks=8,
        shared_memory_gb=64.0,
    )
    assert resources["active_workers"] == 8
    assert resources["threads_per_worker"] == 1
    assert resources["snpsift_memory_gb"] == 4


def test_resource_allocation_memory_constrains_workers(mock_context: PipelineContext) -> None:
    """Verify memory constraint reduces pool concurrency when memory cannot fit 2GB per worker."""
    stage = ParallelCompleteProcessingStage()
    # 8 GB memory -> SnpSift pool heap = 8 * 0.5 = 4 GB -> max 2 workers with 2GB each
    resources = stage._calculate_chunk_worker_resources(
        total_threads=8,
        num_chunks=8,
        shared_memory_gb=8.0,
    )
    assert resources["active_workers"] == 2
    assert resources["threads_per_worker"] == 4  # 8 threads // 2 workers
    assert resources["snpsift_memory_gb"] == 2


def test_resource_allocation_16_workers_64gb() -> None:
    """Verify 16 requested workers under 64-GiB budget allocate 16 workers with 2GB heap each."""
    stage = ParallelCompleteProcessingStage()
    resources = stage._calculate_chunk_worker_resources(
        total_threads=16,
        num_chunks=16,
        shared_memory_gb=64.0,
    )
    assert resources["active_workers"] == 16
    assert resources["threads_per_worker"] == 1
    assert resources["snpsift_memory_gb"] == 2


def test_resource_allocation_16_workers_8gb() -> None:
    """Verify 16 requested workers under 8-GiB budget are bounded to 2 workers with 8 threads."""
    stage = ParallelCompleteProcessingStage()
    resources = stage._calculate_chunk_worker_resources(
        total_threads=16,
        num_chunks=16,
        shared_memory_gb=8.0,
    )
    assert resources["active_workers"] == 2
    assert resources["threads_per_worker"] == 8
    assert resources["snpsift_memory_gb"] == 2


def test_resource_allocation_slurm_mem_override(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify SLURM_MEM_PER_NODE environment variable is honored when budgeting worker pool."""
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "8192")  # 8192 MB = 8 GB
    rm = ResourceManager(config={})
    assert rm.memory_gb == pytest.approx(8.0, rel=1e-2)

    stage = ParallelCompleteProcessingStage()
    resources = stage._calculate_chunk_worker_resources(
        total_threads=16,
        num_chunks=16,
        shared_memory_gb=rm.memory_gb,
    )
    assert resources["active_workers"] == 2
    assert resources["threads_per_worker"] == 8
    assert resources["snpsift_memory_gb"] == 2


def test_snpsift_memory_flag_explicit_and_fallback() -> None:
    """Verify _snpsift_memory_flag respects explicit snpsift_memory_gb config."""
    assert _snpsift_memory_flag({"snpsift_memory_gb": 8}) == "-Xmx8g"
    assert _snpsift_memory_flag({"snpsift_memory_gb": 1}) == "-Xmx1g"
    assert _snpsift_memory_flag({"snpsift_memory_gb": 0}) == "-Xmx1g"
    assert _snpsift_memory_flag({"snpsift_memory_gb": 32}) == "-Xmx16g"  # Bounded to 16GB
    # Fallback to system detection when None
    flag = _snpsift_memory_flag({})
    assert flag.startswith("-Xmx") and flag.endswith("g")


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
    assert resources["snpsift_memory_gb"] == 16


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


def test_worker_command_thread_and_memory_dispatch(tmp_path: Path) -> None:
    """Verify chunk processing passes calculated threads and heap to underlying tools."""
    stage = ParallelCompleteProcessingStage()
    chunk_bed = tmp_path / "chunk_0.bed"
    chunk_bed.write_text("chr1\t1\t100\n")

    worker_cfg = {
        "threads_per_chunk": 4,
        "max_memory_gb": 8.0,
        "snpsift_memory_gb": 4,
        "bcftools_prefilter": 'FILTER="PASS"',
        "filter": '(ANN[*].EFFECT has "missense_variant")',
        "fields_to_extract": "CHROM POS REF ALT",
        "gzip_intermediates": False,
    }

    pp_mod = "variantcentrifuge.stages.parallel_processing_stage"
    with (
        patch(f"{pp_mod}.extract_variants") as mock_extract,
        patch(f"{pp_mod}.apply_snpsift_filter") as mock_filter,
        patch(f"{pp_mod}.extract_fields_bcftools") as mock_fields,
    ):
        stage._process_single_chunk(
            chunk_index=0,
            chunk_bed=chunk_bed,
            vcf_file=str(tmp_path / "in.vcf.gz"),
            base_name="sample",
            intermediate_dir=tmp_path,
            config=worker_cfg,
        )

        mock_extract.assert_called_once()
        extract_cfg = mock_extract.call_args[1]["cfg"]
        assert extract_cfg["threads"] == 4
        assert extract_cfg["bcftools_prefilter"] == 'FILTER="PASS"'

        mock_filter.assert_called_once()
        filter_cfg = mock_filter.call_args[0][2]
        assert filter_cfg["threads"] == 4
        assert filter_cfg["snpsift_memory_gb"] == 4

        mock_fields.assert_called_once()
        fields_cfg = mock_fields.call_args[1]["cfg"]
        assert fields_cfg["threads"] == 4


@pytest.mark.parametrize("sub_4gb", [0.5, 1.0, 2.0, 3.9])
def test_sub_4gb_memory_budget_rejected(sub_4gb: float) -> None:
    """Verify memory budgets under 4 GiB are rejected to protect 2 GiB heap + headroom."""
    stage = ParallelCompleteProcessingStage()
    with pytest.raises(ValueError, match="too small to accommodate"):
        stage._calculate_chunk_worker_resources(
            total_threads=4,
            num_chunks=4,
            shared_memory_gb=sub_4gb,
        )


def test_cli_max_memory_overrides_slurm_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify CLI --max-memory-gb overrides SLURM_MEM_PER_NODE and yields -Xmx2g with 2 workers."""
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "65536")  # 64 GiB in SLURM env
    rm = ResourceManager(config={"max_memory_gb": 8.0})
    assert rm.memory_gb == 8.0

    stage = ParallelCompleteProcessingStage()
    resources = stage._calculate_chunk_worker_resources(
        total_threads=16,
        num_chunks=16,
        shared_memory_gb=rm.memory_gb,
    )
    assert resources["active_workers"] == 2
    assert resources["threads_per_worker"] == 8
    assert resources["snpsift_memory_gb"] == 2

    flag = _snpsift_memory_flag({"snpsift_memory_gb": resources["snpsift_memory_gb"]})
    assert flag == "-Xmx2g"


def test_process_chunks_parallel_dispatches_fields_and_resources_to_single_chunk(
    mock_context: PipelineContext, tmp_path: Path
) -> None:
    """End-to-end parallel dispatch passes fields_to_extract, heap, and threads to workers."""
    from concurrent.futures import Future

    stage = ParallelCompleteProcessingStage()
    vcf = tmp_path / "sample.vcf.gz"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n")
    mock_context.config["vcf_file"] = str(vcf)
    mock_context.config["fields_to_extract"] = "CHROM POS REF ALT"
    mock_context.config["filter"] = "QUAL > 30"
    mock_context.config["threads"] = 8
    mock_context.config["max_memory_gb"] = 16.0

    bed_chunk = tmp_path / "c0.bed"
    bed_chunk.write_text("chr1\t1\t100\n")

    class InlineExecutor:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def submit(self, fn, *args, **kwargs):
            fut = Future()
            try:
                res = fn(*args, **kwargs)
                fut.set_result(res)
            except Exception as exc:
                fut.set_exception(exc)
            return fut

    def mock_extract(vcf_file, bed_file, cfg, output_file):
        Path(output_file).touch()
        return output_file

    def mock_filter(variant_file, filter_string, cfg, output_file):
        Path(output_file).touch()
        return output_file

    def mock_fields(variant_file, fields, cfg, output_file, vcf_samples=None):
        assert fields == "CHROM POS REF ALT"
        Path(output_file).touch()
        return output_file

    pp_mod = "variantcentrifuge.stages.parallel_processing_stage"
    with (
        patch(f"{pp_mod}.ProcessPoolExecutor", InlineExecutor),
        patch.object(stage, "_validate_chunk_tsv", return_value=False),
        patch(f"{pp_mod}.extract_variants", side_effect=mock_extract) as p_extract,
        patch(f"{pp_mod}.apply_snpsift_filter", side_effect=mock_filter) as p_filter,
        patch(f"{pp_mod}.extract_fields_bcftools", side_effect=mock_fields) as p_fields,
    ):
        result = stage._process_chunks_parallel(mock_context, [bed_chunk])
        assert len(result) == 1
        p_extract.assert_called_once()
        p_filter.assert_called_once()
        p_fields.assert_called_once()
        fields_call = p_fields.call_args
        assert fields_call[1]["fields"] == "CHROM POS REF ALT"
        assert fields_call[1]["cfg"]["threads"] == 8
        filter_call = p_filter.call_args
        assert filter_call[0][2]["snpsift_memory_gb"] == 8


def test_filters_thread_allocation_and_command_construction(tmp_path: Path) -> None:
    """Verify bcftools receives max(0, threads-1) and bgzip receives threads."""
    from variantcentrifuge.filters import (
        _resolve_filter_threads,
        apply_bcftools_prefilter,
        apply_snpsift_filter,
        extract_variants,
    )

    # 1. Thread resolution helper
    assert _resolve_filter_threads({"threads": 1}) == (0, 0)
    assert _resolve_filter_threads({"threads": 4}) == (3, 3)
    assert _resolve_filter_threads({"threads": 4, "bcftools_threads": 2}) == (3, 2)
    assert _resolve_filter_threads({}) == (0, 0)

    bed = tmp_path / "dummy.bed"
    bed.write_text("chr1\t1\t10\n")
    vcf = tmp_path / "dummy.vcf.gz"
    vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\n")

    captured_cmds: list[list[str]] = []

    def mock_run_command(cmd, *args, **kwargs):
        captured_cmds.append(list(cmd))
        out_f = kwargs.get("output_file")
        if out_f:
            Path(out_f).touch()
        res = Mock()
        res.stderr = ""
        res.returncode = 0
        return res

    filters_mod = "variantcentrifuge.filters"
    with (
        patch(f"{filters_mod}.run_command", side_effect=mock_run_command),
        patch(f"{filters_mod}._count_vcf_records", return_value=1),
        patch(f"{filters_mod}.validate_bed_vcf_contigs"),
        patch("variantcentrifuge.vcf_eff_one_per_line.restrict_vcf_ann_to_alt", return_value=False),
    ):
        # Case A: 1 thread -> no --threads in bcftools, no -@ in bgzip
        captured_cmds.clear()
        out_a = str(tmp_path / "out_a.vcf.gz")
        apply_bcftools_prefilter(str(vcf), out_a, 'FILTER="PASS"', {"threads": 1})
        for cmd in captured_cmds:
            assert "--threads" not in cmd

        captured_cmds.clear()
        extract_variants(str(vcf), str(bed), {"threads": 1}, out_a)
        for cmd in captured_cmds:
            if "bcftools" in cmd[0]:
                assert "--threads" not in cmd

        captured_cmds.clear()
        apply_snpsift_filter(str(vcf), "QUAL > 30", {"threads": 1, "snpsift_memory_gb": 2}, out_a)
        bgzip_cmd = next(c for c in captured_cmds if c[0] == "bgzip")
        assert "-@" not in bgzip_cmd
        bcf_idx = next(c for c in captured_cmds if c[0] == "bcftools" and c[1] == "index")
        assert "--threads" not in bcf_idx

        # Case B: 4 threads -> --threads 3 in bcftools, -@ 3 in bgzip (1 foreground + 3 pool)
        captured_cmds.clear()
        apply_bcftools_prefilter(str(vcf), out_a, 'FILTER="PASS"', {"threads": 4})
        view_cmd = next(c for c in captured_cmds if c[0] == "bcftools" and c[1] == "view")
        assert view_cmd[view_cmd.index("--threads") + 1] == "3"
        idx_cmd = next(c for c in captured_cmds if c[0] == "bcftools" and c[1] == "index")
        assert idx_cmd[idx_cmd.index("--threads") + 1] == "3"

        captured_cmds.clear()
        extract_variants(str(vcf), str(bed), {"threads": 4}, out_a)
        norm_cmd = next(c for c in captured_cmds if c[0] == "bcftools" and c[1] == "norm")
        assert norm_cmd[norm_cmd.index("--threads") + 1] == "3"

        captured_cmds.clear()
        apply_snpsift_filter(str(vcf), "QUAL > 30", {"threads": 4, "snpsift_memory_gb": 2}, out_a)
        bgzip_4 = next(c for c in captured_cmds if c[0] == "bgzip")
        assert bgzip_4[bgzip_4.index("-@") + 1] == "3"
        idx_4 = next(c for c in captured_cmds if c[0] == "bcftools" and c[1] == "index")
        assert idx_4[idx_4.index("--threads") + 1] == "3"


def test_resource_allocation_no_snpsift_allows_sub_4gb(
    mock_context: PipelineContext, tmp_path: Path
) -> None:
    """Verify runs without SnpSift (e.g. late_filtering) succeed on sub-4GB budgets."""
    from concurrent.futures import Future

    stage = ParallelCompleteProcessingStage()
    res = stage._calculate_chunk_worker_resources(
        total_threads=2,
        num_chunks=2,
        shared_memory_gb=2.0,
        requires_snpsift=False,
    )
    assert res["active_workers"] == 2
    assert res["threads_per_worker"] == 1
    assert res["snpsift_memory_gb"] == 2
    assert res["max_memory_gb"] == 1.0

    vcf = tmp_path / "sample.vcf.gz"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n")
    mock_context.config["vcf_file"] = str(vcf)
    mock_context.config["fields_to_extract"] = "CHROM POS REF ALT"
    mock_context.config["filter"] = "QUAL > 30"
    mock_context.config["late_filtering"] = True
    mock_context.config["threads"] = 2
    mock_context.config["max_memory_gb"] = 2.0

    bed_chunks = [tmp_path / "c0.bed", tmp_path / "c1.bed"]
    for b in bed_chunks:
        b.write_text("chr1\t1\t100\n")

    class InlineExecutor:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def submit(self, fn, *args, **kwargs):
            fut = Future()
            fut.set_result((Path(args[4]) / f"{args[3]}.chunk_{args[0]}.extracted.tsv", {}, {}))
            return fut

    pp_mod = "variantcentrifuge.stages.parallel_processing_stage"
    with (
        patch(f"{pp_mod}.ProcessPoolExecutor", InlineExecutor),
        patch.object(stage, "_validate_chunk_tsv", return_value=False),
        patch.object(stage, "_process_single_chunk", return_value=(tmp_path / "out.tsv", {}, {})),
    ):
        result = stage._process_chunks_parallel(mock_context, bed_chunks)
        assert len(result) == 2

    # Verify no filter expression also passes
    mock_context.config["filter"] = None
    mock_context.config["filters"] = None
    mock_context.config["late_filtering"] = False
    with (
        patch(f"{pp_mod}.ProcessPoolExecutor", InlineExecutor),
        patch.object(stage, "_validate_chunk_tsv", return_value=False),
        patch.object(stage, "_process_single_chunk", return_value=(tmp_path / "out.tsv", {}, {})),
    ):
        result = stage._process_chunks_parallel(mock_context, bed_chunks)
        assert len(result) == 2
