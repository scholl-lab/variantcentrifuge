"""cfg["temp_dir"] must reach every bulk temp writer through the real stage code (#213)."""

import contextlib
import tempfile
from pathlib import Path
from unittest.mock import patch

from tests.unit.stages.test_processing_stages_critical import (  # noqa: F401 (fixtures)
    base_context,
    mock_workspace,
)
from variantcentrifuge.stages.processing_stages import (
    FieldExtractionStage,
    ParallelCompleteProcessingStage,
    SnpSiftFilterStage,
    VariantExtractionStage,
)

TMP = "/scratch/vc_temp_x"


def test_snpsift_filter_stage_passes_temp_dir(base_context, tmp_path):  # noqa: F811
    vcf = tmp_path / "in.vcf.gz"
    vcf.touch()
    base_context.config = {"filter": "(AF > 0.1)", "temp_dir": TMP}
    base_context.data = vcf
    with patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter") as m:
        SnpSiftFilterStage()._process(base_context)
    assert m.call_args.args[2]["temp_dir"] == TMP


def test_field_extraction_stage_passes_temp_dir(base_context, tmp_path):  # noqa: F811
    vcf = tmp_path / "f.vcf.gz"
    vcf.touch()
    base_context.config = {"fields_to_extract": "CHROM POS", "temp_dir": TMP}
    base_context.filtered_vcf = vcf
    with patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools") as m:
        FieldExtractionStage()._process(base_context)
    assert m.call_args.kwargs["cfg"]["temp_dir"] == TMP


def test_variant_extraction_stage_passes_temp_dir(base_context, tmp_path):  # noqa: F811
    vcf = tmp_path / "in.vcf.gz"
    vcf.touch()
    base_context.config = {"vcf_file": str(vcf), "temp_dir": TMP}
    base_context.gene_bed_file = tmp_path / "g.bed"
    with (
        patch("variantcentrifuge.stages.processing_stages.extract_variants") as m,
        contextlib.suppress(Exception),  # only the call arguments matter here
    ):
        VariantExtractionStage()._process(base_context)
    assert m.call_args.kwargs["cfg"]["temp_dir"] == TMP


def test_parallel_chunk_passes_temp_dir(tmp_path: Path):
    bed = tmp_path / "c.bed"
    bed.write_text("1\t0\t100\n")
    cfg = {
        "threads_per_chunk": 1,
        "fields_to_extract": "CHROM POS",
        "filters": "(AF > 0.1)",
        "temp_dir": TMP,
        "vcf_samples": ["S1"],
    }
    p = "variantcentrifuge.stages.processing_stages."
    with (
        patch(p + "extract_variants") as ev,
        patch(p + "apply_snpsift_filter") as af,
        patch(p + "extract_fields_bcftools") as ef,
    ):
        ParallelCompleteProcessingStage()._process_single_chunk(
            0, bed, "in.vcf.gz", "s", tmp_path, cfg
        )
    assert ev.call_args.kwargs["cfg"]["temp_dir"] == TMP
    assert af.call_args.args[2]["temp_dir"] == TMP
    assert ef.call_args.kwargs["cfg"]["temp_dir"] == TMP


def _run_parallel_chunks(tmp_path: Path, monkeypatch) -> tuple[list[str], Path, Path, object]:
    """Run the real chunk-worker path (external tools faked); return mkstemp paths."""
    from argparse import Namespace
    from concurrent.futures import ThreadPoolExecutor
    from unittest.mock import Mock

    from variantcentrifuge.pipeline_core.context import PipelineContext
    from variantcentrifuge.pipeline_core.workspace import Workspace

    sentinel, tdir = tmp_path / "system_tmp", tmp_path / "temp_dir"
    sentinel.mkdir()
    tdir.mkdir()
    monkeypatch.setenv("TMPDIR", str(sentinel))
    monkeypatch.setattr(tempfile, "tempdir", str(sentinel))
    made: list[str] = []
    real_mkstemp = tempfile.mkstemp

    def spy_mkstemp(*a, **kw):
        fd, name = real_mkstemp(*a, **kw)
        made.append(name)
        return fd, name

    def fake_run(cmd, output_file=None, **kw):
        out = output_file or (cmd[cmd.index("-o") + 1] if "-o" in cmd else None)
        if out:
            Path(out).write_text("")
        return Mock(stderr="", returncode=0)

    ws = Mock(spec=Workspace)
    ws.base_name = "s"
    ws.intermediate_dir = tmp_path / "inter"
    ws.intermediate_dir.mkdir()
    cfg = {
        "threads": 2,
        "max_memory_gb": 64.0,
        "vcf_file": str(tmp_path / "in.vcf.gz"),
        "fields_to_extract": "CHROM POS",
        "filters": "(QUAL > 1)",
        "temp_dir": str(tdir),
    }
    ctx = PipelineContext(args=Namespace(), config=cfg, workspace=ws)
    ctx.vcf_samples = ["S1"]
    beds = [tmp_path / "c0.bed", tmp_path / "c1.bed"]
    for b in beds:
        b.write_text("1\t0\t100\n")
    stage = ParallelCompleteProcessingStage()
    pp = "variantcentrifuge.stages.parallel_processing_stage."
    with (
        patch(pp + "ProcessPoolExecutor", ThreadPoolExecutor),
        patch("tempfile.mkstemp", spy_mkstemp),
        patch("variantcentrifuge.filters.run_command", fake_run),
        patch("variantcentrifuge.filters.validate_bed_vcf_contigs"),
        patch("variantcentrifuge.vcf_eff_one_per_line.restrict_vcf_ann_to_alt", return_value=False),
        patch("variantcentrifuge.extractor.subprocess.run", return_value=Mock(returncode=0)),
    ):
        stage._process_chunks_parallel(ctx, beds)
    return made, sentinel, tdir, stage


def test_parallel_chunk_workers_write_temps_under_temp_dir(tmp_path, monkeypatch):
    """#226: worker temps go to --temp-dir, never to the system tmp (TMPDIR)."""
    made, sentinel, tdir, _ = _run_parallel_chunks(tmp_path, monkeypatch)
    assert len(made) == 6  # extract_variants, SnpSift filter, bcftools query x 2 chunks
    assert all(Path(p).parent == tdir for p in made), made
    assert not any(sentinel.iterdir())


def test_parallel_chunks_need_no_manager(tmp_path, monkeypatch):
    """#226: no multiprocessing.Manager (its AF_UNIX socket lives in TMPDIR); times still merged."""
    no_manager = AssertionError("multiprocessing.Manager must not be started")
    with patch("multiprocessing.managers.SyncManager.start", side_effect=no_manager):
        *_, stage = _run_parallel_chunks(tmp_path, monkeypatch)
    assert set(stage._subtask_times) == {
        "variant_extraction",
        "snpsift_filtering",
        "field_extraction",
    }
