"""Reporting of genes dropped by the transcript filter and --transcript-fallback (#231)."""

import gzip
import logging
from argparse import Namespace
from concurrent.futures import Future
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.parallel_processing_stage import ParallelCompleteProcessingStage
from variantcentrifuge.stages.processing_stages import TranscriptFilterStage
from variantcentrifuge.validators import validate_cli_args

TSV_NAME = "transcript_filter_dropped_genes.tsv"
TX_LIST = "NM_100001.2,NM_300001.1,NM_999999.1"
# GENEA: exact match; GENEC: version-skewed match; GENEB: its only transcript is not requested.
RECORDS = [
    (100, [("GENEA", "NM_100001.2"), ("GENEA", "NM_100002.1")]),
    (200, [("GENEC", "NM_300001.3")]),
    (300, [("GENEB", "NM_200001.1")]),
]


def _vcf_text(records) -> str:
    lines = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]
    for pos, anns in records:
        ann = ",".join(
            f"C|missense_variant|MODERATE|{g}|{g}|transcript|{tx}|protein_coding||||||||"
            for g, tx in anns
        )
        lines.append(f"chr1\t{pos}\t.\tT\tC\t.\tPASS\tANN={ann}\tGT\t0/1")
    return "\n".join(lines) + "\n"


def _genes(vcf_gz: Path) -> list[str]:
    with gzip.open(vcf_gz, "rt") as fh:
        return [ln.split("ANN=")[1].split("|")[3] for ln in fh if not ln.startswith("#")]


def _run_stage(tmp_path: Path, records=RECORDS, **cfg) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    vcf = tmp_path / "in.vcf"
    vcf.write_text(_vcf_text(records), encoding="utf-8")
    ctx = PipelineContext(
        args=None,
        config={"transcript_list": TX_LIST, **cfg},
        workspace=Workspace(output_dir=tmp_path / "out", base_name="s"),
    )
    ctx.extracted_vcf = vcf
    TranscriptFilterStage()._process(ctx)
    return ctx.transcript_filtered_vcf


def test_gene_without_listed_transcript_is_warned_listed_and_kept_by_fallback(
    tmp_path: Path, caplog
) -> None:
    diag = tmp_path / "diag"
    with caplog.at_level(logging.INFO):
        out = _run_stage(tmp_path / "default", diagnostics_output=str(diag))
    assert sorted(_genes(out)) == ["GENEA", "GENEC"]
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1 and "GENEB" in warnings[0]
    assert (diag / TSV_NAME).read_text().splitlines() == ["gene", "GENEB"]

    kept = _run_stage(tmp_path / "keep", transcript_fallback="keep-gene")
    assert sorted(_genes(kept)) == ["GENEA", "GENEB", "GENEC"]


def test_info_summary_counts_exact_skewed_and_genes(tmp_path: Path, caplog) -> None:
    with caplog.at_level(logging.INFO):
        _run_stage(tmp_path)
    summary = [
        r.getMessage()
        for r in caplog.records
        if r.levelno == logging.INFO and "Transcript filtering" in r.getMessage()
    ]
    assert len(summary) == 1
    assert "for 3 transcript IDs" in summary[0]
    assert "matched 2/3 (exact 1, version-skewed 1)" in summary[0]
    assert "genes 3 -> 2" in summary[0]


def test_warning_gene_list_is_truncated(tmp_path: Path, caplog) -> None:
    records = [(100, [("GENEA", "NM_100001.2")])] + [
        (200 + i, [(f"GENE{i:02d}", f"NM_5{i:05d}.1")]) for i in range(12)
    ]
    with caplog.at_level(logging.INFO):
        _run_stage(tmp_path, records=records)
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "12 genes" in warnings[0] and "(+2 more)" in warnings[0]
    assert "GENE11" not in warnings[0]


def test_dropped_genes_tsv_only_with_diagnostics_output(tmp_path: Path) -> None:
    _run_stage(tmp_path)
    assert not list(tmp_path.rglob(TSV_NAME))


def test_chunked_summary_equals_single_pass(tmp_path: Path, caplog) -> None:
    diag = tmp_path / "diag"
    with caplog.at_level(logging.INFO):
        _run_stage(tmp_path / "single", diagnostics_output=str(diag))
    single = [
        r.getMessage()
        for r in caplog.records
        if r.levelno >= logging.INFO and "Transcript filtering" in r.getMessage()
    ]
    caplog.clear()

    ws = Mock(spec=Workspace)
    ws.base_name = "s"
    ws.intermediate_dir = tmp_path / "inter"
    ws.intermediate_dir.mkdir()
    cfg = {
        "threads": 2,
        "max_memory_gb": 64.0,
        "vcf_file": str(tmp_path / "input.vcf.gz"),
        "gzip_intermediates": False,
        "fields_to_extract": "CHROM POS REF ALT",
        "transcript_list": TX_LIST,
        "diagnostics_output": str(diag),
    }
    ctx = PipelineContext(args=Namespace(), config=cfg, workspace=ws)
    ctx.vcf_samples = ["S1"]
    chunks = [tmp_path / "c0.bed", tmp_path / "c1.bed"]
    for b in chunks:
        b.write_text("chr1\t1\t1000\n")
    parts = {"c0.bed": RECORDS[:2], "c1.bed": RECORDS[2:]}

    def fake_extract(vcf_file, bed_file, cfg, output_file):
        with gzip.open(output_file, "wt") as fh:
            fh.write(_vcf_text(parts[Path(bed_file).name]))

    def fake_fields(variant_file, fields, cfg, output_file, vcf_samples=None):
        Path(output_file).write_text("CHROM\tPOS\tREF\tALT\n")

    class InlineExecutor:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            pass

        def submit(self, fn, *args, **kwargs):
            fut = Future()
            fut.set_result(fn(*args, **kwargs))
            return fut

    pp_mod = "variantcentrifuge.stages.parallel_processing_stage"
    stage = ParallelCompleteProcessingStage()
    with (
        caplog.at_level(logging.INFO),
        patch(f"{pp_mod}.ProcessPoolExecutor", InlineExecutor),
        patch(f"{pp_mod}.extract_variants", side_effect=fake_extract),
        patch(f"{pp_mod}.extract_fields_bcftools", side_effect=fake_fields),
        patch.object(stage, "_validate_chunk_tsv", return_value=False),
    ):
        stage._process_chunks_parallel(ctx, chunks)
    chunked = [
        r.getMessage()
        for r in caplog.records
        if r.levelno >= logging.INFO and "Transcript filtering" in r.getMessage()
    ]
    assert len(single) == 2  # INFO summary + WARNING
    assert chunked == single
    assert (diag / TSV_NAME).read_text().splitlines() == ["gene", "GENEB"]


def test_keep_gene_requires_transcript_list() -> None:
    args = Namespace(transcript_fallback="keep-gene")
    with pytest.raises(ValueError, match="--transcript-fallback"):
        validate_cli_args(args, {})


def test_transcript_fallback_cli_choices() -> None:
    from variantcentrifuge.cli import create_parser

    parser = create_parser()
    vcf = ["-v", "in.vcf.gz"]
    args = parser.parse_args([*vcf, "--transcript-fallback", "keep-gene"])
    assert args.transcript_fallback == "keep-gene"
    assert parser.parse_args(vcf).transcript_fallback is None
    with pytest.raises(SystemExit):
        parser.parse_args([*vcf, "--transcript-fallback", "all"])
