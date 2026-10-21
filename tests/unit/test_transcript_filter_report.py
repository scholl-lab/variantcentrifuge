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
            f"C|{(a[2:] or ('missense_variant',))[0]}|MODERATE|{a[0]}|{a[0]}|transcript|{a[1]}"
            "|protein_coding||||||||"
            for a in anns
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


def _run_parallel(tmp_path: Path, parts: dict, reused: tuple[int, ...] = (), **extra) -> Path:
    """Run the parallel stage's chunk loop inline; chunk files named in ``parts`` (c0.bed, ...)."""
    ws = Mock(spec=Workspace)
    ws.base_name = "s"
    ws.intermediate_dir = tmp_path / "inter"
    ws.intermediate_dir.mkdir(parents=True)
    cfg = {
        "threads": 2,
        "max_memory_gb": 64.0,
        "vcf_file": str(tmp_path / "input.vcf.gz"),
        "gzip_intermediates": False,
        "fields_to_extract": "CHROM POS REF ALT",
        "transcript_list": TX_LIST,
        **extra,
    }
    ctx = PipelineContext(args=Namespace(), config=cfg, workspace=ws)
    ctx.vcf_samples = ["S1"]
    chunks = [tmp_path / name for name in parts]
    for b in chunks:
        b.write_text("chr1\t1\t1000\n")

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
        patch(f"{pp_mod}.ProcessPoolExecutor", InlineExecutor),
        patch(f"{pp_mod}.extract_variants", side_effect=fake_extract),
        patch(f"{pp_mod}.extract_fields_bcftools", side_effect=fake_fields),
        patch.object(
            stage,
            "_validate_chunk_tsv",
            side_effect=lambda path, **_: any(f"chunk_{i}." in path.name for i in reused),
        ),
    ):
        stage._process_chunks_parallel(ctx, chunks)
    return ws.intermediate_dir


def _rows(vcf_paths) -> list[str]:
    rows = []
    for path in vcf_paths:
        with gzip.open(path, "rt") as fh:
            rows += [ln.split("ANN=")[1].split("|")[6] for ln in fh if not ln.startswith("#")]
    return sorted(rows)


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

    with caplog.at_level(logging.INFO):
        _run_parallel(
            tmp_path, {"c0.bed": RECORDS[:2], "c1.bed": RECORDS[2:]}, diagnostics_output=str(diag)
        )
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
    assert "readthrough" in parser.format_help()


# GENEG matched in chunk 0, where GENEB appears only as an upstream entry; chunk 1 holds GENEG/LOCX
# only as neighbour entries of a GENEH variant, GENEB whose only transcript is unlisted, and a
# conserved_intergenic_variant between two genes.
NB_TX = "NM_100001.1,NM_300003.1"
NB_PARTS = {
    "c0.bed": [
        (100, [("GENEG", "NM_100001.1"), ("GENEB", "NM_200009.1", "upstream_gene_variant")])
    ],
    "c1.bed": [
        (
            5000,
            [
                ("GENEH", "NM_300003.1"),
                ("GENEG", "NM_100002.1", "downstream_gene_variant"),
                ("LOCX", "NR_900001.1", "upstream_gene_variant"),
            ],
        ),
        (6000, [("GENEB", "NM_200001.1")]),
        (7000, [("LOCY-LOCZ", "NR_800001.1", "conserved_intergenic_variant")]),
    ],
}


def test_keep_gene_ignores_neighbour_entries(tmp_path: Path) -> None:
    single = _run_stage(
        tmp_path / "single",
        records=[r for recs in NB_PARTS.values() for r in recs],
        transcript_list=NB_TX,
        transcript_fallback="keep-gene",
    )
    assert _rows([single]) == ["NM_100001.1", "NM_200001.1", "NM_300003.1"]


def test_dropped_list_ignores_neighbour_only_genes(tmp_path: Path, caplog) -> None:
    diag = tmp_path / "diag"
    with caplog.at_level(logging.INFO):
        _run_parallel(tmp_path, NB_PARTS, transcript_list=NB_TX, diagnostics_output=str(diag))
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1 and "GENEB" in warnings[0] and "LOCX" not in warnings[0]
    assert "LOCY-LOCZ" not in warnings[0]
    assert (diag / TSV_NAME).read_text().splitlines() == ["gene", "GENEB"]


def test_resume_with_reused_chunks_skips_summary(tmp_path: Path, caplog) -> None:
    diag = tmp_path / "diag"
    parts = {"c0.bed": RECORDS[:2], "c1.bed": RECORDS[2:]}
    with caplog.at_level(logging.INFO):
        _run_parallel(
            tmp_path, parts, reused=(0,), transcript_list=TX_LIST, diagnostics_output=str(diag)
        )
    msgs = [r.getMessage() for r in caplog.records if "ranscript filter" in r.getMessage()]
    assert not [r for r in caplog.records if r.levelno == logging.WARNING]
    assert not (diag / TSV_NAME).exists()
    assert len(msgs) == 1 and "skipped" in msgs[0] and "1 of 2 chunks reused" in msgs[0]


def test_explicit_fallback_none_keeps_default_fingerprints(tmp_path: Path) -> None:
    from variantcentrifuge.checkpoint_helpers import hash_configuration
    from variantcentrifuge.cli import create_parser
    from variantcentrifuge.cli_config import populate_cli_config, resolve_presets_and_tumor_normal

    parser = create_parser()
    hashes = []
    for extra in ([], ["--transcript-fallback", "none"]):
        args = parser.parse_args(["-v", "in.vcf.gz", "--transcript-list", "NM_1.1", *extra])
        cfg: dict = {}
        filters, tn_vars = resolve_presets_and_tumor_normal(args, cfg, None)
        populate_cli_config(args, cfg, parser, "GRCh38", filters, "CHROM POS", tn_vars)
        assert cfg.get("transcript_fallback") is None
        hashes.append(hash_configuration(cfg))
    assert hashes[0] == hashes[1]


@pytest.mark.parametrize("fallback", ["keep-gene", None])
def test_split_gene_bed_gives_single_pass_rows(tmp_path: Path, fallback) -> None:
    """A regions BED can split one gene into intervals in different chunks."""
    from variantcentrifuge.pipeline import create_stages_from_config
    from variantcentrifuge.transcript_filter import filter_vcf_to_transcripts
    from variantcentrifuge.utils import split_bed_file

    bed = tmp_path / "regions.bed"
    bed.write_text("1\t90\t110\n1\t290\t310\n")
    records = [(100, [("GENEG", "NM_1.1")]), (300, [("GENEG", "NM_2.1")])]
    cfg = {"fields_to_extract": "CHROM POS", "threads": 4, "transcript_list": "NM_1.1"}
    if fallback:
        cfg["transcript_fallback"] = fallback
    names = [s.name for s in create_stages_from_config(cfg)]

    def run(recs, name: str) -> list[str]:
        vcf = tmp_path / f"{name}.vcf"
        vcf.write_text(_vcf_text(recs).replace("chr1", "1"))
        out = tmp_path / f"{name}.out.vcf.gz"
        filter_vcf_to_transcripts(str(vcf), str(out), {"NM_1.1"}, fallback=fallback or "none")
        return _rows([out])

    if "parallel_complete_processing" in names:  # what the chosen stages do with this BED
        (tmp_path / "chunks").mkdir()
        chunks = split_bed_file(str(bed), 2, str(tmp_path / "chunks"))
        assert len(chunks) == 2
        kept = []
        for i, chunk in enumerate(chunks):
            start, end = (int(x) for x in Path(chunk).read_text().split()[1:3])
            kept += run([r for r in records if start <= r[0] < end], f"chunk{i}")
    else:
        assert "transcript_filtering" in names
        kept = run(records, "single")
    assert sorted(kept) == ["NM_1.1"]
    assert ("parallel_complete_processing" in names) == (fallback is None)


def test_keep_gene_routing_reads_args_when_config_is_a_path(tmp_path: Path) -> None:
    """Python API: args.config may be a file path, so transcript options are read from args too."""
    from variantcentrifuge.pipeline import build_pipeline_stages

    cfg_file = tmp_path / "config.json"
    cfg_file.write_text("{}")
    args = Namespace(
        config=str(cfg_file),
        threads=4,
        transcript_list=TX_LIST,
        transcript_fallback="keep-gene",
        fields_to_extract="CHROM POS",
    )
    stages = {s.name: s for s in build_pipeline_stages(args)}
    assert "parallel_complete_processing" not in stages
    assert "transcript_filtering" in stages

    vcf = tmp_path / "in.vcf"
    vcf.write_text(_vcf_text(RECORDS), encoding="utf-8")
    ctx = PipelineContext(args=args, config={}, workspace=Workspace(tmp_path / "out", "s"))
    ctx = stages["configuration_loading"]._process(ctx)
    ctx.extracted_vcf = vcf
    ctx = stages["transcript_filtering"]._process(ctx)
    assert sorted(_genes(ctx.transcript_filtered_vcf)) == ["GENEA", "GENEB", "GENEC"]


def test_default_chunk_fingerprint_equals_v0_24_0(tmp_path: Path) -> None:
    """Chunk TSVs from a v0.24.0 run stay reusable: the default-config fingerprint is unchanged.

    Pinned value computed with v0.24.0 (f751241b) ParallelCompleteProcessingStage._chunk_fingerprint
    on exactly these inputs.
    """
    from types import SimpleNamespace

    bed = tmp_path / "c.bed"
    bed.write_text("chr1\t1\t100\n")
    ctx = SimpleNamespace(config={"vcf_file": "in.vcf.gz", "transcript_list": "NM_1.1"}, data=None)
    fp = ParallelCompleteProcessingStage()._chunk_fingerprint(ctx, bed, "CHROM POS")
    assert fp == "4af64340f714d02f"
