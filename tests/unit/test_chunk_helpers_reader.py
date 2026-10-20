"""Gene-aware chunk reading and sorting helpers used by ChunkedAnalysisStage."""

from __future__ import annotations

import gzip
from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.chunk_helpers import (
    apply_chunk_annotations,
    calculate_inheritance_safe_chunk_size,
    ensure_sorted_by_gene,
    read_tsv_in_gene_chunks,
)


def _context(tmp_path: Path, **config: object) -> PipelineContext:
    return PipelineContext(Namespace(), dict(config), Workspace(tmp_path / "w", "t"))


def _tsv(path: Path, genes: list[str]) -> Path:
    path.write_text("GENE\tREF\tPOS\n" + "".join(f"{g}\tNA\t{i}\n" for i, g in enumerate(genes)))
    return path


def _read_all(chunks) -> pd.DataFrame:
    return pd.concat(list(chunks), ignore_index=True)


@pytest.mark.parametrize(
    ("genes", "sizes"),
    [
        (["A", "A", "B", "B", "C"], [2, 2, 1]),
        (["A", "A", "A", "A", "A", "B"], [5, 1]),  # one gene never split
    ],
)
def test_chunks_respect_gene_boundaries_and_keep_na(tmp_path, genes, sizes):
    tsv = _tsv(tmp_path / "in.tsv", genes)
    chunks = list(read_tsv_in_gene_chunks(tsv, "GENE", 2))
    assert [len(c) for c in chunks] == sizes
    df = _read_all(chunks)
    assert df["GENE"].tolist() == genes
    assert set(df["REF"]) == {"NA"}  # not parsed as missing


@pytest.mark.parametrize(
    ("genes", "sizes"),
    [
        (["A"] * 5 + ["B"], [4, 1, 1]),  # no boundary: buffer yielded whole (gene split)
        (["A", "B", "B", "B", "B", "C"], [1, 4, 1]),  # emergency yield at a boundary
    ],
)
def test_emergency_buffer_yield(tmp_path, genes, sizes):
    tsv = _tsv(tmp_path / "in.tsv", genes)
    context = _context(tmp_path, max_variants_per_gene=3)
    chunks = list(read_tsv_in_gene_chunks(tsv, "GENE", 1, context))
    assert [len(c) for c in chunks] == sizes
    assert _read_all(chunks)["GENE"].tolist() == genes


def test_requested_chunk_size_capped_by_inheritance_safe_size(tmp_path, monkeypatch):
    import variantcentrifuge.stages.chunk_helpers as mod

    monkeypatch.setattr(mod, "calculate_inheritance_safe_chunk_size", lambda n, ctx: 2)
    tsv = _tsv(tmp_path / "in.tsv", ["A", "A", "B", "C"])
    context = _context(tmp_path)
    context.vcf_samples = ["S1"]
    assert [len(c) for c in read_tsv_in_gene_chunks(tsv, "GENE", 100, context)] == [2, 2]
    assert [len(c) for c in read_tsv_in_gene_chunks(tsv, "GENE", 1, context)] == [2, 1, 1]


@pytest.mark.parametrize(
    ("samples", "expected"),
    [(10, 50000), (300, 25000), (800, 15000), (2000, 8000), (4000, 6000), (20000, 2500)],
)
def test_safe_chunk_size_fallback_tiers(samples, expected):
    assert calculate_inheritance_safe_chunk_size(samples, None) == expected


def test_safe_chunk_size_uses_resource_manager(tmp_path):
    class RM:
        def auto_chunk_size(self, total, samples):
            return 1234

    context = _context(tmp_path)
    context.resource_manager = RM()
    assert calculate_inheritance_safe_chunk_size(10, context) == 1234


def test_safe_chunk_size_falls_back_when_resource_manager_fails(tmp_path):
    class RM:
        def auto_chunk_size(self, total, samples):
            raise RuntimeError("no memory info")

    context = _context(tmp_path)
    context.resource_manager = RM()
    assert calculate_inheritance_safe_chunk_size(10, context) == 50000


@pytest.mark.parametrize("gz", [False, True])
def test_ensure_sorted_by_gene_sorts_stably(tmp_path, gz):
    context = _context(tmp_path)
    context.workspace.intermediate_dir.mkdir(parents=True, exist_ok=True)
    text = "POS\tGENE\n1\tB\n2\tA\n3\tB\n4\tA\n"
    src = tmp_path / ("in.tsv.gz" if gz else "in.tsv")
    if gz:
        with gzip.open(src, "wt") as f:
            f.write(text)
    else:
        src.write_text(text)
    out = ensure_sorted_by_gene(src, "GENE", context)
    assert out != src
    assert out.read_text() == "POS\tGENE\n2\tA\n4\tA\n1\tB\n3\tB\n"
    assert ensure_sorted_by_gene(src, "GENE", context) == out  # reused


def test_ensure_sorted_by_gene_passthrough_cases(tmp_path):
    src = _tsv(tmp_path / "in.tsv", ["B", "A"])
    assert ensure_sorted_by_gene(src, "GENE", _context(tmp_path, assume_sorted=True)) == src
    context = _context(tmp_path)
    context.workspace.intermediate_dir.mkdir(parents=True, exist_ok=True)
    assert ensure_sorted_by_gene(src, "MISSING", context) == src


def test_apply_chunk_annotations_with_gene_list(tmp_path):
    genes = tmp_path / "panel.txt"
    genes.write_text("A\n")
    context = _context(tmp_path)
    context.annotation_configs = {"gene_lists": [str(genes)]}
    df = pd.DataFrame({"CHROM": ["1", "1"], "POS": ["1", "2"], "GENE": ["A", "B"]})
    out = apply_chunk_annotations(df, context)
    assert len(out) == 2 and len(out.columns) > len(df.columns)


def test_apply_chunk_annotations_reraises(tmp_path, monkeypatch):
    import variantcentrifuge.annotator as annotator

    def boom(_cfg):
        raise ValueError("bad annotation config")

    monkeypatch.setattr(annotator, "load_custom_features", boom)
    context = _context(tmp_path)
    context.annotation_configs = {"gene_lists": ["x"]}
    with pytest.raises(ValueError, match="bad annotation config"):
        apply_chunk_annotations(pd.DataFrame({"GENE": ["A"]}), context)
