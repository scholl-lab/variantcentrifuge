"""Forced chunking must match the non-chunked loader (#183) and restore no-op runs (#184)."""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.pipeline_core.runner_planner import restore_completed_stage
from variantcentrifuge.stages.chunked_analysis_stages import (
    ChunkedAnalysisStage,
    _find_gene_column,
)
from variantcentrifuge.stages.data_loading_stages import DataFrameLoadingStage

# Genes deliberately not in sorted order so gene-aware chunking has to sort and restore.
_ROWS = [
    ("1", "100", "A", "T", "SIGNAL", "missense_variant", "0/1", "0/0"),
    ("1", "200", "C", "G", "BASELINE", "stop_gained", "0/0", "1/1"),
    ("1", "300", "G", "A", "SIGNAL", "synonymous_variant", "0/1", "0/1"),
    ("2", "400", "T", "C", "PROTECT", "missense_variant", "./.", "0/1"),
    ("2", "500", "A", "G", "BASELINE", "missense_variant", "1/1", "0/0"),
]


def _write_tsv(path: Path, gene_header: str | None) -> Path:
    gene = [gene_header] if gene_header else []
    header = ["CHROM", "POS", "REF", "ALT", *gene, "ANN[0].EFFECT", "GEN[0].GT", "GEN[1].GT"]
    lines = ["\t".join(header)]
    for row in _ROWS:
        values = list(row) if gene_header else [*row[:4], *row[5:]]
        lines.append("\t".join(values))
    path.write_text("\n".join(lines) + "\n")
    return path


def _context(tmp_path: Path, tsv: Path, **config: object) -> PipelineContext:
    workspace = Workspace(tmp_path / "work", "test")
    cfg = {"threads": 1, "chunk_size": 2, "gzip_intermediates": False, **config}
    context = PipelineContext(Namespace(), cfg, workspace)
    context.data = tsv
    context.vcf_samples = ["S1", "S2"]
    return context


def _serial(tmp_path: Path, tsv: Path) -> PipelineContext:
    context = _context(tmp_path / "serial", tsv, no_chunked_processing=True)
    return DataFrameLoadingStage()._process(context)


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        (["CHROM", "GENE", "ANN[0].GENE"], "GENE"),
        (["CHROM", "ANN[0].GENE"], "ANN[0].GENE"),
        (["CHROM", "ANN_0__GENE"], "ANN_0__GENE"),
        (["CHROM", "ANN[*].GENE"], "ANN[*].GENE"),
        (["CHROM", "EFF[*].GENE"], "EFF[*].GENE"),
        (["CHROM", "SYMBOL"], "SYMBOL"),
        (["CHROM", "Gene"], "Gene"),
        (["CHROM", "POS"], None),
        ([], None),
    ],
)
def test_find_gene_column_header_variants(header, expected):
    assert _find_gene_column(header) == expected


@pytest.mark.parametrize("gene_header", ["GENE", "ANN[0].GENE", "ANN_0__GENE", "ANN[*].GENE"])
@pytest.mark.parametrize("threads", [1, 2])
def test_forced_chunking_matches_nonchunked_loader(tmp_path, gene_header, threads):
    """#183: forced chunking yields the same rows, order, columns and rename map as serial."""
    tsv = _write_tsv(tmp_path / "input.tsv", gene_header)
    serial = _serial(tmp_path, tsv)

    context = _context(tmp_path, tsv, use_chunked_processing=True, threads=threads)
    result = ChunkedAnalysisStage()._process(context)

    assert result.current_dataframe is not None
    pd.testing.assert_frame_equal(
        result.current_dataframe.astype(str),
        serial.current_dataframe.astype(str),
    )
    assert result.column_rename_map == serial.column_rename_map
    assert result.config["dataframe_chunked_analysis_complete"] is True
    assert ChunkedAnalysisStage().get_output_files(result) == [result.chunked_analysis_tsv]


def test_missing_gene_column_falls_back_to_full_load(tmp_path):
    """#183: the fallback loads the data (as serial does) instead of leaving no DataFrame."""
    tsv = _write_tsv(tmp_path / "input.tsv", None)
    serial = _serial(tmp_path, tsv)

    context = _context(tmp_path, tsv, use_chunked_processing=True)
    stage = ChunkedAnalysisStage()
    result = stage._process(context)

    pd.testing.assert_frame_equal(result.current_dataframe, serial.current_dataframe)
    assert result.config["use_chunked_processing"] is False
    assert stage.get_output_files(result) == []

    # Resume of the fallback run: dataframe_loading re-enables chunking, no artifact exists.
    resumed = _context(tmp_path, tsv, use_chunked_processing=True)
    resumed = restore_completed_stage(stage, resumed)
    pd.testing.assert_frame_equal(resumed.current_dataframe, serial.current_dataframe)


def test_noop_stage_restores_without_artifact(tmp_path):
    """#184: a completed no-op chunk stage restores without demanding chunk output."""
    tsv = _write_tsv(tmp_path / "input.tsv", "GENE")
    serial = _serial(tmp_path, tsv)
    df = serial.current_dataframe

    stage = ChunkedAnalysisStage()
    assert stage._process(serial).current_dataframe is df
    assert stage.get_output_files(serial) == []

    restored = restore_completed_stage(stage, serial)
    assert restored.current_dataframe is df
    assert "dataframe_chunked_analysis_complete" not in restored.config


def test_real_chunked_stage_restores_artifact_in_input_order(tmp_path):
    """#184: restoring a genuinely chunked run reproduces the in-memory result."""
    tsv = _write_tsv(tmp_path / "input.tsv", "ANN[0].GENE")
    first = ChunkedAnalysisStage()._process(_context(tmp_path, tsv, use_chunked_processing=True))

    resumed = _context(tmp_path, tsv, use_chunked_processing=True)
    resumed = restore_completed_stage(ChunkedAnalysisStage(), resumed)

    pd.testing.assert_frame_equal(
        resumed.current_dataframe.astype(str), first.current_dataframe.astype(str)
    )
    assert resumed.column_rename_map == first.column_rename_map


def test_restore_still_requires_artifact_when_chunking_was_real(tmp_path):
    tsv = _write_tsv(tmp_path / "input.tsv", "GENE")
    context = _context(tmp_path, tsv, use_chunked_processing=True)
    with pytest.raises(RuntimeError, match="Cannot restore chunked_analysis"):
        ChunkedAnalysisStage()._handle_checkpoint_skip(context)


def test_header_only_input_yields_empty_frame_without_helper_column(tmp_path):
    tsv = tmp_path / "empty.tsv"
    tsv.write_text("CHROM\tPOS\tGENE\n")
    result = ChunkedAnalysisStage()._process(_context(tmp_path, tsv, use_chunked_processing=True))
    assert result.current_dataframe.empty
    assert list(result.current_dataframe.columns) == ["CHROM", "POS", "GENE"]


def test_input_file_precedence(tmp_path):
    context = _context(tmp_path, tmp_path / "data.tsv")
    find = ChunkedAnalysisStage._find_input_file
    assert find(context) == tmp_path / "data.tsv"
    for attr in (
        "extracted_tsv",
        "genotype_replaced_tsv",
        "phenotypes_added_tsv",
        "extra_columns_removed_tsv",
    ):
        setattr(context, attr, tmp_path / f"{attr}.tsv")
        assert find(context) == tmp_path / f"{attr}.tsv"
    empty = _context(tmp_path, tmp_path / "x.tsv")
    empty.data = None
    with pytest.raises(RuntimeError, match="No TSV file available"):
        find(empty)
