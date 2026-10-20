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
        (["CHROM", "ANN[0].GENE", "GENE"], "GENE"),
        (["CHROM", "ANN[*].GENE", "GENE"], "GENE"),
        (["CHROM", "ANN[0].GENE"], "ANN[0].GENE"),
        (["CHROM", "ANN_0__GENE"], "ANN_0__GENE"),
        (["CHROM", "ANN[*].GENE"], None),  # wildcard is not a GENE source
        (["CHROM", "EFF[*].GENE"], None),
        (["CHROM", "ANN[*].GENE", "ANN[0].GENE"], "ANN[0].GENE"),
        (["CHROM", "SYMBOL", "ANN[0].GENE"], "ANN[0].GENE"),
        (["CHROM", "SYMBOL"], "SYMBOL"),
        (["CHROM", "Gene"], "Gene"),
        (["CHROM", "POS"], None),
        ([], None),
    ],
)
def test_find_gene_column_header_variants(header, expected):
    assert _find_gene_column(header) == expected


@pytest.mark.parametrize("gene_header", ["GENE", "ANN[0].GENE", "ANN_0__GENE"])
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


_TRIO = {
    s: {
        "family_id": "F",
        "sample_id": s,
        "father_id": f,
        "mother_id": m,
        "sex": x,
        "affected_status": a,
    }
    for s, f, m, x, a in [
        ("Child", "Father", "Mother", "1", "2"),
        ("Father", "0", "0", "1", "1"),
        ("Mother", "0", "0", "2", "1"),
    ]
}


@pytest.mark.parametrize(
    ("gene_cols", "canonical"),
    [
        (["ANN[0].GENE", "GENE"], "GENE"),
        (["GENE", "ANN[0].GENE"], "GENE"),
        (["ANN[*].GENE", "ANN[0].GENE"], "ANN[0].GENE"),
        (["ANN[0].GENE", "ANN[*].GENE"], "ANN[0].GENE"),
    ],
)
def test_chunks_follow_canonical_gene_and_keep_comp_het(
    tmp_path, monkeypatch, gene_cols, canonical
):
    """Review r1 #1 / r2 #2: chunk by the column that becomes GENE; chunking by a decoy
    annotation column split the pair (comp-het -> autosomal_dominant_possible)."""
    from variantcentrifuge.memory import ResourceManager

    monkeypatch.setattr(ResourceManager, "auto_chunk_size", lambda self, *a, **k: 1)
    values = {c: (["TARGET", "TARGET"] if c == canonical else ["X", "Y"]) for c in gene_cols}
    header = ["CHROM", "POS", "REF", "ALT", *gene_cols, "GEN[0].GT", "GEN[1].GT", "GEN[2].GT"]
    gts = [("0/1", "0/1", "0/0"), ("0/1", "0/0", "0/1")]  # one allele from each parent
    lines = ["\t".join(header)]
    for i, pos in enumerate(("100", "200")):
        lines.append("\t".join(["1", pos, "A", "T", *(values[c][i] for c in gene_cols), *gts[i]]))
    tsv = tmp_path / "trio.tsv"
    tsv.write_text("\n".join(lines) + "\n")

    assert _find_gene_column(header) == canonical
    context = _context(tmp_path, tsv, use_chunked_processing=True, inheritance_mode="simple")
    context.vcf_samples = ["Child", "Father", "Mother"]
    context.pedigree_data = _TRIO
    df = ChunkedAnalysisStage()._process(context).current_dataframe
    assert list(df["GENE"]) == ["TARGET", "TARGET"]
    assert list(df["Inheritance_Pattern"]) == ["compound_heterozygous"] * 2


def test_gene_column_follows_shared_alias_priority(monkeypatch):
    """The selector reuses _STANDARD_COLUMN_ALIASES order rather than its own list."""
    import variantcentrifuge.stages.chunked_analysis_stages as chunked_mod
    import variantcentrifuge.stages.genotype_matrix_builder as builder_mod

    aliases = {"SECOND_GENE": "GENE", "ANN_0__GENE": "GENE"}
    monkeypatch.setattr(chunked_mod, "_STANDARD_COLUMN_ALIASES", aliases)
    monkeypatch.setattr(builder_mod, "_STANDARD_COLUMN_ALIASES", aliases)
    header = ["ANN[0].GENE", "SECOND_GENE"]
    df = pd.DataFrame([["a", "b"]], columns=["ANN_0__GENE", "SECOND_GENE"])
    builder_mod._create_standard_column_aliases(df)
    assert df["GENE"].iloc[0] == "b"
    assert _find_gene_column(header) == "SECOND_GENE"


def test_chunk_reader_keeps_na_alleles_and_literal_quotes(tmp_path):
    """Review r1 #2: chunk parsing must use the non-chunked loader's NA/quoting rules."""
    tsv = tmp_path / "na.tsv"
    tsv.write_text(
        "CHROM\tPOS\tREF\tALT\tGENE\tANN[0].EFFECT\tNOTE\tGEN[0].GT\n"
        '1\t100\tNA\tN\tG1\tNA\t"quoted"\t0/1\n'
        "1\t200\tA\tNA\tG2\tmissense\t\t0/1\n"
    )
    serial = _serial(tmp_path, tsv).current_dataframe
    chunked = (
        ChunkedAnalysisStage()
        ._process(_context(tmp_path, tsv, use_chunked_processing=True))
        .current_dataframe
    )
    pd.testing.assert_frame_equal(chunked.astype(str), serial.astype(str))
    assert chunked["REF"].tolist() == ["NA", "A"]
    assert chunked["ALT"].tolist() == ["N", "NA"]
    assert chunked["EFFECT"].tolist() == ["NA", "missense"]
    assert chunked["NOTE"].iloc[0] == '"quoted"'
    assert pd.isna(chunked["NOTE"].iloc[1])


@pytest.mark.parametrize("gene_header", ["ANN[*].GENE", "EFF[*].GENE"])
def test_wildcard_gene_header_is_not_a_gene_and_falls_back(tmp_path, gene_header):
    """Wildcard columns hold raw annotation text from bcftools; never alias them to GENE.
    Chunking then has no gene column and falls back to the full serial load."""
    tsv = tmp_path / "wild.tsv"
    tsv.write_text(f"CHROM\tPOS\t{gene_header}\n1\t100\tC|missense|X\n")
    serial = _serial(tmp_path, tsv).current_dataframe
    result = ChunkedAnalysisStage()._process(_context(tmp_path, tsv, use_chunked_processing=True))
    assert "GENE" not in serial.columns
    pd.testing.assert_frame_equal(result.current_dataframe, serial)
    assert result.config["use_chunked_processing"] is False


def test_chunk_reader_skips_overlong_rows_and_pads_short_rows_like_serial(tmp_path):
    """Review r2 #3: same bad-line policy as the serial loader."""
    tsv = tmp_path / "bad.tsv"
    tsv.write_text("CHROM\tPOS\tGENE\tREF\n1\t1\tA\tC\n1\t2\tB\tG\textra\n1\t3\tC\n1\t4\tA\tT\n")
    serial = _serial(tmp_path, tsv).current_dataframe
    chunked = (
        ChunkedAnalysisStage()
        ._process(_context(tmp_path, tsv, use_chunked_processing=True))
        .current_dataframe
    )
    pd.testing.assert_frame_equal(chunked.astype(str), serial.astype(str))
    assert chunked["POS"].tolist() == ["1", "3", "4"]
    assert pd.isna(chunked["REF"].iloc[1])


@pytest.mark.parametrize(
    ("scoring", "expected"),
    [
        (None, False),
        ({"formulas": [{"s": "Inheritance_Pattern == 'de_novo'"}]}, True),
        ({"formulas": [{"s": "CADD > 20"}, {"t": 3}]}, False),
        ({"variables": {"Inheritance_Details": "d"}}, True),
        ({"variables": {"CADD": "c"}}, False),
    ],
)
def test_scoring_inheritance_detail_detection(tmp_path, scoring, expected):
    context = _context(tmp_path, tmp_path / "x.tsv")
    context.scoring_config = scoring
    stage = ChunkedAnalysisStage()
    assert stage._check_if_scoring_needs_details(context) is expected
    assert stage._should_calculate_inheritance(context) is expected
    stage._setup_inheritance_config(context)
    assert ("inheritance_analysis_config" in context.config) is expected
    if expected:
        assert context.config["inheritance_analysis_config"]["preserve_details_for_scoring"]


def test_chunk_annotations_and_scoring_see_serial_column_names(tmp_path, monkeypatch):
    """Per-chunk annotation/scoring run on the same names as the non-chunked path."""
    import variantcentrifuge.stages.chunked_analysis_stages as mod

    seen: list[list[str]] = []

    def fake_scoring(df, _cfg):
        seen.append(list(df.columns))
        return df.assign(Score="1")

    monkeypatch.setattr(mod, "apply_scoring", fake_scoring)
    monkeypatch.setattr(
        ChunkedAnalysisStage,
        "_apply_chunk_annotations",
        lambda self, df, ctx: df.assign(Custom_Annotation="x"),
    )
    tsv = _write_tsv(tmp_path / "input.tsv", "ANN[0].GENE")
    context = _context(tmp_path, tsv, use_chunked_processing=True)
    context.annotation_configs = {"gene_lists": ["dummy"]}
    context.scoring_config = {"formulas": [{"Score": "1"}]}

    df = ChunkedAnalysisStage()._process(context).current_dataframe

    assert seen and all("GENE" in cols and "GEN_0__GT" in cols for cols in seen)
    assert set(df["Score"]) == {"1"} and set(df["Custom_Annotation"]) == {"x"}
    assert "ANN_0__GENE" in df.columns  # final reload renames like the serial loader


def test_empty_chunk_is_returned_unprocessed(tmp_path):
    empty = pd.DataFrame(columns=["GENE"])
    context = _context(tmp_path, tmp_path / "x.tsv")
    assert ChunkedAnalysisStage()._process_single_chunk(empty, context, 0) is empty


def test_parallel_chunk_error_propagates(tmp_path, monkeypatch):
    from variantcentrifuge.memory import ResourceManager

    monkeypatch.setattr(ResourceManager, "auto_chunk_size", lambda self, *a, **k: 1)
    monkeypatch.setattr(ResourceManager, "auto_workers", lambda self, **k: 2)

    def boom(self, df, ctx, idx):
        raise RuntimeError("chunk failed")

    monkeypatch.setattr(ChunkedAnalysisStage, "_process_single_chunk", boom)
    tsv = _write_tsv(tmp_path / "input.tsv", "GENE")
    with pytest.raises(RuntimeError, match="chunk failed"):
        ChunkedAnalysisStage()._process(
            _context(tmp_path, tsv, use_chunked_processing=True, threads=2)
        )


def test_get_input_files(tmp_path):
    stage = ChunkedAnalysisStage()
    context = _context(tmp_path, tmp_path / "data.tsv")
    assert stage.get_input_files(context) == [tmp_path / "data.tsv"]
    extracted = tmp_path / "extracted.tsv"
    extracted.write_text("GENE\n")
    context.extracted_tsv = extracted
    context.phenotypes_added_tsv = tmp_path / "missing.tsv"
    assert stage.get_input_files(context) == [extracted]
    context.data = None
    context.extracted_tsv = None
    assert stage.get_input_files(context) == []
