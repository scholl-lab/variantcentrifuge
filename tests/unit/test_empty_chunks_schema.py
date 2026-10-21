"""Unit tests for schema preservation on empty BEDs and empty variant chunks (#149)."""

from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from variantcentrifuge.extractor import extract_fields_bcftools
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import DataFrameLoadingStage
from variantcentrifuge.stages.parallel_processing_stage import ParallelCompleteProcessingStage


def _write_minimal_vcf(path: Path, has_variants: bool = True) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=1000>",
        (
            '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
            "Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | "
            "Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | "
            'CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO">'
        ),
        '##INFO=<ID=NMD,Number=.,Type=String,Description="Nonsense-mediated decay">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2",
    ]
    if has_variants:
        info = (
            "ANN=T|missense_variant|MODERATE|BRCA1|ENSG1|transcript|NM_1|protein_coding||||||||;"
            "NMD=(BRCA1|ENSG1|1|0.95)"
        )
        lines.append(f"chr1\t100\t.\tA\tT\t500\tPASS\t{info}\tGT\t0/1\t0/0")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def test_empty_and_populated_extraction_identical_headers(tmp_path: Path) -> None:
    """Verify that empty and populated extractions produce identical ordered headers."""
    empty_vcf = tmp_path / "empty.vcf"
    pop_vcf = tmp_path / "pop.vcf"
    _write_minimal_vcf(empty_vcf, has_variants=False)
    _write_minimal_vcf(pop_vcf, has_variants=True)

    fields = "CHROM POS ANN[0].GENE NMD[0].GENE GEN[*].GT GEN[*].DP"
    empty_tsv = tmp_path / "empty.tsv"
    pop_tsv = tmp_path / "pop.tsv"

    def mock_subprocess_run(cmd, *args, **kwargs):
        variant_file = cmd[-1]
        stdout = kwargs.get("stdout")
        if stdout and "pop.vcf" in variant_file:
            # Output matches build_bcftools_format_string query format
            ann = "T|missense_variant|MODERATE|BRCA1|ENSG1|transcript|NM_1|protein_coding||||||||"
            nmd = "(BRCA1|ENSG1|1|0.95)"
            stdout.write(f"chr1\t100\t{ann}\t{nmd}\t0/1\t30\t0/0\t45\n")
        return Mock(returncode=0)

    with patch("variantcentrifuge.extractor.subprocess.run", side_effect=mock_subprocess_run):
        extract_fields_bcftools(
            variant_file=str(empty_vcf),
            fields=fields,
            cfg={},
            output_file=str(empty_tsv),
            vcf_samples=["SAMPLE1", "SAMPLE2"],
        )
        extract_fields_bcftools(
            variant_file=str(pop_vcf),
            fields=fields,
            cfg={},
            output_file=str(pop_tsv),
            vcf_samples=["SAMPLE1", "SAMPLE2"],
        )

    empty_header = empty_tsv.read_text().splitlines()[0].split("\t")
    pop_header = pop_tsv.read_text().splitlines()[0].split("\t")

    assert empty_header == pop_header, (
        f"Empty header {empty_header} != Populated header {pop_header}"
    )


@pytest.mark.parametrize("gzip_intermediates", [False, True])
def test_zero_region_bed_produces_loadable_artifact(
    tmp_path: Path, gzip_intermediates: bool
) -> None:
    """Verify empty bed_chunks produces schema-valid header artifact loaded into DataFrame."""
    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir(parents=True, exist_ok=True)
    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = inter_dir
    ws.output_dir = tmp_path
    ws.get_intermediate_path = lambda name: inter_dir / name

    empty_bed = inter_dir / "empty.bed"
    empty_bed.write_text("")

    config = {
        "threads": 4,
        "gzip_intermediates": gzip_intermediates,
        "fields_to_extract": "CHROM POS REF ALT ANN[0].GENE GEN[*].GT",
    }
    ctx = PipelineContext(args=Namespace(), config=config, workspace=ws)
    ctx.gene_bed_file = empty_bed
    ctx.vcf_samples = ["SAMPLE1", "SAMPLE2"]

    stage = ParallelCompleteProcessingStage()
    updated_ctx = stage._process(ctx)

    assert updated_ctx.extracted_tsv.exists()

    loading_stage = DataFrameLoadingStage()
    loaded_ctx = loading_stage._process(updated_ctx)

    assert loaded_ctx.current_dataframe is not None
    assert len(loaded_ctx.current_dataframe) == 0
    assert "CHROM" in loaded_ctx.current_dataframe.columns
    cols = loaded_ctx.current_dataframe.columns
    assert "GEN[0].GT" in cols or "GEN_0__GT" in cols


def test_empty_bed_preserves_extra_sample_fields_schema(tmp_path: Path) -> None:
    """Verify empty BED artifact retains extra_sample_fields in header columns."""
    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir(parents=True, exist_ok=True)
    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = inter_dir
    ws.output_dir = tmp_path
    ws.get_intermediate_path = lambda name: inter_dir / name

    empty_bed = inter_dir / "empty.bed"
    empty_bed.write_text("")

    config = {
        "threads": 2,
        "gzip_intermediates": False,
        "fields_to_extract": "CHROM POS REF ALT GEN[*].GT",
        "append_extra_sample_fields": True,
        "extra_sample_fields": ["DP"],
    }
    ctx = PipelineContext(args=Namespace(), config=config, workspace=ws)
    ctx.gene_bed_file = empty_bed
    ctx.vcf_samples = ["S1", "S2"]

    stage = ParallelCompleteProcessingStage()
    updated_ctx = stage._process(ctx)

    header = updated_ctx.extracted_tsv.read_text().splitlines()[0].split("\t")
    assert "GEN[0].DP" in header
    assert "GEN[1].DP" in header


@pytest.mark.parametrize("all_empty", [False, True])
def test_interleaved_empty_and_populated_chunks_merge_cleanly(
    tmp_path: Path, all_empty: bool
) -> None:
    """Verify chunks where first/middle/last or all are empty merge and load into DataFrame."""
    from variantcentrifuge.extractor import get_expected_extraction_columns

    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir(parents=True, exist_ok=True)
    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = inter_dir
    ws.output_dir = tmp_path
    ws.get_intermediate_path = lambda name: inter_dir / name

    fields = "CHROM POS REF ALT ANN[0].GENE GEN[*].GT"
    config = {
        "threads": 4,
        "gzip_intermediates": False,
        "fields_to_extract": fields,
    }
    ctx = PipelineContext(args=Namespace(), config=config, workspace=ws)
    ctx.vcf_samples = ["S1", "S2"]

    cols = get_expected_extraction_columns(fields, ctx.vcf_samples)
    header = "\t".join(cols) + "\n"

    c0 = inter_dir / "chunk_0.extracted.tsv"
    c1 = inter_dir / "chunk_1.extracted.tsv"
    c2 = inter_dir / "chunk_2.extracted.tsv"
    c3 = inter_dir / "chunk_3.extracted.tsv"
    c4 = inter_dir / "chunk_4.extracted.tsv"

    # c0 is empty (first)
    c0.write_text(header)
    # c1 has row (or empty if all_empty)
    row1 = "\t".join(["chr1", "100", "A", "T", "0/1", "0/0", "BRCA1"]) + "\n"
    c1.write_text(header if all_empty else header + row1)
    # c2 is empty (middle)
    c2.write_text(header)
    # c3 has row (or empty if all_empty)
    row2 = "\t".join(["chr1", "200", "G", "C", "0/0", "0/1", "TP53"]) + "\n"
    c3.write_text(header if all_empty else header + row2)
    # c4 is empty (last)
    c4.write_text(header)

    stage = ParallelCompleteProcessingStage()
    merged_tsv = stage._merge_tsv_outputs(ctx, [c0, c1, c2, c3, c4])
    ctx.extracted_tsv = merged_tsv
    ctx.data = merged_tsv

    loading_stage = DataFrameLoadingStage()
    loaded_ctx = loading_stage._process(ctx)

    df = loaded_ctx.current_dataframe
    assert df is not None
    assert "CHROM" in df.columns
    assert "POS" in df.columns
    gt_col_0 = "GEN_0__GT" if "GEN_0__GT" in df.columns else "GEN[0].GT"
    gt_col_1 = "GEN_1__GT" if "GEN_1__GT" in df.columns else "GEN[1].GT"
    assert gt_col_0 in df.columns
    assert gt_col_1 in df.columns

    if all_empty:
        assert len(df) == 0
    else:
        assert len(df) == 2
        assert list(df["CHROM"]) == ["chr1", "chr1"]
        assert list(df["POS"]) == ["100", "200"]
        assert list(df[gt_col_0]) == ["0/1", "0/0"]
        assert list(df[gt_col_1]) == ["0/0", "0/1"]
        gene_col = "GENE" if "GENE" in df.columns else "ANN_0__GENE"
        genes = list(df[gene_col])
        assert "BRCA1" in genes and "TP53" in genes


def test_cached_chunk_rejected_on_header_mismatch(tmp_path: Path) -> None:
    """Verify that cached chunk files on disk are invalidated if their header schema differs."""
    from concurrent.futures import Future

    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir(parents=True, exist_ok=True)
    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = inter_dir
    ws.output_dir = tmp_path
    ws.get_intermediate_path = lambda name: inter_dir / name

    config = {
        "threads": 2,
        "gzip_intermediates": False,
        "fields_to_extract": "CHROM POS REF ALT GEN[*].GT",
        "vcf_file": str(tmp_path / "dummy.vcf.gz"),
    }
    Path(config["vcf_file"]).write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n")

    ctx = PipelineContext(args=Namespace(), config=config, workspace=ws)
    ctx.vcf_samples = ["S1", "S2"]

    chunk_bed0 = tmp_path / "chunk_0.bed"
    chunk_bed0.write_text("chr1\t1\t100\n")

    # Stale chunk TSV has obsolete schema without sample columns
    stale_chunk = inter_dir / "sample.chunk_0.extracted.tsv"
    stale_chunk.write_text("CHROM\tPOS\tREF\tALT\nchr1\t50\tA\tC\n")

    stage = ParallelCompleteProcessingStage()

    # When validated with expected header, stale chunk should be rejected as invalid
    expected_header = "CHROM\tPOS\tREF\tALT\tGEN[0].GT\tGEN[1].GT"
    assert not stage._validate_chunk_tsv(stale_chunk, expected_header=expected_header)

    mock_future = Future()
    mock_future.set_result((stale_chunk, {}, {}))
    mock_executor = Mock()
    mock_executor.submit.return_value = mock_future
    mock_executor.__enter__ = Mock(return_value=mock_executor)
    mock_executor.__exit__ = Mock(return_value=False)

    pp_mod = "variantcentrifuge.stages.parallel_processing_stage"
    with (
        patch(f"{pp_mod}.ProcessPoolExecutor", return_value=mock_executor),
        patch(f"{pp_mod}.as_completed", return_value=[mock_future]),
    ):
        stage._process_chunks_parallel(ctx, [chunk_bed0])
        # Crucial: verify that executor.submit was called (stale chunk was NOT reused)
        assert mock_executor.submit.call_count == 1


def test_incompatible_chunk_header_rejected_during_merge(tmp_path: Path) -> None:
    """Verify that traditional merging rejects chunks with incompatible headers."""
    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir()
    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = inter_dir
    ws.get_intermediate_path = lambda name: inter_dir / name

    chunk0 = inter_dir / "chunk0.tsv"
    chunk1 = inter_dir / "chunk1.tsv"

    chunk0.write_text("CHROM\tPOS\tANN\n")
    chunk1.write_text("CHROM\tPOS\tANN[0].GENE\tGEN[0].GT\nchr1\t100\tBRCA1\t0/1\n")

    stage = ParallelCompleteProcessingStage()
    with pytest.raises(ValueError, match="Incompatible chunk header"):
        stage._merge_tsvs_traditional([chunk0, chunk1], inter_dir / "merged.tsv")


def test_checkpoint_skip_rejects_cached_tsv_with_incompatible_header(tmp_path: Path) -> None:
    """Verify _handle_checkpoint_skip rejects cached TSV if header schema does not match."""
    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir(parents=True, exist_ok=True)
    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = inter_dir

    cached_tsv = inter_dir / "sample.extracted.tsv"
    cached_tsv.write_text("CHROM\tPOS\tANN\nchr1\t100\traw_ann\n")

    config = {
        "fields_to_extract": "CHROM POS ANN[0].GENE GEN[*].GT",
    }
    ctx = PipelineContext(args=Namespace(), config=config, workspace=ws)
    ctx.vcf_samples = ["SAMPLE1"]

    stage = ParallelCompleteProcessingStage()
    with pytest.raises(
        RuntimeError, match="Cannot restore parallel_complete_processing from checkpoint"
    ):
        stage._handle_checkpoint_skip(ctx)


def test_sites_only_vcf_with_empty_sample_list_omits_gen_columns(tmp_path: Path) -> None:
    """Verify sites-only VCF (vcf_samples=[]) omits GEN[*].GT columns and merges cleanly."""
    empty_vcf = tmp_path / "empty.vcf"
    pop_vcf = tmp_path / "pop.vcf"
    _write_minimal_vcf(empty_vcf, has_variants=False)
    _write_minimal_vcf(pop_vcf, has_variants=True)

    fields = "CHROM POS ANN[0].GENE GEN[*].GT"
    empty_tsv = tmp_path / "empty.tsv"
    pop_tsv = tmp_path / "pop.tsv"

    def mock_subprocess_run(cmd, *args, **kwargs):
        variant_file = cmd[-1]
        stdout = kwargs.get("stdout")
        if stdout and "pop.vcf" in variant_file:
            ann = "T|missense_variant|MODERATE|BRCA1|ENSG1|transcript|NM_1|protein_coding||||||||"
            stdout.write(f"chr1\t100\t{ann}\n")
        return Mock(returncode=0)

    with patch("variantcentrifuge.extractor.subprocess.run", side_effect=mock_subprocess_run):
        extract_fields_bcftools(
            variant_file=str(empty_vcf),
            fields=fields,
            cfg={},
            output_file=str(empty_tsv),
            vcf_samples=[],
        )
        extract_fields_bcftools(
            variant_file=str(pop_vcf),
            fields=fields,
            cfg={},
            output_file=str(pop_tsv),
            vcf_samples=[],
        )

    empty_header = empty_tsv.read_text().splitlines()[0].split("\t")
    pop_header = pop_tsv.read_text().splitlines()[0].split("\t")

    assert empty_header == pop_header == ["CHROM", "POS", "ANN[0].GENE"]
    assert "GEN[*].GT" not in empty_header

    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = tmp_path
    ws.get_intermediate_path = lambda name: tmp_path / name

    ctx = PipelineContext(
        args=Namespace(),
        config={"threads": 2, "gzip_intermediates": False, "fields_to_extract": fields},
        workspace=ws,
    )
    ctx.vcf_samples = []

    stage = ParallelCompleteProcessingStage()
    merged = stage._merge_tsv_outputs(ctx, [empty_tsv, pop_tsv])
    merged_lines = merged.read_text().splitlines()
    assert merged_lines[0].split("\t") == ["CHROM", "POS", "ANN[0].GENE"]
    assert len(merged_lines) == 2
