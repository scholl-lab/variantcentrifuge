"""Unit tests for schema preservation on empty BEDs and empty variant chunks (#149)."""

from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock

import pytest

from variantcentrifuge.extractor import extract_fields_bcftools
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import DataFrameLoadingStage
from variantcentrifuge.stages.processing_stages import ParallelCompleteProcessingStage


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

    fields = "CHROM POS REF ALT ANN[0].GENE ANN[0].EFFECT NMD[0].PERC GEN[*].GT"
    empty_tsv = tmp_path / "empty.tsv"
    pop_tsv = tmp_path / "pop.tsv"

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

    # Empty BED file
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

    # DataFrameLoadingStage should load this file without EmptyDataError
    loading_stage = DataFrameLoadingStage()
    loaded_ctx = loading_stage._process(updated_ctx)

    assert loaded_ctx.current_dataframe is not None
    assert len(loaded_ctx.current_dataframe) == 0
    assert "CHROM" in loaded_ctx.current_dataframe.columns
    cols = loaded_ctx.current_dataframe.columns
    assert "GEN[0].GT" in cols or "GEN_0__GT" in cols


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

    # chunk0 has fewer columns
    chunk0.write_text("CHROM\tPOS\tANN\n")
    chunk1.write_text("CHROM\tPOS\tANN[0].GENE\tGEN[0].GT\nchr1\t100\tBRCA1\t0/1\n")

    stage = ParallelCompleteProcessingStage()
    with pytest.raises(ValueError, match="Incompatible chunk header"):
        stage._merge_tsvs_traditional([chunk0, chunk1], inter_dir / "merged.tsv")
