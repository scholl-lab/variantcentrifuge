"""Unit tests for PhenotypeIntegrationStage."""

import gzip
from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.post_processing_stages import PhenotypeIntegrationStage


def _setup_context(
    tmp_path: Path, gzip_intermediates: bool = False
) -> tuple[PipelineContext, Path]:
    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir(parents=True, exist_ok=True)
    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = inter_dir
    ws.get_intermediate_path = lambda name: inter_dir / name

    ext = ".tsv.gz" if gzip_intermediates else ".tsv"
    extracted_tsv = inter_dir / f"sample.extracted.sorted{ext}"

    tsv_data = (
        "CHROM\tPOS\tANN[0].GENE\tGEN[0].GT\tGEN[1].GT\n"
        "chr1\t100\tGENE1\t0/1\t0/0\n"
        "chr2\t200\tGENE2\t0/0\t1/1\n"
    )
    if gzip_intermediates:
        with gzip.open(extracted_tsv, "wt") as f:
            f.write(tsv_data)
    else:
        extracted_tsv.write_text(tsv_data)

    config = {
        "gzip_intermediates": gzip_intermediates,
        "no_replacement": True,
    }
    context = PipelineContext(args=Namespace(no_replacement=True), config=config, workspace=ws)
    context.extracted_tsv = extracted_tsv
    context.data = extracted_tsv
    context.vcf_samples = ["SampleA", "SampleB"]
    context.phenotype_data = {
        "SampleA": {"HPO:0001"},
        "SampleB": {"HPO:0002"},
    }
    return context, extracted_tsv


def test_phenotype_stage_soft_dependencies() -> None:
    """Verify PhenotypeIntegrationStage depends softly on data_sorting."""
    stage = PhenotypeIntegrationStage()
    assert "data_sorting" in stage.soft_dependencies


@pytest.mark.parametrize("gzip_intermediates", [False, True])
def test_phenotype_stage_processes_raw_indexed_genotype_columns(
    tmp_path: Path, gzip_intermediates: bool
) -> None:
    """Verify PhenotypeIntegrationStage recognizes GEN[i].GT columns and adds Phenotypes."""
    context, _ = _setup_context(tmp_path, gzip_intermediates=gzip_intermediates)
    stage = PhenotypeIntegrationStage()

    result_ctx = stage._process(context)

    assert result_ctx.phenotypes_added_tsv is not None
    assert result_ctx.phenotypes_added_tsv.exists()
    assert result_ctx.data == result_ctx.phenotypes_added_tsv

    compression = "gzip" if gzip_intermediates else None
    df = pd.read_csv(result_ctx.phenotypes_added_tsv, sep="\t", compression=compression)

    assert "Phenotypes" in df.columns
    # Row 0: SampleA is 0/1 (variant), SampleB is 0/0 (ref) -> SampleA(HPO:0001)
    assert df.loc[0, "Phenotypes"] == "SampleA(HPO:0001)"
    # Row 1: SampleA is 0/0 (ref), SampleB is 1/1 (variant) -> SampleB(HPO:0002)
    assert df.loc[1, "Phenotypes"] == "SampleB(HPO:0002)"


@pytest.mark.parametrize("gzip_intermediates", [False, True])
def test_pipeline_chain_sorting_phenotype_dataframe_loading(
    tmp_path: Path, gzip_intermediates: bool
) -> None:
    """Verify DataSortingStage -> PhenotypeIntegrationStage -> DataFrameLoadingStage chain."""
    from variantcentrifuge.stages.analysis_stages import DataFrameLoadingStage
    from variantcentrifuge.stages.processing_stages import DataSortingStage

    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir(parents=True, exist_ok=True)
    ws = Mock(spec=Workspace)
    ws.base_name = "sample"
    ws.intermediate_dir = inter_dir
    ws.output_dir = tmp_path
    ws.get_intermediate_path = lambda name: inter_dir / name

    ext = ".tsv.gz" if gzip_intermediates else ".tsv"
    extracted_tsv = inter_dir / f"sample.extracted{ext}"

    # Unsorted rows (ZZZ before AAA)
    tsv_data = (
        "CHROM\tPOS\tANN[0].GENE\tGEN[0].GT\tGEN[1].GT\n"
        "chr2\t200\tZZZ\t0/0\t1/1\n"
        "chr1\t100\tAAA\t0/1\t0/0\n"
    )
    if gzip_intermediates:
        with gzip.open(extracted_tsv, "wt") as f:
            f.write(tsv_data)
    else:
        extracted_tsv.write_text(tsv_data)

    config = {
        "gzip_intermediates": gzip_intermediates,
        "no_replacement": True,
        "sort_gene": True,
        "fields_to_extract": "CHROM POS ANN[0].GENE GEN[*].GT",
    }
    context = PipelineContext(args=Namespace(no_replacement=True), config=config, workspace=ws)
    context.extracted_tsv = extracted_tsv
    context.data = extracted_tsv
    context.vcf_samples = ["SampleA", "SampleB"]
    context.phenotype_data = {
        "SampleA": {"HPO:0001"},
        "SampleB": {"HPO:0002"},
    }

    # 1. Run DataSortingStage
    sorting_stage = DataSortingStage()
    ctx_sorted = sorting_stage._process(context)
    assert ctx_sorted.data.exists()

    # 2. Run PhenotypeIntegrationStage
    phenotype_stage = PhenotypeIntegrationStage()
    ctx_pheno = phenotype_stage._process(ctx_sorted)
    assert ctx_pheno.phenotypes_added_tsv.exists()

    # 3. Run DataFrameLoadingStage
    loading_stage = DataFrameLoadingStage()
    ctx_loaded = loading_stage._process(ctx_pheno)

    df = ctx_loaded.current_dataframe
    assert df is not None
    assert len(df) == 2

    # Assert rows are sorted alphabetically by gene (AAA, then ZZZ)
    gene_col = "GENE" if "GENE" in df.columns else "ANN_0__GENE"
    assert list(df[gene_col]) == ["AAA", "ZZZ"]

    # Assert Phenotypes column exists and matches sorted order
    assert list(df["Phenotypes"]) == ["SampleA(HPO:0001)", "SampleB(HPO:0002)"]
