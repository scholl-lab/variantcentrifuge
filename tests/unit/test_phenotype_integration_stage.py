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
