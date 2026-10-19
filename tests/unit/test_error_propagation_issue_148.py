"""Regression tests for Issue #148: Exception propagation across pipeline stages."""

from __future__ import annotations

import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from variantcentrifuge.cli import main
from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.chunk_helpers import (
    apply_chunk_annotations,
    apply_chunk_inheritance_analysis,
)
from variantcentrifuge.stages.inheritance_and_scoring_stages import InheritanceAnalysisStage
from variantcentrifuge.stages.phenotype_case_control_stages import (
    PhenotypeCaseControlAssignmentStage,
)
from variantcentrifuge.stages.setup_stages import PhenotypeLoadingStage


def test_malformed_csv_raises_during_assignment(tmp_path: Path):
    """Synthetic malformed CSV raises during assignment and does not mark complete."""
    workspace = Workspace(tmp_path / "work", "test")
    bad_csv = tmp_path / "bad_pheno.csv"
    bad_csv.write_text("sample,phenotype\nS1,case\nS2,control,unexpected_extra_field\n")

    config = {
        "phenotype_file": str(bad_csv),
        "phenotype_sample_column": "sample",
        "phenotype_value_column": "phenotype",
        "case_phenotypes": ["case"],
        "control_phenotypes": ["control"],
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.vcf_samples = ["S1", "S2"]

    # Phenotype loading succeeds (records file path), but assignment must raise and stay incomplete
    loading_stage = PhenotypeLoadingStage()
    loading_stage._process(context)
    context.mark_complete(loading_stage.name)
    assert context.is_complete("phenotype_loading")

    stage = PhenotypeCaseControlAssignmentStage()
    with pytest.raises((pd.errors.ParserError, ValueError, OSError)):
        stage._process(context)
    assert not context.is_complete("phenotype_case_control_assignment")
    assert not context.config.get("case_samples")


def test_case_control_file_read_exceptions_propagate(tmp_path: Path):
    """Case/control file read exceptions propagate with original causes."""
    workspace = Workspace(tmp_path / "work", "test")
    non_existent_case = tmp_path / "missing_cases.txt"
    config = {
        "case_samples_file": str(non_existent_case),
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.vcf_samples = ["S1", "S2"]
    stage = PhenotypeCaseControlAssignmentStage()
    with pytest.raises(FileNotFoundError):
        stage._process(context)


def test_inheritance_failure_propagates_in_non_chunked(tmp_path: Path):
    """Sequential/non-chunked inheritance failure propagates and is not swallowed."""
    workspace = Workspace(tmp_path / "work", "test")
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["BRCA1"],
            "GT": ["S1(0/1);S2(0/0)"],
        }
    )
    config = {
        "inheritance_mode": "simple",
        "use_chunked_processing": False,
        "threads": 1,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.current_dataframe = df
    context.vcf_samples = ["S1", "S2"]
    context.pedigree_data = {"S1": {"sample_id": "S1"}}

    stage = InheritanceAnalysisStage()
    with (
        patch(
            "variantcentrifuge.stages.analysis_stages.analyze_inheritance",
            side_effect=RuntimeError("Inheritance computation engine failure"),
        ),
        pytest.raises(RuntimeError, match="Inheritance computation engine failure"),
    ):
        stage._process(context)


def test_inheritance_failure_propagates_in_chunked(tmp_path: Path):
    """Chunked inheritance analysis failure propagates and is not replaced with error rows."""
    workspace = Workspace(tmp_path / "work", "test")
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["BRCA1"],
            "GT": ["S1(0/1);S2(0/0)"],
        }
    )
    config = {
        "inheritance_analysis_config": {
            "inheritance_mode": "simple",
            "vcf_samples": ["S1", "S2"],
            "pedigree_data": {},
        }
    }
    context = PipelineContext(Namespace(), config, workspace)

    with (
        patch(
            "variantcentrifuge.inheritance.analyzer.analyze_inheritance",
            side_effect=RuntimeError("Chunk inheritance failure"),
        ),
        pytest.raises(RuntimeError, match="Chunk inheritance failure"),
    ):
        apply_chunk_inheritance_analysis(df, context)


def test_chunk_annotation_failure_propagates(tmp_path: Path):
    """Chunk custom annotation failure propagates and stops chunk processing."""
    workspace = Workspace(tmp_path / "work", "test")
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["BRCA1"],
        }
    )
    context = PipelineContext(Namespace(), {}, workspace)
    context.annotation_configs = {"bed_files": ["/non/existent.bed"]}

    with (
        patch(
            "variantcentrifuge.annotator.load_custom_features",
            side_effect=ValueError("Feature load failed"),
        ),
        pytest.raises(ValueError, match="Feature load failed"),
    ):
        apply_chunk_annotations(df, context)


def test_parallel_inheritance_failure_propagates_and_not_checkpointed(tmp_path: Path):
    """Parallel inheritance failure (threads > 1) propagates and is not checkpointed."""
    workspace = Workspace(tmp_path / "work", "test")
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["BRCA1"],
            "GT": ["S1(0/1);S2(0/0)"],
        }
    )
    config = {
        "inheritance_mode": "simple",
        "threads": 2,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.current_dataframe = df
    context.vcf_samples = ["S1", "S2"]
    context.pedigree_data = {"S1": {"sample_id": "S1"}}

    stage = InheritanceAnalysisStage()
    with (
        patch(
            "variantcentrifuge.stages.analysis_stages.analyze_inheritance",
            side_effect=RuntimeError("Parallel inheritance worker failure"),
        ),
        pytest.raises(RuntimeError, match="Parallel inheritance worker failure"),
    ):
        stage._process(context)

    assert not context.is_complete("inheritance_analysis")


def test_custom_annotation_failure_propagates_both_modes(tmp_path: Path):
    """Custom annotation failures propagate in both non-chunked and chunked modes."""
    from variantcentrifuge.stages.analysis_stages import CustomAnnotationStage

    workspace = Workspace(tmp_path / "work", "test")
    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["BRCA1"],
        }
    )
    context = PipelineContext(Namespace(), {}, workspace)
    context.current_dataframe = df
    context.annotation_configs = {"bed_files": ["/fake.bed"]}

    stage = CustomAnnotationStage()
    with (
        patch(
            "variantcentrifuge.stages.analysis_stages.annotate_dataframe_with_features",
            side_effect=RuntimeError("Non-chunked annotation crash"),
        ),
        patch(
            "variantcentrifuge.stages.analysis_stages.load_custom_features",
            return_value=[{"type": "bed", "path": "/fake.bed"}],
        ),
        pytest.raises(RuntimeError, match="Non-chunked annotation crash"),
    ):
        stage._process(context)
    assert not context.is_complete("custom_annotation")

    with (
        patch(
            "variantcentrifuge.annotator.annotate_dataframe_with_features",
            side_effect=RuntimeError("Chunk annotation crash"),
        ),
        patch(
            "variantcentrifuge.annotator.load_custom_features",
            return_value=[{"type": "bed", "path": "/fake.bed"}],
        ),
        pytest.raises(RuntimeError, match="Chunk annotation crash"),
    ):
        apply_chunk_annotations(df, context)


def test_cli_returns_nonzero_on_stage_failure(tmp_path: Path):
    """CLI execution returns code 1 when pipeline stage raises an exception."""
    test_vcf = tmp_path / "test.vcf"
    test_vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with (
        patch.object(
            sys,
            "argv",
            [
                "variantcentrifuge",
                "--vcf",
                str(test_vcf),
                "--output-dir",
                str(tmp_path),
                "--filters",
                "PASS",
            ],
        ),
        patch(
            "variantcentrifuge.cli.run_pipeline",
            side_effect=RuntimeError("Injected pipeline crash"),
        ),
    ):
        exit_code = main()
        assert exit_code == 1
