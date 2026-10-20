"""Tests for Issue #150: Missing cohort rejection and TSV fallback fail-on-empty handling."""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.pipeline import _validate_requested_outputs
from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.association_analysis_stages import AssociationAnalysisStage
from variantcentrifuge.stages.gene_burden_stages import GeneBurdenAnalysisStage
from variantcentrifuge.stages.tsv_output_stages import TSVOutputStage


def test_missing_cohorts_rejected_for_gene_burden(tmp_path: Path):
    """Missing cohorts produce explicit ValueError when gene burden is requested."""
    workspace = Workspace(tmp_path / "work", "test")
    stage = GeneBurdenAnalysisStage()

    # Both missing
    ctx_both = PipelineContext(
        Namespace(),
        {"perform_gene_burden": True, "case_samples": [], "control_samples": []},
        workspace,
    )
    with pytest.raises(ValueError, match="both case and control cohorts are missing"):
        stage._process(ctx_both)

    # Missing cases
    ctx_no_cases = PipelineContext(
        Namespace(),
        {"perform_gene_burden": True, "case_samples": [], "control_samples": ["CTRL1"]},
        workspace,
    )
    with pytest.raises(ValueError, match="case cohort is missing"):
        stage._process(ctx_no_cases)

    # Missing controls
    ctx_no_controls = PipelineContext(
        Namespace(),
        {"perform_gene_burden": True, "case_samples": ["CASE1"], "control_samples": []},
        workspace,
    )
    with pytest.raises(ValueError, match="control cohort is missing"):
        stage._process(ctx_no_controls)

    # Disabled burden skips without error
    ctx_disabled = PipelineContext(Namespace(), {"perform_gene_burden": False}, workspace)
    stage._process(ctx_disabled)


def test_missing_cohorts_rejected_for_association(tmp_path: Path):
    """Missing cohorts produce explicit ValueError when association is requested."""
    workspace = Workspace(tmp_path / "work", "test")
    stage = AssociationAnalysisStage()

    # Both missing
    ctx_both = PipelineContext(
        Namespace(),
        {"perform_association": True, "case_samples": [], "control_samples": []},
        workspace,
    )
    with pytest.raises(ValueError, match="both case and control cohorts are missing"):
        stage._process(ctx_both)

    # Missing cases
    ctx_no_cases = PipelineContext(
        Namespace(),
        {"perform_association": True, "case_samples": [], "control_samples": ["CTRL1"]},
        workspace,
    )
    with pytest.raises(ValueError, match="case cohort is missing"):
        stage._process(ctx_no_cases)

    # Missing controls
    ctx_no_controls = PipelineContext(
        Namespace(),
        {"perform_association": True, "case_samples": ["CASE1"], "control_samples": []},
        workspace,
    )
    with pytest.raises(ValueError, match="control cohort is missing"):
        stage._process(ctx_no_controls)

    # Disabled association skips without error
    ctx_disabled = PipelineContext(Namespace(), {"perform_association": False}, workspace)
    stage._process(ctx_disabled)


def test_zero_variant_cohort_produces_sidecars_validating_outputs(tmp_path: Path):
    """Valid zero-variant cohorts produce header-only sidecars that pass output validation."""
    workspace = Workspace(tmp_path / "work", "test")
    workspace.output_dir.mkdir(parents=True, exist_ok=True)

    burden_out = str(workspace.output_dir / "burden.tsv")
    assoc_out = str(workspace.output_dir / "assoc.tsv")
    final_out = workspace.output_dir / "final.tsv"
    final_out.write_text("CHROM\tPOS\tREF\tALT\n")

    config = {
        "output_file": str(final_out),
        "case_samples": ["CASE1"],
        "control_samples": ["CTRL1"],
        "perform_gene_burden": True,
        "gene_burden_output": burden_out,
        "perform_association": True,
        "association_output": assoc_out,
        "gzip_intermediates": False,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.current_dataframe = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "GENE", "GT"])
    context.final_output_path = final_out

    # Process gene burden on empty df -> writes header-only sidecar
    GeneBurdenAnalysisStage()._process(context)
    assert Path(burden_out).exists()
    assert pd.read_csv(burden_out, sep="\t").empty

    # Process association on empty df -> writes header-only sidecar
    AssociationAnalysisStage()._process(context)
    assert Path(assoc_out).exists()
    assert pd.read_csv(assoc_out, sep="\t").empty

    # Output validation must pass without raising RuntimeError
    _validate_requested_outputs(context)


def test_restored_header_only_checkpoint_fail_on_empty_exits_2(tmp_path: Path):
    """Restored header-only checkpoint DataFrame exits 2 under fail_on_empty=True."""
    workspace = Workspace(tmp_path / "work", "test")
    config = {
        "output_file": str(tmp_path / "final.tsv"),
        "fail_on_empty": True,
    }
    context = PipelineContext(Namespace(), config, workspace)
    context.current_dataframe = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT"])

    stage = TSVOutputStage()
    with pytest.raises(SystemExit) as excinfo:
        stage._process(context)
    assert excinfo.value.code == 2
