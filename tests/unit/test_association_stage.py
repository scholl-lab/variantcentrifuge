"""
Unit tests for AssociationAnalysisStage, including coexistence tests with
GeneBurdenAnalysisStage (CORE-05 validation).

Tests stage skip/run behavior, config guards, output storage in context,
and full independence between the two analysis stages.
"""

from __future__ import annotations

from argparse import Namespace
from unittest.mock import Mock

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import (
    AssociationAnalysisStage,
    GeneBurdenAnalysisStage,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_workspace(tmp_path):
    """Create a mock workspace pointing to tmp_path."""
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    workspace.intermediate_dir = tmp_path / "intermediate"
    workspace.intermediate_dir.mkdir()
    workspace.get_intermediate_path = lambda x: workspace.intermediate_dir / x
    workspace.get_output_path = lambda x, ext=".tsv": workspace.output_dir / f"{x}{ext}"
    workspace.base_name = "test"
    return workspace


@pytest.fixture
def minimal_df():
    """A minimal DataFrame with GENE and GT columns for analysis stages."""
    return pd.DataFrame(
        {
            "GENE": ["BRCA1", "BRCA1", "TP53"],
            "GT": [
                "CASE1(0/1);CASE2(0/1);CTRL1(0/0)",
                "CASE1(0/1);CTRL2(0/1)",
                "CASE1(1/1);CTRL1(0/1)",
            ],
        }
    )


@pytest.fixture
def case_control_config(tmp_path):
    """Config dict with case/control samples and analysis flags.

    Uses 10 cases and 10 controls to pass the tiered sample size guard
    (Phase 19: cases < 10 refuses to run).
    """
    cases = [f"CASE{i}" for i in range(1, 11)]
    controls = [f"CTRL{i}" for i in range(1, 11)]
    return {
        "perform_association": True,
        "perform_gene_burden": False,
        "case_samples": cases,
        "control_samples": controls,
        "gene_burden_mode": "samples",
        "correction_method": "fdr",
        "association_tests": ["fisher"],
        "output_dir": str(tmp_path),
        "output_file_base": "test_output",
        "gzip_intermediates": False,
    }


def _make_context(config: dict, df: pd.DataFrame | None, workspace) -> PipelineContext:
    """Construct a PipelineContext with the given config and DataFrame."""
    context = PipelineContext(
        args=Namespace(),
        config=config,
        workspace=workspace,
    )
    if df is not None:
        context.current_dataframe = df
    return context


# ---------------------------------------------------------------------------
# AssociationAnalysisStage: skip/run behavior
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestAssociationAnalysisStageSkipBehavior:
    """AssociationAnalysisStage skips when perform_association is not set or False."""

    def test_skips_when_perform_association_not_in_config(self, mock_workspace, minimal_df):
        """Stage returns context unchanged when perform_association is absent."""
        config = {}  # no perform_association key
        context = _make_context(config, minimal_df, mock_workspace)

        stage = AssociationAnalysisStage()
        result = stage._process(context)

        # Context unchanged: no association_results set
        assert result.association_results is None

    def test_skips_when_perform_association_is_false(self, mock_workspace, minimal_df):
        """Stage returns context unchanged when perform_association=False."""
        config = {"perform_association": False}
        context = _make_context(config, minimal_df, mock_workspace)

        stage = AssociationAnalysisStage()
        result = stage._process(context)

        assert result.association_results is None

    def test_skips_when_no_case_samples(self, mock_workspace, minimal_df):
        """Stage returns context unchanged when case_samples is empty."""
        config = {
            "perform_association": True,
            "case_samples": [],
            "control_samples": ["CTRL1"],
        }
        context = _make_context(config, minimal_df, mock_workspace)

        stage = AssociationAnalysisStage()
        result = stage._process(context)

        assert result.association_results is None

    def test_skips_when_no_control_samples(self, mock_workspace, minimal_df):
        """Stage returns context unchanged when control_samples is empty."""
        config = {
            "perform_association": True,
            "case_samples": ["CASE1"],
            "control_samples": [],
        }
        context = _make_context(config, minimal_df, mock_workspace)

        stage = AssociationAnalysisStage()
        result = stage._process(context)

        assert result.association_results is None

    def test_skips_when_dataframe_is_none(self, mock_workspace):
        """Stage returns context unchanged when no DataFrame is loaded."""
        config = {
            "perform_association": True,
            "case_samples": ["CASE1"],
            "control_samples": ["CTRL1"],
        }
        context = _make_context(config, None, mock_workspace)

        stage = AssociationAnalysisStage()
        result = stage._process(context)

        assert result.association_results is None


# ---------------------------------------------------------------------------
# AssociationAnalysisStage: stage metadata
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestAssociationAnalysisStageMetadata:
    """Stage name, dependencies, and soft_dependencies are correct."""

    def test_stage_name_is_association_analysis(self):
        """Stage.name == 'association_analysis'."""
        stage = AssociationAnalysisStage()
        assert stage.name == "association_analysis"

    def test_dependencies_include_required_stages(self):
        """Stage.dependencies includes dataframe_loading and sample_config_loading."""
        stage = AssociationAnalysisStage()
        assert "dataframe_loading" in stage.dependencies
        assert "sample_config_loading" in stage.dependencies

    def test_soft_dependencies_include_custom_annotation(self):
        """Stage.soft_dependencies includes custom_annotation."""
        stage = AssociationAnalysisStage()
        assert "custom_annotation" in stage.soft_dependencies


# ---------------------------------------------------------------------------
# AssociationAnalysisStage: successful run
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestAssociationAnalysisStageRun:
    """AssociationAnalysisStage produces results when properly configured."""

    def test_sets_association_results_on_context(
        self, mock_workspace, minimal_df, case_control_config
    ):
        """Stage sets context.association_results to a non-empty DataFrame."""
        context = _make_context(case_control_config, minimal_df, mock_workspace)

        stage = AssociationAnalysisStage()
        result = stage._process(context)

        assert result.association_results is not None
        assert isinstance(result.association_results, pd.DataFrame)
        assert len(result.association_results) > 0

    def test_sets_association_output_path_in_config(
        self, mock_workspace, minimal_df, case_control_config
    ):
        """Stage sets context.config['association_output'] after successful run."""
        context = _make_context(case_control_config, minimal_df, mock_workspace)

        stage = AssociationAnalysisStage()
        result = stage._process(context)

        assert "association_output" in result.config
        assert result.config["association_output"] is not None

    def test_association_results_has_expected_columns(
        self, mock_workspace, minimal_df, case_control_config
    ):
        """association_results DataFrame has gene, fisher_p_value, fisher_or columns."""
        context = _make_context(case_control_config, minimal_df, mock_workspace)

        stage = AssociationAnalysisStage()
        result = stage._process(context)

        df = result.association_results
        assert "gene" in df.columns
        assert "fisher_p_value" in df.columns
        assert "fisher_or" in df.columns

    def test_association_config_key_used_not_perform_gene_burden(
        self, mock_workspace, minimal_df, tmp_path
    ):
        """Stage guard reads perform_association, not perform_gene_burden."""
        cases = [f"CASE{i}" for i in range(1, 11)]
        controls = [f"CTRL{i}" for i in range(1, 11)]
        config = {
            "perform_association": True,
            "perform_gene_burden": False,  # explicitly False; should NOT block this stage
            "case_samples": cases,
            "control_samples": controls,
            "gene_burden_mode": "samples",
            "correction_method": "fdr",
            "association_tests": ["fisher"],
            "output_dir": str(tmp_path),
            "output_file_base": "test_output",
            "gzip_intermediates": False,
        }
        context = _make_context(config, minimal_df, mock_workspace)

        stage = AssociationAnalysisStage()
        result = stage._process(context)

        # Association runs even though perform_gene_burden=False
        assert result.association_results is not None


# ---------------------------------------------------------------------------
# CORE-05 Coexistence tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestStageCoexistence:
    """
    CORE-05: GeneBurdenAnalysisStage and AssociationAnalysisStage are fully
    independent. Both can run in the same pipeline invocation.
    """

    def test_both_stages_coexist_no_name_conflict(self):
        """Both stages have distinct names and can be in the same stage list."""
        gb_stage = GeneBurdenAnalysisStage()
        assoc_stage = AssociationAnalysisStage()

        assert gb_stage.name != assoc_stage.name
        assert gb_stage.name == "gene_burden_analysis"
        assert assoc_stage.name == "association_analysis"

    def test_both_stages_independent_dependencies(self):
        """Neither stage declares the other as a dependency."""
        gb_stage = GeneBurdenAnalysisStage()
        assoc_stage = AssociationAnalysisStage()

        # Neither should depend on the other
        assert assoc_stage.name not in gb_stage.dependencies
        assert assoc_stage.name not in gb_stage.soft_dependencies
        assert gb_stage.name not in assoc_stage.dependencies
        assert gb_stage.name not in assoc_stage.soft_dependencies

    def test_gene_burden_stage_checks_perform_gene_burden_not_perform_association(
        self, mock_workspace, minimal_df
    ):
        """GeneBurdenAnalysisStage skips when perform_gene_burden=False.

        The perform_association config key must have no effect on this stage.
        """
        config = {
            "perform_gene_burden": False,
            "perform_association": True,  # should have no effect on GB stage
            "case_samples": ["CASE1"],
            "control_samples": ["CTRL1"],
        }
        context = _make_context(config, minimal_df, mock_workspace)

        gb_stage = GeneBurdenAnalysisStage()
        result = gb_stage._process(context)

        # Gene burden must skip regardless of perform_association
        assert result.gene_burden_results is None

    def test_association_stage_checks_perform_association_not_perform_gene_burden(
        self, mock_workspace, minimal_df
    ):
        """AssociationAnalysisStage skips when perform_association=False.

        The perform_gene_burden config key must have no effect on this stage.
        """
        config = {
            "perform_gene_burden": True,  # should have no effect on association stage
            "perform_association": False,
            "case_samples": ["CASE1"],
            "control_samples": ["CTRL1"],
        }
        context = _make_context(config, minimal_df, mock_workspace)

        assoc_stage = AssociationAnalysisStage()
        result = assoc_stage._process(context)

        # Association must skip regardless of perform_gene_burden
        assert result.association_results is None

    def test_both_stages_run_simultaneously_perform_both_true(
        self, mock_workspace, minimal_df, tmp_path
    ):
        """Both stages produce results when both perform_ flags are True."""
        cases = [f"CASE{i}" for i in range(1, 11)]
        controls = [f"CTRL{i}" for i in range(1, 11)]
        config = {
            "perform_gene_burden": True,
            "perform_association": True,
            "case_samples": cases,
            "control_samples": controls,
            "gene_burden_mode": "samples",
            "correction_method": "fdr",
            "association_tests": ["fisher"],
            "output_dir": str(tmp_path),
            "output_file_base": "test_output",
            "gzip_intermediates": False,
        }
        context = _make_context(config, minimal_df, mock_workspace)

        gb_stage = GeneBurdenAnalysisStage()
        assoc_stage = AssociationAnalysisStage()

        # Run gene burden first, then association
        context = gb_stage._process(context)
        context = assoc_stage._process(context)

        # Both should have results
        assert context.gene_burden_results is not None
        assert context.association_results is not None

    def test_both_stages_do_not_share_output_config_keys(
        self, mock_workspace, minimal_df, tmp_path
    ):
        """Gene burden and association use distinct config keys for output paths."""
        cases = [f"CASE{i}" for i in range(1, 11)]
        controls = [f"CTRL{i}" for i in range(1, 11)]
        config = {
            "perform_gene_burden": True,
            "perform_association": True,
            "case_samples": cases,
            "control_samples": controls,
            "gene_burden_mode": "samples",
            "correction_method": "fdr",
            "association_tests": ["fisher"],
            "output_dir": str(tmp_path),
            "output_file_base": "test_output",
            "gzip_intermediates": False,
        }
        context = _make_context(config, minimal_df, mock_workspace)

        gb_stage = GeneBurdenAnalysisStage()
        assoc_stage = AssociationAnalysisStage()

        context = gb_stage._process(context)
        context = assoc_stage._process(context)

        # Each stage uses its own config key
        gb_output = context.config.get("gene_burden_output")
        assoc_output = context.config.get("association_output")

        assert gb_output is not None
        assert assoc_output is not None
        assert gb_output != assoc_output  # Different output files

    def test_only_gene_burden_true_association_skips(self, mock_workspace, minimal_df, tmp_path):
        """When only perform_gene_burden=True, AssociationAnalysisStage skips cleanly."""
        config = {
            "perform_gene_burden": True,
            "perform_association": False,
            "case_samples": ["CASE1", "CASE2"],
            "control_samples": ["CTRL1", "CTRL2"],
            "gene_burden_mode": "samples",
            "correction_method": "fdr",
            "output_dir": str(tmp_path),
            "output_file_base": "test_output",
            "gzip_intermediates": False,
        }
        context = _make_context(config, minimal_df, mock_workspace)

        gb_stage = GeneBurdenAnalysisStage()
        assoc_stage = AssociationAnalysisStage()

        context = gb_stage._process(context)
        context = assoc_stage._process(context)

        assert context.gene_burden_results is not None
        assert context.association_results is None

    def test_only_association_true_gene_burden_skips(self, mock_workspace, minimal_df, tmp_path):
        """When only perform_association=True, GeneBurdenAnalysisStage skips cleanly."""
        cases = [f"CASE{i}" for i in range(1, 11)]
        controls = [f"CTRL{i}" for i in range(1, 11)]
        config = {
            "perform_gene_burden": False,
            "perform_association": True,
            "case_samples": cases,
            "control_samples": controls,
            "gene_burden_mode": "samples",
            "correction_method": "fdr",
            "association_tests": ["fisher"],
            "output_dir": str(tmp_path),
            "output_file_base": "test_output",
            "gzip_intermediates": False,
        }
        context = _make_context(config, minimal_df, mock_workspace)

        gb_stage = GeneBurdenAnalysisStage()
        assoc_stage = AssociationAnalysisStage()

        context = gb_stage._process(context)
        context = assoc_stage._process(context)

        assert context.gene_burden_results is None
        assert context.association_results is not None

    def test_gene_burden_stage_does_not_inspect_perform_association_key(
        self, mock_workspace, minimal_df
    ):
        """GeneBurdenAnalysisStage source code does not reference perform_association."""
        import inspect

        from variantcentrifuge.stages.analysis_stages import GeneBurdenAnalysisStage

        # Get the _process method source code
        source = inspect.getsource(GeneBurdenAnalysisStage._process)

        # The guard must be perform_gene_burden, not perform_association
        assert "perform_gene_burden" in source
        # It must not be gated on perform_association in the _process logic
        # (it may appear in comments, but not as a config.get() guard)
        assert 'config.get("perform_association")' not in source
