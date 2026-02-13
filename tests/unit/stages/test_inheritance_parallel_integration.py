"""Test InheritanceAnalysisStage integration with parallel analyzer."""

import logging
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import InheritanceAnalysisStage


class TestInheritanceParallelIntegration:
    """Test that InheritanceAnalysisStage correctly uses parallel analyzer."""

    @pytest.fixture
    def context(self, tmp_path):
        """Create test context with multiple threads."""
        workspace = Workspace(tmp_path, "test")
        args = Mock(
            threads=4,  # Multiple threads to trigger parallel
            pedigree_file=tmp_path / "test.ped",
        )

        # Create test pedigree file
        ped_content = (
            "#FamilyID\tIndividualID\tPaternalID\tMaternalID\tSex\tPhenotype\n"
            "FAM1\tproband\tfather\tmother\t1\t2\n"
            "FAM1\tmother\t0\t0\t2\t1\n"
            "FAM1\tfather\t0\t0\t1\t1\n"
        )
        args.pedigree_file.write_text(ped_content)

        config = {
            "threads": 4,
            "pedigree_file": str(args.pedigree_file),
            "inheritance_mode": "columns",
            "min_variants_for_parallel_inheritance": 50,  # Low threshold for testing
        }

        ctx = PipelineContext(args=args, config=config, workspace=workspace)

        # Create sample DataFrame with enough variants
        data = {
            "CHROM": ["chr1"] * 100,
            "POS": list(range(1000, 1100)),
            "REF": ["A"] * 100,
            "ALT": ["T"] * 100,
            "GENE": (
                ["GENE1"] * 20 + ["GENE2"] * 20 + ["GENE3"] * 20 + ["GENE4"] * 20 + ["GENE5"] * 20
            ),
            "IMPACT": ["MODERATE"] * 100,
        }
        ctx.current_dataframe = pd.DataFrame(data)

        # Mark dependencies complete
        ctx.mark_complete("dataframe_loading")
        ctx.mark_complete("pedigree_loading")

        # Set VCF samples
        ctx.vcf_samples = ["proband", "mother", "father"]

        # Load pedigree data into context
        from variantcentrifuge.ped_reader import read_pedigree

        ctx.pedigree_data = read_pedigree(str(args.pedigree_file))

        return ctx

    def test_uses_parallel_analyzer_when_appropriate(self, context, caplog):
        """Test that parallel analyzer is used with multiple threads and large dataset."""
        stage = InheritanceAnalysisStage()

        with caplog.at_level(logging.INFO):
            # Mock the parallel analyzer to verify it's called
            with patch(
                "variantcentrifuge.inheritance.parallel_analyzer.analyze_inheritance_parallel"
            ) as mock_parallel:
                # Make it return a valid DataFrame
                mock_parallel.return_value = context.current_dataframe.copy()
                mock_parallel.return_value["Inheritance_Pattern"] = "autosomal_dominant"
                mock_parallel.return_value["Inheritance_Details"] = "{}"

                stage(context)

                # Verify parallel analyzer was called
                mock_parallel.assert_called_once()
                call_args = mock_parallel.call_args

                # Check arguments
                assert call_args.kwargs["n_workers"] == 4
                assert call_args.kwargs["min_variants_for_parallel"] == 50
                assert call_args.kwargs["use_vectorized_comp_het"] is True  # Default is True

        # Check logs
        assert "Using parallel inheritance analyzer with 4 workers" in caplog.text

    def test_uses_sequential_analyzer_for_small_dataset(self, context, caplog):
        """Test that sequential analyzer is used for small datasets."""
        # Reduce dataset size
        context.current_dataframe = context.current_dataframe.iloc[:10]  # Only 10 variants

        stage = InheritanceAnalysisStage()

        with caplog.at_level(logging.INFO):
            # Mock the sequential analyzer where it's imported
            with patch("variantcentrifuge.stages.analysis_stages.analyze_inheritance") as mock_seq:
                # Make it return a valid DataFrame
                mock_seq.return_value = context.current_dataframe.copy()
                mock_seq.return_value["Inheritance_Pattern"] = "autosomal_dominant"
                mock_seq.return_value["Inheritance_Details"] = "{}"

                stage(context)

                # Verify sequential analyzer was called
                mock_seq.assert_called_once()

        # Should not mention parallel analyzer in logs
        assert "Using parallel inheritance analyzer" not in caplog.text

    def test_uses_sequential_analyzer_with_single_thread(self, context, caplog):
        """Test that sequential analyzer is used with single thread."""
        # Set to single thread
        context.config["threads"] = 1

        stage = InheritanceAnalysisStage()

        with caplog.at_level(logging.INFO):
            # Mock the sequential analyzer where it's imported
            with patch("variantcentrifuge.stages.analysis_stages.analyze_inheritance") as mock_seq:
                # Make it return a valid DataFrame
                mock_seq.return_value = context.current_dataframe.copy()
                mock_seq.return_value["Inheritance_Pattern"] = "autosomal_dominant"
                mock_seq.return_value["Inheritance_Details"] = "{}"

                stage(context)

                # Verify sequential analyzer was called
                mock_seq.assert_called_once()

        # Should not mention parallel analyzer in logs
        assert "Using parallel inheritance analyzer" not in caplog.text

    def test_subtask_timing_for_compound_het(self, context):
        """Test that compound het timing is recorded when significant."""
        stage = InheritanceAnalysisStage()

        # Mock the parallel analyzer to simulate slow compound het analysis
        with patch(
            "variantcentrifuge.inheritance.parallel_analyzer.analyze_inheritance_parallel"
        ) as mock_parallel:
            import time

            def slow_analysis(*args, **kwargs):
                # Simulate 2 second compound het analysis
                time.sleep(2)
                df = context.current_dataframe.copy()
                df["Inheritance_Pattern"] = "compound_heterozygous"
                df["Inheritance_Details"] = "{}"
                return df

            mock_parallel.side_effect = slow_analysis

            stage(context)

        # Check that compound het timing was recorded
        assert hasattr(stage, "_subtask_times")
        assert "compound_het_analysis" in stage._subtask_times
        # Use a tolerance for timing assertions to account for system variations
        assert stage._subtask_times["compound_het_analysis"] >= 1.8

    def test_configuration_options_passed_correctly(self, context):
        """Test that configuration options are passed to parallel analyzer."""
        # Set specific config options
        context.config["no_vectorized_comp_het"] = True
        context.config["min_variants_for_parallel_inheritance"] = 75

        stage = InheritanceAnalysisStage()

        with patch(
            "variantcentrifuge.inheritance.parallel_analyzer.analyze_inheritance_parallel"
        ) as mock_parallel:
            mock_parallel.return_value = context.current_dataframe.copy()
            mock_parallel.return_value["Inheritance_Pattern"] = "autosomal_dominant"
            mock_parallel.return_value["Inheritance_Details"] = "{}"

            stage(context)

            # Check that config was passed correctly
            call_args = mock_parallel.call_args
            assert (
                call_args.kwargs["use_vectorized_comp_het"] is False
            )  # Inverted from no_vectorized_comp_het
            assert call_args.kwargs["min_variants_for_parallel"] == 75

    def test_all_subtasks_timed(self, context):
        """Test that all subtasks are properly timed."""
        stage = InheritanceAnalysisStage()

        with patch(
            "variantcentrifuge.inheritance.parallel_analyzer.analyze_inheritance_parallel"
        ) as mock_parallel:
            mock_parallel.return_value = context.current_dataframe.copy()
            mock_parallel.return_value["Inheritance_Pattern"] = "autosomal_dominant"
            mock_parallel.return_value["Inheritance_Details"] = "{}"

            stage(context)

        # Check all expected subtasks are timed
        # Note: sample_column_preparation is only done if GT column exists and
        # samples not already in DataFrame
        expected_subtasks = [
            "inheritance_calculation",
            "output_processing",
            "column_cleanup",
        ]

        for subtask in expected_subtasks:
            assert subtask in stage._subtask_times
            assert stage._subtask_times[subtask] >= 0
