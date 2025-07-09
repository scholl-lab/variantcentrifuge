"""Unit tests for InheritanceAnalysisStage."""

from unittest.mock import patch

import pandas as pd
import pytest

from tests.mocks.fixtures import create_test_context
from variantcentrifuge.stages.analysis_stages import InheritanceAnalysisStage


class TestInheritanceAnalysisStage:
    """Test InheritanceAnalysisStage."""

    @pytest.fixture
    def context(self):
        """Create test context with DataFrame and sample data."""
        ctx = create_test_context()
        # Create a test DataFrame with GT column
        ctx.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr2", "chr3"],
                "POS": [100, 200, 300, 400],
                "REF": ["A", "C", "G", "T"],
                "ALT": ["G", "T", "A", "C"],
                "QUAL": [100, 90, 85, 95],
                "GENE": ["GENE1", "GENE2", "GENE2", "GENE3"],
                "GT": ["0/1,0/0,0/0", "0/1,0/1,0/0", "0/1,0/0,0/1", "1/1,0/1,0/1"],
            }
        )
        ctx.vcf_samples = ["Child", "Father", "Mother"]
        ctx.pedigree_data = {
            "Father": {
                "family_id": "FAM001",
                "sample_id": "Father",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "1",
            },
            "Mother": {
                "family_id": "FAM001",
                "sample_id": "Mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
            "Child": {
                "family_id": "FAM001",
                "sample_id": "Child",
                "father_id": "Father",
                "mother_id": "Mother",
                "sex": "1",
                "affected_status": "2",
            },
        }
        ctx.config["inheritance_mode"] = "simple"
        ctx.mark_complete("dataframe_loading")
        ctx.mark_complete("sample_config_loading")
        ctx.mark_complete("pedigree_loading")
        return ctx

    def test_sample_columns_created_from_gt(self, context):
        """Test that sample columns are created from GT field."""
        stage = InheritanceAnalysisStage()

        # Mock analyze_inheritance to check inputs
        with patch("variantcentrifuge.stages.analysis_stages.analyze_inheritance") as mock_analyze:
            # Return DataFrame with inheritance columns added
            mock_df = context.current_dataframe.copy()
            mock_df["Inheritance_Pattern"] = "unknown"
            mock_df["Inheritance_Details"] = "{}"
            mock_analyze.return_value = mock_df

            # Mock process_inheritance_output
            with patch(
                "variantcentrifuge.inheritance.analyzer.process_inheritance_output"
            ) as mock_process:
                mock_process.return_value = mock_df

                stage(context)

        # Verify analyze_inheritance was called
        mock_analyze.assert_called_once()
        call_args = mock_analyze.call_args

        # Check DataFrame has sample columns
        df_arg = call_args[1]["df"]
        assert "Child" in df_arg.columns
        assert "Father" in df_arg.columns
        assert "Mother" in df_arg.columns

        # Verify genotypes were extracted correctly
        assert df_arg.loc[0, "Child"] == "0/1"
        assert df_arg.loc[0, "Father"] == "0/0"
        assert df_arg.loc[0, "Mother"] == "0/0"

        # Compound het candidate in GENE2
        assert df_arg.loc[1, "Child"] == "0/1"
        assert df_arg.loc[1, "Father"] == "0/1"
        assert df_arg.loc[1, "Mother"] == "0/0"

        assert df_arg.loc[2, "Child"] == "0/1"
        assert df_arg.loc[2, "Father"] == "0/0"
        assert df_arg.loc[2, "Mother"] == "0/1"

    def test_replaced_genotype_format(self, context):
        """Test handling of replaced genotype format."""
        # Change GT format to replaced format
        context.current_dataframe["GT"] = [
            "Child(0/1:30);Father(0/0:25);Mother(0/0:28)",
            "Child(0/1:25);Father(0/1:30);Mother(0/0:27)",
            "Child(0/1:28);Father(0/0:26);Mother(0/1:29)",
            "Child(1/1:35);Father(0/1:30);Mother(0/1:32)",
        ]

        stage = InheritanceAnalysisStage()

        with patch("variantcentrifuge.stages.analysis_stages.analyze_inheritance") as mock_analyze:
            mock_df = context.current_dataframe.copy()
            mock_df["Inheritance_Pattern"] = "unknown"
            mock_df["Inheritance_Details"] = "{}"
            mock_analyze.return_value = mock_df

            with patch(
                "variantcentrifuge.inheritance.analyzer.process_inheritance_output"
            ) as mock_process:
                mock_process.return_value = mock_df

                stage(context)

        # Check DataFrame has sample columns with correct genotypes
        df_arg = mock_analyze.call_args[1]["df"]
        assert df_arg.loc[0, "Child"] == "0/1"
        assert df_arg.loc[0, "Father"] == "0/0"
        assert df_arg.loc[0, "Mother"] == "0/0"

    def test_inheritance_with_set_samples(self, context):
        """Test that vcf_samples is converted from set to list."""
        # Make vcf_samples a set
        context.vcf_samples = {"Child", "Father", "Mother"}

        stage = InheritanceAnalysisStage()

        with patch("variantcentrifuge.stages.analysis_stages.analyze_inheritance") as mock_analyze:
            mock_df = context.current_dataframe.copy()
            mock_df["Inheritance_Pattern"] = "unknown"
            mock_df["Inheritance_Details"] = "{}"
            mock_analyze.return_value = mock_df

            with patch(
                "variantcentrifuge.inheritance.analyzer.process_inheritance_output"
            ) as mock_process:
                mock_process.return_value = mock_df

                stage(context)

        # Verify sample_list was passed as list
        call_args = mock_analyze.call_args
        sample_list = call_args[1]["sample_list"]
        assert isinstance(sample_list, list)
        assert len(sample_list) == 3
        assert set(sample_list) == {"Child", "Father", "Mother"}

    def test_no_gt_column(self, context):
        """Test behavior when GT column is missing."""
        # Remove GT column
        context.current_dataframe = context.current_dataframe.drop(columns=["GT"])

        stage = InheritanceAnalysisStage()

        with patch("variantcentrifuge.stages.analysis_stages.analyze_inheritance") as mock_analyze:
            mock_df = context.current_dataframe.copy()
            mock_df["Inheritance_Pattern"] = "unknown"
            mock_df["Inheritance_Details"] = "{}"
            mock_analyze.return_value = mock_df

            with patch(
                "variantcentrifuge.inheritance.analyzer.process_inheritance_output"
            ) as mock_process:
                mock_process.return_value = mock_df

                stage(context)

        # Should still call analyze_inheritance but without sample columns
        mock_analyze.assert_called_once()
        df_arg = mock_analyze.call_args[1]["df"]
        assert "Child" not in df_arg.columns
        assert "Father" not in df_arg.columns
        assert "Mother" not in df_arg.columns

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = InheritanceAnalysisStage()
        assert stage.dependencies == {"dataframe_loading"}
        assert stage.soft_dependencies == {"custom_annotation"}

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = InheritanceAnalysisStage()
        assert stage.parallel_safe is False  # Must run in main process to modify DataFrame

    def test_empty_dataframe(self, context):
        """Test with empty DataFrame."""
        context.current_dataframe = pd.DataFrame()

        stage = InheritanceAnalysisStage()

        with patch("variantcentrifuge.stages.analysis_stages.analyze_inheritance") as mock_analyze:
            mock_df = pd.DataFrame()
            mock_df["Inheritance_Pattern"] = []
            mock_df["Inheritance_Details"] = []
            mock_analyze.return_value = mock_df

            with patch(
                "variantcentrifuge.inheritance.analyzer.process_inheritance_output"
            ) as mock_process:
                mock_process.return_value = mock_df

                result = stage(context)

        # Should handle empty DataFrame gracefully
        assert result is context

    def test_missing_gene_column(self, context):
        """Test inheritance analysis without GENE column."""
        # Remove GENE column
        context.current_dataframe = context.current_dataframe.drop(columns=["GENE"])

        stage = InheritanceAnalysisStage()

        with patch("variantcentrifuge.stages.analysis_stages.analyze_inheritance") as mock_analyze:
            mock_df = context.current_dataframe.copy()
            mock_df["Inheritance_Pattern"] = "unknown"
            mock_df["Inheritance_Details"] = "{}"
            mock_analyze.return_value = mock_df

            with patch(
                "variantcentrifuge.inheritance.analyzer.process_inheritance_output"
            ) as mock_process:
                mock_process.return_value = mock_df

                stage(context)

        # Should still work without GENE column
        mock_analyze.assert_called_once()
