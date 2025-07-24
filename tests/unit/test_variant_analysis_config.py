"""Unit tests for VariantAnalysisStage case/control configuration passing."""

# import tempfile  # Unused import
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from tests.mocks.fixtures import create_test_context
from variantcentrifuge.stages.analysis_stages import VariantAnalysisStage


class TestVariantAnalysisConfig:
    """Test VariantAnalysisStage case/control configuration passing."""

    @pytest.fixture
    def context(self):
        """Create test context with case/control assignments."""
        context = create_test_context(
            gene_name="BRCA1",
            vcf_file="/tmp/test.vcf",
            config_overrides={
                "case_samples": ["S001", "S002"],
                "control_samples": ["S003", "S004"],
                "case_phenotypes": ["HP:0000001"],
                "control_phenotypes": [],
            },
        )

        # Mock VCF samples
        context.vcf_samples = ["S001", "S002", "S003", "S004"]

        # Create a mock DataFrame with test data
        test_df = pd.DataFrame(
            {
                "CHROM": ["1", "1", "2"],
                "POS": [1000, 2000, 3000],
                "REF": ["A", "G", "C"],
                "ALT": ["T", "C", "A"],
                "GT": ["S001(0/1)", "S002(1/1)", "S003(0/1)"],
            }
        )
        context.current_dataframe = test_df

        # Mock workspace for temporary files
        mock_workspace = MagicMock()
        mock_workspace.get_intermediate_path.return_value = Path("/tmp/temp_analysis.tsv")
        mock_workspace.base_name = "test_analysis"
        context.workspace = mock_workspace

        return context

    @patch("variantcentrifuge.stages.analysis_stages.analyze_variants")
    def test_case_control_config_passed_to_analyze_variants(self, mock_analyze_variants, context):
        """Test that case/control samples are passed to analyze_variants function."""
        # Mock analyze_variants to return test results
        mock_analyze_variants.return_value = [
            "CHROM\tPOS\tREF\tALT\tproband_count\tcontrol_count",
            "1\t1000\tA\tT\t2\t2",
            "1\t2000\tG\tC\t2\t2",
        ]

        # Mark dependencies as complete
        context.mark_complete("dataframe_loading")

        stage = VariantAnalysisStage()
        _ = stage(context)  # Execute stage, result not needed

        # Verify analyze_variants was called
        assert mock_analyze_variants.call_count == 1

        # Get the analysis_config that was passed
        call_args = mock_analyze_variants.call_args
        analysis_config = call_args[0][1]  # Second argument is the config

        # Verify case/control samples were passed
        assert "case_samples" in analysis_config
        assert "control_samples" in analysis_config
        assert analysis_config["case_samples"] == ["S001", "S002"]
        assert analysis_config["control_samples"] == ["S003", "S004"]

        # Verify phenotype terms were also passed for legacy compatibility
        assert "case_phenotypes" in analysis_config
        assert "control_phenotypes" in analysis_config
        assert analysis_config["case_phenotypes"] == ["HP:0000001"]
        assert analysis_config["control_phenotypes"] == []

    @patch("variantcentrifuge.stages.analysis_stages.analyze_variants")
    def test_empty_case_control_config_passed(self, mock_analyze_variants, context):
        """Test that empty case/control samples are passed correctly."""
        # Set empty case/control samples
        context.config["case_samples"] = []
        context.config["control_samples"] = []

        # Mock analyze_variants to return test results
        mock_analyze_variants.return_value = [
            "CHROM\tPOS\tREF\tALT\tproband_count\tcontrol_count",
            "1\t1000\tA\tT\t0\t4",
        ]

        # Mark dependencies as complete
        context.mark_complete("dataframe_loading")

        stage = VariantAnalysisStage()
        _ = stage(context)  # Execute stage, result not needed

        # Get the analysis_config that was passed
        call_args = mock_analyze_variants.call_args
        analysis_config = call_args[0][1]

        # Verify empty lists were passed
        assert analysis_config["case_samples"] == []
        assert analysis_config["control_samples"] == []

    @patch("variantcentrifuge.stages.analysis_stages.analyze_variants")
    def test_missing_case_control_config_defaults(self, mock_analyze_variants, context):
        """Test that missing case/control samples default to empty lists."""
        # Remove case/control samples from config
        if "case_samples" in context.config:
            del context.config["case_samples"]
        if "control_samples" in context.config:
            del context.config["control_samples"]

        # Mock analyze_variants to return test results
        mock_analyze_variants.return_value = [
            "CHROM\tPOS\tREF\tALT",
            "1\t1000\tA\tT",
        ]

        # Mark dependencies as complete
        context.mark_complete("dataframe_loading")

        stage = VariantAnalysisStage()
        _ = stage(context)  # Execute stage, result not needed

        # Get the analysis_config that was passed
        call_args = mock_analyze_variants.call_args
        analysis_config = call_args[0][1]

        # Verify defaults to empty lists
        assert analysis_config["case_samples"] == []
        assert analysis_config["control_samples"] == []

    @patch("variantcentrifuge.stages.analysis_stages.analyze_variants")
    def test_vcf_samples_config_passed(self, mock_analyze_variants, context):
        """Test that VCF samples and other config are passed correctly."""
        # Mock analyze_variants to return test results
        mock_analyze_variants.return_value = ["CHROM\tPOS", "1\t1000"]

        # Mark dependencies as complete
        context.mark_complete("dataframe_loading")

        stage = VariantAnalysisStage()
        _ = stage(context)  # Execute stage, result not needed

        # Get the analysis_config that was passed
        call_args = mock_analyze_variants.call_args
        analysis_config = call_args[0][1]

        # Verify VCF samples and other required config
        assert analysis_config["vcf_sample_names"] == ["S001", "S002", "S003", "S004"]
        assert analysis_config["sample_list"] == "S001,S002,S003,S004"
        assert analysis_config["base_name"] == "test_analysis"
        assert "reference" in analysis_config
