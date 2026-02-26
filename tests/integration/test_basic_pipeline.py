"""Basic integration test for the stage-based pipeline."""

import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

pytestmark = pytest.mark.integration

from tests.mocks import create_test_context, create_test_phenotype_file, create_test_vcf
from variantcentrifuge.pipeline_core import PipelineRunner
from variantcentrifuge.stages import (
    ConfigurationLoadingStage,
    PhenotypeLoadingStage,
    SampleConfigLoadingStage,
)


class TestBasicPipeline:
    """Integration tests for basic pipeline functionality."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def test_data(self, temp_dir):
        """Create test data files."""
        # Create test VCF
        vcf_path = create_test_vcf(temp_dir / "test.vcf", samples=["Sample1", "Sample2", "Sample3"])

        # Create phenotype file
        phenotype_path = create_test_phenotype_file(
            temp_dir / "phenotypes.tsv",
            phenotypes={
                "Sample1": "Case",
                "Sample2": "Control",
                "Sample3": "Case",
            },
        )

        return {
            "vcf": vcf_path,
            "phenotypes": phenotype_path,
            "output_dir": temp_dir / "output",
        }

    def test_configuration_loading(self, test_data):
        """Test configuration loading stage."""
        # Create context
        context = create_test_context(
            vcf_file=str(test_data["vcf"]),
            output_dir=str(test_data["output_dir"]),
            gene_name="BRCA1",
            config_overrides={
                "test_option": "test_value",
            },
        )

        # Run configuration stage
        stage = ConfigurationLoadingStage()
        result = stage(context)

        # Verify configuration loaded
        assert result.config["vcf_file"] == str(test_data["vcf"])
        assert result.config["gene_name"] == "BRCA1"
        assert result.config["test_option"] == "test_value"
        assert result.is_complete("configuration_loading")

    def test_phenotype_loading(self, test_data):
        """Test phenotype loading stage."""
        # Create context with phenotype file
        context = create_test_context(
            vcf_file=str(test_data["vcf"]),
            output_dir=str(test_data["output_dir"]),
            phenotype_file=str(test_data["phenotypes"]),
            phenotype_sample_column="Sample",
            phenotype_value_column="Phenotype",
        )

        # Load configuration first
        ConfigurationLoadingStage()(context)

        # Run phenotype loading
        stage = PhenotypeLoadingStage()
        result = stage(context)

        # Verify phenotypes loaded
        assert result.phenotype_data is not None
        assert len(result.phenotype_data) == 3
        assert result.phenotype_data["Sample1"] == {"Case"}
        assert result.phenotype_data["Sample2"] == {"Control"}
        assert result.phenotype_data["Sample3"] == {"Case"}
        assert result.is_complete("phenotype_loading")

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_sample_loading(self, mock_get_samples, test_data):
        """Test sample configuration loading."""
        # Mock VCF sample extraction
        mock_get_samples.return_value = ["Sample1", "Sample2", "Sample3"]

        # Create context
        context = create_test_context(
            vcf_file=str(test_data["vcf"]),
            output_dir=str(test_data["output_dir"]),
        )

        # Load configuration first
        ConfigurationLoadingStage()(context)

        # Run sample loading
        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Verify samples loaded
        assert result.vcf_samples == ["Sample1", "Sample2", "Sample3"]
        assert result.is_complete("sample_config_loading")
        mock_get_samples.assert_called_once_with(str(test_data["vcf"]))

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_parallel_configuration_stages(self, mock_get_samples, test_data):
        """Test parallel execution of configuration stages."""
        # Mock VCF sample extraction
        mock_get_samples.return_value = ["Sample1", "Sample2", "Sample3"]

        # Create context
        context = create_test_context(
            vcf_file=str(test_data["vcf"]),
            output_dir=str(test_data["output_dir"]),
            phenotype_file=str(test_data["phenotypes"]),
            phenotype_sample_column="Sample",
            phenotype_value_column="Phenotype",
        )

        # For parallel execution test, ensure config has phenotype settings
        # since stages may run before configuration_loading completes
        context.config["phenotype_file"] = str(test_data["phenotypes"])
        context.config["phenotype_sample_column"] = "Sample"
        context.config["phenotype_value_column"] = "Phenotype"

        # Create stages
        stages = [
            ConfigurationLoadingStage(),
            PhenotypeLoadingStage(),
            SampleConfigLoadingStage(),
        ]

        # Run pipeline
        runner = PipelineRunner(max_workers=3)
        result = runner.run(stages, context)

        # Verify all stages completed
        assert result.is_complete("configuration_loading")
        assert result.is_complete("phenotype_loading")
        assert result.is_complete("sample_config_loading")

        # Verify data loaded
        assert result.config is not None
        assert result.phenotype_data is not None
        assert result.vcf_samples is not None

    def test_execution_plan(self):
        """Test execution plan generation."""
        # Create stages with dependencies
        stages = [
            ConfigurationLoadingStage(),  # No dependencies
            PhenotypeLoadingStage(),  # No dependencies
            SampleConfigLoadingStage(),  # Depends on configuration_loading
        ]

        runner = PipelineRunner()
        plan = runner.dry_run(stages)

        # Configuration and phenotype can run in parallel
        # Sample loading must wait for configuration
        assert len(plan) == 2
        assert set(plan[0]) == {"configuration_loading", "phenotype_loading"}
        assert plan[1] == ["sample_config_loading"]

    def test_workspace_cleanup(self, test_data):
        """Test workspace cleanup after pipeline."""
        context = create_test_context(
            vcf_file=str(test_data["vcf"]),
            output_dir=str(test_data["output_dir"]),
        )

        # Create some temporary files
        temp_file = context.workspace.get_temp_path("test.txt")
        temp_file.write_text("test data")
        assert temp_file.exists()

        # Use context manager for automatic cleanup
        with context.workspace:
            assert temp_file.exists()

        # Temp files should be cleaned up
        assert not temp_file.exists()

    def test_stage_failure_handling(self):
        """Test handling of stage failures."""
        from variantcentrifuge.pipeline_core import Stage

        class FailingStage(Stage):
            @property
            def name(self):
                return "failing_stage"

            @property
            def dependencies(self):
                return {"configuration_loading"}

            def _process(self, context):
                raise ValueError("Intentional failure")

        context = create_test_context()
        stages = [
            ConfigurationLoadingStage(),
            FailingStage(),
        ]

        runner = PipelineRunner()

        with pytest.raises(ValueError) as exc_info:
            runner.run(stages, context)

        assert "Intentional failure" in str(exc_info.value)
        assert context.is_complete("configuration_loading")
        assert not context.is_complete("failing_stage")
