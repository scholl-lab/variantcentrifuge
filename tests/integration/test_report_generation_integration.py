"""Simplified integration tests for IGV and HTML report generation."""

import json
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest
import pandas as pd

from variantcentrifuge.stages.output_stages import IGVReportStage, HTMLReportStage
from tests.mocks.fixtures import create_test_context


class TestReportGenerationIntegration:
    """Simplified integration tests for report generation functionality."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def sample_dataframe(self):
        """Create a sample DataFrame for testing."""
        return pd.DataFrame({
            "CHROM": ["17", "17"],
            "POS": [41223094, 7579472],
            "REF": ["T", "G"],
            "ALT": ["C", "C"],
            "GT": ["Sample1(0/1);Sample2(0/0)", "Sample1(1/1);Sample2(0/1)"],
            "GENE": ["BRCA1", "TP53"],
        })

    @pytest.fixture
    def simple_bam_mapping(self, temp_output_dir):
        """Create a simple BAM mapping file."""
        bam_file = temp_output_dir / "bam_mapping.txt"
        bam_content = """Sample1,/path/to/sample1.bam
Sample2,/path/to/sample2.bam
"""
        bam_file.write_text(bam_content)
        return bam_file

    def test_igv_report_stage_integration(self, temp_output_dir, sample_dataframe, simple_bam_mapping):
        """Test IGV report stage with mocked dependencies."""
        # Create test context
        context = create_test_context(
            output_dir=str(temp_output_dir),
            config_overrides={
                "igv_enabled": True,
                "igv_reference": "hg19",
                "bam_mapping_file": str(simple_bam_mapping),
            }
        )
        context.current_dataframe = sample_dataframe
        context.mark_complete("tsv_output")
        
        # Create a fake TSV output file
        tsv_file = temp_output_dir / "test_output.tsv"
        sample_dataframe.to_csv(tsv_file, sep='\t', index=False)
        context.final_output_path = tsv_file

        # Mock the IGV report generation function
        with patch("variantcentrifuge.stages.output_stages.generate_igv_report") as mock_igv:
            # Execute stage
            stage = IGVReportStage()
            result = stage(context)

            # Verify IGV report generation was called
            mock_igv.assert_called_once()
            call_kwargs = mock_igv.call_args[1]
            assert call_kwargs["variants_tsv"] == str(tsv_file)
            assert call_kwargs["bam_mapping_file"] == str(simple_bam_mapping)
            assert call_kwargs["igv_reference"] == "hg19"

    def test_html_report_stage_integration(self, temp_output_dir, sample_dataframe):
        """Test HTML report stage with mocked dependencies."""
        # Create test context
        context = create_test_context(
            output_dir=str(temp_output_dir),
            config_overrides={"html_report": True}
        )
        context.current_dataframe = sample_dataframe
        context.mark_complete("tsv_output")
        
        # Create a fake TSV output file
        tsv_file = temp_output_dir / "test_output.tsv"
        sample_dataframe.to_csv(tsv_file, sep='\t', index=False)
        context.final_output_path = tsv_file

        # Mock the report generation functions
        with patch("variantcentrifuge.stages.output_stages.produce_report_json") as mock_json, \
             patch("variantcentrifuge.stages.output_stages.generate_html_report") as mock_html:
            
            # Mock produce_report_json to create the expected JSON files
            def create_json_files(*args, **kwargs):
                report_dir = temp_output_dir / "report"
                report_dir.mkdir(parents=True, exist_ok=True)
                
                # Create variants.json
                variants_json = report_dir / "variants.json"
                variants_json.write_text('[{"chrom": "17", "pos": "41223094"}]')
                
                # Create summary.json
                summary_json = report_dir / "summary.json"
                summary_json.write_text('{"num_variants": 2}')

            mock_json.side_effect = create_json_files
            
            # Mock generate_html_report to create index.html
            def create_html_file(*args, **kwargs):
                report_dir = temp_output_dir / "report"
                index_html = report_dir / "index.html"
                index_html.write_text("<html><body>Test Report</body></html>")
            
            mock_html.side_effect = create_html_file

            # Execute stage
            stage = HTMLReportStage()
            result = stage(context)

            # Verify report generation was called
            mock_json.assert_called_once()
            mock_html.assert_called_once()
            
            # Verify report paths were set
            assert "html" in result.report_paths

    def test_combined_reports_integration(self, temp_output_dir, sample_dataframe, simple_bam_mapping):
        """Test both IGV and HTML report generation together."""
        # Create test context for both reports
        context = create_test_context(
            output_dir=str(temp_output_dir),
            config_overrides={
                "igv_enabled": True,
                "igv_reference": "hg19",
                "bam_mapping_file": str(simple_bam_mapping),
                "html_report": True,
            }
        )
        context.current_dataframe = sample_dataframe
        context.mark_complete("tsv_output")
        
        # Create a fake TSV output file
        tsv_file = temp_output_dir / "test_output.tsv"
        sample_dataframe.to_csv(tsv_file, sep='\t', index=False)
        context.final_output_path = tsv_file

        # Mock both report generation functions
        with patch("variantcentrifuge.stages.output_stages.generate_igv_report") as mock_igv, \
             patch("variantcentrifuge.stages.output_stages.produce_report_json") as mock_json, \
             patch("variantcentrifuge.stages.output_stages.generate_html_report") as mock_html:
            
            # Mock JSON creation for HTML report
            def create_json_files(*args, **kwargs):
                report_dir = temp_output_dir / "report"
                report_dir.mkdir(parents=True, exist_ok=True)
                
                variants_json = report_dir / "variants.json"
                variants_json.write_text('[{"chrom": "17", "pos": "41223094"}]')
                
                summary_json = report_dir / "summary.json"
                summary_json.write_text('{"num_variants": 2}')

            mock_json.side_effect = create_json_files
            
            # Mock generate_html_report to create index.html
            def create_html_file(*args, **kwargs):
                report_dir = temp_output_dir / "report"
                index_html = report_dir / "index.html"
                index_html.write_text("<html><body>Test Report</body></html>")
            
            mock_html.side_effect = create_html_file

            # Execute both stages
            igv_stage = IGVReportStage()
            html_stage = HTMLReportStage()
            
            igv_result = igv_stage(context)
            html_result = html_stage(igv_result)

            # Verify both report types were called
            mock_igv.assert_called_once()
            mock_json.assert_called_once()
            mock_html.assert_called_once()

    def test_report_error_handling(self, temp_output_dir, sample_dataframe):
        """Test error handling in report generation."""
        # Create test context
        context = create_test_context(
            output_dir=str(temp_output_dir),
            config_overrides={"html_report": True}
        )
        context.current_dataframe = sample_dataframe
        context.mark_complete("tsv_output")
        
        # Create a fake TSV output file
        tsv_file = temp_output_dir / "test_output.tsv"
        sample_dataframe.to_csv(tsv_file, sep='\t', index=False)
        context.final_output_path = tsv_file

        # Mock report generation to raise an exception
        with patch("variantcentrifuge.stages.output_stages.produce_report_json") as mock_json:
            mock_json.side_effect = Exception("Simulated failure")

            # Execute stage - should raise the exception
            stage = HTMLReportStage()
            with pytest.raises(Exception, match="Simulated failure"):
                stage(context)

    def test_missing_dependencies_handling(self, temp_output_dir):
        """Test handling when required dependencies are missing."""
        # Create test context without TSV output
        context = create_test_context(
            output_dir=str(temp_output_dir),
            config_overrides={"igv_enabled": True, "bam_mapping_file": "/nonexistent/path"}
        )
        # Don't mark tsv_output as complete
        context.final_output_path = None

        # Execute stage with missing dependencies - should raise RuntimeError
        stage = IGVReportStage()
        with pytest.raises(RuntimeError, match="Stage 'igv_report' requires these stages to complete first"):
            stage(context)