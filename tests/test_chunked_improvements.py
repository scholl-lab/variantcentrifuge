"""Tests for the chunked processing improvements.

Tests for the chunked processing improvements including parallel processing,
error handling, timing tracking, and unified cleanup.
"""

import argparse
import logging
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import (
    ChunkedAnalysisStage,
    cleanup_sample_columns,
    handle_inheritance_analysis_error,
)


class TestChunkedProcessingImprovements:
    """Test suite for chunked processing improvements."""

    @pytest.fixture
    def mock_args(self) -> argparse.Namespace:
        """Create mock command-line arguments."""
        return argparse.Namespace(
            vcf_file="test.vcf",
            output_dir="output",
            output_file="test_output.tsv",
            log_level="INFO",
            chunks=100,
            threads=4,  # Enable parallel processing
        )

    @pytest.fixture
    def mock_workspace(self, tmp_path) -> Workspace:
        """Create mock workspace for testing."""
        workspace = Mock(spec=Workspace)
        workspace.output_dir = tmp_path / "output"
        workspace.intermediate_dir = tmp_path / "intermediate"
        workspace.intermediate_dir.mkdir(parents=True, exist_ok=True)
        workspace.output_dir.mkdir(parents=True, exist_ok=True)
        return workspace

    @pytest.fixture
    def sample_config(self) -> dict:
        """Sample configuration for testing."""
        return {
            "separator": ";",
            "extract_fields_separator": ":",
            "genotype_replacement_map": {r"[2-9]": "1"},
            "sample_list": "Sample1,Sample2,Sample3",
            "use_chunked_processing": True,
            "chunk_size": 100,
            "inheritance_mode": "simple",
            "vcf_samples": ["Sample1", "Sample2", "Sample3"],
            "threads": 4,
            "min_variants_for_parallel_inheritance": 10,
        }

    @pytest.fixture
    def sample_dataframe(self) -> pd.DataFrame:
        """Sample DataFrame with genetic variant data."""
        # Create a larger DataFrame to trigger parallel processing
        base_data = {
            "CHROM": ["chr1"] * 20,
            "POS": list(range(100, 120)),
            "ID": [f"rs{i}" for i in range(1, 21)],
            "REF": ["A"] * 20,
            "ALT": ["T"] * 20,
            "GENE": ["BRCA1"] * 20,
            "GT": ["Sample1(0/1);Sample2(0/0);Sample3(1/1)"] * 20,
            "IMPACT": ["HIGH"] * 20,
        }
        return pd.DataFrame(base_data)

    @pytest.fixture
    def sample_pedigree(self) -> dict:
        """Sample pedigree data for inheritance analysis."""
        return {
            "Sample1": {
                "family_id": "FAM001",
                "individual_id": "Sample1",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "2",
            },
            "Sample2": {
                "family_id": "FAM001",
                "individual_id": "Sample2",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
            "Sample3": {
                "family_id": "FAM001",
                "individual_id": "Sample3",
                "father_id": "Sample1",
                "mother_id": "Sample2",
                "sex": "1",
                "affected_status": "2",
            },
        }

    def test_parallel_processing_support(
        self, sample_config, sample_dataframe, sample_pedigree, mock_args, mock_workspace
    ):
        """Test that chunked processing supports parallel inheritance analysis."""
        # Set up context for parallel processing
        context = PipelineContext(mock_args, sample_config, mock_workspace)
        context.vcf_samples = sample_config["vcf_samples"]
        context.pedigree_data = sample_pedigree

        # Set up inheritance config
        context.config["inheritance_analysis_config"] = {
            "inheritance_mode": "simple",
            "vcf_samples": sample_config["vcf_samples"],
            "pedigree_data": sample_pedigree,
            "use_vectorized_comp_het": True,
            "preserve_details_for_scoring": False,
        }

        # Mock the parallel analyzer to verify it's called
        with (
            patch(
                "variantcentrifuge.inheritance.parallel_analyzer.analyze_inheritance_parallel"
            ) as mock_parallel,
            patch("variantcentrifuge.inheritance.analyzer.analyze_inheritance") as mock_sequential,
        ):
            # Make the parallel analyzer return the DataFrame with inheritance patterns
            mock_parallel.return_value = sample_dataframe.copy()
            mock_parallel.return_value["Inheritance_Pattern"] = "autosomal_dominant"

            chunked_stage = ChunkedAnalysisStage()
            chunked_stage._apply_chunk_inheritance_analysis(sample_dataframe, context)

            # Verify parallel analyzer was called (not sequential)
            mock_parallel.assert_called_once()
            mock_sequential.assert_not_called()

            # Verify the parallel analyzer was called with correct parameters
            call_args = mock_parallel.call_args
            assert call_args[1]["n_workers"] == 4
            assert call_args[1]["min_variants_for_parallel"] == 10
            assert call_args[1]["sample_list"] == sample_config["vcf_samples"]
            assert call_args[1]["pedigree_data"] == sample_pedigree
            assert call_args[1]["use_vectorized_comp_het"] is True

    def test_error_handling_consistency(
        self, sample_config, sample_dataframe, mock_args, mock_workspace
    ):
        """Test that error handling works consistently for chunked processing."""
        context = PipelineContext(mock_args, sample_config, mock_workspace)
        context.vcf_samples = sample_config["vcf_samples"]

        # Set up inheritance config
        context.config["inheritance_analysis_config"] = {
            "inheritance_mode": "simple",
            "vcf_samples": sample_config["vcf_samples"],
            "pedigree_data": {},
            "use_vectorized_comp_het": True,
            "preserve_details_for_scoring": True,  # Test details preservation
        }

        # Mock the sample column creation to raise an exception
        with patch(
            "variantcentrifuge.stages.analysis_stages.create_sample_columns_from_gt"
        ) as mock_create_columns:
            mock_create_columns.side_effect = Exception("Test error")

            chunked_stage = ChunkedAnalysisStage()
            result = chunked_stage._apply_chunk_inheritance_analysis(sample_dataframe, context)

            # Verify error recovery
            assert isinstance(result, pd.DataFrame)
            assert len(result) == len(sample_dataframe)
            assert "Inheritance_Pattern" in result.columns
            assert "Inheritance_Details" in result.columns
            assert all(result["Inheritance_Pattern"] == "error")
            assert all(result["Inheritance_Details"] == "{}")

    def test_timing_tracking_functionality(
        self, sample_config, sample_dataframe, sample_pedigree, mock_args, mock_workspace
    ):
        """Test that timing tracking works correctly in chunked processing."""
        context = PipelineContext(mock_args, sample_config, mock_workspace)
        context.vcf_samples = sample_config["vcf_samples"]
        context.pedigree_data = sample_pedigree

        # Set up inheritance config
        context.config["inheritance_analysis_config"] = {
            "inheritance_mode": "simple",
            "vcf_samples": sample_config["vcf_samples"],
            "pedigree_data": sample_pedigree,
            "use_vectorized_comp_het": True,
            "preserve_details_for_scoring": False,
        }

        # Just run the function and verify it executes without error
        # (timing is logged internally but we can verify the function ran)
        chunked_stage = ChunkedAnalysisStage()
        result = chunked_stage._apply_chunk_inheritance_analysis(sample_dataframe, context)

        # Verify the result contains inheritance patterns
        assert "Inheritance_Pattern" in result.columns
        assert isinstance(result, pd.DataFrame)
        assert len(result) == len(sample_dataframe)

    def test_unified_cleanup_function(self):
        """Test the unified cleanup_sample_columns function."""
        # Create test DataFrame with sample columns
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2"],
                "POS": [100, 200],
                "GT": ["0/1", "1/1"],
                "Sample1": ["0/1", "1/1"],
                "Sample2": ["0/0", "0/1"],
                "Sample3": ["1/1", "0/0"],
                "Inheritance_Pattern": ["autosomal_dominant", "autosomal_recessive"],
                "Inheritance_Details": ["{}", "{}"],
            }
        )

        vcf_samples = ["Sample1", "Sample2", "Sample3"]

        # Test default behavior (preserve GT, Inheritance_Pattern, Inheritance_Details)
        result_df = cleanup_sample_columns(df, vcf_samples)

        # Should remove Sample1, Sample2, Sample3 but keep GT, Inheritance_Pattern,
        # Inheritance_Details
        assert "Sample1" not in result_df.columns
        assert "Sample2" not in result_df.columns
        assert "Sample3" not in result_df.columns
        assert "GT" in result_df.columns
        assert "Inheritance_Pattern" in result_df.columns
        assert "Inheritance_Details" in result_df.columns
        assert "CHROM" in result_df.columns
        assert "POS" in result_df.columns

        # Test with empty preserve_columns (remove all sample columns)
        result_empty = cleanup_sample_columns(df, vcf_samples, preserve_columns=[])

        # Should remove all sample columns but keep GT (GT is not a VCF sample name)
        assert "Sample1" not in result_empty.columns
        assert "Sample2" not in result_empty.columns
        assert "Sample3" not in result_empty.columns
        assert "GT" in result_empty.columns  # GT is not a VCF sample name, so it's preserved
        assert "Inheritance_Pattern" in result_empty.columns
        assert "Inheritance_Details" in result_empty.columns

        # Test with custom preserve_columns
        result_custom = cleanup_sample_columns(df, vcf_samples, preserve_columns=["Sample1", "GT"])

        # Should only preserve Sample1 and GT
        assert "Sample1" in result_custom.columns
        assert "Sample2" not in result_custom.columns
        assert "Sample3" not in result_custom.columns
        assert "GT" in result_custom.columns

    def test_error_handling_function(self):
        """Test the unified handle_inheritance_analysis_error function."""
        # Create test DataFrame
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2"],
                "POS": [100, 200],
                "GT": ["0/1", "1/1"],
            }
        )

        test_error = ValueError("Test error message")

        # Test without details preservation
        result = handle_inheritance_analysis_error(
            df, test_error, preserve_details_for_scoring=False
        )

        assert isinstance(result, pd.DataFrame)
        assert len(result) == len(df)
        assert "Inheritance_Pattern" in result.columns
        assert "Inheritance_Details" not in result.columns
        assert all(result["Inheritance_Pattern"] == "error")

        # Test with details preservation
        result_with_details = handle_inheritance_analysis_error(
            df, test_error, preserve_details_for_scoring=True
        )

        assert "Inheritance_Pattern" in result_with_details.columns
        assert "Inheritance_Details" in result_with_details.columns
        assert all(result_with_details["Inheritance_Pattern"] == "error")
        assert all(result_with_details["Inheritance_Details"] == "{}")

    def test_sequential_fallback(
        self, sample_config, sample_dataframe, sample_pedigree, mock_args, mock_workspace
    ):
        """Test that chunked processing falls back to sequential for small datasets."""
        # Configure for sequential processing (low thread count or small dataset)
        sample_config["threads"] = 1

        context = PipelineContext(mock_args, sample_config, mock_workspace)
        context.vcf_samples = sample_config["vcf_samples"]
        context.pedigree_data = sample_pedigree

        # Set up inheritance config
        context.config["inheritance_analysis_config"] = {
            "inheritance_mode": "simple",
            "vcf_samples": sample_config["vcf_samples"],
            "pedigree_data": sample_pedigree,
            "use_vectorized_comp_het": True,
            "preserve_details_for_scoring": False,
        }

        # Mock the analyzers to verify which one is called
        with (
            patch(
                "variantcentrifuge.inheritance.parallel_analyzer.analyze_inheritance_parallel"
            ) as mock_parallel,
            patch("variantcentrifuge.inheritance.analyzer.analyze_inheritance") as mock_sequential,
        ):
            # Make the sequential analyzer return the DataFrame with inheritance patterns
            mock_sequential.return_value = sample_dataframe.copy()
            mock_sequential.return_value["Inheritance_Pattern"] = "autosomal_dominant"

            chunked_stage = ChunkedAnalysisStage()
            chunked_stage._apply_chunk_inheritance_analysis(sample_dataframe, context)

            # Verify sequential analyzer was called (not parallel)
            mock_sequential.assert_called_once()
            mock_parallel.assert_not_called()

            # Verify the sequential analyzer was called with correct parameters
            call_args = mock_sequential.call_args
            assert call_args[1]["sample_list"] == sample_config["vcf_samples"]
            assert call_args[1]["pedigree_data"] == sample_pedigree
            assert call_args[1]["use_vectorized_comp_het"] is True


if __name__ == "__main__":
    # Run tests with logging
    logging.basicConfig(level=logging.DEBUG)
    pytest.main([__file__, "-v"])
