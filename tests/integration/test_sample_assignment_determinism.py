"""
Integration tests for sample assignment determinism.

These tests verify that the entire pipeline produces consistent results
when run multiple times with the same input, preventing regression of
the critical sample assignment randomization bug.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch
from variantcentrifuge.pipeline_refactored import run_refactored_pipeline
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace


class TestPipelineSampleDeterminism:
    """Integration tests for end-to-end sample assignment determinism."""

    @pytest.fixture
    def temp_workspace(self):
        """Create a temporary workspace for testing."""
        temp_dir = Path(tempfile.mkdtemp())
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)

    @pytest.fixture
    def mock_test_data(self):
        """Mock test data for pipeline runs."""
        return {
            "vcf_samples": ["Sample3", "Sample1", "Sample2"],
            "test_tsv": "CHROM\tPOS\tREF\tALT\tGT\tGENE\nchr1\t100\tA\tT\t0/1:1/1:0/0\tBRCA1\n",
            "config": {
                "gene_name": "BRCA1",
                "vcf_file": "/tmp/test.vcf",
                "output_file": "test_output.tsv",
                "use_new_pipeline": True,
                "no_stats": True,
                "no_metadata": True,
                "extract": ["CHROM", "POS", "REF", "ALT", "GT", "GENE"],
                "extract_fields_separator": ":",
            },
        }

    @patch("variantcentrifuge.helpers.get_vcf_names")
    @patch("variantcentrifuge.stages.processing_stages.extract_fields")
    @patch("variantcentrifuge.stages.processing_stages.Path.exists")
    def test_multiple_pipeline_runs_identical_results(
        self, mock_exists, mock_extract, mock_get_names, temp_workspace, mock_test_data
    ):
        """Test that multiple pipeline runs produce identical results."""
        # Setup mocks
        mock_get_names.return_value = mock_test_data["vcf_samples"]
        mock_exists.return_value = True

        def mock_extract_fields(variant_file, fields, cfg, output_file):
            """Mock extract_fields to write test TSV data."""
            with open(output_file, "w") as f:
                f.write(mock_test_data["test_tsv"])

        mock_extract.side_effect = mock_extract_fields

        # Run pipeline multiple times
        results = []
        for run_num in range(3):
            run_workspace = temp_workspace / f"run_{run_num}"
            run_workspace.mkdir()

            config = mock_test_data["config"].copy()
            workspace = Workspace(str(run_workspace), "test")
            from argparse import Namespace

            args = Namespace()
            context = PipelineContext(args, config, workspace)

            with patch("pathlib.Path.exists", return_value=True):
                # Mock the pipeline run to focus on sample assignment
                from variantcentrifuge.stages.setup_stages import SampleConfigLoadingStage
                from variantcentrifuge.stages.processing_stages import (
                    FieldExtractionStage,
                    GenotypeReplacementStage,
                )

                # Load samples
                sample_stage = SampleConfigLoadingStage()
                context = sample_stage._process(context)

                # Extract fields (mocked)
                field_stage = FieldExtractionStage()
                context = field_stage._process(context)

                # Replace genotypes
                genotype_stage = GenotypeReplacementStage()
                context = genotype_stage._process(context)

                # Store sample order and genotype replacement result
                result_data = {
                    "vcf_samples": context.vcf_samples.copy(),
                    "genotype_file_exists": context.genotype_replaced_tsv
                    and context.genotype_replaced_tsv.exists(),
                }

                # Read genotype replacement output if it exists
                if result_data["genotype_file_exists"]:
                    with open(context.genotype_replaced_tsv, "r") as f:
                        result_data["genotype_content"] = f.read()

                results.append(result_data)

        # Verify all runs produced identical results
        first_result = results[0]
        for i, result in enumerate(results[1:], 1):
            assert (
                result["vcf_samples"] == first_result["vcf_samples"]
            ), f"Run {i+1} sample order differs from run 1"

            if first_result["genotype_file_exists"] and result["genotype_file_exists"]:
                assert (
                    result["genotype_content"] == first_result["genotype_content"]
                ), f"Run {i+1} genotype replacement output differs from run 1"

    def test_sample_order_propagation_through_stages(self, temp_workspace, mock_test_data):
        """Test that sample order is consistently propagated through all pipeline stages."""
        from variantcentrifuge.stages.setup_stages import SampleConfigLoadingStage
        from variantcentrifuge.stages.analysis_stages import InheritanceAnalysisStage

        test_samples = mock_test_data["vcf_samples"]

        with patch("variantcentrifuge.helpers.get_vcf_names", return_value=test_samples):
            workspace = Workspace(str(temp_workspace), "test")
            from argparse import Namespace

            args = Namespace()
            context = PipelineContext(args, mock_test_data["config"], workspace)

            # Track sample order through stages
            sample_orders = {}

            # 1. Sample loading stage
            sample_stage = SampleConfigLoadingStage()
            context = sample_stage._process(context)
            sample_orders["after_sample_loading"] = context.vcf_samples.copy()

            # 2. Inheritance analysis stage
            context.pedigree_data = {}  # Empty pedigree for this test
            context.current_dataframe = None
            inheritance_stage = InheritanceAnalysisStage()
            context = inheritance_stage._process(context)
            # Note: vcf_samples should remain unchanged by inheritance stage
            sample_orders["after_inheritance"] = context.vcf_samples.copy()

            # Verify sample order consistency
            original_order = test_samples
            for stage_name, stage_samples in sample_orders.items():
                assert (
                    stage_samples == original_order
                ), f"Sample order changed at {stage_name}: expected {original_order}, got {stage_samples}"

    @patch("variantcentrifuge.helpers.get_vcf_names")
    def test_genotype_replacement_determinism(self, mock_get_names, temp_workspace):
        """Test that genotype replacement produces deterministic output."""
        from variantcentrifuge.stages.processing_stages import GenotypeReplacementStage
        from variantcentrifuge.replacer import replace_genotypes

        test_samples = ["Sample1", "Sample2", "Sample3"]
        mock_get_names.return_value = test_samples

        # Create test TSV input
        test_input = [
            "CHROM\tPOS\tREF\tALT\tGT\tGENE",
            "chr1\t100\tA\tT\t0/1:1/1:0/0\tBRCA1",
            "chr1\t200\tC\tG\t1/1:0/1:0/0\tTP53",
        ]

        config = {
            "sample_list": ",".join(test_samples),
            "separator": ";",
            "extract_fields_separator": ":",
            "append_extra_sample_fields": False,
            "extra_sample_fields": [],
            "genotype_replacement_map": {},
        }

        # Run genotype replacement multiple times
        results = []
        for _ in range(3):
            output_lines = list(replace_genotypes(iter(test_input), config))
            results.append(output_lines)

        # All results should be identical
        first_result = results[0]
        for i, result in enumerate(results[1:], 1):
            assert result == first_result, f"Genotype replacement run {i+1} differs from run 1"

        # Verify sample order in output
        gt_line = first_result[2]  # Second data line
        assert "Sample1(" in gt_line, "Sample1 should appear in GT column"
        assert "Sample2(" in gt_line, "Sample2 should appear in GT column"
        # Sample3 has 0/0 genotype so should not appear (filtered out)

    def test_inheritance_analysis_determinism_with_real_data(self, temp_workspace):
        """Test inheritance analysis determinism with realistic family data."""
        from variantcentrifuge.stages.analysis_stages import InheritanceAnalysisStage
        import pandas as pd

        # Family trio data
        test_samples = ["Child", "Father", "Mother"]
        test_df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2"],
                "POS": [100, 200],
                "REF": ["A", "C"],
                "ALT": ["T", "G"],
                "GENE": ["BRCA1", "TP53"],
                "GT": [
                    "Child(0/1);Father(0/0);Mother(1/1)",  # De novo
                    "Child(1/1);Father(0/1);Mother(0/1)",  # Recessive
                ],
            }
        )

        pedigree_data = {
            "Child": {
                "sample_id": "Child",
                "father": "Father",
                "mother": "Mother",
                "affected": True,
            },
            "Father": {"sample_id": "Father", "affected": False},
            "Mother": {"sample_id": "Mother", "affected": False},
        }

        # Run inheritance analysis multiple times
        results = []
        for _ in range(3):
            workspace = Workspace(str(temp_workspace), "test")
            from argparse import Namespace

            args = Namespace()
            context = PipelineContext(
                args, {"inheritance_mode": "simple", "calculate_inheritance": True}, workspace
            )
            context.vcf_samples = test_samples.copy()
            context.current_dataframe = test_df.copy()
            context.pedigree_data = pedigree_data.copy()
            context.mark_complete("dataframe_loading")

            stage = InheritanceAnalysisStage()
            result_context = stage._process(context)

            if result_context.current_dataframe is not None:
                # Extract inheritance columns for comparison
                inheritance_cols = [
                    col for col in result_context.current_dataframe.columns if "Inheritance" in col
                ]
                if inheritance_cols:
                    inheritance_result = result_context.current_dataframe[inheritance_cols].to_dict(
                        "records"
                    )
                    results.append(inheritance_result)

        # All inheritance analysis results should be identical
        if results:
            first_result = results[0]
            for i, result in enumerate(results[1:], 1):
                assert result == first_result, f"Inheritance analysis run {i+1} differs from run 1"


@pytest.mark.integration
class TestPipelineReproducibility:
    """Integration tests for overall pipeline reproducibility."""

    def test_pipeline_output_hash_consistency(self):
        """Test that pipeline produces identical output hashes across runs."""
        # This test would run the full pipeline multiple times and compare
        # file hashes to ensure complete reproducibility
        pytest.skip("Requires full external tool setup - placeholder for future implementation")

    def test_sample_assignment_regression_prevention(self):
        """Test to prevent regression of sample assignment bug."""
        from variantcentrifuge.helpers import get_vcf_samples

        # Ensure get_vcf_samples never returns a set
        with patch("variantcentrifuge.helpers.get_vcf_names", return_value=["A", "B", "C"]):
            result = get_vcf_samples("dummy.vcf")

            assert isinstance(result, list), "get_vcf_samples must return list"
            assert not isinstance(result, set), "get_vcf_samples must not return set"

            # Multiple calls should return identical results
            results = [get_vcf_samples("dummy.vcf") for _ in range(10)]
            assert all(r == result for r in results), "Results must be deterministic"
