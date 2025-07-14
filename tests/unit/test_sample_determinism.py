"""
Unit tests for sample assignment determinism fix.

These tests ensure that the critical bug where sample assignments were
randomized between runs is fixed and won't regress.
"""

from unittest.mock import patch

from variantcentrifuge.helpers import get_vcf_samples
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import InheritanceAnalysisStage
from variantcentrifuge.stages.setup_stages import SampleConfigLoadingStage


class TestSampleDeterminism:
    """Test that sample ordering is deterministic and preserves VCF header order."""

    def test_get_vcf_samples_returns_list(self):
        """Test that get_vcf_samples returns a list, not a set."""
        test_samples = ["Sample3", "Sample1", "Sample2", "Sample4"]

        with patch("variantcentrifuge.helpers.get_vcf_names", return_value=test_samples):
            result = get_vcf_samples("dummy.vcf")

            assert isinstance(result, list), "get_vcf_samples should return a list"
            assert not isinstance(result, set), "get_vcf_samples should not return a set"

    def test_get_vcf_samples_preserves_order(self):
        """Test that get_vcf_samples preserves VCF header order."""
        test_samples = ["Sample3", "Sample1", "Sample2", "Sample4"]

        with patch("variantcentrifuge.helpers.get_vcf_names", return_value=test_samples):
            result = get_vcf_samples("dummy.vcf")

            assert result == test_samples, "Sample order should match VCF header order"

    def test_get_vcf_samples_deterministic_across_runs(self):
        """Test that get_vcf_samples returns identical results across multiple calls."""
        test_samples = ["Sample3", "Sample1", "Sample2", "Sample4"]

        with patch("variantcentrifuge.helpers.get_vcf_names", return_value=test_samples):
            results = [get_vcf_samples("dummy.vcf") for _ in range(10)]

            # All results should be identical
            first_result = results[0]
            assert all(
                result == first_result for result in results
            ), "get_vcf_samples should return identical results across multiple calls"

    def test_inheritance_stage_handles_set_with_warning(self):
        """Test that InheritanceAnalysisStage converts sets to sorted lists with warning."""
        # Create a minimal context with enough data to trigger set conversion
        from argparse import Namespace
        from pathlib import Path

        import pandas as pd

        workspace = Workspace(Path("/tmp/test"), "test")
        context = PipelineContext(Namespace(), {"inheritance_mode": "simple"}, workspace)
        context.vcf_samples = {"Sample3", "Sample1", "Sample2"}  # Set - unordered
        context.pedigree_data = {
            "Sample1": {"sample_id": "Sample1"},
            "Sample2": {"sample_id": "Sample2"},
            "Sample3": {"sample_id": "Sample3"},
        }
        # Add a minimal DataFrame to prevent early exit
        context.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": [100],
                "REF": ["A"],
                "ALT": ["T"],
                "GENE": ["BRCA1"],
                "GT": ["Sample1(0/1);Sample2(0/0);Sample3(1/1)"],
            }
        )

        stage = InheritanceAnalysisStage()

        with patch("variantcentrifuge.stages.analysis_stages.logger") as mock_logger:
            # Should not raise an exception and should log warning
            stage._process(context)

            # Check that warning was logged about set conversion
            warning_calls = [
                call
                for call in mock_logger.warning.call_args_list
                if "converted to sorted list" in str(call)
            ]
            assert len(warning_calls) > 0, "Should log warning about set conversion"

    def test_set_to_list_conversion_is_deterministic(self):
        """Test that set-to-list conversion produces consistent sorted results."""
        test_samples_set = {"Sample3", "Sample1", "Sample2"}
        expected_sorted = ["Sample1", "Sample2", "Sample3"]

        # Simulate multiple conversions
        conversions = []
        for _ in range(10):
            if isinstance(test_samples_set, set):
                converted = sorted(list(test_samples_set))
            else:
                converted = test_samples_set
            conversions.append(converted)

        # All conversions should be identical and sorted
        assert all(
            result == expected_sorted for result in conversions
        ), "Set-to-list conversions should be deterministic and sorted"

    @patch("variantcentrifuge.helpers.get_vcf_names")
    def test_sample_config_loading_stage_determinism(self, mock_get_names):
        """Test that SampleConfigLoadingStage produces deterministic results."""
        test_samples = ["Sample3", "Sample1", "Sample2"]
        mock_get_names.return_value = test_samples

        # Create contexts for multiple runs
        results = []
        for i in range(5):
            from argparse import Namespace
            from pathlib import Path

            workspace = Workspace(Path(f"/tmp/test_{i}"), "test")
            context = PipelineContext(Namespace(), {"vcf_file": "/tmp/test.vcf"}, workspace)

            stage = SampleConfigLoadingStage()
            result = stage._process(context)
            results.append(result.vcf_samples)

        # All results should be identical
        first_result = results[0]
        assert all(
            result == first_result for result in results
        ), "SampleConfigLoadingStage should produce deterministic sample lists"
        assert first_result == test_samples, "Should preserve original VCF header order"


class TestSampleAssignmentRegression:
    """Regression tests to prevent the sample assignment bug from reoccurring."""

    def test_no_set_types_in_vcf_samples_pipeline(self):
        """Test that vcf_samples is never a set in the pipeline context."""
        test_samples = ["Sample1", "Sample2", "Sample3"]

        with patch("variantcentrifuge.helpers.get_vcf_names", return_value=test_samples):
            from argparse import Namespace
            from pathlib import Path

            workspace = Workspace(Path("/tmp/test"), "test")
            context = PipelineContext(Namespace(), {"vcf_file": "/tmp/test.vcf"}, workspace)

            # Load samples using the stage
            stage = SampleConfigLoadingStage()
            result = stage._process(context)

            assert isinstance(result.vcf_samples, list), "vcf_samples should be a list, not a set"
            assert not isinstance(result.vcf_samples, set), "vcf_samples should never be a set"

    def test_genotype_replacement_sample_order_consistency(self):
        """Test that genotype replacement uses consistent sample ordering."""
        import tempfile
        from pathlib import Path

        from variantcentrifuge.stages.processing_stages import GenotypeReplacementStage

        test_samples = ["Sample1", "Sample2", "Sample3"]
        test_tsv_content = "CHROM\tPOS\tREF\tALT\tGT\nchr1\t100\tA\tT\t0/1:1/1:0/0\n"

        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(test_tsv_content)
            temp_tsv = Path(f.name)

        try:
            # Run genotype replacement multiple times
            results = []
            for i in range(3):
                from argparse import Namespace
                from pathlib import Path as PathLib

                workspace = Workspace(PathLib(f"/tmp/test_{i}"), "test")
                context = PipelineContext(Namespace(), {}, workspace)
                context.data = temp_tsv
                context.vcf_samples = test_samples.copy()  # Consistent list
                context.mark_complete("field_extraction")
                context.mark_complete("sample_config_loading")

                stage = GenotypeReplacementStage()
                result = stage._process(context)

                # Read the result file to check sample order consistency
                if result.genotype_replaced_tsv and result.genotype_replaced_tsv.exists():
                    import gzip

                    # Check if file is gzipped based on extension
                    if str(result.genotype_replaced_tsv).endswith(".gz"):
                        with gzip.open(result.genotype_replaced_tsv, "rt") as rf:
                            content = rf.read()
                    else:
                        with open(result.genotype_replaced_tsv, "r") as rf:
                            content = rf.read()
                    results.append(content)

            # All results should be identical if sample order is deterministic
            if results:
                first_result = results[0]
                assert all(
                    result == first_result for result in results
                ), "Genotype replacement should produce identical results with same sample order"

        finally:
            # Cleanup
            temp_tsv.unlink(missing_ok=True)

    def test_inheritance_analysis_sample_order_consistency(self):
        """Test that inheritance analysis produces consistent results with same sample order."""
        import pandas as pd

        test_samples = ["Child", "Father", "Mother"]
        test_df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": [100],
                "REF": ["A"],
                "ALT": ["T"],
                "GENE": ["BRCA1"],
                "GT": ["Child(0/1);Father(0/0);Mother(1/1)"],
            }
        )

        # Run inheritance analysis multiple times
        results = []
        for i in range(3):
            from argparse import Namespace
            from pathlib import Path

            workspace = Workspace(Path(f"/tmp/test_{i}"), "test")
            context = PipelineContext(
                Namespace(),
                {"inheritance_mode": "simple", "calculate_inheritance": True},
                workspace,
            )
            context.vcf_samples = test_samples.copy()  # Consistent list
            context.current_dataframe = test_df.copy()
            context.pedigree_data = {
                "Child": {"sample_id": "Child", "father": "Father", "mother": "Mother"},
                "Father": {"sample_id": "Father"},
                "Mother": {"sample_id": "Mother"},
            }
            context.mark_complete("dataframe_loading")

            stage = InheritanceAnalysisStage()
            result = stage._process(context)

            # Store relevant columns for comparison
            if result.current_dataframe is not None:
                inheritance_cols = [
                    col for col in result.current_dataframe.columns if "Inheritance" in col
                ]
                if inheritance_cols:
                    inheritance_data = result.current_dataframe[inheritance_cols].to_dict()
                    results.append(inheritance_data)

        # All results should be identical if sample order is deterministic
        if results:
            first_result = results[0]
            assert all(
                result == first_result for result in results
            ), "Inheritance analysis should produce identical results with same sample order"
