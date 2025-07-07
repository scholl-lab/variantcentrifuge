"""Test PipelineContext pickling for ProcessPoolExecutor support."""

import pickle
from argparse import Namespace

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace


def process_context_in_subprocess(ctx):
    """Function to run in subprocess for testing pickling."""
    ctx.mark_complete("subprocess_stage")
    return ctx.is_complete("subprocess_stage")


class TestPipelineContextPickling:
    """Test that PipelineContext can be pickled and unpickled correctly."""

    @pytest.fixture
    def workspace(self, tmp_path):
        """Create a real workspace for testing."""
        return Workspace(output_dir=tmp_path, base_name="test")

    def test_basic_pickling(self, workspace):
        """Test basic pickling and unpickling of PipelineContext."""
        # Create context
        context = PipelineContext(
            args=Namespace(test=True),
            config={"key": "value", "threads": 4},
            workspace=workspace,
        )

        # Add some state
        context.mark_complete("test_stage", result={"data": "result"})
        context.vcf_samples = ["sample1", "sample2"]
        context.config["extra"] = "data"

        # Pickle and unpickle
        pickled = pickle.dumps(context)
        restored = pickle.loads(pickled)

        # Verify state was preserved
        assert restored.config == context.config
        assert restored.vcf_samples == context.vcf_samples
        assert restored.is_complete("test_stage")
        assert restored.get_result("test_stage") == {"data": "result"}

        # Verify lock was recreated
        assert hasattr(restored, "_lock")
        assert restored._lock is not None

    def test_thread_safety_after_unpickling(self, workspace):
        """Test that thread safety is maintained after unpickling."""
        context = PipelineContext(
            args=Namespace(),
            config={},
            workspace=workspace,
        )

        # Pickle and unpickle
        pickled = pickle.dumps(context)
        restored = pickle.loads(pickled)

        # Test thread-safe operations work
        restored.mark_complete("stage1")
        restored.mark_complete("stage2", result="test")

        assert restored.is_complete("stage1")
        assert restored.is_complete("stage2")
        assert restored.get_result("stage2") == "test"

    def test_complex_state_pickling(self, workspace):
        """Test pickling with complex state including DataFrames."""
        import pandas as pd

        context = PipelineContext(
            args=Namespace(complex=True),
            config={"nested": {"data": [1, 2, 3]}},
            workspace=workspace,
        )

        # Add complex state
        context.current_dataframe = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
        context.statistics = {"total_variants": 1000, "filtered": 500, "genes": ["GENE1", "GENE2"]}
        context.report_paths = {
            "html": workspace.output_dir / "report.html",
            "excel": workspace.output_dir / "report.xlsx",
        }

        # Pickle and unpickle
        pickled = pickle.dumps(context)
        restored = pickle.loads(pickled)

        # Verify complex state
        assert restored.config == context.config
        pd.testing.assert_frame_equal(restored.current_dataframe, context.current_dataframe)
        assert restored.statistics == context.statistics
        assert restored.report_paths == context.report_paths

    def test_concurrent_access_after_unpickling(self, workspace):
        """Test concurrent access patterns after unpickling."""
        from concurrent.futures import ThreadPoolExecutor

        context = PipelineContext(
            args=Namespace(),
            config={},
            workspace=workspace,
        )

        # Pickle and unpickle
        pickled = pickle.dumps(context)
        restored = pickle.loads(pickled)

        # Test concurrent marking of stages
        def mark_stage(ctx, stage_name):
            ctx.mark_complete(stage_name, result=f"result_{stage_name}")
            return ctx.is_complete(stage_name)

        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = []
            for i in range(10):
                future = executor.submit(mark_stage, restored, f"stage_{i}")
                futures.append(future)

            # All should complete successfully
            results = [f.result() for f in futures]
            assert all(results)

        # Verify all stages marked
        assert len(restored.completed_stages) == 10

    def test_pickle_with_process_pool(self, workspace):
        """Test that context can be used with ProcessPoolExecutor."""
        from concurrent.futures import ProcessPoolExecutor

        context = PipelineContext(
            args=Namespace(),
            config={"test": "data"},
            workspace=workspace,
        )

        # This should now work without pickling errors
        with ProcessPoolExecutor(max_workers=1) as executor:
            future = executor.submit(process_context_in_subprocess, context)
            result = future.result()

        assert result is True
