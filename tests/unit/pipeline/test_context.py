"""Unit tests for PipelineContext."""

import argparse
from datetime import datetime
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core import PipelineContext, Workspace


class TestPipelineContext:
    """Test suite for PipelineContext."""

    @pytest.fixture
    def mock_args(self):
        """Create mock command-line arguments."""
        args = Mock(spec=argparse.Namespace)
        args.vcf_file = "test.vcf"
        args.output_dir = "/tmp/test"
        args.gene_name = "BRCA1"
        return args

    @pytest.fixture
    def mock_workspace(self):
        """Create mock workspace."""
        return Mock(spec=Workspace)

    @pytest.fixture
    def context(self, mock_args, mock_workspace):
        """Create a test PipelineContext."""
        return PipelineContext(args=mock_args, config={"test": "value"}, workspace=mock_workspace)

    def test_initialization(self, mock_args, mock_workspace):
        """Test PipelineContext initialization."""
        context = PipelineContext(args=mock_args, config={"key": "value"}, workspace=mock_workspace)

        assert context.args == mock_args
        assert context.config == {"key": "value"}
        assert context.workspace == mock_workspace
        assert isinstance(context.start_time, datetime)
        assert context.data is None
        assert len(context.completed_stages) == 0
        assert len(context.stage_results) == 0

    def test_mark_complete(self, context):
        """Test marking a stage as complete."""
        context.mark_complete("test_stage")
        assert context.is_complete("test_stage")
        assert "test_stage" in context.completed_stages

        # With result
        context.mark_complete("stage_with_result", result={"data": 123})
        assert context.is_complete("stage_with_result")
        assert context.get_result("stage_with_result") == {"data": 123}

    def test_is_complete(self, context):
        """Test checking stage completion."""
        assert not context.is_complete("unknown_stage")

        context.mark_complete("known_stage")
        assert context.is_complete("known_stage")

    def test_get_result(self, context):
        """Test getting stage results."""
        assert context.get_result("unknown_stage") is None

        context.mark_complete("stage1", result="result1")
        assert context.get_result("stage1") == "result1"

    def test_update_data(self, context):
        """Test updating primary data artifact."""
        assert context.data is None

        context.update_data("new_data.tsv")
        assert context.data == "new_data.tsv"

        df = pd.DataFrame({"col": [1, 2, 3]})
        context.update_data(df)
        assert isinstance(context.data, pd.DataFrame)

    def test_add_report_path(self, context):
        """Test adding report paths."""
        assert len(context.report_paths) == 0

        context.add_report_path("html", Path("/tmp/report.html"))
        assert context.report_paths["html"] == Path("/tmp/report.html")

        context.add_report_path("excel", Path("/tmp/report.xlsx"))
        assert len(context.report_paths) == 2

    def test_get_execution_time(self, context):
        """Test execution time calculation."""
        # Mock the start time to be in the past
        past_time = datetime.now()
        context.start_time = past_time

        # Small delay to ensure time difference
        import time

        time.sleep(0.1)

        exec_time = context.get_execution_time()
        assert exec_time > 0
        assert isinstance(exec_time, float)

    def test_thread_safety(self, context):
        """Test thread-safe operations."""
        import threading

        def mark_stages():
            for i in range(100):
                context.mark_complete(f"stage_{i}", result=i)

        # Run multiple threads
        threads = []
        for _ in range(5):
            t = threading.Thread(target=mark_stages)
            threads.append(t)
            t.start()

        for t in threads:
            t.join()

        # Verify all stages marked
        assert len(context.completed_stages) == 100
        for i in range(100):
            assert context.is_complete(f"stage_{i}")
            assert context.get_result(f"stage_{i}") == i

    def test_repr(self, context):
        """Test string representation."""
        repr_str = repr(context)
        assert "PipelineContext" in repr_str
        assert "stages_completed=0" in repr_str
        assert "has_data=False" in repr_str

        context.mark_complete("stage1")
        context.update_data("some_data")
        repr_str = repr(context)
        assert "stages_completed=1" in repr_str
        assert "has_data=True" in repr_str
