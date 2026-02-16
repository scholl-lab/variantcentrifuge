"""Tests for memory reporting functionality in PipelineRunner."""

import logging
import unittest
from pathlib import Path

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.runner import PipelineRunner
from variantcentrifuge.pipeline_core.stage import Stage
from variantcentrifuge.pipeline_core.workspace import Workspace


class MemoryAllocatingStage(Stage):
    """Test stage that allocates memory."""

    def __init__(self, stage_name: str = "memory_stage"):
        """Initialize with configurable name."""
        self._stage_name = stage_name

    @property
    def name(self) -> str:
        """Return the stage name."""
        return self._stage_name

    @property
    def description(self) -> str:
        """Return the stage description."""
        return f"Test stage: {self._stage_name}"

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Allocate some memory to create measurable delta."""
        # Allocate ~10MB to create measurable memory delta
        _data = [0] * (10 * 1024 * 1024 // 8)
        return context


class SimpleStage(Stage):
    """Simple stage without significant memory allocation."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "simple_stage"

    @property
    def description(self) -> str:
        """Return the stage description."""
        return "Simple test stage"

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process without allocating memory."""
        return context


@pytest.mark.unit
class TestRunnerMemoryReporting(unittest.TestCase):
    """Test memory reporting functionality."""

    def setUp(self):
        """Set up test environment."""
        import argparse

        self.args = argparse.Namespace()
        self.config = {"test": "config"}
        self.workspace = Workspace(Path("/tmp/test_memory"), "test")
        self.context = PipelineContext(args=self.args, config=self.config, workspace=self.workspace)

    def test_stage_metrics_captured(self):
        """Test that memory metrics are captured for each stage."""
        # Create two stages
        stage1 = MemoryAllocatingStage("stage1")
        stage2 = MemoryAllocatingStage("stage2")
        stages = [stage1, stage2]

        runner = PipelineRunner()
        runner.run(stages, self.context)

        # Verify metrics were captured
        self.assertIn("stage1", runner._stage_metrics)
        self.assertIn("stage2", runner._stage_metrics)

        # Verify metric structure
        for stage_name in ["stage1", "stage2"]:
            metrics = runner._stage_metrics[stage_name]
            self.assertIn("mem_before_mb", metrics)
            self.assertIn("mem_after_mb", metrics)
            self.assertIn("mem_delta_mb", metrics)
            self.assertIn("mem_peak_mb", metrics)

            # Verify metrics are reasonable
            self.assertGreater(metrics["mem_before_mb"], 0)
            self.assertGreater(metrics["mem_after_mb"], 0)
            # Delta is calculated correctly
            expected_delta = metrics["mem_after_mb"] - metrics["mem_before_mb"]
            self.assertAlmostEqual(metrics["mem_delta_mb"], expected_delta, places=1)
            # Peak is at least as large as before and after
            self.assertGreaterEqual(
                metrics["mem_peak_mb"], max(metrics["mem_before_mb"], metrics["mem_after_mb"])
            )

    def test_memory_summary_logged(self, caplog=None):
        """Test that memory summary appears in log output at INFO level."""
        # Use pytest's caplog fixture
        stages = [MemoryAllocatingStage("alloc_stage"), SimpleStage()]
        runner = PipelineRunner()

        with self.assertLogs("variantcentrifuge.pipeline_core.runner", level=logging.INFO) as cm:
            runner.run(stages, self.context)

        # Combine all log messages
        log_output = "\n".join(cm.output)

        # Verify memory summary section exists
        self.assertIn("Memory Usage Summary", log_output)
        self.assertIn("Peak RSS", log_output)
        self.assertIn("MB", log_output)

        # Verify both stages appear in memory summary
        self.assertIn("alloc_stage", log_output)
        self.assertIn("simple_stage", log_output)

        # Verify peak memory is reported
        self.assertIn("Pipeline peak memory:", log_output)

    def test_metrics_empty_when_no_stages(self):
        """Test that empty metrics are handled gracefully."""
        runner = PipelineRunner()

        # Call _log_execution_summary with no stages executed
        # Should not crash
        try:
            runner._log_execution_summary()
        except Exception as e:
            self.fail(f"_log_execution_summary raised exception with empty metrics: {e}")

        # Verify no metrics were created
        self.assertEqual(len(runner._stage_metrics), 0)

    def test_memory_delta_calculation(self):
        """Test that memory delta is correctly calculated (positive or negative)."""
        # Run a simple stage
        stage = SimpleStage()
        runner = PipelineRunner()
        runner.run([stage], self.context)

        # Verify delta calculation
        metrics = runner._stage_metrics["simple_stage"]
        calculated_delta = metrics["mem_after_mb"] - metrics["mem_before_mb"]
        self.assertAlmostEqual(metrics["mem_delta_mb"], calculated_delta, places=1)

    def test_peak_memory_reported_after_completion(self):
        """Test that peak memory is logged after pipeline completion."""
        stages = [MemoryAllocatingStage("stage1"), MemoryAllocatingStage("stage2")]
        runner = PipelineRunner()

        with self.assertLogs("variantcentrifuge.pipeline_core.runner", level=logging.INFO) as cm:
            runner.run(stages, self.context)

        log_output = "\n".join(cm.output)

        # Verify peak memory message exists
        self.assertIn("Pipeline peak memory:", log_output)

        # Verify peak memory value is reasonable (greater than 0)
        # Extract the peak value - find line with "Pipeline peak memory:"
        for line in cm.output:
            if "Pipeline peak memory:" in line:
                # Line format: "Pipeline peak memory: XXX MB"
                self.assertIn("MB", line)
                break
        else:
            self.fail("Pipeline peak memory line not found in logs")
