"""Tests for subtask timing functionality."""

import time
import unittest
from pathlib import Path
from typing import Set

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.runner import PipelineRunner
from variantcentrifuge.pipeline_core.stage import Stage
from variantcentrifuge.pipeline_core.workspace import Workspace


class StageWithSubtasks(Stage):
    """Test stage that uses subtask timing."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "test_stage_with_subtasks"

    @property
    def description(self) -> str:
        """Return the stage description."""
        return "Test stage with subtasks"

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process with subtask timing."""
        # Subtask 1
        t1 = self._start_subtask("subtask1")
        time.sleep(0.1)
        self._end_subtask("subtask1", t1)

        # Subtask 2
        t2 = self._start_subtask("subtask2")
        time.sleep(0.2)
        self._end_subtask("subtask2", t2)

        return context


class SimpleStage(Stage):
    """Simple stage without subtasks."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "simple_stage"

    @property
    def description(self) -> str:
        """Return the stage description."""
        return "Simple stage"

    @property
    def dependencies(self) -> Set[str]:
        """Return the stage dependencies."""
        return {"test_stage_with_subtasks"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process the stage."""
        time.sleep(0.1)
        return context


class TestSubtaskTiming(unittest.TestCase):
    """Test subtask timing functionality."""

    def setUp(self):
        """Set up test environment."""
        import argparse

        self.args = argparse.Namespace()
        self.config = {"test": "config"}
        self.workspace = Workspace(Path("/tmp/test_subtask"), "test")
        self.context = PipelineContext(args=self.args, config=self.config, workspace=self.workspace)

    def test_subtask_timing_recorded(self):
        """Test that subtask times are recorded."""
        stage = StageWithSubtasks()
        stage(self.context)

        # Check subtask times were recorded
        subtask_times = stage.subtask_times
        self.assertIn("subtask1", subtask_times)
        self.assertIn("subtask2", subtask_times)

        # Verify approximate timing
        self.assertGreater(subtask_times["subtask1"], 0.09)
        self.assertLess(subtask_times["subtask1"], 0.15)
        self.assertGreater(subtask_times["subtask2"], 0.19)
        self.assertLess(subtask_times["subtask2"], 0.25)

    def test_stage_without_subtasks(self):
        """Test stage without subtasks has empty subtask times."""
        stage = SimpleStage()
        # Mark dependency as complete
        self.context.mark_complete("test_stage_with_subtasks")
        stage(self.context)

        # Check no subtask times recorded
        self.assertEqual(stage.subtask_times, {})

    def test_runner_captures_subtask_times(self):
        """Test that runner captures subtask times from stages."""
        stages = [StageWithSubtasks(), SimpleStage()]
        runner = PipelineRunner()

        # Run pipeline
        runner.run(stages, self.context)

        # Check runner captured subtask times
        self.assertIn("test_stage_with_subtasks", runner._subtask_times)
        self.assertNotIn("simple_stage", runner._subtask_times)

        # Verify subtask times content
        subtasks = runner._subtask_times["test_stage_with_subtasks"]
        self.assertIn("subtask1", subtasks)
        self.assertIn("subtask2", subtasks)
