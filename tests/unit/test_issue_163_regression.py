"""Regression tests for Issue #163:
- Module-level logging.getLogger(__name__) in runner
- Guard psutil gracefully for memory tracking in runner (import failure and runtime errors)
"""

from __future__ import annotations

import logging
from unittest.mock import Mock, patch

import pytest

import variantcentrifuge.pipeline_core.runner as runner_module
from variantcentrifuge.pipeline_core import PipelineContext, PipelineRunner, Stage


class SimpleTestStage(Stage):
    """Simple stage for runner execution tests."""

    @property
    def name(self) -> str:
        return "simple_stage"

    @property
    def dependencies(self) -> set[str]:
        return set()

    def _process(self, context: PipelineContext) -> PipelineContext:
        context.mark_complete(self.name)
        return context


@pytest.mark.unit
class TestRunnerLogger:
    """Verify module-level logger conventions in runner.py."""

    def test_runner_module_logger(self) -> None:
        """runner.py must have a module-level logger named after its module (__name__)."""
        assert hasattr(runner_module, "logger")
        assert isinstance(runner_module.logger, logging.Logger)
        assert runner_module.logger.name == "variantcentrifuge.pipeline_core.runner"


@pytest.mark.unit
class TestRunnerPsutilGracefulGuarding:
    """Verify psutil is guarded gracefully when unavailable or failing."""

    def test_runner_executes_when_psutil_none(self) -> None:
        """PipelineRunner should complete execution without error when psutil is None."""
        runner = PipelineRunner(enable_checkpoints=False, max_workers=1)
        stage = SimpleTestStage()
        context = PipelineContext(args=Mock(), config={}, workspace=Mock())

        with patch.object(runner_module, "psutil", None):
            result_context = runner.run([stage], context)

        assert result_context.is_complete("simple_stage")
        # Memory metrics should either be empty or contain graceful defaults
        assert stage.name not in runner._stage_metrics or runner._stage_metrics[stage.name] == {}

    def test_runner_handles_psutil_runtime_exception(self) -> None:
        """PipelineRunner should gracefully catch runtime exceptions during memory measurement."""
        runner = PipelineRunner(enable_checkpoints=False, max_workers=1)
        stage = SimpleTestStage()
        context = PipelineContext(args=Mock(), config={}, workspace=Mock())

        mock_psutil = Mock()
        mock_psutil.Process.side_effect = RuntimeError("Process memory info unavailable")

        with patch.object(runner_module, "psutil", mock_psutil):
            result_context = runner.run([stage], context)

        assert result_context.is_complete("simple_stage")

    def test_runner_tracks_memory_when_psutil_available(self) -> None:
        """PipelineRunner tracks memory when psutil works normally."""
        runner = PipelineRunner(enable_checkpoints=False, max_workers=1)
        stage = SimpleTestStage()
        context = PipelineContext(args=Mock(), config={}, workspace=Mock())

        mock_process = Mock()
        mock_process.memory_info.return_value = Mock(rss=100 * 1024 * 1024)  # 100 MB
        mock_psutil = Mock()
        mock_psutil.Process.return_value = mock_process

        with patch.object(runner_module, "psutil", mock_psutil):
            result_context = runner.run([stage], context)

        assert result_context.is_complete("simple_stage")
        assert stage.name in runner._stage_metrics
        assert runner._stage_metrics[stage.name]["mem_before_mb"] == pytest.approx(100.0)
        assert runner._stage_metrics[stage.name]["mem_after_mb"] == pytest.approx(100.0)
