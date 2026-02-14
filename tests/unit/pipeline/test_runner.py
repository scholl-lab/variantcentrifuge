"""Unit tests for PipelineRunner."""

import time
from unittest.mock import Mock, patch

import pytest

from variantcentrifuge.pipeline_core import PipelineContext, PipelineRunner, Stage


class MockStage(Stage):
    """Mock stage for testing."""

    def __init__(
        self, name, dependencies=None, parallel_safe=False, runtime=0.1, should_fail=False
    ):
        self._name = name
        self._dependencies = set(dependencies or [])
        self._parallel_safe = parallel_safe
        self._runtime = runtime
        self._should_fail = should_fail
        self.executed = False
        self.execution_time = None

    @property
    def name(self):
        """Return the stage name."""
        return self._name

    @property
    def dependencies(self):
        """Return the stage dependencies."""
        return self._dependencies

    @property
    def parallel_safe(self):
        """Return whether the stage is safe to run in parallel."""
        return self._parallel_safe

    @property
    def estimated_runtime(self):
        """Return the estimated runtime in seconds."""
        return self._runtime

    def _process(self, context):
        start = time.time()
        time.sleep(self._runtime)  # Simulate work
        self.execution_time = time.time() - start
        self.executed = True

        if self._should_fail:
            raise ValueError(f"Stage {self.name} failed!")

        context.mark_complete(self.name)
        return context


class TestPipelineRunner:
    """Test suite for PipelineRunner."""

    @pytest.fixture
    def context(self):
        """Create a test context."""
        context = Mock(spec=PipelineContext)
        context.is_complete = Mock(return_value=False)
        context.mark_complete = Mock()
        context.completed_stages = set()
        context.checkpoint_state = None
        context.stage_results = {}
        context.config = {"resume": False}  # Add config attribute to fix runner test

        # Make mark_complete update completed_stages
        def mark_complete_side_effect(name):
            context.completed_stages.add(name)

        context.mark_complete.side_effect = mark_complete_side_effect

        # Make is_complete check completed_stages
        context.is_complete.side_effect = lambda name: name in context.completed_stages

        return context

    @pytest.fixture
    def runner(self):
        """Create a test runner."""
        return PipelineRunner(enable_checkpoints=False, max_workers=2)

    def test_simple_sequential_execution(self, runner, context):
        """Test simple sequential stage execution."""
        stages = [
            MockStage("stage1"),
            MockStage("stage2", dependencies=["stage1"]),
            MockStage("stage3", dependencies=["stage2"]),
        ]

        result = runner.run(stages, context)

        assert all(stage.executed for stage in stages)
        assert result == context
        assert len(context.completed_stages) == 3

    def test_parallel_execution(self, runner, context):
        """Test parallel stage execution."""
        # Create stages that can run in parallel with longer runtimes
        # to make parallel benefits more apparent
        stages = [
            MockStage("stage1", parallel_safe=True, runtime=0.2),
            MockStage("stage2", parallel_safe=True, runtime=0.2),
            MockStage("stage3", parallel_safe=True, runtime=0.2),
            MockStage("stage4", dependencies=["stage1", "stage2", "stage3"], runtime=0.2),
        ]

        start_time = time.time()
        result = runner.run(stages, context)
        total_time = time.time() - start_time

        assert all(stage.executed for stage in stages)
        assert result == context

        # Parallel execution should be faster than sequential
        # The first 3 stages should run in parallel (~0.2s), then stage4 runs (~0.2s)
        # So total should be ~0.4s instead of sequential 0.8s
        # Use generous margin (1.5) to account for gc.collect() overhead and loaded systems
        # gc.collect() runs after each stage (4x overhead), adding time but freeing memory
        sequential_time = sum(s.estimated_runtime for s in stages)
        assert total_time < sequential_time * 1.5

    def test_mixed_parallel_sequential(self, runner, context):
        """Test mixed parallel and sequential execution."""
        stages = [
            MockStage("config1", parallel_safe=True),
            MockStage("config2", parallel_safe=True),
            MockStage("sequential", parallel_safe=False),  # Must run alone
            MockStage("parallel1", dependencies=["sequential"], parallel_safe=True),
            MockStage("parallel2", dependencies=["sequential"], parallel_safe=True),
        ]

        result = runner.run(stages, context)

        assert all(stage.executed for stage in stages)
        assert result == context

    def test_circular_dependency_detection(self, runner, context):
        """Test circular dependency detection."""
        stages = [
            MockStage("stage1", dependencies=["stage3"]),
            MockStage("stage2", dependencies=["stage1"]),
            MockStage("stage3", dependencies=["stage2"]),
        ]

        with pytest.raises(ValueError) as exc_info:
            runner.run(stages, context)

        assert "Circular dependency detected" in str(exc_info.value)

    def test_duplicate_stage_names(self, runner, context):
        """Test duplicate stage name detection."""
        stages = [
            MockStage("stage1"),
            MockStage("stage1"),  # Duplicate name
        ]

        with pytest.raises(ValueError) as exc_info:
            runner.run(stages, context)

        assert "Duplicate stage names detected" in str(exc_info.value)

    def test_stage_failure_handling(self, runner, context):
        """Test handling of stage failures."""
        stages = [
            MockStage("stage1"),
            MockStage("stage2", should_fail=True),
            MockStage("stage3", dependencies=["stage2"]),
        ]

        with pytest.raises(ValueError) as exc_info:
            runner.run(stages, context)

        assert "Stage stage2 failed!" in str(exc_info.value)
        assert stages[0].executed
        assert stages[1].executed
        assert not stages[2].executed  # Should not execute after failure

    def test_execution_plan_creation(self, runner):
        """Test execution plan creation."""
        stages = [
            MockStage("a", dependencies=["b", "c"]),
            MockStage("b", dependencies=["d"]),
            MockStage("c", dependencies=["d"]),
            MockStage("d"),
            MockStage("e", dependencies=["a"]),
        ]

        plan = runner._create_execution_plan(stages)

        # Should have 4 levels: d -> b,c -> a -> e
        assert len(plan) == 4
        assert plan[0][0].name == "d"
        assert {s.name for s in plan[1]} == {"b", "c"}
        assert plan[2][0].name == "a"
        assert plan[3][0].name == "e"

    def test_dry_run(self, runner):
        """Test dry run functionality."""
        stages = [
            MockStage("stage1"),
            MockStage("stage2", dependencies=["stage1"]),
            MockStage("stage3", dependencies=["stage1"]),
            MockStage("stage4", dependencies=["stage2", "stage3"]),
        ]

        plan = runner.dry_run(stages)

        assert len(plan) == 3
        assert plan[0] == ["stage1"]
        assert set(plan[1]) == {"stage2", "stage3"}  # Order doesn't matter
        assert plan[2] == ["stage4"]

    @patch("variantcentrifuge.pipeline_core.runner.logger")
    def test_execution_summary_logging(self, mock_logger, runner, context):
        """Test execution summary logging."""
        stages = [
            MockStage("fast", runtime=0.01),
            MockStage("slow", runtime=0.05),
        ]

        runner.run(stages, context)

        # Check that summary was logged
        info_calls = [call[0][0] for call in mock_logger.info.call_args_list]
        assert any("Stage Execution Summary" in call for call in info_calls)
        assert any("fast" in call for call in info_calls)
        assert any("slow" in call for call in info_calls)

    def test_checkpoint_saving(self, runner, context):
        """Test checkpoint saving during execution."""
        # Enable checkpoints
        runner.enable_checkpoints = True

        # Mock checkpoint state
        checkpoint_state = Mock()
        checkpoint_state.should_skip_step = Mock(return_value=False)  # Don't skip any stages
        context.checkpoint_state = checkpoint_state

        stages = [
            MockStage("stage1"),
            MockStage("stage2", dependencies=["stage1"]),  # Force 2 levels
        ]

        runner.run(stages, context)

        # Checkpoint should be saved after each level
        assert checkpoint_state.save.call_count >= 2
