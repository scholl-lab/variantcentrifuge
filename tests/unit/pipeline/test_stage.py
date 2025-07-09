"""Unit tests for Stage abstract base class."""

from unittest.mock import Mock, patch

import pytest

from variantcentrifuge.pipeline_core import PipelineContext, Stage


class ConcreteStage(Stage):
    """Concrete implementation for testing."""

    def __init__(self, name="test_stage", dependencies=None, parallel_safe=False):
        self._name = name
        self._dependencies = dependencies or set()
        self._parallel_safe = parallel_safe
        self.process_called = False
        self.pre_execute_called = False
        self.post_execute_called = False

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

    def _process(self, context):
        self.process_called = True
        return context

    def _pre_execute(self, context):
        self.pre_execute_called = True

    def _post_execute(self, context):
        self.post_execute_called = True


class FailingStage(Stage):
    """Stage that fails during processing."""

    @property
    def name(self):
        """Return the stage name."""
        return "failing_stage"

    def _process(self, context):
        raise ValueError("Stage failed!")


class TestStage:
    """Test suite for Stage base class."""

    @pytest.fixture
    def context(self):
        """Create a mock PipelineContext."""
        context = Mock(spec=PipelineContext)
        context.is_complete = Mock(return_value=False)
        context.mark_complete = Mock()
        context.checkpoint_state = None  # Explicitly set to None to avoid checkpoint skipping
        return context

    def test_stage_properties(self):
        """Test stage property defaults."""
        stage = ConcreteStage()
        assert stage.name == "test_stage"
        assert stage.dependencies == set()
        assert stage.description == "Stage: test_stage"
        assert stage.parallel_safe is False
        assert stage.estimated_runtime == 1.0

    def test_stage_with_dependencies(self):
        """Test stage with dependencies."""
        stage = ConcreteStage(dependencies={"dep1", "dep2"})
        assert stage.dependencies == {"dep1", "dep2"}

    def test_successful_execution(self, context):
        """Test successful stage execution."""
        stage = ConcreteStage()

        result = stage(context)

        assert stage.pre_execute_called
        assert stage.process_called
        assert stage.post_execute_called
        assert result == context
        context.mark_complete.assert_called_once_with("test_stage")

    def test_dependency_validation(self, context):
        """Test dependency validation."""
        stage = ConcreteStage(dependencies={"dep1", "dep2"})

        # Mock missing dependencies
        context.is_complete.side_effect = lambda x: x == "dep1"

        with pytest.raises(RuntimeError) as exc_info:
            stage(context)

        assert "requires these stages to complete first: dep2" in str(exc_info.value)

    def test_skip_if_already_complete(self, context):
        """Test skipping execution if already complete."""
        stage = ConcreteStage()

        # Mock stage as already complete
        context.is_complete.side_effect = lambda x: x == "test_stage"

        result = stage(context)

        # Process should not be called
        assert not stage.process_called
        assert result == context
        context.mark_complete.assert_not_called()

    def test_stage_failure(self, context):
        """Test stage failure handling."""
        stage = FailingStage()

        with pytest.raises(ValueError) as exc_info:
            stage(context)

        assert "Stage failed!" in str(exc_info.value)
        context.mark_complete.assert_not_called()

    @patch("variantcentrifuge.pipeline_core.stage.logger")
    def test_logging(self, mock_logger, context):
        """Test stage execution logging."""
        stage = ConcreteStage()
        stage(context)

        # Check logging calls
        assert mock_logger.info.called
        info_calls = [call[0][0] for call in mock_logger.info.call_args_list]
        assert any("Executing Stage: test_stage" in call for call in info_calls)
        assert any("completed successfully" in call for call in info_calls)

    def test_parallel_safe_property(self):
        """Test parallel_safe property."""
        safe_stage = ConcreteStage(parallel_safe=True)
        unsafe_stage = ConcreteStage(parallel_safe=False)

        assert safe_stage.parallel_safe is True
        assert unsafe_stage.parallel_safe is False

    def test_input_output_files(self, context):
        """Test default input/output file methods."""
        stage = ConcreteStage()

        assert stage.get_input_files(context) == []
        assert stage.get_output_files(context) == []

    def test_validate_prerequisites(self, context):
        """Test prerequisite validation."""
        stage = ConcreteStage()

        # Should not raise by default
        stage.validate_prerequisites(context)

    def test_repr(self):
        """Test string representation."""
        stage1 = ConcreteStage()
        assert repr(stage1) == "ConcreteStage(name='test_stage')"

        stage2 = ConcreteStage(dependencies={"dep1"})
        assert "depends_on={'dep1'}" in repr(stage2)
