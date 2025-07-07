"""Unit tests for error handling utilities."""

import logging
import pytest
import time
from pathlib import Path
from unittest.mock import Mock, patch

from variantcentrifuge.pipeline_core.error_handling import (
    PipelineError,
    ToolNotFoundError,
    FileFormatError,
    DataValidationError,
    StageExecutionError,
    retry_on_failure,
    graceful_error_handling,
    FileRecoveryStrategy,
    ToolRecoveryStrategy,
    ErrorRecoveryManager,
    handle_stage_error,
    validate_file_exists,
    validate_output_directory,
)


class TestExceptionClasses:
    """Test custom exception classes."""

    def test_pipeline_error(self):
        """Test PipelineError initialization."""
        error = PipelineError("Test error", stage="test_stage", details={"key": "value"})
        assert str(error) == "Test error"
        assert error.stage == "test_stage"
        assert error.details == {"key": "value"}

    def test_tool_not_found_error(self):
        """Test ToolNotFoundError."""
        error = ToolNotFoundError("bcftools", stage="extraction")
        assert "bcftools" in str(error)
        assert "not found in PATH" in str(error)
        assert error.stage == "extraction"
        assert error.details["tool"] == "bcftools"

    def test_file_format_error(self):
        """Test FileFormatError."""
        error = FileFormatError("/path/to/file.txt", "VCF", stage="loading")
        assert "/path/to/file.txt" in str(error)
        assert "VCF" in str(error)
        assert error.details["file"] == "/path/to/file.txt"
        assert error.details["expected_format"] == "VCF"

    def test_data_validation_error(self):
        """Test DataValidationError."""
        error = DataValidationError("Invalid value", "QUAL", stage="filtering")
        assert "Invalid value" in str(error)
        assert error.details["field"] == "QUAL"

    def test_stage_execution_error(self):
        """Test StageExecutionError."""
        original = ValueError("Original error")
        error = StageExecutionError("test_stage", original, recovery_attempted=True)
        assert "test_stage" in str(error)
        assert "Original error" in str(error)
        assert error.details["recovery_attempted"] is True
        assert error.original_error is original


class TestRetryDecorator:
    """Test retry_on_failure decorator."""

    def test_successful_execution(self):
        """Test successful execution without retries."""
        call_count = 0

        @retry_on_failure(max_attempts=3)
        def success_func():
            nonlocal call_count
            call_count += 1
            return "success"

        result = success_func()
        assert result == "success"
        assert call_count == 1

    def test_retry_then_success(self):
        """Test retry mechanism with eventual success."""
        call_count = 0

        @retry_on_failure(max_attempts=3, delay=0.1)
        def eventual_success():
            nonlocal call_count
            call_count += 1
            if call_count < 3:
                raise ValueError("Temporary failure")
            return "success"

        start_time = time.time()
        result = eventual_success()
        elapsed = time.time() - start_time

        assert result == "success"
        assert call_count == 3
        assert elapsed >= 0.2  # At least 2 delays of 0.1s

    def test_max_attempts_exceeded(self):
        """Test failure after max attempts."""
        call_count = 0

        @retry_on_failure(max_attempts=3, delay=0.01)
        def always_fails():
            nonlocal call_count
            call_count += 1
            raise ValueError("Persistent failure")

        with pytest.raises(ValueError, match="Persistent failure"):
            always_fails()

        assert call_count == 3

    def test_specific_exceptions(self):
        """Test retry only on specific exceptions."""

        @retry_on_failure(max_attempts=3, exceptions=(ValueError,))
        def raises_type_error():
            raise TypeError("Wrong type")

        # Should not retry TypeError
        with pytest.raises(TypeError):
            raises_type_error()

    def test_exponential_backoff(self):
        """Test exponential backoff."""
        delays = []

        @retry_on_failure(max_attempts=3, delay=0.1, backoff=2.0)
        def track_delays():
            nonlocal delays
            if delays:
                delays.append(time.time())
            else:
                delays.append(time.time())
            raise ValueError("Force retry")

        with pytest.raises(ValueError):
            track_delays()

        # Check delays increase
        assert len(delays) == 3
        if len(delays) > 2:
            delay1 = delays[1] - delays[0]
            delay2 = delays[2] - delays[1]
            assert delay2 > delay1 * 1.5  # Allow some tolerance


class TestGracefulErrorHandling:
    """Test graceful_error_handling context manager."""

    def test_successful_execution(self):
        """Test successful execution."""
        with graceful_error_handling("test_stage"):
            result = 1 + 1
        assert result == 2

    def test_pipeline_error_passthrough(self):
        """Test PipelineError passes through."""
        with pytest.raises(PipelineError):
            with graceful_error_handling("test_stage"):
                raise PipelineError("Test error")

    def test_file_not_found_conversion(self):
        """Test FileNotFoundError conversion."""
        with pytest.raises(PipelineError) as exc_info:
            with graceful_error_handling("test_stage"):
                raise FileNotFoundError("missing.txt")

        assert "Required file not found" in str(exc_info.value)
        assert exc_info.value.stage == "test_stage"

    def test_permission_error_conversion(self):
        """Test PermissionError conversion."""
        with pytest.raises(PipelineError) as exc_info:
            with graceful_error_handling("test_stage"):
                raise PermissionError("Cannot write")

        assert "Permission denied" in str(exc_info.value)

    def test_generic_error_conversion(self):
        """Test generic error conversion."""
        with pytest.raises(StageExecutionError) as exc_info:
            with graceful_error_handling("test_stage"):
                raise RuntimeError("Unexpected")

        assert exc_info.value.stage == "test_stage"
        assert isinstance(exc_info.value.original_error, RuntimeError)

    def test_no_reraise(self):
        """Test with reraise=False."""
        result = None
        with graceful_error_handling("test_stage", reraise=False, fallback_value="fallback"):
            raise ValueError("Error")
            result = "should not reach"

        # Should not set result due to error
        assert result is None


class TestRecoveryStrategies:
    """Test error recovery strategies."""

    def test_file_recovery_strategy(self):
        """Test FileRecoveryStrategy."""
        strategy = FileRecoveryStrategy()

        # Test FileNotFoundError recovery with alternatives
        error = FileNotFoundError("missing.txt")
        context = {
            "original_path": "/path/to/missing.txt",
            "alternative_paths": ["/nonexistent.txt", "/tmp"],
        }

        with patch("pathlib.Path.exists") as mock_exists:
            mock_exists.side_effect = [False, True]
            assert strategy.can_recover(error, context)
            result = strategy.recover(error, context)
            assert result == "/tmp"

    def test_file_recovery_create(self, tmp_path):
        """Test file creation recovery."""
        strategy = FileRecoveryStrategy()
        missing_file = tmp_path / "subdir" / "missing.txt"

        error = FileNotFoundError(str(missing_file))
        context = {
            "original_path": str(missing_file),
            "can_create": True,
        }

        assert strategy.can_recover(error, context)
        result = strategy.recover(error, context)
        assert result == str(missing_file)
        assert missing_file.exists()

    def test_tool_recovery_strategy(self):
        """Test ToolRecoveryStrategy."""
        strategy = ToolRecoveryStrategy()

        # Test tool not found with alternatives
        error = ToolNotFoundError("bcftools")
        context = {"alternative_tools": ["samtools", "vcftools"]}

        with patch("shutil.which") as mock_which:
            mock_which.side_effect = [None, "vcftools"]
            assert strategy.can_recover(error, context)
            result = strategy.recover(error, context)
            assert result == "vcftools"


class TestErrorRecoveryManager:
    """Test ErrorRecoveryManager."""

    def test_recovery_manager(self):
        """Test error recovery manager."""
        manager = ErrorRecoveryManager()

        # Test with FileNotFoundError
        error = FileNotFoundError("missing.txt")
        context = {
            "alternative_paths": ["/exists.txt"],
        }

        with patch("pathlib.Path.exists", return_value=True):
            result = manager.attempt_recovery(error, context)
            assert result == "/exists.txt"

    def test_custom_strategy(self):
        """Test adding custom strategy."""
        manager = ErrorRecoveryManager()

        # Create custom strategy
        class CustomStrategy:
            def can_recover(self, error, context):
                return isinstance(error, ValueError) and "custom" in context

            def recover(self, error, context):
                return "custom_recovery"

        manager.add_strategy(CustomStrategy())

        error = ValueError("test")
        result = manager.attempt_recovery(error, {"custom": True})
        assert result == "custom_recovery"


class TestValidationFunctions:
    """Test validation helper functions."""

    def test_validate_file_exists(self, tmp_path):
        """Test file validation."""
        # Create test file
        test_file = tmp_path / "test.txt"
        test_file.write_text("content")

        # Test valid file
        result = validate_file_exists(test_file, "test_stage")
        assert result == test_file

        # Test missing file
        with pytest.raises(FileNotFoundError):
            validate_file_exists(tmp_path / "missing.txt", "test_stage")

        # Test directory instead of file
        with pytest.raises(FileFormatError):
            validate_file_exists(tmp_path, "test_stage")

    def test_validate_output_directory(self, tmp_path):
        """Test output directory validation."""
        # Test existing directory
        result = validate_output_directory(tmp_path, "test_stage")
        assert result == tmp_path

        # Test creating directory
        new_dir = tmp_path / "new" / "nested"
        result = validate_output_directory(new_dir, "test_stage", create=True)
        assert result == new_dir
        assert new_dir.exists()

        # Test no create flag
        missing_dir = tmp_path / "missing"
        with pytest.raises(FileNotFoundError):
            validate_output_directory(missing_dir, "test_stage", create=False)

        # Test file instead of directory
        file_path = tmp_path / "file.txt"
        file_path.write_text("content")
        with pytest.raises(FileFormatError):
            validate_output_directory(file_path, "test_stage")


class TestHandleStageError:
    """Test handle_stage_error function."""

    def test_successful_recovery(self):
        """Test successful error recovery."""
        error = FileNotFoundError("missing.txt")
        context = {"alternative_paths": ["/tmp"]}

        with patch("pathlib.Path.exists", return_value=True):
            result = handle_stage_error("test_stage", error, context, reraise=False)
            assert result == "/tmp"

    def test_failed_recovery_with_reraise(self):
        """Test failed recovery with reraise."""
        error = ValueError("Unrecoverable")

        with pytest.raises(StageExecutionError) as exc_info:
            handle_stage_error("test_stage", error, reraise=True)

        assert exc_info.value.stage == "test_stage"
        assert exc_info.value.details["recovery_attempted"] is True

    def test_failed_recovery_no_reraise(self):
        """Test failed recovery without reraise."""
        error = ValueError("Unrecoverable")
        result = handle_stage_error("test_stage", error, reraise=False)
        assert result is None
