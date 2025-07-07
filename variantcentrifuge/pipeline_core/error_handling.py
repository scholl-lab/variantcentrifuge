"""
Enhanced error handling utilities for the refactored pipeline.

This module provides error handling utilities including:
- Custom exception classes for different error types
- Retry decorators for transient failures
- Context managers for graceful error handling
- Error recovery strategies
"""

import logging
import time
from contextlib import contextmanager
from functools import wraps
from pathlib import Path
from typing import Callable, Dict, List, Optional, Type, Union

logger = logging.getLogger(__name__)


class PipelineError(Exception):
    """Base exception for all pipeline errors."""

    def __init__(self, message: str, stage: Optional[str] = None, details: Optional[Dict] = None):
        """Initialize pipeline error.

        Parameters
        ----------
        message : str
            Error message
        stage : str, optional
            Stage where error occurred
        details : dict, optional
            Additional error details
        """
        super().__init__(message)
        self.stage = stage
        self.details = details or {}


class ToolNotFoundError(PipelineError):
    """Raised when a required external tool is not found."""

    def __init__(self, tool: str, stage: Optional[str] = None):
        """Initialize tool not found error."""
        message = f"Required tool '{tool}' not found in PATH"
        super().__init__(message, stage, {"tool": tool})


class FileFormatError(PipelineError):
    """Raised when a file has an invalid format."""

    def __init__(self, file_path: str, expected_format: str, stage: Optional[str] = None):
        """Initialize file format error."""
        message = f"Invalid file format for {file_path}. Expected: {expected_format}"
        super().__init__(message, stage, {"file": file_path, "expected_format": expected_format})


class DataValidationError(PipelineError):
    """Raised when data validation fails."""

    def __init__(self, message: str, field: str, stage: Optional[str] = None):
        """Initialize data validation error."""
        super().__init__(message, stage, {"field": field})


class StageExecutionError(PipelineError):
    """Raised when a stage fails to execute properly."""

    def __init__(
        self, stage_name: str, original_error: Exception, recovery_attempted: bool = False
    ):
        """Initialize stage execution error."""
        message = f"Stage '{stage_name}' failed: {str(original_error)}"
        super().__init__(
            message,
            stage_name,
            {
                "original_error": str(original_error),
                "error_type": type(original_error).__name__,
                "recovery_attempted": recovery_attempted,
            },
        )
        self.original_error = original_error
        self.recovery_attempted = recovery_attempted

    def __reduce__(self):
        """Custom pickling to handle multiprocessing correctly."""
        return (
            self.__class__,
            (self.stage, self.original_error, self.recovery_attempted),
            self.__dict__,
        )


def retry_on_failure(
    max_attempts: int = 3,
    delay: float = 1.0,
    backoff: float = 2.0,
    exceptions: tuple = (Exception,),
    logger: Optional[logging.Logger] = None,
) -> Callable:
    """Decorator to retry function on failure with exponential backoff.

    Parameters
    ----------
    max_attempts : int
        Maximum number of attempts
    delay : float
        Initial delay between attempts in seconds
    backoff : float
        Backoff multiplier for delay
    exceptions : tuple
        Tuple of exceptions to catch
    logger : logging.Logger, optional
        Logger for retry messages

    Returns
    -------
    Callable
        Decorated function
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            _logger = logger or logging.getLogger(func.__module__)
            current_delay = delay

            for attempt in range(max_attempts):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    if attempt == max_attempts - 1:
                        _logger.error(f"Failed after {max_attempts} attempts: {e}")
                        raise

                    _logger.warning(
                        f"Attempt {attempt + 1} failed: {e}. "
                        f"Retrying in {current_delay:.1f} seconds..."
                    )
                    time.sleep(current_delay)
                    current_delay *= backoff

        return wrapper

    return decorator


@contextmanager
def graceful_error_handling(
    stage_name: str,
    reraise: bool = True,
    fallback_value=None,
    logger: Optional[logging.Logger] = None,
):
    """Context manager for graceful error handling in stages.

    Parameters
    ----------
    stage_name : str
        Name of the stage for error reporting
    reraise : bool
        Whether to reraise exceptions after logging
    fallback_value : any
        Value to return on error if not reraising
    logger : logging.Logger, optional
        Logger instance

    Yields
    ------
    None

    Examples
    --------
    >>> with graceful_error_handling("my_stage"):
    ...     # Stage processing code
    ...     pass
    """
    _logger = logger or logging.getLogger(__name__)

    try:
        yield
    except PipelineError:
        # Pipeline errors are already formatted nicely
        raise
    except FileNotFoundError as e:
        _logger.error(f"File not found in {stage_name}: {e}")
        if reraise:
            raise PipelineError(f"Required file not found: {e}", stage=stage_name)
        return fallback_value
    except PermissionError as e:
        _logger.error(f"Permission denied in {stage_name}: {e}")
        if reraise:
            raise PipelineError(f"Permission denied: {e}", stage=stage_name)
        return fallback_value
    except Exception as e:
        _logger.error(f"Unexpected error in {stage_name}: {e}", exc_info=True)
        if reraise:
            raise StageExecutionError(stage_name, e)
        return fallback_value


class ErrorRecoveryStrategy:
    """Base class for error recovery strategies."""

    def can_recover(self, error: Exception, context: Dict) -> bool:
        """Check if recovery is possible for this error.

        Parameters
        ----------
        error : Exception
            The error that occurred
        context : dict
            Context information about the error

        Returns
        -------
        bool
            True if recovery is possible
        """
        return False

    def recover(self, error: Exception, context: Dict) -> any:
        """Attempt to recover from the error.

        Parameters
        ----------
        error : Exception
            The error that occurred
        context : dict
            Context information about the error

        Returns
        -------
        any
            Recovery result or None
        """
        return None


class FileRecoveryStrategy(ErrorRecoveryStrategy):
    """Recovery strategy for file-related errors."""

    def can_recover(self, error: Exception, context: Dict) -> bool:
        """Check if file error can be recovered."""
        if isinstance(error, FileNotFoundError):
            # Can try alternative paths or create file
            return "alternative_paths" in context or "can_create" in context
        if isinstance(error, PermissionError):
            # Can try alternative location
            return "alternative_location" in context
        return False

    def recover(self, error: Exception, context: Dict) -> any:
        """Attempt file error recovery."""
        if isinstance(error, FileNotFoundError):
            if "alternative_paths" in context:
                for alt_path in context["alternative_paths"]:
                    if Path(alt_path).exists():
                        logger.info(f"Using alternative path: {alt_path}")
                        return alt_path

            if context.get("can_create"):
                path = Path(context["original_path"])
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()
                logger.info(f"Created missing file: {path}")
                return str(path)

        if isinstance(error, PermissionError) and "alternative_location" in context:
            alt_loc = context["alternative_location"]
            logger.info(f"Using alternative location: {alt_loc}")
            return alt_loc

        return None


class ToolRecoveryStrategy(ErrorRecoveryStrategy):
    """Recovery strategy for external tool errors."""

    def can_recover(self, error: Exception, context: Dict) -> bool:
        """Check if tool error can be recovered."""
        if isinstance(error, ToolNotFoundError):
            return "alternative_tools" in context
        if hasattr(error, "returncode") and error.returncode != 0:
            return "retry_with_different_params" in context
        return False

    def recover(self, error: Exception, context: Dict) -> any:
        """Attempt tool error recovery."""
        if isinstance(error, ToolNotFoundError) and "alternative_tools" in context:
            for alt_tool in context["alternative_tools"]:
                # Check if alternative tool exists
                from shutil import which

                if which(alt_tool):
                    logger.info(f"Using alternative tool: {alt_tool}")
                    return alt_tool

        if hasattr(error, "returncode") and "retry_with_different_params" in context:
            # Return modified parameters for retry
            return context["retry_with_different_params"]

        return None


class ErrorRecoveryManager:
    """Manages error recovery strategies."""

    def __init__(self):
        """Initialize error recovery manager."""
        self.strategies: List[ErrorRecoveryStrategy] = [
            FileRecoveryStrategy(),
            ToolRecoveryStrategy(),
        ]

    def add_strategy(self, strategy: ErrorRecoveryStrategy) -> None:
        """Add a recovery strategy.

        Parameters
        ----------
        strategy : ErrorRecoveryStrategy
            Strategy to add
        """
        self.strategies.append(strategy)

    def attempt_recovery(self, error: Exception, context: Dict) -> Optional[any]:
        """Attempt to recover from an error.

        Parameters
        ----------
        error : Exception
            The error that occurred
        context : dict
            Context information about the error

        Returns
        -------
        any or None
            Recovery result if successful, None otherwise
        """
        for strategy in self.strategies:
            if strategy.can_recover(error, context):
                try:
                    result = strategy.recover(error, context)
                    if result is not None:
                        return result
                except Exception as e:
                    logger.warning(f"Recovery strategy failed: {e}")

        return None


# Global error recovery manager
error_recovery_manager = ErrorRecoveryManager()


def handle_stage_error(
    stage_name: str,
    error: Exception,
    context: Optional[Dict] = None,
    reraise: bool = True,
) -> Optional[any]:
    """Handle errors that occur during stage execution.

    Parameters
    ----------
    stage_name : str
        Name of the stage where error occurred
    error : Exception
        The error that occurred
    context : dict, optional
        Additional context for error handling
    reraise : bool
        Whether to reraise the error after handling

    Returns
    -------
    any or None
        Recovery result if successful, None otherwise

    Raises
    ------
    StageExecutionError
        If reraise is True and recovery fails
    """
    context = context or {}
    context["stage"] = stage_name

    # Log the error
    logger.error(f"Error in stage '{stage_name}': {error}")

    # Attempt recovery
    recovery_result = error_recovery_manager.attempt_recovery(error, context)

    if recovery_result is not None:
        logger.info(f"Successfully recovered from error in '{stage_name}'")
        return recovery_result

    if reraise:
        raise StageExecutionError(stage_name, error, recovery_attempted=True)

    return None


def validate_file_exists(file_path: Union[str, Path], stage_name: str) -> Path:
    """Validate that a file exists and is readable.

    Parameters
    ----------
    file_path : str or Path
        Path to validate
    stage_name : str
        Stage name for error reporting

    Returns
    -------
    Path
        Validated path object

    Raises
    ------
    FileNotFoundError
        If file doesn't exist
    PermissionError
        If file isn't readable
    """
    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")

    if not path.is_file():
        raise FileFormatError(str(path), "file", stage_name)

    # Check if readable
    try:
        with open(path, "r") as f:
            pass
    except PermissionError:
        raise PermissionError(f"Cannot read file: {path}")

    return path


def validate_output_directory(
    output_dir: Union[str, Path], stage_name: str, create: bool = True
) -> Path:
    """Validate output directory.

    Parameters
    ----------
    output_dir : str or Path
        Output directory path
    stage_name : str
        Stage name for error reporting
    create : bool
        Whether to create directory if it doesn't exist

    Returns
    -------
    Path
        Validated directory path

    Raises
    ------
    PermissionError
        If directory cannot be created or written to
    """
    path = Path(output_dir)

    if path.exists():
        if not path.is_dir():
            raise FileFormatError(str(path), "directory", stage_name)
    elif create:
        try:
            path.mkdir(parents=True, exist_ok=True)
        except PermissionError:
            raise PermissionError(f"Cannot create directory: {path}")
    else:
        raise FileNotFoundError(f"Output directory does not exist: {path}")

    # Check if writable
    test_file = path / ".write_test"
    try:
        test_file.touch()
        test_file.unlink()
    except PermissionError:
        raise PermissionError(f"Cannot write to directory: {path}")

    return path
