"""
Stage - Abstract base class for all pipeline stages.

This module provides the unified Stage abstraction that all pipeline components
inherit from, ensuring consistent behavior and interface across the pipeline.
"""

import logging
import time
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Set

from .context import PipelineContext

logger = logging.getLogger(__name__)


class Stage(ABC):
    """Abstract base class for all pipeline stages - unified abstraction.

    All pipeline components (configuration loaders, processors, analyzers,
    reporters) inherit from this class. Each stage declares its dependencies
    and implements the _process method to perform its work.

    The stage execution is handled by __call__, which validates dependencies,
    logs execution, handles errors, and tracks timing.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Unique identifier for the stage.

        Returns
        -------
        str
            The stage name used for dependency tracking and logging
        """
        pass

    @property
    def dependencies(self) -> Set[str]:
        """Stage names that must complete before this stage.

        Returns
        -------
        Set[str]
            Set of stage names this stage depends on
        """
        return set()

    @property
    def soft_dependencies(self) -> Set[str]:
        """Stage names that should run before this stage if present.

        These are optional dependencies - if the stage is in the pipeline,
        it should run before this stage, but if it's not present, that's OK.

        Returns
        -------
        Set[str]
            Set of stage names this stage prefers to run after
        """
        return set()

    @property
    def description(self) -> str:
        """Human-readable description for logging.

        Returns
        -------
        str
            Description of what this stage does
        """
        return f"Stage: {self.name}"

    @property
    def parallel_safe(self) -> bool:
        """Whether this stage can run in parallel with other parallel-safe stages.

        Returns
        -------
        bool
            True if stage can run concurrently with other parallel-safe stages
        """
        return False

    @property
    def estimated_runtime(self) -> float:
        """Estimated runtime in seconds for scheduling optimization.

        Returns
        -------
        float
            Estimated seconds this stage takes to run
        """
        return 1.0

    def __call__(self, context: PipelineContext) -> PipelineContext:
        """Execute the stage with pre/post processing.

        This method handles:
        - Dependency validation
        - Execution logging
        - Error handling
        - Timing tracking
        - Checkpoint updates

        Parameters
        ----------
        context : PipelineContext
            The pipeline context

        Returns
        -------
        PipelineContext
            Updated context after stage execution

        Raises
        ------
        RuntimeError
            If dependencies are not satisfied
        Exception
            If stage execution fails
        """
        # Validate dependencies
        missing_deps = []
        for dep in self.dependencies:
            if not context.is_complete(dep):
                missing_deps.append(dep)

        if missing_deps:
            raise RuntimeError(
                f"Stage '{self.name}' requires these stages to complete first: "
                f"{', '.join(missing_deps)}"
            )

        # Check if already complete (for resume scenarios)
        if context.is_complete(self.name):
            logger.info(f"Stage '{self.name}' already complete, skipping")
            return context

        # Log execution start
        logger.info(f"Executing {self.description}")
        start_time = time.time()

        # Execute stage
        try:
            # Pre-execution hook
            self._pre_execute(context)

            # Main processing
            updated_context = self._process(context)

            # Post-execution hook
            self._post_execute(updated_context)

            # Mark complete
            elapsed = time.time() - start_time
            updated_context.mark_complete(self.name)

            logger.info(f"Stage '{self.name}' completed successfully in {elapsed:.1f}s")

            return updated_context

        except Exception as e:
            elapsed = time.time() - start_time
            logger.error(f"Stage '{self.name}' failed after {elapsed:.1f}s: {e}")
            raise

    @abstractmethod
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Core processing logic - must be implemented by subclasses.

        Parameters
        ----------
        context : PipelineContext
            The pipeline context

        Returns
        -------
        PipelineContext
            Updated context after processing
        """
        pass

    def _pre_execute(self, context: PipelineContext) -> None:
        """Execute hook called before stage execution.

        Override in subclasses to perform setup tasks.

        Parameters
        ----------
        context : PipelineContext
            The pipeline context
        """
        pass

    def _post_execute(self, context: PipelineContext) -> None:
        """Execute hook called after successful stage execution.

        Override in subclasses to perform cleanup tasks.

        Parameters
        ----------
        context : PipelineContext
            The pipeline context
        """
        pass

    def get_input_files(self, context: PipelineContext) -> List[Path]:
        """Return input files for checkpoint tracking.

        Override in subclasses to specify input files.

        Parameters
        ----------
        context : PipelineContext
            The pipeline context

        Returns
        -------
        List[Path]
            List of input file paths
        """
        return []

    def get_output_files(self, context: PipelineContext) -> List[Path]:
        """Return output files for checkpoint tracking.

        Override in subclasses to specify output files.

        Parameters
        ----------
        context : PipelineContext
            The pipeline context

        Returns
        -------
        List[Path]
            List of output file paths
        """
        return []

    def validate_prerequisites(self, context: PipelineContext) -> None:
        """Validate stage-specific prerequisites beyond dependencies.

        Override in subclasses to add custom validation.
        Raise exceptions if prerequisites are not met.

        Parameters
        ----------
        context : PipelineContext
            The pipeline context

        Raises
        ------
        ValueError
            If prerequisites are not satisfied
        """
        pass

    def __repr__(self) -> str:
        """Return string representation of the stage."""
        deps = f", depends_on={self.dependencies}" if self.dependencies else ""
        return f"{self.__class__.__name__}(name='{self.name}'{deps})"
