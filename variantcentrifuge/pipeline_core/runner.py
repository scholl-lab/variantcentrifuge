"""
PipelineRunner - Executes stages in dependency order with parallel support.

This module provides the PipelineRunner class that orchestrates stage execution,
handling dependencies, parallel execution, and error recovery.
"""

import logging
import time
from collections import defaultdict, deque
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional

from .context import PipelineContext
from .stage import Stage

logger = logging.getLogger(__name__)


class PipelineRunner:
    """Executes stages in dependency order with support for parallel execution.

    The runner analyzes stage dependencies, determines execution order, and
    runs stages in parallel where possible. It handles:
    - Dependency graph analysis
    - Topological sorting
    - Parallel execution of independent stages
    - Error handling and propagation
    - Checkpoint integration

    Attributes
    ----------
    enable_checkpoints : bool
        Whether to enable checkpoint tracking
    max_workers : int
        Maximum number of parallel workers
    """

    def __init__(self, enable_checkpoints: bool = False, max_workers: Optional[int] = None):
        """Initialize the pipeline runner.

        Parameters
        ----------
        enable_checkpoints : bool
            Enable checkpoint tracking for resume capability
        max_workers : int, optional
            Maximum parallel workers (default: CPU count)
        """
        self.enable_checkpoints = enable_checkpoints
        self.max_workers = max_workers
        self._execution_times: Dict[str, float] = {}

    def run(self, stages: List[Stage], context: PipelineContext) -> PipelineContext:
        """Execute all stages in dependency order with parallelization.

        Parameters
        ----------
        stages : List[Stage]
            List of stages to execute
        context : PipelineContext
            Initial pipeline context

        Returns
        -------
        PipelineContext
            Final context after all stages complete

        Raises
        ------
        ValueError
            If circular dependencies are detected
        Exception
            If any stage fails
        """
        start_time = time.time()
        logger.info(f"Starting pipeline execution with {len(stages)} stages")

        # Build stage map
        stage_map = {stage.name: stage for stage in stages}

        # Validate no duplicate stage names
        if len(stage_map) != len(stages):
            raise ValueError("Duplicate stage names detected")

        # Analyze dependencies and get execution plan
        execution_plan = self._create_execution_plan(stages)

        logger.info(f"Execution plan has {len(execution_plan)} levels")
        for level, level_stages in enumerate(execution_plan):
            stage_names = [s.name for s in level_stages]
            logger.debug(f"Level {level}: {stage_names}")

        # Execute stages level by level
        for level, level_stages in enumerate(execution_plan):
            context = self._execute_level(level_stages, context, level)

            # Save checkpoint if enabled
            if self.enable_checkpoints and context.checkpoint_state:
                context.checkpoint_state.save()

        total_time = time.time() - start_time
        logger.info(f"Pipeline execution completed in {total_time:.1f}s")

        # Log execution time summary
        self._log_execution_summary()

        return context

    def _create_execution_plan(self, stages: List[Stage]) -> List[List[Stage]]:
        """Create execution plan with stages grouped by dependency levels.

        Stages at the same level can run in parallel.

        Parameters
        ----------
        stages : List[Stage]
            All stages to execute

        Returns
        -------
        List[List[Stage]]
            Stages grouped by execution level

        Raises
        ------
        ValueError
            If circular dependencies are detected
        """
        # Build dependency graph
        graph = {stage.name: stage for stage in stages}
        dependencies = {stage.name: stage.dependencies for stage in stages}

        # Build reverse dependency graph (who depends on each stage)
        dependents = defaultdict(set)
        for stage_name, deps in dependencies.items():
            for dep in deps:
                if dep in graph:  # Only track dependencies that exist
                    dependents[dep].add(stage_name)

        # Calculate in-degree (number of dependencies each stage has)
        in_degree = {}
        for stage_name in graph:
            in_degree[stage_name] = len(dependencies[stage_name])

        # Find stages with no dependencies (in-degree 0)
        queue = deque([name for name, degree in in_degree.items() if degree == 0])

        # Group stages by level for parallel execution
        execution_plan = []
        processed = set()

        while queue:
            # Get all stages that can run at this level
            current_level = []
            level_size = len(queue)

            for _ in range(level_size):
                stage_name = queue.popleft()
                current_level.append(graph[stage_name])
                processed.add(stage_name)

                # Reduce in-degree for stages that depend on this one
                for dependent in dependents[stage_name]:
                    in_degree[dependent] -= 1
                    if in_degree[dependent] == 0:
                        queue.append(dependent)

            # Sort stages in level by estimated runtime (longest first)
            # This helps with better work distribution
            current_level.sort(key=lambda s: s.estimated_runtime, reverse=True)
            execution_plan.append(current_level)

        # Check for circular dependencies
        if len(processed) != len(stages):
            unprocessed = set(graph.keys()) - processed
            raise ValueError(f"Circular dependency detected involving stages: {unprocessed}")

        return execution_plan

    def _execute_level(
        self, stages: List[Stage], context: PipelineContext, level: int
    ) -> PipelineContext:
        """Execute all stages at a given dependency level.

        Parameters
        ----------
        stages : List[Stage]
            Stages to execute (can run in parallel)
        context : PipelineContext
            Current pipeline context
        level : int
            Dependency level (for logging)

        Returns
        -------
        PipelineContext
            Updated context after all stages complete
        """
        if not stages:
            return context

        # Single stage - no parallelization needed
        if len(stages) == 1:
            stage = stages[0]
            logger.info(f"Level {level}: Executing single stage '{stage.name}'")
            return self._execute_stage(stage, context)

        # Multiple stages - check if parallel execution is beneficial
        parallel_safe_stages = [s for s in stages if s.parallel_safe]
        sequential_stages = [s for s in stages if not s.parallel_safe]

        logger.info(
            f"Level {level}: {len(stages)} stages "
            f"({len(parallel_safe_stages)} parallel-safe, "
            f"{len(sequential_stages)} sequential)"
        )

        # Execute sequential stages first
        for stage in sequential_stages:
            context = self._execute_stage(stage, context)

        # Execute parallel-safe stages concurrently
        if parallel_safe_stages:
            context = self._execute_parallel_stages(parallel_safe_stages, context)

        return context

    def _execute_stage(self, stage: Stage, context: PipelineContext) -> PipelineContext:
        """Execute a single stage and track timing.

        Parameters
        ----------
        stage : Stage
            Stage to execute
        context : PipelineContext
            Current context

        Returns
        -------
        PipelineContext
            Updated context
        """
        start_time = time.time()
        result = stage(context)
        elapsed = time.time() - start_time
        self._execution_times[stage.name] = elapsed
        return result

    def _execute_parallel_stages(
        self, stages: List[Stage], context: PipelineContext
    ) -> PipelineContext:
        """Execute multiple stages in parallel.

        Parameters
        ----------
        stages : List[Stage]
            Parallel-safe stages to execute
        context : PipelineContext
            Current context

        Returns
        -------
        PipelineContext
            Updated context with all stage results merged
        """
        logger.info(f"Executing {len(stages)} stages in parallel")

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all stages
            future_to_stage: Dict[Future, Stage] = {}
            for stage in stages:
                # Each parallel stage gets a copy of context to avoid conflicts
                # The context has thread-safe methods for updating shared state
                future = executor.submit(self._execute_stage, stage, context)
                future_to_stage[future] = stage

            # Wait for all to complete
            for future in as_completed(future_to_stage):
                stage = future_to_stage[future]
                try:
                    # Get result (this will raise any exceptions)
                    future.result()
                    logger.debug(f"Parallel stage '{stage.name}' completed")
                except Exception as e:
                    logger.error(f"Parallel stage '{stage.name}' failed: {e}")
                    # Cancel remaining futures
                    for f in future_to_stage:
                        if not f.done():
                            f.cancel()
                    raise

        return context

    def _log_execution_summary(self) -> None:
        """Log summary of stage execution times."""
        if not self._execution_times:
            return

        logger.info("=" * 60)
        logger.info("Stage Execution Summary")
        logger.info("=" * 60)

        # Sort by execution time
        sorted_times = sorted(self._execution_times.items(), key=lambda x: x[1], reverse=True)

        total_time = sum(self._execution_times.values())

        for stage_name, elapsed in sorted_times:
            percentage = (elapsed / total_time) * 100 if total_time > 0 else 0
            logger.info(f"{stage_name:30s} {elapsed:6.1f}s ({percentage:4.1f}%)")

        logger.info("-" * 60)
        logger.info(f"{'Total stage time:':30s} {total_time:6.1f}s")
        logger.info("=" * 60)

    def dry_run(self, stages: List[Stage]) -> List[List[str]]:
        """Perform a dry run to show execution plan without running stages.

        Parameters
        ----------
        stages : List[Stage]
            Stages to analyze

        Returns
        -------
        List[List[str]]
            Stage names grouped by execution level
        """
        execution_plan = self._create_execution_plan(stages)
        return [[stage.name for stage in level] for level in execution_plan]
