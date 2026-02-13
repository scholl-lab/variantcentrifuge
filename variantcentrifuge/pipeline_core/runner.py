"""
PipelineRunner - Executes stages in dependency order with parallel support.

This module provides the PipelineRunner class that orchestrates stage execution,
handling dependencies, parallel execution, and error recovery.
"""

import logging
import multiprocessing
import time
from collections import defaultdict, deque
from concurrent.futures import (
    Future,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    as_completed,
)
from typing import Literal

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

    def __init__(
        self,
        enable_checkpoints: bool = False,
        max_workers: int | None = None,
        executor_type: Literal["thread", "process"] = "thread",
        enable_stage_batching: bool = True,
    ):
        """Initialize the pipeline runner.

        Parameters
        ----------
        enable_checkpoints : bool
            Enable checkpoint tracking for resume capability
        max_workers : int, optional
            Maximum parallel workers (default: CPU count)
        executor_type : {"thread", "process"}, optional
            Type of executor for parallel execution (default: "thread")
        enable_stage_batching : bool
            Enable batching of lightweight stages (default: True)
        """
        self.enable_checkpoints = enable_checkpoints
        self.max_workers = max_workers or multiprocessing.cpu_count()
        self.executor_type = executor_type
        self.enable_stage_batching = enable_stage_batching
        self._execution_times: dict[str, float] = {}
        self._stage_metrics: dict[str, dict[str, float]] = {}
        self._subtask_times: dict[str, dict[str, float]] = {}

    def run(self, stages: list[Stage], context: PipelineContext) -> PipelineContext:
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

        # Handle resume logic if checkpoint is enabled
        if self.enable_checkpoints and context.checkpoint_state:
            # Handle selective resume (--resume-from)
            if context.config.get("resume_from"):
                stages = self._handle_selective_resume(stages, context)
            # Handle traditional resume (--resume)
            elif context.config.get("resume", False):
                # Checkpoint state should already be loaded during initialization
                # Clean up any stale running stages first
                context.checkpoint_state.cleanup_stale_stages()

                pipeline_version = context.config.get("pipeline_version", "refactored_pipeline")
                if context.checkpoint_state.can_resume(context.config, pipeline_version):
                    logger.info("Resume mode: Checking pipeline state...")
                    logger.info(context.checkpoint_state.get_summary())

                    # Mark completed stages as complete in context
                    for stage in stages:
                        if context.checkpoint_state.should_skip_step(stage.name):
                            logger.info(f"Skipping completed stage: {stage.name}")
                            context.mark_complete(stage.name)

                            # Handle checkpoint skip logic for composite stages
                            if hasattr(stage, "_handle_checkpoint_skip"):
                                context = stage._handle_checkpoint_skip(context)

                    resume_point = context.checkpoint_state.get_resume_point()
                    if resume_point:
                        logger.info(f"Resuming from step: {resume_point}")
                else:
                    logger.warning("Cannot resume - configuration or version mismatch")
                    # Start fresh but don't clear existing state

        # Analyze dependencies and get execution plan (after resume logic)
        execution_plan = self._create_execution_plan(stages)

        logger.info(f"Execution plan has {len(execution_plan)} levels")
        for level, level_stages in enumerate(execution_plan):
            stage_names = [s.name for s in level_stages]
            logger.debug(f"Level {level}: {stage_names}")

        # Execute stages level by level
        for level, level_stages in enumerate(execution_plan):
            # Filter out already completed stages
            remaining_stages = [
                stage for stage in level_stages if not context.is_complete(stage.name)
            ]

            if remaining_stages:
                context = self._execute_level(remaining_stages, context, level)
            else:
                logger.info(f"Level {level}: All stages already complete, skipping")

            # Save checkpoint if enabled
            if self.enable_checkpoints and context.checkpoint_state:
                context.checkpoint_state.save()

        total_time = time.time() - start_time
        logger.info(f"Pipeline execution completed in {total_time:.1f}s")

        # Log execution time summary
        self._log_execution_summary()

        return context

    def _handle_selective_resume(
        self, stages: list[Stage], context: PipelineContext
    ) -> list[Stage]:
        """Handle restart from a specific stage.

        This method implements restart behavior: the specified stage and all
        subsequent stages will be re-executed, ignoring their previous completion
        status. Only stages that come before the restart point remain marked as complete.

        Parameters
        ----------
        stages : List[Stage]
            Original list of stages
        context : PipelineContext
            Pipeline context with configuration

        Returns
        -------
        List[Stage]
            Filtered list of stages to execute (from restart point onwards)
        """
        resume_from = context.config.get("resume_from")
        if not resume_from:
            return stages

        # Load checkpoint state
        if not context.checkpoint_state.load():
            raise ValueError("Cannot resume: No checkpoint file found")

        # Validate resume point
        stage_names = [stage.name for stage in stages]
        is_valid, error_msg = context.checkpoint_state.validate_resume_from_stage(
            resume_from, stage_names
        )

        if not is_valid:
            raise ValueError(f"Cannot resume from '{resume_from}': {error_msg}")

        # Validate dependencies and get stages to execute
        stages_to_execute = self._get_stages_to_execute_from(resume_from, stages, context)

        # Restart from specified stage - clear completion status for resume stage and all
        # subsequent stages
        completed_stages = context.checkpoint_state.get_available_resume_points()
        current_stage_names = {stage.name for stage in stages}
        stage_map = {stage.name: stage for stage in stages}

        # Get all stages in execution order to determine which come before/after resume point
        execution_plan = self._create_execution_plan(stages)
        all_stages_ordered = []
        for level in execution_plan:
            all_stages_ordered.extend(level)

        # Find the index of the resume stage
        resume_index = None
        for i, stage in enumerate(all_stages_ordered):
            if stage.name == resume_from:
                resume_index = i
                break

        if resume_index is None:
            raise ValueError(f"Could not find resume stage '{resume_from}' in execution plan")

        # Determine stages that come before the resume point (these can stay complete)
        stages_before_resume = {stage.name for stage in all_stages_ordered[:resume_index]}

        # Clear completion status in checkpoint for restart stage and all subsequent stages
        stages_to_clear = {stage.name for stage in all_stages_ordered[resume_index:]}
        for stage_name in stages_to_clear:
            context.checkpoint_state.clear_step_completion(stage_name)

        logger.info(f"Cleared completion status for {len(stages_to_clear)} stages to force restart")

        # Add prerequisite stages that should be considered complete for restart
        # These are fundamental setup stages that must have completed for pipeline to reach
        # this point
        prerequisite_stages = {
            "configuration_loading",
            "sample_config_loading",
            "gene_bed_creation",
            "pedigree_loading",
            "phenotype_loading",
            "scoring_config_loading",
            "annotation_config_loading",
        }

        # Mark prerequisite stages as complete if they exist in the pipeline and come before
        # resume point
        for prereq_stage in prerequisite_stages:
            if prereq_stage in stage_map and prereq_stage in stages_before_resume:
                logger.debug(f"Marking prerequisite stage as complete for restart: {prereq_stage}")
                context.mark_complete(prereq_stage)

                # Call checkpoint skip logic to restore context data (VCF samples, configs, etc.)
                stage_instance = stage_map[prereq_stage]
                if hasattr(stage_instance, "_handle_checkpoint_skip"):
                    logger.debug(
                        f"Calling checkpoint skip logic for prerequisite stage: {prereq_stage}"
                    )
                    context = stage_instance._handle_checkpoint_skip(context)

        # Mark only stages that come BEFORE the resume point as complete
        for completed_stage in completed_stages:
            if completed_stage in stages_before_resume:
                # Mark stage as complete in context
                context.mark_complete(completed_stage)

                # For stages that are in the original stages list (before filtering),
                # call their checkpoint skip logic to restore virtual dependencies
                if completed_stage in stage_map:
                    stage_instance = stage_map[completed_stage]
                    if hasattr(stage_instance, "_handle_checkpoint_skip"):
                        logger.debug(
                            f"Calling checkpoint skip logic for completed stage: {completed_stage}"
                        )
                        context = stage_instance._handle_checkpoint_skip(context)

                if completed_stage in current_stage_names:
                    logger.info(f"Marking completed stage as done: {completed_stage}")
                else:
                    logger.debug(f"Marking virtual completed stage as done: {completed_stage}")
            else:
                # Stage comes at or after resume point - do not mark as complete (force re-run)
                logger.debug(f"Stage '{completed_stage}' will be re-run (at or after resume point)")

        logger.info(
            f"Restart mode: Starting from '{resume_from}' with {len(stages_to_execute)} "
            f"stages to execute"
        )
        return stages_to_execute

    def _get_stages_to_execute_from(
        self, resume_from: str, stages: list[Stage], context: PipelineContext
    ) -> list[Stage]:
        """Get the filtered list of stages to execute when resuming from a specific stage.

        Parameters
        ----------
        resume_from : str
            Name of the stage to resume from
        stages : List[Stage]
            Original list of stages
        context : PipelineContext
            Pipeline context

        Returns
        -------
        List[Stage]
            Filtered stages to execute
        """
        # Build stage map
        stage_map = {stage.name: stage for stage in stages}

        # Find the resume stage
        if resume_from not in stage_map:
            raise ValueError(f"Resume stage '{resume_from}' not found in pipeline")

        # Get all stages in execution order
        execution_plan = self._create_execution_plan(stages)
        all_stages_ordered = []
        for level in execution_plan:
            all_stages_ordered.extend(level)

        # Find the index of the resume stage
        resume_index = None
        for i, stage in enumerate(all_stages_ordered):
            if stage.name == resume_from:
                resume_index = i
                break

        if resume_index is None:
            raise ValueError(f"Could not find resume stage '{resume_from}' in execution plan")

        # Return stages from resume point onwards
        return all_stages_ordered[resume_index:]

    def get_execution_plan(self, stages: list[Stage], resume_from: str | None = None) -> list[str]:
        """Get the planned execution order, optionally starting from a specific stage.

        Parameters
        ----------
        stages : List[Stage]
            List of stages
        resume_from : Optional[str]
            Stage name to resume from

        Returns
        -------
        List[str]
            Ordered list of stage names that will execute
        """
        execution_plan = self._create_execution_plan(stages)
        all_stages = []
        for level in execution_plan:
            all_stages.extend([stage.name for stage in level])

        if resume_from:
            if resume_from in all_stages:
                resume_index = all_stages.index(resume_from)
                return all_stages[resume_index:]
            else:
                return []

        return all_stages

    def validate_resume_point(self, stage_name: str, stages: list[Stage]) -> tuple:
        """Validate that a resume point is valid and return any issues.

        Parameters
        ----------
        stage_name : str
            Name of the stage to resume from
        stages : List[Stage]
            List of available stages

        Returns
        -------
        tuple
            (is_valid: bool, error_message: str)
        """
        stage_names = [stage.name for stage in stages]

        # Check if stage exists
        if stage_name not in stage_names:
            return False, f"Stage '{stage_name}' does not exist in current configuration"

        # Build dependency graph to check for issues
        try:
            _ = self._create_execution_plan(stages)

            # Check if stage has dependencies that might be skipped
            stage_map = {stage.name: stage for stage in stages}
            target_stage = stage_map[stage_name]

            # Get all stages that will execute from this point
            planned_stages = self.get_execution_plan(stages, stage_name)
            planned_set = set(planned_stages)

            # Check dependencies
            missing_deps = target_stage.dependencies - planned_set
            if missing_deps:
                return (
                    False,
                    f"Stage '{stage_name}' depends on {missing_deps} which won't be executed",
                )

            return True, ""

        except Exception as e:
            return False, f"Error validating resume point: {e!s}"

    def _create_execution_plan(self, stages: list[Stage]) -> list[list[Stage]]:
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
        stage_names = set(graph.keys())

        # Build effective dependencies including soft dependencies that exist
        dependencies = {}
        for stage in stages:
            # Start with hard dependencies
            deps = set(stage.dependencies)

            # Add soft dependencies if they exist in the pipeline
            if hasattr(stage, "soft_dependencies"):
                for soft_dep in stage.soft_dependencies:
                    if soft_dep in stage_names:
                        deps.add(soft_dep)

            dependencies[stage.name] = deps

        # Build reverse dependency graph (who depends on each stage)
        dependents = defaultdict(set)
        for stage_name, deps in dependencies.items():
            for dep in deps:
                if dep in graph:  # Only track dependencies that exist
                    dependents[dep].add(stage_name)

        # Calculate in-degree (number of dependencies each stage has)
        in_degree = {}
        for stage_name in graph:
            # Only count dependencies that exist in the graph
            valid_deps = [dep for dep in dependencies[stage_name] if dep in graph]
            in_degree[stage_name] = len(valid_deps)

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

            # Log level composition for debugging
            parallel_count = sum(1 for s in current_level if s.parallel_safe)
            logger.debug(
                f"Level {len(execution_plan)}: {len(current_level)} stages "
                f"({parallel_count} parallel-safe)"
            )

            execution_plan.append(current_level)

        # Check for circular dependencies
        if len(processed) != len(stages):
            unprocessed = set(graph.keys()) - processed
            logger.debug(f"Processed stages: {processed}")
            logger.debug(f"Unprocessed stages: {unprocessed}")
            logger.debug("Dependencies for unprocessed stages:")
            for stage_name in unprocessed:
                logger.debug(f"  {stage_name}: {dependencies[stage_name]}")
            raise ValueError(f"Circular dependency detected involving stages: {unprocessed}")

        return execution_plan

    def _execute_level(
        self, stages: list[Stage], context: PipelineContext, level: int
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

        # Capture subtask times if any were recorded
        if hasattr(stage, "subtask_times") and stage.subtask_times:
            self._subtask_times[stage.name] = stage.subtask_times

        return result

    def _execute_parallel_stages(
        self, stages: list[Stage], context: PipelineContext
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
        logger.info(
            f"Executing {len(stages)} stages in parallel using {self.executor_type} executor"
        )

        # Batch lightweight stages if enabled
        if self.enable_stage_batching:
            stages = self._batch_lightweight_stages(stages)

        # Choose executor based on type
        executor_class = (
            ProcessPoolExecutor if self.executor_type == "process" else ThreadPoolExecutor
        )

        # Adjust max_workers based on stage count
        effective_workers = min(self.max_workers, len(stages))

        with executor_class(max_workers=effective_workers) as executor:
            # Submit all stages
            future_to_stage: dict[Future, Stage] = {}
            for stage in stages:
                # Each parallel stage gets a copy of context to avoid conflicts
                # The context has thread-safe methods for updating shared state
                future = executor.submit(self._execute_stage, stage, context)
                future_to_stage[future] = stage

            # Wait for all to complete and collect results
            completed_count = 0
            updated_contexts = []
            for future in as_completed(future_to_stage):
                stage = future_to_stage[future]
                completed_count += 1
                try:
                    # Get the updated context from the stage execution
                    updated_context = future.result()
                    updated_contexts.append(updated_context)
                    logger.debug(
                        f"Parallel stage '{stage.name}' completed ({completed_count}/{len(stages)})"
                    )
                except Exception as e:
                    logger.error(f"Parallel stage '{stage.name}' failed: {e}")
                    # Cancel remaining futures
                    for f in future_to_stage:
                        if not f.done():
                            f.cancel()
                    raise

        # Merge all updates back into the main context
        for i, updated_context in enumerate(updated_contexts):
            # Log completed stages and important state from each context
            logger.debug(
                f"Merging context updates: {len(updated_context.completed_stages)} "
                f"completed stages, {len(updated_context.stage_results)} stage results"
            )
            context.merge_from(updated_context)

        logger.debug(f"Merged updates from {len(updated_contexts)} parallel stages")

        return context

    def _batch_lightweight_stages(self, stages: list[Stage]) -> list[Stage]:
        """Batch lightweight stages together for more efficient execution.

        Parameters
        ----------
        stages : List[Stage]
            Stages to potentially batch

        Returns
        -------
        List[Stage]
            Stages with lightweight ones potentially batched
        """
        # For now, return stages as-is
        # TODO: Implement intelligent batching based on estimated runtime
        return stages

    def _get_executor_for_stage(self, stage: Stage):
        """Get the appropriate executor type for a stage.

        Parameters
        ----------
        stage : Stage
            Stage to execute

        Returns
        -------
        type
            Executor class (ThreadPoolExecutor or ProcessPoolExecutor)
        """
        # Use ProcessPoolExecutor for CPU-intensive stages
        cpu_intensive_stages = {
            "variant_scoring",
            "inheritance_analysis",
            "gene_burden_analysis",
            "statistics_generation",
        }

        if stage.name in cpu_intensive_stages and self.executor_type == "process":
            return ProcessPoolExecutor
        return ThreadPoolExecutor

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

            # Log subtask times if available
            if stage_name in self._subtask_times:
                subtasks = self._subtask_times[stage_name]
                if subtasks:
                    # Sort subtasks by time
                    sorted_subtasks = sorted(subtasks.items(), key=lambda x: x[1], reverse=True)
                    for subtask_name, subtask_elapsed in sorted_subtasks:
                        subtask_percentage = (subtask_elapsed / elapsed) * 100 if elapsed > 0 else 0
                        # Add "(avg)" suffix for averaged chunk processing times
                        display_name = subtask_name
                        if stage_name == "parallel_complete_processing" and subtask_name in [
                            "variant_extraction",
                            "snpsift_filtering",
                            "field_extraction",
                        ]:
                            display_name += " (avg)"
                        logger.info(
                            f"  └─ {display_name:32s} {subtask_elapsed:6.1f}s "
                            f"({subtask_percentage:4.1f}%)"
                        )

        logger.info("-" * 60)
        logger.info(f"{'Total stage time:':30s} {total_time:6.1f}s")
        logger.info("=" * 60)

    def dry_run(self, stages: list[Stage]) -> list[list[str]]:
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
