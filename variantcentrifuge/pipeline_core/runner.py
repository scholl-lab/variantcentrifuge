"""PipelineRunner - Executes stages in dependency order with parallel support."""

import gc
import logging
import multiprocessing
import time
from concurrent.futures import (
    Future,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    as_completed,
)
from typing import Any, Literal

try:
    import psutil
except ImportError:
    psutil = None  # type: ignore[assignment]

from .context import PipelineContext
from .runner_planner import (
    clear_stage_and_downstream,
    create_execution_plan,
    get_execution_plan,
    get_stages_to_execute_from,
    handle_selective_resume,
    log_runner_summary,
    ordered_stages,
    restore_completed_stage,
    restore_or_recompute_completed_stages,
    validate_resume_point,
)
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

                pipeline_version = context.config.get("pipeline_version", "pipeline")
                if context.checkpoint_state.can_resume(context.config, pipeline_version):
                    logger.info("Resume mode: Checking pipeline state...")
                    logger.info(context.checkpoint_state.get_summary())

                    ordered_stages = self._ordered_stages(stages)
                    context = self._restore_or_recompute_completed_stages(ordered_stages, context)

                    assert context.checkpoint_state is not None
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
        all_ordered = self._ordered_stages(stages)

        # Execute stages level by level
        for level, level_stages in enumerate(execution_plan):
            remaining_stages: list[Stage] = []
            for stage in level_stages:
                if self._should_skip_stage(stage, context):
                    if not context.is_complete(stage.name):
                        context.mark_complete(stage.name)
                        if hasattr(stage, "_handle_checkpoint_skip"):
                            context = stage._handle_checkpoint_skip(context)
                else:
                    remaining_stages.append(stage)
                    if self.enable_checkpoints and context.checkpoint_state:
                        self._clear_stage_and_downstream(
                            stage.name, all_ordered, context.checkpoint_state
                        )

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

        # Log peak memory usage
        if self._stage_metrics:
            mem_vals = [
                m["mem_after_mb"]
                for m in self._stage_metrics.values()
                if m.get("mem_after_mb") is not None
            ]
            if mem_vals:
                peak_mb = max(mem_vals)
                logger.info(f"Pipeline peak memory: {peak_mb:.0f} MB")

        return context

    def _handle_selective_resume(
        self, stages: list[Stage], context: PipelineContext
    ) -> list[Stage]:
        """Handle restart from a specific stage."""
        return handle_selective_resume(stages, context)

    def _should_skip_stage(self, stage: Stage, context: PipelineContext) -> bool:
        """Determine if a stage should be skipped due to checkpointing or completion."""
        if context.is_complete(stage.name):
            return True
        if not (self.enable_checkpoints and context.checkpoint_state):
            return False
        if not getattr(stage, "enable_checkpoint", True):
            return False
        return bool(context.checkpoint_state.should_skip_step(stage.name))

    def _ordered_stages(self, stages: list[Stage]) -> list[Stage]:
        """Return stages flattened in dependency execution order."""
        return ordered_stages(stages)

    def _clear_stage_and_downstream(
        self, stage_name: str, ordered: list[Stage], checkpoint_state: Any
    ) -> None:
        """Clear persisted completion for a stage and every downstream ordered stage."""
        clear_stage_and_downstream(stage_name, ordered, checkpoint_state)

    def _restore_completed_stage(self, stage: Stage, context: PipelineContext) -> PipelineContext:
        """Restore or skip a completed checkpoint stage according to its policy."""
        return restore_completed_stage(stage, context)

    def _restore_or_recompute_completed_stages(
        self, ordered: list[Stage], context: PipelineContext
    ) -> PipelineContext:
        """Restore safe completed stages, then clear volatile state for recomputation."""
        return restore_or_recompute_completed_stages(ordered, context)

    def _get_stages_to_execute_from(
        self, resume_from: str, stages: list[Stage], context: PipelineContext
    ) -> list[Stage]:
        """Get the filtered list of stages to execute when resuming from a specific stage."""
        return get_stages_to_execute_from(resume_from, stages, context)

    def get_execution_plan(self, stages: list[Stage], resume_from: str | None = None) -> list[str]:
        """Get the planned execution order, optionally starting from a specific stage."""
        return get_execution_plan(stages, resume_from=resume_from)

    def validate_resume_point(self, stage_name: str, stages: list[Stage]) -> tuple:
        """Validate that a resume point is valid and return any issues."""
        return validate_resume_point(stage_name, stages)

    def _create_execution_plan(self, stages: list[Stage]) -> list[list[Stage]]:
        """Create execution plan with stages grouped by dependency levels."""
        return create_execution_plan(stages)

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

    def _get_process_memory_mb(self) -> float | None:
        """Return current process RSS in MB, or None if psutil is unavailable."""
        if psutil is None:
            return None
        try:
            return float(psutil.Process().memory_info().rss / (1024 * 1024))
        except Exception as e:
            logger.debug(f"psutil memory tracking failed: {e}")
            return None

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

        mem_before_mb = self._get_process_memory_mb()
        if mem_before_mb is not None and logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"[{stage.name}] Starting - Memory: {mem_before_mb:.1f} MB")

        result = stage(context)
        elapsed = time.time() - start_time
        self._execution_times[stage.name] = elapsed

        # Capture subtask times if any were recorded
        if hasattr(stage, "subtask_times") and stage.subtask_times:
            self._subtask_times[stage.name] = stage.subtask_times

        gc.collect()

        mem_after_mb = self._get_process_memory_mb()
        if mem_before_mb is not None and mem_after_mb is not None:
            mem_delta_mb = mem_after_mb - mem_before_mb
            self._stage_metrics[stage.name] = {
                "mem_before_mb": mem_before_mb,
                "mem_after_mb": mem_after_mb,
                "mem_delta_mb": mem_delta_mb,
                "mem_peak_mb": max(mem_before_mb, mem_after_mb),
            }
            if logger.isEnabledFor(logging.DEBUG):
                freed_mb = mem_before_mb - mem_after_mb
                logger.debug(
                    f"[{stage.name}] Complete in {elapsed:.1f}s - "
                    f"Memory: {mem_after_mb:.1f} MB (freed {freed_mb:.1f} MB)"
                )

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
            updated_contexts = []
            for completed_count, future in enumerate(as_completed(future_to_stage), 1):
                stage = future_to_stage[future]
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
        for _i, updated_context in enumerate(updated_contexts):
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
        # NOTE: Intelligent batching deferred — current fixed-size batching is sufficient
        # for observed workloads. Revisit if stage count exceeds ~100.
        return stages

    def _get_executor_for_stage(self, stage: Stage) -> type:
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
        log_runner_summary(
            self._execution_times,
            self._subtask_times,
            self._stage_metrics,
            log=logger,
        )

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
