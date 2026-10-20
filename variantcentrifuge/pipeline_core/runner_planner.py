"""Execution planning and dependency analysis for pipeline stages."""

from __future__ import annotations

import logging
from collections import defaultdict, deque
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .context import PipelineContext
    from .stage import Stage

logger = logging.getLogger(__name__)


def create_execution_plan(stages: list[Stage]) -> list[list[Stage]]:
    """Create execution plan with stages grouped by dependency levels.

    Stages at the same level can run in parallel.
    """
    graph = {stage.name: stage for stage in stages}
    stage_names = set(graph.keys())

    dependencies: dict[str, set[str]] = {}
    for stage in stages:
        deps = set(stage.dependencies)
        if hasattr(stage, "soft_dependencies"):
            for soft_dep in stage.soft_dependencies:
                if soft_dep in stage_names:
                    deps.add(soft_dep)
        dependencies[stage.name] = deps

    dependents: dict[str, set[str]] = defaultdict(set)
    for stage_name, deps in dependencies.items():
        for dep in deps:
            if dep in graph:
                dependents[dep].add(stage_name)

    in_degree: dict[str, int] = {}
    for stage_name in graph:
        valid_deps = [dep for dep in dependencies[stage_name] if dep in graph]
        in_degree[stage_name] = len(valid_deps)

    queue = deque([name for name, degree in in_degree.items() if degree == 0])
    execution_plan: list[list[Stage]] = []
    processed: set[str] = set()

    while queue:
        current_level: list[Stage] = []
        level_size = len(queue)

        for _ in range(level_size):
            stage_name = queue.popleft()
            current_level.append(graph[stage_name])
            processed.add(stage_name)

            for dependent in sorted(dependents[stage_name]):
                in_degree[dependent] -= 1
                if in_degree[dependent] == 0:
                    queue.append(dependent)

        # Name breaks runtime ties so the order never depends on PYTHONHASHSEED (#219).
        current_level.sort(key=lambda s: (-s.estimated_runtime, s.name))
        parallel_count = sum(1 for s in current_level if s.parallel_safe)
        logger.debug(
            f"Level {len(execution_plan)}: {len(current_level)} stages "
            f"({parallel_count} parallel-safe)"
        )
        execution_plan.append(current_level)

    if len(processed) != len(stages):
        unprocessed = sorted(set(graph.keys()) - processed)
        raise ValueError(f"Circular dependency detected involving stages: {unprocessed}")

    return execution_plan


def ordered_stages(stages: list[Stage]) -> list[Stage]:
    """Return stages flattened in dependency execution order."""
    ordered: list[Stage] = []
    for level in create_execution_plan(stages):
        ordered.extend(level)
    return ordered


def get_execution_plan(stages: list[Stage], resume_from: str | None = None) -> list[str]:
    """Get the planned execution order, optionally starting from a specific stage."""
    execution_plan = create_execution_plan(stages)
    all_stages = []
    for level in execution_plan:
        all_stages.extend([stage.name for stage in level])

    if resume_from:
        if resume_from in all_stages:
            resume_index = all_stages.index(resume_from)
            return all_stages[resume_index:]
        return []

    return all_stages


def validate_resume_point(stage_name: str, stages: list[Stage]) -> tuple[bool, str]:
    """Validate that a resume point is valid and return any issues."""
    stage_names = [stage.name for stage in stages]
    if stage_name not in stage_names:
        return False, f"Stage '{stage_name}' does not exist in current configuration"

    try:
        _ = create_execution_plan(stages)
        stage_map = {stage.name: stage for stage in stages}
        target_stage = stage_map[stage_name]

        planned_stages = get_execution_plan(stages, stage_name)
        planned_set = set(planned_stages)

        missing_deps = target_stage.dependencies - planned_set
        if missing_deps:
            return (
                False,
                f"Stage '{stage_name}' depends on {missing_deps} which won't be executed",
            )

        return True, ""
    except Exception as e:
        return False, f"Error validating resume point: {e!s}"


def clear_stage_and_downstream(
    stage_name: str, ordered: list[Stage], checkpoint_state: Any
) -> None:
    """Clear persisted completion for a stage and every downstream ordered stage."""
    names = [stage.name for stage in ordered]
    start = names.index(stage_name)
    for name in names[start:]:
        checkpoint_state.clear_step_completion(name, save=False)
    checkpoint_state.save()


def restore_completed_stage(stage: Stage, context: PipelineContext) -> PipelineContext:
    """Restore or skip a completed checkpoint stage according to its policy."""
    policy = stage.checkpoint_resume_policy
    if policy == "restore":
        if not hasattr(stage, "_handle_checkpoint_skip"):
            raise RuntimeError(
                f"Stage '{stage.name}' is restorable but has no checkpoint skip handler"
            )
        context = stage._handle_checkpoint_skip(context)
        context.mark_complete(stage.name)
    elif policy == "skip":
        context.mark_complete(stage.name)
    else:
        raise RuntimeError(f"Stage '{stage.name}' must be recomputed")
    return context


def restore_or_recompute_completed_stages(
    ordered: list[Stage], context: PipelineContext
) -> PipelineContext:
    """Restore safe completed stages, then clear volatile state for recomputation."""
    assert context.checkpoint_state is not None
    for stage in ordered:
        if not context.checkpoint_state.should_skip_step(stage.name):
            logger.info(
                "Recomputing checkpoint-completed stage '%s' and downstream stages",
                stage.name,
            )
            clear_stage_and_downstream(stage.name, ordered, context.checkpoint_state)
            break

        policy = stage.checkpoint_resume_policy
        if policy in {"restore", "skip"}:
            logger.info(f"Restoring completed stage from checkpoint: {stage.name}")
            context = restore_completed_stage(stage, context)
            continue

        logger.info(
            "Recomputing checkpoint-completed stage '%s' and downstream stages",
            stage.name,
        )
        clear_stage_and_downstream(stage.name, ordered, context.checkpoint_state)
        break
    return context


def handle_selective_resume(stages: list[Stage], context: PipelineContext) -> list[Stage]:
    """Handle restart from a specific stage."""
    assert context.checkpoint_state is not None, "checkpoint_state required for selective resume"
    checkpoint_state = context.checkpoint_state

    resume_from = context.config.get("resume_from")
    if not resume_from:
        return stages

    if not checkpoint_state.load():
        raise ValueError("Cannot resume: No checkpoint file found")

    stage_names = [stage.name for stage in stages]
    is_valid, error_msg = checkpoint_state.validate_resume_from_stage(resume_from, stage_names)

    if not is_valid:
        raise ValueError(f"Cannot resume from '{resume_from}': {error_msg}")

    all_stages_ordered = ordered_stages(stages)
    stage_names_ordered = [stage.name for stage in all_stages_ordered]
    resume_index = stage_names_ordered.index(resume_from)

    safe_restart = None
    for prior_stage in all_stages_ordered[:resume_index]:
        if prior_stage.checkpoint_resume_policy in {"restore", "skip"}:
            safe_restart = prior_stage.name
        if not checkpoint_state.should_skip_step(prior_stage.name):
            if context.is_complete(prior_stage.name):
                continue
            raise ValueError(
                f"Cannot resume from '{resume_from}' safely because prerequisite "
                f"stage '{prior_stage.name}' was not completed in the checkpoint. "
                f"Resume from '{prior_stage.name}' or an earlier stage."
            )
        if prior_stage.checkpoint_resume_policy == "recompute":
            suggestion = safe_restart or stage_names_ordered[0]
            raise ValueError(
                f"Cannot resume from '{resume_from}' safely because required "
                f"DataFrame state from '{prior_stage.name}' is not durable. "
                f"Resume from '{suggestion}' to restore durable inputs and "
                "re-run in-memory analysis stages."
            )
        context = restore_completed_stage(prior_stage, context)

    clear_stage_and_downstream(resume_from, all_stages_ordered, checkpoint_state)
    stages_to_execute = all_stages_ordered[resume_index:]

    logger.info(
        f"Restart mode: Starting from '{resume_from}' with {len(stages_to_execute)} "
        f"stages to execute"
    )
    return stages_to_execute


def get_stages_to_execute_from(
    resume_from: str, stages: list[Stage], context: PipelineContext
) -> list[Stage]:
    """Get the filtered list of stages to execute when resuming from a specific stage."""
    stage_map = {stage.name for stage in stages}
    if resume_from not in stage_map:
        raise ValueError(f"Resume stage '{resume_from}' not found in pipeline")

    all_stages_ordered = ordered_stages(stages)
    names = [s.name for s in all_stages_ordered]
    if resume_from not in names:
        raise ValueError(f"Could not find resume stage '{resume_from}' in execution plan")

    resume_index = names.index(resume_from)
    return all_stages_ordered[resume_index:]


def log_runner_summary(
    execution_times: dict[str, float],
    subtask_times: dict[str, dict[str, float]],
    stage_metrics: dict[str, dict[str, float]],
    log: logging.Logger | None = None,
) -> None:
    """Log summary of stage execution times and memory usage."""
    if not execution_times:
        return

    active_logger = log or logger
    active_logger.info("=" * 60)
    active_logger.info("Stage Execution Summary")
    active_logger.info("=" * 60)

    sorted_times = sorted(execution_times.items(), key=lambda x: x[1], reverse=True)
    total_time = sum(execution_times.values())

    for stage_name, elapsed in sorted_times:
        percentage = (elapsed / total_time) * 100 if total_time > 0 else 0
        active_logger.info(f"{stage_name:30s} {elapsed:6.1f}s ({percentage:4.1f}%)")

        subtasks = subtask_times.get(stage_name)
        if subtasks:
            sorted_subtasks = sorted(subtasks.items(), key=lambda x: x[1], reverse=True)
            for subtask_name, subtask_elapsed in sorted_subtasks:
                subtask_percentage = (subtask_elapsed / elapsed) * 100 if elapsed > 0 else 0
                display_name = subtask_name
                if stage_name == "parallel_complete_processing" and subtask_name in [
                    "variant_extraction",
                    "snpsift_filtering",
                    "field_extraction",
                ]:
                    display_name += " (avg)"
                active_logger.info(
                    f"  └─ {display_name:32s} {subtask_elapsed:6.1f}s ({subtask_percentage:4.1f}%)"
                )

    active_logger.info("-" * 60)
    active_logger.info(f"{'Total stage time:':30s} {total_time:6.1f}s")
    active_logger.info("=" * 60)

    valid_stage_metrics = {
        name: m
        for name, m in stage_metrics.items()
        if m.get("mem_before_mb") is not None
        and m.get("mem_after_mb") is not None
        and m.get("mem_delta_mb") is not None
    }
    if valid_stage_metrics:
        active_logger.info("")
        active_logger.info("=" * 60)
        active_logger.info("Memory Usage Summary")
        active_logger.info("=" * 60)

        sorted_metrics = sorted(
            valid_stage_metrics.items(),
            key=lambda x: abs(x[1]["mem_delta_mb"]),
            reverse=True,
        )

        active_logger.info(f"{'Stage':<30} {'Before':>10} {'After':>10} {'Delta':>10}")
        active_logger.info("-" * 60)

        for stage_name, metrics in sorted_metrics:
            before_mb = metrics["mem_before_mb"]
            after_mb = metrics["mem_after_mb"]
            delta_mb = metrics["mem_delta_mb"]
            delta_str = f"{delta_mb:+.1f} MB"
            active_logger.info(
                f"{stage_name:<30} {before_mb:>8.1f} MB {after_mb:>8.1f} MB {delta_str:>10}"
            )

        peak_mb = max(m["mem_after_mb"] for m in valid_stage_metrics.values())
        active_logger.info("-" * 60)
        active_logger.info(f"Peak RSS: {peak_mb:.1f} MB")
        active_logger.info("=" * 60)
