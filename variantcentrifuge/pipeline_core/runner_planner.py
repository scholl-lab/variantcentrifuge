"""Execution planning and dependency analysis for pipeline stages."""

from __future__ import annotations

import logging
from collections import defaultdict, deque
from typing import TYPE_CHECKING

if TYPE_CHECKING:
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

            for dependent in dependents[stage_name]:
                in_degree[dependent] -= 1
                if in_degree[dependent] == 0:
                    queue.append(dependent)

        current_level.sort(key=lambda s: s.estimated_runtime, reverse=True)
        parallel_count = sum(1 for s in current_level if s.parallel_safe)
        logger.debug(
            f"Level {len(execution_plan)}: {len(current_level)} stages "
            f"({parallel_count} parallel-safe)"
        )
        execution_plan.append(current_level)

    if len(processed) != len(stages):
        unprocessed = set(graph.keys()) - processed
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

    if stage_metrics:
        active_logger.info("")
        active_logger.info("=" * 60)
        active_logger.info("Memory Usage Summary")
        active_logger.info("=" * 60)

        sorted_metrics = sorted(
            stage_metrics.items(),
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

        peak_mb = max(m["mem_after_mb"] for m in stage_metrics.values())
        active_logger.info("-" * 60)
        active_logger.info(f"Peak RSS: {peak_mb:.1f} MB")
        active_logger.info("=" * 60)
