"""
Display Utilities - Enhanced status and information display for VariantCentrifuge.

This module provides comprehensive display functions for checkpoint status,
stage information, execution timelines, and resume recommendations.
"""

import logging
from datetime import datetime
from typing import Any

from .checkpoint import PipelineState, StepInfo
from .pipeline_core.stage import Stage
from .stages.stage_registry import StageInfo, get_registry

logger = logging.getLogger(__name__)


def display_enhanced_status(pipeline_state: PipelineState, available_stages: list[str]) -> None:
    """Display comprehensive checkpoint status with resume suggestions.

    Parameters
    ----------
    pipeline_state : PipelineState
        Current pipeline state with checkpoint information
    available_stages : List[str]
        List of available stage names for current configuration
    """
    status = pipeline_state.get_detailed_status()

    print("\n" + "=" * 80)
    print("VARIANTCENTRIFUGE CHECKPOINT STATUS")
    print("=" * 80)

    if not status.get("has_checkpoint"):
        print("❌ No checkpoint file found.")
        print("   Start a new pipeline run with --enable-checkpoint to create checkpoints.")
        return

    # Basic status information
    print(f"📄 Checkpoint file: {pipeline_state.state_file}")
    print(f"🏭 Pipeline version: {status.get('pipeline_version', 'unknown')}")

    if status.get("start_time"):
        start_time = datetime.fromtimestamp(status["start_time"])
        print(f"⏰ Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    completed_stages = pipeline_state.get_completed_stages()
    total_stages = len(status.get("stages", []))

    print(f"📊 Progress: {len(completed_stages)}/{total_stages} stages completed")

    if status.get("total_runtime", 0) > 0:
        print(f"⏱️  Total runtime: {status['total_runtime']:.1f}s")

    if status.get("last_completed_time"):
        last_time = datetime.fromtimestamp(status["last_completed_time"])
        print(f"🏁 Last completed: {last_time.strftime('%Y-%m-%d %H:%M:%S')}")

    # Display execution timeline
    if completed_stages:
        print("\n📈 Execution Timeline:")
        display_stage_timeline(completed_stages)

    # Display resume suggestions
    suggestions = status.get("resume_suggestions", [])
    if suggestions:
        print("\n💡 Resume Suggestions:")
        for suggestion in suggestions:
            print(f"   • {suggestion['stage']} - {suggestion['reason']}")
            print(f"     Command: variantcentrifuge {suggestion['command']} ...")

    # Display performance summary
    if completed_stages:
        print("\n⚡ Performance Summary:")
        display_performance_summary(completed_stages)

    # Display available resume points
    resume_points = [stage_name for stage_name, _ in completed_stages]
    if resume_points:
        print("\n🎯 Available Resume Points:")
        print(f"   {', '.join(resume_points)}")
        print("\n💻 Usage Examples:")
        print(f"   variantcentrifuge --resume-from {resume_points[-1]} ...")
        if len(resume_points) > 1:
            mid_point = len(resume_points) // 2
            print(f"   variantcentrifuge --resume-from {resume_points[mid_point]} ...")
        print("   variantcentrifuge --interactive-resume ...")

    print("\n" + "=" * 80)


def display_available_stages(stages: list[Stage], config: dict[str, Any]) -> None:
    """Display all available stages for current configuration.

    Parameters
    ----------
    stages : List[Stage]
        List of stage instances
    config : Dict[str, Any]
        Configuration dictionary
    """
    registry = get_registry()

    print("\n" + "=" * 80)
    print("AVAILABLE PIPELINE STAGES")
    print("=" * 80)

    print(f"📋 Configuration: {len(stages)} stages active")

    # Group stages by category
    categories: dict[str, list[tuple[Stage, StageInfo | None]]] = {}
    for stage in stages:
        stage_info = registry.get_stage_info(stage.name)
        category = stage_info.category if stage_info else "unknown"

        if category not in categories:
            categories[category] = []
        categories[category].append((stage, stage_info))

    # Display each category
    for category, stage_list in categories.items():
        print(f"\n🏷️  {category.upper()} ({len(stage_list)} stages):")
        print("-" * 60)

        for stage, stage_info in sorted(stage_list, key=lambda x: x[0].name):
            # Stage name and aliases
            name = stage.name
            aliases = stage_info.aliases if stage_info else []
            aliases_str = f" (aliases: {', '.join(aliases)})" if aliases else ""

            # Description
            description = stage_info.description if stage_info else stage.description
            description = description or "No description available"

            # Dependencies
            deps = stage.dependencies
            soft_deps: set[str] = getattr(stage, "soft_dependencies", set())

            print(f"   📌 {name}{aliases_str}")
            print(f"      {description}")

            if deps:
                print(f"      Dependencies: {', '.join(sorted(deps))}")
            if soft_deps:
                print(f"      Soft deps: {', '.join(sorted(soft_deps))}")

            if stage_info and stage_info.estimated_runtime > 0:
                print(f"      Est. runtime: {stage_info.estimated_runtime:.1f}s")

            print()

    print("=" * 80)


def display_stage_timeline(completed_stages: list[tuple[str, StepInfo]]) -> None:
    """Display execution timeline with durations.

    Parameters
    ----------
    completed_stages : List[Tuple[str, StepInfo]]
        List of completed stages with their step information
    """
    if not completed_stages:
        print("   No completed stages to display.")
        return

    print("-" * 70)
    print(f"{'#':<3} {'Stage':<25} {'Duration':<12} {'Start Time':<15} {'Status'}")
    print("-" * 70)

    total_duration: float = 0
    for i, (stage_name, step_info) in enumerate(completed_stages, 1):
        # Duration
        duration = step_info.duration or 0
        total_duration += duration
        duration_str = f"{duration:.1f}s" if duration > 0 else "N/A"

        # Start time
        start_time_str = ""
        if step_info.start_time:
            start_time = datetime.fromtimestamp(step_info.start_time)
            start_time_str = start_time.strftime("%H:%M:%S")

        # Status indicator
        status_icon = {"completed": "✅", "failed": "❌", "running": "⏳"}.get(
            step_info.status, "❓"
        )

        print(f"{i:<3} {stage_name:<25} {duration_str:<12} {start_time_str:<15} {status_icon}")

        # Show error if failed
        if step_info.status == "failed" and step_info.error:
            print(f"    ❌ Error: {step_info.error}")

    print("-" * 70)
    print(f"{'Total':<28} {total_duration:.1f}s")


def display_performance_summary(completed_stages: list[tuple[str, StepInfo]]) -> None:
    """Display performance summary for completed stages.

    Parameters
    ----------
    completed_stages : List[Tuple[str, StepInfo]]
        List of completed stages with their step information
    """
    # Calculate statistics
    durations = [step_info.duration for _, step_info in completed_stages if step_info.duration]

    if not durations:
        print("   No timing data available.")
        return

    total_time = sum(durations)
    avg_time = total_time / len(durations)

    # Sort by duration for top performers
    stages_by_duration = sorted(
        [(name, info) for name, info in completed_stages if info.duration],
        key=lambda x: x[1].duration or 0.0,
        reverse=True,
    )

    print(f"   Total execution time: {total_time:.1f}s")
    print(f"   Average stage time: {avg_time:.1f}s")
    print(
        f"   Fastest stage: {stages_by_duration[-1][0]} ({stages_by_duration[-1][1].duration:.1f}s)"
    )
    print(
        f"   Slowest stage: {stages_by_duration[0][0]} ({stages_by_duration[0][1].duration:.1f}s)"
    )

    # Show top 3 time consumers
    print("\n   🐌 Top Time Consumers:")
    for i, (stage_name, step_info) in enumerate(stages_by_duration[:3], 1):
        dur = step_info.duration or 0.0
        percentage = (dur / total_time) * 100
        print(f"      {i}. {stage_name}: {dur:.1f}s ({percentage:.1f}%)")
