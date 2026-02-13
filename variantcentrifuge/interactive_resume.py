"""
Interactive Resume Interface - User-friendly stage selection for selective resume.

This module provides interactive command-line interfaces for selecting resume points,
displaying checkpoint status, and guiding users through safe resume operations.
"""

import logging
from datetime import datetime

from .checkpoint import PipelineState
from .stages.stage_registry import get_registry

logger = logging.getLogger(__name__)


def interactive_resume_selection(
    pipeline_state: PipelineState, available_stages: list[str]
) -> str | None:
    """Interactive CLI for selecting resume point.

    Parameters
    ----------
    pipeline_state : PipelineState
        Current pipeline state with checkpoint information
    available_stages : List[str]
        List of all available stage names in current configuration

    Returns
    -------
    Optional[str]
        Selected stage name to resume from, or None if cancelled
    """
    print("\n" + "=" * 70)
    print("INTERACTIVE RESUME POINT SELECTION")
    print("=" * 70)

    # Get detailed status information
    status = pipeline_state.get_detailed_status()

    if not status.get("has_checkpoint"):
        print("❌ No checkpoint file found. Cannot resume.")
        return None

    completed_stages = pipeline_state.get_completed_stages()

    if not completed_stages:
        print("❌ No completed stages found in checkpoint. Cannot resume.")
        return None

    print("\n📊 Pipeline Status:")
    print(f"   Total stages completed: {len(completed_stages)}")
    print(f"   Total runtime: {status.get('total_runtime', 0):.1f}s")

    if status.get("last_completed_time"):
        last_time = datetime.fromtimestamp(status["last_completed_time"])
        print(f"   Last completed: {last_time.strftime('%Y-%m-%d %H:%M:%S')}")

    # Display completed stages in table format
    print("\n📋 Completed Stages:")
    print("-" * 70)
    print(f"{'#':<3} {'Stage Name':<25} {'Duration':<10} {'Completed':<20}")
    print("-" * 70)

    for i, (stage_name, step_info) in enumerate(completed_stages, 1):
        duration = f"{step_info.duration:.1f}s" if step_info.duration else "N/A"
        completed_time = ""
        if step_info.end_time:
            completed_time = datetime.fromtimestamp(step_info.end_time).strftime("%H:%M:%S")

        print(f"{i:<3} {stage_name:<25} {duration:<10} {completed_time:<20}")

    print("-" * 70)

    # Show resume suggestions
    suggestions = status.get("resume_suggestions", [])
    if suggestions:
        print("\n💡 Resume Suggestions:")
        for i, suggestion in enumerate(suggestions, 1):
            print(f"   {i}. {suggestion['stage']} - {suggestion['reason']}")

    # Get available resume points
    resume_points = [stage_name for stage_name, _ in completed_stages]

    while True:
        print("\n🎯 Select Resume Point:")
        print("   Enter stage name, number from list above, or:")
        print("   'list' - show all available stages")
        print("   'info <stage>' - show stage details")
        print("   'quit' - cancel and exit")

        try:
            user_input = input("\nYour choice: ").strip()
        except (KeyboardInterrupt, EOFError):
            print("\n\n❌ Cancelled by user.")
            return None

        if not user_input:
            continue

        if user_input.lower() in ["quit", "q", "exit"]:
            print("❌ Cancelled by user.")
            return None

        if user_input.lower() == "list":
            _display_all_available_stages(available_stages)
            continue

        if user_input.lower().startswith("info "):
            stage_name = user_input[5:].strip()
            _display_stage_info(stage_name, available_stages, pipeline_state)
            continue

        # Try to resolve user input to a stage name
        selected_stage = _resolve_user_input(user_input, resume_points, completed_stages)

        if not selected_stage:
            print(f"❌ Invalid selection: '{user_input}'")
            print("   Please enter a valid stage name or number.")
            continue

        # Validate the resume point
        is_valid, error_msg = pipeline_state.validate_resume_from_stage(
            selected_stage, available_stages
        )

        if not is_valid:
            print(f"❌ Cannot resume from '{selected_stage}': {error_msg}")
            continue

        # Show confirmation with impact summary
        if _confirm_resume_selection(selected_stage, available_stages, pipeline_state):
            return selected_stage

        # User declined, continue the loop


def _resolve_user_input(
    user_input: str, resume_points: list[str], completed_stages: list[tuple[str, object]]
) -> str | None:
    """Resolve user input to a stage name.

    Parameters
    ----------
    user_input : str
        User's input (stage name or number)
    resume_points : List[str]
        List of available resume point names
    completed_stages : List[Tuple[str, object]]
        List of completed stages with their info

    Returns
    -------
    Optional[str]
        Resolved stage name or None if invalid
    """
    # Try as direct stage name first
    if user_input in resume_points:
        return user_input

    # Try as number (1-indexed)
    try:
        index = int(user_input) - 1
        if 0 <= index < len(completed_stages):
            return completed_stages[index][0]
    except ValueError:
        pass

    # Try partial matching
    matches = [stage for stage in resume_points if stage.lower().startswith(user_input.lower())]
    if len(matches) == 1:
        return matches[0]
    elif len(matches) > 1:
        print(f"❌ Ambiguous selection '{user_input}'. Matches: {', '.join(matches)}")
        return None

    # Try alias resolution through registry
    registry = get_registry()
    resolved = registry.resolve_stage_name(user_input)
    if resolved and resolved in resume_points:
        return resolved

    return None


def _display_all_available_stages(available_stages: list[str]) -> None:
    """Display all available stages in current configuration."""
    registry = get_registry()

    print(f"\n📚 All Available Stages ({len(available_stages)}):")
    print("-" * 50)

    # Group by category
    categories = {}
    for stage_name in available_stages:
        stage_info = registry.get_stage_info(stage_name)
        if stage_info:
            category = stage_info.category
            if category not in categories:
                categories[category] = []
            categories[category].append((stage_name, stage_info))

    for category, stages in categories.items():
        print(f"\n{category.upper()}:")
        for stage_name, stage_info in sorted(stages):
            aliases_str = (
                f" (aliases: {', '.join(stage_info.aliases)})" if stage_info.aliases else ""
            )
            print(f"  • {stage_name}{aliases_str}")


def _display_stage_info(
    stage_name: str, available_stages: list[str], pipeline_state: PipelineState
) -> None:
    """Display detailed information about a specific stage."""
    registry = get_registry()

    # Resolve stage name
    resolved_name = registry.resolve_stage_name(stage_name)
    if not resolved_name:
        print(f"❌ Stage '{stage_name}' not found.")
        return

    stage_info = registry.get_stage_info(resolved_name)
    if not stage_info:
        print(f"❌ No information available for stage '{resolved_name}'.")
        return

    print(f"\n📋 Stage Information: {resolved_name}")
    print("-" * 50)
    print(f"Category: {stage_info.category}")
    print(f"Description: {stage_info.description}")

    if stage_info.aliases:
        print(f"Aliases: {', '.join(stage_info.aliases)}")

    if stage_info.estimated_runtime > 0:
        print(f"Estimated runtime: {stage_info.estimated_runtime:.1f}s")

    # Check if stage was completed
    completed_stages = pipeline_state.get_available_resume_points()
    if resolved_name in completed_stages:
        print("✅ Status: Completed")

        # Get step info
        step_info = pipeline_state.state["steps"].get(resolved_name)
        if step_info and step_info.duration:
            print(f"Actual runtime: {step_info.duration:.1f}s")
    else:
        print("⏳ Status: Not completed")

    # Show dependencies if available
    if resolved_name in available_stages:
        try:
            # Create temporary stage instance to get dependencies
            stage_class = stage_info.class_ref
            temp_stage = stage_class()

            if temp_stage.dependencies:
                print(f"Dependencies: {', '.join(temp_stage.dependencies)}")
            if hasattr(temp_stage, "soft_dependencies") and temp_stage.soft_dependencies:
                print(f"Soft dependencies: {', '.join(temp_stage.soft_dependencies)}")

        except Exception:
            pass  # Skip dependency info if we can't instantiate


def _confirm_resume_selection(
    selected_stage: str, available_stages: list[str], pipeline_state: PipelineState
) -> bool:
    """Show confirmation dialog with impact summary."""
    print(f"\n🎯 Resume Point Selected: {selected_stage}")
    print("=" * 50)

    # Show what stages will be executed
    completed_stages = pipeline_state.get_available_resume_points()
    stages_to_execute = pipeline_state.get_stages_to_execute(selected_stage, available_stages)

    print("📊 Impact Summary:")
    print(f"   Stages to skip: {len(completed_stages)} completed stages")
    print(f"   Stages to execute: {len(stages_to_execute)} stages")

    if len(stages_to_execute) <= 5:
        print(f"   Execution plan: {' → '.join(stages_to_execute)}")
    else:
        print(
            f"   Execution plan: {' → '.join(stages_to_execute[:3])} ... → {stages_to_execute[-1]}"
        )

    # Calculate estimated time savings
    total_saved_time = 0
    for stage_name in completed_stages:
        step_info = pipeline_state.state["steps"].get(stage_name)
        if step_info and step_info.duration:
            total_saved_time += step_info.duration

    if total_saved_time > 0:
        print(f"   Estimated time saved: {total_saved_time:.1f}s")

    # Safety check
    registry = get_registry()
    stage_info = registry.get_stage_info(selected_stage)

    # Warn about potentially risky resume points
    if stage_info and stage_info.category in ["setup", "processing"]:
        print("\n⚠️  Safety Notice:")
        print(f"   Resuming from '{selected_stage}' will skip data processing stages.")
        print("   Ensure input data and configuration haven't changed.")

    # Get user confirmation
    while True:
        try:
            confirm = input(f"\n✅ Resume from '{selected_stage}'? [y/N]: ").strip().lower()
        except (KeyboardInterrupt, EOFError):
            print("\n❌ Cancelled by user.")
            return False

        if confirm in ["y", "yes"]:
            return True
        elif confirm in ["n", "no", ""]:
            return False
        else:
            print("Please enter 'y' for yes or 'n' for no.")


def handle_interactive_resume(
    args, pipeline_state: PipelineState, available_stages: list[str]
) -> str | None:
    """Handle the --interactive-resume CLI option.

    Parameters
    ----------
    args : argparse.Namespace
        CLI arguments
    pipeline_state : PipelineState
        Current pipeline state
    available_stages : List[str]
        List of available stage names

    Returns
    -------
    Optional[str]
        Selected stage name or None if cancelled
    """
    try:
        return interactive_resume_selection(pipeline_state, available_stages)
    except Exception as e:
        logger.error(f"Interactive resume failed: {e}")
        print(f"\n❌ Interactive resume failed: {e}")
        return None
