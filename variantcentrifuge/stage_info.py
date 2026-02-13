"""
Stage Information Display - Detailed stage metadata and visualization utilities.

This module provides functions for displaying stage information, dependency trees,
safety indicators, and other metadata useful for understanding pipeline structure.
"""

import logging
from typing import Any

from .checkpoint import PipelineState
from .pipeline_core.stage import Stage
from .stages.stage_registry import StageInfo, get_registry

logger = logging.getLogger(__name__)


def display_stage_details(stage_name: str, pipeline_state: PipelineState | None = None) -> None:
    """Display comprehensive information about a specific stage.

    Parameters
    ----------
    stage_name : str
        Name or alias of the stage
    pipeline_state : Optional[PipelineState]
        Pipeline state for execution history (optional)
    """
    registry = get_registry()

    # Resolve stage name and get info
    resolved_name = registry.resolve_stage_name(stage_name)
    if not resolved_name:
        print(f"❌ Stage '{stage_name}' not found.")
        _suggest_similar_stages(stage_name)
        return

    stage_info = registry.get_stage_info(resolved_name)
    if not stage_info:
        print(f"❌ No information available for stage '{resolved_name}'.")
        return

    print("\n" + "=" * 60)
    print(f"STAGE DETAILS: {resolved_name}")
    print("=" * 60)

    # Basic information
    print(f"📋 Category: {stage_info.category}")
    print(f"📝 Description: {stage_info.description}")

    if stage_info.aliases:
        print(f"🏷️  Aliases: {', '.join(stage_info.aliases)}")

    if stage_info.estimated_runtime > 0:
        print(f"⏱️  Estimated runtime: {stage_info.estimated_runtime:.1f}s")

    # Dependencies
    try:
        temp_stage = stage_info.class_ref()

        if temp_stage.dependencies:
            print(f"📌 Hard dependencies: {', '.join(sorted(temp_stage.dependencies))}")

        if hasattr(temp_stage, "soft_dependencies") and temp_stage.soft_dependencies:
            print(f"📌 Soft dependencies: {', '.join(sorted(temp_stage.soft_dependencies))}")

        if not temp_stage.dependencies and not getattr(temp_stage, "soft_dependencies", set()):
            print("📌 Dependencies: None (can run independently)")

        # Parallel safety
        if hasattr(temp_stage, "parallel_safe"):
            safety_icon = "✅" if temp_stage.parallel_safe else "⚠️"
            print(f"⚡ Parallel safe: {safety_icon} {temp_stage.parallel_safe}")

    except Exception as e:
        logger.debug(f"Could not instantiate stage {resolved_name}: {e}")
        print("📌 Dependencies: Unable to determine")

    # Execution history from checkpoint
    if pipeline_state:
        _display_stage_execution_history(resolved_name, pipeline_state)

    # Resume safety information
    if pipeline_state:
        _display_stage_resume_safety(resolved_name, pipeline_state)

    print("=" * 60)


def _display_stage_execution_history(stage_name: str, pipeline_state: PipelineState) -> None:
    """Display execution history for a stage from checkpoint data."""
    step_info = pipeline_state.state["steps"].get(stage_name)

    print("\n📊 Execution History:")

    if not step_info:
        print("   No execution history found")
        return

    print(f"   Status: {step_info.status}")

    if step_info.start_time:
        from datetime import datetime

        start_time = datetime.fromtimestamp(step_info.start_time)
        print(f"   Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    if step_info.end_time:
        from datetime import datetime

        end_time = datetime.fromtimestamp(step_info.end_time)
        print(f"   Completed: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")

    if step_info.duration:
        print(f"   Duration: {step_info.duration:.1f}s")

    if step_info.error:
        print(f"   Error: {step_info.error}")

    # Input/output files
    if step_info.input_files:
        print(f"   Input files: {len(step_info.input_files)} files")
        for file_info in step_info.input_files[:3]:  # Show first 3
            print(f"     • {file_info.path}")
        if len(step_info.input_files) > 3:
            print(f"     ... and {len(step_info.input_files) - 3} more")

    if step_info.output_files:
        print(f"   Output files: {len(step_info.output_files)} files")
        for file_info in step_info.output_files[:3]:  # Show first 3
            print(f"     • {file_info.path}")
        if len(step_info.output_files) > 3:
            print(f"     ... and {len(step_info.output_files) - 3} more")


def _display_stage_resume_safety(stage_name: str, pipeline_state: PipelineState) -> None:
    """Display resume safety information for a stage."""
    print("\n🔒 Resume Safety:")

    completed_stages = pipeline_state.get_available_resume_points()

    if stage_name in completed_stages:
        print("   ✅ Safe resume point (stage completed successfully)")
    else:
        print("   ❌ Not a valid resume point (stage not completed)")
        return

    # Check for potential issues
    registry = get_registry()
    stage_info = registry.get_stage_info(stage_name)

    if stage_info:
        if stage_info.category == "setup":
            print("   ⚠️  Setup stage - resuming here may skip configuration changes")
        elif stage_info.category == "processing":
            print("   ⚠️  Processing stage - ensure input data hasn't changed")
        elif stage_info.category == "analysis":
            print("   ✅ Analysis stage - generally safe for resume")
        elif stage_info.category == "output":
            print("   ✅ Output stage - safe for regenerating reports")


def _suggest_similar_stages(stage_name: str) -> None:
    """Suggest similar stage names when an exact match is not found."""
    registry = get_registry()
    all_stages = registry.get_all_stages()

    # Find stages that contain the search term
    partial_matches = []
    for name, stage_info in all_stages.items():
        if stage_name.lower() in name.lower():
            partial_matches.append(name)
        # Also check aliases
        for alias in stage_info.aliases:
            if stage_name.lower() in alias.lower():
                partial_matches.append(name)
                break

    if partial_matches:
        print("\n💡 Did you mean one of these?")
        for match in sorted(set(partial_matches))[:5]:  # Show top 5
            print(f"   • {match}")
    else:
        print("\n💡 Use 'variantcentrifuge --list-stages' to see all available stages.")


def display_dependency_tree(stage_name: str, all_stages: list[Stage], max_depth: int = 3) -> None:
    """Display dependency tree for a specific stage.

    Parameters
    ----------
    stage_name : str
        Name of the stage to analyze
    all_stages : List[Stage]
        List of all available stages
    max_depth : int
        Maximum depth to display in the tree
    """
    stage_map = {stage.name: stage for stage in all_stages}

    if stage_name not in stage_map:
        print(f"❌ Stage '{stage_name}' not found in current pipeline.")
        return

    print(f"\n🌳 Dependency Tree: {stage_name}")
    print("-" * 50)

    visited: set[str] = set()
    _display_dependencies_recursive(stage_name, stage_map, visited, 0, max_depth, "")


def _display_dependencies_recursive(
    stage_name: str,
    stage_map: dict[str, Stage],
    visited: set[str],
    depth: int,
    max_depth: int,
    prefix: str,
) -> None:
    """Recursively display stage dependencies."""
    if depth > max_depth or stage_name in visited:
        if stage_name in visited:
            print(f"{prefix}🔄 {stage_name} (circular reference)")
        else:
            print(f"{prefix}... (max depth reached)")
        return

    visited.add(stage_name)

    stage = stage_map.get(stage_name)
    if not stage:
        print(f"{prefix}❓ {stage_name} (not in current pipeline)")
        return

    # Stage icon based on status
    icon = "📌"
    print(f"{prefix}{icon} {stage_name}")

    # Get dependencies
    deps = stage.dependencies
    soft_deps: set[str] = getattr(stage, "soft_dependencies", set())

    # Filter to only existing stages
    active_deps = deps & set(stage_map.keys())
    active_soft_deps = soft_deps & set(stage_map.keys())

    all_deps = list(active_deps) + list(active_soft_deps)

    if not all_deps:
        return

    # Display dependencies
    for i, dep in enumerate(sorted(all_deps)):
        is_last = i == len(all_deps) - 1
        new_prefix = prefix + ("    " if is_last else "│   ")
        connector = "└── " if is_last else "├── "

        dep_type = "(soft)" if dep in active_soft_deps else ""
        print(f"{prefix}{connector}{dep} {dep_type}".strip())

        _display_dependencies_recursive(
            dep, stage_map, visited.copy(), depth + 1, max_depth, new_prefix
        )


def display_stage_categories_detailed() -> None:
    """Display detailed information about all stage categories."""
    registry = get_registry()
    categories = registry.get_categories()

    print("\n" + "=" * 80)
    print("STAGE CATEGORIES (DETAILED)")
    print("=" * 80)

    category_info = {
        "setup": {
            "icon": "🏗️",
            "description": "Configuration loading and pipeline initialization",
            "purpose": "Prepare the pipeline environment, load configurations, and validate inputs",
            "resume_safety": "⚠️ Caution - may skip configuration changes",
        },
        "processing": {
            "icon": "⚙️",
            "description": "Data extraction, filtering, and transformation",
            "purpose": "Extract variants from VCF files, apply filters, and transform data",
            "resume_safety": "⚠️ Caution - ensure input data hasn't changed",
        },
        "analysis": {
            "icon": "🔬",
            "description": "Variant analysis, scoring, and statistical computations",
            "purpose": "Perform inheritance analysis, variant scoring, and statistical tests",
            "resume_safety": "✅ Generally safe for resume",
        },
        "output": {
            "icon": "📄",
            "description": "Report generation and result formatting",
            "purpose": "Generate Excel, HTML, and IGV reports from analyzed data",
            "resume_safety": "✅ Safe for regenerating reports",
        },
    }

    for category in sorted(categories):
        stages = registry.get_stages_by_category(category)
        info = category_info.get(
            category,
            {
                "icon": "📦",
                "description": f"Stages in {category} category",
                "purpose": "No detailed information available",
                "resume_safety": "❓ Safety unknown",
            },
        )

        print(f"\n{info['icon']} {category.upper()}")
        print(f"   Description: {info['description']}")
        print(f"   Purpose: {info['purpose']}")
        print(f"   Resume safety: {info['resume_safety']}")
        print(f"   Stages ({len(stages)}):")

        for stage_name in sorted(stages.keys()):
            stage_info = stages[stage_name]
            runtime_str = (
                f" ({stage_info.estimated_runtime:.1f}s)"
                if stage_info.estimated_runtime > 0
                else ""
            )
            print(f"     • {stage_name}{runtime_str}")
            if stage_info.aliases:
                print(f"       Aliases: {', '.join(stage_info.aliases)}")

    print("\n" + "=" * 80)


def analyze_stage_impact(stage_name: str, all_stages: list[Stage]) -> dict[str, Any]:
    """Analyze the impact of running or skipping a specific stage.

    Parameters
    ----------
    stage_name : str
        Name of the stage to analyze
    all_stages : List[Stage]
        List of all available stages

    Returns
    -------
    Dict[str, Any]
        Analysis results including dependencies and dependents
    """
    stage_map = {stage.name: stage for stage in all_stages}

    if stage_name not in stage_map:
        return {"error": f"Stage '{stage_name}' not found"}

    target_stage = stage_map[stage_name]

    # Find what this stage depends on
    dependencies = target_stage.dependencies
    soft_dependencies: set[str] = getattr(target_stage, "soft_dependencies", set())

    # Find what depends on this stage
    dependents = set()
    soft_dependents = set()

    for stage in all_stages:
        if stage_name in stage.dependencies:
            dependents.add(stage.name)
        if hasattr(stage, "soft_dependencies") and stage_name in stage.soft_dependencies:
            soft_dependents.add(stage.name)

    return {
        "stage_name": stage_name,
        "dependencies": dependencies,
        "soft_dependencies": soft_dependencies,
        "dependents": dependents,
        "soft_dependents": soft_dependents,
        "can_run_independently": len(dependencies) == 0,
        "critical_path": len(dependents) > 0,
        "total_impact": len(dependents) + len(soft_dependents),
    }


def display_stage_impact_analysis(stage_name: str, all_stages: list[Stage]) -> None:
    """Display comprehensive impact analysis for a stage."""
    analysis = analyze_stage_impact(stage_name, all_stages)

    if "error" in analysis:
        print(f"❌ {analysis['error']}")
        return

    print(f"\n🎯 Impact Analysis: {stage_name}")
    print("-" * 50)

    # Dependencies
    if analysis["dependencies"]:
        print(f"📥 Requires: {', '.join(sorted(analysis['dependencies']))}")

    if analysis["soft_dependencies"]:
        print(f"📥 Prefers: {', '.join(sorted(analysis['soft_dependencies']))} (soft)")

    if analysis["can_run_independently"]:
        print("🆓 Can run independently: Yes")
    else:
        print("🔗 Can run independently: No")

    # Dependents
    if analysis["dependents"]:
        print(f"📤 Required by: {', '.join(sorted(analysis['dependents']))}")

    if analysis["soft_dependents"]:
        print(f"📤 Preferred by: {', '.join(sorted(analysis['soft_dependents']))} (soft)")

    # Impact summary
    total_impact = analysis["total_impact"]
    if total_impact == 0:
        print("💥 Impact: Minimal (no other stages depend on this)")
    elif total_impact <= 3:
        print(f"💥 Impact: Low ({total_impact} stages affected)")
    elif total_impact <= 7:
        print(f"💥 Impact: Medium ({total_impact} stages affected)")
    else:
        print(f"💥 Impact: High ({total_impact} stages affected)")

    if analysis["critical_path"]:
        print("⚠️  Critical: Yes (other stages depend on this)")
    else:
        print("✅ Critical: No (safe to skip/modify)")


def get_stage_safety_indicator(
    stage_name: str, category: str, pipeline_state: PipelineState | None = None
) -> str:
    """Get a safety indicator for resuming from a stage.

    Parameters
    ----------
    stage_name : str
        Name of the stage
    category : str
        Category of the stage
    pipeline_state : Optional[PipelineState]
        Pipeline state for validation

    Returns
    -------
    str
        Safety indicator emoji and text
    """
    if pipeline_state:
        completed_stages = pipeline_state.get_available_resume_points()
        if stage_name not in completed_stages:
            return "❌ Not completed"

    safety_map = {
        "setup": "⚠️ Caution",
        "processing": "⚠️ Caution",
        "analysis": "✅ Safe",
        "output": "✅ Safe",
    }

    return safety_map.get(category, "❓ Unknown")


def format_stage_summary(
    stage_name: str, stage_info: StageInfo, pipeline_state: PipelineState | None = None
) -> str:
    """Format a one-line summary of a stage.

    Parameters
    ----------
    stage_name : str
        Name of the stage
    stage_info : StageInfo
        Stage information
    pipeline_state : Optional[PipelineState]
        Pipeline state for execution info

    Returns
    -------
    str
        Formatted summary line
    """
    # Basic info
    category_icon = {"setup": "🏗️", "processing": "⚙️", "analysis": "🔬", "output": "📄"}.get(
        stage_info.category, "📦"
    )

    # Runtime info
    runtime_str = (
        f"{stage_info.estimated_runtime:.1f}s" if stage_info.estimated_runtime > 0 else "~"
    )

    # Execution status
    status_str = ""
    if pipeline_state:
        step_info = pipeline_state.state["steps"].get(stage_name)
        if step_info:
            if step_info.status == "completed":
                actual_runtime = step_info.duration or 0
                status_str = f" ✅({actual_runtime:.1f}s)"
            elif step_info.status == "failed":
                status_str = " ❌"
            elif step_info.status == "running":
                status_str = " ⏳"

    # Safety indicator
    safety = get_stage_safety_indicator(stage_name, stage_info.category, pipeline_state)
    safety_icon = safety.split()[0]  # Just the emoji

    return f"{category_icon} {stage_name:<20} {runtime_str:>6} {safety_icon}{status_str}"
