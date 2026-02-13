"""
Utilities for stage dependency management.

This module provides helper functions to manage conditional dependencies
between stages, allowing stages to declare optional dependencies that
are only enforced when the dependent stages are actually in the pipeline.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from variantcentrifuge.pipeline_core import Stage


def filter_active_dependencies(
    declared_dependencies: set[str], active_stages: list["Stage"]
) -> set[str]:
    """Filter dependencies to only include stages that are active in the pipeline.

    This allows stages to declare dependencies on optional stages without
    causing errors when those stages aren't included.

    Parameters
    ----------
    declared_dependencies : Set[str]
        All dependencies declared by a stage
    active_stages : List[Stage]
        List of stages actually in the pipeline

    Returns
    -------
    Set[str]
        Filtered set of dependencies that are actually present
    """
    active_stage_names = {stage.name for stage in active_stages}
    return declared_dependencies & active_stage_names


def get_conditional_dependencies(
    required: set[str], optional: set[str], active_stages: list["Stage"]
) -> set[str]:
    """Get dependencies considering both required and optional dependencies.

    Parameters
    ----------
    required : Set[str]
        Dependencies that must be present (will error if missing)
    optional : Set[str]
        Dependencies that are only included if present in pipeline
    active_stages : List[Stage]
        List of stages actually in the pipeline

    Returns
    -------
    Set[str]
        Combined set of required and active optional dependencies
    """
    active_stage_names = {stage.name for stage in active_stages}

    # All required dependencies must be present
    missing_required = required - active_stage_names
    if missing_required:
        raise ValueError(f"Missing required dependencies: {missing_required}")

    # Only include optional dependencies that are present
    active_optional = optional & active_stage_names

    return required | active_optional


class ConditionalDependencyMixin:
    """Mixin to add conditional dependency support to stages.

    Stages can use this mixin to declare both required and optional
    dependencies, making it easier to handle variable pipeline configurations.
    """

    @property
    def required_dependencies(self) -> set[str]:
        """Dependencies that must always be present."""
        return set()

    @property
    def optional_dependencies(self) -> set[str]:
        """Dependencies that are only used if present in pipeline."""
        return set()

    def get_active_dependencies(self, active_stages: list["Stage"]) -> set[str]:
        """Get the actual dependencies based on active stages.

        Parameters
        ----------
        active_stages : List[Stage]
            List of stages in the pipeline

        Returns
        -------
        Set[str]
            Active dependencies for this stage
        """
        return get_conditional_dependencies(
            self.required_dependencies, self.optional_dependencies, active_stages
        )
