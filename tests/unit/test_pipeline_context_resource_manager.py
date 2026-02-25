"""
Tests for shared ResourceManager on PipelineContext (37-01 PERF-04).

Verifies that:
- PipelineContext accepts resource_manager=None (default)
- PipelineContext accepts a ResourceManager instance
- merge_from does NOT overwrite resource_manager on the parent
- Fallback pattern: when resource_manager is None, stages can create a local one
"""

from argparse import Namespace
from unittest.mock import MagicMock

import pytest

from variantcentrifuge.memory.resource_manager import ResourceManager
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace


@pytest.fixture
def workspace(tmp_path):
    """Create a minimal Workspace for test contexts."""
    return Workspace(output_dir=tmp_path, base_name="test")


@pytest.fixture
def minimal_context(workspace):
    """Create a minimal PipelineContext with defaults."""
    return PipelineContext(
        args=Namespace(),
        config={"threads": 1},
        workspace=workspace,
    )


@pytest.mark.unit
class TestPipelineContextResourceManager:
    """Tests for PipelineContext.resource_manager field."""

    def test_resource_manager_field_exists(self):
        """PipelineContext must have a resource_manager dataclass field."""
        assert "resource_manager" in PipelineContext.__dataclass_fields__

    def test_resource_manager_defaults_to_none(self, minimal_context):
        """resource_manager defaults to None when not explicitly set."""
        assert minimal_context.resource_manager is None

    def test_resource_manager_accepts_instance(self, workspace):
        """PipelineContext accepts a ResourceManager instance."""
        rm = ResourceManager(config={})
        context = PipelineContext(
            args=Namespace(),
            config={},
            workspace=workspace,
            resource_manager=rm,
        )
        assert context.resource_manager is rm

    def test_resource_manager_can_be_set_after_creation(self, minimal_context):
        """resource_manager can be assigned after context creation."""
        rm = ResourceManager(config={})
        minimal_context.resource_manager = rm
        assert minimal_context.resource_manager is rm

    def test_merge_from_does_not_overwrite_resource_manager(self, workspace):
        """merge_from must NOT overwrite the parent's resource_manager.

        Parallel contexts share the parent's ResourceManager; they must not
        replace it with their own (or with None).
        """
        rm_parent = ResourceManager(config={})
        parent = PipelineContext(
            args=Namespace(),
            config={"threads": 2},
            workspace=workspace,
            resource_manager=rm_parent,
        )

        # Child context has no resource_manager (as created in parallel execution)
        child = PipelineContext(
            args=Namespace(),
            config={"threads": 2},
            workspace=workspace,
        )
        child.mark_complete("child_stage")

        parent.merge_from(child)

        # Parent's resource_manager must still be the original instance
        assert parent.resource_manager is rm_parent

    def test_merge_from_does_not_overwrite_with_different_instance(self, workspace):
        """merge_from must not replace parent's rm even if child has one."""
        rm_parent = ResourceManager(config={})
        rm_child = ResourceManager(config={})
        parent = PipelineContext(
            args=Namespace(),
            config={},
            workspace=workspace,
            resource_manager=rm_parent,
        )
        child = PipelineContext(
            args=Namespace(),
            config={},
            workspace=workspace,
            resource_manager=rm_child,
        )

        parent.merge_from(child)

        # Parent retains its own resource_manager; child's is ignored
        assert parent.resource_manager is rm_parent

    def test_fallback_pattern_when_resource_manager_is_none(self, minimal_context):
        """Stages can create a local ResourceManager when context.resource_manager is None.

        This models the fallback pattern used in analysis_stages.py.
        """
        assert minimal_context.resource_manager is None

        # Simulate the fallback pattern used in stages
        rm = minimal_context.resource_manager
        if rm is None:
            rm = ResourceManager(config=minimal_context.config)

        # Fallback should produce a valid ResourceManager
        assert rm is not None
        assert isinstance(rm, ResourceManager)
        # And the context field is still None (not mutated by stages)
        assert minimal_context.resource_manager is None

    def test_resource_manager_with_mock(self, workspace):
        """Tests that use PipelineContext can mock resource_manager."""
        mock_rm = MagicMock(spec=ResourceManager)
        mock_rm.auto_chunk_size.return_value = 5000
        mock_rm.memory_gb = 8.0
        mock_rm.memory_safety_factor = 0.8

        context = PipelineContext(
            args=Namespace(),
            config={},
            workspace=workspace,
            resource_manager=mock_rm,
        )

        assert context.resource_manager is mock_rm
        # Mock behaves correctly
        assert context.resource_manager.auto_chunk_size(100000, 100) == 5000
