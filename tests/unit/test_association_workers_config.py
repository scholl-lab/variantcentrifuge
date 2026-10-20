# File: tests/unit/test_association_workers_config.py
"""Unit tests for association_workers config plumbing (Phase 27, Plan 02).

Validates that AssociationConfig.association_workers has the correct default,
accepts custom values, and is correctly read by _build_assoc_config_from_context
from a PipelineContext.config dict.
"""

from unittest.mock import MagicMock

import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.stages.analysis_stages import _build_assoc_config_from_context


@pytest.mark.unit
def test_association_config_default_workers():
    """AssociationConfig.association_workers defaults to 1 (sequential)."""
    config = AssociationConfig()
    assert config.association_workers == 1


@pytest.mark.unit
def test_association_config_custom_workers():
    """AssociationConfig.association_workers accepts a custom positive value."""
    config = AssociationConfig(association_workers=4)
    assert config.association_workers == 4


@pytest.mark.unit
def test_association_config_auto_workers():
    """AssociationConfig.association_workers accepts -1 (auto/os.cpu_count()).

    Resolution from -1 to actual cpu_count happens in the engine at execution
    time, not in the config dataclass. The config merely stores the sentinel.
    """
    config = AssociationConfig(association_workers=-1)
    assert config.association_workers == -1


@pytest.mark.unit
def test_build_assoc_config_reads_workers():
    """_build_assoc_config_from_context reads association_workers from context.config.

    Verifies that when context.config contains association_workers=4, the
    resulting AssociationConfig has association_workers == 4.
    """
    context = MagicMock()
    context.config = {"association_workers": 4}

    result = _build_assoc_config_from_context(context)

    assert result.association_workers == 4


@pytest.mark.unit
@pytest.mark.parametrize(
    ("cfg", "expected"),
    [
        ({"association_workers": 0, "threads": 6}, 6),  # default: one worker per --threads
        ({"threads": 3}, 3),
        ({"association_workers": 2, "threads": 6}, 2),  # explicit --association-workers wins
    ],
)
def test_workers_default_from_threads(cfg, expected):
    """--threads N alone gives N association workers; an explicit worker count wins (#212)."""
    context = MagicMock()
    context.config = cfg
    assert _build_assoc_config_from_context(context).association_workers == expected
