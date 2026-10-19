"""Unit and regression tests for annotation-consistent consequence filtering guard (#135)."""

from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.processing_stages import (
    SnpSiftFilterStage,
    validate_annotation_consequence_filtering,
)


def _make_context(config: dict, tmp_path: Path) -> PipelineContext:
    inter_dir = tmp_path / "intermediate"
    inter_dir.mkdir(parents=True, exist_ok=True)
    ws = Mock(spec=Workspace)
    ws.base_name = "test"
    ws.intermediate_dir = inter_dir
    ws.get_intermediate_path = lambda name: inter_dir / name
    ctx = PipelineContext(args=Namespace(), config=config, workspace=ws)
    return ctx


@pytest.mark.parametrize("split_mode", [None, "after_filters"])
@pytest.mark.parametrize("analysis_key", ["perform_gene_burden", "perform_association"])
def test_unsafe_burden_fails_without_before_filters(
    tmp_path: Path, split_mode: str | None, analysis_key: str
) -> None:
    """Verify that consequence filtering for burden/association fails if not before_filters."""
    config = {
        analysis_key: True,
        "filter": "(ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE')",
        "snpeff_splitting_mode": split_mode,
    }
    with pytest.raises(ValueError, match="--split-snpeff-lines before_filters"):
        validate_annotation_consequence_filtering(config)

    ctx = _make_context(config, tmp_path)
    stage = SnpSiftFilterStage()
    with pytest.raises(ValueError, match="--split-snpeff-lines before_filters"):
        stage._process(ctx)


def test_preset_triggers_consequence_guard(tmp_path: Path) -> None:
    """Verify that consequence presets like 'coding' trigger the validation guard."""
    config = {
        "perform_gene_burden": True,
        "presets": ["coding"],
        "snpeff_splitting_mode": None,
    }
    with pytest.raises(ValueError, match="--split-snpeff-lines before_filters"):
        validate_annotation_consequence_filtering(config)


def test_safe_configuration_passes_validation() -> None:
    """Verify that before_filters mode passes validation."""
    config = {
        "perform_gene_burden": True,
        "filter": "(ANN[ANY].IMPACT has 'HIGH')",
        "snpeff_splitting_mode": "before_filters",
    }
    # Should not raise
    validate_annotation_consequence_filtering(config)


def test_non_burden_run_passes_validation() -> None:
    """Verify that runs without burden/association pass validation even without splitting."""
    config = {
        "perform_gene_burden": False,
        "perform_association": False,
        "filter": "(ANN[ANY].IMPACT has 'HIGH')",
        "snpeff_splitting_mode": None,
    }
    # Should not raise
    validate_annotation_consequence_filtering(config)
