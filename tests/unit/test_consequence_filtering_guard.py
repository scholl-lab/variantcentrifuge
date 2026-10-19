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


def test_shipped_config_json_presets_dict_does_not_crash() -> None:
    """Verify that shipped config.json containing presets dict does not cause TypeError."""
    import json

    config_path = Path("config/config.json")
    if config_path.exists():
        with open(config_path) as f:
            config = json.load(f)
    else:
        config = {"presets": {"coding": "ANN[ANY].EFFECT has 'missense'", "rare": "AF < 0.01"}}

    # Enable burden with a benign filter
    config["perform_gene_burden"] = True
    config["filter"] = "QUAL > 30"
    config["snpeff_splitting_mode"] = None

    # Must NOT raise TypeError: unhashable type: 'dict'
    validate_annotation_consequence_filtering(config)


@pytest.mark.parametrize(
    "config_preset_field",
    [
        {"preset": "coding"},
        {"presets_list": ["high"]},
        {"preset": ["missense"]},
        {"presets": ["damaging"]},
    ],
)
def test_preset_variations_trigger_guard(config_preset_field: dict) -> None:
    """Verify that various ways to specify consequence presets trigger the guard."""
    config = {
        "perform_gene_burden": True,
        "snpeff_splitting_mode": None,
        **config_preset_field,
    }
    with pytest.raises(ValueError, match="--split-snpeff-lines before_filters"):
        validate_annotation_consequence_filtering(config)


def test_variant_extraction_stage_runs_consequence_guard(tmp_path: Path) -> None:
    """Verify VariantExtractionStage enforces consequence guard before extraction."""
    from variantcentrifuge.stages.processing_stages import VariantExtractionStage

    config = {
        "perform_gene_burden": True,
        "filter": "(ANN[ANY].IMPACT has 'HIGH')",
        "snpeff_splitting_mode": None,
        "vcf_file": "/tmp/dummy.vcf",
    }
    ctx = _make_context(config, tmp_path)
    stage = VariantExtractionStage()
    with pytest.raises(Exception, match="--split-snpeff-lines before_filters"):
        stage._process(ctx)
