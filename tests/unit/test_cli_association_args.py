"""
Unit tests for CLI argument propagation to AssociationConfig.

Tests that the four new CLI arguments added in Phase 28 propagate correctly
through the cfg dict to _build_assoc_config_from_context() and AssociationConfig:

- --skat-method   -> cfg["skat_method"]           -> AssociationConfig.skat_method
- --min-cases     -> cfg["association_min_cases"]  -> AssociationConfig.min_cases
- --max-case-control-ratio -> cfg["association_max_case_control_ratio"]
                                                   -> AssociationConfig.max_case_control_ratio
- --min-case-carriers -> cfg["association_min_case_carriers"]
                                                   -> AssociationConfig.min_case_carriers

Coverage:
- CLI value from each of the four new args
- Default value when CLI arg not provided (no key in cfg)
- JSON association section values work when CLI args absent
- CLI keys override JSON association section (CLI wins)
"""

from __future__ import annotations

from argparse import Namespace
from unittest.mock import Mock

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import _build_assoc_config_from_context

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_context(config: dict) -> PipelineContext:
    """Build a minimal PipelineContext with the given config dict."""
    workspace = Mock(spec=Workspace)
    return PipelineContext(
        args=Namespace(),
        config=config,
        workspace=workspace,
    )


# ---------------------------------------------------------------------------
# Tests: --skat-method
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_skat_method_from_cli():
    """--skat-method SKATO propagates as cfg["skat_method"] = "SKATO"."""
    ctx = _make_context({"skat_method": "SKATO"})
    assoc_config = _build_assoc_config_from_context(ctx)
    assert assoc_config.skat_method == "SKATO"


@pytest.mark.unit
def test_skat_method_burden_from_cli():
    """--skat-method Burden propagates as cfg["skat_method"] = "Burden"."""
    ctx = _make_context({"skat_method": "Burden"})
    assoc_config = _build_assoc_config_from_context(ctx)
    assert assoc_config.skat_method == "Burden"


@pytest.mark.unit
def test_skat_method_default():
    """When skat_method not in cfg, AssociationConfig uses default 'SKAT'."""
    ctx = _make_context({})
    assoc_config = _build_assoc_config_from_context(ctx)
    assert assoc_config.skat_method == "SKAT"


# ---------------------------------------------------------------------------
# Tests: --min-cases
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_min_cases_from_cli():
    """--min-cases 100 propagates as cfg["association_min_cases"] = 100."""
    ctx = _make_context({"association_min_cases": 100})
    assoc_config = _build_assoc_config_from_context(ctx)
    assert assoc_config.min_cases == 100


@pytest.mark.unit
def test_min_cases_default():
    """When association_min_cases not in cfg, AssociationConfig uses default 200."""
    ctx = _make_context({})
    assoc_config = _build_assoc_config_from_context(ctx)
    assert assoc_config.min_cases == 200


# ---------------------------------------------------------------------------
# Tests: --max-case-control-ratio
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_max_case_control_ratio_from_cli():
    """--max-case-control-ratio 10 propagates as cfg["association_max_case_control_ratio"]."""
    ctx = _make_context({"association_max_case_control_ratio": 10.0})
    assoc_config = _build_assoc_config_from_context(ctx)
    assert assoc_config.max_case_control_ratio == 10.0


@pytest.mark.unit
def test_max_case_control_ratio_default():
    """When association_max_case_control_ratio not in cfg, default is 20.0."""
    ctx = _make_context({})
    assoc_config = _build_assoc_config_from_context(ctx)
    assert assoc_config.max_case_control_ratio == 20.0


# ---------------------------------------------------------------------------
# Tests: --min-case-carriers
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_min_case_carriers_from_cli():
    """--min-case-carriers 5 propagates as cfg["association_min_case_carriers"] = 5."""
    ctx = _make_context({"association_min_case_carriers": 5})
    assoc_config = _build_assoc_config_from_context(ctx)
    assert assoc_config.min_case_carriers == 5


@pytest.mark.unit
def test_min_case_carriers_default():
    """When association_min_case_carriers not in cfg, default is 10."""
    ctx = _make_context({})
    assoc_config = _build_assoc_config_from_context(ctx)
    assert assoc_config.min_case_carriers == 10


# ---------------------------------------------------------------------------
# Tests: JSON config path (association section)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_json_overrides_for_thresholds():
    """JSON association section values propagate when CLI keys are absent."""
    ctx = _make_context(
        {
            "association": {
                "min_cases": 50,
                "skat_method": "Burden",
                "max_case_control_ratio": 5.0,
                "min_case_carriers": 3,
            }
        }
    )
    assoc_config = _build_assoc_config_from_context(ctx)
    assert assoc_config.min_cases == 50
    assert assoc_config.skat_method == "Burden"
    assert assoc_config.max_case_control_ratio == 5.0
    assert assoc_config.min_case_carriers == 3


# ---------------------------------------------------------------------------
# Tests: CLI overrides JSON (precedence)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_cli_overrides_json():
    """CLI cfg key association_min_cases overrides JSON association section min_cases."""
    ctx = _make_context(
        {
            # CLI-set key (higher precedence for non-nullable fields)
            "association_min_cases": 300,
            # JSON association section (lower precedence)
            "association": {
                "min_cases": 50,
            },
        }
    )
    assoc_config = _build_assoc_config_from_context(ctx)
    # CLI wins: 300, not 50
    assert assoc_config.min_cases == 300


@pytest.mark.unit
def test_cli_skat_method_overrides_json():
    """CLI skat_method key overrides JSON association section skat_method."""
    ctx = _make_context(
        {
            # CLI-set key
            "skat_method": "SKATO",
            # JSON association section
            "association": {
                "skat_method": "Burden",
            },
        }
    )
    assoc_config = _build_assoc_config_from_context(ctx)
    # CLI wins: SKATO, not Burden
    assert assoc_config.skat_method == "SKATO"
