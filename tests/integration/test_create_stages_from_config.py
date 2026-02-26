"""Integration tests for create_stages_from_config() config key mappings.

Verifies that config keys correctly activate the corresponding pipeline stages.
This tests the TD-02 fix: previously perform_association and perform_gene_burden
were not mapped, so passing them in config silently had no effect.
"""

import pytest

from variantcentrifuge.pipeline import create_stages_from_config
from variantcentrifuge.stages.analysis_stages import (
    AssociationAnalysisStage,
    GeneBurdenAnalysisStage,
)
from variantcentrifuge.stages.output_stages import (
    ExcelReportStage,
    HTMLReportStage,
    IGVReportStage,
)


def _stage_types(config: dict) -> list[type]:
    """Return list of stage classes produced by create_stages_from_config."""
    return [type(s) for s in create_stages_from_config(config)]


@pytest.mark.integration
def test_perform_association_true_activates_association_stage():
    """Config with perform_association=True must include AssociationAnalysisStage."""
    stage_types = _stage_types({"perform_association": True})
    assert AssociationAnalysisStage in stage_types, (
        "AssociationAnalysisStage missing when perform_association=True"
    )


@pytest.mark.integration
def test_perform_gene_burden_true_activates_gene_burden_stage():
    """Config with perform_gene_burden=True must include GeneBurdenAnalysisStage."""
    stage_types = _stage_types({"perform_gene_burden": True})
    assert GeneBurdenAnalysisStage in stage_types, (
        "GeneBurdenAnalysisStage missing when perform_gene_burden=True"
    )


@pytest.mark.integration
def test_empty_config_excludes_association_and_burden_stages():
    """Empty config must NOT include AssociationAnalysisStage or GeneBurdenAnalysisStage."""
    stage_types = _stage_types({})
    assert AssociationAnalysisStage not in stage_types, (
        "AssociationAnalysisStage should not appear with empty config"
    )
    assert GeneBurdenAnalysisStage not in stage_types, (
        "GeneBurdenAnalysisStage should not appear with empty config"
    )


@pytest.mark.integration
def test_perform_association_false_excludes_association_stage():
    """Explicit perform_association=False must not activate the association stage."""
    stage_types = _stage_types({"perform_association": False})
    assert AssociationAnalysisStage not in stage_types


@pytest.mark.integration
def test_perform_gene_burden_false_excludes_gene_burden_stage():
    """Explicit perform_gene_burden=False must not activate the gene burden stage."""
    stage_types = _stage_types({"perform_gene_burden": False})
    assert GeneBurdenAnalysisStage not in stage_types


@pytest.mark.integration
def test_xlsx_true_activates_excel_stage():
    """Config with xlsx=True must include ExcelReportStage."""
    stage_types = _stage_types({"xlsx": True})
    assert ExcelReportStage in stage_types


@pytest.mark.integration
def test_html_report_true_activates_html_stage():
    """Config with html_report=True must include HTMLReportStage."""
    stage_types = _stage_types({"html_report": True})
    assert HTMLReportStage in stage_types


@pytest.mark.integration
def test_igv_true_activates_igv_stage():
    """Config with igv=True must include IGVReportStage."""
    stage_types = _stage_types({"igv": True})
    assert IGVReportStage in stage_types
