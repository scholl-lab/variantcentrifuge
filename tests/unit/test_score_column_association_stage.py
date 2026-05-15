"""Stage-level tests for score-column variant weights."""

from __future__ import annotations

from argparse import Namespace
from unittest.mock import Mock

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import AssociationAnalysisStage


def _workspace(tmp_path):
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    workspace.intermediate_dir = tmp_path / "intermediate"
    workspace.intermediate_dir.mkdir()
    workspace.get_intermediate_path = lambda name: workspace.intermediate_dir / name
    workspace.get_output_path = lambda name, ext=".tsv": workspace.output_dir / f"{name}{ext}"
    workspace.base_name = "score_column_stage"
    return workspace


def _context(tmp_path, df, extra_config):
    cases = [f"CASE{i}" for i in range(10)]
    controls = [f"CTRL{i}" for i in range(10)]
    config = {
        "perform_association": True,
        "case_samples": cases,
        "control_samples": controls,
        "association_tests": ["logistic_burden"],
        "variant_weights": "score_column",
        "variant_weight_column": "nephro_candidate_score",
        "output_dir": str(tmp_path),
        "output_file_base": "score_column_stage",
        "gzip_intermediates": False,
    }
    config.update(extra_config)
    ctx = PipelineContext(args=Namespace(), config=config, workspace=_workspace(tmp_path))
    ctx.current_dataframe = df
    ctx.vcf_samples = cases + controls
    return ctx


def _df_with_score():
    sample_cols = {}
    for i in range(20):
        sample_cols[f"GEN_{i}__GT"] = ["0/1", "0/1", "0/0"]
    return pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1", "GENE1"],
            "nephro_candidate_score": [2.0, 5.0, 9.0],
            "dbNSFP_CADD_phred": [10.0, 20.0, 30.0],
            "dbNSFP_REVEL_score": [0.1, 0.2, 0.3],
            "ANN_0__EFFECT": ["missense_variant", "stop_gained", "frameshift_variant"],
            **sample_cols,
        }
    )


def test_stage_fails_when_score_column_missing(tmp_path):
    df = _df_with_score().drop(columns=["nephro_candidate_score"])
    context = _context(tmp_path, df, {})

    with pytest.raises(ValueError, match="variant weight column 'nephro_candidate_score'"):
        AssociationAnalysisStage()._process(context)


def test_stage_passes_score_column_to_lazy_builder(monkeypatch, tmp_path):
    from variantcentrifuge.association.engine import AssociationEngine

    captured_gene_data = []

    def fake_run_all(self, gene_burden_data):
        captured_gene_data.extend(gene_burden_data)
        return pd.DataFrame(
            {
                "gene": ["GENE1"],
                "n_cases": [10],
                "n_controls": [10],
                "n_variants": [3],
                "logistic_burden_pvalue": [1.0],
            }
        )

    monkeypatch.setattr(AssociationEngine, "run_all", fake_run_all)
    context = _context(tmp_path, _df_with_score(), {})

    AssociationAnalysisStage()._process(context)

    builder = captured_gene_data[0]["_genotype_matrix_builder"]
    assert builder.score_column == "nephro_candidate_score"
    assert builder.cadd_column is None
    assert builder.revel_column is None
    assert builder.effect_column == "ANN_0__EFFECT"


def test_stage_inline_column_spec_wins_over_variant_weight_column(monkeypatch, tmp_path):
    from variantcentrifuge.association.engine import AssociationEngine

    captured_gene_data = []

    def fake_run_all(self, gene_burden_data):
        captured_gene_data.extend(gene_burden_data)
        return pd.DataFrame(
            {
                "gene": ["GENE1"],
                "n_cases": [10],
                "n_controls": [10],
                "n_variants": [3],
                "logistic_burden_pvalue": [1.0],
            }
        )

    monkeypatch.setattr(AssociationEngine, "run_all", fake_run_all)
    df = _df_with_score()
    df["inline_score"] = [1.0, 2.0, 3.0]
    context = _context(
        tmp_path,
        df,
        {
            "variant_weights": "column:inline_score",
            "variant_weight_column": "nephro_candidate_score",
        },
    )

    AssociationAnalysisStage()._process(context)

    builder = captured_gene_data[0]["_genotype_matrix_builder"]
    assert builder.score_column == "inline_score"
