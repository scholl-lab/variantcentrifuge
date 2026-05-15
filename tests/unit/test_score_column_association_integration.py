"""Small association regression for score-column weights across burden and Python SKAT."""

from __future__ import annotations

from argparse import Namespace
from unittest.mock import Mock

import pandas as pd

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import AssociationAnalysisStage


def test_score_column_weights_run_logistic_burden_and_python_skat(tmp_path):
    cases = [f"CASE{i}" for i in range(10)]
    controls = [f"CTRL{i}" for i in range(10)]
    samples = cases + controls

    rows = []
    for variant_idx, score in enumerate([2.0, 5.0, 9.0]):
        row = {
            "GENE": "GENE1",
            "nephro_candidate_score": score,
            "ANN_0__EFFECT": "missense_variant",
        }
        for sample_idx, _sample in enumerate(samples):
            if variant_idx == 0:
                gt = "0/1" if sample_idx < 8 else "0/0"
            elif variant_idx == 1:
                gt = "0/1" if 4 <= sample_idx < 14 else "0/0"
            else:
                gt = "1/1" if sample_idx in (0, 1, 10) else "0/0"
            row[f"GEN_{sample_idx}__GT"] = gt
        rows.append(row)

    df = pd.DataFrame(rows)
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    workspace.intermediate_dir = tmp_path / "intermediate"
    workspace.intermediate_dir.mkdir()
    workspace.get_intermediate_path = lambda name: workspace.intermediate_dir / name
    workspace.get_output_path = lambda name, ext=".tsv": workspace.output_dir / f"{name}{ext}"
    workspace.base_name = "score_column_integration"

    context = PipelineContext(
        args=Namespace(),
        config={
            "perform_association": True,
            "case_samples": cases,
            "control_samples": controls,
            "association_tests": ["logistic_burden", "skat"],
            "skat_backend": "python",
            "variant_weights": "score_column",
            "variant_weight_column": "nephro_candidate_score",
            "variant_weight_params": {
                "score_min": 0,
                "score_max": 10,
                "floor": 0.1,
                "combine_with_beta": True,
            },
            "output_dir": str(tmp_path),
            "output_file_base": "score_column_integration",
            "gzip_intermediates": False,
        },
        workspace=workspace,
    )
    context.current_dataframe = df
    context.vcf_samples = samples

    result_context = AssociationAnalysisStage()._process(context)

    result_df = result_context.association_results
    assert result_df is not None
    assert not result_df.empty
    assert "logistic_burden_pvalue" in result_df.columns
    assert "skat_pvalue" in result_df.columns
    assert "acat_o_pvalue" in result_df.columns
