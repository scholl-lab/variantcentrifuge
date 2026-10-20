"""Tests for SKAT null model invalidation and unassigned sample exclusion (#140)."""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.engine import AssociationEngine
from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest
from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.association_analysis_stages import AssociationAnalysisStage


def test_consecutive_genes_with_different_sample_counts():
    """Verify that consecutive genes with different sample counts refit null model without error."""
    test = PurePythonSKATTest()
    test.check_dependencies()
    config = AssociationConfig(trait_type="binary")

    # Gene 1: 20 samples
    y1 = np.array([1.0] * 10 + [0.0] * 10)
    g1 = np.random.default_rng(42).binomial(2, 0.2, (20, 4)).astype(float)
    data1 = {
        "phenotype_vector": y1,
        "genotype_matrix": g1,
        "covariate_matrix": None,
    }
    res1 = test.run("GENE1", data1, config)
    assert res1.p_value is not None

    # Gene 2: sample 0 dropped (19 samples)
    y2 = y1[1:]
    g2 = g1[1:, :]
    data2 = {
        "phenotype_vector": y2,
        "genotype_matrix": g2,
        "covariate_matrix": None,
    }
    res2 = test.run("GENE2", data2, config)
    assert res2.p_value is not None

    # Gene 2 with a fresh test instance should give the identical p-value
    test_fresh = PurePythonSKATTest()
    test_fresh.check_dependencies()
    res2_fresh = test_fresh.run("GENE2", data2, config)
    assert res2.p_value == pytest.approx(res2_fresh.p_value, rel=1e-6)

    # Gene 3: same cohort as Gene 2 (19 samples) reuses valid null model
    saved_null = test._null_model
    data3 = {
        "phenotype_vector": y2,
        "genotype_matrix": np.roll(g2, 1, axis=1),
        "covariate_matrix": None,
    }
    res3 = test.run("GENE3", data3, config)
    assert res3.p_value is not None
    assert test._null_model is saved_null, "Unchanged cohort must reuse existing null model"
    res3_fresh = test_fresh.run("GENE3", data3, config)
    assert res3.p_value == pytest.approx(res3_fresh.p_value, rel=1e-6)


def test_equal_sized_masks_retaining_different_samples_match_fresh_fit():
    """Verify equal-sized masks match fresh-fit SKAT and ACAT-V results."""
    test = PurePythonSKATTest()
    test.check_dependencies()
    config = AssociationConfig(trait_type="binary")

    rng = np.random.default_rng(123)
    y = np.array([1.0] * 10 + [0.0] * 10)
    cov = rng.normal(0, 1, (20, 2))
    g = rng.binomial(2, 0.25, (20, 5)).astype(float)

    # Drop sample 0 (samples 1..19 retained)
    mask_a = np.ones(20, dtype=bool)
    mask_a[0] = False
    data_a = {
        "phenotype_vector": y[mask_a],
        "genotype_matrix": g[mask_a],
        "covariate_matrix": cov[mask_a],
    }

    # Drop sample 19 (samples 0..18 retained)
    mask_b = np.ones(20, dtype=bool)
    mask_b[19] = False
    data_b = {
        "phenotype_vector": y[mask_b],
        "genotype_matrix": g[mask_b],
        "covariate_matrix": cov[mask_b],
    }

    res_a = test.run("GENE_A", data_a, config)
    res_b = test.run("GENE_B", data_b, config)

    # Compare with fresh fits
    fresh_a = PurePythonSKATTest()
    fresh_a.check_dependencies()
    res_fresh_a = fresh_a.run("GENE_A", data_a, config)

    fresh_b = PurePythonSKATTest()
    fresh_b.check_dependencies()
    res_fresh_b = fresh_b.run("GENE_B", data_b, config)

    assert res_a.p_value == pytest.approx(res_fresh_a.p_value, rel=1e-6)
    assert res_b.p_value == pytest.approx(res_fresh_b.p_value, rel=1e-6)
    if "acat_v_p" in res_a.extra and res_a.extra["acat_v_p"] is not None:
        assert res_a.extra["acat_v_p"] == pytest.approx(res_fresh_a.extra["acat_v_p"], rel=1e-6)
    if "acat_v_p" in res_b.extra and res_b.extra["acat_v_p"] is not None:
        assert res_b.extra["acat_v_p"] == pytest.approx(res_fresh_b.extra["acat_v_p"], rel=1e-6)


def test_unassigned_samples_excluded_from_regression_inputs(tmp_path: Path):
    """Verify unassigned samples are excluded from regression inputs."""
    workspace = Workspace(tmp_path / "work", "test")
    # 10 cases, 10 controls, and 2 unassigned samples interleaved
    cases = [f"CASE_{i}" for i in range(10)]
    controls = [f"CTRL_{i}" for i in range(10)]
    unassigned = ["EXCLUDED_1", "EXCLUDED_2"]
    all_vcf_samples = [
        *cases[:5],
        unassigned[0],
        *cases[5:],
        *controls[:5],
        unassigned[1],
        *controls[5:],
    ]

    config = {
        "perform_association": True,
        "output_dir": str(tmp_path),
        "association": {
            "association_tests": ["skat_python"],
            "trait_type": "binary",
        },
        "case_samples": cases,
        "control_samples": controls,
    }

    # Create dummy dataframe with per-sample genotype columns
    data = {
        "CHROM": ["chr1"] * 5,
        "POS": [100, 200, 300, 400, 500],
        "REF": ["A"] * 5,
        "ALT": ["T"] * 5,
        "GENE": ["TEST_GENE"] * 5,
    }
    for idx, s in enumerate(all_vcf_samples):
        col = f"GEN_{idx}__GT"
        if "CASE" in s:
            data[col] = ["0/1", "0/1", "0/1", "0/0", "0/1"]
        elif "CTRL" in s:
            data[col] = ["0/0", "0/1", "0/0", "0/0", "0/0"]
        else:
            # Unassigned samples
            data[col] = ["1/1", "1/1", "1/1", "1/1", "1/1"]

    df = pd.DataFrame(data)
    context = PipelineContext(Namespace(), config, workspace)
    context.current_dataframe = df
    context.vcf_samples = all_vcf_samples

    stage = AssociationAnalysisStage()
    stage._process(context)

    # In results, verify association test ran and unassigned samples were not treated as controls
    assert context.association_results is not None
    assert not context.association_results.empty

    # Explicitly pre-subsetted cohort: only cases and controls in identical relative order
    subset_samples = [s for s in all_vcf_samples if s not in unassigned]
    subset_data = {
        "CHROM": ["chr1"] * 5,
        "POS": [100, 200, 300, 400, 500],
        "REF": ["A"] * 5,
        "ALT": ["T"] * 5,
        "GENE": ["TEST_GENE"] * 5,
    }
    for idx, s in enumerate(subset_samples):
        orig_col = f"GEN_{all_vcf_samples.index(s)}__GT"
        subset_data[f"GEN_{idx}__GT"] = data[orig_col]

    workspace_sub = Workspace(tmp_path / "work_sub", "test_sub")
    context_sub = PipelineContext(Namespace(), config, workspace_sub)
    context_sub.current_dataframe = pd.DataFrame(subset_data)
    context_sub.vcf_samples = subset_samples

    stage_sub = AssociationAnalysisStage()
    stage_sub._process(context_sub)

    assert context_sub.association_results is not None
    p_interleaved = context.association_results["skat_pvalue"].values[0]
    p_subsetted = context_sub.association_results["skat_pvalue"].values[0]
    assert p_interleaved == pytest.approx(p_subsetted, rel=1e-6)


def test_serial_and_parallel_workers_produce_identical_results_on_masked_cohorts():
    """Verify that serial and parallel worker execution match on masked cohorts."""
    rng = np.random.default_rng(42)
    y = np.array([1.0] * 10 + [0.0] * 10)
    g = rng.binomial(2, 0.25, (20, 5)).astype(float)
    cov = rng.normal(0, 1, (20, 2))

    def make_data():
        data_list = []
        for i in range(5):
            m = np.ones(20, dtype=bool)
            if i > 0:
                m[i] = False
            data_list.append(
                {
                    "GENE": f"GENE{i + 1}",
                    "phenotype_vector": y[m].copy(),
                    "genotype_matrix": np.roll(g, i, axis=1)[m].copy(),
                    "covariate_matrix": cov[m].copy(),
                }
            )
        return data_list

    cfg_serial = AssociationConfig(trait_type="binary", association_workers=1)
    engine_serial = AssociationEngine.from_names(["skat_python"], cfg_serial)
    df_serial = engine_serial.run_all(make_data())

    cfg_parallel = AssociationConfig(trait_type="binary", association_workers=2)
    engine_parallel = AssociationEngine.from_names(["skat_python"], cfg_parallel)
    df_parallel = engine_parallel.run_all(make_data())

    assert not df_serial.empty and not df_parallel.empty
    assert list(df_serial["gene"]) == list(df_parallel["gene"])
    for p_ser, p_par in zip(df_serial["skat_pvalue"], df_parallel["skat_pvalue"], strict=True):
        assert p_ser == pytest.approx(p_par, rel=1e-6)
