"""Score-artifact export from the association stage's matrix builders (--score-artifact)."""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd
import pytest

from tests.unit.meta_helpers import make_stratum
from variantcentrifuge.association.meta.artifact import (
    emit_stratum_artifact,
    fit_stratum_null,
    load_artifact,
    save_artifact,
)
from variantcentrifuge.association.meta.manifest import load_manifest_artifacts, read_study_manifest
from variantcentrifuge.association.meta.runner import run_meta
from variantcentrifuge.stages.genotype_matrix_builder import _GenotypeMatrixBuilder

GT = {0: "0/0", 1: "0/1", 2: "1/1"}


def _builder(n=200, seed=3):
    rng = np.random.default_rng(seed)
    pheno = (np.arange(n) % 3 == 0).astype(float)
    # v1 rare, v2 common ALT (AF>0.5 -> minor-allele recoded by build_genotype_matrix),
    # v3 is v1 repeated (transcript-split row).
    d1 = rng.binomial(2, 0.03, n)
    d2 = rng.binomial(2, 0.8, n)
    dosages = [d1, d2, d1]
    samples = [f"S{i}" for i in range(n)]
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr1", "chr1"],
            "POS": [100, 200, 100],
            "REF": ["A", "C", "A"],
            "ALT": ["G", "T", "G"],
            **{s: [GT[int(d[i])] for d in dosages] for i, s in enumerate(samples)},
        }
    )
    b = _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=samples,
        gt_columns=samples,
        is_binary=True,
        missing_site_threshold=0.1,
        missing_sample_threshold=0.8,
        phenotype_vector=pheno,
        covariate_matrix=None,
    )
    return b, np.column_stack([d1, d2]).astype(float), pheno


def test_emit_keeps_alt_coding_dedupes_and_hashes(tmp_path):
    b, alt, pheno = _builder()
    art = emit_stratum_artifact(
        str(tmp_path / "s1"),
        {"G1": b},
        stratum_id="s1",
        genome_build="GRCh38",
        annotation_version="GRCh38.99",
        mask_id="m",
        mask_hash="h",
    )
    gs = art.genes["G1"]
    assert gs.variant_keys == ["1:100:A:G", "1:200:C:T"]
    residuals = fit_stratum_null(pheno, None, "binary").extra["residuals"]
    # U must follow the ALT allele (the key), not the minor allele.
    np.testing.assert_allclose(gs.u, alt.T @ residuals, atol=1e-9)
    np.testing.assert_allclose(gs.ac, alt.sum(axis=0))
    assert art.n_cases == int(pheno.sum())
    assert art.sample_id_hashes[0] == hashlib.sha256(b"S0").hexdigest()
    assert load_artifact(str(tmp_path / "s1")).genes["G1"].variant_keys == gs.variant_keys


def test_manifest_relative_paths_and_names(tmp_path):
    keys = {"G": [f"1:{100 + j}:A:G" for j in range(5)]}
    for sid, seed in (("a", 1), ("b", 2)):
        save_artifact(make_stratum(sid, seed=seed, gene_keys=keys), str(tmp_path / sid))
    manifest = tmp_path / "study.tsv"
    manifest.write_text("stratum_id\tartifact_path\nadult\ta\nchild\tb\n")
    arts = load_manifest_artifacts(read_study_manifest(str(manifest)))
    assert [a.stratum_id for a in arts] == ["adult", "child"]
    df = run_meta(str(manifest), str(tmp_path / "out.tsv"))
    assert df.loc[0, "per_stratum_directions"].split(",")[0].startswith("adult:")


def test_unknown_schema_version_refused(tmp_path):
    keys = {"G": ["1:100:A:G"]}
    path = str(tmp_path / "a")
    save_artifact(make_stratum("a", seed=1, gene_keys=keys), path)
    meta = tmp_path / "a.meta.json"
    meta.write_text(meta.read_text().replace('"schema_version": 1', '"schema_version": 2'))
    with pytest.raises(ValueError, match="schema_version"):
        load_artifact(path)


def test_stage_writes_score_artifact_for_fisher_only_run(monkeypatch, tmp_path):
    from argparse import Namespace
    from unittest.mock import Mock

    from variantcentrifuge.association.engine import AssociationEngine
    from variantcentrifuge.pipeline_core.context import PipelineContext
    from variantcentrifuge.pipeline_core.workspace import Workspace
    from variantcentrifuge.stages.analysis_stages import AssociationAnalysisStage

    captured: list[dict] = []

    def fake_run_all(self, gene_burden_data):
        captured.extend(gene_burden_data)
        return pd.DataFrame({"gene": ["GENE1"], "fisher_pvalue": [1.0]})

    monkeypatch.setattr(AssociationEngine, "run_all", fake_run_all)
    cases = [f"CASE{i}" for i in range(10)]
    controls = [f"CTRL{i}" for i in range(10)]
    gts = {
        f"GEN_{i}__GT": ["0/1" if i % 4 == 0 else "0/0", "0/0" if i % 3 else "1/1"]
        for i in range(20)
    }
    df = pd.DataFrame(
        {
            "CHROM": ["1", "1"],
            "POS": [10, 20],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
            "GENE": ["GENE1", "GENE1"],
            **gts,
        }
    )
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    workspace.get_output_path = lambda name, ext=".tsv": tmp_path / f"{name}{ext}"
    workspace.base_name = "s"
    artifact = str(tmp_path / "cohortA.npz")
    config = {
        "perform_association": True,
        "case_samples": cases,
        "control_samples": controls,
        "association_tests": ["fisher"],
        "reference": "GRCh38.99",
        "filters": "(ANN[0].IMPACT has 'HIGH')",
        "score_artifact": artifact,
        "gzip_intermediates": False,
    }
    ctx = PipelineContext(args=Namespace(), config=config, workspace=workspace)
    ctx.current_dataframe = df
    ctx.vcf_samples = cases + controls

    AssociationAnalysisStage()._process(ctx)

    art = load_artifact(artifact)
    assert art.stratum_id == "cohortA" and art.genome_build == "GRCh38"
    assert art.genes["GENE1"].variant_keys == ["1:10:A:G", "1:20:C:T"]
    assert (art.n_cases, art.n_controls) == (10, 10)
    assert all("_genotype_matrix_builder" not in d for d in captured)
