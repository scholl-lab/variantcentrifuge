"""Genotype missing-data policy: --missing-as-ref, visible site QC, artifact v2 (#200)."""

from __future__ import annotations

import argparse
import json
import logging
from argparse import Namespace
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from tests.unit.meta_helpers import make_stratum
from variantcentrifuge.association.genotype_matrix import build_genotype_matrix
from variantcentrifuge.association.meta.artifact import (
    V1_MISSING_POLICY,
    emit_stratum_artifact,
    load_artifact,
    save_artifact,
)
from variantcentrifuge.association.meta.manifest import (
    MetaCompatibilityError,
    validate_comparability,
)
from variantcentrifuge.cli_association_config import (
    populate_association_config,
    validate_association_args,
)
from variantcentrifuge.cli_association_groups import add_association_argument_group
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import _build_assoc_config_from_context
from variantcentrifuge.stages.association_stage_helpers import _report_missingness_qc
from variantcentrifuge.stages.genotype_matrix_builder import _GenotypeMatrixBuilder

N = 50


def _sparse_gene(n_sites: int = 4, carriers: int = 4) -> tuple[pd.DataFrame, list[str]]:
    """Merged single-sample VCFs: each site called only in its het carriers, rest ./.."""
    samples = [f"S{i}" for i in range(N)]
    rows = []
    for j in range(n_sites):
        hets = {(j * carriers + k) % N for k in range(carriers)}
        rows.append({s: "0/1" if i in hets else "./." for i, s in enumerate(samples)})
    df = pd.DataFrame(rows)
    df.insert(0, "CHROM", "chr1")
    df.insert(1, "POS", [100 + j for j in range(n_sites)])
    df.insert(2, "REF", "A")
    df.insert(3, "ALT", "G")
    return df, samples


def _builder(df, samples, *, missing_as_ref=False, site=0.10) -> _GenotypeMatrixBuilder:
    pheno = (np.arange(N) % 2).astype(float)
    return _GenotypeMatrixBuilder(
        df, samples, samples, True, site, 0.80, pheno, None, missing_as_ref=missing_as_ref
    )


# -- genotype matrix ----------------------------------------------------------------


def test_sparse_sites_dropped_by_default_and_kept_as_ref():
    df, samples = _sparse_gene()
    assert _builder(df, samples)()["n_variants_qc"] == 0  # all 92% missing -> dropped
    res = _builder(df, samples, missing_as_ref=True)()
    assert res["n_variants_qc"] == 4
    assert res["genotype_matrix"].sum() == 16  # exactly the observed het carriers


def test_raised_threshold_inflates_carriers_only_without_missing_as_ref():
    df, samples = _sparse_gene(n_sites=1)
    gt_cols = samples
    # Default policy (SKAT bestguess): AF_obs = 0.5 -> every ./. becomes a het carrier.
    geno, *_ = build_genotype_matrix(df, samples, gt_cols, missing_site_threshold=1.0)
    assert int((geno[:, 0] > 0).sum()) == N
    geno, *_ = build_genotype_matrix(
        df, samples, gt_cols, missing_site_threshold=1.0, missing_as_ref=True
    )
    assert int((geno[:, 0] > 0).sum()) == 4


def test_missing_as_ref_fills_partial_calls_and_keeps_haploid_ploidy():
    samples = ["a", "b", "c", "d"]
    df = pd.DataFrame({"a": ["./1", "1"], "b": ["0/1", "."], "c": ["0/0", "0"], "d": ["0/0", ""]})
    geno, mafs, *_ = build_genotype_matrix(df, samples, samples, missing_as_ref=True)
    np.testing.assert_array_equal(geno[:, 0], [0, 1, 0, 0])
    np.testing.assert_array_equal(geno[:, 1], [1, 0, 0, 0])
    assert mafs[1] == pytest.approx(0.25)  # haploid row: AN = 4, not 8


def test_empty_gene_reports_zero_qc_variants():
    df, samples = _sparse_gene()
    assert _builder(df.iloc[0:0], samples)()["n_variants_qc"] == 0


def test_variant_keep_mask_aligns_with_gene_df_rows_and_is_discarded():
    from variantcentrifuge.association.engine_workers import (
        _apply_builder_result_to_gene_data,
        _discard_per_gene_matrix_payload,
    )

    df, samples = _sparse_gene()
    df.loc[2, samples] = ["0/1" if i < 6 else "0/0" for i in range(N)]  # fully called row
    res = _builder(df, samples)()
    mask = res["variant_keep_mask"]
    assert mask.dtype == bool and mask.tolist() == [False, False, True, False]
    assert res["n_variants_qc"] == int(mask.sum())
    assert _builder(df.iloc[0:0], samples)()["variant_keep_mask"].shape == (0,)
    gene_data: dict = {}
    _apply_builder_result_to_gene_data(gene_data, res, "G1")
    assert gene_data["variant_keep_mask"] is mask
    _discard_per_gene_matrix_payload(gene_data)
    assert "variant_keep_mask" not in gene_data and gene_data["n_variants_qc"] == 1


# -- stage reporting ----------------------------------------------------------------


def test_report_adds_column_and_one_aggregate_warning(caplog):
    df, _samples = _sparse_gene()
    gene_data = {"G1": {"gene_df": df, "n_variants_qc": 0}, "G2": {"gene_df": df.iloc[:2]}}
    gene_data["G2"]["n_variants_qc"] = 2
    results = pd.DataFrame({"gene": ["G1", "G2", "G3"], "n_variants": [4, 2, 1], "p": [1, 1, 1]})
    cfg = Mock(missing_site_threshold=0.1)
    with caplog.at_level(logging.WARNING):
        _report_missingness_qc(results, gene_data, cfg)
    assert list(results.columns) == ["gene", "n_variants", "n_variants_qc", "p"]
    assert results["n_variants_qc"].tolist() == [0, 2, pd.NA]
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "removed 4 of 6 variants" in warnings[0].getMessage()
    assert "in 1 of 2 genes" in warnings[0].getMessage()


def test_report_is_silent_without_drops_or_matrix_tests(caplog):
    df, _ = _sparse_gene()
    results = pd.DataFrame({"gene": ["G1"], "n_variants": [4]})
    _report_missingness_qc(results, {"G1": {}}, Mock(missing_site_threshold=0.1))
    assert list(results.columns) == ["gene", "n_variants"]
    with caplog.at_level(logging.WARNING):
        _report_missingness_qc(
            results, {"G1": {"gene_df": df, "n_variants_qc": 4}}, Mock(missing_site_threshold=0.1)
        )
    assert results["n_variants_qc"].tolist() == [4]
    assert not caplog.records


# -- score artifact -----------------------------------------------------------------


def _emit(tmp_path, name, builder):
    return emit_stratum_artifact(
        str(tmp_path / name),
        {"G1": builder},
        stratum_id=name,
        genome_build="GRCh38",
        annotation_version="x",
        mask_id="m",
        mask_hash="h",
    )


def test_artifact_records_policy_and_full_an_with_missing_as_ref(tmp_path):
    df, samples = _sparse_gene()
    art = _emit(tmp_path, "s1", _builder(df, samples, missing_as_ref=True, site=0.2))
    gs = art.genes["G1"]
    np.testing.assert_array_equal(gs.an, np.full(4, 2.0 * N))
    np.testing.assert_array_equal(gs.ac, np.full(4, 4.0))
    policy = {"imputation": "missing_as_ref", "site_threshold": 0.2, "sample_threshold": 0.8}
    assert art.missing_policy == policy
    meta = json.loads((tmp_path / "s1.meta.json").read_text())
    assert meta["schema_version"] == 2
    assert meta["missing_policy"] == policy
    assert load_artifact(str(tmp_path / "s1")).missing_policy == policy


def test_v1_artifact_loads_as_bestguess(tmp_path):
    path = str(tmp_path / "old")
    save_artifact(make_stratum("old", seed=1, gene_keys={"G": ["1:100:A:G"]}), path)
    meta_file = tmp_path / "old.meta.json"
    meta = json.loads(meta_file.read_text())
    meta["schema_version"] = 1
    del meta["missing_policy"]
    meta_file.write_text(json.dumps(meta))
    art = load_artifact(path)
    assert art.schema_version == 1
    assert art.missing_policy == V1_MISSING_POLICY


def test_unknown_future_schema_refused(tmp_path):
    path = str(tmp_path / "new")
    save_artifact(make_stratum("new", seed=1, gene_keys={"G": ["1:100:A:G"]}), path)
    meta_file = tmp_path / "new.meta.json"
    text = meta_file.read_text().replace('"schema_version": 2', '"schema_version": 3')
    meta_file.write_text(text)
    with pytest.raises(ValueError, match="schema_version 3"):
        load_artifact(path)


# -- meta comparability ---------------------------------------------------------------


def _pair():
    keys = {"G": [f"1:{100 + j}:A:G" for j in range(3)]}
    return make_stratum("a", seed=1, gene_keys=keys), make_stratum("b", seed=2, gene_keys=keys)


def test_meta_refuses_mixed_imputation():
    a, b = _pair()
    b.missing_policy = {**b.missing_policy, "imputation": "missing_as_ref"}
    with pytest.raises(MetaCompatibilityError, match="a=bestguess, b=missing_as_ref"):
        validate_comparability([a, b])


def test_meta_warns_on_threshold_mismatch(caplog):
    a, b = _pair()
    a.missing_policy = {**a.missing_policy, "site_threshold": 0.1}
    b.missing_policy = {**b.missing_policy, "site_threshold": 0.2}
    with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
        validate_comparability([a, b])
    assert "site_threshold (a=0.1, b=0.2)" in caplog.text


# -- CLI / config ---------------------------------------------------------------------


def _parse(argv: list[str]) -> Namespace:
    parser = argparse.ArgumentParser()
    add_association_argument_group(parser)
    return parser.parse_args(argv)


def _assoc_config(cfg: dict):
    ctx = PipelineContext(args=Namespace(), config=cfg, workspace=Mock(spec=Workspace))
    return _build_assoc_config_from_context(ctx)


def test_cli_flags_map_to_config_and_win_over_json():
    args = _parse(
        ["--missing-site-threshold", "0.5", "--missing-sample-threshold", "1", "--missing-as-ref"]
    )
    args.perform_association = True
    cfg: dict = {"association": {"missing_site_threshold": 0.3, "missing_as_ref": False}}
    populate_association_config(args, cfg, pytest.fail)
    ac = _assoc_config(cfg)
    assert (ac.missing_site_threshold, ac.missing_sample_threshold, ac.missing_as_ref) == (
        0.5,
        1.0,
        True,
    )


def test_unset_cli_flags_leave_json_and_defaults():
    args = _parse([])
    args.perform_association = True
    cfg: dict = {"association": {"missing_site_threshold": 0.3, "missing_as_ref": True}}
    populate_association_config(args, cfg, pytest.fail)
    ac = _assoc_config(cfg)
    assert (ac.missing_site_threshold, ac.missing_sample_threshold, ac.missing_as_ref) == (
        0.3,
        0.8,
        True,
    )
    assert _assoc_config({}).missing_as_ref is False


@pytest.mark.parametrize("value", ["1.5", "-0.1", "abc"])
def test_threshold_flags_reject_values_outside_unit_interval(value):
    with pytest.raises(SystemExit):
        _parse(["--missing-site-threshold", value])


def test_missing_flags_require_perform_association():
    args = _parse(["--missing-as-ref"])
    args.perform_association = False
    args.association_tests = None
    fail = Mock(side_effect=SystemExit)
    with pytest.raises(SystemExit):
        validate_association_args(args, {}, fail)
    assert "--missing-as-ref requires --perform-association" in fail.call_args[0][0]


def test_json_missing_as_ref_must_be_boolean():
    with pytest.raises(ValueError, match="'missing_as_ref' must be a boolean"):
        _assoc_config({"association": {"missing_as_ref": "yes"}})


# -- review round 1 -------------------------------------------------------------------


def test_fisher_aggregators_honour_missing_as_ref():
    from variantcentrifuge.gene_burden import (
        _aggregate_gene_burden_from_columns,
        _aggregate_gene_burden_from_gt,
    )

    cases, ctrls = {"A", "B"}, {"C", "D"}
    cols = pd.DataFrame(
        {
            "GENE": ["G", "G"],
            "a": ["0/1", "./."],
            "b": ["./.", "./1"],
            "c": ["0/0", "."],
            "d": ["0/0", "0/0"],
        }
    )
    order, gt_cols = ["A", "B", "C", "D"], ["a", "b", "c", "d"]
    packed = pd.DataFrame({"GENE": ["G", "G"], "GT": ["A(0/1);B(./.)", "A(./.);B(./1);C(.)"]})
    for missing_as_ref, called, carriers in ((False, (0, 1), (0, 0)), (True, (2, 2), (1, 0))):
        for rec in (
            _aggregate_gene_burden_from_columns(
                cols, cases, ctrls, order, gt_cols, missing_as_ref=missing_as_ref
            )[0],
            _aggregate_gene_burden_from_gt(packed, cases, ctrls, missing_as_ref=missing_as_ref)[0],
        ):
            assert (rec["proband_called_count"], rec["control_called_count"]) == called
            assert (rec["proband_carrier_count"], rec["control_carrier_count"]) == carriers
    rec = _aggregate_gene_burden_from_columns(cols, cases, ctrls, order, gt_cols, True)[0]
    assert rec["control_called_chromosomes"] == 4  # C: "0/0" + "." -> diploid


def _stage_ctx(tmp_path, config: dict, df: pd.DataFrame, samples: list[str]):
    from variantcentrifuge.stages.analysis_stages import AssociationAnalysisStage  # noqa: F401

    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    workspace.get_output_path = lambda name, ext=".tsv": tmp_path / f"{name}{ext}"
    workspace.base_name = "s"
    ctx = PipelineContext(args=Namespace(), config=config, workspace=workspace)
    ctx.current_dataframe = df
    ctx.vcf_samples = samples
    return ctx


def _stage_df(n: int = 20) -> tuple[pd.DataFrame, list[str], list[str], list[str]]:
    cases = [f"CASE{i}" for i in range(n // 2)]
    ctrls = [f"CTRL{i}" for i in range(n // 2)]
    # Site 1 called in every sample; site 2 only in its two carriers (merged-VCF style).
    gts = {
        f"GEN_{i}__GT": ["0/1" if i % 5 == 0 else "0/0", "0/1" if i in (0, 11) else "./."]
        for i in range(n)
    }
    df = pd.DataFrame(
        {
            "CHROM": ["1", "1"],
            "POS": [10, 20],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
            "GENE": ["GENE1"] * 2,
            **gts,
        }
    )
    return df, cases, ctrls, cases + ctrls


@pytest.mark.parametrize("missing_as_ref", [False, True])
def test_stage_passes_missing_as_ref_to_fisher_and_builders(monkeypatch, tmp_path, missing_as_ref):
    from variantcentrifuge.association.engine import AssociationEngine
    from variantcentrifuge.stages.analysis_stages import AssociationAnalysisStage

    captured: list[dict] = []

    def fake_run_all(self, gene_burden_data, gene_hook=None):
        captured.extend(dict(d) for d in gene_burden_data)
        return pd.DataFrame({"gene": ["GENE1"], "n_variants": [2], "fisher_pvalue": [1.0]})

    monkeypatch.setattr(AssociationEngine, "run_all", fake_run_all)
    df, cases, ctrls, samples = _stage_df()
    config = {
        "perform_association": True,
        "case_samples": cases,
        "control_samples": ctrls,
        "association_tests": ["fisher", "logistic_burden"],
        "gzip_intermediates": False,
        "missing_as_ref": missing_as_ref or None,
    }
    AssociationAnalysisStage()._process(_stage_ctx(tmp_path, config, df, samples))
    (gene,) = captured
    assert gene["_genotype_matrix_builder"].missing_as_ref is missing_as_ref
    called = 10 if missing_as_ref else 1  # complete-case: only CASE0 is called at site 2
    assert gene["proband_called_count"] == called


def test_artifact_only_run_warns_about_dropped_sites(monkeypatch, tmp_path, caplog):
    from variantcentrifuge.association.engine import AssociationEngine
    from variantcentrifuge.stages.analysis_stages import AssociationAnalysisStage

    monkeypatch.setattr(
        AssociationEngine,
        "run_all",
        lambda self, d, gene_hook=None: pd.DataFrame({"gene": ["GENE1"]}),
    )
    df, cases, ctrls, samples = _stage_df()
    config = {
        "perform_association": True,
        "case_samples": cases,
        "control_samples": ctrls,
        "association_tests": ["fisher"],
        "gzip_intermediates": False,
        "score_artifact": str(tmp_path / "s.npz"),
        "reference": "GRCh38",
    }
    with caplog.at_level(logging.WARNING):
        AssociationAnalysisStage()._process(_stage_ctx(tmp_path, config, df, samples))
    assert "Score artifact: site missingness QC removed 1 of 2 variants" in caplog.text


@pytest.mark.parametrize("policy", ["drop", None, {}, {"imputation": "mean"}])
def test_v2_artifact_without_valid_policy_refused(tmp_path, policy):
    path = str(tmp_path / "x")
    save_artifact(make_stratum("x", seed=1, gene_keys={"G": ["1:100:A:G"]}), path)
    meta_file = tmp_path / "x.meta.json"
    meta = json.loads(meta_file.read_text())
    if policy == "drop":
        del meta["missing_policy"]
    else:
        meta["missing_policy"] = policy
    meta_file.write_text(json.dumps(meta))
    with pytest.raises(ValueError, match="missing_policy"):
        load_artifact(path)


@pytest.mark.parametrize(
    "value", [-0.1, 1.5, True, float("nan")], ids=["neg", "gt1", "bool", "nan"]
)
@pytest.mark.parametrize("key", ["missing_site_threshold", "missing_sample_threshold"])
def test_json_thresholds_range_checked(key, value):
    with pytest.raises(ValueError, match=f"'{key}' .*between 0 and 1"):
        _assoc_config({"association": {key: value}})
