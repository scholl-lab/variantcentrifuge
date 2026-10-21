"""Quantitative association uses the numeric phenotype column (#187); results are
invariant to repeated identical transcript annotations (#188).

Synthetic data follow the public reproduction recipe in issue #187 (240 samples, three
genes x four variants, age covariate, outcome ``20 + (i*13) % 91 + i/100``).
"""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.association_analysis_stages import AssociationAnalysisStage
from variantcentrifuge.stages.association_outcome import build_quantitative_outcome

N = 240
SAMPLES = [f"S{i:03}" for i in range(N)]
GENES = [("SIGNAL", 20, 4), ("NULLGENE", 4, 12), ("PROTECT", 4, 48)]

# Reference: R 4.x `summary(lm(y ~ burden + age))`, burden = G %*% dbeta(MAF, 1, 25),
# MAF = colSums(G) / (2 * 240), on the recipe above (script: lm187.R in the C workstream
# scratch dir). SIGNAL also matches the issue's public LinearBurdenTest/R lm values.
R_LM = {
    "SIGNAL": (-0.045038660426, 0.305674712252, 0.882987512898),
    "NULLGENE": (-0.235649065651, 0.331832631112, 0.478312789204),
    "PROTECT": (-0.163229508309, 0.317794686975, 0.607987841519),
}


def _genotypes() -> dict[str, np.ndarray]:
    out = {}
    for gene, n_case, n_ctrl in GENES:
        g = np.zeros((N, 4))
        for v in range(4):
            g[v * n_case // 4 : (v + 1) * n_case // 4, v] = 1
            g[60 + v * n_ctrl // 4 : 60 + (v + 1) * n_ctrl // 4, v] = 1
        out[gene] = g
    return out


def _outcome(i: int, scale: float = 1.0, offset: float = 0.0) -> float:
    return scale * (20 + (i * 13) % 91 + i / 100) + offset


def _age(i: int) -> int:
    return 30 + (i * 17) % 51


def _variants_df(vcf_order: list[str]) -> pd.DataFrame:
    """Per-sample GT columns GEN_j__GT hold the genotype of ``vcf_order[j]``."""
    index = {s: i for i, s in enumerate(SAMPLES)}
    rows = []
    for gi, (gene, g) in enumerate(_genotypes().items()):
        for v in range(4):
            row = {
                "CHROM": "1",
                "POS": str(1000 + gi * 1000 + v * 10),
                "REF": "A",
                "ALT": "C",
                "GENE": gene,
                "FEATUREID": f"TX_{gene}",
                "EFFECT": "missense_variant",
            }
            for j, sample in enumerate(vcf_order):
                row[f"GEN_{j}__GT"] = "0/1" if g[index[sample], v] else "0/0"
            rows.append(row)
    return pd.DataFrame(rows)


def _context(tmp_path: Path, config: dict, vcf_order: list[str] | None = None) -> PipelineContext:
    vcf_order = vcf_order or SAMPLES
    workspace = Mock(spec=Workspace)
    tmp_path.mkdir(parents=True, exist_ok=True)
    workspace.output_dir = tmp_path
    covariates = tmp_path / "covariates.tsv"
    covariates.write_text(
        "sample\tage\n"
        + "".join(f"{s}\t{_age(i)}\n" for i, s in reversed(list(enumerate(SAMPLES))))
    )
    cfg = {
        "perform_association": True,
        "trait_type": "quantitative",
        "association_tests": ["linear_burden"],
        "covariate_file": str(covariates),
        "phenotype_value_column": "gfr",
        "output_dir": str(tmp_path),
        "output_file_base": "t",
        "gzip_intermediates": False,
        **config,
    }
    context = PipelineContext(args=Namespace(), config=cfg, workspace=workspace)
    context.current_dataframe = _variants_df(vcf_order)
    context.vcf_samples = list(vcf_order)
    context.phenotype_data = {s: {f"{_outcome(i):.10f}"} for i, s in enumerate(SAMPLES)}
    return context


def _run(context: PipelineContext) -> pd.DataFrame:
    result = AssociationAnalysisStage()._process(context)
    assert result.association_results is not None
    return result.association_results.set_index("gene")


def _assert_matches(results: pd.DataFrame, expected: dict, rel: float = 1e-9) -> None:
    for gene, (beta, se, p) in expected.items():
        row = results.loc[gene]
        assert row["linear_burden_beta"] == pytest.approx(beta, rel=rel)
        assert row["linear_burden_se"] == pytest.approx(se, rel=rel)
        assert row["linear_burden_pvalue"] == pytest.approx(p, rel=rel)


def _numpy_ols(keep: list[int]) -> dict[str, tuple[float, float, float]]:
    """Independent OLS (normal equations) of outcome ~ burden + age on samples ``keep``."""
    y = np.array([_outcome(i) for i in keep])
    age = np.array([_age(i) for i in keep], dtype=float)
    out = {}
    for gene, g in _genotypes().items():
        geno = g[keep]
        weights = stats.beta.pdf(geno.sum(axis=0) / (2 * len(keep)), 1, 25)
        x = np.column_stack([np.ones(len(keep)), geno @ weights, age])
        xtx_inv = np.linalg.inv(x.T @ x)
        coef = xtx_inv @ x.T @ y
        dof = len(keep) - x.shape[1]
        resid = y - x @ coef
        se = np.sqrt(resid @ resid / dof * xtx_inv[1, 1])
        out[gene] = (coef[1], se, 2 * stats.t.sf(abs(coef[1] / se), dof))
    return out


def test_quantitative_outcome_matches_r_lm_without_cohorts(tmp_path):
    results = _run(_context(tmp_path, {}))
    _assert_matches(results, R_LM)
    assert (results["n_cases"] == 0).all() and (results["n_controls"] == 0).all()
    assert (tmp_path / "t.association.tsv").exists()


def test_affine_outcome_scales_beta_and_se_keeps_p(tmp_path):
    context = _context(tmp_path, {})
    context.phenotype_data = {s: {f"{_outcome(i, 3, 17):.10f}"} for i, s in enumerate(SAMPLES)}
    expected = {g: (3 * b, 3 * se, p) for g, (b, se, p) in R_LM.items()}
    _assert_matches(_run(context), expected)


def test_shuffled_vcf_sample_order_is_aligned(tmp_path):
    order = list(np.random.default_rng(187).permutation(SAMPLES))
    _assert_matches(_run(_context(tmp_path, {}, vcf_order=order)), R_LM)


def test_missing_outcomes_excluded_from_matrix_phenotype_and_covariates(tmp_path, caplog):
    context = _context(tmp_path, {})
    missing = {0, 3, 61, 100, 239}
    context.phenotype_data[SAMPLES[0]] = set()  # empty cell
    context.phenotype_data[SAMPLES[3]] = {"NA"}
    context.phenotype_data[SAMPLES[61]] = {"."}
    context.phenotype_data[SAMPLES[100]] = {"nan"}
    del context.phenotype_data[SAMPLES[239]]  # not in phenotype file

    with caplog.at_level("WARNING"):
        results = _run(context)

    _assert_matches(results, _numpy_ols([i for i in range(N) if i not in missing]))
    assert "5 of 240 VCF sample(s) have no 'gfr' outcome" in caplog.text


def test_cohort_labels_do_not_replace_quantitative_outcome(tmp_path):
    cohorts = {"case_samples": SAMPLES[:60], "control_samples": ",".join(SAMPLES[60:])}
    _assert_matches(_run(_context(tmp_path, cohorts)), R_LM)


@pytest.mark.parametrize(
    ("config", "message"),
    [
        ({"association_tests": ["fisher", "linear_burden"]}, "need a binary trait"),
        ({"score_artifact": "x.npz"}, "binary traits only"),
    ],
)
def test_quantitative_rejects_binary_only_requests(tmp_path, config, message):
    with pytest.raises(ValueError, match=message):
        AssociationAnalysisStage()._process(_context(tmp_path, config))


def test_quantitative_without_phenotype_values_fails_clearly(tmp_path):
    context = _context(tmp_path, {})
    context.phenotype_data = None
    with pytest.raises(ValueError, match="needs the numeric outcome"):
        AssociationAnalysisStage()._process(context)


def test_quantitative_requires_vcf_samples(tmp_path):
    context = _context(tmp_path, {})
    context.vcf_samples = []
    with pytest.raises(ValueError, match="requires VCF sample IDs"):
        AssociationAnalysisStage()._process(context)


@pytest.mark.parametrize("bad", ["abc", "inf", "-Infinity", "1,5"])
def test_build_outcome_rejects_non_numeric_or_non_finite(bad):
    data = {"A": {"1.5"}, "B": {bad}}
    with pytest.raises(ValueError, match="non-numeric or non-finite"):
        build_quantitative_outcome(data, ["A", "B"], "gfr")


def test_build_outcome_rejects_conflicting_values():
    with pytest.raises(ValueError, match="conflicting"):
        build_quantitative_outcome({"A": {"1", "2"}}, ["A"], "gfr")


def test_build_outcome_requires_overlap_with_vcf():
    with pytest.raises(ValueError, match="No VCF sample has a quantitative outcome"):
        build_quantitative_outcome({"X": {"1"}}, ["A", "B"], "gfr")


def test_build_outcome_aligns_to_vcf_order():
    samples, values = build_quantitative_outcome(
        {"C": {"3"}, "A": {"1.0"}, "B": {" NA "}, "Z": {"9"}}, ["A", "B", "C"], "gfr"
    )
    assert samples == ["A", "C"]
    np.testing.assert_array_equal(values, [1.0, 3.0])


# --- binary path unchanged ---------------------------------------------------------------

_BINARY = {
    "trait_type": "binary",
    "association_tests": ["fisher", "logistic_burden", "skat"],
    "skat_backend": "python",
    "case_samples": SAMPLES[:60],
    "control_samples": SAMPLES[60:],
}


def test_binary_still_requires_cohorts(tmp_path):
    config = {**_BINARY, "case_samples": [], "control_samples": []}
    with pytest.raises(ValueError, match="both case and control cohorts are missing"):
        AssociationAnalysisStage()._process(_context(tmp_path, config))


def test_binary_outcome_is_cohort_membership_not_phenotype_values(tmp_path):
    with_values = _run(_context(tmp_path / "a", _BINARY))
    context = _context(tmp_path / "b", _BINARY)
    context.phenotype_data = None
    pd.testing.assert_frame_equal(with_values, _run(context))
    assert (with_values["n_cases"] == 60).all() and (with_values["n_controls"] == 180).all()


@pytest.mark.parametrize("trait", ["binary", "quantitative"])
def test_repeated_identical_transcript_annotation_is_invariant(tmp_path, trait):
    """#188: an extra qualifying transcript row for an existing allele changes nothing."""
    config = _BINARY if trait == "binary" else {"association_tests": ["linear_burden", "skat"]}
    base = _run(_context(tmp_path / "a", config))

    context = _context(tmp_path / "b", config)
    df = context.current_dataframe
    duplicate = df[df["GENE"] == "SIGNAL"].iloc[[0]].assign(FEATUREID="TX_SIGNAL_SECOND")
    context.current_dataframe = pd.concat([df, duplicate], ignore_index=True)
    assert len(context.current_dataframe) == 13

    pd.testing.assert_frame_equal(base, _run(context))
    assert base.loc["SIGNAL", "n_variants"] == 4


@pytest.mark.parametrize(
    ("n_cases", "n_controls", "ok", "message"),
    [
        (9, 100, False, "fewer than 10 cases"),
        (20, 30, True, "no practical power"),
        (60, 180, True, "underpowered for SKAT"),
        (300, 7000, True, "exceeding the 20:1"),
    ],
)
def test_binary_cohort_size_guard(caplog, n_cases, n_controls, ok, message):
    from variantcentrifuge.stages.association_analysis_stages import _cohort_size_ok

    with caplog.at_level("WARNING"):
        assert _cohort_size_ok(n_cases, n_controls) is ok
    assert message in caplog.text


def test_remove_sample_substring_applies_to_phenotype_ids(tmp_path):
    """Like binary phenotype matching, phenotype-file IDs get the same substring removal
    that SampleConfigLoadingStage applied to the VCF IDs (e.g. 'S000.bam' -> 'S000')."""
    context = _context(tmp_path, {"remove_sample_substring": ".bam"})
    context.phenotype_data = {f"{s}.bam": v for s, v in context.phenotype_data.items()}
    _assert_matches(_run(context), R_LM)


def test_build_outcome_substring_removal_merges_and_checks_conflicts():
    samples, values = build_quantitative_outcome(
        {"A.bam": {"1"}, "B": {"2"}}, ["A", "B"], "gfr", remove_substring=".bam"
    )
    assert samples == ["A", "B"]
    np.testing.assert_array_equal(values, [1.0, 2.0])
    with pytest.raises(ValueError, match="conflicting"):
        build_quantitative_outcome({"A.bam": {"1"}, "A": {"2"}}, ["A"], "gfr", ".bam")


# --- remaining stage responsibilities ------------------------------------------------------


@pytest.mark.parametrize("existing", ["gz", "plain", None])
def test_checkpoint_skip_restores_output_path(tmp_path, existing, caplog):
    context = _context(tmp_path, {})
    if existing:
        name = "t.association.tsv" + (".gz" if existing == "gz" else "")
        (tmp_path / name).write_text("")
    with caplog.at_level("WARNING"):
        result = AssociationAnalysisStage()._handle_checkpoint_skip(context)
    if existing:
        assert result.config["association_output"].endswith(name)
    else:
        assert "association_output" not in result.config
        assert "not found during checkpoint skip" in caplog.text
    off = _context(tmp_path, {"perform_association": False})
    assert (
        "association_output" not in AssociationAnalysisStage()._handle_checkpoint_skip(off).config
    )


def test_empty_dataframe_writes_header_only_output(tmp_path):
    context = _context(tmp_path, {"gzip_intermediates": True})
    context.current_dataframe = context.current_dataframe.iloc[0:0]
    result = AssociationAnalysisStage()._process(context)
    assert result.association_results.empty
    out = pd.read_csv(tmp_path / "t.association.tsv.gz", sep="\t")
    assert list(out.columns)[:4] == ["gene", "n_cases", "n_controls", "n_variants"]
    assert "linear_burden_pvalue" in out.columns
    assert AssociationAnalysisStage().get_output_files(result) == [
        tmp_path / "t.association.tsv.gz"
    ]


def test_missing_gene_column_and_missing_weight_column(tmp_path):
    context = _context(tmp_path, {})
    context.current_dataframe = context.current_dataframe.drop(columns="GENE")
    assert AssociationAnalysisStage()._process(context).association_results is None

    context = _context(tmp_path, {"variant_weights": "column:NOPE"})
    with pytest.raises(ValueError, match="'NOPE' is missing"):
        AssociationAnalysisStage()._process(context)


def test_missing_dataframe_raises(tmp_path):
    context = _context(tmp_path, {})
    context.current_dataframe = None
    with pytest.raises(RuntimeError, match="No DataFrame loaded"):
        AssociationAnalysisStage()._process(context)


@pytest.mark.parametrize("model", [None, "cadd"])
def test_coast_without_required_columns_is_skipped(tmp_path, model):
    only_coast = _context(tmp_path, {"association_tests": ["coast"], "coast_classification": model})
    assert AssociationAnalysisStage()._process(only_coast).association_results is None

    both = _context(tmp_path, {"association_tests": ["coast", "linear_burden"]})
    results = _run(both)
    assert "coast_pvalue" not in results.columns
    _assert_matches(results, R_LM)


def test_pca_file_from_pipeline_stage_is_used(tmp_path):
    pca = tmp_path / "pcs.tsv"
    pca.write_text(
        "sample\tPC1\tPC2\n"
        + "".join(f"{s}\t{np.sin(i):.6f}\t{np.cos(i * 0.7):.6f}\n" for i, s in enumerate(SAMPLES))
    )
    context = _context(tmp_path, {"pca_components": 2})
    context.stage_results["pca_computation"] = {"pca_file": str(pca)}
    results = _run(context)
    assert context.config["pca_file"] == str(pca)
    assert set(results.index) == set(R_LM)


@pytest.mark.parametrize("trait", ["binary", "quantitative"])
def test_diagnostics_output_written(tmp_path, trait):
    config = {**(_BINARY if trait == "binary" else {}), "diagnostics_output": str(tmp_path / "d")}
    _run(_context(tmp_path, config))
    assert any((tmp_path / "d").iterdir())


def test_binary_regression_from_packed_gt_column(tmp_path):
    """Packed GT strings are unpacked per sample for regression tests."""
    context = _context(tmp_path, {**_BINARY, "association_tests": ["logistic_burden"]})
    df = context.current_dataframe
    gt_cols = [f"GEN_{j}__GT" for j in range(N)]
    df["GT"] = [
        ";".join(f"{s}({row[c]})" for s, c in zip(SAMPLES, gt_cols, strict=True) if row[c] != "0/0")
        for _, row in df.iterrows()
    ]
    context.current_dataframe = df.drop(columns=gt_cols)
    packed = _run(context)
    per_sample = _run(
        _context(tmp_path / "ref", {**_BINARY, "association_tests": ["logistic_burden"]})
    )
    pd.testing.assert_series_equal(
        packed["logistic_burden_pvalue"], per_sample["logistic_burden_pvalue"]
    )


def test_cohort_size_abort_returns_without_results(tmp_path):
    config = {**_BINARY, "case_samples": SAMPLES[:5], "control_samples": SAMPLES[5:]}
    assert (
        AssociationAnalysisStage()._process(_context(tmp_path, config)).association_results is None
    )


def test_no_ratio_warning_without_controls(caplog):
    from variantcentrifuge.stages.association_analysis_stages import _cohort_size_ok

    with caplog.at_level("WARNING"):
        assert _cohort_size_ok(300, 0) is True
    assert "ratio" not in caplog.text


def test_stage_metadata():
    stage = AssociationAnalysisStage()
    assert stage.name == "association_analysis"
    assert "dataframe_loading" in stage.dependencies
    assert "chunked_analysis" in stage.soft_dependencies
    assert stage.parallel_safe is False
    assert stage.get_input_files(Mock()) == []


@pytest.mark.parametrize(
    ("cases", "controls", "message"),
    [([], SAMPLES[60:], "case cohort is missing"), (SAMPLES[:60], [], "control cohort is missing")],
)
def test_binary_single_missing_cohort_raises(tmp_path, cases, controls, message):
    config = {**_BINARY, "case_samples": cases, "control_samples": controls}
    with pytest.raises(ValueError, match=message):
        AssociationAnalysisStage()._process(_context(tmp_path, config))


def test_binary_comma_string_cohorts_are_parsed(tmp_path):
    config = {**_BINARY, "case_samples": ",".join(SAMPLES[:60]), "association_tests": ["fisher"]}
    context = _context(tmp_path, config)
    results = _run(context)
    assert context.config["case_samples"] == SAMPLES[:60]
    assert (results["n_cases"] == 60).all()


@pytest.mark.parametrize("tests", [["linear_burden"], ["linear_burden", "skat"]])
def test_untestable_gene_keeps_row_and_skip_code(tmp_path, tests):
    """A gene linear burden cannot test stays in the output with its warning code."""
    ctx = _context(tmp_path, {"association_tests": tests, "skat_backend": "python"})
    df = ctx.current_dataframe
    gcols = [c for c in df.columns if c.startswith("GEN_")]
    df.loc[df.GENE == "NULLGENE", gcols] = "0/1"  # every sample carries every variant
    res = _run(ctx)
    assert "NULLGENE" in res.index
    assert "BURDEN_NOT_ESTIMABLE" in res.loc["NULLGENE", "warnings"]
    assert pd.isna(res.loc["NULLGENE", "linear_burden_pvalue"])


# --- duplicate sample IDs after substring removal (final review item 1) -------------------


@pytest.mark.parametrize("restore", [False, True], ids=["fresh", "restore"])
@pytest.mark.parametrize("source", ["vcf", "from_tsv"])
def test_colliding_normalized_sample_ids_are_rejected(monkeypatch, source, restore):
    """'DUP.bam' and 'DUP' both become 'DUP' with --remove-sample-substring .bam; with equal
    outcomes the stage used to analyse one sample's genotypes twice (SIGNAL p 0.779->0.911)."""
    from variantcentrifuge import from_tsv
    from variantcentrifuge.stages import sample_and_pedigree_stages as sps

    raw = ["DUP.bam", *SAMPLES[1:80], "DUP"]
    config = {"remove_sample_substring": ".bam"}
    if source == "vcf":
        config["vcf_file"] = "in.vcf.gz"
        monkeypatch.setattr(sps, "get_vcf_samples", lambda _f: list(raw))
    else:
        config["from_tsv"] = "in.tsv"

        def prepare(ctx):  # mirrors from_tsv's own substring removal
            ctx.vcf_samples = [s.replace(".bam", "") for s in raw]

        monkeypatch.setattr(from_tsv, "prepare_from_tsv_context", prepare)
    context = PipelineContext(args=Namespace(), config=config, workspace=None)  # type: ignore[arg-type]
    stage = sps.SampleConfigLoadingStage()
    run = stage._handle_checkpoint_skip if restore else stage._process
    with pytest.raises(ValueError, match=r"not unique after removing substring '\.bam': DUP\."):
        run(context)


def test_non_colliding_substring_removal_still_loads(monkeypatch):
    from variantcentrifuge.stages import sample_and_pedigree_stages as sps

    monkeypatch.setattr(sps, "get_vcf_samples", lambda _f: ["A.bam", "B.bam", "C"])
    context = PipelineContext(
        args=Namespace(),
        config={"vcf_file": "in.vcf.gz", "remove_sample_substring": ".bam"},
        workspace=None,  # type: ignore[arg-type]
    )
    assert sps.SampleConfigLoadingStage()._process(context).vcf_samples == ["A", "B", "C"]


@pytest.mark.parametrize("tests", [["skat"], ["linear_burden"]])
def test_outcome_explained_by_covariates_reports_every_gene(tmp_path, tests):
    """An outcome equal to a covariate is untestable; genes stay listed with the reason."""
    ctx = _context(tmp_path, {"association_tests": tests, "skat_backend": "python"})
    ctx.phenotype_data = {s: {str(_age(i))} for i, s in enumerate(SAMPLES)}
    res = _run(ctx)
    assert sorted(res.index) == sorted(_genotypes())
    reason = "skat_skip_reason" if tests == ["skat"] else "warnings"
    assert res[reason].astype(str).str.contains("OUTCOME_EXPLAINED_BY_COVARIATES").all()


def test_skato_null_fit_skips_keep_skato_column(tmp_path):
    """SKAT-O rows skipped at the null fit keep the SKAT-O output column."""
    ctx = _context(
        tmp_path, {"association_tests": ["skat"], "skat_backend": "python", "skat_method": "SKATO"}
    )
    ctx.phenotype_data = {s: {str(_age(i))} for i, s in enumerate(SAMPLES)}
    res = _run(ctx)
    assert "skat_o_pvalue" in res.columns and "skat_pvalue" not in res.columns
    assert (res["skat_method"] == "SKATO").all()


def test_cohort_size_ok_imbalance_warning_is_two_sided(caplog):
    """Review C-r1 #1: stage warning also fires when cases outnumber controls."""
    from variantcentrifuge.stages.association_analysis_stages import _cohort_size_ok

    with caplog.at_level("WARNING", logger="variantcentrifuge"):
        assert _cohort_size_ok(1000, 10) is True
    assert "cases outnumber controls 100.0:1" in caplog.text


def test_covariate_case_control_table_written_by_stage(tmp_path, caplog):
    """Review C-r1 #2: stage writes covariate_case_control.tsv with raw batch labels and warns
    about a level with no cases (quasi-separation is left to Firth, #222)."""
    config = {
        **_BINARY,
        "diagnostics_output": str(tmp_path / "d"),
        "categorical_covariates": ["batch"],
    }
    context = _context(tmp_path, config)
    level = dict.fromkeys(SAMPLES[:10] + SAMPLES[100:220], "Y")
    level.update(dict.fromkeys(SAMPLES[220:], "Z"))  # 20 controls, 0 cases
    Path(context.config["covariate_file"]).write_text(
        "sample\tbatch\n" + "".join(f"{s}\t{level.get(s, 'X')}\n" for s in SAMPLES)
    )
    with caplog.at_level("WARNING", logger="variantcentrifuge"):
        _run(context)
    table = pd.read_csv(tmp_path / "d" / "covariate_case_control.tsv", sep="\t")
    assert table.to_dict("records") == [
        {"covariate": "batch", "level": "X", "n_cases": 50, "n_controls": 40},
        {"covariate": "batch", "level": "Y", "n_cases": 10, "n_controls": 120},
        {"covariate": "batch", "level": "Z", "n_cases": 0, "n_controls": 20},
    ]
    assert "level 'Z' has 0 cases and 20 controls; possible confounding" in caplog.text


def test_cohort_size_ok_uses_configured_ratio(caplog):
    """Ratio 7.8 warns when max_case_control_ratio=5, and stays silent when imbalance is None."""
    from variantcentrifuge.stages.association_analysis_stages import _cohort_size_ok

    with caplog.at_level("WARNING", logger="variantcentrifuge"):
        assert _cohort_size_ok(788, 101, 5.0) is True
    assert "cases outnumber controls 7.8:1" in caplog.text
    caplog.clear()
    with caplog.at_level("WARNING", logger="variantcentrifuge"):
        _cohort_size_ok(788, 101, None)  # diagnostics emit it instead (no duplicate)
    assert "outnumber" not in caplog.text
