"""Covariate pre-flight checks and logistic convergence-failure summary (#222).

All data are synthetic.  The issue's shape (a near-constant 0/1 indicator next to
unscaled PCs) is reproduced at 1,000 samples: one control has ``age_imputed = 1`` and
ten heavy-tailed PCs span roughly -600..830.  On v0.23.0 every logistic burden p-value
of that cohort is NA ("Both Logit and Firth failed to converge").
"""

from __future__ import annotations

import logging
from argparse import Namespace
from pathlib import Path
from typing import ClassVar
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.covariate_checks import (
    check_convergence_failures,
    preflight_covariates,
)
from variantcentrifuge.association.tests import logistic_burden
from variantcentrifuge.association.tests.logistic_burden import LogisticBurdenTest
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.association_analysis_stages import AssociationAnalysisStage

N = 1000
SAMPLES = [f"S{i:04}" for i in range(N)]
PCS = [f"PC{k}" for k in range(1, 11)]


def _cohort() -> tuple[np.ndarray, pd.DataFrame]:
    """Binary phenotype (10 % cases) and covariates age, sex, age_imputed, PC1..PC10."""
    rng = np.random.default_rng(0)
    y = np.zeros(N)
    y[: N // 10] = 1
    rng.shuffle(y)
    cov = pd.DataFrame(
        {
            "age": np.round(rng.uniform(18, 90, N)),
            "sex": rng.integers(0, 2, N).astype(float),
            "age_imputed": 0.0,
        },
        index=SAMPLES,
    )
    cov.iloc[np.flatnonzero(y == 0)[0], 2] = 1.0  # 1 of 1,000, a control
    cov[PCS] = rng.standard_t(2, (N, 10)) * 5  # unscaled, heavy-tailed
    return y, cov


def _context(tmp_path: Path, y: np.ndarray, cov: pd.DataFrame, n_genes: int = 10):
    rng = np.random.default_rng(1)
    rows = []
    for g in range(n_genes):
        geno = rng.random((N, 8)) < 0.01
        for v in range(8):
            rows.append(
                {
                    "CHROM": "1",
                    "POS": str(1000 + 100 * g + v),
                    "REF": "A",
                    "ALT": "C",
                    "GENE": f"G{g:02}",
                    "FEATUREID": f"TX{g}",
                    "EFFECT": "missense_variant",
                    **{f"GEN_{j}__GT": "0/1" if geno[j, v] else "0/0" for j in range(N)},
                }
            )
    tmp_path.mkdir(parents=True, exist_ok=True)
    cov_path = tmp_path / "covariates.tsv"
    cov.rename_axis("sample").to_csv(cov_path, sep="\t")
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    cfg = {
        "perform_association": True,
        "trait_type": "binary",
        "association_tests": ["logistic_burden"],
        "covariate_file": str(cov_path),
        "case_samples": [s for s, v in zip(SAMPLES, y, strict=True) if v == 1],
        "control_samples": [s for s, v in zip(SAMPLES, y, strict=True) if v == 0],
        "output_dir": str(tmp_path),
        "output_file_base": "t",
        "gzip_intermediates": False,
    }
    context = PipelineContext(args=Namespace(), config=cfg, workspace=workspace)
    context.current_dataframe = pd.DataFrame(rows)
    context.vcf_samples = list(SAMPLES)
    return context


def _run(context: PipelineContext) -> pd.DataFrame:
    return AssociationAnalysisStage()._process(context).association_results.set_index("gene")


# --- acceptance 1: near-constant indicator + unscaled PCs -----------------------------------


def test_near_constant_indicator_with_unscaled_pcs_converges(tmp_path, caplog):
    """v0.23.0: 10/10 genes NA.  Now age_imputed is dropped by name and every gene converges."""
    y, cov = _cohort()
    with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
        res = _run(_context(tmp_path, y, cov))
    assert res["logistic_burden_pvalue"].notna().all()
    assert "age_imputed (1)" in caplog.text
    assert not res["warnings"].str.contains("FIRTH_CONVERGE_FAIL").any()


def test_dropped_indicator_equals_run_without_it(tmp_path):
    """Dropping is exactly the run without the column (no other side effect)."""
    y, cov = _cohort()
    dropped = _run(_context(tmp_path / "a", y, cov))
    without = _run(_context(tmp_path / "b", y, cov.drop(columns="age_imputed")))
    pd.testing.assert_frame_equal(dropped, without)


# --- acceptance 3: separation and collinearity fail up front, naming columns ----------------


def test_separating_covariate_fails_naming_it(tmp_path):
    y, cov = _cohort()
    cov = cov.drop(columns="age_imputed")
    cov["screen"] = y + 0.5 * cov["sex"]  # cases all >= 1, controls all <= 0.5
    with pytest.raises(ValueError, match=r"separate cases from controls.*: screen\."):
        _run(_context(tmp_path, y, cov, n_genes=1))


def test_collinear_covariates_fail_naming_them(tmp_path):
    y, cov = _cohort()
    cov["PC1_copy"] = cov["PC1"] * 2.0 + 3.0
    with pytest.raises(ValueError, match=r"rank-deficient or collinear.*: PC1, PC1_copy\."):
        _run(_context(tmp_path, y, cov, n_genes=1))


def _cov(*cols: np.ndarray) -> np.ndarray:
    return np.column_stack(cols).astype(float)


def test_preflight_well_conditioned_returns_same_object():
    rng = np.random.default_rng(3)
    cov = _cov(rng.normal(50, 15, 200), rng.integers(0, 2, 200), rng.normal(0, 1, 200))
    y = (rng.random(200) < 0.3).astype(float)
    out, names = preflight_covariates(cov, ["age", "sex", "PC1"], y, binary=True)
    assert out is cov and names == ["age", "sex", "PC1"]
    assert preflight_covariates(None, [], y, binary=True) == (None, [])


def test_preflight_drop_boundary_and_all_dropped():
    y = np.r_[np.ones(10), np.zeros(90)]
    four, five = np.zeros(100), np.zeros(100)
    four[[0, 20, 30, 40]] = 1.0  # 4 non-modal values: dropped
    five[[0, 20, 30, 40, 50]] = 1.0  # 5: kept
    out, names = preflight_covariates(_cov(four, five), ["four", "five"], y, binary=True)
    assert names == ["five"] and np.array_equal(out, five[:, None])
    assert preflight_covariates(_cov(four), ["four"], y, binary=True) == (None, [])


def test_preflight_near_collinear_above_limit():
    rng = np.random.default_rng(4)
    a, noise = rng.normal(0, 1, 500), rng.normal(0, 1, 500)
    y = (rng.random(500) < 0.5).astype(float)
    ok = _cov(a, a + 1e-3 * noise, rng.normal(0, 1, 500))  # condition ~1e3
    assert preflight_covariates(ok, ["a", "b", "c"], y, binary=True)[0] is ok
    bad = _cov(a, a + 1e-9 * noise, rng.normal(0, 1, 500))  # condition ~1e9
    with pytest.raises(ValueError, match=r"condition number .*: a, b\."):
        preflight_covariates(bad, ["a", "b", "c"], y, binary=True)


def test_preflight_quasi_separating_indicator_and_quantitative_skip():
    y = np.r_[np.ones(20), np.zeros(80)]
    ind = np.zeros(100)
    ind[20:30] = 1.0  # 10 non-modal values, all controls: quasi-complete separation
    with pytest.raises(ValueError, match="ind"):
        preflight_covariates(_cov(ind), ["ind"], y, binary=True)
    assert preflight_covariates(_cov(ind), ["ind"], y, binary=False)[1] == ["ind"]
    # A single-class phenotype is left to the cohort checks.
    assert preflight_covariates(_cov(ind), ["ind"], np.zeros(100), binary=True)[1] == ["ind"]


# --- standardisation of ill-conditioned covariates for the logistic fit ---------------------


def _logistic(cov: np.ndarray, y: np.ndarray, geno: np.ndarray):
    data = {
        "genotype_matrix": geno,
        "variant_mafs": geno.mean(axis=0) / 2,
        "phenotype_vector": y,
        "covariate_matrix": cov,
    }
    return LogisticBurdenTest().run("G", data, AssociationConfig())


def test_ill_conditioned_units_are_scaled_beta_unchanged():
    """Reference by closed form: the burden coefficient of a logistic model with an
    intercept is invariant to affine covariate transforms, so age in 1e-7 units must
    give the beta/SE/p of age in years (condition ~1e2, fitted unscaled)."""
    rng = np.random.default_rng(5)
    n = 800
    age = rng.normal(50, 15, n)
    geno = (rng.random((n, 6)) < 0.05).astype(float)
    eta = -1.5 + 0.02 * (age - 50) + 0.4 * geno.sum(axis=1)
    y = (rng.random(n) < 1 / (1 + np.exp(-eta))).astype(float)
    years = _cov(age, rng.normal(0, 1, n))
    tiny = years * np.array([1e7, 1.0]) + np.array([1e9, 0.0])
    assert np.linalg.cond(np.column_stack([np.ones(n), tiny])) > logistic_burden._SCALE_CONDITION
    ref, got = _logistic(years, y, geno), _logistic(tiny, y, geno)
    assert got.effect_size == pytest.approx(ref.effect_size, rel=1e-6)
    assert got.se == pytest.approx(ref.se, rel=1e-6)
    assert got.p_value == pytest.approx(ref.p_value, rel=1e-6)


def test_well_conditioned_covariates_are_not_rescaled():
    cov = _cov(np.arange(10.0), np.arange(10.0) % 3)
    assert logistic_burden._fit_covariates(cov) is cov


# --- acceptance 2: widespread convergence failure fails the stage ---------------------------


def test_most_genes_failing_fails_stage_with_summary(tmp_path, monkeypatch):
    y, cov = _cohort()
    cov = cov.drop(columns="age_imputed")
    monkeypatch.setattr(logistic_burden, "_SEPARATION_BSE_THRESHOLD", 0.0)  # every fit "separates"
    monkeypatch.setattr(logistic_burden, "_firth_logistic", lambda *a, **k: None)
    with pytest.raises(RuntimeError, match=r"10 of 10 tested genes \(CONVERGENCE_SEPARATION=10\)"):
        _run(_context(tmp_path, y, cov))


def _results(pvalues: list, warnings: list[str]) -> pd.DataFrame:
    return pd.DataFrame({"logistic_burden_pvalue": pvalues, "warnings": warnings})


def test_convergence_failure_threshold(caplog):
    y = np.tile([1.0, 0.0], 10)
    cov = np.linspace(0, 1, 20)[:, None]
    half = _results([0.5, None], ["", "FIRTH_CONVERGE_FAIL;CONVERGENCE_MAX_ITER"])
    with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
        check_convergence_failures(half, ["x"], y, cov)  # 1 of 2 = 50 %: logged only
    assert "1 of 2 tested genes failed to converge (CONVERGENCE_MAX_ITER=1)" in caplog.text
    most = _results([0.5, None, None], ["", "CONVERGENCE_MAX_ITER", "CONVERGENCE_MAX_ITER"])
    with pytest.raises(
        RuntimeError, match=r"suspect covariates: none identified \(covariates: x\)"
    ):
        check_convergence_failures(most, ["x"], y, cov)
    with pytest.raises(RuntimeError, match="none \\(no covariates\\)"):
        check_convergence_failures(most, [], y, None)
    check_convergence_failures(
        _results([None], ["CONVERGENCE_MAX_ITER"]).drop(columns="logistic_burden_pvalue"),
        [],
        y,
        None,
    )  # no logistic burden: nothing to check
    check_convergence_failures(_results([0.1], [""]), ["x"], y, cov)


def test_suspect_covariates_named_from_null_fit():
    y = np.tile([1.0, 0.0, 0.0, 0.0, 0.0], 20)
    sep = np.zeros(100)
    sep[51] = 1.0  # one control: huge null-model SE
    cov = _cov(np.linspace(0, 1, 100), sep)
    most = _results([None], ["CONVERGENCE_SEPARATION"])
    with pytest.raises(RuntimeError, match=r"suspect covariates: sep\."):
        check_convergence_failures(most, ["lin", "sep"], y, cov)
    dup = _cov(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
    with pytest.raises(RuntimeError, match=r"suspect covariates: a, b\."):
        check_convergence_failures(most, ["a", "b"], y, dup)


@pytest.mark.parametrize(
    ("fit", "code"),
    [
        (np.linalg.LinAlgError("Singular matrix"), "CONVERGENCE_SINGULAR_HESSIAN"),
        (ValueError("boom"), "CONVERGENCE_FIT_ERROR"),
        ("not_converged", "CONVERGENCE_MAX_ITER"),
    ],
)
def test_na_gene_records_logit_failure_reason(monkeypatch, fit, code):
    import statsmodels.api as sm

    class _Fit:
        mle_retvals: ClassVar[dict] = {"converged": False}
        bse = np.array([0.1, 0.1])

    class _Logit:
        def __init__(self, *_a):
            pass

        def fit(self, **_k):
            if isinstance(fit, Exception):
                raise fit
            return _Fit()

    monkeypatch.setattr(sm, "Logit", _Logit)
    monkeypatch.setattr(logistic_burden, "_firth_logistic", lambda *a, **k: None)
    geno = np.zeros((20, 1))
    geno[:5] = 1
    res = _logistic(None, np.r_[np.ones(10), np.zeros(10)], geno)
    assert res.p_value is None
    assert res.extra["warnings"][-2:] == ["FIRTH_CONVERGE_FAIL", code]
