"""--permutations (#195): empirical p-values from a fixed-null parametric bootstrap.

References are closed form: Besag & Clifford (1991) stopping counts, and the exact
null law of their sequential p-value (P(p <= s) = s on its support), which makes the
randomized p exactly Uniform(0, 1) for the KS test.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest
import scipy.stats

from tests.unit.test_association_primary_test import _parse, _validate
from variantcentrifuge.association import permutation as perm
from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.engine import AssociationEngine


@pytest.mark.unit
class TestValidation:
    def test_unset_by_default_and_mapped(self):
        from variantcentrifuge.cli_association_config import populate_association_config

        args = _parse([])
        assert (args.permutations, args.permute_below_p, args.seed) == (None, None, None)
        cfg: dict = {"perform_association": True}
        argv = ["--permutations", "100", "--permute-below-p", "0.01", "--seed", "7"]
        populate_association_config(_parse(argv), cfg, pytest.fail)
        assert (cfg["permutations"], cfg["permute_below_p"], cfg["seed"]) == (100, 0.01, 7)

    @pytest.mark.parametrize(
        ("argv", "tests", "file_assoc"),
        [
            (["--permutations", "1000"], "logistic_burden,skat", None),
            (["--permutations", "10", "--permute-below-p", "1"], ["skat"], None),
            (["--permutations", "0"], ["skat"], None),
            (["--primary-test", "burden"], "fisher,logistic_burden", {"permutations": 5}),
        ],
    )
    def test_valid(self, argv, tests, file_assoc):
        assert _validate(argv, tests, file_assoc) is None

    @pytest.mark.parametrize(
        ("argv", "tests", "file_assoc", "fragment"),
        [
            (["--permutations", "9", "--trait-type", "quantitative"], ["skat"], None, "binary"),
            (["--permutations", "9"], ["fisher"], None, "'skat' or 'logistic_burden'"),
            (["--permutations", "-1"], ["skat"], None, ">= 0"),
            (["--permutations", "9", "--permute-below-p", "0"], ["skat"], None, "(0, 1]"),
            (["--permute-below-p", "0.01"], ["skat"], None, "require --permutations"),
            (["--seed", "3"], ["skat"], None, "require --permutations"),
            ([], ["fisher"], {"permutations": 9}, "'skat' or 'logistic_burden'"),
        ],
    )
    def test_invalid_fails_before_extraction(self, argv, tests, file_assoc, fragment):
        msg = _validate(argv, tests, file_assoc)
        assert msg is not None and fragment in msg

    def test_requires_perform_association(self):
        from variantcentrifuge.cli_association_config import validate_association_args

        args = _parse(["--permutations", "10"])
        args.perform_association = False
        args.association_tests = None
        with pytest.raises(SystemExit):
            validate_association_args(args, {}, MagicMock(side_effect=SystemExit))

    def test_config_section_keys_and_merge(self):
        from variantcentrifuge.stages.analysis_stage_helpers import (
            _build_assoc_config_from_context,
            _validate_association_config_dict,
        )

        _validate_association_config_dict({"permutations": 10, "permute_below_p": 0.5, "seed": 1})
        with pytest.raises(ValueError, match="permute_below_p"):
            _validate_association_config_dict({"permute_below_p": "x"})
        ctx = SimpleNamespace(config={"association": {"permutations": 10, "seed": 4}})
        cfg = _build_assoc_config_from_context(ctx)
        assert (cfg.permutations, cfg.permute_below_p, cfg.seed) == (10, 1e-4, 4)

    def test_statistic_choice(self):
        stat = perm.permutation_statistic
        assert stat(["logistic_burden", "skat"], "logistic_burden", "SKATO", "binary") == "burden"
        assert stat(["logistic_burden", "skat"], "acat_o", "SKATO", "binary") == "skat"
        assert stat(["skat"], "skat", "Burden", "binary") == "burden"
        assert stat(["fisher", "logistic_burden"], "acat_o", "SKAT", "binary") == "burden"
        assert perm.is_primary_statistic("skat", "skat", "SKAT")
        assert not perm.is_primary_statistic("skat", "skat", "SKATO")
        assert not perm.is_primary_statistic("skat", "acat_o", "SKAT")
        assert perm.is_primary_statistic("burden", "logistic_burden", "SKATO")


class _Scripted:
    """Statistic returning a scripted sequence: first the observed value, then replicates."""

    def __init__(self, t_obs: float, replicates: np.ndarray):
        self.t_obs, self.rest, self.batches = t_obs, replicates, []

    def __call__(self, s: np.ndarray) -> np.ndarray:
        if s.ndim == 1:
            return np.float64(self.t_obs)
        out, self.rest = self.rest[: s.shape[1]], self.rest[s.shape[1] :]
        self.batches.append(s.shape[1])
        return out


def _run(stat, n_perm: int, n: int = 50):
    mu = np.full(n, 0.3)
    return perm.empirical_p(
        stat, np.ones((n, 2)), mu, np.zeros(n), np.random.default_rng(0), n_perm
    )


@pytest.mark.unit
class TestBesagClifford:
    def test_stops_at_hth_exceedance_across_batches(self):
        reps = np.zeros(5000)
        hits = [3, 50, 127, 128, 129, 200, 300, 381, 383, 500]  # 0-based; 10th at index 500
        reps[hits] = 5.0
        stat = _Scripted(4.0, reps)
        assert _run(stat, 5000) == (10 / 501, 501)
        assert stat.batches == [128, 256, 512]  # one GEMM per batch, doubling

    def test_ties_count_as_exceedances(self):
        stat = _Scripted(4.0, np.full(100, 4.0))
        assert _run(stat, 100) == (1.0, perm.H_STOP)

    def test_n_reached_without_exceedance_never_zero(self):
        stat = _Scripted(4.0, np.zeros(300))
        assert _run(stat, 300) == (1 / 301, 300)
        assert stat.batches == [128, 172]

    def test_fewer_than_h_exceedances(self):
        reps = np.zeros(300)
        reps[[5, 6, 299]] = 9.0
        assert _run(_Scripted(4.0, reps), 300) == (4 / 301, 300)

    def test_real_gemm_extremes(self):
        n = 200
        mu = np.full(n, 0.2)
        rng = np.random.default_rng(1)
        a = np.ones((n, 1))
        stat = perm._STATISTICS["burden"]
        # observed score far beyond any replicate: N used, p = 1/(N+1)
        assert perm.empirical_p(stat, a, mu, np.full(n, 50.0), rng, 400) == (1 / 401, 400)
        # zero observed statistic: every replicate exceeds, stop at h
        assert perm.empirical_p(stat, a, mu, np.zeros(n), rng, 400) == (1.0, perm.H_STOP)


def _gene_dicts(n_genes: int, n: int, case_frac: float, seed: int, signal: str | None = None):
    rng = np.random.default_rng(seed)
    cov = np.column_stack([rng.normal(size=n), rng.integers(0, 2, n)])
    p_case = 1 / (1 + np.exp(-(np.log(case_frac / (1 - case_frac)) + 0.7 * cov[:, 0])))
    y = (rng.random(n) < p_case).astype(float)
    genes = []
    for i in range(n_genes):
        maf = rng.uniform(0.005, 0.04, 6)
        g = (rng.random((n, 6)) < maf).astype(float) + (rng.random((n, 6)) < maf)
        name = f"G{i:03d}"
        if name == signal:
            g[y == 1, 0] = 1.0  # every case carries variant 0
        built = {
            "genotype_matrix": g,
            "variant_mafs": g.mean(0) / 2,
            "phenotype_vector": y,
            "covariate_matrix": cov,
        }
        genes.append(
            {
                "GENE": name,
                "proband_count": int(y.sum()),
                "control_count": int(n - y.sum()),
                "n_qualifying_variants": 6,
                "_genotype_matrix_builder": lambda built=built: dict(built),
            }
        )
    return genes


def _engine_run(workers: int, seed: int = 0, **kw):
    config = AssociationConfig(
        association_workers=workers, permutations=2000, permute_below_p=0.5, seed=seed, **kw
    )
    engine = AssociationEngine.from_names(["logistic_burden", "skat"], config)
    return engine.run_all(_gene_dicts(60, 400, 0.2, seed=5))


@pytest.mark.unit
class TestEngine:
    def test_default_off_adds_no_columns(self):
        engine = AssociationEngine.from_names(["skat"], AssociationConfig())
        df = engine.run_all(_gene_dicts(3, 200, 0.3, seed=1))
        assert not {"empirical_p", "n_permutations", "permutation_method"} & set(df.columns)

    def test_workers_and_reruns_identical(self, caplog):
        with caplog.at_level(logging.WARNING, logger="variantcentrifuge"):
            one = _engine_run(1)
        assert "empirical_p is for the bootstrap_skat_q statistic" in caplog.text
        pd.testing.assert_frame_equal(one, _engine_run(4))
        pd.testing.assert_frame_equal(one, _engine_run(1))
        other_seed = _engine_run(1, seed=1)
        assert not np.array_equal(one["empirical_p"], other_seed["empirical_p"], equal_nan=True)
        below = pd.to_numeric(one["primary_pvalue"]) < 0.5
        assert one.loc[below, "empirical_p"].notna().all()
        assert one.loc[~below, "empirical_p"].isna().all()
        assert one.loc[~below, "n_permutations"].isna().all()
        assert (one.loc[below, "empirical_p"] > 0).all()
        assert set(one.loc[below, "permutation_method"]) == {"bootstrap_skat_q"}
        assert str(one["n_permutations"].dtype) == "Int64"

    def test_burden_primary_and_strong_signal(self):
        config = AssociationConfig(
            permutations=500,
            permute_below_p=1.0,
            primary_test="burden",
            association_workers=1,
            variant_weights="uniform",  # beta(1,25) would down-weight the common signal variant
        )
        engine = AssociationEngine.from_names(["logistic_burden", "skat"], config)
        df = engine.run_all(_gene_dicts(4, 400, 0.2, seed=2, signal="G001")).set_index("gene")
        assert set(df["permutation_method"]) == {"bootstrap_burden_score"}
        assert df.loc["G001", "empirical_p"] == 1 / 501
        assert df.loc["G001", "n_permutations"] == 500

    def test_untestable_genes_get_no_empirical_p(self):
        config = SimpleNamespace(
            variant_weights="beta:1,25", variant_weight_params=None, seed=0, permutations=10
        )
        empty = {"genotype_matrix": np.zeros((4, 0)), "phenotype_vector": np.ones(4)}
        assert perm._gene_empirical("A", lambda: empty, config, "skat", {}) is None
        tiny = {  # 3 samples, 3 design columns: no residual df
            "genotype_matrix": np.ones((3, 2)),
            "phenotype_vector": np.array([0.0, 1.0, 0.0]),
            "covariate_matrix": np.array([[1.0, 0.0], [0.0, 1.0], [2.0, 5.0]]),
            "variant_mafs": np.full(2, 0.5),
        }
        nulls: dict = {}
        assert perm._gene_empirical("B", lambda: tiny, config, "skat", nulls) is None
        assert list(nulls.values()) == [None]  # the failed fit is cached


@pytest.mark.unit
def test_null_calibration_imbalanced_with_covariates():
    """8:1 cases:controls, covariate-dependent outcome, no genetic effect: randomized p ~ U(0,1)."""
    n_perm = 300
    config = SimpleNamespace(
        variant_weights="beta:1,25", variant_weight_params=None, seed=3, permutations=n_perm
    )
    nulls: dict = {}
    ps = []
    for d in _gene_dicts(300, 450, 8 / 9, seed=9):
        r = perm._gene_empirical(d["GENE"], d["_genotype_matrix_builder"], config, "skat", nulls)
        if r is not None:
            ps.append(r[0])
    p = np.array(ps)
    h = perm.H_STOP
    support = np.unique(np.r_[np.arange(1, h + 1) / (n_perm + 1), h / np.arange(h, n_perm + 1)])
    idx = np.searchsorted(support, p - 1e-15)
    assert np.allclose(support[idx], p)
    prev = np.where(idx > 0, support[np.maximum(idx - 1, 0)], 0.0)
    u = prev + np.random.default_rng(0).random(p.size) * (p - prev)
    assert len(nulls) == 1 and p.size > 250
    assert scipy.stats.kstest(u, "uniform").pvalue > 0.01
