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
            (["--permutations", "9", "--seed", "-1"], ["skat"], None, "--seed must be >= 0"),
            ([], ["skat"], {"permutations": 9, "seed": -2}, "--seed must be >= 0"),
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
        for key in ("permutations", "seed", "permute_below_p"):  # bool is an int in Python
            with pytest.raises(ValueError, match=key):
                _validate_association_config_dict({key: True})
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
        # selection is on the permuted statistic's own p (SKAT), not the ACAT-O primary
        below = pd.to_numeric(one["skat_pvalue"]) < 0.5
        assert (below != (pd.to_numeric(one["primary_pvalue"]) < 0.5)).any()
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

    def test_gene_without_permuted_test_p_is_not_permuted(self):
        """A gene SKAT-O skips (one variant: rank_deficient) gets no SKAT-Q empirical p."""
        genes = _gene_dicts(3, 400, 0.2, seed=4)
        one_var = genes[0]["_genotype_matrix_builder"]()
        one_var["genotype_matrix"] = one_var["genotype_matrix"][:, :1]
        one_var["variant_mafs"] = one_var["variant_mafs"][:1]
        genes[0]["_genotype_matrix_builder"] = lambda: dict(one_var)
        config = AssociationConfig(
            permutations=200, permute_below_p=1.0, skat_method="SKATO", association_workers=1
        )
        df = AssociationEngine.from_names(["logistic_burden", "skat"], config).run_all(genes)
        row = df.set_index("gene").loc["G000"]
        assert row["skat_skip_reason"] == "rank_deficient"
        assert np.isnan(row["empirical_p"]) and pd.isna(row["n_permutations"])
        assert row["permutation_method"] == "skipped_rank_deficient"
        assert df.set_index("gene").loc["G001", "permutation_method"] == "bootstrap_skat_q"

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


@pytest.mark.unit
def test_skato_selects_on_skat_q_own_p():
    """Under SKAT-O, genes are selected on SKAT Q's p (the permuted statistic), not SKAT-O's."""
    genes = _gene_dicts(40, 2000, 0.2, seed=11)
    for k, gd in enumerate(genes[:8]):  # burden-driven: SKAT-O p well below SKAT Q p
        b = gd["_genotype_matrix_builder"]()
        y, g, rng = (
            b["phenotype_vector"],
            b["genotype_matrix"].copy(),
            np.random.default_rng(100 + k),
        )
        for j in range(g.shape[1]):
            add = (y == 1) & (rng.random(y.size) < 0.012)
            g[add, j] = np.maximum(g[add, j], 1)
        b.update(genotype_matrix=g, variant_mafs=g.mean(0) / 2)
        gd["_genotype_matrix_builder"] = lambda b=b: dict(b)

    def run(method):
        config = AssociationConfig(
            permutations=2000, permute_below_p=1e-3, seed=1, skat_method=method
        )
        engine = AssociationEngine.from_names(["logistic_burden", "skat"], config)
        return engine.run_all([dict(d) for d in genes]).set_index("gene")

    skato, skat = run("SKATO"), run("SKAT")
    assert ((skato["skat_o_pvalue"] < 1e-3) & (skat["skat_pvalue"] >= 1e-3)).any()  # repro shape
    np.testing.assert_array_equal(skato["empirical_p"], skat["empirical_p"])
    assert skato["empirical_p"].notna().sum() == (skat["skat_pvalue"] < 1e-3).sum() > 0


@pytest.mark.unit
def test_skato_selection_ignores_reported_skato_p():
    """A conservative reported SKAT-O p (> 7P) must not hide a gene whose SKAT Q p is < P."""
    genes = _gene_dicts(3, 400, 0.2, seed=2, signal="G001")
    builders = {d["GENE"]: d["_genotype_matrix_builder"] for d in genes}
    frame = pd.DataFrame({"gene": list(builders), "skat_o_pvalue": [0.9, 0.9, np.nan]})
    config = AssociationConfig(
        permutations=300, permute_below_p=1e-2, skat_method="SKATO", association_workers=1
    )
    out = perm.add_empirical_p(frame, builders, ["skat"], config, "acat_o").set_index("gene")
    assert out.loc["G001", "empirical_p"] > 0  # SKAT Q p 5e-4 < P: selected and permuted
    assert np.isnan(out.loc["G000", "empirical_p"])  # null gene: SKAT Q p >= P
    assert out.loc["G002", "permutation_method"] == "skipped_no_pvalue"


@pytest.mark.unit
def test_fitted_probabilities_of_one_do_not_give_the_floor():
    """Five outlier cases get mu_hat == 1.0 exactly (no separation): v0.24 dev gave 1/(N+1)."""
    rng = np.random.default_rng(0)
    n = 1000
    x = rng.normal(size=n)
    y = (rng.random(n) < 1 / (1 + np.exp(-x))).astype(float)
    x[:5], y[:5] = 200.0, 1.0
    g = (rng.random((n, 4)) < 0.02).astype(float)
    built = {
        "genotype_matrix": g,
        "variant_mafs": g.mean(0) / 2,
        "phenotype_vector": y,
        "covariate_matrix": x[:, None],
    }
    config = SimpleNamespace(
        variant_weights="uniform", variant_weight_params=None, seed=0, permutations=1000
    )
    nulls: dict = {}
    p, used = perm._gene_empirical("A", lambda: dict(built), config, "skat", nulls)
    assert (next(iter(nulls.values())).extra["mu_hat"] == 1.0).sum() == 5  # repro shape
    assert p > 0.05 and used < 1000


@pytest.mark.unit
def test_nonfinite_statistic_is_skipped_not_floored(monkeypatch):
    n = 50
    stat, a, mu = perm._STATISTICS["skat"], np.ones((n, 2)), np.full(n, 0.3)
    e_obs = np.zeros(n)
    e_obs[0] = np.nan
    p, used = perm.empirical_p(stat, a, mu, e_obs, np.random.default_rng(0), 100)
    assert np.isnan(p) and used == 0
    mu_bad = mu.copy()
    mu_bad[0] = np.nan  # a replicate statistic turns NaN
    p, used = perm.empirical_p(stat, a, mu_bad, np.zeros(n), np.random.default_rng(0), 100)
    assert np.isnan(p) and used == 0

    genes = _gene_dicts(2, 400, 0.2, seed=4)
    builders = {d["GENE"]: d["_genotype_matrix_builder"] for d in genes}
    frame = pd.DataFrame({"gene": list(builders), "skat_pvalue": [0.01, 0.9]})
    config = AssociationConfig(permutations=100, permute_below_p=0.5, association_workers=1)
    monkeypatch.setattr(perm, "empirical_p", lambda *args, **kw: (np.nan, 0))
    out = perm.add_empirical_p(frame, builders, ["skat"], config, "skat").set_index("gene")
    assert np.isnan(out.loc["G000", "empirical_p"]) and pd.isna(out.loc["G000", "n_permutations"])
    assert out.loc["G000", "permutation_method"] == "skipped_nonfinite"
    assert pd.isna(out.loc["G001", "permutation_method"])


@pytest.mark.unit
def test_null_cache_is_bounded():
    """One null fit per distinct sample mask: the cache keeps at most 8 (Codex: 226 MiB at 64)."""
    from variantcentrifuge.association.meta.artifact import _MAX_CACHED_NULLS

    config = SimpleNamespace(
        variant_weights="beta:1,25", variant_weight_params=None, seed=0, permutations=20
    )
    genes = _gene_dicts(1, 200, 0.3, seed=1)
    base = genes[0]["_genotype_matrix_builder"]()
    nulls: dict = {}
    for k in range(3 * _MAX_CACHED_NULLS):
        b = {key: val[k:] if key != "variant_mafs" else val for key, val in base.items()}
        perm._gene_empirical(f"G{k}", lambda b=b: b, config, "skat", nulls)
        assert len(nulls) <= _MAX_CACHED_NULLS
    assert len(nulls) == _MAX_CACHED_NULLS


@pytest.mark.unit
def test_robust_skat_and_skato_select_on_the_same_skat_q_p():
    """Fixture r19/le10: robust SKAT p 1.6e-4, SKAT Q's analytic p 1.1e-6 (Codex review).

    v0.24 dev selected robust SKAT on its calibrated p (skipped at P = 1e-5) but robust
    SKAT-O on SKAT Q's analytic p (permuted): the same bootstrap, different genes.
    """
    from tests.unit.test_skat_robust_parity import _load

    geno, _w, y, _mu, design = _load("r19", "le10")
    built = {
        "genotype_matrix": geno,
        "variant_mafs": geno.mean(axis=0) / 2,
        "phenotype_vector": y,
        "covariate_matrix": design[:, 1:],
    }
    gene = {
        "GENE": "le10",
        "proband_count": int(y.sum()),
        "control_count": int(y.size - y.sum()),
        "n_qualifying_variants": geno.shape[1],
        "_genotype_matrix_builder": lambda: dict(built),
    }

    def run(method, robust):
        config = AssociationConfig(
            permutations=200, permute_below_p=1e-5, skat_method=method, skat_robust=robust
        )
        return AssociationEngine.from_names(["skat"], config).run_all([dict(gene)]).iloc[0]

    skat, rskat, rskato = run("SKAT", False), run("SKAT", True), run("SKATO", True)
    assert skat["skat_pvalue"] < 1e-5 <= rskat["skat_pvalue"]  # repro shape
    assert skat["empirical_p"] > 0
    assert rskat["empirical_p"] == rskato["empirical_p"] == skat["empirical_p"]
