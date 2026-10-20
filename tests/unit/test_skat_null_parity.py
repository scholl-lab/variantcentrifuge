"""Python SKAT null model / score tests vs R SKAT 2.2.5 (#179, #180, #182, #186).

Fixtures and references come from ``scripts/generate_skat_null_parity.R`` (R 4.5.3,
SKAT 2.2.5; provenance line at the top of ``r_reference.csv``). R used identical
CSV inputs and explicit weights with ``SKAT_Null_Model(..., Adjustment=FALSE)``,
``method="davies"`` (SKAT ``r.corr=0``, Burden ``r.corr=1``) and
``method="optimal.adj"`` (SKAT-O).

Tolerances: SKAT/Burden ``rel=1e-6`` -- Davies is run at ``acc=1e-6`` in both
implementations and both GLM fits stop at ~1e-8 deviance change (observed max
3e-8). Quantitative SKAT-O uses the repository's SKAT-O convention ``rel=1e-3``;
binary SKAT-O integration accuracy is tracked separately (#181).
"""

from __future__ import annotations

import csv
import pathlib

import numpy as np
import pytest
import scipy.stats

from variantcentrifuge.association.backends.null_model import full_rank_design
from variantcentrifuge.association.backends.python_backend import PythonSKATBackend
from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.tests.linear_burden import LinearBurdenTest
from variantcentrifuge.association.tests.skat_python import PurePythonSKATTest

FIXTURES = pathlib.Path(__file__).parent.parent / "fixtures" / "skat_null_parity"


def _references() -> dict[tuple[str, str], float]:
    with open(FIXTURES / "r_reference.csv") as fh:
        rows = csv.DictReader(line for line in fh if not line.startswith("#"))
        return {(r["name"], r["method"]): float(r["p_value"]) for r in rows}


R_REF = _references()
SPECS = sorted({name for name, _method in R_REF})


def _load(name: str) -> tuple[np.ndarray, np.ndarray | None, np.ndarray, np.ndarray]:
    data = np.genfromtxt(FIXTURES / f"{name}.csv", delimiter=",", names=True)
    cols = data.dtype.names
    assert cols is not None
    geno = np.column_stack([data[c] for c in cols if c.startswith("g")])
    xcols = [c for c in cols if c.startswith("x")]
    cov = np.column_stack([data[c] for c in xcols]) if xcols else None
    with open(FIXTURES / "weights.csv") as fh:
        weights = np.array([float(r["w"]) for r in csv.DictReader(fh) if r["name"] == name])
    return data["y"], cov, geno, weights


@pytest.fixture(scope="module")
def backend() -> PythonSKATBackend:
    b = PythonSKATBackend()
    b.detect_environment()
    return b


def _p(backend, y, cov, geno, w, trait, method) -> float:
    null = backend.fit_null_model(y, cov, trait)
    p = backend.test_gene("G", geno, null, method, weights=w)["p_value"]
    assert p is not None
    return float(p)


@pytest.mark.parametrize("method", ["SKAT", "Burden"])
@pytest.mark.parametrize("name", SPECS)
def test_skat_and_burden_match_r(backend, name, method):
    y, cov, geno, w = _load(name)
    p = _p(backend, y, cov, geno, w, name.split("_")[0], method)
    assert p == pytest.approx(R_REF[(name, method)], rel=1e-6, abs=0)


@pytest.mark.parametrize("name", [s for s in SPECS if s.startswith("quantitative")])
def test_quantitative_skato_matches_r(backend, name):
    y, cov, geno, w = _load(name)
    p = _p(backend, y, cov, geno, w, "quantitative", "SKATO")
    assert p == pytest.approx(R_REF[(name, "SKATO")], rel=1e-3, abs=0)


@pytest.mark.parametrize("a,b", [(0.1, 0.0), (10.0, -3.0), (-2.5, 100.0), (1e4, 0.0)])
@pytest.mark.parametrize("name", ["quantitative_n60_cov3_sig1", "quantitative_n300_cov0_sig1"])
def test_quantitative_affine_invariance(backend, name, a, b):
    y, cov, geno, w = _load(name)
    for method in ("SKAT", "Burden", "SKATO"):
        base = _p(backend, y, cov, geno, w, "quantitative", method)
        moved = _p(backend, a * y + b, cov, geno, w, "quantitative", method)
        assert moved == pytest.approx(base, rel=1e-9, abs=0), method

    data = {"genotype_matrix": geno, "variant_mafs": geno.mean(axis=0) / 2.0}
    if cov is not None:
        data["covariate_matrix"] = cov
    cfg = AssociationConfig(trait_type="quantitative")
    lb = LinearBurdenTest()
    p0 = lb.run("G", {**data, "phenotype_vector": y}, cfg).p_value
    p1 = lb.run("G", {**data, "phenotype_vector": a * y + b}, cfg).p_value
    assert p0 is not None and p1 == pytest.approx(p0, rel=1e-9, abs=0)


def _redundant_designs(cov: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "duplicate": np.column_stack([cov, cov[:, 0]]),
        "constant": np.column_stack([cov, np.full(len(cov), 7.0)]),
        "collinear": np.column_stack([cov, cov[:, 0] + 2.0 * cov[:, 1]]),
        "zero": np.column_stack([np.zeros(len(cov)), cov]),
    }


@pytest.mark.parametrize("trait", ["binary", "quantitative"])
def test_redundant_covariates_match_reduced_design(backend, trait):
    # R SKAT 2.2.5 replaces a rank-deficient X1 by svd(X1)$u (including null-space
    # columns), so its p-values change with redundancy; we use the reduced design.
    y, cov, geno, w = _load(f"{trait}_n300_cov3_sig1")
    assert cov is not None
    for method in ("SKAT", "Burden", "SKATO"):
        reduced = _p(backend, y, cov, geno, w, trait, method)
        assert reduced == pytest.approx(R_REF[(f"{trait}_n300_cov3_sig1", method)], rel=1e-3, abs=0)
        for label, design in _redundant_designs(cov).items():
            p = _p(backend, y, design, geno, w, trait, method)
            assert p == pytest.approx(reduced, rel=1e-9, abs=0), (method, label)


def test_full_rank_design_drops_dependent_columns_and_logs(caplog):
    x = np.random.default_rng(0).normal(size=(50, 2))
    design, r_sel = full_rank_design(np.column_stack([x, x[:, 1] * 1e6, np.ones(50)]), 50)
    assert design.shape == (50, 3)
    np.testing.assert_allclose(design[:, 1:].T @ design[:, 1:], np.eye(2), atol=1e-12)
    np.testing.assert_allclose(design[:, 1:].sum(axis=0), 0.0, atol=1e-12)
    assert 1.0 <= np.linalg.cond(r_sel) < 10.0
    assert "dropped 2 covariate column(s)" in caplog.text
    # Covariates in tiny units are not mistaken for dependence
    assert full_rank_design(x * 1e-9, 50)[0].shape == (50, 3)
    assert full_rank_design(None, 50)[0].shape == (50, 1)
    with pytest.raises(ValueError, match="finite"):
        full_rank_design(np.full((50, 1), np.nan), 50)


def test_duplicate_covariate_does_not_crash_public_api():
    rng = np.random.default_rng(7654)
    geno = rng.binomial(2, np.linspace(0.01, 0.15, 9), (300, 9)).astype(float)
    x = rng.normal(size=(300, 1))
    y = np.zeros(300)
    y[:10] = 1.0
    rng.shuffle(y)
    test = PurePythonSKATTest()
    test.check_dependencies()
    data = {
        "genotype_matrix": geno,
        "phenotype_vector": y,
        "n_qualifying_variants": 9,
        "proband_count": 10,
        "control_count": 290,
    }
    for method in ("SKAT", "SKATO"):
        cfg = AssociationConfig(skat_method=method)
        dup = test.run("G", {**data, "covariate_matrix": np.column_stack([x, x])}, cfg)
        single = test.run("G", {**data, "covariate_matrix": x}, cfg)
        assert dup.p_value is not None
        assert dup.p_value == pytest.approx(single.p_value, rel=1e-9, abs=0)


DEGENERATE = {
    "constant": (np.full(60, 5.0), "CONSTANT_PHENOTYPE"),
    "all_zero": (np.zeros(60), "CONSTANT_PHENOTYPE"),
    "all_nan": (np.full(60, np.nan), "NONFINITE_PHENOTYPE"),
    "with_inf": (np.r_[np.inf, np.arange(59.0)], "NONFINITE_PHENOTYPE"),
    "with_nan": (np.r_[np.nan, np.arange(59.0)], "NONFINITE_PHENOTYPE"),
}


@pytest.mark.parametrize("case", list(DEGENERATE))
def test_degenerate_outcomes_rejected(backend, case):
    y, code = DEGENERATE[case]
    rng = np.random.default_rng(331)
    geno = rng.binomial(2, np.linspace(0.03, 0.22, 12), size=(60, 12)).astype(float)
    cov = rng.normal(size=(60, 3))
    for trait in ("binary", "quantitative"):
        with pytest.raises(ValueError, match=code):
            backend.fit_null_model(y, cov, trait)

    data = {
        "genotype_matrix": geno,
        "variant_mafs": geno.mean(axis=0) / 2.0,
        "phenotype_vector": y,
        "covariate_matrix": cov,
        "n_qualifying_variants": 12,
    }
    skat = PurePythonSKATTest()
    skat.check_dependencies()
    res = skat.run("G", data, AssociationConfig(trait_type="quantitative"))
    assert res.p_value is None and res.extra["skat_warnings"] == code

    res = LinearBurdenTest().run("G", data, AssociationConfig(trait_type="quantitative"))
    assert res.p_value is None and res.effect_size is None and res.se is None
    assert res.extra["warnings"] == [code]


def test_null_model_requires_residual_degrees_of_freedom(backend):
    y = np.array([0.0, 1.0, 2.0])
    with pytest.raises(ValueError, match="INSUFFICIENT_RESIDUAL_DF"):
        backend.fit_null_model(y, np.array([[1.0, 0.0], [0.0, 1.0], [2.0, 5.0]]), "quantitative")


def test_linear_burden_rejects_nonfinite_fit_statistics(monkeypatch):
    import statsmodels.api as sm

    class _Fit:
        params = np.array([0.0, np.nan])
        bse = np.array([1.0, np.nan])
        pvalues = np.array([1.0, np.nan])

        def conf_int(self):
            return np.array([[0.0, 0.0], [np.nan, np.nan]])

    class _OLS:
        def __init__(self, *args):
            pass

        def fit(self):
            return _Fit()

    monkeypatch.setattr(sm, "OLS", _OLS)
    geno = np.eye(4)[:, :2]
    res = LinearBurdenTest().run(
        "G",
        {
            "genotype_matrix": geno,
            "variant_mafs": np.full(2, 0.1),
            "phenotype_vector": np.arange(4.0),
        },
        AssociationConfig(trait_type="quantitative"),
    )
    assert res.p_value is None
    assert res.extra["warnings"] == ["NONFINITE_STATISTICS"]


def test_variants_explained_by_covariates_have_no_kernel(backend):
    y, cov, _geno, _w = _load("quantitative_n60_cov3_sig1")
    assert cov is not None
    geno = np.column_stack([np.ones(len(y)), cov[:, 0]])
    null = backend.fit_null_model(y, cov, "quantitative")
    for method in ("SKAT", "Burden", "SKATO"):
        res = backend.test_gene("G", geno, null, method, weights=np.ones(2))
        assert res["p_value"] is None and res["skip_reason"] == "no_informative_variants"
        assert res["n_variants"] == 2 and res["n_marker_test"] == 0

    # Neither variant is explained alone, but their weighted sum is a covariate:
    # the burden lies in the nuisance space (round-off floor), SKAT still tests.
    rng = np.random.default_rng(3)
    pair = rng.binomial(2, 0.2, (len(y), 2)).astype(float)
    null = backend.fit_null_model(y, np.column_stack([cov, pair.sum(axis=1)]), "quantitative")
    burden = backend.test_gene("G", pair, null, "Burden", weights=np.ones(2))
    assert burden["p_value"] is None and burden["skip_reason"] == "zero_variance"
    skat = backend.test_gene("G", pair, null, "SKAT", weights=np.ones(2))
    assert skat["p_value"] is not None and skat["n_marker_test"] == 2


def _confounded():
    """Review B1-r1 design: covariate x strongly confounds variant 0 and y."""
    rng = np.random.default_rng(123)
    geno = rng.binomial(2, [0.1, 0.2], (300, 2)).astype(float)
    x = geno[:, 0] + 0.1 * rng.normal(size=300)
    y = 2 * x + rng.normal(size=300)
    return rng, geno, x, y


@pytest.mark.parametrize("trait", ["binary", "quantitative"])
def test_covariate_shift_and_scale_invariance(backend, trait):
    rng, geno, x, y = _confounded()
    if trait == "binary":
        y = (y > np.median(y)).astype(float)
    z = rng.normal(size=300)
    w = np.ones(2)
    for method in ("SKAT", "Burden", "SKATO"):
        base = _p(backend, y, x[:, None], geno, w, trait, method)
        # x + 1e7 loses ~1e7 * eps / sd(x) ~ 1e-8 relative input precision
        assert _p(backend, y, (x + 1e7)[:, None], geno, w, trait, method) == pytest.approx(
            base, rel=1e-6, abs=0
        )
        for cov in (x[:, None] * 1e8, x[:, None] * 1e-8, np.column_stack([x, np.full(300, 1e16)])):
            assert _p(backend, y, cov, geno, w, trait, method) == pytest.approx(
                base, rel=1e-9, abs=0
            )
        with_z = _p(backend, y, np.column_stack([x, z]), geno, w, trait, method)
        tiny_z = _p(backend, y, np.column_stack([x, 1e-16 * z]), geno, w, trait, method)
        assert tiny_z == pytest.approx(with_z, rel=1e-9, abs=0), method


def test_small_projected_burden_is_tested(backend):
    _rng, geno, _x, y = _confounded()
    cov = geno[:, :1]
    w = np.array([1.0, 1e-8])
    p = _p(backend, y, cov, geno, w, "quantitative", "Burden")
    # Independent OLS score test with explicit hat matrix
    design = np.column_stack([np.ones(300), cov])
    hat = design @ np.linalg.solve(design.T @ design, design.T)
    resid = y - hat @ y
    burden = geno @ w
    pb = burden - hat @ burden
    z = (resid @ burden) / np.sqrt(resid @ resid / 298 * (pb @ pb))
    # |Pb| ~ 1e-8 |b|: cancellation leaves ~8 significant digits in either computation
    assert p == pytest.approx(2 * scipy.stats.norm.sf(abs(z)), rel=1e-5, abs=0)


def test_exact_translation_of_covariate_is_bit_invariant(backend):
    """Review B1-r2 item 1: +1e13 on a 1/16-grid covariate is exactly representable."""
    _rng, geno, x, _y = _confounded()
    x = np.round(x * 16) / 16
    y = 2 * x + np.random.default_rng(42).normal(size=300)
    for method in ("SKAT", "Burden", "SKATO"):
        base = _p(backend, y, x[:, None], geno, np.ones(2), "quantitative", method)
        moved = _p(backend, y, (x + 1e13)[:, None], geno, np.ones(2), "quantitative", method)
        assert moved == base, method


def test_near_collinear_covariates_explain_variants_in_all_consumers(backend):
    """Review B1-r2 item 3: both variants lie in span(1, g0, g0 + 2^-21 g1)."""
    from variantcentrifuge.association.meta.artifact import compute_gene_score
    from variantcentrifuge.association.meta.combine import meta_burden, meta_skat

    rng = np.random.default_rng(123)
    geno = rng.binomial(2, [0.1, 0.2], (300, 2)).astype(float)
    y = rng.binomial(1, 0.4, 300).astype(float)
    cov = np.column_stack([geno[:, 0], geno[:, 0] + 2.0**-21 * geno[:, 1]])
    null = backend.fit_null_model(y, cov, "binary")
    w = np.ones(2)
    for method in ("SKAT", "Burden", "SKATO"):
        res = backend.test_gene("G", geno, null, method, weights=w)
        assert res["p_value"] is None and res["skip_reason"] == "no_informative_variants"
    score = compute_gene_score(geno, geno, null, ["1:1:A:G", "1:2:A:G"])
    assert not score.u.any() and not score.v.any()
    assert meta_skat([score.u], [score.v], w)[0] == 1.0
    assert meta_burden([score.u], [score.v], w)[0] is None

    from variantcentrifuge.association.backends.coast_python import PythonCOASTBackend
    from variantcentrifuge.association.tests.acat import compute_acat_v

    assert PythonCOASTBackend()._run_allelic_skat(geno, w, null) == 1.0
    acat = compute_acat_v(
        geno, null.extra["residuals"], "binary", 1.0, null.extra["mu_hat"], w, null_model=null
    )
    assert acat is None


@pytest.mark.parametrize("trait", ["binary", "quantitative"])
def test_explained_variant_cannot_mask_small_independent_variant(backend, trait):
    """Review B1-r2 item 4: weights [1, 1e-10] must equal [0, 1e-10] (g0 is a covariate)."""
    rng = np.random.default_rng(5)
    geno = rng.binomial(2, [0.2, 0.001], (100000, 2)).astype(float)
    eta = 0.5 * geno[:, 1] + rng.normal(size=100000)
    y = eta if trait == "quantitative" else (eta > 0.5).astype(float)
    null = backend.fit_null_model(y, geno[:, :1], trait)
    big = backend.test_gene("G", geno, null, "Burden", weights=np.array([1.0, 1e-10]))
    ref = backend.test_gene("G", geno, null, "Burden", weights=np.array([0.0, 1e-10]))
    assert ref["p_value"] is not None and ref["p_value"] < 1e-3
    assert big["p_value"] == pytest.approx(ref["p_value"], rel=1e-9, abs=0)
    # One retained variant: SKAT follows the single-variant (rank < 2) policy either way
    big = backend.test_gene("G", geno, null, "SKAT", weights=np.array([1.0, 1e-10]))
    ref = backend.test_gene("G", geno, null, "SKAT", weights=np.array([0.0, 1e-10]))
    assert big["skip_reason"] == ref["skip_reason"] == "rank_deficient"


@pytest.mark.parametrize("name", [s for s in SPECS if s.startswith("quantitative")])
def test_linear_burden_matches_r_lm(name):
    y, cov, geno, _w = _load(name)
    data = {"genotype_matrix": geno, "variant_mafs": geno.mean(axis=0) / 2.0, "phenotype_vector": y}
    if cov is not None:
        data["covariate_matrix"] = cov
    res = LinearBurdenTest().run("G", data, AssociationConfig(trait_type="quantitative"))
    # Both are double-precision least squares on a well-conditioned design
    assert res.p_value == pytest.approx(R_REF[(name, "LinearBurden")], rel=1e-9, abs=0)
    assert res.effect_size == pytest.approx(R_REF[(name, "LinearBurden_beta")], rel=1e-9, abs=0)
    assert res.se == pytest.approx(R_REF[(name, "LinearBurden_se")], rel=1e-9, abs=0)


@pytest.mark.parametrize("case", ["constant_with_covariate", "constant_no_covariate", "collinear"])
def test_linear_burden_rejects_non_estimable_burden(case):
    rng = np.random.default_rng(1)
    age = rng.normal(50, 10, 200)
    y = 0.1 * age + rng.normal(size=200)
    carrier = rng.binomial(2, 0.3, (200, 1)).astype(float)
    geno, cov = {
        "constant_with_covariate": (np.ones((200, 1)), age[:, None]),
        "constant_no_covariate": (np.ones((200, 1)), None),
        "collinear": (carrier, np.column_stack([age, 3.0 * carrier[:, 0] + 1.0])),
    }[case]
    data = {"genotype_matrix": geno, "variant_mafs": np.array([0.3]), "phenotype_vector": y}
    if cov is not None:
        data["covariate_matrix"] = cov
    res = LinearBurdenTest().run("G", data, AssociationConfig(trait_type="quantitative"))
    assert res.p_value is None and res.effect_size is None
    assert res.extra["warnings"] == ["BURDEN_NOT_ESTIMABLE"]


def test_linear_burden_rejects_nonfinite_covariates():
    geno = np.eye(6)[:, :2]
    data = {
        "genotype_matrix": geno,
        "variant_mafs": np.full(2, 0.1),
        "phenotype_vector": np.arange(6.0),
        "covariate_matrix": np.r_[np.nan, np.arange(5.0)][:, None],
    }
    res = LinearBurdenTest().run("G", data, AssociationConfig(trait_type="quantitative"))
    assert res.p_value is None and res.extra["warnings"] == ["NONFINITE_COVARIATES"]


def _near_dependent_rare(trait: str):
    """Review B1-r3 item 1: n=200k, rare variants; [x, z] vs [x, x + 2^-23 z]."""
    rng = np.random.default_rng(42)
    n = 200000
    x = rng.choice([-1.0, 1.0], n)
    z = rng.choice([-1.0, 1.0], n)
    geno = rng.binomial(2, [0.001, 0.002], (n, 2)).astype(float)
    if trait == "binary":
        pr = np.where(x > 0, 0.1, 0.0001)
        odds = pr / (1 - pr) * 5 ** geno[:, 0]
        y = rng.binomial(1, odds / (1 + odds)).astype(float)
    else:
        y = x + 0.5 * geno[:, 0] + rng.normal(size=n)
    return geno, x, z, y


@pytest.mark.parametrize("trait", ["binary", "quantitative"])
def test_near_dependent_equivalent_design_is_never_silently_lost(backend, trait):
    """Review B1-r3 item 1: [x, x + 2^-23 z] must not zero informative variants.

    ponytail: with the conservative round-off bound (kappa ~ 1.7e7, n = 200k) the
    projection is indeterminate, so the gene is skipped explicitly instead of
    matching the [x, z] result; a certified tighter bound would recover the test.
    """
    from variantcentrifuge.association.backends.null_model import UntestableGeneError
    from variantcentrifuge.association.meta.artifact import compute_gene_score

    geno, x, z, y = _near_dependent_rare(trait)
    w = np.ones(2)
    keys = ["1:1:A:G", "1:2:A:G"]
    plain = backend.fit_null_model(y, np.column_stack([x, z]), trait)
    near = backend.fit_null_model(y, np.column_stack([x, x + 2.0**-23 * z]), trait)
    for method in ("SKAT", "Burden"):
        ref = backend.test_gene("G", geno, plain, method, weights=w)
        assert ref["p_value"] < 1e-6 and ref["n_marker_test"] == 2
        res = backend.test_gene("G", geno, near, method, weights=w)
        if res["p_value"] is None:
            assert res["skip_reason"] == "ill_conditioned_covariates"
        else:
            assert res["p_value"] == pytest.approx(ref["p_value"], rel=1e-6, abs=0)
    if trait == "binary":
        assert np.all(np.abs(compute_gene_score(geno, geno, plain, keys).u) > 1.0)
        with pytest.raises(UntestableGeneError, match="ILL_CONDITIONED_COVARIATES"):
            compute_gene_score(geno, geno, near, keys)


def test_indeterminate_projection_skips_instead_of_zeroing(backend, monkeypatch):
    """When the round-off bound exceeds 1e-6 nothing is zeroed; the gene is skipped."""
    from variantcentrifuge.association.backends import null_model as nm_mod
    from variantcentrifuge.association.tests.acat import compute_acat_v

    y, cov, geno, w = _load("quantitative_n60_cov3_sig1")
    null = backend.fit_null_model(y, cov, "quantitative")
    # near-parallel covariates (cond ~1e9 after normalization)
    null.extra["covariate_r"] = np.array([[1.0, 1.0, 0.0], [0.0, 1e-9, 0.0], [0.0, 0.0, 1.0]])
    assert nm_mod.explained_columns(geno, null) is None
    for method in ("SKAT", "Burden", "SKATO"):
        res = backend.test_gene("G", geno, null, method, weights=w)
        assert res["p_value"] is None and res["skip_reason"] == "ill_conditioned_covariates"
    resid = null.extra["residuals"]
    assert compute_acat_v(geno, resid, "quantitative", 1.0, None, w, null_model=null) is None

    from variantcentrifuge.association.backends.coast_python import PythonCOASTBackend
    from variantcentrifuge.association.backends.null_model import UntestableGeneError
    from variantcentrifuge.association.meta.artifact import compute_gene_score

    with pytest.raises(UntestableGeneError):
        PythonCOASTBackend()._run_allelic_skat(geno, w, null)
    with pytest.raises(UntestableGeneError):  # never exported or combined
        compute_gene_score(geno, geno, null, [f"1:{j}:A:G" for j in range(geno.shape[1])])

    monkeypatch.setattr(
        "variantcentrifuge.association.tests.linear_burden.explained_columns", lambda *a: None
    )
    data = {"genotype_matrix": geno, "variant_mafs": geno.mean(axis=0) / 2, "phenotype_vector": y}
    res = LinearBurdenTest().run("G", data, AssociationConfig(trait_type="quantitative"))
    assert res.p_value is None and res.extra["warnings"] == ["ILL_CONDITIONED_COVARIATES"]


@pytest.mark.parametrize("extra", ["duplicate", "zero"])
def test_acat_v_raw_design_is_rank_reduced(backend, extra):
    """Review B1-r3 item 4: redundant raw design columns give the full-rank result."""
    from variantcentrifuge.association.tests.acat import compute_acat_v

    rng = np.random.default_rng(123)
    x = rng.normal(size=300)
    geno = rng.binomial(2, 0.05, (300, 3)).astype(float)
    y = rng.binomial(1, 0.4, 300).astype(float)
    null = backend.fit_null_model(y, x[:, None], "binary")
    args = (geno, null.extra["residuals"], "binary", 1.0, null.extra["mu_hat"], np.ones(3))
    base = compute_acat_v(*args, design_matrix=np.column_stack([np.ones(300), x]))
    col = x if extra == "duplicate" else np.zeros(300)
    redundant = compute_acat_v(*args, design_matrix=np.column_stack([np.ones(300), x, col]))
    assert base is not None
    assert redundant == pytest.approx(base, rel=1e-12, abs=0)
    assert compute_acat_v(*args, null_model=null) == pytest.approx(base, rel=1e-9, abs=0)


def test_linear_burden_near_dependent_covariates_use_reduced_basis():
    """Review B1-r3 item 3: [x, x + 2^-25 g] drops the second covariate (lm tol 1e-7)."""
    rng = np.random.default_rng(123)
    g = rng.binomial(2, 0.1, (300, 1)).astype(float)
    x = rng.normal(size=300)
    y = g[:, 0] + rng.normal(size=300)
    cfg = AssociationConfig(trait_type="quantitative")
    base = {"genotype_matrix": g, "variant_mafs": np.array([0.1]), "phenotype_vector": y}
    near = LinearBurdenTest().run(
        "G", {**base, "covariate_matrix": np.column_stack([x, x + 2.0**-25 * g[:, 0]])}, cfg
    )
    only_x = LinearBurdenTest().run("G", {**base, "covariate_matrix": x[:, None]}, cfg)
    assert near.p_value == pytest.approx(only_x.p_value, rel=1e-12, abs=0)
    assert near.effect_size == pytest.approx(only_x.effect_size, rel=1e-12, abs=0)
    exact = LinearBurdenTest().run(
        "G", {**base, "covariate_matrix": np.column_stack([x, g[:, 0]])}, cfg
    )
    assert exact.p_value is None and exact.extra["warnings"] == ["BURDEN_NOT_ESTIMABLE"]


def test_backend_input_guards(backend):
    y, cov, geno, _w = _load("quantitative_n60_cov3_sig1")
    null = backend.fit_null_model(y, cov, "quantitative")
    with pytest.raises(ValueError, match="concrete variant weights"):
        backend.test_gene("G", geno, null, "SKAT", weights=None)
    with pytest.raises(ValueError, match="finite"):
        backend.test_gene("G", geno, null, "SKAT", weights=np.full(geno.shape[1], np.nan))
    empty = np.empty((len(y), 0))
    for method in ("SKAT", "SKATO"):
        res = backend.test_gene("G", empty, null, method, weights=np.empty(0))
        assert res["p_value"] is None
    # 1-D covariate input is ignored by linear burden (no 2-D covariate matrix)
    data = {
        "genotype_matrix": geno,
        "variant_mafs": geno.mean(axis=0) / 2,
        "phenotype_vector": y,
        "covariate_matrix": cov[:, 0],
    }
    plain = {k: v for k, v in data.items() if k != "covariate_matrix"}
    cfg = AssociationConfig(trait_type="quantitative")
    assert LinearBurdenTest().run("G", data, cfg).p_value == pytest.approx(
        LinearBurdenTest().run("G", plain, cfg).p_value, rel=1e-12, abs=0
    )


def test_meta_score_of_indeterminate_gene_is_not_exported(backend):
    """Review B1-r4 item 1: all variants in span of [g0, g0 + 2^-16 g_j]; meta must not give p=0."""
    from variantcentrifuge.association.backends.null_model import UntestableGeneError
    from variantcentrifuge.association.meta.artifact import compute_gene_score

    rng = np.random.default_rng(0)
    geno = rng.binomial(2, np.linspace(0.02, 0.3, 8), (300, 8)).astype(float)
    y = rng.binomial(1, 0.05, 300).astype(float)
    cov = geno.copy()
    cov[:, 1:] = geno[:, :1] + 2.0**-16 * geno[:, 1:]
    null = backend.fit_null_model(y, cov, "binary")
    keys = [f"1:{j}:A:G" for j in range(8)]
    try:
        score = compute_gene_score(geno, geno, null, keys)
    except UntestableGeneError as exc:
        assert exc.reason == "ILL_CONDITIONED_COVARIATES"
    else:  # classified as explained: nothing left to combine
        assert not score.u.any() and not score.v.any()


@pytest.mark.parametrize("trait", ["binary", "quantitative"])
def test_burden_in_nuisance_span_is_not_tested(backend, trait):
    """Review B1-r4 item 5: each variant informative, their sum exactly in the nuisance span."""
    rng = np.random.default_rng(25)
    geno = rng.binomial(2, [0.1, 0.2], (10000, 2)).astype(float)
    a = rng.binomial(2, 0.2, (10000, 3)).astype(float)
    y = rng.binomial(1, 0.05, 10000).astype(float) if trait == "binary" else rng.normal(size=10000)
    cov = a.copy()
    cov[:, 1:] = a[:, :1] + 2.0**-16 * a[:, 1:]
    cov = np.column_stack([cov, a[:, 0] + 2.0**-16 * geno.sum(axis=1)])
    null = backend.fit_null_model(y, cov, trait)
    res = backend.test_gene("G", geno, null, "Burden", weights=np.ones(2))
    assert res["p_value"] is None
    assert res["skip_reason"] in ("zero_variance", "ill_conditioned_covariates")


def test_zero_weight_variant_cannot_veto_gene(backend):
    """Review B1-r4 item 6: weights [0, 1, 1] equal deleting the zero-weight column."""
    from variantcentrifuge.association.tests.acat import compute_acat_v

    rng = np.random.default_rng(123)
    geno = rng.binomial(2, 0.1, (300, 3)).astype(float)
    x = geno[:, 0] + 1e-8 * rng.normal(size=300)
    y = geno[:, 1] + rng.normal(size=300)
    null = backend.fit_null_model(y, x[:, None], "quantitative")
    for method in ("SKAT", "Burden", "SKATO"):
        full = backend.test_gene("G", geno, null, method, weights=np.array([0.0, 1.0, 1.0]))
        sub = backend.test_gene("G", geno[:, 1:], null, method, weights=np.ones(2))
        assert full["p_value"] == pytest.approx(sub["p_value"], rel=1e-12, abs=0), method
        assert full["n_variants"] == 3 and full["n_marker_test"] == 2
    resid, s2 = null.extra["residuals"], null.extra["sigma2"]
    full_acat = compute_acat_v(
        geno, resid, "quantitative", s2, None, np.array([0.0, 1, 1]), null_model=null
    )
    sub_acat = compute_acat_v(
        geno[:, 1:], resid, "quantitative", s2, None, np.ones(2), null_model=null
    )
    assert full_acat == pytest.approx(sub_acat, rel=1e-12, abs=0)


# Seeded subset of the in-span scans (1,794 constructions in two runs, 0 exceedances,
# max ratio/bound 0.054): x = [g0, g0 + sep * g_j] puts every g_j exactly in span.
@pytest.mark.parametrize("trait", ["binary", "quantitative"])
@pytest.mark.parametrize("sep_exp", [8, 16, 23, 30])
@pytest.mark.parametrize("n,k", [(300, 4), (2000, 2), (2000, 8)])
def test_in_span_columns_never_informative_and_within_bound(trait, sep_exp, n, k):
    from variantcentrifuge.association.backends.null_model import (
        _project,
        explained_columns,
        fit_glm_null_model,
    )

    rng = np.random.default_rng(1000 * sep_exp + n + k)
    geno = rng.binomial(2, np.linspace(0.02, 0.3, k), (n, k)).astype(float)
    cov = geno.copy()
    cov[:, 1:] = geno[:, :1] + 2.0**-sep_exp * geno[:, 1:]
    y = rng.binomial(1, 0.3, n).astype(float) if trait == "binary" else rng.normal(size=n)
    null = fit_glm_null_model(y, cov, trait)
    if null.extra["design_matrix"].shape[1] != k + 1:
        pytest.skip("covariates rank-reduced: genotypes no longer exactly in span")
    adjusted, z, bound = _project(geno, null)
    ratio = np.linalg.norm(adjusted, axis=0) / np.linalg.norm(z, axis=0)
    assert np.all(ratio <= bound)
    explained = explained_columns(geno, null)
    assert explained is None or explained.all()


@pytest.mark.parametrize("seed", [0, 6, 17, 31])
def test_skato_with_one_retained_variant_follows_single_variant_policy(backend, seed):
    """Final review item 1: covariate explains g0; one informative variant remains."""
    rng = np.random.default_rng(seed)
    geno = rng.binomial(2, [0.1, 0.2], (300, 2)).astype(float)
    y = rng.normal(size=300)
    null = backend.fit_null_model(y, geno[:, :1], "quantitative")
    single = backend.test_gene("G", geno[:, 1:], null, "SKATO", weights=np.ones(1))
    for method in ("SKAT", "SKATO"):
        res = backend.test_gene("G", geno, null, method, weights=np.ones(2))
        assert res["p_value"] is None and res["skip_reason"] == single["skip_reason"]
    burden = backend.test_gene("G", geno, null, "Burden", weights=np.ones(2))
    ref = backend.test_gene("G", geno[:, 1:], null, "Burden", weights=np.ones(1))
    assert burden["p_value"] == pytest.approx(ref["p_value"], rel=1e-12, abs=0)


@pytest.mark.parametrize("scale", [1.0, 1e-6, 1e6])
def test_outcome_explained_by_covariates_is_untestable(backend, scale):
    """Final review item 2: y = x exactly; no residual variation (unit invariant)."""
    from variantcentrifuge.association.backends.null_model import UntestableGeneError

    rng = np.random.default_rng(0)
    geno = rng.binomial(2, 0.1, (300, 8)).astype(float)
    x = rng.integers(0, 10, 300).astype(float)[:, None] * scale
    y = x[:, 0]
    with pytest.raises(UntestableGeneError, match="OUTCOME_EXPLAINED_BY_COVARIATES"):
        backend.fit_null_model(y, x, "quantitative")
    data = {
        "genotype_matrix": geno,
        "variant_mafs": geno.mean(axis=0) / 2,
        "phenotype_vector": y,
        "covariate_matrix": x,
        "n_qualifying_variants": 8,
    }
    cfg = AssociationConfig(trait_type="quantitative")
    res = LinearBurdenTest().run("G", data, cfg)
    assert res.p_value is None and res.extra["warnings"] == ["OUTCOME_EXPLAINED_BY_COVARIATES"]
    skat = PurePythonSKATTest()
    skat.check_dependencies()
    for _ in range(2):  # second call uses the cached failure
        res = skat.run("G", data, cfg)
        assert (
            res.p_value is None and res.extra["skat_warnings"] == "OUTCOME_EXPLAINED_BY_COVARIATES"
        )
    # A genuine residual still tests normally
    data["phenotype_vector"] = y + scale * rng.normal(size=300)
    assert skat.run("G", data, cfg).p_value is not None


def test_binary_outcome_separated_by_covariates_is_untestable(backend):
    from variantcentrifuge.association.backends.null_model import UntestableGeneError

    x = np.arange(100, dtype=float)[:, None]
    y = (x[:, 0] >= 50).astype(float)
    with pytest.raises(UntestableGeneError, match="OUTCOME_EXPLAINED_BY_COVARIATES"):
        backend.fit_null_model(y, x, "binary")


def test_insufficient_residual_df_skips_skat_gene():
    rng = np.random.default_rng(1)
    geno = rng.binomial(2, 0.2, (6, 3)).astype(float)
    data = {
        "genotype_matrix": geno,
        "phenotype_vector": np.array([0.0, 1, 0, 1, 0, 1]),
        "covariate_matrix": rng.normal(size=(6, 5)),
        "n_qualifying_variants": 3,
    }
    skat = PurePythonSKATTest()
    skat.check_dependencies()
    res = skat.run("G", data, AssociationConfig(skat_method="SKAT"))
    assert res.p_value is None and res.extra["skat_warnings"] == "INSUFFICIENT_RESIDUAL_DF"


def test_meta_burden_skips_cancelled_collapsed_variance(backend):
    """Final review item 4: both variants informative, their sum is the covariate."""
    from variantcentrifuge.association.meta.artifact import compute_gene_score
    from variantcentrifuge.association.meta.combine import meta_burden

    rng = np.random.default_rng(21)
    geno = rng.binomial(2, 0.1, (60, 2)).astype(float)
    x = geno.sum(axis=1)[:, None]
    y = rng.binomial(1, 1 / (1 + np.exp(4 - 3 * x[:, 0]))).astype(float)
    null = backend.fit_null_model(y, x, "binary")
    assert backend.test_gene("G", geno, null, "Burden", weights=np.ones(2))["p_value"] is None
    score = compute_gene_score(geno, geno, null, ["1:1:A:G", "1:2:A:G"])
    assert meta_burden([score.u], [score.v], np.ones(2)) == (None, 0.0, 0)
