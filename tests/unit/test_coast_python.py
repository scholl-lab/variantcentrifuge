"""
Unit tests for PythonCOASTBackend components and PurePythonCOASTTest lifecycle.

Tests each layer of the pure Python COAST implementation independently:
- Layer 1: _aggregate_by_category() — genotype aggregation by annotation category
- Layer 2: _run_burden_test() — regression burden tests (OLS and logit LRT)
- Layer 3: _compute_allelic_skat_weights() — annotation-aware SKAT weights
- Layer 4: PythonCOASTBackend.test_gene() — full 7-component omnibus pipeline
- Layer 5: PurePythonCOASTTest lifecycle — engine integration (check_deps -> run -> finalize)

Covers requirements: COAST-PY-01 through COAST-PY-05
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.backends.coast_python import (
    PythonCOASTBackend,
    _aggregate_by_category,
    _compute_allelic_skat_weights,
    _run_burden_test,
)
from variantcentrifuge.association.tests.allelic_series_python import PurePythonCOASTTest

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _make_gene_df(
    effects: list[str],
    impacts: list[str],
    sift: list[str] | None = None,
    polyphen: list[str] | None = None,
) -> pd.DataFrame:
    """Build a per-variant DataFrame for COAST annotation testing."""
    data: dict[str, list] = {"EFFECT": effects, "IMPACT": impacts}
    if sift is not None:
        data["dbNSFP_SIFT_pred"] = sift
    if polyphen is not None:
        data["dbNSFP_Polyphen2_HDIV_pred"] = polyphen
    return pd.DataFrame(data)


def _make_3cat_gene_df(n_bmv: int = 3, n_dmv: int = 3, n_ptv: int = 3) -> pd.DataFrame:
    """Build a gene_df with all 3 COAST categories: BMV, DMV, PTV."""
    effects = (
        ["missense_variant"] * n_bmv  # BMV: missense + tolerated/benign
        + ["missense_variant"] * n_dmv  # DMV: missense + deleterious
        + ["frameshift_variant"] * n_ptv  # PTV: frameshift + HIGH
    )
    impacts = ["MODERATE"] * n_bmv + ["MODERATE"] * n_dmv + ["HIGH"] * n_ptv
    sift = ["tolerated"] * n_bmv + ["deleterious"] * n_dmv + ["."] * n_ptv
    polyphen = ["benign"] * n_bmv + ["probably_damaging"] * n_dmv + ["."] * n_ptv
    return _make_gene_df(effects, impacts, sift, polyphen)


def _make_coast_contingency_data(
    n_samples: int = 100,
    n_bmv: int = 3,
    n_dmv: int = 3,
    n_ptv: int = 3,
    seed: int = 42,
    trait_type: str = "binary",
) -> tuple[dict[str, Any], Any]:
    """
    Build contingency_data dict and AssociationConfig for PurePythonCOASTTest testing.

    Returns (contingency_data, config) matching the lifecycle test pattern.
    """
    from variantcentrifuge.association.base import AssociationConfig

    n_variants = n_bmv + n_dmv + n_ptv
    rng = np.random.default_rng(seed)

    geno = rng.choice([0, 1, 2], size=(n_samples, n_variants), p=[0.6, 0.3, 0.1]).astype(np.float64)

    if trait_type == "binary":
        n_cases = n_samples // 2
        n_controls = n_samples - n_cases
        phenotype = np.array([1.0] * n_cases + [0.0] * n_controls)
    else:
        n_cases = 0
        n_controls = n_samples
        phenotype = rng.standard_normal(n_samples)

    gene_df = _make_3cat_gene_df(n_bmv=n_bmv, n_dmv=n_dmv, n_ptv=n_ptv)

    contingency_data: dict[str, Any] = {
        "genotype_matrix": geno,
        "phenotype_vector": phenotype,
        "gene_df": gene_df,
        "proband_count": int(phenotype.sum()) if trait_type == "binary" else 0,
        "control_count": int((1 - phenotype).sum()) if trait_type == "binary" else n_samples,
        "n_qualifying_variants": n_variants,
        "covariate_matrix": None,
    }
    config = AssociationConfig(
        trait_type=trait_type,
        coast_weights=[1.0, 2.0, 3.0],
    )
    return contingency_data, config


# ---------------------------------------------------------------------------
# TestCOASTBurdenAggregation
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestCOASTBurdenAggregation:
    """Test _aggregate_by_category() for count/indicator and none/sum/max methods."""

    def _make_geno_and_anno(
        self, n_samples: int = 20, seed: int = 0
    ) -> tuple[np.ndarray, np.ndarray]:
        """Build (geno, anno_codes) with 2 BMV + 2 DMV + 2 PTV."""
        rng = np.random.default_rng(seed)
        geno = rng.choice([0, 1, 2], size=(n_samples, 6), p=[0.6, 0.3, 0.1]).astype(np.float64)
        anno_codes = np.array([1, 1, 2, 2, 3, 3])  # 2 BMV + 2 DMV + 2 PTV
        return geno, anno_codes

    def test_count_none_returns_3_columns(self):
        """method='none': output shape is (n_samples, 3) — one column per category."""
        n = 30
        geno, anno = self._make_geno_and_anno(n_samples=n)
        result = _aggregate_by_category(geno, anno, weights=[1.0, 1.0, 1.0], method="none")
        assert result.shape == (n, 3), f"Expected (n, 3), got {result.shape}"

    def test_count_none_column_sums_match_manual(self):
        """method='none': each column equals the per-sample sum for that category."""
        n = 20
        geno, anno = self._make_geno_and_anno(n_samples=n)
        result = _aggregate_by_category(geno, anno, weights=[1.0, 1.0, 1.0], method="none")
        # BMV column (index 0) = sum of variants with anno==1 (columns 0 and 1)
        expected_bmv = geno[:, 0] + geno[:, 1]
        expected_dmv = geno[:, 2] + geno[:, 3]
        expected_ptv = geno[:, 4] + geno[:, 5]
        np.testing.assert_allclose(result[:, 0], expected_bmv, err_msg="BMV column mismatch")
        np.testing.assert_allclose(result[:, 1], expected_dmv, err_msg="DMV column mismatch")
        np.testing.assert_allclose(result[:, 2], expected_ptv, err_msg="PTV column mismatch")

    def test_indicator_converts_to_binary(self):
        """indicator=True: all values are 0 or 1 (carrier indicator)."""
        geno, anno = self._make_geno_and_anno(n_samples=30)
        result = _aggregate_by_category(
            geno, anno, weights=[1.0, 1.0, 1.0], indicator=True, method="none"
        )
        unique_vals = np.unique(result)
        assert set(unique_vals).issubset({0.0, 1.0}), (
            f"indicator=True produced non-binary values: {unique_vals}"
        )

    def test_indicator_detects_any_alt_allele(self):
        """indicator=True: a sample carrying >=1 alt allele gets indicator=1."""
        n = 5
        # All samples carry alt alleles in all categories
        geno = np.ones((n, 6), dtype=np.float64) * 2.0
        anno = np.array([1, 1, 2, 2, 3, 3])
        result = _aggregate_by_category(
            geno, anno, weights=[1.0, 1.0, 1.0], indicator=True, method="none"
        )
        assert np.all(result == 1.0), f"All samples should be carriers; got {result}"

    def test_indicator_zero_for_no_alt_allele(self):
        """indicator=True: sample with all ref genotypes gets indicator=0."""
        n = 5
        geno = np.zeros((n, 6), dtype=np.float64)
        anno = np.array([1, 1, 2, 2, 3, 3])
        result = _aggregate_by_category(
            geno, anno, weights=[1.0, 1.0, 1.0], indicator=True, method="none"
        )
        assert np.all(result == 0.0), f"All ref samples should have indicator=0; got {result}"

    def test_sum_method_returns_1_column(self):
        """method='sum': output shape is (n_samples, 1)."""
        geno, anno = self._make_geno_and_anno(n_samples=20)
        result = _aggregate_by_category(geno, anno, weights=[1.0, 2.0, 3.0], method="sum")
        assert result.shape == (20, 1), f"Expected (n, 1), got {result.shape}"

    def test_sum_method_returns_weighted_sum(self):
        """method='sum': result = 1*N_bmv + 2*N_dmv + 3*N_ptv per sample."""
        n = 20
        geno, anno = self._make_geno_and_anno(n_samples=n)
        weights = [1.0, 2.0, 3.0]
        result = _aggregate_by_category(geno, anno, weights=weights, method="sum")
        # Manual computation
        n_bmv = geno[:, 0] + geno[:, 1]
        n_dmv = geno[:, 2] + geno[:, 3]
        n_ptv = geno[:, 4] + geno[:, 5]
        expected = (1.0 * n_bmv + 2.0 * n_dmv + 3.0 * n_ptv)[:, np.newaxis]
        np.testing.assert_allclose(result, expected, err_msg="Weighted sum mismatch")

    def test_max_method_returns_1_column(self):
        """method='max': output shape is (n_samples, 1)."""
        geno, anno = self._make_geno_and_anno(n_samples=20)
        result = _aggregate_by_category(geno, anno, weights=[1.0, 2.0, 3.0], method="max")
        assert result.shape == (20, 1), f"Expected (n, 1), got {result.shape}"

    def test_max_method_returns_max_weighted_value(self):
        """method='max': result = max(1*N_bmv, 2*N_dmv, 3*N_ptv) per sample."""
        n = 20
        geno, anno = self._make_geno_and_anno(n_samples=n)
        weights = [1.0, 2.0, 3.0]
        result = _aggregate_by_category(geno, anno, weights=weights, method="max")
        n_bmv = geno[:, 0] + geno[:, 1]
        n_dmv = geno[:, 2] + geno[:, 3]
        n_ptv = geno[:, 4] + geno[:, 5]
        expected = np.maximum.reduce([1.0 * n_bmv, 2.0 * n_dmv, 3.0 * n_ptv])[:, np.newaxis]
        np.testing.assert_allclose(result, expected, err_msg="Weighted max mismatch")

    def test_empty_category_yields_zero_column(self):
        """All variants in one category: other category columns are zero."""
        n = 20
        rng = np.random.default_rng(99)
        geno = rng.choice([0, 1, 2], size=(n, 4), p=[0.6, 0.3, 0.1]).astype(np.float64)
        # All 4 variants are PTV (code 3); no BMV or DMV
        anno = np.array([3, 3, 3, 3])
        result = _aggregate_by_category(geno, anno, weights=[1.0, 1.0, 1.0], method="none")
        # BMV column (0) and DMV column (1) should be all zeros
        assert np.all(result[:, 0] == 0.0), "BMV column should be zero when no BMV variants"
        assert np.all(result[:, 1] == 0.0), "DMV column should be zero when no DMV variants"
        # PTV column (2) should match sum of all 4 variant columns
        np.testing.assert_allclose(result[:, 2], geno.sum(axis=1))

    def test_unknown_method_raises(self):
        """Unknown method raises ValueError."""
        geno, anno = self._make_geno_and_anno()
        with pytest.raises(ValueError, match="Unknown aggregation method"):
            _aggregate_by_category(geno, anno, weights=[1.0, 1.0, 1.0], method="invalid")


# ---------------------------------------------------------------------------
# TestCOASTBurdenTests
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestCOASTBurdenTests:
    """Test _run_burden_test() for quantitative and binary traits."""

    def _make_predictor(self, n: int, k: int, seed: int = 1) -> np.ndarray:
        """Generate random predictor matrix (n_samples, k)."""
        rng = np.random.default_rng(seed)
        return rng.standard_normal((n, k))

    def _make_binary_phenotype(self, n: int, seed: int = 1) -> np.ndarray:
        """Generate binary 0/1 phenotype."""
        rng = np.random.default_rng(seed)
        return rng.integers(0, 2, n).astype(np.float64)

    def _make_continuous_phenotype(self, n: int, seed: int = 1) -> np.ndarray:
        """Generate continuous phenotype."""
        rng = np.random.default_rng(seed)
        return rng.standard_normal(n)

    def test_quantitative_1df_returns_valid_p(self):
        """OLS 1-df test returns p in [0, 1)."""
        n = 100
        predictor = self._make_predictor(n, k=1, seed=10)
        phenotype = self._make_continuous_phenotype(n, seed=20)  # different seed
        p = _run_burden_test(predictor, phenotype, covariates=None, trait_type="quantitative")
        assert p is not None, "Expected a p-value for 1-df OLS"
        assert 0.0 <= p <= 1.0, f"p={p} out of [0, 1]"

    def test_quantitative_3df_returns_valid_p(self):
        """OLS 3-df F-test returns p in [0, 1)."""
        n = 100
        predictor = self._make_predictor(n, k=3, seed=11)
        phenotype = self._make_continuous_phenotype(n, seed=21)  # different seed
        p = _run_burden_test(predictor, phenotype, covariates=None, trait_type="quantitative")
        assert p is not None, "Expected a p-value for 3-df OLS F-test"
        assert 0.0 <= p <= 1.0, f"p={p} out of [0, 1]"

    def test_binary_1df_lrt_returns_valid_p(self):
        """Logit LRT 1-df returns p in [0, 1)."""
        n = 100
        predictor = self._make_predictor(n, k=1, seed=12)
        phenotype = self._make_binary_phenotype(n, seed=22)  # different seed
        p = _run_burden_test(predictor, phenotype, covariates=None, trait_type="binary")
        assert p is not None, "Expected a p-value for binary 1-df LRT"
        assert 0.0 <= p <= 1.0, f"p={p} out of [0, 1]"

    def test_binary_3df_lrt_returns_valid_p(self):
        """Logit LRT 3-df returns p in [0, 1)."""
        n = 100
        predictor = self._make_predictor(n, k=3, seed=13)
        phenotype = self._make_binary_phenotype(n, seed=23)  # different seed
        p = _run_burden_test(predictor, phenotype, covariates=None, trait_type="binary")
        assert p is not None, "Expected a p-value for binary 3-df LRT"
        assert 0.0 <= p <= 1.0, f"p={p} out of [0, 1]"

    def test_strong_signal_gives_small_p(self):
        """Predictor strongly correlated with phenotype: p < 0.05."""
        n = 200
        # Phenotype is binary; predictor is nearly perfectly correlated
        phenotype = np.array([1.0] * (n // 2) + [0.0] * (n // 2))
        # Cases have high predictor values; controls have low values
        rng = np.random.default_rng(77)
        predictor = np.zeros((n, 1))
        predictor[: n // 2, 0] = 1.0 + rng.standard_normal(n // 2) * 0.1
        predictor[n // 2 :, 0] = 0.0 + rng.standard_normal(n // 2) * 0.1
        p = _run_burden_test(predictor, phenotype, covariates=None, trait_type="binary")
        assert p is not None, "Expected a p-value for strong signal"
        assert p < 0.05, f"Expected p < 0.05 for strong signal, got {p:.6e}"

    def test_null_predictor_gives_large_mean_p(self):
        """
        Null predictor (uncorrelated with phenotype) gives moderate mean p over 20 trials.

        Mean p over independent trials should be in [0.2, 0.8] (roughly uniform).
        """
        rng = np.random.default_rng(42)
        n = 200
        p_values = []
        for _ in range(20):
            phenotype = rng.integers(0, 2, n).astype(np.float64)
            predictor = rng.standard_normal((n, 1))  # completely random, no correlation
            p = _run_burden_test(predictor, phenotype, covariates=None, trait_type="binary")
            if p is not None:
                p_values.append(p)

        assert len(p_values) >= 15, f"Too many None p-values: {20 - len(p_values)}"
        mean_p = float(np.mean(p_values))
        assert 0.15 < mean_p < 0.85, (
            f"Mean p={mean_p:.3f} outside [0.15, 0.85] — null distribution appears biased"
        )


# ---------------------------------------------------------------------------
# TestCOASTAllelicSKATWeights
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestCOASTAllelicSKATWeights:
    """Test _compute_allelic_skat_weights() and allelic SKAT component."""

    def test_weights_use_aaf_not_2aaf(self):
        """
        Verify w_j = sqrt(w_anno / (aaf*(1-aaf))), NOT sqrt(w_anno / (2*aaf*(1-aaf))).

        Decision IMPL-56: variance denominator uses aaf*(1-aaf) to match AllelicSeries R.
        """
        n = 100
        # Create genotype with known allele frequency
        geno = np.zeros((n, 3), dtype=np.float64)
        # BMV (code 1): 20% AAF => geno[:,0] has ~20 alt alleles
        geno[:20, 0] = 1.0
        # DMV (code 2): 30% AAF => geno[:,1] has ~30 alt alleles
        geno[:30, 1] = 1.0
        # PTV (code 3): 10% AAF => geno[:,2] has ~10 alt alleles
        geno[:10, 2] = 1.0

        anno_codes = np.array([1, 2, 3])
        coast_weights = [1.0, 2.0, 3.0]

        result = _compute_allelic_skat_weights(geno, anno_codes, coast_weights)

        # Manual computation: aaf = mean(geno[:,j]) / 2
        aaf_bmv = geno[:, 0].mean() / 2.0  # 20/100/2 = 0.1
        aaf_dmv = geno[:, 1].mean() / 2.0  # 30/100/2 = 0.15
        aaf_ptv = geno[:, 2].mean() / 2.0  # 10/100/2 = 0.05

        # w_j = sqrt(coast_weight / (aaf*(1-aaf)))
        expected_bmv = np.sqrt(1.0 / (aaf_bmv * (1.0 - aaf_bmv)))
        expected_dmv = np.sqrt(2.0 / (aaf_dmv * (1.0 - aaf_dmv)))
        expected_ptv = np.sqrt(3.0 / (aaf_ptv * (1.0 - aaf_ptv)))

        np.testing.assert_allclose(result[0], expected_bmv, rtol=1e-6, err_msg="BMV weight")
        np.testing.assert_allclose(result[1], expected_dmv, rtol=1e-6, err_msg="DMV weight")
        np.testing.assert_allclose(result[2], expected_ptv, rtol=1e-6, err_msg="PTV weight")

        # Verify NOT using 2*aaf*(1-aaf)
        wrong_bmv = np.sqrt(1.0 / (2.0 * aaf_bmv * (1.0 - aaf_bmv)))
        assert not np.isclose(result[0], wrong_bmv, rtol=1e-6), (
            "Weight matches 2*aaf*(1-aaf) denominator — expected aaf*(1-aaf)"
        )

    def test_monomorphic_variant_clamped(self):
        """Variant with aaf=0 is clamped to 1e-8; no division by zero."""
        n = 100
        geno = np.zeros((n, 3), dtype=np.float64)
        # All variants are monomorphic (no alt alleles)
        anno_codes = np.array([1, 2, 3])
        coast_weights = [1.0, 2.0, 3.0]

        # Should not raise
        result = _compute_allelic_skat_weights(geno, anno_codes, coast_weights)

        assert result.shape == (3,), f"Expected 3 weights, got {result.shape}"
        assert np.all(np.isfinite(result)), f"Expected finite weights, got {result}"
        assert np.all(result > 0), f"Expected positive weights, got {result}"

    def test_weights_shape_matches_n_variants(self):
        """Output weight array has shape (n_variants,)."""
        rng = np.random.default_rng(42)
        n, k = 50, 9
        geno = rng.choice([0, 1, 2], size=(n, k), p=[0.6, 0.3, 0.1]).astype(np.float64)
        anno_codes = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3])
        result = _compute_allelic_skat_weights(geno, anno_codes, [1.0, 2.0, 3.0])
        assert result.shape == (k,), f"Expected ({k},), got {result.shape}"

    def test_higher_weight_gives_larger_weight(self):
        """PTV weight (coast_weight=3) > BMV weight (coast_weight=1) for same AAF."""
        n = 100
        # All three variants have same AAF (10% alt)
        geno = np.zeros((n, 3), dtype=np.float64)
        geno[:10, :] = 1.0  # 10% AAF for all three
        anno_codes = np.array([1, 2, 3])  # BMV, DMV, PTV
        coast_weights = [1.0, 2.0, 3.0]

        result = _compute_allelic_skat_weights(geno, anno_codes, coast_weights)

        # Same AAF => ratio determined by coast_weight: w3 > w2 > w1
        assert result[2] > result[1] > result[0], (
            f"Expected w3 > w2 > w1 for equal AAF, got {result}"
        )

    def test_allelic_skat_via_backend_returns_valid_p(self):
        """Full allelic SKAT via PythonCOASTBackend._run_allelic_skat() returns p in (0, 1]."""
        from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

        rng = np.random.default_rng(55)
        n, k = 100, 9
        geno = rng.choice([0, 1, 2], size=(n, k), p=[0.6, 0.3, 0.1]).astype(np.float64)
        phenotype = rng.integers(0, 2, n).astype(np.float64)
        anno_codes = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3])
        coast_weights = [1.0, 2.0, 3.0]

        skat_backend = PythonSKATBackend()
        skat_backend.detect_environment()
        null_model = skat_backend.fit_null_model(phenotype, None, "binary")

        coast_backend = PythonCOASTBackend()
        skat_weights = _compute_allelic_skat_weights(geno, anno_codes, coast_weights)
        p = coast_backend._run_allelic_skat(geno, skat_weights, null_model)

        assert 0.0 < p <= 1.0, f"Allelic SKAT p-value {p} out of (0, 1]"


# ---------------------------------------------------------------------------
# TestCOASTOmnibus
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestCOASTOmnibus:
    """Test full test_gene() pipeline: 6 burden + 1 SKAT + Cauchy combination."""

    @pytest.fixture(scope="class")
    def backend_and_null(self):
        """Class-scoped backend + null model to avoid refitting for each test."""
        from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

        rng = np.random.default_rng(100)
        n = 100
        phenotype = rng.integers(0, 2, n).astype(np.float64)
        skat_backend = PythonSKATBackend()
        skat_backend.detect_environment()
        null_model = skat_backend.fit_null_model(phenotype, None, "binary")
        coast_backend = PythonCOASTBackend()
        return coast_backend, null_model, phenotype

    def _make_geno_and_anno(self, n: int = 100, seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
        """Build (geno, anno_codes) with 3 BMV + 3 DMV + 3 PTV variants."""
        rng = np.random.default_rng(seed)
        geno = rng.choice([0, 1, 2], size=(n, 9), p=[0.6, 0.3, 0.1]).astype(np.float64)
        anno_codes = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3])
        return geno, anno_codes

    def test_omnibus_returns_valid_p(self, backend_and_null):
        """test_gene() returns p_value in (0, 1]."""
        backend, null_model, phenotype = backend_and_null
        geno, anno = self._make_geno_and_anno()
        result = backend.test_gene(
            gene="TEST_GENE",
            geno_filtered=geno,
            anno_codes_filtered=anno,
            phenotype=phenotype,
            covariates=None,
            coast_weights=[1.0, 2.0, 3.0],
            trait_type="binary",
            null_model=null_model,
        )
        p = result["p_value"]
        assert p is not None, "Expected a p_value from test_gene()"
        assert 0.0 < p <= 1.0, f"p_value={p} out of (0, 1]"

    def test_result_has_6_burden_pvalues(self, backend_and_null):
        """test_gene() result contains exactly 6 burden p-values."""
        backend, null_model, phenotype = backend_and_null
        geno, anno = self._make_geno_and_anno(seed=43)
        result = backend.test_gene(
            gene="TEST_GENE",
            geno_filtered=geno,
            anno_codes_filtered=anno,
            phenotype=phenotype,
            covariates=None,
            coast_weights=[1.0, 2.0, 3.0],
            trait_type="binary",
            null_model=null_model,
        )
        burden_pvals = result["burden_p_values"]
        assert len(burden_pvals) == 6, (
            f"Expected 6 burden p-values, got {len(burden_pvals)}: {burden_pvals}"
        )

    def test_burden_labels_match_expected(self, backend_and_null):
        """Burden labels are: base_count, base_ind, sum_count, sum_ind, max_count, max_ind."""
        backend, null_model, phenotype = backend_and_null
        geno, anno = self._make_geno_and_anno(seed=44)
        result = backend.test_gene(
            gene="TEST_GENE",
            geno_filtered=geno,
            anno_codes_filtered=anno,
            phenotype=phenotype,
            covariates=None,
            coast_weights=[1.0, 2.0, 3.0],
            trait_type="binary",
            null_model=null_model,
        )
        expected_labels = ["base_count", "base_ind", "sum_count", "sum_ind", "max_count", "max_ind"]
        assert result["burden_labels"] == expected_labels, (
            f"Expected labels {expected_labels}, got {result['burden_labels']}"
        )

    def test_result_has_skat_pvalue(self, backend_and_null):
        """test_gene() result contains skat_p_value in (0, 1]."""
        backend, null_model, phenotype = backend_and_null
        geno, anno = self._make_geno_and_anno(seed=45)
        result = backend.test_gene(
            gene="TEST_GENE",
            geno_filtered=geno,
            anno_codes_filtered=anno,
            phenotype=phenotype,
            covariates=None,
            coast_weights=[1.0, 2.0, 3.0],
            trait_type="binary",
            null_model=null_model,
        )
        p_skat = result["skat_p_value"]
        if p_skat is not None:
            assert 0.0 < p_skat <= 1.0, f"skat_p_value={p_skat} out of (0, 1]"

    def test_deterministic_across_calls(self, backend_and_null):
        """Same input produces identical output across two calls (no random state)."""
        backend, null_model, phenotype = backend_and_null
        geno, anno = self._make_geno_and_anno(seed=46)
        kwargs = {
            "gene": "TEST_GENE",
            "geno_filtered": geno,
            "anno_codes_filtered": anno,
            "phenotype": phenotype,
            "covariates": None,
            "coast_weights": [1.0, 2.0, 3.0],
            "trait_type": "binary",
            "null_model": null_model,
        }
        r1 = backend.test_gene(**kwargs)
        r2 = backend.test_gene(**kwargs)

        assert r1["p_value"] == r2["p_value"], (
            f"Non-deterministic omnibus p: {r1['p_value']} != {r2['p_value']}"
        )
        assert r1["burden_p_values"] == r2["burden_p_values"], "Non-deterministic burden p-values"

    def test_null_phenotype_moderate_mean_p(self):
        """
        Under null (no association), mean omnibus p over 20 genes should be in [0.15, 0.85].

        This is a distribution sanity check: COAST should not be systematically biased.
        """
        from variantcentrifuge.association.backends.python_backend import PythonSKATBackend

        rng = np.random.default_rng(7777)
        n = 100
        phenotype = rng.integers(0, 2, n).astype(np.float64)

        skat_backend = PythonSKATBackend()
        skat_backend.detect_environment()
        null_model = skat_backend.fit_null_model(phenotype, None, "binary")

        coast = PythonCOASTBackend()
        anno_codes = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3])

        p_values = []
        for i in range(20):
            geno = rng.choice([0, 1, 2], size=(n, 9), p=[0.6, 0.3, 0.1]).astype(np.float64)
            result = coast.test_gene(
                gene=f"GENE_{i}",
                geno_filtered=geno,
                anno_codes_filtered=anno_codes,
                phenotype=phenotype,
                covariates=None,
                coast_weights=[1.0, 2.0, 3.0],
                trait_type="binary",
                null_model=null_model,
            )
            if result["p_value"] is not None:
                p_values.append(result["p_value"])

        assert len(p_values) >= 15, f"Too many None p-values: {20 - len(p_values)}"
        mean_p = float(np.mean(p_values))
        assert 0.15 < mean_p < 0.85, (
            f"Mean omnibus p={mean_p:.3f} outside [0.15, 0.85] — distribution appears biased"
        )


# ---------------------------------------------------------------------------
# TestPurePythonCOASTTestLifecycle
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestPurePythonCOASTTestLifecycle:
    """Integration tests for PurePythonCOASTTest engine-pattern lifecycle."""

    def test_parallel_safe_is_true(self):
        """PurePythonCOASTTest.parallel_safe is True (no rpy2 restriction)."""
        assert PurePythonCOASTTest.parallel_safe is True

    def test_name_is_coast(self):
        """PurePythonCOASTTest().name == 'coast' (same as COASTTest for registry swap)."""
        assert PurePythonCOASTTest().name == "coast"

    def test_effect_column_names_all_none(self):
        """effect_column_names() returns all None values (no effect size for COAST)."""
        test = PurePythonCOASTTest()
        col_names = test.effect_column_names()
        assert col_names["effect"] is None
        assert col_names["se"] is None
        assert col_names["ci_lower"] is None
        assert col_names["ci_upper"] is None

    def test_effect_column_names_has_all_four_keys(self):
        """effect_column_names() returns dict with exactly 4 keys."""
        test = PurePythonCOASTTest()
        col_names = test.effect_column_names()
        assert set(col_names.keys()) == {"effect", "se", "ci_lower", "ci_upper"}

    def test_full_lifecycle(self):
        """Full lifecycle: check_dependencies -> prepare -> run -> finalize."""
        test = PurePythonCOASTTest()
        test.check_dependencies()
        test.prepare(1)

        contingency_data, config = _make_coast_contingency_data()
        result = test.run("TEST_GENE", contingency_data, config)
        test.finalize()

        # Verify TestResult fields
        assert result.gene == "TEST_GENE"
        assert result.test_name == "coast"
        assert result.effect_size is None
        assert result.se is None
        assert result.ci_lower is None
        assert result.ci_upper is None
        # p_value can be None or float
        if result.p_value is not None:
            assert 0.0 < result.p_value <= 1.0, f"p_value={result.p_value} out of (0, 1]"

    def test_extra_keys_match_coast_test(self):
        """TestResult.extra contains expected COAST sub-result keys."""
        test = PurePythonCOASTTest()
        test.check_dependencies()
        test.prepare(1)

        contingency_data, config = _make_coast_contingency_data()
        result = test.run("TEST_GENE", contingency_data, config)
        test.finalize()

        if result.p_value is not None:
            # Full run: check all expected keys
            expected_keys = {
                "coast_burden_p_value",
                "coast_skat_p_value",
                "coast_n_bmv",
                "coast_n_dmv",
                "coast_n_ptv",
            }
            missing = expected_keys - set(result.extra.keys())
            assert not missing, (
                f"Missing expected extra keys: {missing}. Got: {set(result.extra.keys())}"
            )
            # Verify counts
            assert result.extra["coast_n_bmv"] == 3
            assert result.extra["coast_n_dmv"] == 3
            assert result.extra["coast_n_ptv"] == 3

    def test_no_genotype_matrix_returns_none(self):
        """Missing genotype_matrix: p_value=None with coast_skip_reason."""
        from variantcentrifuge.association.base import AssociationConfig

        test = PurePythonCOASTTest()
        test.check_dependencies()
        test.prepare(1)

        # No genotype_matrix key
        contingency_data = {
            "proband_count": 50,
            "control_count": 50,
            "n_qualifying_variants": 9,
        }
        config = AssociationConfig(trait_type="binary")
        result = test.run("NO_GENO", contingency_data, config)
        test.finalize()

        assert result.p_value is None
        assert "coast_skip_reason" in result.extra
        assert "NO_GENOTYPE_MATRIX" in result.extra["coast_skip_reason"]

    def test_partial_categories_proceeds(self):
        """
        Only PTV variants (no BMV/DMV): COAST-03 partial-category fallback.

        With the COAST-03 fix, COAST proceeds when at least 1 category is present.
        p_value is numeric (not None), coast_status == 'partial', and
        coast_missing_categories reports BMV and DMV.
        """
        from variantcentrifuge.association.base import AssociationConfig

        test = PurePythonCOASTTest()
        test.check_dependencies()
        test.prepare(1)

        n = 60
        # Build gene_df with ONLY PTV variants
        gene_df = _make_gene_df(
            effects=["frameshift_variant", "stop_gained", "splice_donor_variant"],
            impacts=["HIGH", "HIGH", "HIGH"],
            sift=[".", ".", "."],
            polyphen=[".", ".", "."],
        )
        rng = np.random.default_rng(42)
        geno = rng.choice([0, 1, 2], size=(n, 3), p=[0.6, 0.3, 0.1]).astype(np.float64)
        phenotype = np.array([1.0] * 30 + [0.0] * 30)

        contingency_data = {
            "genotype_matrix": geno,
            "phenotype_vector": phenotype,
            "gene_df": gene_df,
            "proband_count": 30,
            "control_count": 30,
            "n_qualifying_variants": 3,
            "covariate_matrix": None,
        }
        config = AssociationConfig(trait_type="binary")
        result = test.run("PTV_ONLY_GENE", contingency_data, config)
        test.finalize()

        # COAST-03: partial-category fallback -- test should proceed, not skip
        assert result.p_value is not None, (
            "COAST-03: expected numeric p_value for partial COAST (PTV-only gene)"
        )
        assert result.extra.get("coast_status") == "partial"
        missing_cats = result.extra.get("coast_missing_categories", "")
        assert "BMV" in missing_cats
        assert "DMV" in missing_cats
        assert result.extra.get("coast_n_bmv") == 0
        assert result.extra.get("coast_n_dmv") == 0
        assert result.extra.get("coast_n_ptv") == 3

    def test_engine_registration_swap_with_python_backend(self):
        """coast_backend='python': from_names(['coast'], config) returns PurePythonCOASTTest."""
        from variantcentrifuge.association.base import AssociationConfig
        from variantcentrifuge.association.engine import AssociationEngine

        config = AssociationConfig(trait_type="binary", coast_backend="python")
        engine = AssociationEngine.from_names(["coast"], config)

        # _tests is a dict {test_name: test_instance}
        assert len(engine._tests) == 1
        assert "coast" in engine._tests, f"Expected 'coast' key in engine._tests: {engine._tests}"
        assert isinstance(engine._tests["coast"], PurePythonCOASTTest), (
            f"Expected PurePythonCOASTTest with coast_backend='python', "
            f"got {type(engine._tests['coast'])}"
        )

    def test_multiple_genes_sequential(self):
        """Process 3 genes sequentially; null model is reused across genes."""
        test = PurePythonCOASTTest()
        test.check_dependencies()
        test.prepare(3)

        results = []
        for i in range(3):
            contingency_data, config = _make_coast_contingency_data(seed=i + 10)
            result = test.run(f"GENE_{i}", contingency_data, config)
            results.append(result)

        test.finalize()

        assert len(results) == 3
        for i, result in enumerate(results):
            assert result.gene == f"GENE_{i}"
            assert result.test_name == "coast"

    def test_quantitative_trait_lifecycle(self):
        """Quantitative phenotype runs without error and returns valid p."""
        test = PurePythonCOASTTest()
        test.check_dependencies()
        test.prepare(1)

        contingency_data, config = _make_coast_contingency_data(trait_type="quantitative", seed=200)
        result = test.run("QUANT_GENE", contingency_data, config)
        test.finalize()

        assert result.test_name == "coast"
        if result.p_value is not None:
            assert 0.0 < result.p_value <= 1.0

    def test_genes_processed_counter_increments(self):
        """_genes_processed increments for each successful run."""
        test = PurePythonCOASTTest()
        test.check_dependencies()
        test.prepare(5)

        initial_count = test._genes_processed
        contingency_data, config = _make_coast_contingency_data()
        result = test.run("GENE_1", contingency_data, config)
        test.finalize()

        if result.p_value is not None:
            # Counter incremented for successful (non-skipped) runs
            assert test._genes_processed == initial_count + 1

    def test_no_gene_df_returns_none(self):
        """Missing gene_df: p_value=None with NO_GENE_DF skip reason."""
        from variantcentrifuge.association.base import AssociationConfig

        test = PurePythonCOASTTest()
        test.check_dependencies()
        test.prepare(1)

        rng = np.random.default_rng(42)
        n = 60
        geno = rng.choice([0, 1, 2], size=(n, 9), p=[0.6, 0.3, 0.1]).astype(np.float64)
        phenotype = np.array([1.0] * 30 + [0.0] * 30)

        contingency_data = {
            "genotype_matrix": geno,
            "phenotype_vector": phenotype,
            # No gene_df
            "proband_count": 30,
            "control_count": 30,
            "n_qualifying_variants": 9,
        }
        config = AssociationConfig(trait_type="binary")
        result = test.run("NO_GENE_DF", contingency_data, config)
        test.finalize()

        assert result.p_value is None
        assert "coast_skip_reason" in result.extra
        assert "NO_GENE_DF" in result.extra["coast_skip_reason"]
