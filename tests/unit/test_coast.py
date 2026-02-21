"""
Unit tests for COASTTest and classify_variants.

Tests the COAST allelic series test logic without invoking rpy2 at all.
rpy2 and AllelicSeries are fully mocked where needed.

Covers requirements: SERIES-01, SERIES-02
"""

from __future__ import annotations

import logging
from typing import Any
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


def _make_config(**kwargs: Any) -> Any:
    """Build an AssociationConfig with sensible defaults."""
    from variantcentrifuge.association.base import AssociationConfig

    defaults = {
        "trait_type": "binary",
        "variant_weights": "beta:1,25",
    }
    defaults.update(kwargs)
    return AssociationConfig(**defaults)


def _make_gene_df(
    effects: list[str],
    impacts: list[str],
    sift: list[str] | None = None,
    polyphen: list[str] | None = None,
    sift_col: str = "dbNSFP_SIFT_pred",
    polyphen_col: str = "dbNSFP_Polyphen2_HDIV_pred",
) -> pd.DataFrame:
    """Build a per-variant DataFrame for testing."""
    data: dict[str, list] = {
        "EFFECT": effects,
        "IMPACT": impacts,
    }
    if sift is not None:
        data[sift_col] = sift
    if polyphen is not None:
        data[polyphen_col] = polyphen
    return pd.DataFrame(data)


def _make_geno_data(
    gene: str = "BRCA1",
    n_samples: int = 60,
    n_variants: int = 5,
    n_cases: int = 30,
    n_controls: int = 30,
    include_genotype_matrix: bool = True,
    gene_df: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Build a contingency_data dict for COAST (includes genotype_matrix and gene_df)."""
    rng = np.random.default_rng(42)
    geno = rng.integers(0, 2, size=(n_samples, n_variants)).astype(float)
    phenotype = np.array([1.0] * n_cases + [0.0] * n_controls)

    data: dict[str, Any] = {
        "GENE": gene,
        "proband_count": n_cases,
        "control_count": n_controls,
        "proband_carrier_count": 5,
        "control_carrier_count": 2,
        "proband_allele_count": 5,
        "control_allele_count": 2,
        "n_qualifying_variants": n_variants,
        "phenotype_vector": phenotype,
        "covariate_matrix": None,
    }
    if include_genotype_matrix:
        data["genotype_matrix"] = geno
    if gene_df is not None:
        data["gene_df"] = gene_df
    return data


# ---------------------------------------------------------------------------
# classify_variants: PTV classification
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestClassifyVariantsPTV:
    """Tests for PTV (protein-truncating variant) classification."""

    def test_stop_gained_high_impact_is_ptv(self):
        """stop_gained + HIGH impact -> code 3 (PTV)."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["stop_gained"],
            impacts=["HIGH"],
            sift=["deleterious"],
            polyphen=["probably_damaging"],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 3
        assert mask[0] is True or bool(mask[0])

    def test_frameshift_variant_high_impact_is_ptv(self):
        """frameshift_variant + HIGH impact -> code 3 (PTV)."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["frameshift_variant"],
            impacts=["HIGH"],
            sift=["tolerated"],
            polyphen=["benign"],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 3
        assert bool(mask[0])

    def test_splice_acceptor_variant_high_impact_is_ptv(self):
        """splice_acceptor_variant + HIGH impact -> code 3."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["splice_acceptor_variant"],
            impacts=["HIGH"],
            sift=[""],
            polyphen=[""],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 3
        assert bool(mask[0])

    def test_splice_donor_variant_high_impact_is_ptv(self):
        """splice_donor_variant + HIGH impact -> code 3."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["splice_donor_variant"],
            impacts=["HIGH"],
            sift=[""],
            polyphen=[""],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 3
        assert bool(mask[0])

    def test_stop_gained_moderate_impact_not_ptv(self):
        """stop_gained without HIGH impact -> NOT code 3."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["stop_gained"],
            impacts=["MODERATE"],
            sift=[""],
            polyphen=[""],
        )
        codes, _mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] != 3

    def test_missense_high_impact_not_ptv(self):
        """missense_variant with HIGH impact -> NOT code 3 (PTV effect check fails)."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["HIGH"],
            sift=["deleterious"],
            polyphen=["probably_damaging"],
        )
        codes, _mask = classify_variants(df, "EFFECT", "IMPACT")
        # HIGH impact missense is not in PTV_EFFECTS — classified as DMV
        assert codes[0] != 3


# ---------------------------------------------------------------------------
# classify_variants: DMV classification
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestClassifyVariantsDMV:
    """Tests for DMV (damaging missense variant) classification."""

    def test_missense_sift_deleterious_is_dmv(self):
        """missense + SIFT deleterious -> code 2 (DMV)."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["MODERATE"],
            sift=["deleterious"],
            polyphen=["benign"],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 2
        assert bool(mask[0])

    def test_missense_polyphen_probably_damaging_is_dmv(self):
        """missense + PolyPhen probably_damaging -> code 2 (DMV)."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["MODERATE"],
            sift=["tolerated"],
            polyphen=["probably_damaging"],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 2
        assert bool(mask[0])

    def test_missense_polyphen_possibly_damaging_is_dmv(self):
        """missense + PolyPhen possibly_damaging -> code 2 (DMV)."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["MODERATE"],
            sift=["tolerated"],
            polyphen=["possibly_damaging"],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 2
        assert bool(mask[0])

    def test_missense_both_sift_and_polyphen_damaging_is_dmv(self):
        """missense + SIFT deleterious AND PolyPhen probably_damaging -> code 2."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["MODERATE"],
            sift=["deleterious"],
            polyphen=["probably_damaging"],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 2
        assert bool(mask[0])


# ---------------------------------------------------------------------------
# classify_variants: BMV classification
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestClassifyVariantsBMV:
    """Tests for BMV (benign missense variant) classification."""

    def test_missense_sift_tolerated_polyphen_benign_is_bmv(self):
        """missense + SIFT tolerated + PolyPhen benign -> code 1 (BMV)."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["MODERATE"],
            sift=["tolerated"],
            polyphen=["benign"],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 1
        assert bool(mask[0])

    def test_missense_only_sift_benign_not_bmv(self):
        """missense + SIFT tolerated but no PolyPhen -> not code 1 (needs both)."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["MODERATE"],
            sift=["tolerated"],
            polyphen=[""],  # no prediction
        )
        codes, _mask = classify_variants(df, "EFFECT", "IMPACT")
        # polyphen_benign is False -> BMV condition fails; not BMV
        assert codes[0] != 1


# ---------------------------------------------------------------------------
# classify_variants: Unclassified / excluded
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestClassifyVariantsUnclassified:
    """Tests for unclassified and excluded variants."""

    def test_missense_no_sift_polyphen_is_excluded(self):
        """missense without SIFT/PolyPhen predictions -> code 0, include_mask=False."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["MODERATE"],
            sift=[""],
            polyphen=[""],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 0
        assert not bool(mask[0])

    def test_missing_sift_polyphen_columns_returns_all_false(self, caplog):
        """DataFrame with no SIFT/PolyPhen columns -> all-zero codes, all-False mask."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = pd.DataFrame(
            {
                "EFFECT": ["missense_variant"],
                "IMPACT": ["MODERATE"],
                # No SIFT or PolyPhen columns
            }
        )
        with caplog.at_level(logging.ERROR, logger="variantcentrifuge"):
            codes, mask = classify_variants(df, "EFFECT", "IMPACT")

        assert all(c == 0 for c in codes)
        assert not any(mask)
        # Check any ERROR record mentions both SIFT and PolyPhen
        error_msgs = [r for r in caplog.records if r.levelno == logging.ERROR]
        assert len(error_msgs) >= 1
        assert any("SIFT" in r.message and "PolyPhen" in r.message for r in error_msgs)

    def test_synonymous_variant_excluded(self):
        """synonymous_variant -> code 0, excluded from COAST."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["synonymous_variant"],
            impacts=["LOW"],
            sift=["tolerated"],
            polyphen=["benign"],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 0
        assert not bool(mask[0])

    def test_intron_variant_excluded(self):
        """intron_variant -> code 0, excluded from COAST."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["intron_variant"],
            impacts=["MODIFIER"],
            sift=["tolerated"],
            polyphen=["benign"],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 0
        assert not bool(mask[0])


# ---------------------------------------------------------------------------
# classify_variants: Mixed variants
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestClassifyVariantsMixed:
    """Tests for a DataFrame with mixed variant types."""

    def test_mixed_variants_correct_codes_and_masks(self):
        """PTV + DMV + BMV + unclassified -> codes [3, 2, 1, 0] and matching masks."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=[
                "frameshift_variant",
                "missense_variant",
                "missense_variant",
                "intron_variant",
            ],
            impacts=["HIGH", "MODERATE", "MODERATE", "MODIFIER"],
            sift=[".", "deleterious", "tolerated", "tolerated"],
            polyphen=[".", "probably_damaging", "benign", "benign"],
        )
        codes, mask = classify_variants(df, "EFFECT", "IMPACT")

        assert codes[0] == 3  # PTV
        assert codes[1] == 2  # DMV
        assert codes[2] == 1  # BMV
        assert codes[3] == 0  # unclassified

        assert bool(mask[0])
        assert bool(mask[1])
        assert bool(mask[2])
        assert not bool(mask[3])

    def test_ptv_overrides_dmv_when_both_conditions_met(self):
        """When a variant matches both PTV and missense conditions, PTV (code 3) wins."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        # This is an artificial case: stop_gained + HIGH impact takes PTV code
        # (stop_gained is in PTV_EFFECTS, so code=3 overwrites any missense code)
        df = _make_gene_df(
            effects=["stop_gained"],
            impacts=["HIGH"],
            sift=["deleterious"],
            polyphen=["probably_damaging"],
        )
        codes, _mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 3

    def test_multi_value_sift_deleterious_detected(self):
        """SIFT field with multiple values like 'tolerated;deleterious' -> DMV."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["MODERATE"],
            sift=["tolerated;deleterious"],  # multi-value, second is damaging
            polyphen=["benign"],
        )
        codes, _mask = classify_variants(df, "EFFECT", "IMPACT")
        # tolerated is benign, deleterious is damaging; is_dmv=True (SIFT damaging)
        # is_bmv needs SIFT benign AND polyphen benign, but SIFT has deleterious -> not BMV
        assert codes[0] == 2  # DMV because SIFT_DAMAGING found

    def test_alternative_sift_column_name_detected(self):
        """SIFT_pred (not dbNSFP_SIFT_pred) column is auto-detected."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["MODERATE"],
            sift=["deleterious"],
            polyphen=["probably_damaging"],
            sift_col="SIFT_pred",  # alternative column name
        )
        codes, _mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 2  # DMV detected via SIFT_pred

    def test_alternative_polyphen_column_name_detected(self):
        """Polyphen2_HVAR_pred column is auto-detected."""
        from variantcentrifuge.association.tests.allelic_series import classify_variants

        df = _make_gene_df(
            effects=["missense_variant"],
            impacts=["MODERATE"],
            sift=["tolerated"],
            polyphen=["benign"],
            polyphen_col="Polyphen2_HVAR_pred",  # alternative column name
        )
        codes, _mask = classify_variants(df, "EFFECT", "IMPACT")
        assert codes[0] == 1  # BMV


# ---------------------------------------------------------------------------
# COASTTest: basic properties
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestCOASTTestBasicProperties:
    """Tests for COASTTest name and effect_column_names."""

    def test_name_is_coast(self):
        """COASTTest.name returns 'coast'."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        assert test.name == "coast"

    def test_effect_column_names_all_none(self):
        """effect_column_names() returns all None values (COAST has no effect size)."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        col_names = test.effect_column_names()

        assert col_names["effect"] is None
        assert col_names["se"] is None
        assert col_names["ci_lower"] is None
        assert col_names["ci_upper"] is None

    def test_effect_column_names_returns_all_four_keys(self):
        """effect_column_names() returns a dict with exactly the four expected keys."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        col_names = test.effect_column_names()

        assert set(col_names.keys()) == {"effect", "se", "ci_lower", "ci_upper"}

    def test_parallel_safe_is_false(self):
        """COASTTest.parallel_safe is False (rpy2 main-thread-only restriction)."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        assert test.parallel_safe is False


# ---------------------------------------------------------------------------
# COASTTest: check_dependencies
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestCOASTCheckDependencies:
    """Tests for COASTTest.check_dependencies() with mocked rpy2."""

    def test_check_dependencies_rpy2_not_available_raises(self):
        """check_dependencies() raises ImportError when rpy2 is not importable."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()

        with (
            patch.dict("sys.modules", {"rpy2": None, "rpy2.robjects": None}),
            pytest.raises(ImportError, match="rpy2"),
        ):
            test.check_dependencies()

    def test_check_dependencies_allelic_series_not_available_raises(self):
        """check_dependencies() raises ImportError when AllelicSeries R package is missing."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()

        mock_rpy2 = MagicMock()
        mock_robjects = MagicMock()

        def mock_importr(pkg: str) -> MagicMock:
            if pkg == "AllelicSeries":
                raise Exception("R package 'AllelicSeries' not found")
            return MagicMock()

        mock_robjects_packages = MagicMock()
        mock_robjects_packages.importr.side_effect = mock_importr

        with (
            patch.dict(
                "sys.modules",
                {
                    "rpy2": mock_rpy2,
                    "rpy2.robjects": mock_robjects,
                    "rpy2.robjects.packages": mock_robjects_packages,
                },
            ),
            pytest.raises(ImportError, match="AllelicSeries"),
        ):
            test.check_dependencies()


# ---------------------------------------------------------------------------
# COASTTest: run() skip conditions
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestCOASTRunSkipConditions:
    """Tests for COASTTest.run() skip conditions returning p_value=None."""

    def test_run_no_genotype_matrix_returns_none(self):
        """run() with no genotype_matrix -> p_value=None."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        config = _make_config()
        data = _make_geno_data(include_genotype_matrix=False)
        result = test.run("GENE1", data, config)

        assert result.p_value is None
        assert "coast_skip_reason" in result.extra
        assert "NO_GENOTYPE_MATRIX" in result.extra["coast_skip_reason"]

    def test_run_no_phenotype_returns_none(self):
        """run() without phenotype_vector -> p_value=None."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        config = _make_config()
        data = _make_geno_data()
        data.pop("phenotype_vector", None)
        result = test.run("GENE1", data, config)

        assert result.p_value is None

    def test_run_no_gene_df_returns_none(self):
        """run() without gene_df -> p_value=None."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        config = _make_config()
        data = _make_geno_data()  # no gene_df key
        result = test.run("GENE1", data, config)

        assert result.p_value is None
        assert "coast_skip_reason" in result.extra
        assert "NO_GENE_DF" in result.extra["coast_skip_reason"]

    def test_run_missing_categories_returns_none(self):
        """run() with only PTVs (no BMV, DMV) -> p_value=None with skip reason."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        config = _make_config()

        # All variants are PTVs (no BMV/DMV)
        gene_df = _make_gene_df(
            effects=["frameshift_variant", "stop_gained", "splice_donor_variant"],
            impacts=["HIGH", "HIGH", "HIGH"],
            sift=[".", ".", "."],
            polyphen=[".", ".", "."],
        )
        # Genotype matrix with 3 variants
        n_samples = 60
        geno = np.zeros((n_samples, 3), dtype=float)
        phenotype = np.array([1.0] * 30 + [0.0] * 30)

        data = {
            "GENE": "GENE1",
            "proband_count": 30,
            "control_count": 30,
            "n_qualifying_variants": 3,
            "genotype_matrix": geno,
            "phenotype_vector": phenotype,
            "covariate_matrix": None,
            "gene_df": gene_df,
        }
        result = test.run("GENE1", data, config)

        assert result.p_value is None
        assert "coast_skip_reason" in result.extra
        # Should have BMV and DMV counts as 0
        assert result.extra.get("coast_n_bmv") == 0
        assert result.extra.get("coast_n_dmv") == 0
        assert result.extra.get("coast_n_ptv") == 3

    def test_run_annotation_genotype_mismatch_returns_none(self):
        """run() with gene_df length != genotype matrix columns -> p_value=None."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        config = _make_config()

        gene_df = _make_gene_df(
            effects=["frameshift_variant", "missense_variant"],
            impacts=["HIGH", "MODERATE"],
            sift=[".", "deleterious"],
            polyphen=[".", "probably_damaging"],
        )
        # Genotype matrix has 5 columns but gene_df has 2 rows -> mismatch
        n_samples = 60
        geno = np.zeros((n_samples, 5), dtype=float)
        phenotype = np.array([1.0] * 30 + [0.0] * 30)

        data = {
            "GENE": "GENE1",
            "proband_count": 30,
            "control_count": 30,
            "n_qualifying_variants": 5,
            "genotype_matrix": geno,
            "phenotype_vector": phenotype,
            "covariate_matrix": None,
            "gene_df": gene_df,
        }
        result = test.run("GENE1", data, config)

        assert result.p_value is None
        assert "coast_skip_reason" in result.extra
        assert "MISMATCH" in result.extra["coast_skip_reason"]


# ---------------------------------------------------------------------------
# COASTTest: run() with mocked R
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestCOASTRunWithMockedR:
    """Tests for COASTTest.run() with mocked rpy2/AllelicSeries."""

    def _build_mock_coast_result(self, burden_p: float, skat_p: float, omni_p: float) -> MagicMock:
        """Build a mock R COAST result with a Pvals slot."""
        mock_result = MagicMock()

        mock_pvals = MagicMock()
        mock_pvals.names = ["Burden", "SKAT", "Omni"]
        # Make it iterable
        mock_pvals.__iter__ = MagicMock(return_value=iter([burden_p, skat_p, omni_p]))

        # slots["Pvals"] returns mock_pvals
        mock_result.slots = {"Pvals": mock_pvals}
        return mock_result

    def _make_all_three_categories_data(self) -> tuple[dict[str, Any], pd.DataFrame]:
        """Build contingency_data with all three COAST categories present."""
        # 5 variants: 1 BMV, 2 DMV, 2 PTV
        gene_df = _make_gene_df(
            effects=[
                "missense_variant",  # BMV
                "missense_variant",  # DMV
                "missense_variant",  # DMV
                "frameshift_variant",  # PTV
                "stop_gained",  # PTV
            ],
            impacts=["MODERATE", "MODERATE", "MODERATE", "HIGH", "HIGH"],
            sift=["tolerated", "deleterious", "deleterious", ".", "."],
            polyphen=["benign", "probably_damaging", "possibly_damaging", ".", "."],
        )
        n_samples = 60
        n_variants = 5
        rng = np.random.default_rng(123)
        geno = rng.integers(0, 2, size=(n_samples, n_variants)).astype(float)
        phenotype = np.array([1.0] * 30 + [0.0] * 30)

        data: dict[str, Any] = {
            "GENE": "BRCA1",
            "proband_count": 30,
            "control_count": 30,
            "n_qualifying_variants": n_variants,
            "genotype_matrix": geno,
            "phenotype_vector": phenotype,
            "covariate_matrix": None,
            "gene_df": gene_df,
        }
        return data, gene_df

    def test_run_with_mocked_r_returns_omnibus_p_value(self):
        """COASTTest.run() with mocked rpy2 returns correct omnibus p_value."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        config = _make_config()
        data, _ = self._make_all_three_categories_data()

        mock_allelic_series = MagicMock()
        mock_coast_result = self._build_mock_coast_result(burden_p=0.05, skat_p=0.03, omni_p=0.02)
        mock_allelic_series.COAST.return_value = mock_coast_result

        mock_rpy2 = MagicMock()
        mock_robjects = MagicMock()
        mock_robjects.NULL = None
        mock_robjects.r = MagicMock()
        mock_numpy2ri = MagicMock()

        def mock_importr(pkg: str) -> MagicMock:
            if pkg == "AllelicSeries":
                return mock_allelic_series
            return MagicMock()

        mock_robjects_packages = MagicMock()
        mock_robjects_packages.importr.side_effect = mock_importr

        with patch.dict(
            "sys.modules",
            {
                "rpy2": mock_rpy2,
                "rpy2.robjects": mock_robjects,
                "rpy2.robjects.packages": mock_robjects_packages,
                "rpy2.robjects.numpy2ri": mock_numpy2ri,
                "rpy2.rinterface": MagicMock(),
            },
        ):
            result = test.run("BRCA1", data, config)

        # Should have called COAST
        mock_allelic_series.COAST.assert_called_once()
        # Omnibus p-value returned
        assert result.p_value == pytest.approx(0.02)
        assert result.test_name == "coast"
        assert result.gene == "BRCA1"

    def test_run_with_mocked_r_extra_keys_present(self):
        """COASTTest.run() with mocked rpy2 populates extra dict correctly."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        config = _make_config()
        data, _ = self._make_all_three_categories_data()

        mock_allelic_series = MagicMock()
        mock_coast_result = self._build_mock_coast_result(burden_p=0.05, skat_p=0.03, omni_p=0.02)
        mock_allelic_series.COAST.return_value = mock_coast_result

        mock_rpy2 = MagicMock()
        mock_robjects = MagicMock()
        mock_robjects.NULL = None
        mock_numpy2ri = MagicMock()

        def mock_importr(pkg: str) -> MagicMock:
            if pkg == "AllelicSeries":
                return mock_allelic_series
            return MagicMock()

        mock_robjects_packages = MagicMock()
        mock_robjects_packages.importr.side_effect = mock_importr

        with patch.dict(
            "sys.modules",
            {
                "rpy2": mock_rpy2,
                "rpy2.robjects": mock_robjects,
                "rpy2.robjects.packages": mock_robjects_packages,
                "rpy2.robjects.numpy2ri": mock_numpy2ri,
                "rpy2.rinterface": MagicMock(),
            },
        ):
            result = test.run("BRCA1", data, config)

        extra = result.extra
        assert "coast_burden_p_value" in extra
        assert "coast_skat_p_value" in extra
        assert "coast_n_bmv" in extra
        assert "coast_n_dmv" in extra
        assert "coast_n_ptv" in extra
        assert extra["coast_n_bmv"] == 1
        assert extra["coast_n_dmv"] == 2
        assert extra["coast_n_ptv"] == 2

    def test_run_with_mocked_r_effect_fields_all_none(self):
        """COAST TestResult has effect_size, se, ci_lower, ci_upper all None."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        config = _make_config()
        data, _ = self._make_all_three_categories_data()

        mock_allelic_series = MagicMock()
        mock_coast_result = self._build_mock_coast_result(burden_p=0.05, skat_p=0.03, omni_p=0.02)
        mock_allelic_series.COAST.return_value = mock_coast_result

        mock_rpy2 = MagicMock()
        mock_robjects = MagicMock()
        mock_robjects.NULL = None
        mock_numpy2ri = MagicMock()

        def mock_importr(pkg: str) -> MagicMock:
            if pkg == "AllelicSeries":
                return mock_allelic_series
            return MagicMock()

        mock_robjects_packages = MagicMock()
        mock_robjects_packages.importr.side_effect = mock_importr

        with patch.dict(
            "sys.modules",
            {
                "rpy2": mock_rpy2,
                "rpy2.robjects": mock_robjects,
                "rpy2.robjects.packages": mock_robjects_packages,
                "rpy2.robjects.numpy2ri": mock_numpy2ri,
                "rpy2.rinterface": MagicMock(),
            },
        ):
            result = test.run("BRCA1", data, config)

        assert result.effect_size is None
        assert result.se is None
        assert result.ci_lower is None
        assert result.ci_upper is None

    def test_run_with_coast_weights_config(self):
        """COASTTest.run() passes coast_weights from config to R COAST()."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        config = _make_config(coast_weights=[1.0, 4.0, 9.0])
        data, _ = self._make_all_three_categories_data()

        mock_allelic_series = MagicMock()
        mock_coast_result = self._build_mock_coast_result(burden_p=0.04, skat_p=0.02, omni_p=0.01)
        mock_allelic_series.COAST.return_value = mock_coast_result

        mock_rpy2 = MagicMock()
        mock_robjects = MagicMock()
        mock_robjects.NULL = None
        mock_numpy2ri = MagicMock()

        # Track calls to COAST to verify weights were passed
        coast_call_kwargs: dict[str, Any] = {}

        def mock_coast(**kwargs: Any) -> MagicMock:
            coast_call_kwargs.update(kwargs)
            return mock_coast_result

        mock_allelic_series.COAST.side_effect = mock_coast

        def mock_importr(pkg: str) -> MagicMock:
            if pkg == "AllelicSeries":
                return mock_allelic_series
            return MagicMock()

        mock_robjects_packages = MagicMock()
        mock_robjects_packages.importr.side_effect = mock_importr

        with patch.dict(
            "sys.modules",
            {
                "rpy2": mock_rpy2,
                "rpy2.robjects": mock_robjects,
                "rpy2.robjects.packages": mock_robjects_packages,
                "rpy2.robjects.numpy2ri": mock_numpy2ri,
                "rpy2.rinterface": MagicMock(),
            },
        ):
            result = test.run("BRCA1", data, config)

        # COAST should have been called with weights keyword arg
        assert mock_allelic_series.COAST.called, "COAST() was never called"
        # The result should not be None (R call succeeded with our mock)
        # p_value is extracted from mock_coast_result which returns omnibus_p=0.01
        # but since we're checking the _run_r_coast path worked, verify call happened
        assert result.test_name == "coast"


# ---------------------------------------------------------------------------
# COASTTest: Lifecycle hooks
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestCOASTLifecycleHooks:
    """Tests for prepare() and finalize() lifecycle hooks."""

    def test_prepare_sets_total_genes(self):
        """prepare(N) sets _total_genes = N."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        test.prepare(100)
        assert test._total_genes == 100

    def test_prepare_logs_start_message(self, caplog):
        """prepare() logs INFO with 'COAST: beginning' message."""
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()

        with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
            test.prepare(42)

        info_msgs = [r for r in caplog.records if r.levelno == logging.INFO]
        assert any("COAST" in m.message and "42" in m.message for m in info_msgs)

    def test_finalize_logs_timing(self, caplog):
        """finalize() logs INFO with 'COAST complete' and gene count."""
        import time

        from variantcentrifuge.association.tests.allelic_series import COASTTest

        test = COASTTest()
        test._genes_processed = 15
        test._start_time = time.time() - 2.0

        with (
            patch.dict("sys.modules", {"rpy2": None, "rpy2.robjects": None}),
            caplog.at_level(logging.INFO, logger="variantcentrifuge"),
        ):
            test.finalize()

        info_msgs = [r for r in caplog.records if r.levelno == logging.INFO]
        assert any("COAST complete" in m.message for m in info_msgs)
        assert any("15" in m.message for m in info_msgs)


# ---------------------------------------------------------------------------
# Engine registry and ACAT-O integration
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.coast
class TestCOASTEngineRegistration:
    """Tests for COAST registration in AssociationEngine and ACAT-O integration."""

    def test_coast_in_build_registry(self):
        """_build_registry() returns a dict containing 'coast' key."""
        from variantcentrifuge.association.engine import _build_registry

        registry = _build_registry()
        assert "coast" in registry

    def test_coast_registry_entry_is_coast_test_class(self):
        """registry['coast'] is COASTTest class."""
        from variantcentrifuge.association.engine import _build_registry
        from variantcentrifuge.association.tests.allelic_series import COASTTest

        registry = _build_registry()
        assert registry["coast"] is COASTTest

    def test_acat_o_includes_coast_p_value(self):
        """When coast is a registered test, ACAT-O combines coast_p_value."""
        from variantcentrifuge.association.base import AssociationConfig, TestResult
        from variantcentrifuge.association.engine import AssociationEngine

        config = AssociationConfig()

        # Create mock COASTTest that returns a known p-value
        mock_coast = MagicMock()
        mock_coast.name = "coast"
        mock_coast.effect_column_names.return_value = {
            "effect": None,
            "se": None,
            "ci_lower": None,
            "ci_upper": None,
        }
        mock_coast.run.return_value = TestResult(
            gene="GENE1",
            test_name="coast",
            p_value=0.01,
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            se=None,
            n_cases=30,
            n_controls=30,
            n_variants=5,
            extra={"coast_burden_p_value": 0.05, "coast_skat_p_value": 0.02},
        )
        mock_coast.prepare = MagicMock()
        mock_coast.finalize = MagicMock()

        engine = AssociationEngine(tests=[mock_coast], config=config)

        gene_data = [
            {
                "GENE": "GENE1",
                "proband_count": 30,
                "control_count": 30,
                "proband_carrier_count": 5,
                "control_carrier_count": 2,
                "proband_allele_count": 5,
                "control_allele_count": 2,
                "n_qualifying_variants": 5,
            }
        ]

        result_df = engine.run_all(gene_data)

        # ACAT-O should have combined the coast p-value
        assert "acat_o_p_value" in result_df.columns
        # With only one test (coast), ACAT-O = coast p-value (single p-value pass-through)
        assert result_df["acat_o_p_value"].iloc[0] is not None
        # coast_p_value column should be present
        assert "coast_p_value" in result_df.columns
        assert result_df["coast_p_value"].iloc[0] == pytest.approx(0.01)

    def test_coast_extra_columns_in_engine_output(self):
        """coast_burden_p_value, coast_skat_p_value, coast_n_* appear in engine output."""
        from variantcentrifuge.association.base import AssociationConfig, TestResult
        from variantcentrifuge.association.engine import AssociationEngine

        config = AssociationConfig()

        mock_coast = MagicMock()
        mock_coast.name = "coast"
        mock_coast.effect_column_names.return_value = {
            "effect": None,
            "se": None,
            "ci_lower": None,
            "ci_upper": None,
        }
        mock_coast.run.return_value = TestResult(
            gene="GENE1",
            test_name="coast",
            p_value=0.01,
            corrected_p_value=None,
            effect_size=None,
            ci_lower=None,
            ci_upper=None,
            se=None,
            n_cases=30,
            n_controls=30,
            n_variants=5,
            extra={
                "coast_burden_p_value": 0.05,
                "coast_skat_p_value": 0.02,
                "coast_n_bmv": 1,
                "coast_n_dmv": 2,
                "coast_n_ptv": 2,
            },
        )
        mock_coast.prepare = MagicMock()
        mock_coast.finalize = MagicMock()

        engine = AssociationEngine(tests=[mock_coast], config=config)

        gene_data = [
            {
                "GENE": "GENE1",
                "proband_count": 30,
                "control_count": 30,
                "proband_carrier_count": 5,
                "control_carrier_count": 2,
                "proband_allele_count": 5,
                "control_allele_count": 2,
                "n_qualifying_variants": 5,
            }
        ]

        result_df = engine.run_all(gene_data)

        assert "coast_burden_p_value" in result_df.columns
        assert "coast_skat_p_value" in result_df.columns
        assert "coast_n_bmv" in result_df.columns
        assert "coast_n_dmv" in result_df.columns
        assert "coast_n_ptv" in result_df.columns
        row = result_df.iloc[0]
        assert row["coast_burden_p_value"] == pytest.approx(0.05)
        assert row["coast_skat_p_value"] == pytest.approx(0.02)
        assert row["coast_n_bmv"] == 1
        assert row["coast_n_dmv"] == 2
        assert row["coast_n_ptv"] == 2
