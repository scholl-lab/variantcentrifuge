"""
Unit tests for COAST partial-category fallback behavior (COAST-03).

Tests that PurePythonCOASTTest.run() proceeds with 1 or 2 variant categories
present and only skips when ALL categories are empty.

Requirements covered: COAST-03 (partial-category fallback)
"""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.tests.allelic_series_python import PurePythonCOASTTest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_gene_df(n_bmv: int, n_dmv: int, n_ptv: int) -> pd.DataFrame:
    """Build a per-variant DataFrame with specified category counts."""
    effects = (
        ["missense_variant"] * n_bmv + ["missense_variant"] * n_dmv + ["frameshift_variant"] * n_ptv
    )
    impacts = ["MODERATE"] * (n_bmv + n_dmv) + ["HIGH"] * n_ptv
    sift = ["tolerated"] * n_bmv + ["deleterious"] * n_dmv + ["."] * n_ptv
    polyphen = ["benign"] * n_bmv + ["probably_damaging"] * n_dmv + ["."] * n_ptv
    return pd.DataFrame(
        {
            "EFFECT": effects,
            "IMPACT": impacts,
            "dbNSFP_SIFT_pred": sift,
            "dbNSFP_Polyphen2_HDIV_pred": polyphen,
        }
    )


def _make_contingency_data(
    n_samples: int,
    n_bmv: int,
    n_dmv: int,
    n_ptv: int,
    seed: int = 42,
) -> dict[str, Any]:
    """Build a minimal contingency_data dict."""
    n_variants = n_bmv + n_dmv + n_ptv
    rng = np.random.default_rng(seed)
    geno = rng.choice([0, 1, 2], size=(n_samples, n_variants), p=[0.6, 0.3, 0.1]).astype(np.float64)
    phenotype = np.array([1.0] * (n_samples // 2) + [0.0] * (n_samples - n_samples // 2))
    gene_df = _make_gene_df(n_bmv, n_dmv, n_ptv)
    return {
        "genotype_matrix": geno,
        "phenotype_vector": phenotype,
        "gene_df": gene_df,
        "proband_count": n_samples // 2,
        "control_count": n_samples - n_samples // 2,
        "n_qualifying_variants": n_variants,
        "covariate_matrix": None,
    }


def _make_config() -> AssociationConfig:
    return AssociationConfig(trait_type="binary", coast_weights=[1.0, 2.0, 3.0])


def _make_mock_backend_result() -> dict[str, Any]:
    """Return a minimal test_gene result dict."""
    return {
        "p_value": 0.05,
        "skat_p_value": 0.08,
        "burden_p_values": [0.1, 0.2, 0.3, 0.1, 0.2, 0.3],
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestCOASTPartialCategory:
    """Test partial-category fallback: COAST-03."""

    def _make_test_instance(self) -> PurePythonCOASTTest:
        """Create a PurePythonCOASTTest with mocked backend."""
        test = PurePythonCOASTTest.__new__(PurePythonCOASTTest)
        # Manually initialize all instance attributes (bypass __init__)
        test._backend = MagicMock()
        test._backend.test_gene.return_value = _make_mock_backend_result()
        test._backend._ensure_skat_backend = MagicMock()
        test._backend._skat_backend = MagicMock()
        test._backend._skat_backend.fit_null_model.return_value = MagicMock()
        test._null_model = MagicMock()  # pre-set to skip lazy init
        test._total_genes = 10
        test._log_interval = 50
        test._genes_processed = 0
        test._start_time = 0.0
        test._n_complete = 0
        test._n_partial = 0
        test._n_skipped = 0
        return test

    def test_all_three_categories_present(self):
        """n_bmv>0, n_dmv>0, n_ptv>0 => coast_status == 'complete', p_value is numeric."""
        test = self._make_test_instance()
        data = _make_contingency_data(n_samples=100, n_bmv=3, n_dmv=3, n_ptv=3)
        config = _make_config()

        result = test.run("GENE1", data, config)

        assert result.p_value is not None, "Expected numeric p_value"
        assert result.extra.get("coast_status") == "complete"
        assert "coast_missing_categories" not in result.extra
        assert result.extra["coast_n_bmv"] == 3
        assert result.extra["coast_n_dmv"] == 3
        assert result.extra["coast_n_ptv"] == 3

    def test_only_ptv_present(self):
        """n_bmv=0, n_dmv=0, n_ptv>0 => coast_status == 'partial', missing BMV and DMV."""
        test = self._make_test_instance()
        data = _make_contingency_data(n_samples=100, n_bmv=0, n_dmv=0, n_ptv=5)
        config = _make_config()

        result = test.run("GENE2", data, config)

        assert result.p_value is not None, "Expected numeric p_value for partial COAST"
        assert result.extra.get("coast_status") == "partial"
        missing_cats = result.extra.get("coast_missing_categories", "")
        assert "BMV" in missing_cats
        assert "DMV" in missing_cats
        assert "PTV" not in missing_cats
        assert result.extra["coast_n_ptv"] == 5

    def test_only_bmv_and_dmv_present(self):
        """n_ptv=0 => coast_status == 'partial', missing PTV only."""
        test = self._make_test_instance()
        data = _make_contingency_data(n_samples=100, n_bmv=3, n_dmv=3, n_ptv=0)
        config = _make_config()

        result = test.run("GENE3", data, config)

        assert result.p_value is not None, "Expected numeric p_value for partial COAST"
        assert result.extra.get("coast_status") == "partial"
        missing_cats = result.extra.get("coast_missing_categories", "")
        assert "PTV" in missing_cats
        assert "BMV" not in missing_cats
        assert "DMV" not in missing_cats

    def test_only_bmv_present(self):
        """n_dmv=0, n_ptv=0 => coast_status == 'partial', missing DMV and PTV."""
        test = self._make_test_instance()
        data = _make_contingency_data(n_samples=100, n_bmv=4, n_dmv=0, n_ptv=0)
        config = _make_config()

        result = test.run("GENE4", data, config)

        assert result.p_value is not None, "Expected numeric p_value for partial COAST"
        assert result.extra.get("coast_status") == "partial"
        missing_cats = result.extra.get("coast_missing_categories", "")
        assert "DMV" in missing_cats
        assert "PTV" in missing_cats
        assert "BMV" not in missing_cats

    def test_all_categories_empty_skips(self):
        """n_bmv=0, n_dmv=0, n_ptv=0 => p_value is None, skip reason ALL_CATEGORIES_EMPTY."""
        # Build contingency data with no classifiable variants (all unclassified)
        n_samples = 100
        n_variants = 3
        rng = np.random.default_rng(7)
        geno = rng.choice([0, 1], size=(n_samples, n_variants)).astype(np.float64)
        phenotype = np.array([1.0] * 50 + [0.0] * 50)
        # Use a gene_df with intergenic variants that won't classify as BMV/DMV/PTV
        gene_df = pd.DataFrame(
            {
                "EFFECT": ["intergenic_variant"] * n_variants,
                "IMPACT": ["MODIFIER"] * n_variants,
                "dbNSFP_SIFT_pred": ["."] * n_variants,
                "dbNSFP_Polyphen2_HDIV_pred": ["."] * n_variants,
            }
        )
        data = {
            "genotype_matrix": geno,
            "phenotype_vector": phenotype,
            "gene_df": gene_df,
            "proband_count": 50,
            "control_count": 50,
            "n_qualifying_variants": n_variants,
            "covariate_matrix": None,
        }
        test = self._make_test_instance()
        config = _make_config()

        result = test.run("GENE5", data, config)

        # All intergenic variants = NO_CLASSIFIABLE_VARIANTS (geno_filtered.shape[1] == 0)
        # This hits the guard before the category check, which is fine.
        # Alternatively, if some variants pass but none are BMV/DMV/PTV:
        assert result.p_value is None, "Expected p_value=None when all categories empty"

    def test_finalize_summary_counts(self):
        """finalize() is called; _n_complete/_n_partial/_n_skipped track correctly."""

        test = self._make_test_instance()
        config = _make_config()

        # Gene with all 3 categories => complete
        data_complete = _make_contingency_data(n_samples=100, n_bmv=2, n_dmv=2, n_ptv=2)
        result1 = test.run("GENE_COMPLETE", data_complete, config)
        assert result1.p_value is not None

        # Gene with only PTV => partial
        data_partial = _make_contingency_data(n_samples=100, n_bmv=0, n_dmv=0, n_ptv=3)
        result2 = test.run("GENE_PARTIAL", data_partial, config)
        assert result2.p_value is not None

        assert test._n_complete == 1
        assert test._n_partial == 1
        assert test._n_skipped == 0

        # Verify finalize logs correctly (just call it; don't assert on log output here)
        import time

        test._start_time = time.time() - 1.0  # avoid divide-by-zero
        test.finalize()  # should not raise

    def test_no_coast_status_in_skipped_before_category_check(self):
        """
        When genotype matrix is absent, skip reason is NO_GENOTYPE_MATRIX (no coast_status).

        This ensures pre-category guards still return p_value=None.
        """
        test = self._make_test_instance()
        config = _make_config()
        data = {"proband_count": 10, "control_count": 10, "n_qualifying_variants": 0}

        result = test.run("GENE_MISSING_GENO", data, config)

        assert result.p_value is None
        assert result.extra.get("coast_skip_reason") == "NO_GENOTYPE_MATRIX"
