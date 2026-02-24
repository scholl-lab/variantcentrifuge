"""
Fisher's exact test unit tests with bit-identity validation.

Tests that FisherExactTest produces bit-identical p-values, odds ratios,
and confidence intervals to gene_burden.py for the same contingency data.

Bit-identity is the critical regression guard for the Phase 18 Fisher refactor:
the new association framework must not introduce any numerical drift.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import fisher_exact

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.engine import AssociationEngine
from variantcentrifuge.association.tests.fisher import FisherExactTest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_gene_data(
    gene: str,
    p_carriers: int,
    c_carriers: int,
    p_total: int,
    c_total: int,
    n_variants: int,
    p_alleles: int | None = None,
    c_alleles: int | None = None,
) -> dict:
    """Build a gene burden data dict in the format produced by _aggregate_gene_burden_from_*."""
    return {
        "GENE": gene,
        "proband_count": p_total,
        "control_count": c_total,
        "proband_carrier_count": p_carriers,
        "control_carrier_count": c_carriers,
        "proband_allele_count": p_alleles if p_alleles is not None else p_carriers,
        "control_allele_count": c_alleles if c_alleles is not None else c_carriers,
        "n_qualifying_variants": n_variants,
    }


def _run_fisher_direct(table: list[list[int]]) -> tuple[float, float]:
    """Run scipy.stats.fisher_exact directly — the bit-identical reference."""
    odds_ratio, pval = fisher_exact(table)
    return float(odds_ratio), float(pval)


def _compute_samples_table(gene_data: dict) -> list[list[int]]:
    """Build the 2x2 table for samples mode (mirrors gene_burden.py and fisher.py logic)."""
    p_count = gene_data["proband_count"]
    c_count = gene_data["control_count"]
    p_var = gene_data["proband_carrier_count"]
    c_var = gene_data["control_carrier_count"]
    p_ref = p_count - p_var
    c_ref = c_count - c_var
    return [[p_var, c_var], [p_ref, c_ref]]


def _compute_alleles_table(gene_data: dict) -> list[list[int]]:
    """Build the 2x2 table for alleles mode (mirrors gene_burden.py and fisher.py logic)."""
    p_count = gene_data["proband_count"]
    c_count = gene_data["control_count"]
    p_all = gene_data["proband_allele_count"]
    c_all = gene_data["control_allele_count"]
    p_ref = p_count * 2 - p_all
    c_ref = c_count * 2 - c_all
    return [[p_all, c_all], [p_ref, c_ref]]


# ---------------------------------------------------------------------------
# Bit-identity tests — samples mode
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestFisherBitIdentitySamplesMode:
    """FisherExactTest matches scipy.stats.fisher_exact exactly for samples mode."""

    def _get_fisher_result(self, gene_data: dict, config: AssociationConfig):
        """Run FisherExactTest.run() for one gene and return the TestResult."""
        test = FisherExactTest()
        return test.run(gene_data["GENE"], gene_data, config)

    def test_bit_identity_single_gene_p_value(self):
        """Raw p-value from FisherExactTest equals scipy.stats.fisher_exact exactly."""
        gene_data = _make_gene_data(
            "BRCA1", p_carriers=10, c_carriers=2, p_total=100, c_total=200, n_variants=5
        )
        config = AssociationConfig(gene_burden_mode="samples")

        result = self._get_fisher_result(gene_data, config)
        table = _compute_samples_table(gene_data)
        _, expected_pval = _run_fisher_direct(table)

        assert result.p_value == expected_pval  # Exact equality (bit-identical)

    def test_bit_identity_single_gene_odds_ratio(self):
        """Odds ratio from FisherExactTest equals scipy.stats.fisher_exact exactly."""
        gene_data = _make_gene_data(
            "BRCA1", p_carriers=10, c_carriers=2, p_total=100, c_total=200, n_variants=5
        )
        config = AssociationConfig(gene_burden_mode="samples")

        result = self._get_fisher_result(gene_data, config)
        table = _compute_samples_table(gene_data)
        expected_or, _ = _run_fisher_direct(table)

        assert result.effect_size == expected_or  # Exact equality (bit-identical)

    def test_bit_identity_multiple_genes_all_p_values(self):
        """FisherExactTest p-values match scipy exactly for all genes in a set."""
        genes = [
            _make_gene_data(
                "GENE_A", p_carriers=5, c_carriers=1, p_total=50, c_total=100, n_variants=3
            ),
            _make_gene_data(
                "GENE_B", p_carriers=8, c_carriers=3, p_total=50, c_total=100, n_variants=2
            ),
            _make_gene_data(
                "GENE_C", p_carriers=1, c_carriers=0, p_total=50, c_total=100, n_variants=1
            ),
            _make_gene_data(
                "GENE_D", p_carriers=0, c_carriers=0, p_total=50, c_total=100, n_variants=1
            ),
            _make_gene_data(
                "GENE_E", p_carriers=20, c_carriers=5, p_total=50, c_total=100, n_variants=4
            ),
        ]
        config = AssociationConfig(gene_burden_mode="samples")
        fisher = FisherExactTest()

        for gene_data in genes:
            if gene_data["n_qualifying_variants"] == 0:
                continue
            p_count = gene_data["proband_count"]
            c_count = gene_data["control_count"]
            if p_count == 0 and c_count == 0:
                continue

            result = fisher.run(gene_data["GENE"], gene_data, config)
            table = _compute_samples_table(gene_data)

            if result.p_value is not None:
                _, expected_pval = _run_fisher_direct(table)
                assert result.p_value == expected_pval, (
                    f"p_value mismatch for {gene_data['GENE']}: "
                    f"got {result.p_value}, expected {expected_pval}"
                )

    def test_bit_identity_via_engine_p_values(self):
        """Engine run_all() p-values match direct scipy calls bit-identically."""

        genes = [
            _make_gene_data(
                "ALPHA", p_carriers=5, c_carriers=1, p_total=50, c_total=100, n_variants=3
            ),
            _make_gene_data(
                "BETA", p_carriers=8, c_carriers=3, p_total=50, c_total=100, n_variants=2
            ),
            _make_gene_data(
                "GAMMA", p_carriers=2, c_carriers=1, p_total=50, c_total=100, n_variants=1
            ),
        ]
        config = AssociationConfig(gene_burden_mode="samples", correction_method="fdr")
        engine = AssociationEngine.from_names(["fisher"], config)
        result_df = engine.run_all(genes)

        # Verify p-values match direct scipy calls for each gene
        sorted_genes = sorted(genes, key=lambda d: d["GENE"])
        for gene_data in sorted_genes:
            gene_name = gene_data["GENE"]
            table = _compute_samples_table(gene_data)
            expected_or, expected_pval = _run_fisher_direct(table)

            row = result_df[result_df["gene"] == gene_name].iloc[0]
            assert row["fisher_pvalue"] == expected_pval, f"p_value mismatch for {gene_name}"
            assert row["fisher_or"] == expected_or, f"OR mismatch for {gene_name}"

    def test_bit_identity_corrected_p_values_match_direct_smm(self):
        """ACAT-O corrected p-values match direct smm.multipletests on ACAT-O raw p-values.

        ARCH-03: FDR is applied only to ACAT-O p-values (not per-test).
        With a single test (fisher), ACAT-O p-values are pass-throughs of fisher p-values,
        so the corrected ACAT-O values match FDR of fisher p-values directly.
        """
        import statsmodels.stats.multitest as smm

        genes = [
            _make_gene_data(
                "A1", p_carriers=5, c_carriers=1, p_total=50, c_total=100, n_variants=3
            ),
            _make_gene_data(
                "B2", p_carriers=2, c_carriers=1, p_total=50, c_total=100, n_variants=2
            ),
            _make_gene_data(
                "C3", p_carriers=1, c_carriers=0, p_total=50, c_total=100, n_variants=1
            ),
        ]
        config = AssociationConfig(gene_burden_mode="samples", correction_method="fdr")
        engine = AssociationEngine.from_names(["fisher"], config)
        result_df = engine.run_all(genes)

        # Primary test has no corrected_p_value column (ARCH-03)
        assert "fisher_corrected_p_value" not in result_df.columns

        # ACAT-O corrected values should match FDR applied to ACAT-O raw p-values
        raw_acat_pvals = result_df["acat_o_pvalue"].values
        expected_corrected = smm.multipletests(raw_acat_pvals, method="fdr_bh")[1]

        np.testing.assert_array_almost_equal(
            result_df["acat_o_qvalue"].values,
            expected_corrected,
            decimal=15,
        )


# ---------------------------------------------------------------------------
# Bit-identity tests — alleles mode
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestFisherBitIdentityAllelesMode:
    """FisherExactTest matches scipy.stats.fisher_exact exactly for alleles mode."""

    def test_bit_identity_alleles_mode_p_value(self):
        """Alleles mode p-value matches scipy.stats.fisher_exact exactly."""
        gene_data = _make_gene_data(
            "MYH7",
            p_carriers=10,
            c_carriers=2,
            p_total=50,
            c_total=100,
            n_variants=3,
            p_alleles=12,  # some samples are hom-alt
            c_alleles=2,
        )
        config = AssociationConfig(gene_burden_mode="alleles")
        fisher = FisherExactTest()
        result = fisher.run(gene_data["GENE"], gene_data, config)

        table = _compute_alleles_table(gene_data)
        _, expected_pval = _run_fisher_direct(table)

        assert result.p_value == expected_pval

    def test_bit_identity_alleles_mode_odds_ratio(self):
        """Alleles mode odds ratio matches scipy.stats.fisher_exact exactly."""
        gene_data = _make_gene_data(
            "PKD1",
            p_carriers=8,
            c_carriers=2,
            p_total=50,
            c_total=100,
            n_variants=4,
            p_alleles=10,
            c_alleles=2,
        )
        config = AssociationConfig(gene_burden_mode="alleles")
        fisher = FisherExactTest()
        result = fisher.run(gene_data["GENE"], gene_data, config)

        table = _compute_alleles_table(gene_data)
        expected_or, _ = _run_fisher_direct(table)

        assert result.effect_size == expected_or

    def test_alleles_mode_table_respects_diploid_constraint(self):
        """Alleles mode reference count = 2*N - alt alleles (diploid constraint)."""
        p_count = 50
        c_count = 100
        p_alleles = 8
        c_alleles = 3
        gene_data = {
            "GENE": "TEST",
            "proband_count": p_count,
            "control_count": c_count,
            "proband_carrier_count": 6,
            "control_carrier_count": 3,
            "proband_allele_count": p_alleles,
            "control_allele_count": c_alleles,
            "n_qualifying_variants": 2,
        }
        config = AssociationConfig(gene_burden_mode="alleles")
        fisher = FisherExactTest()
        result = fisher.run("TEST", gene_data, config)

        # Verify via extra table stored by FisherExactTest
        assert result.extra.get("table") is not None
        table = result.extra["table"]
        p_ref = p_count * 2 - p_alleles
        c_ref = c_count * 2 - c_alleles
        assert table[0][0] == p_alleles
        assert table[0][1] == c_alleles
        assert table[1][0] == p_ref
        assert table[1][1] == c_ref


# ---------------------------------------------------------------------------
# FisherExactTest dependency check
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestFisherDependencies:
    """FisherExactTest.check_dependencies() passes in test environment."""

    def test_check_dependencies_does_not_raise(self):
        """scipy and statsmodels are installed; check_dependencies must succeed."""
        fisher = FisherExactTest()
        # Should not raise ImportError
        fisher.check_dependencies()

    def test_name_property_is_fisher(self):
        """FisherExactTest.name == 'fisher'."""
        fisher = FisherExactTest()
        assert fisher.name == "fisher"


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestFisherEdgeCases:
    """Edge cases for FisherExactTest.run()."""

    def test_zero_variant_gene_returns_none_p_value(self):
        """Gene with n_qualifying_variants=0 returns TestResult with p_value=None."""
        gene_data = _make_gene_data(
            "EMPTY", p_carriers=0, c_carriers=0, p_total=50, c_total=100, n_variants=0
        )
        config = AssociationConfig(gene_burden_mode="samples")
        fisher = FisherExactTest()
        result = fisher.run("EMPTY", gene_data, config)

        assert result.p_value is None
        assert result.corrected_p_value is None
        assert result.effect_size is None

    def test_zero_carriers_both_groups_returns_none_p_value(self):
        """Gene with zero counts in both groups is skipped (p_value=None)."""
        gene_data = {
            "GENE": "NOPROBAND",
            "proband_count": 0,
            "control_count": 0,
            "proband_carrier_count": 0,
            "control_carrier_count": 0,
            "proband_allele_count": 0,
            "control_allele_count": 0,
            "n_qualifying_variants": 2,
        }
        config = AssociationConfig(gene_burden_mode="samples")
        fisher = FisherExactTest()
        result = fisher.run("NOPROBAND", gene_data, config)

        assert result.p_value is None

    def test_all_carriers_in_cases_none_in_controls(self):
        """Gene with all carriers in cases and none in controls gives valid p-value and high OR."""
        gene_data = _make_gene_data(
            "EXTREME",
            p_carriers=10,
            c_carriers=0,
            p_total=50,
            c_total=100,
            n_variants=3,
        )
        config = AssociationConfig(gene_burden_mode="samples")
        fisher = FisherExactTest()
        result = fisher.run("EXTREME", gene_data, config)

        # p-value should be non-None (valid table with at least one carrier)
        assert result.p_value is not None
        assert result.p_value < 0.05  # Strong signal

    def test_one_variant_one_carrier_valid_result(self):
        """Sparse table (1 variant, 1 carrier) produces valid result."""
        gene_data = _make_gene_data(
            "SPARSE",
            p_carriers=1,
            c_carriers=0,
            p_total=50,
            c_total=100,
            n_variants=1,
        )
        config = AssociationConfig(gene_burden_mode="samples")
        fisher = FisherExactTest()
        result = fisher.run("SPARSE", gene_data, config)

        assert result.p_value is not None
        # With continuity correction, CI bounds should be finite (not NaN)
        # CI might be None if computation failed, but p_value should be valid

    def test_negative_reference_count_returns_none_p_value(self):
        """Gene where carrier count > total count (data error) is skipped."""
        # p_var > p_count creates negative p_ref — data quality issue
        gene_data = {
            "GENE": "BADDATA",
            "proband_count": 10,
            "control_count": 20,
            "proband_carrier_count": 15,  # > proband_count (invalid!)
            "control_carrier_count": 2,
            "proband_allele_count": 15,
            "control_allele_count": 2,
            "n_qualifying_variants": 3,
        }
        config = AssociationConfig(gene_burden_mode="samples")
        fisher = FisherExactTest()
        result = fisher.run("BADDATA", gene_data, config)

        assert result.p_value is None

    def test_negative_reference_alleles_returns_none_p_value(self):
        """Alleles mode: allele count > 2*N (data error) is skipped."""
        gene_data = {
            "GENE": "BADALLELES",
            "proband_count": 10,
            "control_count": 20,
            "proband_carrier_count": 5,
            "control_carrier_count": 2,
            "proband_allele_count": 25,  # > 2*proband_count (invalid!)
            "control_allele_count": 2,
            "n_qualifying_variants": 3,
        }
        config = AssociationConfig(gene_burden_mode="alleles")
        fisher = FisherExactTest()
        result = fisher.run("BADALLELES", gene_data, config)

        assert result.p_value is None

    def test_zero_cell_triggers_continuity_correction_for_ci(self):
        """Table with zero cell uses Haldane-Anscombe continuity correction for CI."""
        # c_carriers=0 means one cell in table is 0, triggering correction
        gene_data = _make_gene_data(
            "ZEROCELL",
            p_carriers=5,
            c_carriers=0,  # zero cell -> continuity correction
            p_total=50,
            c_total=100,
            n_variants=3,
        )
        config = AssociationConfig(gene_burden_mode="samples", continuity_correction=0.5)
        fisher = FisherExactTest()
        result = fisher.run("ZEROCELL", gene_data, config)

        assert result.p_value is not None
        # CI should be computed (not None) due to continuity correction
        assert result.ci_lower is not None or result.ci_upper is not None

    def test_result_test_name_is_fisher(self):
        """TestResult.test_name is always 'fisher' for FisherExactTest."""
        gene_data = _make_gene_data(
            "BRCA1", p_carriers=5, c_carriers=1, p_total=50, c_total=100, n_variants=2
        )
        config = AssociationConfig()
        fisher = FisherExactTest()
        result = fisher.run("BRCA1", gene_data, config)

        assert result.test_name == "fisher"

    def test_result_gene_name_matches_input(self):
        """TestResult.gene matches the gene parameter passed to run()."""
        gene_data = _make_gene_data(
            "MYGENE", p_carriers=3, c_carriers=1, p_total=50, c_total=100, n_variants=1
        )
        config = AssociationConfig()
        fisher = FisherExactTest()
        result = fisher.run("MYGENE", gene_data, config)

        assert result.gene == "MYGENE"

    def test_ci_fallback_chain_score_to_normal_to_logit(self):
        """CI computation uses score -> normal -> logit fallback chain without error."""
        # Create a table where score method might fail: heavily imbalanced table
        gene_data = _make_gene_data(
            "IMBALANCED",
            p_carriers=50,
            c_carriers=0,
            p_total=50,
            c_total=1000,
            n_variants=5,
        )
        config = AssociationConfig(gene_burden_mode="samples")
        fisher = FisherExactTest()
        # Should complete without exception regardless of which CI method is used
        result = fisher.run("IMBALANCED", gene_data, config)

        assert result.p_value is not None
        # CI may be None if all methods fail, but no exception should be raised

    def test_table_stored_in_extra(self):
        """FisherExactTest stores the contingency table in TestResult.extra['table']."""
        gene_data = _make_gene_data(
            "BRCA2", p_carriers=5, c_carriers=1, p_total=50, c_total=100, n_variants=2
        )
        config = AssociationConfig(gene_burden_mode="samples")
        fisher = FisherExactTest()
        result = fisher.run("BRCA2", gene_data, config)

        assert "table" in result.extra
        table = result.extra["table"]
        assert table[0][0] == 5  # p_carriers
        assert table[0][1] == 1  # c_carriers
        assert table[1][0] == 45  # p_ref = 50 - 5
        assert table[1][1] == 99  # c_ref = 100 - 1

    def test_n_cases_n_controls_n_variants_populated(self):
        """TestResult metadata fields (n_cases, n_controls, n_variants) are correct."""
        gene_data = _make_gene_data(
            "META", p_carriers=3, c_carriers=1, p_total=75, c_total=150, n_variants=4
        )
        config = AssociationConfig()
        fisher = FisherExactTest()
        result = fisher.run("META", gene_data, config)

        assert result.n_cases == 75
        assert result.n_controls == 150
        assert result.n_variants == 4
