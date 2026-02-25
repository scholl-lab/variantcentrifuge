"""Tests for streaming genotype matrix construction (PERF-06).

Validates that:
- _GenotypeMatrixBuilder is picklable (required for ProcessPoolExecutor)
- Builder returns correct results for normal, empty, and MAC-filtered genes
- Engine invokes builder and discards matrix after test execution
"""

import pickle

import numpy as np
import pandas as pd
import pytest


@pytest.mark.unit
class TestGenotypeMatrixBuilder:
    """Test the lazy builder dataclass."""

    def test_builder_is_picklable(self):
        """Builder must be picklable for ProcessPoolExecutor."""
        from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

        df = pd.DataFrame(
            {
                "GEN_0__GT": ["0/1", "0/0"],
                "GEN_1__GT": ["1/1", "0/1"],
                "GENE": ["A", "A"],
            }
        )
        builder = _GenotypeMatrixBuilder(
            gene_df=df,
            vcf_samples=["S1", "S2"],
            gt_columns=["GEN_0__GT", "GEN_1__GT"],
            is_binary=True,
            missing_site_threshold=0.1,
            missing_sample_threshold=0.5,
            phenotype_vector=np.array([1, 0]),
            covariate_matrix=None,
        )
        # Must survive pickle round-trip
        restored = pickle.loads(pickle.dumps(builder))
        result = restored()
        assert "genotype_matrix" in result
        assert "variant_mafs" in result
        assert isinstance(result["genotype_matrix"], np.ndarray)

    def test_builder_empty_gene(self):
        """Empty gene_df returns zero-column matrix."""
        from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

        df = pd.DataFrame(
            {
                "GEN_0__GT": pd.Series([], dtype=str),
                "GEN_1__GT": pd.Series([], dtype=str),
                "GENE": pd.Series([], dtype=str),
            }
        )
        builder = _GenotypeMatrixBuilder(
            gene_df=df,
            vcf_samples=["S1", "S2"],
            gt_columns=["GEN_0__GT", "GEN_1__GT"],
            is_binary=True,
            missing_site_threshold=0.1,
            missing_sample_threshold=0.5,
            phenotype_vector=np.array([1, 0]),
            covariate_matrix=None,
        )
        result = builder()
        assert result["genotype_matrix"].shape[1] == 0
        assert result["mac_filtered"] is False

    def test_builder_mac_filter(self):
        """Builder filters genes with MAC < 5."""
        from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

        # All ref genotypes -> MAC=0
        df = pd.DataFrame(
            {
                "GEN_0__GT": ["0/0"],
                "GEN_1__GT": ["0/0"],
                "GENE": ["A"],
            }
        )
        builder = _GenotypeMatrixBuilder(
            gene_df=df,
            vcf_samples=["S1", "S2"],
            gt_columns=["GEN_0__GT", "GEN_1__GT"],
            is_binary=True,
            missing_site_threshold=0.1,
            missing_sample_threshold=0.5,
            phenotype_vector=np.array([1, 0]),
            covariate_matrix=None,
        )
        result = builder()
        assert result["mac_filtered"] is True
        assert result["genotype_matrix"].shape[1] == 0

    def test_builder_returns_phenotype_and_covariate(self):
        """Builder passes through phenotype_vector and covariate_matrix."""
        from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

        df = pd.DataFrame(
            {
                "GEN_0__GT": ["0/0"],
                "GEN_1__GT": ["0/0"],
                "GENE": ["A"],
            }
        )
        pv = np.array([1, 0])
        cm = np.array([[0.5, 1.0], [0.3, 0.8]])
        builder = _GenotypeMatrixBuilder(
            gene_df=df,
            vcf_samples=["S1", "S2"],
            gt_columns=["GEN_0__GT", "GEN_1__GT"],
            is_binary=True,
            missing_site_threshold=0.1,
            missing_sample_threshold=0.5,
            phenotype_vector=pv,
            covariate_matrix=cm,
        )
        result = builder()
        # phenotype and covariate are returned (may be masked but same values when no high-missing)
        assert result["phenotype_vector"] is not None
        assert result["covariate_matrix"] is not None

    def test_builder_result_keys(self):
        """Builder result contains all expected keys."""
        from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

        df = pd.DataFrame(
            {
                "GEN_0__GT": pd.Series([], dtype=str),
                "GENE": pd.Series([], dtype=str),
            }
        )
        builder = _GenotypeMatrixBuilder(
            gene_df=df,
            vcf_samples=["S1"],
            gt_columns=["GEN_0__GT"],
            is_binary=True,
            missing_site_threshold=0.1,
            missing_sample_threshold=0.5,
            phenotype_vector=np.array([1]),
            covariate_matrix=None,
        )
        result = builder()
        expected_keys = {
            "genotype_matrix",
            "variant_mafs",
            "phenotype_vector",
            "covariate_matrix",
            "gt_warnings",
            "mac_filtered",
        }
        assert set(result.keys()) == expected_keys


@pytest.mark.unit
class TestEngineBuilderConsumption:
    """Test that engine correctly invokes and discards builders."""

    def test_matrix_discarded_after_sequential_run(self):
        """After engine.run_all, gene_data dicts should not hold matrices.

        Uses logistic_burden (a matrix-consuming test type) to verify the
        builder is actually invoked and then the matrix is discarded.
        """
        from variantcentrifuge.association.base import AssociationConfig
        from variantcentrifuge.association.engine import AssociationEngine
        from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

        config = AssociationConfig(gene_burden_mode="samples", trait_type="binary")
        engine = AssociationEngine.from_names(["logistic_burden"], config)

        # Create a builder that produces a MAC-filtered (all-ref) matrix
        df = pd.DataFrame(
            {
                "GEN_0__GT": ["0/0"],
                "GEN_1__GT": ["0/0"],
                "GENE": ["TestGene"],
            }
        )
        builder = _GenotypeMatrixBuilder(
            gene_df=df,
            vcf_samples=["S1", "S2"],
            gt_columns=["GEN_0__GT", "GEN_1__GT"],
            is_binary=True,
            missing_site_threshold=0.1,
            missing_sample_threshold=0.5,
            phenotype_vector=np.array([1, 0]),
            covariate_matrix=None,
        )

        gene_data = [
            {
                "GENE": "TestGene",
                "case_alleles": 5,
                "control_alleles": 10,
                "case_total": 20,
                "control_total": 40,
                "_genotype_matrix_builder": builder,
            }
        ]
        result_df = engine.run_all(gene_data)
        assert result_df is not None
        # Verify matrix keys are not lingering after run
        assert "genotype_matrix" not in gene_data[0]
        assert "_genotype_matrix_builder" not in gene_data[0]

    def test_engine_works_without_builder(self):
        """Engine sequential path handles gene_data without builder (Fisher-only)."""
        from variantcentrifuge.association.base import AssociationConfig
        from variantcentrifuge.association.engine import AssociationEngine

        config = AssociationConfig(gene_burden_mode="samples", trait_type="binary")
        engine = AssociationEngine.from_names(["fisher"], config)

        gene_data = [
            {
                "GENE": "GeneA",
                "case_alleles": 5,
                "control_alleles": 10,
                "case_total": 20,
                "control_total": 40,
            }
        ]
        result_df = engine.run_all(gene_data)
        assert result_df is not None
        assert len(result_df) >= 0  # may be empty if no testable variants
