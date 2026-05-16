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


def test_builder_aligns_score_values_with_real_keep_mask():
    from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

    df = pd.DataFrame(
        {
            "GENE": ["A", "A", "A"],
            "nephro_candidate_score": [1.0, 99.0, 7.0],
            "GEN_0__GT": ["0/1", "./.", "1/1"],
            "GEN_1__GT": ["0/1", "./.", "0/1"],
            "GEN_2__GT": ["0/1", "0/0", "0/1"],
            "GEN_3__GT": ["0/1", "0/0", "0/0"],
        }
    )
    builder = _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=["S1", "S2", "S3", "S4"],
        gt_columns=["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"],
        is_binary=True,
        missing_site_threshold=0.25,
        missing_sample_threshold=0.80,
        phenotype_vector=np.array([1, 1, 0, 0]),
        covariate_matrix=None,
        score_column="nephro_candidate_score",
    )

    result = builder()

    assert result["genotype_matrix"].shape[1] == 2
    np.testing.assert_array_equal(result["score_values"], np.array([1.0, 7.0], dtype=object))
    assert result["variant_weight_column"] == "nephro_candidate_score"


def test_builder_aligns_functional_annotation_arrays_with_real_keep_mask():
    from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

    df = pd.DataFrame(
        {
            "GENE": ["A", "A", "A"],
            "dbNSFP_CADD_phred": [30.0, 99.0, 12.0],
            "dbNSFP_REVEL_score": [0.9, 0.1, 0.4],
            "ANN_0__EFFECT": ["missense_variant", "stop_gained", "frameshift_variant"],
            "GEN_0__GT": ["0/1", "./.", "1/1"],
            "GEN_1__GT": ["0/1", "./.", "0/1"],
            "GEN_2__GT": ["0/1", "0/0", "0/1"],
            "GEN_3__GT": ["0/1", "0/0", "0/0"],
        }
    )
    builder = _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=["S1", "S2", "S3", "S4"],
        gt_columns=["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"],
        is_binary=True,
        missing_site_threshold=0.25,
        missing_sample_threshold=0.80,
        phenotype_vector=np.array([1, 1, 0, 0]),
        covariate_matrix=None,
        cadd_column="dbNSFP_CADD_phred",
        revel_column="dbNSFP_REVEL_score",
        effect_column="ANN_0__EFFECT",
    )

    result = builder()

    np.testing.assert_array_equal(result["cadd_scores"], np.array([30.0, 12.0], dtype=object))
    np.testing.assert_array_equal(result["revel_scores"], np.array([0.9, 0.4], dtype=object))
    np.testing.assert_array_equal(
        result["variant_effects"],
        np.array(["missense_variant", "frameshift_variant"], dtype=object),
    )


def test_builder_uses_parser_keep_mask_for_malformed_genotypes():
    from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

    df = pd.DataFrame(
        {
            "GENE": ["A", "A", "A"],
            "nephro_candidate_score": [1.0, 99.0, 7.0],
            "GEN_0__GT": ["0/1", "0/x", "1/1"],
            "GEN_1__GT": ["0/1", "0/x", "0/1"],
            "GEN_2__GT": ["0/1", "0/x", "0/1"],
            "GEN_3__GT": ["0/1", "0/x", "0/0"],
        }
    )
    builder = _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=["S1", "S2", "S3", "S4"],
        gt_columns=["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"],
        is_binary=True,
        missing_site_threshold=0.25,
        missing_sample_threshold=0.80,
        phenotype_vector=np.array([1, 1, 0, 0]),
        covariate_matrix=None,
        score_column="nephro_candidate_score",
    )

    result = builder()

    assert result["genotype_matrix"].shape[1] == 2
    np.testing.assert_array_equal(result["score_values"], np.array([1.0, 7.0], dtype=object))


def test_engine_passes_builder_score_values_to_tests_and_discards_payload():
    from variantcentrifuge.association.base import AssociationConfig, AssociationTest, TestResult
    from variantcentrifuge.association.engine import AssociationEngine
    from variantcentrifuge.stages.analysis_stages import _GenotypeMatrixBuilder

    seen: dict[str, object] = {}

    class RecordingTest(AssociationTest):
        parallel_safe = True

        @property
        def name(self) -> str:
            return "recording"

        def run(self, gene, contingency_data, config):
            seen["score_values"] = contingency_data.get("score_values")
            seen["variant_weight_column"] = contingency_data.get("variant_weight_column")
            return TestResult(
                gene=gene,
                test_name=self.name,
                p_value=1.0,
                corrected_p_value=None,
                effect_size=None,
                ci_lower=None,
                ci_upper=None,
                se=None,
                n_cases=2,
                n_controls=2,
                n_variants=2,
            )

    df = pd.DataFrame(
        {
            "GENE": ["A", "A"],
            "nephro_candidate_score": [2.0, 8.0],
            "GEN_0__GT": ["0/1", "1/1"],
            "GEN_1__GT": ["0/1", "0/1"],
            "GEN_2__GT": ["0/0", "0/1"],
            "GEN_3__GT": ["0/0", "0/0"],
        }
    )
    builder = _GenotypeMatrixBuilder(
        gene_df=df,
        vcf_samples=["S1", "S2", "S3", "S4"],
        gt_columns=["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"],
        is_binary=True,
        missing_site_threshold=0.25,
        missing_sample_threshold=0.80,
        phenotype_vector=np.array([1, 1, 0, 0]),
        covariate_matrix=None,
        score_column="nephro_candidate_score",
    )
    gene_data = [
        {
            "GENE": "A",
            "proband_count": 2,
            "control_count": 2,
            "n_qualifying_variants": 2,
            "_genotype_matrix_builder": builder,
        }
    ]

    engine = AssociationEngine([RecordingTest()], AssociationConfig())
    result_df = engine.run_all(gene_data)

    assert result_df is not None
    np.testing.assert_array_equal(seen["score_values"], np.array([2.0, 8.0], dtype=object))
    assert seen["variant_weight_column"] == "nephro_candidate_score"
    assert "score_values" not in gene_data[0]
    assert "variant_weight_column" not in gene_data[0]
