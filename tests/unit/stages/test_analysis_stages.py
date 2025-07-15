"""Unit tests for analysis stages."""

from argparse import Namespace
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import (
    ChunkedAnalysisStage,
    CustomAnnotationStage,
    DataFrameLoadingStage,
    GeneBurdenAnalysisStage,
    InheritanceAnalysisStage,
    ParallelAnalysisOrchestrator,
    StatisticsGenerationStage,
    VariantScoringStage,
)


@pytest.fixture
def mock_workspace(tmp_path):
    """Create a mock workspace."""
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    workspace.intermediate_dir = tmp_path / "intermediate"
    workspace.intermediate_dir.mkdir()
    workspace.get_intermediate_path = lambda x: workspace.intermediate_dir / x
    workspace.get_output_path = lambda x, ext=".tsv": workspace.output_dir / f"{x}{ext}"
    return workspace


@pytest.fixture
def base_context(mock_workspace):
    """Create a base pipeline context."""
    return PipelineContext(
        args=Namespace(),
        config={},
        workspace=mock_workspace,
    )


class TestDataFrameLoadingStage:
    """Test the DataFrameLoadingStage."""

    def test_basic_loading(self, base_context, tmp_path):
        """Test basic DataFrame loading from TSV."""
        # Create test file
        test_file = tmp_path / "test.tsv"
        test_file.write_text("col1\tcol2\tcol3\nval1\tval2\tval3\n")

        base_context.config = {"chunk_size": None}
        base_context.data = test_file

        stage = DataFrameLoadingStage()
        result = stage._process(base_context)

        assert result.current_dataframe is not None
        assert isinstance(result.current_dataframe, pd.DataFrame)
        assert len(result.current_dataframe) == 1
        assert list(result.current_dataframe.columns) == ["col1", "col2", "col3"]

    def test_chunked_loading(self, base_context, tmp_path):
        """Test chunked loading for large files."""
        # Create test file
        test_file = tmp_path / "test.tsv"
        test_file.write_text("col1\tcol2\n" + "\n".join(f"val{i}\tdata{i}" for i in range(100)))

        base_context.config = {"chunks": 10}
        base_context.data = test_file

        stage = DataFrameLoadingStage()
        result = stage._process(base_context)

        # Should set chunked processing flag
        assert result.config.get("use_chunked_processing") is True
        assert result.current_dataframe is None

    def test_compressed_file_loading(self, base_context, tmp_path):
        """Test loading gzipped files."""
        import gzip

        test_file = tmp_path / "test.tsv.gz"
        with gzip.open(test_file, "wt") as f:
            f.write("col1\tcol2\nval1\tval2\n")

        base_context.config = {"chunk_size": None}
        base_context.data = test_file

        stage = DataFrameLoadingStage()
        result = stage._process(base_context)

        assert result.current_dataframe is not None
        assert len(result.current_dataframe) == 1


class TestCustomAnnotationStage:
    """Test the CustomAnnotationStage."""

    def test_no_annotations(self, base_context):
        """Test when no annotations are configured."""
        base_context.annotation_configs = {}
        base_context.current_dataframe = pd.DataFrame(
            {"Gene": ["BRCA1", "TP53"], "Variant": ["c.100A>T", "c.200G>A"]}
        )

        stage = CustomAnnotationStage()
        result = stage._process(base_context)

        # Should not modify DataFrame when no annotations
        assert result.current_dataframe.equals(base_context.current_dataframe)

    @patch("variantcentrifuge.stages.analysis_stages.annotate_dataframe_with_features")
    @patch("variantcentrifuge.stages.analysis_stages.load_custom_features")
    def test_with_annotations(self, mock_load, mock_annotate, base_context):
        """Test with annotation configurations."""
        base_context.annotation_configs = {
            "bed_files": ["regions.bed"],
            "gene_lists": ["genes.txt"],
        }
        base_context.current_dataframe = pd.DataFrame(
            {"Gene": ["BRCA1", "TP53"], "Variant": ["c.100A>T", "c.200G>A"]}
        )

        # Mock the load and annotation functions
        mock_load.return_value = [{"type": "bed", "data": "mock_data"}]
        mock_result_df = base_context.current_dataframe.copy()
        mock_result_df["Custom_Annotation"] = ["Region1", "InGeneList"]
        mock_annotate.return_value = mock_result_df

        stage = CustomAnnotationStage()
        result = stage._process(base_context)

        assert mock_load.called
        assert mock_annotate.called
        assert "Custom_Annotation" in result.current_dataframe.columns


class TestInheritanceAnalysisStage:
    """Test the InheritanceAnalysisStage."""

    def test_no_pedigree(self, base_context):
        """Test when no pedigree is loaded."""
        base_context.pedigree_data = None
        base_context.vcf_samples = []  # No samples
        base_context.current_dataframe = pd.DataFrame({"Gene": ["BRCA1"], "Variant": ["c.100A>T"]})

        stage = InheritanceAnalysisStage()
        result = stage._process(base_context)

        # Should skip inheritance analysis due to no samples
        assert result.current_dataframe.equals(base_context.current_dataframe)

    @patch("variantcentrifuge.stages.analysis_stages.analyze_inheritance")
    def test_with_pedigree(self, mock_analyze, base_context):
        """Test inheritance analysis with pedigree."""
        base_context.pedigree_data = {"FAM1": {"proband": "SAMPLE1"}}
        base_context.vcf_samples = ["SAMPLE1", "MOTHER", "FATHER"]
        base_context.config = {"inheritance_mode": "simple"}
        base_context.current_dataframe = pd.DataFrame(
            {"Gene": ["BRCA1"], "Variant": ["c.100A>T"], "SAMPLE1": ["0/1"]}
        )

        # Mock the analysis function
        mock_result = base_context.current_dataframe.copy()
        mock_result["Inheritance_Pattern"] = ["de_novo"]
        mock_analyze.return_value = mock_result

        stage = InheritanceAnalysisStage()
        result = stage._process(base_context)

        assert mock_analyze.called
        assert "Inheritance_Pattern" in result.current_dataframe.columns


class TestVariantScoringStage:
    """Test the VariantScoringStage."""

    def test_no_scoring_config(self, base_context):
        """Test when no scoring config is loaded."""
        base_context.scoring_config = None
        base_context.current_dataframe = pd.DataFrame({"Gene": ["BRCA1"], "Variant": ["c.100A>T"]})

        stage = VariantScoringStage()
        result = stage._process(base_context)

        # Should skip scoring
        assert result.current_dataframe.equals(base_context.current_dataframe)

    @patch("variantcentrifuge.stages.analysis_stages.apply_scoring")
    def test_with_scoring(self, mock_scoring, base_context):
        """Test variant scoring with configuration."""
        base_context.scoring_config = {
            "variable_config": {"AF": "allele_frequency"},
            "formula_config": {"Score": "1 - AF"},
        }
        base_context.current_dataframe = pd.DataFrame(
            {"Gene": ["BRCA1"], "allele_frequency": [0.01]}
        )

        # Mock the scoring function
        mock_result = base_context.current_dataframe.copy()
        mock_result["Score"] = [0.99]
        mock_scoring.return_value = mock_result

        stage = VariantScoringStage()
        result = stage._process(base_context)

        assert mock_scoring.called
        assert "Score" in result.current_dataframe.columns


class TestStatisticsGenerationStage:
    """Test the StatisticsGenerationStage."""

    def test_no_stats_flag(self, base_context):
        """Test when no_stats is True."""
        base_context.config = {"no_stats": True}
        base_context.current_dataframe = pd.DataFrame(
            {"Gene": ["BRCA1", "TP53"], "Impact": ["HIGH", "MODERATE"]}
        )

        stage = StatisticsGenerationStage()
        result = stage._process(base_context)

        # Should skip statistics
        assert result.statistics == {}

    @patch("variantcentrifuge.stages.analysis_stages.StatsEngine")
    def test_generate_stats(self, mock_stats_engine, base_context):
        """Test statistics generation."""
        base_context.config = {"no_stats": False}
        base_context.current_dataframe = pd.DataFrame(
            {"Gene": ["BRCA1", "TP53", "BRCA1"], "Impact": ["HIGH", "MODERATE", "HIGH"]}
        )

        # Mock the StatsEngine
        mock_engine = Mock()
        mock_engine.compute.return_value = {
            "total_variants": 3,
            "unique_genes": 2,
            "impact_counts": {"HIGH": 2, "MODERATE": 1},
        }
        mock_stats_engine.return_value = mock_engine

        stage = StatisticsGenerationStage()
        result = stage._process(base_context)

        assert mock_stats_engine.called
        assert result.statistics == mock_engine.compute.return_value


class TestGeneBurdenAnalysisStage:
    """Test the GeneBurdenAnalysisStage."""

    def test_no_burden_analysis(self, base_context):
        """Test when gene burden analysis is not requested."""
        base_context.config = {"perform_gene_burden": False}
        base_context.current_dataframe = pd.DataFrame({"Gene": ["BRCA1"], "Variant": ["c.100A>T"]})

        stage = GeneBurdenAnalysisStage()
        result = stage._process(base_context)

        # Should skip analysis
        assert result.gene_burden_results is None

    @patch("variantcentrifuge.stages.analysis_stages.perform_gene_burden_analysis")
    @patch("variantcentrifuge.helpers.assign_case_control_counts")
    def test_gene_burden_analysis(self, mock_assign_counts, mock_burden, base_context):
        """Test gene burden calculation."""
        base_context.config = {
            "perform_gene_burden": True,
            "case_samples": ["CASE1", "CASE2"],
            "control_samples": ["CTRL1", "CTRL2"],
            "gene_column": "Gene",
        }
        base_context.vcf_samples = ["CASE1", "CASE2", "CTRL1", "CTRL2"]
        base_context.current_dataframe = pd.DataFrame(
            {
                "GENE": ["BRCA1", "BRCA1", "TP53"],
                "CASE1": ["0/1", "0/0", "0/1"],
                "CASE2": ["0/1", "0/1", "0/0"],
                "CTRL1": ["0/0", "0/0", "0/0"],
                "CTRL2": ["0/0", "0/1", "0/0"],
            }
        )

        # Mock assign_case_control_counts to add the necessary columns
        df_with_counts = base_context.current_dataframe.copy()
        df_with_counts["proband_count"] = [2, 2, 2]
        df_with_counts["control_count"] = [2, 2, 2]
        df_with_counts["proband_variant_count"] = [2, 1, 1]
        df_with_counts["control_variant_count"] = [0, 1, 0]
        df_with_counts["proband_allele_count"] = [3, 1, 1]
        df_with_counts["control_allele_count"] = [0, 1, 0]
        df_with_counts["proband_homozygous_count"] = [0, 0, 0]
        df_with_counts["control_homozygous_count"] = [0, 0, 0]
        mock_assign_counts.return_value = df_with_counts

        # Mock the burden calculation
        mock_result = pd.DataFrame(
            {
                "GENE": ["BRCA1", "TP53"],
                "case_carriers": [2, 1],
                "control_carriers": [1, 0],
                "p_value": [0.05, 0.5],
            }
        )
        mock_burden.return_value = mock_result

        stage = GeneBurdenAnalysisStage()
        result = stage._process(base_context)

        assert mock_assign_counts.called
        assert mock_burden.called
        assert result.gene_burden_results is not None


class TestChunkedAnalysisStage:
    """Test the ChunkedAnalysisStage."""

    def test_dataframe_mode(self, base_context):
        """Test when data is already in DataFrame."""
        base_context.current_dataframe = pd.DataFrame(
            {"Gene": ["BRCA1", "TP53"], "Variant": ["c.100A>T", "c.200G>A"]}
        )
        base_context.config = {"use_chunked_processing": False}

        stage = ChunkedAnalysisStage()

        result = stage._process(base_context)

        # Should pass through without chunking
        assert result == base_context  # No modification when not chunked

    def test_chunked_processing(self, base_context, tmp_path):
        """Test chunked file processing."""
        test_file = tmp_path / "large.tsv"
        test_file.write_text("Gene\tVariant\n" + "\n".join(f"GENE{i}\tvar{i}" for i in range(100)))

        base_context.data = test_file
        base_context.current_dataframe = None
        base_context.config = {"use_chunked_processing": True, "chunks": 10}
        base_context.annotation_configs = {}
        base_context.scoring_config = None
        base_context.pedigree_data = None
        base_context.vcf_samples = []

        stage = ChunkedAnalysisStage()

        # The chunked processing stage is complex and processes in chunks
        # For this test, we verify it runs without error
        result = stage._process(base_context)

        # In the actual implementation, chunked processing happens
        # and results are written to file, not stored in DataFrame
        assert result == base_context  # Context is returned as-is


class TestParallelAnalysisOrchestrator:
    """Test the ParallelAnalysisOrchestrator."""

    def test_single_threaded(self, base_context):
        """Test single-threaded execution."""
        base_context.config = {"threads": 1}
        base_context.current_dataframe = pd.DataFrame(
            {"Gene": ["BRCA1", "TP53"], "Variant": ["c.100A>T", "c.200G>A"]}
        )

        stage = ParallelAnalysisOrchestrator()

        result = stage._process(base_context)

        # Should process without error (simplified implementation)
        assert result == base_context

    def test_parallel_execution(self, base_context):
        """Test parallel analysis execution."""
        base_context.config = {"threads": 4}
        base_context.current_dataframe = pd.DataFrame(
            {"Gene": ["BRCA1", "TP53", "KRAS", "EGFR"], "Variant": ["var1", "var2", "var3", "var4"]}
        )

        stage = ParallelAnalysisOrchestrator()

        result = stage._process(base_context)

        # Should process without error (simplified implementation)
        assert result == base_context
