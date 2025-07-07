"""Tests for all analysis stages."""

import json
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pandas as pd
import pytest

from variantcentrifuge.pipeline_core import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import (
    DataFrameLoadingStage,
    CustomAnnotationStage,
    InheritanceAnalysisStage,
    VariantScoringStage,
    StatisticsGenerationStage,
    GeneBurdenAnalysisStage,
    ChunkedAnalysisStage,
)


class TestDataFrameLoadingStage:
    """Test DataFrameLoadingStage."""

    def test_basic_properties(self):
        """Test stage basic properties."""
        stage = DataFrameLoadingStage()
        assert stage.name == "dataframe_loading"
        assert "Load data into DataFrame" in stage.description
        assert stage.dependencies == {"field_extraction"}
        assert stage.parallel_safe is False

    def test_process_with_tsv_file(self, tmp_path):
        """Test processing TSV file into DataFrame."""
        # Create test TSV
        tsv_file = tmp_path / "test.tsv"
        tsv_content = "CHROM\tPOS\tREF\tALT\tQUAL\n1\t100\tA\tT\t30\n2\t200\tG\tC\t40\n"
        tsv_file.write_text(tsv_content)

        # Create context
        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={}, workspace=workspace)
        context.data = tsv_file

        # Process
        stage = DataFrameLoadingStage()
        result = stage._process(context)

        # Verify
        assert result.current_dataframe is not None
        assert len(result.current_dataframe) == 2
        assert list(result.current_dataframe.columns) == ["CHROM", "POS", "REF", "ALT", "QUAL"]

    def test_process_with_gzipped_tsv(self, tmp_path):
        """Test processing gzipped TSV file."""
        import gzip

        # Create gzipped TSV
        tsv_file = tmp_path / "test.tsv.gz"
        tsv_content = b"CHROM\tPOS\tREF\tALT\n1\t100\tA\tT\n"
        with gzip.open(tsv_file, "wb") as f:
            f.write(tsv_content)

        # Create context
        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={}, workspace=workspace)
        context.data = tsv_file

        # Process
        stage = DataFrameLoadingStage()
        result = stage._process(context)

        # Verify
        assert result.current_dataframe is not None
        assert len(result.current_dataframe) == 1

    def test_process_with_chunked_processing(self, tmp_path):
        """Test behavior when chunked processing is enabled."""
        # Create a dummy file
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text("CHROM\tPOS\n1\t100\n")

        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={"chunks": 1000}, workspace=workspace)
        context.data = tsv_file

        stage = DataFrameLoadingStage()
        result = stage._process(context)

        # When chunks is set, it should enable chunked processing
        assert result.config.get("use_chunked_processing") is True
        assert result.current_dataframe is None


class TestCustomAnnotationStage:
    """Test CustomAnnotationStage."""

    def test_basic_properties(self):
        """Test stage basic properties."""
        stage = CustomAnnotationStage()
        assert stage.name == "custom_annotation"
        assert "Apply custom annotations" in stage.description
        assert stage.dependencies == {"dataframe_loading"}
        assert stage.parallel_safe is False

    def test_process_with_bed_annotation(self, tmp_path):
        """Test BED file annotation."""
        # Create test DataFrame
        df = pd.DataFrame(
            {
                "CHROM": ["1", "1", "2"],
                "POS": [100, 200, 300],
                "REF": ["A", "G", "C"],
                "ALT": ["T", "C", "G"],
            }
        )

        # Create context
        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(
            args=Mock(), config={"annotate_bed": ["exon.bed"]}, workspace=workspace
        )
        context.current_dataframe = df
        context.annotation_configs = None  # No configs means it returns early

        # Process
        stage = CustomAnnotationStage()
        result = stage._process(context)

        # Should return context unchanged when no annotation configs
        assert result.current_dataframe.equals(df)

    def test_process_with_gene_list_annotation(self, tmp_path):
        """Test gene list annotation returns early without configs."""
        # Create test DataFrame
        df = pd.DataFrame(
            {
                "CHROM": ["1", "2", "3"],
                "POS": [100, 200, 300],
                "GENE": ["BRCA1", "BRCA2", "TP53"],
            }
        )

        # Create context
        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(
            args=Mock(), config={"annotate_gene_list": ["cancer_genes.txt"]}, workspace=workspace
        )
        context.current_dataframe = df
        context.annotation_configs = None  # No configs

        # Process
        stage = CustomAnnotationStage()
        result = stage._process(context)

        # Should return unchanged
        assert result.current_dataframe.equals(df)

    def test_process_without_annotations(self, tmp_path):
        """Test processing without annotations."""
        df = pd.DataFrame({"CHROM": ["1"], "POS": [100]})

        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={}, workspace=workspace)
        context.current_dataframe = df
        context.annotation_data = {}

        stage = CustomAnnotationStage()
        result = stage._process(context)

        assert result.current_dataframe.equals(df)


class TestInheritanceAnalysisStage:
    """Test InheritanceAnalysisStage."""

    def test_basic_properties(self):
        """Test stage basic properties."""
        stage = InheritanceAnalysisStage()
        assert stage.name == "inheritance_analysis"
        assert "Calculate inheritance patterns" in stage.description
        assert stage.dependencies == {"dataframe_loading", "custom_annotation", "pedigree_loading"}
        assert stage.parallel_safe is False

    @patch("variantcentrifuge.stages.analysis_stages.analyze_inheritance")
    def test_process_with_pedigree(self, mock_analyze, tmp_path):
        """Test inheritance analysis with pedigree."""
        # Create test DataFrame
        df = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": [100, 200],
                "GENE": ["GENE1", "GENE1"],
                "Sample1": ["0/1", "0/0"],
                "Sample2": ["0/0", "0/1"],
                "Sample3": ["0/1", "0/1"],
            }
        )

        # Mock pedigree data
        pedigree_data = {
            "Sample1": {
                "family": "FAM1",
                "individual": "Sample1",
                "father": "Sample2",
                "mother": "Sample3",
                "sex": "M",
                "phenotype": "affected",
            },
            "Sample2": {
                "family": "FAM1",
                "individual": "Sample2",
                "father": "0",
                "mother": "0",
                "sex": "M",
                "phenotype": "unaffected",
            },
            "Sample3": {
                "family": "FAM1",
                "individual": "Sample3",
                "father": "0",
                "mother": "0",
                "sex": "F",
                "phenotype": "unaffected",
            },
        }

        # Mock inheritance analysis result
        mock_analyze.return_value = df.copy()
        mock_analyze.return_value["Inheritance_Pattern"] = ["de_novo", "autosomal_recessive"]

        # Create context
        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(
            args=Mock(), config={"inheritance_mode": "simple"}, workspace=workspace
        )
        context.current_dataframe = df
        context.pedigree_data = pedigree_data
        context.vcf_samples = ["Sample1", "Sample2", "Sample3"]

        # Process
        stage = InheritanceAnalysisStage()
        result = stage._process(context)

        # Verify
        mock_analyze.assert_called_once()
        assert "Inheritance_Pattern" in result.current_dataframe.columns

    def test_process_without_pedigree(self, tmp_path):
        """Test skipping when no pedigree data."""
        df = pd.DataFrame({"CHROM": ["1"], "POS": [100]})

        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={}, workspace=workspace)
        context.current_dataframe = df
        context.pedigree_data = None

        stage = InheritanceAnalysisStage()
        result = stage._process(context)

        assert result.current_dataframe.equals(df)


class TestVariantScoringStage:
    """Test VariantScoringStage."""

    def test_basic_properties(self):
        """Test stage basic properties."""
        stage = VariantScoringStage()
        assert stage.name == "variant_scoring"
        assert "Calculate variant scores" in stage.description
        assert stage.dependencies == {
            "dataframe_loading",
            "custom_annotation",
            "inheritance_analysis",
            "scoring_config_loading",
        }
        assert stage.parallel_safe is False

    @patch("variantcentrifuge.stages.analysis_stages.apply_scoring")
    def test_process_with_scoring_config(self, mock_scoring, tmp_path):
        """Test variant scoring with config."""
        # Create test DataFrame
        df = pd.DataFrame(
            {
                "CHROM": ["1", "2"],
                "POS": [100, 200],
                "QUAL": [30, 40],
                "IMPACT": ["HIGH", "MODERATE"],
            }
        )

        # Mock scoring config
        scoring_config = {
            "variable_assignment": {"qual": "QUAL", "impact": "IMPACT"},
            "formulas": [{"name": "TestScore", "formula": "qual * 2"}],
        }

        # Mock scoring result
        mock_scoring.return_value = df.copy()
        mock_scoring.return_value["TestScore"] = [60, 80]

        # Create context
        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={}, workspace=workspace)
        context.current_dataframe = df
        context.scoring_config = scoring_config

        # Process
        stage = VariantScoringStage()
        result = stage._process(context)

        # Verify
        mock_scoring.assert_called_once()
        assert result.current_dataframe is not None

    def test_process_without_scoring_config(self, tmp_path):
        """Test skipping when no scoring config."""
        df = pd.DataFrame({"CHROM": ["1"], "POS": [100]})

        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={}, workspace=workspace)
        context.current_dataframe = df
        context.scoring_config = None

        stage = VariantScoringStage()
        result = stage._process(context)

        assert result.current_dataframe.equals(df)


class TestStatisticsGenerationStage:
    """Test StatisticsGenerationStage."""

    def test_basic_properties(self):
        """Test stage basic properties."""
        stage = StatisticsGenerationStage()
        assert stage.name == "statistics_generation"
        assert "Generate summary statistics" in stage.description
        assert stage.dependencies == {"dataframe_loading"}
        assert stage.parallel_safe is False

    @patch("variantcentrifuge.stages.analysis_stages.StatsEngine")
    def test_process_generates_statistics(self, mock_stats_engine, tmp_path):
        """Test statistics generation."""
        # Create test DataFrame
        df = pd.DataFrame(
            {
                "CHROM": ["1", "1", "2"],
                "POS": [100, 200, 300],
                "GENE": ["GENE1", "GENE1", "GENE2"],
                "IMPACT": ["HIGH", "MODERATE", "HIGH"],
            }
        )

        # Mock statistics engine and result
        mock_engine = Mock()
        mock_engine.compute.return_value = {
            "total_variants": 3,
            "variants_per_gene": {"GENE1": 2, "GENE2": 1},
            "impact_distribution": {"HIGH": 2, "MODERATE": 1},
        }
        mock_stats_engine.return_value = mock_engine

        # Create context
        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={}, workspace=workspace)
        context.current_dataframe = df

        # Process
        stage = StatisticsGenerationStage()
        result = stage._process(context)

        # Verify
        mock_stats_engine.assert_called_once()
        mock_engine.compute.assert_called_once_with(df)
        assert result.statistics is not None
        assert result.statistics["total_variants"] == 3

    def test_process_with_empty_dataframe(self, tmp_path):
        """Test statistics with empty DataFrame."""
        df = pd.DataFrame()

        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={}, workspace=workspace)
        context.current_dataframe = df

        stage = StatisticsGenerationStage()
        result = stage._process(context)

        # Should still process but with empty stats
        assert result.statistics is not None


class TestGeneBurdenAnalysisStage:
    """Test GeneBurdenAnalysisStage."""

    def test_basic_properties(self):
        """Test stage basic properties."""
        stage = GeneBurdenAnalysisStage()
        assert stage.name == "gene_burden_analysis"
        assert "Perform gene burden analysis" in stage.description
        # Check dependencies - may vary based on configuration
        assert "dataframe_loading" in stage.dependencies
        assert stage.parallel_safe is False
        assert stage.estimated_runtime == 1.0

    @patch("variantcentrifuge.stages.analysis_stages.perform_gene_burden_analysis")
    @patch("variantcentrifuge.stages.analysis_stages.analyze_variants")
    def test_process_with_case_control(self, mock_analyze, mock_burden, tmp_path):
        """Test gene burden analysis with cases and controls."""
        # Create test DataFrame
        df = pd.DataFrame(
            {
                "GENE": ["GENE1", "GENE1", "GENE2"],
                "Sample1": ["0/1", "0/0", "0/1"],
                "Sample2": ["0/0", "0/1", "0/0"],
                "Sample3": ["0/1", "0/1", "0/0"],
                "Sample4": ["0/0", "0/0", "0/1"],
            }
        )

        # Mock analyze variants result
        mock_analyze.return_value = df.copy()

        # Mock burden test result
        burden_result = pd.DataFrame(
            {
                "gene": ["GENE1", "GENE2"],
                "case_count": [2, 1],
                "control_count": [1, 1],
                "p_value": [0.5, 0.8],
                "odds_ratio": [2.0, 1.0],
            }
        )
        mock_burden.return_value = burden_result

        # Create context
        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(
            args=Mock(),
            config={
                "perform_gene_burden": True,
                "gene_burden_mode": "samples",
                "correction_method": "fdr",
                "case_samples": ["Sample1", "Sample3"],
                "control_samples": ["Sample2", "Sample4"],
            },
            workspace=workspace,
        )
        context.current_dataframe = df

        # Process
        stage = GeneBurdenAnalysisStage()
        result = stage._process(context)

        # Verify
        mock_burden.assert_called_once()
        assert result.gene_burden_results is not None
        assert len(result.gene_burden_results) == 2

    def test_process_without_gene_burden(self, tmp_path):
        """Test skipping when gene burden not requested."""
        df = pd.DataFrame({"GENE": ["GENE1"], "Sample1": ["0/1"]})

        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={}, workspace=workspace)
        context.current_dataframe = df

        stage = GeneBurdenAnalysisStage()
        result = stage._process(context)

        assert result.gene_burden_results is None

    def test_process_without_samples(self, tmp_path):
        """Test warning when no case/control samples."""
        df = pd.DataFrame({"GENE": ["GENE1"], "Sample1": ["0/1"]})

        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(
            args=Mock(), config={"perform_gene_burden": True}, workspace=workspace
        )
        context.current_dataframe = df

        stage = GeneBurdenAnalysisStage()
        result = stage._process(context)

        assert result.gene_burden_results is None


class TestChunkedAnalysisStage:
    """Test ChunkedAnalysisStage."""

    def test_basic_properties(self):
        """Test stage basic properties."""
        stage = ChunkedAnalysisStage()
        assert stage.name == "chunked_analysis"
        assert "Process large files in memory-efficient chunks" in stage.description
        # ChunkedAnalysisStage has more dependencies
        expected_deps = {
            "field_extraction",
            "phenotype_integration",
            "annotation_config_loading",
            "scoring_config_loading",
            "pedigree_loading",
        }
        assert stage.dependencies == expected_deps
        assert stage.parallel_safe is False

    def test_process_with_chunked_mode(self, tmp_path):
        """Test chunked processing."""
        # Create test TSV
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text("CHROM\tPOS\n1\t100\n1\t200\n")

        # Create context
        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(
            args=Mock(),
            config={
                "use_chunked_processing": True,
                "chunk_size": 1,
                "output_file": "output.tsv",
            },
            workspace=workspace,
        )
        context.data = tsv_file

        # Process
        stage = ChunkedAnalysisStage()
        result = stage._process(context)

        # Verify
        assert result.config.get("chunked_processing_complete") is True

    def test_process_without_chunked_mode(self, tmp_path):
        """Test skipping when chunked mode disabled."""
        workspace = Workspace(tmp_path, "test")
        context = PipelineContext(args=Mock(), config={}, workspace=workspace)

        stage = ChunkedAnalysisStage()
        result = stage._process(context)

        assert result.config.get("chunked_processing_complete") is None
