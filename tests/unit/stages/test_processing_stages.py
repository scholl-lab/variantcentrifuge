"""Unit tests for processing stages."""

from pathlib import Path
from unittest.mock import Mock, patch, ANY

import pandas as pd
import pytest
from variantcentrifuge.stages.processing_stages import (
    GeneBedCreationStage,
    VariantExtractionStage,
    ParallelVariantExtractionStage,
    FieldExtractionStage,
    GenotypeReplacementStage,
    PhenotypeIntegrationStage,
)
from tests.mocks.fixtures import create_test_context


class TestGeneBedCreationStage:
    """Test GeneBedCreationStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context(
            gene_name="BRCA1 BRCA2",
            config_overrides={
                "reference": "GRCh37.75",
                "interval_expand": 0,
                "add_chr": True,
            },
        )
        # Mark configuration loading as complete
        ctx.mark_complete("configuration_loading")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.get_gene_bed")
    @patch("variantcentrifuge.stages.processing_stages.normalize_genes")
    def test_gene_bed_creation(self, mock_normalize, mock_get_bed, context):
        """Test BED file creation."""
        # Setup mocks
        mock_normalize.return_value = "BRCA1 BRCA2"
        mock_get_bed.return_value = "/tmp/genes.bed"

        stage = GeneBedCreationStage()
        result = stage(context)

        # Verify normalization
        mock_normalize.assert_called_once_with("BRCA1 BRCA2", None, ANY)

        # Verify BED creation
        mock_get_bed.assert_called_once_with(
            reference="GRCh37.75",
            gene_name="BRCA1 BRCA2",
            interval_expand=0,
            add_chr=True,
            output_dir=str(context.workspace.output_dir),
        )

        # Verify context updates
        assert result.config["normalized_genes"] == "BRCA1 BRCA2"
        assert result.gene_bed_file == Path("/tmp/genes.bed")

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = GeneBedCreationStage()
        assert stage.dependencies == {"configuration_loading"}


class TestVariantExtractionStage:
    """Test VariantExtractionStage."""

    @pytest.fixture
    def context(self):
        """Create test context with gene bed."""
        ctx = create_test_context(vcf_file="/tmp/test.vcf", config_overrides={"threads": 2})
        ctx.gene_bed_file = Path("/tmp/genes.bed")
        ctx.mark_complete("gene_bed_creation")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.extract_variants")
    def test_variant_extraction(self, mock_extract, context):
        """Test variant extraction."""
        stage = VariantExtractionStage()
        result = stage(context)

        # Verify extraction call
        mock_extract.assert_called_once()
        call_args = mock_extract.call_args[1]
        assert call_args["vcf_file"] == "/tmp/test.vcf"
        assert call_args["bed_file"] == "/tmp/genes.bed"
        assert call_args["threads"] == 2
        assert ".variants.vcf.gz" in call_args["output_file"]

        # Verify context updates
        assert result.extracted_vcf is not None
        assert result.data == result.extracted_vcf

    @patch("variantcentrifuge.stages.processing_stages.extract_variants")
    def test_with_bcftools_filter(self, mock_extract, context):
        """Test extraction with bcftools filter."""
        context.config["bcftools_filter"] = "QUAL >= 30"

        stage = VariantExtractionStage()
        stage(context)

        # Verify filter passed
        call_args = mock_extract.call_args[1]
        assert call_args["bcftools_filter"] == "QUAL >= 30"


class TestParallelVariantExtractionStage:
    """Test ParallelVariantExtractionStage."""

    @pytest.fixture
    def context(self):
        """Create test context for parallel processing."""
        ctx = create_test_context(vcf_file="/tmp/test.vcf", config_overrides={"threads": 4})
        ctx.gene_bed_file = Path("/tmp/genes.bed")
        ctx.mark_complete("gene_bed_creation")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.extract_variants")
    def test_fallback_to_single_thread(self, mock_extract, context):
        """Test fallback when threads <= 1."""
        context.config["threads"] = 1

        stage = ParallelVariantExtractionStage()
        stage(context)

        # Should call extract_variants once (not parallel)
        assert mock_extract.call_count == 1

    @patch("variantcentrifuge.stages.processing_stages.ProcessPoolExecutor")
    @patch("variantcentrifuge.stages.processing_stages.run_command")
    @patch("variantcentrifuge.stages.processing_stages.extract_variants")
    @patch("variantcentrifuge.stages.processing_stages.split_bed_file")
    def test_parallel_extraction(self, mock_split, mock_extract, mock_run, mock_executor, context):
        """Test parallel extraction with multiple chunks."""
        # Setup mocks
        mock_split.return_value = ["/tmp/chunk0.bed", "/tmp/chunk1.bed"]

        # Mock ProcessPoolExecutor to avoid pickling issues
        mock_executor_instance = Mock()
        mock_executor.return_value.__enter__ = Mock(return_value=mock_executor_instance)
        mock_executor.return_value.__exit__ = Mock(return_value=None)

        # Mock futures
        future1 = Mock()
        future1.result.return_value = None
        future2 = Mock()
        future2.result.return_value = None
        mock_executor_instance.submit.side_effect = [future1, future2]

        with patch("variantcentrifuge.stages.processing_stages.shutil.move"):
            with patch(
                "variantcentrifuge.stages.processing_stages.as_completed",
                return_value=[future1, future2],
            ):
                stage = ParallelVariantExtractionStage()
                result = stage(context)

        # Verify BED splitting
        mock_split.assert_called_once()

        # Verify parallel extraction submissions
        assert mock_executor_instance.submit.call_count == 2

        # Verify output
        assert result.extracted_vcf is not None
        assert result.data == result.extracted_vcf

    def test_stage_properties(self):
        """Test stage properties."""
        stage = ParallelVariantExtractionStage()
        assert not stage.parallel_safe  # Manages own parallelism
        assert stage.estimated_runtime == 30.0


class TestFieldExtractionStage:
    """Test FieldExtractionStage."""

    @pytest.fixture
    def context(self):
        """Create test context with extracted VCF."""
        ctx = create_test_context(
            config_overrides={
                "extract": ["CHROM", "POS", "REF", "ALT", "QUAL"],
                "gzip_intermediates": False,
            }
        )
        ctx.data = Path("/tmp/variants.vcf.gz")
        # Mark all potential dependencies as complete
        ctx.mark_complete("variant_extraction")
        ctx.mark_complete("parallel_variant_extraction")
        ctx.mark_complete("multiallelic_split")
        ctx.mark_complete("snpsift_filtering")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.extract_fields")
    def test_field_extraction(self, mock_extract, context):
        """Test field extraction to TSV."""
        stage = FieldExtractionStage()
        result = stage(context)

        # Verify extraction
        mock_extract.assert_called_once()
        call_args = mock_extract.call_args[0]
        assert call_args[0] == "/tmp/variants.vcf.gz"
        assert call_args[1] == ["CHROM", "POS", "REF", "ALT", "QUAL"]
        assert ".extracted.tsv" in call_args[2]

        # Verify split_multi_sample_columns is True
        assert mock_extract.call_args[1]["split_multi_sample_columns"] is True

        # Verify context updates
        assert result.extracted_tsv is not None
        assert result.data == result.extracted_tsv

    @patch("variantcentrifuge.stages.processing_stages.extract_fields")
    def test_with_gzip_compression(self, mock_extract, context):
        """Test extraction with gzip compression."""
        context.config["gzip_intermediates"] = True

        stage = FieldExtractionStage()
        stage(context)

        # Verify gzipped output
        output_path = mock_extract.call_args[0][2]
        assert output_path.endswith(".gz")

    def test_no_fields_error(self, context):
        """Test error when no fields specified."""
        context.config["extract"] = []

        stage = FieldExtractionStage()
        with pytest.raises(ValueError, match="No fields specified"):
            stage(context)


class TestGenotypeReplacementStage:
    """Test GenotypeReplacementStage."""

    @pytest.fixture
    def context(self):
        """Create test context with samples."""
        ctx = create_test_context(
            config_overrides={"replace_genotypes": True, "missing_string": "NA"}
        )
        ctx.data = Path("/tmp/extracted.tsv")
        ctx.vcf_samples = ["Sample1", "Sample2", "Sample3"]
        ctx.mark_complete("field_extraction")
        ctx.mark_complete("sample_config_loading")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.replace_genotypes")
    def test_genotype_replacement(self, mock_replace, context):
        """Test genotype replacement."""
        stage = GenotypeReplacementStage()
        result = stage(context)

        # Verify replacement
        mock_replace.assert_called_once()
        call_args = mock_replace.call_args[0]
        assert call_args[0] == "/tmp/extracted.tsv"
        assert call_args[1] == ["Sample1", "Sample2", "Sample3"]
        assert ".genotype_replaced.tsv" in call_args[2]

        # Verify options
        kwargs = mock_replace.call_args[1]
        assert kwargs["missing_string"] == "NA"
        assert kwargs["output_to_stdout"] is False

        # Verify context update
        assert result.data.name.endswith(".genotype_replaced.tsv")

    def test_skip_if_disabled(self, context):
        """Test skipping when replacement disabled."""
        context.config["replace_genotypes"] = False

        stage = GenotypeReplacementStage()
        result = stage(context)

        # Data should not change
        assert result.data == Path("/tmp/extracted.tsv")

    def test_skip_if_no_samples(self, context):
        """Test skipping when no samples available."""
        context.vcf_samples = []

        with patch("variantcentrifuge.stages.processing_stages.replace_genotypes") as mock:
            stage = GenotypeReplacementStage()
            stage(context)

            # Should not call replace_genotypes
            mock.assert_not_called()


class TestPhenotypeIntegrationStage:
    """Test PhenotypeIntegrationStage."""

    @pytest.fixture
    def context(self):
        """Create test context with phenotype data."""
        ctx = create_test_context()
        ctx.data = Path("/tmp/genotypes.tsv")
        ctx.phenotype_data = {"Sample1": "Affected", "Sample2": "Unaffected", "Sample3": "Affected"}
        ctx.vcf_samples = ["Sample1", "Sample2", "Sample3"]
        ctx.mark_complete("field_extraction")
        ctx.mark_complete("phenotype_loading")
        ctx.mark_complete("genotype_replacement")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.aggregate_phenotypes_for_samples")
    @patch("variantcentrifuge.stages.processing_stages.pd.read_csv")
    def test_phenotype_integration(self, mock_read, mock_aggregate, context):
        """Test phenotype integration."""
        # Setup mocks
        mock_df = Mock(spec=pd.DataFrame)
        mock_df.__setitem__ = Mock()  # Support item assignment
        mock_df.to_csv = Mock()  # Add to_csv method
        mock_read.return_value = mock_df
        mock_aggregate.return_value = ["Affected", "Unaffected", "Affected"]

        stage = PhenotypeIntegrationStage()
        stage(context)

        # Verify aggregation
        mock_aggregate.assert_called_once_with(
            phenotypes=context.phenotype_data,
            samples=["Sample1", "Sample2", "Sample3"],
            missing_string="",
        )

        # Verify phenotype column added
        assert mock_df.__setitem__.called
        call_args = mock_df.__setitem__.call_args[0]
        assert call_args[0] == "Phenotypes"
        assert call_args[1] == ["Affected", "Unaffected", "Affected"]

        # Verify output written
        mock_df.to_csv.assert_called_once()

    def test_skip_if_no_phenotypes(self, context):
        """Test skipping when no phenotype data."""
        context.phenotype_data = None

        stage = PhenotypeIntegrationStage()
        result = stage(context)

        # Data should not change
        assert result.data == Path("/tmp/genotypes.tsv")
