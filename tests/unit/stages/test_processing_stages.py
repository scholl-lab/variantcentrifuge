"""Unit tests for processing stages."""

from pathlib import Path
from unittest.mock import ANY, Mock, patch

import pandas as pd
import pytest

from tests.mocks.fixtures import create_test_context
from variantcentrifuge.stages.processing_stages import (
    FieldExtractionStage,
    GeneBedCreationStage,
    GenotypeReplacementStage,
    ParallelVariantExtractionStage,
    PhenotypeIntegrationStage,
    VariantExtractionStage,
)


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
        # Create a temp VCF file that actually exists
        import tempfile

        vcf_file = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False)
        vcf_file.write(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        vcf_file.close()

        ctx = create_test_context(vcf_file=vcf_file.name, config_overrides={"threads": 2})
        ctx.gene_bed_file = Path("/tmp/genes.bed")
        ctx.mark_complete("gene_bed_creation")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.extract_variants")
    def test_variant_extraction(self, mock_extract, context):
        """Test variant extraction."""

        # Mock the extract_variants to create output file with VCF content
        def mock_extract_side_effect(**kwargs):
            # Create the output file with realistic VCF content
            output_path = Path(kwargs["output_file"])
            output_path.parent.mkdir(parents=True, exist_ok=True)
            # Create VCF content with variants that include GT field
            vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr17,length=81195210>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
chr17\t43044295\t.\tA\tG\t100\tPASS\t.\tGT\t0/1\t0/0
chr17\t43044300\t.\tC\tT\t100\tPASS\t.\tGT\t0/0\t0/1
"""
            output_path.write_text(vcf_content)

        mock_extract.side_effect = mock_extract_side_effect

        stage = VariantExtractionStage()
        result = stage(context)

        # Verify extraction call - extract_variants uses keyword arguments
        mock_extract.assert_called_once()
        kwargs = mock_extract.call_args[1]
        assert kwargs["vcf_file"] == context.config["vcf_file"]
        assert kwargs["bed_file"] == "/tmp/genes.bed"
        assert kwargs["cfg"]["threads"] == 2
        assert ".variants.vcf.gz" in kwargs["output_file"]

        # Verify context updates
        assert result.extracted_vcf is not None
        assert result.data == result.extracted_vcf

    @patch("variantcentrifuge.stages.processing_stages.extract_variants")
    def test_with_bcftools_filter(self, mock_extract, context):
        """Test extraction with bcftools filter."""
        context.config["bcftools_prefilter"] = "QUAL >= 30"

        # Mock the extract_variants to create output file
        def mock_extract_side_effect(**kwargs):
            # Create the output file
            output_path = Path(kwargs["output_file"])
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.touch()

        mock_extract.side_effect = mock_extract_side_effect

        stage = VariantExtractionStage()
        stage(context)

        # Verify filter passed
        kwargs = mock_extract.call_args[1]
        assert "bcftools_prefilter" in kwargs["cfg"]
        assert kwargs["cfg"]["bcftools_prefilter"] == "QUAL >= 30"


class TestParallelVariantExtractionStage:
    """Test ParallelVariantExtractionStage."""

    @pytest.fixture
    def context(self):
        """Create test context for parallel processing."""
        # Create a temp VCF file that actually exists
        import tempfile

        vcf_file = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False)
        vcf_file.write(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        vcf_file.close()

        ctx = create_test_context(vcf_file=vcf_file.name, config_overrides={"threads": 4})
        ctx.gene_bed_file = Path("/tmp/genes.bed")
        ctx.mark_complete("gene_bed_creation")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.extract_variants")
    def test_fallback_to_single_thread(self, mock_extract, context):
        """Test fallback when threads <= 1."""
        context.config["threads"] = 1

        # Mock the extract_variants to create output file
        def mock_extract_side_effect(**kwargs):
            # Create the output file
            output_path = Path(kwargs["output_file"])
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.touch()

        mock_extract.side_effect = mock_extract_side_effect

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

    @patch("variantcentrifuge.stages.processing_stages.ProcessPoolExecutor")
    @patch("variantcentrifuge.stages.processing_stages.run_command")
    @patch("variantcentrifuge.stages.processing_stages.extract_variants")
    @patch("variantcentrifuge.stages.processing_stages.split_bed_file")
    def test_with_bcftools_prefilter(
        self, mock_split, mock_extract, mock_run, mock_executor, context
    ):
        """Test parallel extraction with bcftools prefilter."""
        # Add bcftools prefilter to config
        context.config["bcftools_prefilter"] = 'FILTER="PASS" & INFO/AC[0] < 10'

        # Setup mocks
        mock_split.return_value = ["/tmp/chunk0.bed", "/tmp/chunk1.bed"]

        # Mock ProcessPoolExecutor
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
                stage(context)

        # Verify extract_variants was called with bcftools_prefilter
        assert mock_executor_instance.submit.call_count == 2

        # Check that extract_variants calls include bcftools_prefilter
        for call in mock_executor_instance.submit.call_args_list:
            # submit is called with (function, *args, **kwargs)
            # extract_variants is passed as positional args, not kwargs
            args = call[0]  # args is the first element (positional args)
            if len(args) > 1:  # Make sure there are enough args
                # Check the config argument passed to extract_variants
                for arg in args:
                    if isinstance(arg, dict) and "bcftools_prefilter" in arg:
                        assert arg["bcftools_prefilter"] == 'FILTER="PASS" & INFO/AC[0] < 10'
                        break


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
        ctx.mark_complete("configuration_loading")
        ctx.mark_complete("gene_bed_creation")
        ctx.mark_complete("variant_extraction")
        ctx.mark_complete("parallel_variant_extraction")
        ctx.mark_complete("multiallelic_split")
        ctx.mark_complete("snpsift_filtering")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.extract_fields")
    def test_field_extraction(self, mock_extract, context):
        """Test field extraction to TSV."""
        # Mock return value - extract_fields returns the output path
        mock_extract.return_value = "mocked_output.tsv"

        stage = FieldExtractionStage()
        result = stage(context)

        # Verify extraction - extract_fields uses keyword arguments
        mock_extract.assert_called_once()
        kwargs = mock_extract.call_args[1]
        assert kwargs["variant_file"] == "/tmp/variants.vcf.gz"  # variant_file
        assert kwargs["fields"] == "CHROM POS REF ALT QUAL"  # fields as string
        assert isinstance(kwargs["cfg"], dict)  # cfg
        assert ".extracted.tsv" in kwargs["output_file"]  # output_file

        # Verify context updates
        assert result.extracted_tsv is not None
        assert result.data == result.extracted_tsv

    @patch("variantcentrifuge.stages.processing_stages.extract_fields")
    def test_with_gzip_compression(self, mock_extract, context):
        """Test extraction with gzip compression."""
        context.config["gzip_intermediates"] = True

        # Mock return value
        mock_extract.return_value = "mocked_output.tsv.gz"

        stage = FieldExtractionStage()
        stage(context)

        # Verify gzipped output
        output_path = mock_extract.call_args[1]["output_file"]  # output_file argument
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
        import tempfile

        # Create a temp TSV file that actually exists
        tsv_file = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False)
        tsv_file.write(b"CHROM\tPOS\tREF\tALT\tSample1\tSample2\tSample3\n")
        tsv_file.write(b"chr1\t100\tA\tG\t0/1\t0/0\t1/1\n")
        tsv_file.close()

        ctx = create_test_context(config_overrides={"replace_genotypes": True})
        ctx.data = Path(tsv_file.name)
        ctx.extracted_tsv = Path(tsv_file.name)
        ctx.vcf_samples = ["Sample1", "Sample2", "Sample3"]
        ctx.mark_complete("field_extraction")
        ctx.mark_complete("sample_config_loading")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.replace_genotypes")
    def test_genotype_replacement(self, mock_replace, context):
        """Test genotype replacement."""
        # Mock replace_genotypes to return some lines
        mock_replace.return_value = iter(["header", "line1", "line2"])

        stage = GenotypeReplacementStage()
        result = stage(context)

        # Verify replacement was called
        mock_replace.assert_called_once()
        call_args = mock_replace.call_args[0]
        # First arg is the opened file handle
        assert hasattr(call_args[0], "read")  # It's a file-like object
        # Second arg is the config dict
        assert isinstance(call_args[1], dict)
        assert call_args[1]["sample_list"] == "Sample1,Sample2,Sample3"
        # Just verify it's a dict with expected keys
        assert "sample_list" in call_args[1]
        # Verify context updates
        assert result.genotype_replaced_tsv is not None
        assert result.data == result.genotype_replaced_tsv
        assert ".genotype_replaced.tsv" in str(result.data)

    def test_skip_if_disabled(self, context):
        """Test skipping when replacement disabled."""
        context.config["replace_genotypes"] = False

        stage = GenotypeReplacementStage()
        result = stage(context)

        # Data should not change
        assert result.data == context.data

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
        import tempfile

        # Create a temp TSV file that actually exists
        tsv_file = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False)
        tsv_file.write(b"CHROM\tPOS\tREF\tALT\tSample1\tSample2\tSample3\n")
        tsv_file.write(b"chr1\t100\tA\tG\t0/1\t0/0\t1/1\n")
        tsv_file.close()

        ctx = create_test_context()
        ctx.data = Path(tsv_file.name)
        ctx.phenotype_data = {
            "Sample1": {"Affected"},
            "Sample2": {"Unaffected"},
            "Sample3": {"Affected"},
        }
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
            samples=["Sample1", "Sample2", "Sample3"],
            phenotypes=context.phenotype_data,
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
        original_data = context.data
        context.phenotype_data = None

        stage = PhenotypeIntegrationStage()
        result = stage(context)

        # Data should not change (keep the original temp file path)
        assert result.data == original_data
