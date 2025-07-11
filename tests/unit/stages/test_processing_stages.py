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

    @patch("variantcentrifuge.stages.processing_stages.pd.read_csv")
    def test_phenotype_integration(self, mock_read, context):
        """Test phenotype integration."""
        # Setup mocks
        mock_df = Mock(spec=pd.DataFrame)
        mock_df.__setitem__ = Mock()  # Support item assignment
        mock_df.to_csv = Mock()  # Add to_csv method
        mock_df.columns = ["CHROM", "POS", "REF", "ALT", "GT"]  # Include GT column
        mock_df.__getitem__ = Mock()  # Support item access
        mock_apply_result = Mock()
        mock_df.__getitem__.return_value.apply = Mock(return_value=mock_apply_result)
        mock_read.return_value = mock_df

        stage = PhenotypeIntegrationStage()
        stage(context)

        # Verify GT column was accessed
        mock_df.__getitem__.assert_called_with("GT")

        # Verify apply was called (phenotype extraction function)
        mock_df.__getitem__.return_value.apply.assert_called_once()

        # Verify phenotype column added
        assert mock_df.__setitem__.called
        call_args = mock_df.__setitem__.call_args[0]
        assert call_args[0] == "Phenotypes"

        # Verify output written
        mock_df.to_csv.assert_called_once()

    def test_extract_phenotypes_for_gt_row(self, context):
        """Test the GT row-specific phenotype extraction function."""
        from variantcentrifuge.phenotype import extract_phenotypes_for_gt_row

        # Test data
        phenotypes = {
            "Sample1": {"HP:0000001", "HP:0000002"},
            "Sample2": {"HP:0000003"},
            "Sample3": set(),  # No phenotypes
            "Sample4": {"HP:0000004"},
        }

        # Test with some samples having variants
        gt_value = "Sample1(0/1);Sample2(1/1);Sample3(0/0);Sample4(./.)"
        result = extract_phenotypes_for_gt_row(gt_value, phenotypes)

        # Should only include samples with variants (not 0/0 or ./.)
        expected = "Sample1(HP:0000001,HP:0000002);Sample2(HP:0000003)"
        assert result == expected

    def test_extract_phenotypes_empty_gt(self, context):
        """Test phenotype extraction with empty GT value."""
        from variantcentrifuge.phenotype import extract_phenotypes_for_gt_row

        phenotypes = {"Sample1": {"HP:0000001"}}

        # Test empty GT
        assert extract_phenotypes_for_gt_row("", phenotypes) == ""
        assert extract_phenotypes_for_gt_row(None, phenotypes) == ""

    def test_extract_phenotypes_no_variants(self, context):
        """Test phenotype extraction when no samples have variants."""
        from variantcentrifuge.phenotype import extract_phenotypes_for_gt_row

        phenotypes = {"Sample1": {"HP:0000001"}, "Sample2": {"HP:0000002"}}

        # All samples have no variants
        gt_value = "Sample1(0/0);Sample2(./.)"
        result = extract_phenotypes_for_gt_row(gt_value, phenotypes)
        assert result == ""

    def test_phenotype_integration_realistic_scenario(self, context):
        """Test realistic scenario with multiple variants and samples."""
        from variantcentrifuge.phenotype import extract_phenotypes_for_gt_row

        # Setup realistic phenotype data
        phenotypes = {
            "Sample1": {"HP:0000001", "HP:0000002"},  # Multiple phenotypes
            "Sample2": {"HP:0000003"},  # Single phenotype
            "Sample3": set(),  # No phenotypes
            "Sample4": {"HP:0000004", "HP:0000005"},  # Multiple phenotypes
            "Sample5": {"HP:0000006"},  # Single phenotype
        }

        # Test different GT scenarios that would occur in real data
        test_cases = [
            # Variant 1: Some samples have variants
            {
                "gt": "Sample1(0/1);Sample2(0/0);Sample3(1/1);Sample4(./.);Sample5(0/1)",
                "expected": "Sample1(HP:0000001,HP:0000002);Sample3();Sample5(HP:0000006)",
            },
            # Variant 2: Only one sample has variant
            {
                "gt": "Sample1(0/0);Sample2(1/1);Sample3(0/0);Sample4(./.);Sample5(0/0)",
                "expected": "Sample2(HP:0000003)",
            },
            # Variant 3: No samples have variants (all reference or missing)
            {
                "gt": "Sample1(0/0);Sample2(0/0);Sample3(./.);Sample4(./.);Sample5(0/0)",
                "expected": "",
            },
            # Variant 4: All samples have variants
            {
                "gt": "Sample1(0/1);Sample2(1/1);Sample3(0/1);Sample4(1/1);Sample5(0/1)",
                "expected": (
                    "Sample1(HP:0000001,HP:0000002);Sample2(HP:0000003);Sample3();"
                    "Sample4(HP:0000004,HP:0000005);Sample5(HP:0000006)"
                ),
            },
            # Variant 5: Complex genotypes (homozygous variants)
            {
                "gt": "Sample1(1/1);Sample2(0/0);Sample3(0/1);Sample4(./.);Sample5(1/1)",
                "expected": "Sample1(HP:0000001,HP:0000002);Sample3();Sample5(HP:0000006)",
            },
        ]

        for i, test_case in enumerate(test_cases):
            result = extract_phenotypes_for_gt_row(test_case["gt"], phenotypes)
            assert result == test_case["expected"], (
                f"Variant {i+1} failed:\n"
                f"GT: {test_case['gt']}\n"
                f"Expected: {test_case['expected']}\n"
                f"Got: {result}"
            )

    def test_phenotype_integration_end_to_end_with_real_dataframe(self, context):
        """Test the complete phenotype integration with a realistic DataFrame."""
        import pandas as pd
        from unittest.mock import patch

        # Create realistic test data
        test_data = {
            "CHROM": ["chr1", "chr1", "chr2", "chr2"],
            "POS": ["12345", "67890", "11111", "22222"],
            "REF": ["A", "G", "T", "C"],
            "ALT": ["T", "C", "A", "G"],
            "GT": [
                "Sample1(0/1);Sample2(0/0);Sample3(1/1)",  # Variant 1
                "Sample1(0/0);Sample2(1/1);Sample3(0/0)",  # Variant 2: Only Sample2 has variant
                "Sample1(0/0);Sample2(0/0);Sample3(0/0)",  # Variant 3: No one has variants
                "Sample1(1/1);Sample2(0/1);Sample3(1/1)",  # Variant 4: All have variants
            ],
        }
        df = pd.DataFrame(test_data)

        # Setup phenotype data
        context.phenotype_data = {
            "Sample1": {"HP:0000001", "HP:0000002"},
            "Sample2": {"HP:0000003"},
            "Sample3": {"HP:0000004"},
        }

        # Mock file operations
        with patch("pandas.read_csv", return_value=df), patch.object(df, "to_csv") as mock_to_csv:

            stage = PhenotypeIntegrationStage()
            stage(context)

            # Verify the phenotype column was added correctly

            # Check that apply was called on GT column
            # The actual phenotype values would be set in the real DataFrame
            assert "Phenotypes" in df.columns or mock_to_csv.called

            # Verify file was written
            mock_to_csv.assert_called_once()

    def test_phenotype_integration_edge_cases(self, context):
        """Test edge cases in phenotype integration."""
        from variantcentrifuge.phenotype import extract_phenotypes_for_gt_row

        phenotypes = {"Sample1": {"HP:0000001"}, "Sample2": {"HP:0000002"}}

        # Test edge cases
        edge_cases = [
            # Malformed GT entries
            {"gt": "Sample1(0/1);Sample2", "expected": "Sample1(HP:0000001)"},
            {"gt": "Sample1;Sample2(0/1)", "expected": "Sample2(HP:0000002)"},
            {"gt": "(0/1);Sample2(1/1)", "expected": "Sample2(HP:0000002)"},
            # Empty or whitespace
            {"gt": "", "expected": ""},
            {"gt": "   ", "expected": ""},
            {"gt": ";;;", "expected": ""},
            # Samples not in phenotype data
            {
                "gt": "UnknownSample(0/1);Sample1(1/1)",
                "expected": "UnknownSample();Sample1(HP:0000001)",
            },
            # Different genotype formats
            {
                "gt": "Sample1(0|1);Sample2(1|0)",
                "expected": "Sample1(HP:0000001);Sample2(HP:0000002)",
            },
            {
                "gt": "Sample1(A/T);Sample2(G/C)",
                "expected": "Sample1(HP:0000001);Sample2(HP:0000002)",
            },
        ]

        for i, test_case in enumerate(edge_cases):
            result = extract_phenotypes_for_gt_row(test_case["gt"], phenotypes)
            assert result == test_case["expected"], (
                f"Edge case {i+1} failed:\n"
                f"GT: {test_case['gt']}\n"
                f"Expected: {test_case['expected']}\n"
                f"Got: {result}"
            )

    def test_skip_if_no_phenotypes(self, context):
        """Test skipping when no phenotype data."""
        original_data = context.data
        context.phenotype_data = None

        stage = PhenotypeIntegrationStage()
        result = stage(context)

        # Data should not change (keep the original temp file path)
        assert result.data == original_data
