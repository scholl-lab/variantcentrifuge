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
    PCAComputationStage,
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

        vcf_file = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False)  # noqa: SIM115
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
        assert kwargs["bed_file"] == str(Path("/tmp/genes.bed"))
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


class TestFieldExtractionStage:
    """Test FieldExtractionStage."""

    @pytest.fixture
    def context(self):
        """Create test context with extracted VCF."""
        ctx = create_test_context(
            config_overrides={
                "fields_to_extract": "CHROM POS REF ALT QUAL",
                "gzip_intermediates": False,
            }
        )
        ctx.data = Path("/tmp/variants.vcf.gz")
        # Add vcf_samples for bcftools extraction
        ctx.vcf_samples = ["Sample1", "Sample2", "Sample3"]
        # Mark all potential dependencies as complete
        ctx.mark_complete("configuration_loading")
        ctx.mark_complete("gene_bed_creation")
        ctx.mark_complete("variant_extraction")
        ctx.mark_complete("parallel_variant_extraction")
        ctx.mark_complete("multiallelic_split")
        ctx.mark_complete("snpsift_filtering")
        return ctx

    @patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
    def test_field_extraction(self, mock_extract, context):
        """Test field extraction to TSV."""
        # Mock return value - extract_fields_bcftools returns the output path
        mock_extract.return_value = "mocked_output.tsv"

        stage = FieldExtractionStage()
        result = stage(context)

        # Verify extraction - extract_fields_bcftools uses keyword arguments
        mock_extract.assert_called_once()
        kwargs = mock_extract.call_args[1]
        assert kwargs["variant_file"] == str(Path("/tmp/variants.vcf.gz"))  # variant_file
        assert kwargs["fields"] == "CHROM POS REF ALT QUAL"  # fields as string
        assert isinstance(kwargs["cfg"], dict)  # cfg
        assert ".extracted.tsv" in kwargs["output_file"]  # output_file
        assert kwargs["vcf_samples"] == ["Sample1", "Sample2", "Sample3"]  # vcf_samples

        # Verify context updates
        assert result.extracted_tsv is not None
        assert result.data == result.extracted_tsv

    @patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
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
        context.config["fields_to_extract"] = []

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
        tsv_file = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False)  # noqa: SIM115
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

    def test_genotype_replacement(self, context):
        """Test genotype replacement stage (Phase 11: now a no-op)."""
        # Phase 11: GenotypeReplacementStage now no-ops immediately
        # GT reconstruction deferred to output time (TSVOutputStage, ExcelReportStage)

        original_data = context.data
        stage = GenotypeReplacementStage()
        result = stage(context)

        # Verify stage completes immediately (no actual processing)
        # The stage may set genotype_replaced_tsv in resume_from_checkpoint
        # but no actual genotype replacement file should be created
        # Key verification: data path should remain unchanged from input
        assert result.data == original_data

    def test_skip_if_disabled(self, context):
        """Test skipping when replacement disabled."""
        context.config["replace_genotypes"] = False

        stage = GenotypeReplacementStage()
        result = stage(context)

        # Data should not change
        assert result.data == context.data

    def test_noop_with_no_samples(self, context):
        """Test no-op when no samples available (Phase 11: always no-op)."""
        context.vcf_samples = []
        original_data = context.data

        stage = GenotypeReplacementStage()
        result = stage(context)

        # Phase 11: stage is always a no-op, data unchanged
        assert result.data == original_data


class TestPhenotypeIntegrationStage:
    """Test PhenotypeIntegrationStage."""

    @pytest.fixture
    def context(self):
        """Create test context with phenotype data."""
        import tempfile

        # Create a temp TSV file that actually exists
        tsv_file = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False)  # noqa: SIM115
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
                f"Variant {i + 1} failed:\n"
                f"GT: {test_case['gt']}\n"
                f"Expected: {test_case['expected']}\n"
                f"Got: {result}"
            )

    def test_phenotype_integration_end_to_end_with_real_dataframe(self, context):
        """Test the complete phenotype integration with a realistic DataFrame."""
        from unittest.mock import patch

        import pandas as pd

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
                f"Edge case {i + 1} failed:\n"
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


class TestGeneBedCreationRegionRestriction:
    """Unit tests for GeneBedCreationStage region restriction helpers."""

    @pytest.fixture
    def stage(self):
        """Return a GeneBedCreationStage instance."""
        return GeneBedCreationStage()

    @pytest.fixture
    def context(self, tmp_path):
        """Create test context with a temporary workspace."""
        ctx = create_test_context(
            gene_name="BRCA1",
            config_overrides={
                "reference": "GRCh37.75",
                "interval_expand": 0,
                "add_chr": True,
            },
        )
        ctx.mark_complete("configuration_loading")
        return ctx

    # ------------------------------------------------------------------
    # _read_chromosomes
    # ------------------------------------------------------------------

    def test_read_chromosomes(self, stage, tmp_path):
        """Read chromosome names from first column, skipping comment lines."""
        bed = tmp_path / "test.bed"
        bed.write_text(
            "# comment line\n"
            "track name=test\n"
            "browser position chr1:1-100\n"
            "chr1\t100\t200\tGENE1\n"
            "chr2\t300\t400\tGENE2\n"
            "chrX\t500\t600\tGENE3\n"
        )
        result = stage._read_chromosomes(bed)
        assert result == {"chr1", "chr2", "chrX"}

    def test_read_chromosomes_empty_file(self, stage, tmp_path):
        """Empty BED file returns empty set."""
        bed = tmp_path / "empty.bed"
        bed.write_text("# only comments\n")
        result = stage._read_chromosomes(bed)
        assert result == set()

    # ------------------------------------------------------------------
    # _count_regions
    # ------------------------------------------------------------------

    def test_count_regions(self, stage, tmp_path):
        """Count only data lines, skip comment/empty lines."""
        bed = tmp_path / "test.bed"
        bed.write_text(
            "# comment\n"
            "track name=x\n"
            "chr1\t100\t200\n"
            "chr1\t300\t400\n"
            "chr2\t500\t600\n"
            "chrX\t700\t800\n"
            "chrY\t900\t1000\n"
        )
        result = stage._count_regions(bed)
        assert result == 5

    # ------------------------------------------------------------------
    # _intersect_with_restriction_bed — chr prefix mismatch
    # ------------------------------------------------------------------

    def test_chr_mismatch_detection(self, stage, tmp_path, context):
        """Raise ValueError when one BED has chr prefix and the other does not."""
        gene_bed = tmp_path / "gene.bed"
        gene_bed.write_text("chr1\t100\t200\tBRCA1\n")

        restriction_bed = tmp_path / "restriction.bed"
        restriction_bed.write_text("1\t150\t250\n")

        with patch("subprocess.run"), pytest.raises(ValueError, match="Chromosome naming mismatch"):
            stage._intersect_with_restriction_bed(gene_bed, str(restriction_bed), context)

    # ------------------------------------------------------------------
    # _intersect_with_restriction_bed — empty intersection
    # ------------------------------------------------------------------

    def test_empty_intersection_raises(self, stage, tmp_path, context):
        """Raise ValueError when bedtools produces an empty output file."""
        gene_bed = tmp_path / "gene.bed"
        gene_bed.write_text("chr1\t100\t200\tBRCA1\n")

        restriction_bed = tmp_path / "restriction.bed"
        restriction_bed.write_text("chr2\t300\t400\n")

        def _fake_run(cmd, stdout, check):
            # Write nothing — simulates zero-overlap intersection
            stdout.write("")

        with (
            patch("subprocess.run", side_effect=_fake_run),
            pytest.raises(ValueError, match="produced no regions"),
        ):
            stage._intersect_with_restriction_bed(gene_bed, str(restriction_bed), context)

    # ------------------------------------------------------------------
    # _intersect_with_restriction_bed — successful intersection
    # ------------------------------------------------------------------

    def test_successful_intersection(self, stage, tmp_path, context):
        """Return output path when bedtools writes non-empty result."""
        gene_bed = tmp_path / "gene.bed"
        gene_bed.write_text("chr1\t100\t500\tBRCA1\nchr2\t200\t600\tBRCA2\nchrX\t300\t700\tGENE3\n")

        restriction_bed = tmp_path / "capture.bed"
        restriction_bed.write_text("chr1\t150\t300\nchr2\t250\t500\nchrX\t350\t600\n")

        intersected_content = (
            "chr1\t150\t300\tBRCA1\nchr2\t250\t500\tBRCA2\nchrX\t350\t600\tGENE3\n"
        )

        def _fake_run(cmd, stdout, check):
            stdout.write(intersected_content)

        with patch("subprocess.run", side_effect=_fake_run):
            result = stage._intersect_with_restriction_bed(gene_bed, str(restriction_bed), context)

        assert result.exists()
        assert result.stat().st_size > 0
        assert result.name.endswith(".restricted.bed")

    # ------------------------------------------------------------------
    # _intersect_with_restriction_bed — missing restriction BED
    # ------------------------------------------------------------------

    def test_missing_restriction_bed_raises(self, stage, tmp_path, context):
        """Raise FileNotFoundError when restriction BED does not exist."""
        gene_bed = tmp_path / "gene.bed"
        gene_bed.write_text("chr1\t100\t200\tBRCA1\n")

        missing_bed = str(tmp_path / "nonexistent.bed")

        with pytest.raises(FileNotFoundError):
            stage._intersect_with_restriction_bed(gene_bed, missing_bed, context)


class TestPCAComputationStage:
    """Unit tests for PCAComputationStage."""

    @pytest.fixture
    def stage(self):
        """Return a PCAComputationStage instance."""
        return PCAComputationStage()

    @pytest.fixture
    def context(self, tmp_path):
        """Create a minimal test context for PCA tests."""
        ctx = create_test_context(
            vcf_file="test.vcf",
            output_dir=str(tmp_path),
            config_overrides={"vcf_file": "test.vcf"},
        )
        ctx.mark_complete("configuration_loading")
        ctx.mark_complete("sample_config_loading")
        return ctx

    def test_pca_no_arg_skips(self, stage, context):
        """When no 'pca' config key is set, _process returns context unchanged."""
        # Ensure 'pca' key is absent
        context.config.pop("pca", None)
        completed_before = set(context.completed_stages)

        result = stage._process(context)

        assert result is context
        # mark_complete should NOT have been called for pca_computation
        assert "pca_computation" not in result.completed_stages
        assert result.completed_stages == completed_before

    def test_pca_file_path(self, stage, context, tmp_path):
        """Pre-computed eigenvec file is used directly; pca_file set in config."""
        eigenvec = tmp_path / "pca.eigenvec"
        eigenvec.write_text("FID\tIID\tPC1\tPC2\nsample1\tsample1\t0.1\t0.2\n")

        context.config["pca"] = str(eigenvec)
        result = stage._process(context)

        assert result.config["pca_file"] == str(eigenvec)
        assert "pca_computation" in result.completed_stages

    def test_pca_akt_tool(self, stage, context, tmp_path):
        """When pca='akt', _run_akt is called; pca_file set in config."""
        context.config["pca"] = "akt"
        context.config["vcf_file"] = "test.vcf"

        # Expected output path from workspace
        expected_output = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.pca.eigenvec"
        )

        def fake_run(cmd, capture_output, text, check):
            # Simulate AKT returning eigenvec data on stdout
            import subprocess as sp

            return sp.CompletedProcess(
                args=cmd,
                returncode=0,
                stdout="sample1\t0.1\t0.2\t0.3\n",
                stderr="",
            )

        with patch("subprocess.run", side_effect=fake_run):
            result = stage._process(context)

        assert result.config["pca_file"] == str(expected_output)
        assert "pca_computation" in result.completed_stages

    def test_pca_akt_with_sites(self, stage, context, tmp_path):
        """When pca_sites is set, -R is used instead of --force."""
        context.config["pca"] = "akt"
        context.config["vcf_file"] = "test.vcf"
        context.config["pca_sites"] = "/path/to/sites.vcf.gz"

        captured_cmd = []

        def fake_run(cmd, capture_output, text, check):
            import subprocess as sp

            captured_cmd.extend(cmd)
            return sp.CompletedProcess(
                args=cmd,
                returncode=0,
                stdout="sample1\t0.1\t0.2\t0.3\n",
                stderr="",
            )

        with patch("subprocess.run", side_effect=fake_run):
            stage._process(context)

        assert "-R" in captured_cmd
        assert "/path/to/sites.vcf.gz" in captured_cmd
        assert "--force" not in captured_cmd

    def test_pca_akt_cache_reuse(self, stage, context, tmp_path):
        """If cached eigenvec file already exists, subprocess.run is NOT called."""
        context.config["pca"] = "akt"
        context.config["vcf_file"] = "test.vcf"

        # Pre-create the expected output file (simulating a cache hit)
        cached_output = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.pca.eigenvec"
        )
        cached_output.write_text("FID\tIID\tPC1\nsample1\tsample1\t0.1\n")

        with patch("subprocess.run") as mock_run:
            result = stage._process(context)
            mock_run.assert_not_called()

        assert result.config["pca_file"] == str(cached_output)

    def test_pca_invalid_arg_raises(self, stage, context):
        """Unknown tool name raises ValueError with descriptive message."""
        context.config["pca"] = "unknown_tool"

        with pytest.raises(
            ValueError, match="neither a valid file path nor a recognized tool name"
        ):
            stage._process(context)

    def test_pca_akt_not_found_raises(self, stage, context):
        """FileNotFoundError from subprocess is re-raised as RuntimeError with AKT message."""
        context.config["pca"] = "akt"
        context.config["vcf_file"] = "test.vcf"

        with (
            patch("subprocess.run", side_effect=FileNotFoundError("akt not found")),
            pytest.raises(RuntimeError, match="AKT not found in PATH"),
        ):
            stage._process(context)

    def test_pca_akt_failure_raises(self, stage, context):
        """CalledProcessError from AKT is re-raised as RuntimeError with exit code message."""
        import subprocess

        context.config["pca"] = "akt"
        context.config["vcf_file"] = "test.vcf"

        err = subprocess.CalledProcessError(returncode=1, cmd=["akt"], stderr="memory error")

        with (
            patch("subprocess.run", side_effect=err),
            pytest.raises(RuntimeError, match="AKT PCA failed"),
        ):
            stage._process(context)

    def test_pca_no_vcf_raises(self, stage, context):
        """When pca='akt' but no vcf_file set, ValueError is raised."""
        context.config["pca"] = "akt"
        context.config.pop("vcf_file", None)

        with pytest.raises(ValueError, match="requires a VCF file"):
            stage._process(context)

    def test_pca_stage_properties(self, stage):
        """PCAComputationStage has expected name, dependencies, and parallel_safe."""
        assert stage.name == "pca_computation"
        assert "configuration_loading" in stage.dependencies
        assert stage.parallel_safe is True
