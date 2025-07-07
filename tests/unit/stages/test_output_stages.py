"""Unit tests for output stages."""

from pathlib import Path
from unittest.mock import Mock, patch, MagicMock, ANY
import json
import pandas as pd
import pytest

from variantcentrifuge.stages.output_stages import (
    VariantIdentifierStage,
    FinalFilteringStage,
    PseudonymizationStage,
    TSVOutputStage,
    ExcelReportStage,
    HTMLReportStage,
    IGVReportStage,
    MetadataGenerationStage,
    ArchiveCreationStage,
)
from tests.mocks.fixtures import create_test_context


class TestVariantIdentifierStage:
    """Test VariantIdentifierStage."""

    @pytest.fixture
    def context(self):
        """Create test context with DataFrame."""
        ctx = create_test_context()
        # Create a test DataFrame
        ctx.current_dataframe = pd.DataFrame({
            "CHROM": ["chr1", "chr1", "chr2"],
            "POS": [100, 200, 300],
            "REF": ["A", "G", "C"],
            "ALT": ["T", "C", "G"],
            "QUAL": [30, 40, 50]
        })
        ctx.mark_complete("dataframe_loading")
        return ctx

    def test_add_variant_identifier(self, context):
        """Test adding variant identifiers."""
        stage = VariantIdentifierStage()
        result = stage(context)

        # Check that Variant_ID column was added
        assert "Variant_ID" in result.current_dataframe.columns
        
        # Check format of variant IDs
        variant_ids = result.current_dataframe["Variant_ID"].tolist()
        assert variant_ids[0] == "chr1:100:A>T"
        assert variant_ids[1] == "chr1:200:G>C"
        assert variant_ids[2] == "chr2:300:C>G"

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = VariantIdentifierStage()
        assert stage.dependencies == {"dataframe_loading"}

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = VariantIdentifierStage()
        assert stage.parallel_safe is True


class TestFinalFilteringStage:
    """Test FinalFilteringStage."""

    @pytest.fixture
    def context(self):
        """Create test context with DataFrame."""
        ctx = create_test_context()
        ctx.current_dataframe = pd.DataFrame({
            "CHROM": ["chr1", "chr1", "chr2", "chr3"],
            "POS": [100, 200, 300, 400],
            "QUAL": [30, 40, 20, 50],
            "AF": [0.01, 0.1, 0.001, 0.5]
        })
        ctx.mark_complete("dataframe_loading")
        return ctx

    def test_final_filtering(self, context):
        """Test final filtering with pandas query."""
        context.config["final_filter"] = "QUAL >= 30 & AF < 0.1"
        
        stage = FinalFilteringStage()
        result = stage(context)

        # Check filtering worked (only chr1:100 should remain)
        assert len(result.current_dataframe) == 1
        assert all(result.current_dataframe["QUAL"] >= 30)
        assert all(result.current_dataframe["AF"] < 0.1)

    def test_late_filtering(self, context):
        """Test late filtering mode."""
        context.config["late_filtering"] = True
        context.config["filters"] = ["QUAL >= 30", "AF < 0.1"]
        
        with patch("variantcentrifuge.filters.apply_filters") as mock_apply:
            mock_apply.return_value = context.current_dataframe[
                (context.current_dataframe["QUAL"] >= 30) & 
                (context.current_dataframe["AF"] < 0.1)
            ]
            
            stage = FinalFilteringStage()
            result = stage(context)
            
            # Verify filters were applied
            assert mock_apply.called
            assert len(result.current_dataframe) == 2

    def test_no_filtering(self, context):
        """Test when no filtering is needed."""
        # No filters configured
        stage = FinalFilteringStage()
        result = stage(context)

        # DataFrame should be unchanged
        assert len(result.current_dataframe) == 4

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = FinalFilteringStage()
        assert stage.parallel_safe is True


class TestPseudonymizationStage:
    """Test PseudonymizationStage."""

    @pytest.fixture
    def context(self):
        """Create test context with sample data."""
        ctx = create_test_context(
            config_overrides={
                "pseudonymize": True,
                "pseudonymize_schema": "sequential"
            }
        )
        ctx.current_dataframe = pd.DataFrame({
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "Sample1": ["0/1", "1/1"],
            "Sample2": ["0/0", "0/1"],
            "Sample3": ["0/1", "1/1"]
        })
        ctx.vcf_samples = ["Sample1", "Sample2", "Sample3"]
        ctx.mark_complete("dataframe_loading")
        ctx.mark_complete("sample_config_loading")
        return ctx

    @patch("variantcentrifuge.pseudonymize.create_pseudonymization_mapping")
    def test_sequential_pseudonymization(self, mock_create_mapping, context):
        """Test sequential pseudonymization schema."""
        # Mock the mapping function
        mock_mapping = {
            "Sample1": "SAMPLE_001",
            "Sample2": "SAMPLE_002", 
            "Sample3": "SAMPLE_003"
        }
        mock_create_mapping.return_value = mock_mapping
        
        stage = PseudonymizationStage()
        result = stage(context)

        # Check mapping was created
        mock_create_mapping.assert_called_once_with(
            ["Sample1", "Sample2", "Sample3"],
            "sequential",
            None
        )

        # Check DataFrame columns were renamed
        assert "SAMPLE_001" in result.current_dataframe.columns
        assert "SAMPLE_002" in result.current_dataframe.columns
        assert "SAMPLE_003" in result.current_dataframe.columns
        assert "Sample1" not in result.current_dataframe.columns

        # Check context updated
        assert result.pseudonymization_mapping == mock_mapping

    @patch("variantcentrifuge.stages.output_stages.json.dump")
    @patch("variantcentrifuge.pseudonymize.create_pseudonymization_mapping")
    def test_mapping_file_creation(self, mock_create_mapping, mock_json_dump, context):
        """Test pseudonymization mapping file creation."""
        mock_mapping = {"Sample1": "SAMPLE_001"}
        mock_create_mapping.return_value = mock_mapping
        
        stage = PseudonymizationStage()
        result = stage(context)

        # Check JSON was written
        assert mock_json_dump.called
        call_args = mock_json_dump.call_args[0]
        assert call_args[0] == mock_mapping

    def test_skip_if_disabled(self, context):
        """Test skipping when pseudonymization disabled."""
        context.config["pseudonymize"] = False
        
        stage = PseudonymizationStage()
        result = stage(context)

        # Columns should be unchanged
        assert "Sample1" in result.current_dataframe.columns
        assert "SAMPLE_001" not in result.current_dataframe.columns


class TestTSVOutputStage:
    """Test TSVOutputStage."""

    @pytest.fixture
    def context(self):
        """Create test context with DataFrame."""
        ctx = create_test_context()
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.current_dataframe = pd.DataFrame({
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "QUAL": [30, 40]
        })
        ctx.mark_complete("dataframe_loading")
        return ctx

    def test_tsv_output(self, context):
        """Test TSV file output."""
        # Mock the DataFrame's to_csv method
        with patch.object(context.current_dataframe, 'to_csv') as mock_to_csv:
            stage = TSVOutputStage()
            result = stage(context)

            # Check to_csv was called correctly
            mock_to_csv.assert_called_once_with(
                "/tmp/output.tsv",
                sep="\t",
                index=False
            )

            # Check context updated
            assert result.final_output_file == Path("/tmp/output.tsv")

    @patch("variantcentrifuge.stages.output_stages.sys.stdout")
    def test_stdout_output(self, mock_stdout, context):
        """Test output to stdout."""
        context.config["output_file"] = "-"
        
        # Mock the DataFrame's to_csv method
        with patch.object(context.current_dataframe, 'to_csv') as mock_to_csv:
            stage = TSVOutputStage()
            result = stage(context)

            # Check output went to stdout
            mock_to_csv.assert_called_once_with(
                mock_stdout,
                sep="\t",
                index=False
            )

    def test_missing_dataframe_error(self, context):
        """Test error when no DataFrame available."""
        context.current_dataframe = None
        
        stage = TSVOutputStage()
        with pytest.raises(ValueError, match="No data available"):
            stage(context)


class TestExcelReportStage:
    """Test ExcelReportStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context(
            config_overrides={"xlsx": True}
        )
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.current_dataframe = pd.DataFrame({
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 200]
        })
        ctx.variant_stats = {"total_variants": 2}
        ctx.gene_stats = pd.DataFrame({"gene": ["BRCA1"], "variants": [2]})
        ctx.mark_complete("dataframe_loading")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.convert_to_excel")
    def test_excel_generation(self, mock_convert, context):
        """Test Excel report generation."""
        stage = ExcelReportStage()
        result = stage(context)

        # Check convert_to_excel was called
        mock_convert.assert_called_once()
        call_kwargs = mock_convert.call_args[1]
        assert call_kwargs["output_file"] == "/tmp/output.tsv"
        assert isinstance(call_kwargs["variant_df"], pd.DataFrame)

    def test_skip_if_disabled(self, context):
        """Test skipping when Excel disabled."""
        context.config["xlsx"] = False
        
        with patch("variantcentrifuge.stages.output_stages.convert_to_excel") as mock:
            stage = ExcelReportStage()
            stage(context)
            
            # Should not generate Excel
            mock.assert_not_called()

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = ExcelReportStage()
        assert stage.parallel_safe is True


class TestHTMLReportStage:
    """Test HTMLReportStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context(
            config_overrides={
                "html_report": True,
                "html_report_output": "/tmp/report.html"
            }
        )
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.current_dataframe = pd.DataFrame({
            "CHROM": ["chr1"], 
            "POS": [100]
        })
        ctx.mark_complete("dataframe_loading")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.generate_html_report")
    def test_html_generation(self, mock_generate, context):
        """Test HTML report generation."""
        stage = HTMLReportStage()
        result = stage(context)

        # Check report generation
        mock_generate.assert_called_once_with(
            dataframe=ANY,
            output_filename="/tmp/report.html",
            title=ANY,
            metadata=ANY
        )

    def test_skip_if_disabled(self, context):
        """Test skipping when HTML disabled."""
        context.config["html_report"] = False
        
        with patch("variantcentrifuge.stages.output_stages.generate_html_report") as mock:
            stage = HTMLReportStage()
            stage(context)
            
            mock.assert_not_called()

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = HTMLReportStage()
        assert stage.parallel_safe is True


class TestIGVReportStage:
    """Test IGVReportStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context(
            vcf_file="/tmp/input.vcf",
            config_overrides={
                "igv_report": True,
                "reference": "hg38"
            }
        )
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.current_dataframe = pd.DataFrame({
            "CHROM": ["chr1"],
            "POS": [100]
        })
        ctx.mark_complete("dataframe_loading")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.generate_igv_report")
    def test_igv_generation(self, mock_generate, context):
        """Test IGV report generation."""
        stage = IGVReportStage()
        result = stage(context)

        # Check report generation
        mock_generate.assert_called_once()
        call_kwargs = mock_generate.call_args[1]
        assert call_kwargs["vcf_file"] == "/tmp/input.vcf"
        assert call_kwargs["reference"] == "hg38"

    def test_skip_if_disabled(self, context):
        """Test skipping when IGV disabled."""
        context.config["igv_report"] = False
        
        with patch("variantcentrifuge.stages.output_stages.generate_igv_report") as mock:
            stage = IGVReportStage()
            stage(context)
            
            mock.assert_not_called()

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = IGVReportStage()
        assert stage.parallel_safe is True


class TestMetadataGenerationStage:
    """Test MetadataGenerationStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context(
            gene_name="BRCA1",
            vcf_file="/tmp/input.vcf"
        )
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.mark_complete("dataframe_loading")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.generate_metadata")
    def test_metadata_generation(self, mock_generate, context):
        """Test metadata file generation."""
        stage = MetadataGenerationStage()
        result = stage(context)

        # Check metadata generation
        mock_generate.assert_called_once()
        call_args = mock_generate.call_args[0]
        metadata = call_args[0]
        
        # Check metadata contents
        assert metadata["gene_name"] == "BRCA1"
        assert metadata["vcf_file"] == "/tmp/input.vcf"
        assert "timestamp" in metadata

    def test_stage_properties(self):
        """Test stage properties."""
        stage = MetadataGenerationStage()
        assert stage.dependencies == {"tsv_output"}
        assert stage.parallel_safe is True


class TestArchiveCreationStage:
    """Test ArchiveCreationStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context(
            output_dir="/tmp/results",
            config_overrides={"archive_results": True}
        )
        ctx.mark_complete("tsv_output")
        return ctx

    @patch("variantcentrifuge.helpers.create_archive")
    def test_archive_creation(self, mock_create, context):
        """Test archive creation."""
        mock_create.return_value = "/tmp/results_archive.tar.gz"
        
        stage = ArchiveCreationStage()
        result = stage(context)

        # Check archive creation
        mock_create.assert_called_once_with("/tmp/results")
        assert result.archive_path == Path("/tmp/results_archive.tar.gz")

    def test_skip_if_disabled(self, context):
        """Test skipping when archiving disabled."""
        context.config["archive_results"] = False
        
        with patch("variantcentrifuge.stages.output_stages.create_archive") as mock:
            stage = ArchiveCreationStage()
            result = stage(context)
            
            mock.assert_not_called()
            assert result.archive_path is None

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = ArchiveCreationStage()
        assert stage.parallel_safe is True