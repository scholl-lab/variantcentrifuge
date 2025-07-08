"""Unit tests for output stages."""

from pathlib import Path
from unittest.mock import Mock, patch, MagicMock, mock_open
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
        ctx.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1", "chr2"],
                "POS": [100, 200, 300],
                "REF": ["A", "G", "C"],
                "ALT": ["T", "C", "G"],
                "QUAL": [30, 40, 50],
            }
        )
        ctx.mark_complete("dataframe_loading")
        return ctx

    def test_add_variant_identifier(self, context):
        """Test adding variant identifiers."""
        stage = VariantIdentifierStage()
        result = stage(context)

        # Check that VAR_ID column was added
        assert "VAR_ID" in result.current_dataframe.columns

        # Check format of variant IDs
        variant_ids = result.current_dataframe["VAR_ID"].tolist()
        assert variant_ids[0].startswith("var_0001_")
        assert variant_ids[1].startswith("var_0002_")
        assert variant_ids[2].startswith("var_0003_")

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = VariantIdentifierStage()
        assert stage.dependencies == {"dataframe_loading"}

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = VariantIdentifierStage()
        assert stage.parallel_safe is False  # Must run in main process to modify DataFrame

    def test_streaming_mode(self, tmp_path):
        """Test streaming mode for large files."""
        import gzip

        # Create test input file
        input_file = tmp_path / "test_variants.tsv.gz"
        with gzip.open(input_file, "wt") as f:
            f.write("CHROM\tPOS\tREF\tALT\tQUAL\n")
            f.write("chr1\t100\tA\tT\t30\n")
            f.write("chr1\t200\tG\tC\t40\n")
            f.write("chr2\t300\tC\tG\t50\n")

        # Create context for streaming
        ctx = create_test_context()
        ctx.config["stream_variant_ids"] = True
        ctx.data = str(input_file)
        ctx.current_dataframe = None  # No DataFrame in streaming mode
        ctx.workspace = Mock()
        ctx.workspace.get_intermediate_path = lambda x: tmp_path / x
        ctx.mark_complete("dataframe_loading")  # Mark required dependency

        # Run stage
        stage = VariantIdentifierStage()
        result = stage(ctx)

        # Check output file was created
        output_file = tmp_path / "with_variant_ids.tsv.gz"
        assert output_file.exists()
        assert result.data == output_file

        # Verify content
        with gzip.open(output_file, "rt") as f:
            lines = f.readlines()
            assert lines[0].strip() == "VAR_ID\tCHROM\tPOS\tREF\tALT\tQUAL"
            # Check ID format: var_XXXX_yyyy where XXXX is sequential and yyyy is hash
            assert lines[1].startswith("var_0001_") and "\tchr1\t100\tA\tT\t30" in lines[1]
            assert lines[2].startswith("var_0002_") and "\tchr1\t200\tG\tC\t40" in lines[2]
            assert lines[3].startswith("var_0003_") and "\tchr2\t300\tC\tG\t50" in lines[3]

    def test_streaming_mode_missing_key_fields(self, tmp_path):
        """Test streaming mode when key fields are missing."""
        import gzip

        # Create test input file without standard fields
        input_file = tmp_path / "test_variants.tsv.gz"
        with gzip.open(input_file, "wt") as f:
            f.write("Gene\tImpact\tScore\n")
            f.write("BRCA1\tHIGH\t0.95\n")
            f.write("TP53\tMODERATE\t0.80\n")

        # Create context for streaming
        ctx = create_test_context()
        ctx.config["stream_variant_ids"] = True
        ctx.data = str(input_file)
        ctx.current_dataframe = None
        ctx.workspace = Mock()
        ctx.workspace.get_intermediate_path = lambda x: tmp_path / x
        ctx.mark_complete("dataframe_loading")

        # Run stage
        stage = VariantIdentifierStage()
        stage(ctx)

        # Verify fallback to index-based IDs
        output_file = tmp_path / "with_variant_ids.tsv.gz"
        with gzip.open(output_file, "rt") as f:
            lines = f.readlines()
            assert lines[0].strip() == "VAR_ID\tGene\tImpact\tScore"
            assert lines[1].startswith("var_0001_0000\t") and "\tBRCA1\tHIGH\t0.95" in lines[1]
            assert lines[2].startswith("var_0002_0000\t") and "\tTP53\tMODERATE\t0.80" in lines[2]


class TestFinalFilteringStage:
    """Test FinalFilteringStage."""

    @pytest.fixture
    def context(self):
        """Create test context with DataFrame."""
        ctx = create_test_context()
        ctx.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1", "chr2", "chr3"],
                "POS": [100, 200, 300, 400],
                "QUAL": [30, 40, 20, 50],
                "AF": [0.01, 0.1, 0.001, 0.5],
            }
        )
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
        # Use "filter" (singular) as expected by the implementation
        context.config["filter"] = "(QUAL >= 30) & (AF < 0.1)"

        with patch(
            "variantcentrifuge.stages.output_stages.filter_dataframe_with_query"
        ) as mock_apply:
            mock_apply.return_value = context.current_dataframe[
                (context.current_dataframe["QUAL"] >= 30) & (context.current_dataframe["AF"] < 0.1)
            ]

            stage = FinalFilteringStage()
            result = stage(context)

            # Verify filters were applied
            assert mock_apply.called
            # The return value from the mock has 1 row (filtered result)
            assert len(result.current_dataframe) == 1

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
            config_overrides={"pseudonymize": True, "pseudonymize_schema": "sequential"}
        )
        ctx.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1"],
                "POS": [100, 200],
                "Sample1": ["0/1", "1/1"],
                "Sample2": ["0/0", "0/1"],
                "Sample3": ["0/1", "1/1"],
            }
        )
        ctx.vcf_samples = ["Sample1", "Sample2", "Sample3"]
        ctx.mark_complete("dataframe_loading")
        ctx.mark_complete("sample_config_loading")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.create_pseudonymizer")
    def test_sequential_pseudonymization(self, mock_create_pseudonymizer, context):
        """Test sequential pseudonymization schema."""
        # Mock the pseudonymizer object
        mock_pseudonymizer = Mock()
        mock_mapping = {"Sample1": "SAMPLE_001", "Sample2": "SAMPLE_002", "Sample3": "SAMPLE_003"}
        mock_pseudonymizer.create_mapping.return_value = mock_mapping
        mock_pseudonymizer.pseudonymize_dataframe.return_value = context.current_dataframe.rename(
            columns={"Sample1": "SAMPLE_001", "Sample2": "SAMPLE_002", "Sample3": "SAMPLE_003"}
        )
        mock_create_pseudonymizer.return_value = mock_pseudonymizer

        stage = PseudonymizationStage()
        result = stage(context)

        # Check pseudonymizer was created
        mock_create_pseudonymizer.assert_called_once_with("sequential", prefix="SAMPLE")

        # Check mapping was created
        mock_pseudonymizer.create_mapping.assert_called_once_with(
            ["Sample1", "Sample2", "Sample3"], {}
        )

        # Check DataFrame columns were renamed
        assert "SAMPLE_001" in result.current_dataframe.columns
        assert "SAMPLE_002" in result.current_dataframe.columns
        assert "SAMPLE_003" in result.current_dataframe.columns
        assert "Sample1" not in result.current_dataframe.columns

        # Check context updated
        assert result.config["pseudonymization_mapping"] == mock_mapping

    @patch("variantcentrifuge.stages.output_stages.json.dump")
    @patch("variantcentrifuge.stages.output_stages.create_pseudonymizer")
    def test_mapping_file_creation(self, mock_create_pseudonymizer, mock_json_dump, context):
        """Test pseudonymization mapping file creation."""
        # Mock the pseudonymizer object
        mock_pseudonymizer = Mock()
        mock_mapping = {"Sample1": "SAMPLE_001"}
        mock_pseudonymizer.create_mapping.return_value = mock_mapping
        mock_pseudonymizer.pseudonymize_dataframe.return_value = context.current_dataframe
        mock_create_pseudonymizer.return_value = mock_pseudonymizer

        stage = PseudonymizationStage()
        stage(context)

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
        ctx.current_dataframe = pd.DataFrame(
            {"CHROM": ["chr1", "chr1"], "POS": [100, 200], "QUAL": [30, 40]}
        )
        ctx.mark_complete("dataframe_loading")
        ctx.mark_complete("variant_identifier")  # TSVOutputStage depends on this
        return ctx

    def test_tsv_output(self, context):
        """Test TSV file output."""
        # Mock the DataFrame's to_csv method
        with patch.object(context.current_dataframe, "to_csv") as mock_to_csv:
            stage = TSVOutputStage()
            result = stage(context)

            # Check to_csv was called correctly
            # Implementation uses additional parameters na_rep and compression
            assert mock_to_csv.call_count == 1
            call_args = mock_to_csv.call_args
            # First argument should be a Path object in the output directory
            assert isinstance(call_args[0][0], Path)
            assert str(call_args[0][0]).endswith("output.tsv")
            # Check other arguments
            assert call_args[1]["sep"] == "\t"
            assert call_args[1]["index"] is False
            assert call_args[1]["na_rep"] == ""
            assert call_args[1]["compression"] is None

            # Check context updated
            assert result.final_output_path == call_args[0][0]

    def test_stdout_output(self, context):
        """Test output to stdout."""
        context.config["output_file"] = "-"

        # Mock the DataFrame's to_csv method
        with patch.object(context.current_dataframe, "to_csv") as mock_to_csv:
            stage = TSVOutputStage()
            stage(context)

            # Check output went to stdout with proper arguments
            assert mock_to_csv.call_count == 1
            call_args = mock_to_csv.call_args
            # First argument should be sys.stdout (actual stdout, not a mock)
            import sys

            assert call_args[0][0] == sys.stdout
            # Check other arguments
            assert call_args[1]["sep"] == "\t"
            assert call_args[1]["index"] is False
            assert call_args[1]["na_rep"] == ""

    def test_missing_dataframe(self, context):
        """Test handling when no DataFrame available."""
        context.current_dataframe = None

        stage = TSVOutputStage()
        result = stage(context)

        # Should return context without error
        assert result is context
        # final_output_path should not be set when there's no data
        assert result.final_output_path is None

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = TSVOutputStage()
        assert stage.dependencies == {"dataframe_loading", "variant_identifier"}
        assert stage.soft_dependencies == {"variant_scoring", "final_filtering", "pseudonymization"}


class TestExcelReportStage:
    """Test ExcelReportStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context(config_overrides={"xlsx": True})
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.current_dataframe = pd.DataFrame({"CHROM": ["chr1", "chr1"], "POS": [100, 200]})
        ctx.variant_stats = {"total_variants": 2}
        ctx.gene_stats = pd.DataFrame({"gene": ["BRCA1"], "variants": [2]})
        ctx.mark_complete("dataframe_loading")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.convert_to_excel")
    def test_excel_generation(self, mock_convert, context):
        """Test Excel report generation."""
        # Mark TSV output as complete since Excel depends on it
        context.mark_complete("tsv_output")
        # Set final output path which TSV output would have created
        context.final_output_path = context.workspace.output_dir / "output.tsv"
        context.final_output_path.touch()  # Create the file

        stage = ExcelReportStage()
        stage(context)

        # Check convert_to_excel was called
        mock_convert.assert_called_once()
        # Implementation calls with positional args
        call_args = mock_convert.call_args[0]
        assert str(context.final_output_path) in call_args[0]
        assert str(call_args[1]).endswith(".xlsx")

    def test_skip_if_disabled(self, context):
        """Test skipping when Excel disabled."""
        context.config["xlsx"] = False
        # Mark TSV output as complete since Excel depends on it
        context.mark_complete("tsv_output")
        context.final_output_path = Path("/tmp/output.tsv")

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
            config_overrides={"html_report": True, "html_report_output": "/tmp/report.html"}
        )
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.current_dataframe = pd.DataFrame({"CHROM": ["chr1"], "POS": [100]})
        ctx.mark_complete("dataframe_loading")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.generate_html_report")
    @patch("variantcentrifuge.stages.output_stages.produce_report_json")
    def test_html_generation(self, mock_produce_json, mock_generate, context):
        """Test HTML report generation."""
        # Mark TSV output as complete since HTML depends on it
        context.mark_complete("tsv_output")
        # Set final output path which TSV output would have created
        context.final_output_path = context.workspace.output_dir / "output.tsv"
        context.final_output_path.touch()  # Create the file

        stage = HTMLReportStage()
        stage(context)

        # Check JSON was produced first
        mock_produce_json.assert_called_once()
        # Check report generation - implementation uses different parameters
        mock_generate.assert_called_once()
        call_kwargs = mock_generate.call_args[1]
        assert "json_file" in call_kwargs
        assert "output_file" in call_kwargs
        assert "title" in call_kwargs

    def test_skip_if_disabled(self, context):
        """Test skipping when HTML disabled."""
        context.config["html_report"] = False
        # Mark TSV output as complete since HTML depends on it
        context.mark_complete("tsv_output")
        context.final_output_path = Path("/tmp/output.tsv")

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
            vcf_file="/tmp/input.vcf", config_overrides={"igv_report": True, "reference": "hg38"}
        )
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.current_dataframe = pd.DataFrame({"CHROM": ["chr1"], "POS": [100]})
        ctx.mark_complete("dataframe_loading")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.generate_igv_report")
    @patch("variantcentrifuge.stages.output_stages.pd.read_csv")
    @patch("variantcentrifuge.stages.output_stages.match_IGV_link_columns")
    def test_igv_generation(self, mock_match_igv, mock_read_csv, mock_generate, context):
        """Test IGV report generation."""
        # Mark TSV output as complete since IGV depends on it
        context.mark_complete("tsv_output")
        # Set final output path which TSV output would have created
        context.final_output_path = context.workspace.output_dir / "output.tsv"
        context.final_output_path.touch()  # Create the file

        # Mock the CSV reading
        mock_df = pd.DataFrame({"CHROM": ["chr1"], "POS": [100]})
        mock_read_csv.return_value = mock_df

        # Mock IGV column matching
        mock_match_igv.return_value = {"chr": "CHROM", "start": "POS", "end": "POS"}

        stage = IGVReportStage()
        stage(context)

        # Check report generation
        mock_generate.assert_called_once()
        call_kwargs = mock_generate.call_args[1]
        assert call_kwargs["vcf_file"] == "/tmp/input.vcf"
        assert call_kwargs["tsv_file"] == str(context.final_output_path)
        assert call_kwargs["reference"] == "hg19"  # Default from implementation
        assert "chr" in call_kwargs
        assert "start" in call_kwargs
        assert "end" in call_kwargs

    def test_skip_if_disabled(self, context):
        """Test skipping when IGV disabled."""
        context.config["igv_report"] = False
        # Mark TSV output as complete since IGV depends on it
        context.mark_complete("tsv_output")
        context.final_output_path = Path("/tmp/output.tsv")

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
        ctx = create_test_context(gene_name="BRCA1", vcf_file="/tmp/input.vcf")
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.mark_complete("dataframe_loading")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.sanitize_metadata_field")
    @patch("variantcentrifuge.stages.output_stages.get_tool_version")
    def test_metadata_generation(self, mock_get_version, mock_sanitize, context):
        """Test metadata file generation."""
        # Mark TSV output as complete since metadata depends on it
        context.mark_complete("tsv_output")

        # Add normalized_genes to config
        context.config["normalized_genes"] = ["BRCA1"]

        # Mock tool versions
        mock_get_version.return_value = "1.0.0"

        # Mock sanitize to pass through the metadata unchanged
        mock_sanitize.side_effect = lambda x: x

        # Mock json.dump to capture the metadata
        metadata_written = None

        def capture_metadata(data, file, **kwargs):
            nonlocal metadata_written
            metadata_written = data

        with patch(
            "variantcentrifuge.stages.output_stages.json.dump", side_effect=capture_metadata
        ):
            with patch("builtins.open", mock_open()):
                stage = MetadataGenerationStage()
                stage(context)

        # Check metadata contents
        assert metadata_written is not None
        assert metadata_written["input_files"]["genes"] == ["BRCA1"]
        assert metadata_written["input_files"]["vcf"] == "/tmp/input.vcf"
        assert "run_date" in metadata_written
        assert metadata_written["pipeline_version"] == "1.0.0-test"  # From context fixture

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
            output_dir="/tmp/results", config_overrides={"archive_results": True}
        )
        ctx.mark_complete("tsv_output")
        return ctx

    @patch("tarfile.open")
    def test_archive_creation(self, mock_tarfile_open, context):
        """Test archive creation."""
        # Mock the tarfile context manager
        mock_tar = MagicMock()
        mock_tarfile_open.return_value.__enter__.return_value = mock_tar

        # Create the archive path to avoid stat() error
        archive_path = context.workspace.get_archive_path()

        # Mock Path.stat() to avoid FileNotFoundError
        with patch.object(Path, "stat") as mock_stat:
            mock_stat.return_value.st_size = 1024 * 1024  # 1MB

            stage = ArchiveCreationStage()
            result = stage(context)

        # Check tarfile was opened correctly
        assert mock_tarfile_open.called
        call_args = mock_tarfile_open.call_args
        # First argument should be the archive path
        assert call_args[0][0] == archive_path
        assert call_args[0][1] == "w:gz"

        # Check that output_dir was added to archive
        mock_tar.add.assert_called_once()
        add_args = mock_tar.add.call_args[0]
        assert add_args[0] == context.workspace.output_dir

        # Check result has archive path
        assert result.report_paths["archive"] == archive_path

    def test_skip_if_disabled(self, context):
        """Test skipping when archiving disabled."""
        context.config["archive_results"] = False

        with patch("tarfile.open") as mock:
            stage = ArchiveCreationStage()
            result = stage(context)

            mock.assert_not_called()
            # archive_path is not a direct attribute but should be in report_paths
            assert "archive" not in result.report_paths

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = ArchiveCreationStage()
        assert stage.parallel_safe is True
