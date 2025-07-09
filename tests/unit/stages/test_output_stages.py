"""Unit tests for output stages."""

from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pandas as pd
import pytest

from tests.mocks.fixtures import create_test_context
from variantcentrifuge.stages.output_stages import (
    ArchiveCreationStage,
    ExcelReportStage,
    FinalFilteringStage,
    HTMLReportStage,
    IGVReportStage,
    MetadataGenerationStage,
    PseudonymizationStage,
    TSVOutputStage,
    VariantIdentifierStage,
)


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

    def test_soft_dependencies(self):
        """Test soft dependencies are correctly defined."""
        stage = FinalFilteringStage()
        # Should have soft dependencies on stages that compute columns
        assert "variant_scoring" in stage.soft_dependencies
        assert "inheritance_analysis" in stage.soft_dependencies
        assert "custom_annotation" in stage.soft_dependencies

    def test_empty_filter_expression(self, context):
        """Test that empty filter expression returns unmodified DataFrame.

        Migrated from test_filters.py.
        """
        context.config["final_filter"] = ""

        stage = FinalFilteringStage()
        result = stage(context)

        # DataFrame should be unchanged
        assert len(result.current_dataframe) == 4
        pd.testing.assert_frame_equal(result.current_dataframe, context.current_dataframe)

    def test_simple_numeric_filter(self, context):
        """Test simple numeric comparison filter.

        Migrated from test_filters.py.
        """
        # Add string numeric values that should be converted
        context.current_dataframe["score"] = ["0.5", "0.8", "0.3", "0.9"]
        context.config["final_filter"] = "score > 0.4"

        stage = FinalFilteringStage()
        result = stage(context)

        # Check filtering worked
        assert len(result.current_dataframe) == 3
        assert result.current_dataframe["CHROM"].tolist() == ["chr1", "chr1", "chr3"]

    def test_string_equality_filter(self, context):
        """Test string equality filter.

        Migrated from test_filters.py.
        """
        context.current_dataframe["IMPACT"] = ["HIGH", "MODERATE", "HIGH", "LOW"]
        context.config["final_filter"] = 'IMPACT == "HIGH"'

        stage = FinalFilteringStage()
        result = stage(context)

        assert len(result.current_dataframe) == 2
        assert result.current_dataframe["CHROM"].tolist() == ["chr1", "chr2"]

    def test_complex_filter_with_and(self, context):
        """Test complex filter with AND condition.

        Migrated from test_filters.py.
        """
        context.current_dataframe["IMPACT"] = ["HIGH", "HIGH", "MODERATE", "HIGH"]
        context.current_dataframe["score"] = ["0.5", "0.8", "0.3", "0.9"]
        context.config["final_filter"] = 'IMPACT == "HIGH" and score > 0.6'

        stage = FinalFilteringStage()
        result = stage(context)

        assert len(result.current_dataframe) == 2
        assert result.current_dataframe["CHROM"].tolist() == ["chr1", "chr3"]

    def test_complex_filter_with_or(self, context):
        """Test complex filter with OR condition.

        Migrated from test_filters.py.
        """
        context.current_dataframe["IMPACT"] = ["HIGH", "MODERATE", "LOW", "LOW"]
        context.current_dataframe["score"] = ["0.9", "0.4", "0.2", "0.85"]
        context.config["final_filter"] = 'IMPACT == "HIGH" or score > 0.8'

        stage = FinalFilteringStage()
        result = stage(context)

        assert len(result.current_dataframe) == 2
        assert result.current_dataframe["CHROM"].tolist() == ["chr1", "chr3"]

    def test_filter_with_in_operator(self, context):
        """Test filter using 'in' operator.

        Migrated from test_filters.py.
        """
        context.current_dataframe["Inheritance_Pattern"] = [
            "de_novo",
            "autosomal_recessive",
            "compound_heterozygous",
            "unknown",
        ]
        context.config["final_filter"] = (
            'Inheritance_Pattern in ["de_novo", "compound_heterozygous"]'
        )

        stage = FinalFilteringStage()
        result = stage(context)

        assert len(result.current_dataframe) == 2
        assert set(result.current_dataframe["Inheritance_Pattern"]) == {
            "de_novo",
            "compound_heterozygous",
        }

    def test_filter_with_string_contains(self, context):
        """Test filter using string contains method.

        Migrated from test_filters.py.
        """
        context.current_dataframe["Custom_Annotation"] = [
            "InGeneList=cancer_panel",
            "Region=exon",
            "InGeneList=cardiac_panel",
            "Region=intron",
        ]
        context.config["final_filter"] = 'Custom_Annotation.str.contains("cancer")'

        stage = FinalFilteringStage()
        result = stage(context)

        assert len(result.current_dataframe) == 1
        assert result.current_dataframe["CHROM"].tolist() == ["chr1"]

    def test_numeric_conversion_with_mixed_types(self, context):
        """Test that numeric conversion handles mixed types correctly.

        Migrated from test_filters.py.
        """
        context.current_dataframe["dbNSFP_CADD_phred"] = ["25.5", "NA", "30.2", "22.1"]
        context.current_dataframe["gene"] = ["BRCA1", "TP53", "EGFR", "KRAS"]
        context.config["final_filter"] = "dbNSFP_CADD_phred > 26"

        stage = FinalFilteringStage()
        result = stage(context)

        assert len(result.current_dataframe) == 1
        assert result.current_dataframe["gene"].tolist() == ["EGFR"]

    def test_invalid_filter_expression(self, context, caplog):
        """Test that invalid filter expression returns unfiltered data with error log.

        Migrated from test_filters.py.
        """
        # Invalid syntax
        context.config["final_filter"] = 'CHROM === "chr1"'  # Triple equals invalid

        stage = FinalFilteringStage()
        result = stage(context)

        # Should return original dataframe
        pd.testing.assert_frame_equal(result.current_dataframe, context.current_dataframe)

        # Check error was logged
        assert "Failed to apply final filter expression" in caplog.text

    def test_filter_on_computed_columns(self, context):
        """Test filtering on columns that might be computed during analysis.

        Migrated from test_filters.py.
        """
        context.current_dataframe["gene"] = ["BRCA1", "TP53", "EGFR", "KRAS"]
        context.current_dataframe["inheritance_score"] = ["0.95", "0.45", "0.78", "0.82"]
        context.current_dataframe["Inheritance_Confidence"] = ["0.9", "0.3", "0.8", "0.85"]
        context.config["final_filter"] = "inheritance_score > 0.7 and Inheritance_Confidence > 0.75"

        stage = FinalFilteringStage()
        result = stage(context)

        assert len(result.current_dataframe) == 3
        assert set(result.current_dataframe["gene"]) == {"BRCA1", "EGFR", "KRAS"}


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

    @patch("variantcentrifuge.stages.output_stages.finalize_excel_file")
    @patch("variantcentrifuge.stages.output_stages.convert_to_excel")
    def test_excel_generation(self, mock_convert, mock_finalize, context):
        """Test Excel report generation."""
        # Mark dependencies as complete since Excel depends on them
        context.mark_complete("tsv_output")
        context.mark_complete("metadata_generation")
        context.mark_complete("statistics_generation")
        # Set final output path which TSV output would have created
        context.final_output_path = context.workspace.output_dir / "output.tsv"
        context.final_output_path.touch()  # Create the file

        # Configure the mock to return a specific Excel file path
        expected_excel_path = str(context.workspace.output_dir / "output.xlsx")
        mock_convert.return_value = expected_excel_path

        # Create the mock Excel file that convert_to_excel "creates"
        from pathlib import Path

        Path(expected_excel_path).touch()

        stage = ExcelReportStage()
        stage(context)

        # Check convert_to_excel was called
        mock_convert.assert_called_once()
        # Implementation calls with positional args: convert_to_excel(input_file, config)
        call_args = mock_convert.call_args[0]
        assert str(context.final_output_path) == call_args[0]
        assert isinstance(call_args[1], dict)  # Second argument is config

        # Check finalize_excel_file was called
        mock_finalize.assert_called_once()

    def test_skip_if_disabled(self, context):
        """Test skipping when Excel disabled."""
        context.config["xlsx"] = False
        # Mark all dependencies as complete since Excel depends on them
        context.mark_complete("tsv_output")
        context.mark_complete("metadata_generation")
        context.mark_complete("statistics_generation")
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

        # Mock produce_report_json to create the expected JSON files
        def create_json_files(*args, **kwargs):
            report_dir = context.workspace.output_dir / "report"
            report_dir.mkdir(parents=True, exist_ok=True)

            # Create variants.json
            variants_json = report_dir / "variants.json"
            variants_json.write_text('[{"chrom": "chr1", "pos": "100"}]')

            # Create summary.json
            summary_json = report_dir / "summary.json"
            summary_json.write_text('{"num_variants": 1, "num_genes": 1}')

        mock_produce_json.side_effect = create_json_files

        # Mock generate_html_report to create index.html
        def create_html_file(*args, **kwargs):
            report_dir = context.workspace.output_dir / "report"
            index_html = report_dir / "index.html"
            index_html.write_text("<html><body>Test Report</body></html>")

        mock_generate.side_effect = create_html_file

        stage = HTMLReportStage()
        stage(context)

        # Check JSON was produced first
        mock_produce_json.assert_called_once()
        # Check report generation - implementation uses different parameters
        mock_generate.assert_called_once()
        call_kwargs = mock_generate.call_args[1]
        assert "variants_json" in call_kwargs
        assert "summary_json" in call_kwargs
        assert "output_dir" in call_kwargs
        assert "cfg" in call_kwargs

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
            vcf_file="/tmp/input.vcf",
            config_overrides={
                "igv_enabled": True,
                "reference": "hg38",
                "bam_mapping_file": "/tmp/bam_mapping.txt",
                "igv_reference": "hg19",
            },
        )
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.current_dataframe = pd.DataFrame({"CHROM": ["chr1"], "POS": [100]})
        ctx.mark_complete("dataframe_loading")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.generate_igv_report")
    def test_igv_generation(self, mock_generate, context):
        """Test IGV report generation."""
        # Mark TSV output as complete since IGV depends on it
        context.mark_complete("tsv_output")
        # Set final output path which TSV output would have created
        context.final_output_path = context.workspace.output_dir / "output.tsv"
        context.final_output_path.touch()  # Create the file

        # Create a mock BAM mapping file
        bam_mapping_file = context.workspace.output_dir / "bam_mapping.txt"
        bam_mapping_file.write_text("Sample1,/path/to/sample1.bam\n")
        context.config["bam_mapping_file"] = str(bam_mapping_file)

        stage = IGVReportStage()
        stage(context)

        # Check report generation
        mock_generate.assert_called_once()
        call_kwargs = mock_generate.call_args[1]
        assert call_kwargs["variants_tsv"] == str(context.final_output_path)
        assert call_kwargs["bam_mapping_file"] == str(bam_mapping_file)
        assert call_kwargs["igv_reference"] == "hg19"
        assert call_kwargs["integrate_into_main"] is True
        assert "max_workers" in call_kwargs
        assert call_kwargs["igv_flanking"] == 50

    def test_skip_if_disabled(self, context):
        """Test skipping when IGV disabled."""
        context.config["igv_enabled"] = False
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

        stage = MetadataGenerationStage()
        result = stage(context)

        # Check metadata file was created
        assert "metadata" in result.report_paths
        metadata_path = result.report_paths["metadata"]
        assert metadata_path.exists()
        assert metadata_path.suffix == ".tsv"  # Implementation creates TSV, not JSON

        # Read and check contents
        with open(metadata_path) as f:
            content = f.read()

        # Check that TSV format is correct
        assert content.startswith("Parameter\tValue\n")
        assert "Tool\tvariantcentrifuge" in content
        assert "Version\t1.0.0-test" in content

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
