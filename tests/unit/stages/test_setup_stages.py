"""Unit tests for setup stages."""

import json
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.mocks.fixtures import create_test_context
from variantcentrifuge.stages.setup_stages import (
    AnnotationConfigLoadingStage,
    ConfigurationLoadingStage,
    PedigreeLoadingStage,
    PhenotypeCaseControlAssignmentStage,
    PhenotypeLoadingStage,
    SampleConfigLoadingStage,
    ScoringConfigLoadingStage,
)


class TestConfigurationLoadingStage:
    """Test ConfigurationLoadingStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        return create_test_context(
            gene_name="BRCA1", vcf_file="/tmp/test.vcf", config_overrides={"reference": "hg38"}
        )

    @patch("variantcentrifuge.stages.setup_stages.load_config")
    def test_config_loading(self, mock_load_config, context):
        """Test configuration loading."""
        # Mock config loading - return a full config
        mock_load_config.return_value = {
            "filters": {"rare": "AF < 0.01"},
            "fields": ["CHROM", "POS", "REF", "ALT"],
            "reference": "hg37",  # Different from context
        }

        stage = ConfigurationLoadingStage()
        result = stage(context)

        # Check config was loaded
        mock_load_config.assert_called_once()

        # Check config was merged - loaded config takes precedence for existing keys
        assert result.config["reference"] == "hg37"  # From loaded config (not context)
        assert "filters" in result.config  # New added
        assert "fields" in result.config  # New added

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = ConfigurationLoadingStage()
        assert stage.dependencies == set()  # No dependencies

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = ConfigurationLoadingStage()
        assert stage.parallel_safe is True


class TestPhenotypeLoadingStage:
    """Test PhenotypeLoadingStage."""

    @pytest.fixture
    def context(self):
        """Create test context with phenotype file."""
        ctx = create_test_context()
        # Create a temp phenotype file
        phenotype_file = Path(tempfile.mktemp(suffix=".tsv"))
        phenotype_file.write_text("sample\tphenotype\nSample1\tAffected\nSample2\tUnaffected\n")
        ctx.config["phenotype_file"] = str(phenotype_file)
        ctx.config["phenotype_sample_column"] = "sample"
        ctx.config["phenotype_value_column"] = "phenotype"
        return ctx

    def test_phenotype_loading(self, context):
        """Test phenotype file loading."""
        stage = PhenotypeLoadingStage()
        result = stage(context)

        # Check phenotypes were loaded
        assert result.phenotype_data == {"Sample1": {"Affected"}, "Sample2": {"Unaffected"}}

    def test_no_phenotype_file(self, context):
        """Test when no phenotype file specified."""
        context.config["phenotype_file"] = None

        stage = PhenotypeLoadingStage()
        result = stage(context)

        # Should return unchanged
        assert result.phenotype_data is None

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = PhenotypeLoadingStage()
        assert stage.parallel_safe is True


class TestScoringConfigLoadingStage:
    """Test ScoringConfigLoadingStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context()
        ctx.config["scoring_config_path"] = "/tmp/scoring"
        return ctx

    @patch("variantcentrifuge.stages.setup_stages.read_scoring_config")
    def test_scoring_config_loading(self, mock_read_scoring_config, context):
        """Test scoring configuration loading."""
        # Mock scoring config
        mock_scoring_config = {"formulas": {"score": "AF * 10"}}
        mock_read_scoring_config.return_value = mock_scoring_config

        stage = ScoringConfigLoadingStage()
        result = stage(context)

        # Check config was loaded
        mock_read_scoring_config.assert_called_once_with("/tmp/scoring")
        assert result.scoring_config == mock_scoring_config

    def test_no_scoring_config(self, context):
        """Test when no scoring config specified."""
        context.config["scoring_config_path"] = None

        stage = ScoringConfigLoadingStage()
        result = stage(context)

        # Should return unchanged
        assert result.scoring_config is None

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = ScoringConfigLoadingStage()
        assert stage.parallel_safe is True


class TestPedigreeLoadingStage:
    """Test PedigreeLoadingStage."""

    @pytest.fixture
    def context(self):
        """Create test context with PED file."""
        ctx = create_test_context()
        # Create a temp PED file
        ped_file = Path(tempfile.mktemp(suffix=".ped"))
        ped_file.write_text("FAM1\tSample1\t0\t0\t1\t2\n")  # Affected male
        ctx.config["ped_file"] = str(ped_file)
        return ctx

    @patch("variantcentrifuge.stages.setup_stages.read_pedigree")
    def test_pedigree_loading(self, mock_read_pedigree, context):
        """Test pedigree file loading."""
        # Mock pedigree loading with correct structure
        mock_samples = {
            "Sample1": {
                "family_id": "FAM1",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",  # 1=male, 2=female
                "affected_status": "2",  # 2=affected, 1=unaffected
            }
        }
        mock_read_pedigree.return_value = mock_samples

        stage = PedigreeLoadingStage()
        result = stage(context)

        # Check pedigree was loaded
        mock_read_pedigree.assert_called_once_with(context.config["ped_file"])
        assert result.pedigree_data == mock_samples

    def test_no_pedigree_file(self, context):
        """Test when no pedigree file specified."""
        context.config["ped_file"] = None

        stage = PedigreeLoadingStage()
        result = stage(context)

        # Should return unchanged
        assert result.pedigree_data is None

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = PedigreeLoadingStage()
        assert stage.parallel_safe is True


class TestAnnotationConfigLoadingStage:
    """Test AnnotationConfigLoadingStage."""

    @pytest.fixture
    def context(self):
        """Create test context with annotations."""
        ctx = create_test_context()
        # Create temp BED file
        bed_file = Path(tempfile.mktemp(suffix=".bed"))
        bed_file.write_text("chr1\t100\t200\tRegion1\n")

        # Create temp gene list
        gene_list = Path(tempfile.mktemp(suffix=".txt"))
        gene_list.write_text("BRCA1\nBRCA2\n")

        ctx.config["annotate_bed"] = [str(bed_file)]
        ctx.config["annotate_gene_list"] = [str(gene_list)]
        return ctx

    def test_annotation_loading(self, context):
        """Test annotation file loading."""
        stage = AnnotationConfigLoadingStage()
        result = stage(context)

        # Check annotation configs were created
        assert hasattr(result, "annotation_configs")
        assert "bed_files" in result.annotation_configs
        assert "gene_lists" in result.annotation_configs
        assert len(result.annotation_configs["bed_files"]) == 1
        assert len(result.annotation_configs["gene_lists"]) == 1

    def test_json_gene_annotation(self, context):
        """Test JSON gene annotation loading."""
        # Create temp JSON file
        json_file = Path(tempfile.mktemp(suffix=".json"))
        json_file.write_text(
            json.dumps(
                [{"symbol": "BRCA1", "panel": "Cancer"}, {"symbol": "BRCA2", "panel": "Cancer"}]
            )
        )

        context.config["annotate_bed"] = []
        context.config["annotate_gene_list"] = []
        context.config["annotate_json_genes"] = [str(json_file)]
        context.config["json_gene_mapping"] = {"identifier": "symbol", "dataFields": ["panel"]}

        stage = AnnotationConfigLoadingStage()
        result = stage(context)

        # Check JSON annotation was loaded
        assert hasattr(result, "annotation_configs")
        assert "json_genes" in result.annotation_configs
        assert result.annotation_configs["json_genes"] == [str(json_file)]
        assert "json_mapping" in result.annotation_configs

    def test_no_annotations(self, context):
        """Test when no annotations specified."""
        context.config["annotate_bed"] = []
        context.config["annotate_gene_list"] = []
        context.config["annotate_json_genes"] = []

        stage = AnnotationConfigLoadingStage()
        result = stage(context)

        # Should have empty dict
        assert result.annotation_configs == {}

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = AnnotationConfigLoadingStage()
        assert stage.parallel_safe is True


class TestSampleConfigLoadingStage:
    """Test SampleConfigLoadingStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context(vcf_file="/tmp/test.vcf")
        # Mark configuration_loading as complete since SampleConfigLoadingStage depends on it
        ctx.mark_complete("configuration_loading")
        return ctx

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_sample_loading(self, mock_get_samples, context):
        """Test sample loading from VCF."""
        # Mock sample extraction
        mock_get_samples.return_value = ["Sample1", "Sample2", "Sample3"]

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check samples were loaded
        mock_get_samples.assert_called_once_with("/tmp/test.vcf")
        assert result.vcf_samples == ["Sample1", "Sample2", "Sample3"]

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_sample_loading_with_subsetting(self, mock_get_samples, context):
        """Test sample loading with subsetting."""
        # Set up sample files
        case_file = Path(tempfile.mktemp())
        case_file.write_text("Sample1\nSample2\n")
        control_file = Path(tempfile.mktemp())
        control_file.write_text("Sample3\n")

        context.config["case_samples_file"] = str(case_file)
        context.config["control_samples_file"] = str(control_file)

        # Mock full sample list
        mock_get_samples.return_value = ["Sample1", "Sample2", "Sample3", "Sample4"]

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check sample subsetting
        assert result.config["case_samples"] == ["Sample1", "Sample2"]
        assert result.config["control_samples"] == ["Sample3"]
        # VCF samples should be all samples from VCF, not filtered
        assert result.vcf_samples == ["Sample1", "Sample2", "Sample3", "Sample4"]

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = SampleConfigLoadingStage()
        assert stage.parallel_safe is True

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = SampleConfigLoadingStage()
        assert stage.dependencies == {"configuration_loading"}

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_remove_sample_substring_basic(self, mock_get_samples, context):
        """Test basic sample substring removal functionality."""
        # Mock sample extraction with suffixes
        mock_get_samples.return_value = ["Sample1_suffix", "Sample2_suffix", "Sample3_suffix"]

        # Configure substring removal
        context.config["remove_sample_substring"] = "_suffix"

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check substring was removed
        expected_samples = ["Sample1", "Sample2", "Sample3"]
        assert result.vcf_samples == expected_samples
        mock_get_samples.assert_called_once_with("/tmp/test.vcf")

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_remove_sample_substring_partial_match(self, mock_get_samples, context):
        """Test substring removal when only some samples contain the substring."""
        # Mock sample extraction with mixed suffixes
        mock_get_samples.return_value = ["Sample1_suffix", "Sample2", "Sample3_suffix", "Sample4"]

        # Configure substring removal
        context.config["remove_sample_substring"] = "_suffix"

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check substring was removed only where present
        expected_samples = ["Sample1", "Sample2", "Sample3", "Sample4"]
        assert result.vcf_samples == expected_samples

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_remove_sample_substring_no_match(self, mock_get_samples, context):
        """Test substring removal when no samples contain the substring."""
        # Mock sample extraction without target substring
        mock_get_samples.return_value = ["Sample1", "Sample2", "Sample3"]

        # Configure substring removal
        context.config["remove_sample_substring"] = "_nonexistent"

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check samples remain unchanged
        expected_samples = ["Sample1", "Sample2", "Sample3"]
        assert result.vcf_samples == expected_samples

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_remove_sample_substring_empty_string(self, mock_get_samples, context):
        """Test substring removal with empty string (should be ignored)."""
        # Mock sample extraction
        mock_get_samples.return_value = ["Sample1", "Sample2", "Sample3"]

        # Configure empty substring removal
        context.config["remove_sample_substring"] = ""

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check samples remain unchanged
        expected_samples = ["Sample1", "Sample2", "Sample3"]
        assert result.vcf_samples == expected_samples

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_remove_sample_substring_whitespace_only(self, mock_get_samples, context):
        """Test substring removal with whitespace-only string (should be ignored)."""
        # Mock sample extraction
        mock_get_samples.return_value = ["Sample1", "Sample2", "Sample3"]

        # Configure whitespace-only substring removal
        context.config["remove_sample_substring"] = "   "

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check samples remain unchanged
        expected_samples = ["Sample1", "Sample2", "Sample3"]
        assert result.vcf_samples == expected_samples

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_remove_sample_substring_complex_pattern(self, mock_get_samples, context):
        """Test substring removal with complex patterns."""
        # Mock sample extraction with complex patterns
        mock_get_samples.return_value = [
            "SAMPLE_001_processed",
            "SAMPLE_002_processed",
            "SAMPLE_003_raw",
        ]

        # Configure substring removal for processing suffix
        context.config["remove_sample_substring"] = "_processed"

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check substring was removed appropriately
        expected_samples = ["SAMPLE_001", "SAMPLE_002", "SAMPLE_003_raw"]
        assert result.vcf_samples == expected_samples

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_remove_sample_substring_multiple_occurrences(self, mock_get_samples, context):
        """Test substring removal when substring appears multiple times in a sample name."""
        # Mock sample extraction with repeated substrings
        mock_get_samples.return_value = ["test_sample_test", "another_test_sample", "final_test"]

        # Configure substring removal
        context.config["remove_sample_substring"] = "_test"

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check all occurrences were removed (Python's replace removes all)
        expected_samples = ["test_sample", "another_sample", "final"]
        assert result.vcf_samples == expected_samples

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_remove_sample_substring_with_case_samples(self, mock_get_samples, context):
        """Test substring removal works with case/control sample loading."""
        # Mock sample extraction with suffixes
        mock_get_samples.return_value = ["Sample1_suffix", "Sample2_suffix", "Sample3_suffix"]

        # Set up case samples file
        case_file = Path(tempfile.mktemp())
        case_file.write_text("Sample1_suffix\nSample2_suffix\n")
        context.config["case_samples_file"] = str(case_file)

        # Configure substring removal
        context.config["remove_sample_substring"] = "_suffix"

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check VCF samples had substring removed
        expected_samples = ["Sample1", "Sample2", "Sample3"]
        assert result.vcf_samples == expected_samples

        # Check case samples were loaded correctly (before substring removal)
        assert result.config["case_samples"] == ["Sample1_suffix", "Sample2_suffix"]

    @patch("variantcentrifuge.stages.setup_stages.get_vcf_samples")
    def test_remove_sample_substring_none_config(self, mock_get_samples, context):
        """Test behavior when remove_sample_substring is None."""
        # Mock sample extraction
        mock_get_samples.return_value = ["Sample1", "Sample2", "Sample3"]

        # Configure None substring removal
        context.config["remove_sample_substring"] = None

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check samples remain unchanged
        expected_samples = ["Sample1", "Sample2", "Sample3"]
        assert result.vcf_samples == expected_samples


class TestPhenotypeCaseControlAssignmentStage:
    """Test PhenotypeCaseControlAssignmentStage."""

    @pytest.fixture
    def context(self):
        """Create test context with phenotype configuration."""
        context = create_test_context(
            gene_name="BRCA1",
            vcf_file="/tmp/test.vcf",
            config_overrides={
                "phenotype_file": "/tmp/phenotypes.csv",
                "phenotype_sample_column": "sample_id",
                "phenotype_value_column": "phenotype",
                "case_phenotypes": ["HP:0000001", "HP:0000002"],
                "control_phenotypes": [],
            },
        )
        # Mock VCF samples
        context.vcf_samples = ["S001", "S002", "S003", "S004"]
        return context

    @patch(
        "variantcentrifuge.stages.setup_stages."
        "PhenotypeCaseControlAssignmentStage._compute_samples_from_phenotype_file"
    )
    def test_phenotype_based_assignment(self, mock_compute, context):
        """Test phenotype-based case/control assignment."""
        # Mock the computation result
        mock_compute.return_value = (["S001", "S002"], ["S003", "S004"])

        # Mark dependencies as complete
        context.mark_complete("phenotype_loading")
        context.mark_complete("sample_config_loading")

        stage = PhenotypeCaseControlAssignmentStage()
        result = stage(context)

        # Check assignment was applied
        assert result.config["case_samples"] == ["S001", "S002"]
        assert result.config["control_samples"] == ["S003", "S004"]

        # Check method was called with correct parameters
        mock_compute.assert_called_once_with(
            "/tmp/phenotypes.csv",
            "sample_id",
            "phenotype",
            ["HP:0000001", "HP:0000002"],
            [],
            ["S001", "S002", "S003", "S004"],
            "",
        )

    def test_skip_if_explicit_assignments_exist(self, context):
        """Test that assignment is skipped if explicit case/control samples exist."""
        # Set explicit assignments
        context.config["case_samples"] = ["S001"]
        context.config["control_samples"] = ["S002"]

        # Mark dependencies as complete
        context.mark_complete("phenotype_loading")
        context.mark_complete("sample_config_loading")

        stage = PhenotypeCaseControlAssignmentStage()
        result = stage(context)

        # Should remain unchanged
        assert result.config["case_samples"] == ["S001"]
        assert result.config["control_samples"] == ["S002"]

    def test_skip_if_no_phenotype_criteria(self, context):
        """Test that assignment is skipped if no phenotype criteria provided."""
        # Remove phenotype criteria
        context.config["case_phenotypes"] = []
        context.config["control_phenotypes"] = []

        # Mark dependencies as complete
        context.mark_complete("phenotype_loading")
        context.mark_complete("sample_config_loading")

        stage = PhenotypeCaseControlAssignmentStage()
        result = stage(context)

        # Should not have case/control assignments
        assert "case_samples" not in result.config
        assert "control_samples" not in result.config

    @patch("pandas.read_csv")
    def test_compute_samples_from_phenotype_file(self, mock_read_csv, context):
        """Test the phenotype file processing logic."""
        import pandas as pd

        # Mock phenotype data
        mock_df = pd.DataFrame(
            {
                "sample_id": [1001, 1002, 1003, 1004],
                "phenotype": ["HP:0000001", "HP:0000002", "HP:0000003", "HP:0000001"],
            }
        )
        mock_read_csv.return_value = mock_df

        # Mock VCF samples as strings
        vcf_samples = ["1001", "1002", "1003", "1004"]

        stage = PhenotypeCaseControlAssignmentStage()
        case_samples, control_samples = stage._compute_samples_from_phenotype_file(
            "/tmp/phenotypes.csv",
            "sample_id",
            "phenotype",
            ["HP:0000001", "HP:0000002"],
            [],
            vcf_samples,
            "",
        )

        # Should find samples with matching phenotypes
        assert set(case_samples) == {
            "1001",
            "1002",
            "1004",
        }  # 1001 and 1004 have HP:0000001, 1002 has HP:0000002
        assert set(control_samples) == {"1003"}  # Only non-case sample

    @patch("pandas.read_csv")
    def test_compute_samples_with_substring_removal(self, mock_read_csv, context):
        """Test phenotype assignment with sample substring removal."""
        import pandas as pd

        # Mock phenotype data with suffixes
        mock_df = pd.DataFrame(
            {
                "sample_id": ["1001_suffix", "1002_suffix", "1003_suffix", "1004_suffix"],
                "phenotype": ["HP:0000001", "HP:0000002", "HP:0000003", "HP:0000001"],
            }
        )
        mock_read_csv.return_value = mock_df

        # Mock VCF samples without suffixes (already processed)
        vcf_samples = ["1001", "1002", "1003", "1004"]

        stage = PhenotypeCaseControlAssignmentStage()
        case_samples, control_samples = stage._compute_samples_from_phenotype_file(
            "/tmp/phenotypes.csv",
            "sample_id",
            "phenotype",
            ["HP:0000001", "HP:0000002"],
            [],
            vcf_samples,
            "_suffix",
        )

        # Should correctly match after substring removal
        assert set(case_samples) == {"1001", "1002", "1004"}
        assert set(control_samples) == {"1003"}

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = PhenotypeCaseControlAssignmentStage()
        assert stage.dependencies == {"phenotype_loading", "sample_config_loading"}

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = PhenotypeCaseControlAssignmentStage()
        assert stage.parallel_safe is True
