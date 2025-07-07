"""Unit tests for setup stages."""

import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pandas as pd
import pytest

from variantcentrifuge.stages.setup_stages import (
    ConfigurationLoadingStage,
    PhenotypeLoadingStage,
    ScoringConfigLoadingStage,
    PedigreeLoadingStage,
    AnnotationConfigLoadingStage,
    SampleConfigLoadingStage,
)
from tests.mocks.fixtures import create_test_context


class TestConfigurationLoadingStage:
    """Test ConfigurationLoadingStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        return create_test_context(
            gene_name="BRCA1", vcf_file="/tmp/test.vcf", config_overrides={"reference": "hg38"}
        )

    @patch("variantcentrifuge.config.load_config")
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

        # Check config was merged - CLI args take precedence
        assert result.config["reference"] == "hg38"  # Original preserved
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
        assert result.phenotype_data == {"Sample1": "Affected", "Sample2": "Unaffected"}

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

    @patch("variantcentrifuge.scoring.read_scoring_config")
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

    @patch("variantcentrifuge.ped_reader.read_pedigree")
    def test_pedigree_loading(self, mock_read_pedigree, context):
        """Test pedigree file loading."""
        # Mock pedigree loading
        mock_samples = {"Sample1": {"affected": True, "sex": "male"}}
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
        return ctx

    @patch("variantcentrifuge.helpers.get_vcf_samples")
    def test_sample_loading(self, mock_get_samples, context):
        """Test sample loading from VCF."""
        # Mock sample extraction
        mock_get_samples.return_value = ["Sample1", "Sample2", "Sample3"]

        stage = SampleConfigLoadingStage()
        result = stage(context)

        # Check samples were loaded
        mock_get_samples.assert_called_once_with("/tmp/test.vcf")
        assert result.vcf_samples == ["Sample1", "Sample2", "Sample3"]

    @patch("variantcentrifuge.helpers.get_vcf_samples")
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
        assert result.case_samples == ["Sample1", "Sample2"]
        assert result.control_samples == ["Sample3"]
        assert set(result.vcf_samples) == {
            "Sample1",
            "Sample2",
            "Sample3",
        }  # Only specified samples

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = SampleConfigLoadingStage()
        assert stage.parallel_safe is True

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = SampleConfigLoadingStage()
        assert stage.dependencies == {"configuration_loading"}
