"""Unit tests for file-based case/control assignment functionality."""

import tempfile
from pathlib import Path

import pytest

from tests.mocks.fixtures import create_test_context
from variantcentrifuge.stages.setup_stages import PhenotypeCaseControlAssignmentStage

# from unittest.mock import patch  # Unused import


class TestFileBasedAssignment:
    """Test file-based case/control assignment."""

    @pytest.fixture
    def context(self):
        """Create test context for file-based assignment."""
        context = create_test_context(
            gene_name="BRCA1",
            vcf_file="/tmp/test.vcf",
            config_overrides={},
        )
        # Mock VCF samples (already processed with substring removal)
        context.vcf_samples = ["S001", "S002", "S003", "S004", "S005"]
        return context

    def test_file_based_assignment_both_files(self, context):
        """Test file-based assignment with both case and control files."""
        # Create temporary files
        case_file = Path(tempfile.mktemp())
        control_file = Path(tempfile.mktemp())

        try:
            # Write case samples (with suffix that will be removed)
            case_file.write_text("S001_suffix\nS002_suffix\n")
            # Write control samples
            control_file.write_text("S003_suffix\nS004_suffix\n")

            context.config["case_samples_file"] = str(case_file)
            context.config["control_samples_file"] = str(control_file)
            context.config["remove_sample_substring"] = "_suffix"

            # Mark dependencies as complete
            context.mark_complete("phenotype_loading")
            context.mark_complete("sample_config_loading")

            stage = PhenotypeCaseControlAssignmentStage()
            result = stage(context)

            # Check assignments
            assert set(result.config["case_samples"]) == {"S001", "S002"}
            assert set(result.config["control_samples"]) == {"S003", "S004"}

        finally:
            case_file.unlink(missing_ok=True)
            control_file.unlink(missing_ok=True)

    def test_file_based_assignment_case_only(self, context):
        """Test file-based assignment with only case file (controls auto-assigned)."""
        case_file = Path(tempfile.mktemp())

        try:
            case_file.write_text("S001\nS002\n")

            context.config["case_samples_file"] = str(case_file)
            context.config["remove_sample_substring"] = ""

            # Mark dependencies as complete
            context.mark_complete("phenotype_loading")
            context.mark_complete("sample_config_loading")

            stage = PhenotypeCaseControlAssignmentStage()
            result = stage(context)

            # Cases from file, controls are all other VCF samples
            assert set(result.config["case_samples"]) == {"S001", "S002"}
            assert set(result.config["control_samples"]) == {"S003", "S004", "S005"}

        finally:
            case_file.unlink(missing_ok=True)

    def test_file_based_assignment_control_only(self, context):
        """Test file-based assignment with only control file (cases auto-assigned)."""
        control_file = Path(tempfile.mktemp())

        try:
            control_file.write_text("S004\nS005\n")

            context.config["control_samples_file"] = str(control_file)
            context.config["remove_sample_substring"] = ""

            # Mark dependencies as complete
            context.mark_complete("phenotype_loading")
            context.mark_complete("sample_config_loading")

            stage = PhenotypeCaseControlAssignmentStage()
            result = stage(context)

            # Controls from file, cases are all other VCF samples
            assert set(result.config["case_samples"]) == {"S001", "S002", "S003"}
            assert set(result.config["control_samples"]) == {"S004", "S005"}

        finally:
            control_file.unlink(missing_ok=True)

    def test_file_based_assignment_with_substring_removal(self, context):
        """Test file-based assignment with substring removal applied."""
        case_file = Path(tempfile.mktemp())

        try:
            # File contains samples with suffix
            case_file.write_text("S001-N1-DNA1-WES1\nS002-N1-DNA1-WES1\nS999-N1-DNA1-WES1\n")

            context.config["case_samples_file"] = str(case_file)
            context.config["remove_sample_substring"] = "-N1-DNA1-WES1"

            # Mark dependencies as complete
            context.mark_complete("phenotype_loading")
            context.mark_complete("sample_config_loading")

            stage = PhenotypeCaseControlAssignmentStage()
            result = stage(context)

            # Should only include samples that exist in VCF after substring removal
            # S999 should be excluded because it's not in VCF samples
            assert set(result.config["case_samples"]) == {"S001", "S002"}
            assert set(result.config["control_samples"]) == {"S003", "S004", "S005"}

        finally:
            case_file.unlink(missing_ok=True)

    def test_file_based_assignment_missing_file(self, context):
        """Test file-based assignment with missing file."""
        context.config["case_samples_file"] = "/nonexistent/file.txt"

        # Mark dependencies as complete
        context.mark_complete("phenotype_loading")
        context.mark_complete("sample_config_loading")

        stage = PhenotypeCaseControlAssignmentStage()
        result = stage(context)

        # Should return original context without modifications
        assert "case_samples" not in result.config
        assert "control_samples" not in result.config

    def test_file_based_assignment_empty_file(self, context):
        """Test file-based assignment with empty file."""
        case_file = Path(tempfile.mktemp())

        try:
            case_file.write_text("")  # Empty file

            context.config["case_samples_file"] = str(case_file)

            # Mark dependencies as complete
            context.mark_complete("phenotype_loading")
            context.mark_complete("sample_config_loading")

            stage = PhenotypeCaseControlAssignmentStage()
            result = stage(context)

            # All VCF samples should become controls
            assert set(result.config["case_samples"]) == set()
            assert set(result.config["control_samples"]) == {"S001", "S002", "S003", "S004", "S005"}

        finally:
            case_file.unlink(missing_ok=True)

    def test_file_based_assignment_no_vcf_samples(self, context):
        """Test file-based assignment when no VCF samples available."""
        case_file = Path(tempfile.mktemp())

        try:
            case_file.write_text("S001\nS002\n")

            context.config["case_samples_file"] = str(case_file)
            context.vcf_samples = []  # No VCF samples

            # Mark dependencies as complete
            context.mark_complete("phenotype_loading")
            context.mark_complete("sample_config_loading")

            stage = PhenotypeCaseControlAssignmentStage()
            result = stage(context)

            # Should return original context without modifications
            assert "case_samples" not in result.config
            assert "control_samples" not in result.config

        finally:
            case_file.unlink(missing_ok=True)

    def test_file_based_takes_priority_over_phenotype(self, context):
        """Test that file-based assignment takes priority over phenotype-based."""
        case_file = Path(tempfile.mktemp())

        try:
            case_file.write_text("S001\nS002\n")

            # Set both file-based and phenotype-based config
            context.config["case_samples_file"] = str(case_file)
            context.config["case_phenotypes"] = ["HP:0000001"]
            context.config["phenotype_file"] = "/tmp/phenotypes.csv"
            context.config["phenotype_sample_column"] = "sample_id"
            context.config["phenotype_value_column"] = "phenotype"

            # Mark dependencies as complete
            context.mark_complete("phenotype_loading")
            context.mark_complete("sample_config_loading")

            stage = PhenotypeCaseControlAssignmentStage()
            result = stage(context)

            # Should use file-based assignment, not phenotype-based
            assert set(result.config["case_samples"]) == {"S001", "S002"}
            assert set(result.config["control_samples"]) == {"S003", "S004", "S005"}

        finally:
            case_file.unlink(missing_ok=True)
