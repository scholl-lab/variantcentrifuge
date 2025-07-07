"""Unit tests for Workspace."""

import tempfile
from pathlib import Path

import pytest

from variantcentrifuge.pipeline_core import Workspace


class TestWorkspace:
    """Test suite for Workspace."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def workspace(self, temp_dir):
        """Create a test Workspace."""
        return Workspace(temp_dir, "test_base")

    def test_initialization(self, temp_dir):
        """Test Workspace initialization."""
        workspace = Workspace(temp_dir, "test_run")

        assert workspace.output_dir == temp_dir
        assert workspace.base_name == "test_run"
        assert workspace.intermediate_dir == temp_dir / "intermediate"
        assert workspace.intermediate_dir.exists()
        assert workspace.temp_dir.exists()
        assert "vc_temp_" in workspace.temp_dir.name

    def test_get_output_path(self, workspace):
        """Test output path generation."""
        path1 = workspace.get_output_path(".final")
        assert path1 == workspace.output_dir / "test_base.final.tsv"

        path2 = workspace.get_output_path(".filtered", ".vcf")
        assert path2 == workspace.output_dir / "test_base.filtered.vcf"

        path3 = workspace.get_output_path("_processed", ".xlsx")
        assert path3 == workspace.output_dir / "test_base_processed.xlsx"

    def test_get_intermediate_path(self, workspace):
        """Test intermediate path generation."""
        path = workspace.get_intermediate_path("variants.vcf.gz")
        assert path == workspace.intermediate_dir / "variants.vcf.gz"
        assert path.parent.exists()

    def test_get_temp_path(self, workspace):
        """Test temporary path generation."""
        path = workspace.get_temp_path("chunk_0.tsv")
        assert path == workspace.temp_dir / "chunk_0.tsv"
        assert path.parent.exists()

    def test_get_timestamped_path(self, workspace):
        """Test timestamped path generation."""
        # Test with default directory
        path1 = workspace.get_timestamped_path("report.html")
        assert path1.parent == workspace.output_dir
        assert workspace.timestamp in str(path1)
        assert path1.name.startswith("report_")
        assert path1.name.endswith(".html")

        # Test with custom directory
        custom_dir = workspace.output_dir / "custom"
        custom_dir.mkdir()
        path2 = workspace.get_timestamped_path("data.tsv", custom_dir)
        assert path2.parent == custom_dir
        assert workspace.timestamp in str(path2)

    def test_create_subdirectory(self, workspace):
        """Test subdirectory creation."""
        # In output directory
        subdir1 = workspace.create_subdirectory("reports")
        assert subdir1 == workspace.output_dir / "reports"
        assert subdir1.exists()

        # In intermediate directory
        subdir2 = workspace.create_subdirectory("temp_data", in_intermediate=True)
        assert subdir2 == workspace.intermediate_dir / "temp_data"
        assert subdir2.exists()

        # Test idempotency
        subdir3 = workspace.create_subdirectory("reports")
        assert subdir3 == subdir1

    def test_cleanup(self, workspace):
        """Test cleanup functionality."""
        # Create some files
        temp_file = workspace.get_temp_path("test.txt")
        temp_file.write_text("test")

        intermediate_file = workspace.get_intermediate_path("data.tsv")
        intermediate_file.write_text("data")

        output_file = workspace.get_output_path(".final")
        output_file.write_text("output")

        # Cleanup keeping intermediates
        workspace.cleanup(keep_intermediates=True)
        assert not workspace.temp_dir.exists()
        assert workspace.intermediate_dir.exists()
        assert intermediate_file.exists()
        assert output_file.exists()

        # Cleanup everything
        workspace.cleanup(keep_intermediates=False)
        assert not workspace.temp_dir.exists()
        assert not workspace.intermediate_dir.exists()
        assert output_file.exists()  # Output files are never deleted

    def test_get_archive_path(self, workspace):
        """Test archive path generation."""
        archive_path = workspace.get_archive_path()
        assert archive_path.parent == workspace.output_dir.parent
        assert archive_path.name.startswith(workspace.output_dir.name)
        assert archive_path.name.endswith(".tar.gz")
        assert workspace.timestamp in archive_path.name

    def test_list_outputs(self, workspace):
        """Test listing output files."""
        # Initially empty
        assert len(workspace.list_outputs()) == 1  # intermediate dir

        # Create some files
        file1 = workspace.get_output_path(".final")
        file1.write_text("data1")

        file2 = workspace.get_output_path(".stats", ".txt")
        file2.write_text("data2")

        outputs = workspace.list_outputs()
        assert len(outputs) == 3  # 2 files + intermediate dir
        assert file1 in outputs
        assert file2 in outputs

    def test_list_intermediates(self, workspace):
        """Test listing intermediate files."""
        # Initially has temp directory
        intermediates = workspace.list_intermediates()
        assert len(intermediates) == 1  # temp dir

        # Create some files
        file1 = workspace.get_intermediate_path("step1.tsv")
        file1.write_text("data1")

        file2 = workspace.get_intermediate_path("step2.vcf")
        file2.write_text("data2")

        intermediates = workspace.list_intermediates()
        assert len(intermediates) == 3  # 2 files + temp dir
        assert file1 in intermediates
        assert file2 in intermediates

    def test_context_manager(self, temp_dir):
        """Test workspace as context manager."""
        with Workspace(temp_dir, "test_context") as workspace:
            # Create temp file
            temp_file = workspace.get_temp_path("test.txt")
            temp_file.write_text("test")
            assert temp_file.exists()

            temp_dir_path = workspace.temp_dir

        # Temp directory should be cleaned up
        assert not temp_dir_path.exists()

    def test_repr(self, workspace):
        """Test string representation."""
        repr_str = repr(workspace)
        assert "Workspace" in repr_str
        assert str(workspace.output_dir) in repr_str
        assert workspace.base_name in repr_str
