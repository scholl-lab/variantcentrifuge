"""
Workspace - Centralized file path management for pipeline runs.

This module provides the Workspace class that manages all file paths
during a pipeline run, including output, intermediate, and temporary files.
"""

import logging
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class Workspace:
    """Manages all file paths for a pipeline run.

    The Workspace class provides centralized management of file paths,
    ensuring consistent naming and organization of output files. It handles:
    - Output directory creation
    - Intermediate file paths
    - Temporary file management
    - Automatic cleanup

    Attributes
    ----------
    output_dir : Path
        Main output directory
    base_name : str
        Base name for generated files
    timestamp : str
        Timestamp string for unique file naming
    intermediate_dir : Path
        Directory for intermediate files
    temp_dir : Path
        Directory for temporary files (auto-cleaned)
    """

    def __init__(self, output_dir: Path, base_name: str):
        """Initialize workspace with output directory and base name.

        Parameters
        ----------
        output_dir : Path
            Main output directory path
        base_name : str
            Base name for generated files
        """
        self.output_dir = Path(output_dir)
        self.base_name = base_name
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Create directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.intermediate_dir = self.output_dir / "intermediate"
        self.intermediate_dir.mkdir(exist_ok=True)

        # Create temp directory
        self.temp_dir = Path(tempfile.mkdtemp(prefix="vc_temp_", dir=self.intermediate_dir))

        logger.debug(f"Workspace initialized: output_dir={self.output_dir}")
        logger.debug(f"Temporary directory: {self.temp_dir}")

    def get_output_path(self, suffix: str, extension: str = ".tsv") -> Path:
        """Generate output file path with suffix and extension.

        Parameters
        ----------
        suffix : str
            Suffix to append to base name (e.g., ".final", ".filtered")
        extension : str
            File extension including dot (default: ".tsv")

        Returns
        -------
        Path
            Full path for the output file
        """
        filename = f"{self.base_name}{suffix}{extension}"
        return self.output_dir / filename

    def get_intermediate_path(self, name: str) -> Path:
        """Generate intermediate file path.

        Parameters
        ----------
        name : str
            Name of the intermediate file

        Returns
        -------
        Path
            Full path in the intermediate directory
        """
        return self.intermediate_dir / name

    def get_temp_path(self, name: str) -> Path:
        """Generate temporary file path.

        Parameters
        ----------
        name : str
            Name of the temporary file

        Returns
        -------
        Path
            Full path in the temporary directory
        """
        return self.temp_dir / name

    def get_timestamped_path(self, name: str, directory: Optional[Path] = None) -> Path:
        """Generate a timestamped file path.

        Parameters
        ----------
        name : str
            Base name for the file
        directory : Path, optional
            Directory to place file in (default: output_dir)

        Returns
        -------
        Path
            Path with timestamp inserted before extension
        """
        if directory is None:
            directory = self.output_dir

        # Split name and extension
        path = Path(name)
        stem = path.stem
        suffix = path.suffix

        # Insert timestamp
        timestamped_name = f"{stem}_{self.timestamp}{suffix}"
        return directory / timestamped_name

    def create_subdirectory(self, name: str, in_intermediate: bool = False) -> Path:
        """Create a subdirectory in output or intermediate directory.

        Parameters
        ----------
        name : str
            Name of the subdirectory
        in_intermediate : bool
            If True, create in intermediate dir, else in output dir

        Returns
        -------
        Path
            Path to the created subdirectory
        """
        parent = self.intermediate_dir if in_intermediate else self.output_dir
        subdir = parent / name
        subdir.mkdir(exist_ok=True)
        return subdir

    def cleanup(self, keep_intermediates: bool = False) -> None:
        """Clean up temporary files and optionally intermediate files.

        Parameters
        ----------
        keep_intermediates : bool
            If False, also remove intermediate directory
        """
        try:
            # Always clean temp directory
            if self.temp_dir.exists():
                shutil.rmtree(self.temp_dir)
                logger.debug(f"Cleaned up temporary directory: {self.temp_dir}")

            # Optionally clean intermediate directory
            if not keep_intermediates and self.intermediate_dir.exists():
                # Only remove if it's actually in our output directory
                if self.intermediate_dir.parent == self.output_dir:
                    shutil.rmtree(self.intermediate_dir)
                    logger.debug(f"Cleaned up intermediate directory: {self.intermediate_dir}")

        except Exception as e:
            logger.warning(f"Error during cleanup: {e}")

    def get_archive_path(self) -> Path:
        """Get path for results archive file.

        Returns
        -------
        Path
            Path for archive file in parent directory
        """
        archive_name = f"{self.output_dir.name}_{self.timestamp}.tar.gz"
        return self.output_dir.parent / archive_name

    def list_outputs(self) -> list:
        """List all files in the output directory.

        Returns
        -------
        list
            List of Path objects for all output files
        """
        return list(self.output_dir.glob("*"))

    def list_intermediates(self) -> list:
        """List all files in the intermediate directory.

        Returns
        -------
        list
            List of Path objects for all intermediate files
        """
        if self.intermediate_dir.exists():
            return list(self.intermediate_dir.glob("*"))
        return []

    def __repr__(self) -> str:
        """Return string representation of the workspace."""
        return f"Workspace(output_dir='{self.output_dir}', " f"base_name='{self.base_name}')"

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with automatic cleanup."""
        self.cleanup(keep_intermediates=True)  # Keep intermediates by default
