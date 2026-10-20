"""Unit tests for release version synchronization and metadata consistency."""

import json
from pathlib import Path

import pytest

from scripts.check_version_sync import (
    check_version_sync,
    extract_cff_version,
    extract_python_version,
    extract_zenodo_version,
)
from variantcentrifuge.version import __version__

ROOT = Path(__file__).resolve().parents[2]


def test_extract_python_version(tmp_path: Path) -> None:
    """Test extracting version from version.py."""
    version_file = tmp_path / "version.py"
    version_file.write_text('__version__ = "1.2.3"\n', encoding="utf-8")
    assert extract_python_version(version_file) == "1.2.3"


def test_extract_cff_version(tmp_path: Path) -> None:
    """Test extracting version from CITATION.cff."""
    cff_file = tmp_path / "CITATION.cff"
    cff_file.write_text("cff-version: 1.2.0\nversion: 1.2.3\ntitle: test\n", encoding="utf-8")
    assert extract_cff_version(cff_file) == "1.2.3"


def test_extract_zenodo_version(tmp_path: Path) -> None:
    """Test extracting version from .zenodo.json."""
    zenodo_file = tmp_path / ".zenodo.json"
    zenodo_file.write_text(json.dumps({"title": "test", "version": "1.2.3"}), encoding="utf-8")
    assert extract_zenodo_version(zenodo_file) == "1.2.3"


def test_current_repository_versions_are_in_sync() -> None:
    """All repository metadata files must agree on the package version."""
    versions = check_version_sync(ROOT)
    assert versions["version.py"] == __version__
    assert versions["CITATION.cff"] == __version__
    assert versions[".zenodo.json"] == __version__


def test_check_version_sync_with_matching_tag() -> None:
    """Tags with or without 'v' prefix must match the repository version."""
    check_version_sync(ROOT, tag=f"v{__version__}")
    check_version_sync(ROOT, tag=__version__)


def test_check_version_sync_rejects_tag_mismatch() -> None:
    """Mismatched tag must raise ValueError."""
    with pytest.raises(ValueError, match=r"Tag version.*does not match"):
        check_version_sync(ROOT, tag="v9.9.9")


def test_check_version_sync_detects_mismatch(tmp_path: Path) -> None:
    """Discrepancies across version files must raise ValueError."""
    pkg_dir = tmp_path / "variantcentrifuge"
    pkg_dir.mkdir()
    (pkg_dir / "version.py").write_text('__version__ = "1.0.0"\n', encoding="utf-8")
    (tmp_path / "CITATION.cff").write_text("version: 1.0.1\n", encoding="utf-8")
    (tmp_path / ".zenodo.json").write_text(json.dumps({"version": "1.0.0"}), encoding="utf-8")

    with pytest.raises(ValueError, match="Version mismatch detected"):
        check_version_sync(tmp_path)


def test_check_version_sync_detects_missing_file(tmp_path: Path) -> None:
    """Missing version file must raise FileNotFoundError."""
    pkg_dir = tmp_path / "variantcentrifuge"
    pkg_dir.mkdir()
    (pkg_dir / "version.py").write_text('__version__ = "1.0.0"\n', encoding="utf-8")

    with pytest.raises(FileNotFoundError):
        check_version_sync(tmp_path)
