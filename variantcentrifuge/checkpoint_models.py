"""Data models and atomic file utilities for pipeline checkpoints."""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import logging
import os
import uuid
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger("variantcentrifuge")


class AtomicFileOperation:
    """Context manager for atomic file operations with automatic temp file cleanup.

    Ensures that output files are only visible in their final location after
    being completely written, using unique sibling temporary files and atomic os.replace.
    """

    def __init__(self, final_path: str | Path, suffix: str | None = None):
        self.final_path = Path(final_path)
        if suffix is not None:
            self.temp_path = self.final_path.with_name(f"{self.final_path.name}{suffix}")
        else:
            full_suffix = "".join(self.final_path.suffixes)
            stem = (
                self.final_path.name[: -len(full_suffix)] if full_suffix else self.final_path.name
            )
            unique_tag = f".tmp.{os.getpid()}.{uuid.uuid4().hex[:8]}"
            self.temp_path = self.final_path.with_name(f"{stem}{unique_tag}{full_suffix}")
        self.success = False

    def __enter__(self) -> Path:
        self.final_path.parent.mkdir(parents=True, exist_ok=True)
        with contextlib.suppress(FileNotFoundError, OSError):
            self.temp_path.unlink(missing_ok=True)
        return self.temp_path

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        if exc_type is None:
            try:
                os.replace(self.temp_path, self.final_path)
                self.success = True
                logger.debug(f"Atomically published file: {self.final_path}")
            except Exception as err:
                with contextlib.suppress(FileNotFoundError, OSError):
                    self.temp_path.unlink(missing_ok=True)
                logger.error(f"Failed to publish {self.temp_path} to {self.final_path}: {err}")
                raise
        else:
            try:
                self.temp_path.unlink(missing_ok=True)
                logger.debug(f"Cleaned up temp file after error: {self.temp_path}")
            except (FileNotFoundError, OSError):
                pass
            except Exception as cleanup_error:
                logger.warning(f"Failed to clean up temp file {self.temp_path}: {cleanup_error}")
        return None


@dataclass
class FileInfo:
    """Information about a file for validation."""

    path: str
    size: int
    mtime: float
    checksum: str | None = None

    @classmethod
    def from_file(cls, filepath: str, calculate_checksum: bool = False) -> FileInfo:
        """Create FileInfo from an existing file."""
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        stat = os.stat(filepath)
        checksum = cls._calculate_checksum(filepath) if calculate_checksum else None
        return cls(path=filepath, size=stat.st_size, mtime=stat.st_mtime, checksum=checksum)

    @staticmethod
    def _calculate_checksum(filepath: str, chunk_size: int = 8192) -> str:
        """Calculate SHA256 checksum of a file."""
        sha256 = hashlib.sha256()
        open_func: Callable[..., Any] = gzip.open if filepath.endswith(".gz") else open

        with open_func(filepath, "rb") as f:
            while chunk := f.read(chunk_size):
                sha256.update(chunk)  # type: ignore[arg-type]

        return sha256.hexdigest()

    def validate(self, calculate_checksum: bool = False) -> bool:
        """Validate that the file matches the stored info."""
        if not os.path.exists(self.path):
            return False

        stat = os.stat(self.path)
        if stat.st_size != self.size:
            return False

        if abs(stat.st_mtime - self.mtime) > 1.0:
            return False

        if calculate_checksum and self.checksum:
            current_checksum = self._calculate_checksum(self.path)
            if current_checksum != self.checksum:
                return False

        return True


@dataclass
class StepInfo:
    """Information about a pipeline step."""

    name: str
    status: str  # "pending", "running", "finalizing", "completed", "failed"
    start_time: float | None = None
    end_time: float | None = None
    command_hash: str | None = None
    input_files: list[FileInfo] = field(default_factory=list)
    output_files: list[FileInfo] = field(default_factory=list)
    parameters: dict[str, Any] = field(default_factory=dict)
    error: str | None = None

    @property
    def duration(self) -> float | None:
        """Calculate step duration in seconds."""
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        return None

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        data = asdict(self)
        data["input_files"] = [asdict(f) for f in self.input_files]
        data["output_files"] = [asdict(f) for f in self.output_files]
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> StepInfo:
        """Create StepInfo from dictionary."""
        data["input_files"] = [FileInfo(**f) for f in data.get("input_files", [])]
        data["output_files"] = [FileInfo(**f) for f in data.get("output_files", [])]
        return cls(**data)
