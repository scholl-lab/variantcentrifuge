"""Pipeline checkpoint and resume system for VariantCentrifuge.

This module provides a robust checkpoint system that tracks pipeline execution state,
allowing resumption from the last successful step in case of interruption or failure.
"""

import json
import hashlib
import os
import time
import logging
from typing import Dict, Any, Optional, List, Callable, Union
from functools import wraps
from dataclasses import dataclass, field, asdict
from datetime import datetime
import gzip

logger = logging.getLogger("variantcentrifuge")


@dataclass
class FileInfo:
    """Information about a file for validation."""

    path: str
    size: int
    mtime: float
    checksum: Optional[str] = None  # Optional for performance

    @classmethod
    def from_file(cls, filepath: str, calculate_checksum: bool = False) -> "FileInfo":
        """Create FileInfo from an existing file."""
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        stat = os.stat(filepath)
        checksum = None

        if calculate_checksum:
            checksum = cls._calculate_checksum(filepath)

        return cls(path=filepath, size=stat.st_size, mtime=stat.st_mtime, checksum=checksum)

    @staticmethod
    def _calculate_checksum(filepath: str, chunk_size: int = 8192) -> str:
        """Calculate SHA256 checksum of a file."""
        sha256 = hashlib.sha256()

        # Handle both gzipped and regular files
        if filepath.endswith(".gz"):
            open_func = gzip.open
            mode = "rb"
        else:
            open_func = open
            mode = "rb"

        with open_func(filepath, mode) as f:
            while chunk := f.read(chunk_size):
                sha256.update(chunk)

        return sha256.hexdigest()

    def validate(self, calculate_checksum: bool = False) -> bool:
        """Validate that the file matches the stored info."""
        if not os.path.exists(self.path):
            return False

        stat = os.stat(self.path)
        if stat.st_size != self.size:
            return False

        # Allow some tolerance for mtime (filesystem precision)
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
    status: str  # "pending", "running", "completed", "failed"
    start_time: Optional[float] = None
    end_time: Optional[float] = None
    command_hash: Optional[str] = None
    input_files: List[FileInfo] = field(default_factory=list)
    output_files: List[FileInfo] = field(default_factory=list)
    parameters: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None

    @property
    def duration(self) -> Optional[float]:
        """Calculate step duration in seconds."""
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        return None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        data = asdict(self)
        # Convert FileInfo objects to dicts
        data["input_files"] = [asdict(f) for f in self.input_files]
        data["output_files"] = [asdict(f) for f in self.output_files]
        return data

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "StepInfo":
        """Create StepInfo from dictionary."""
        # Convert file dicts back to FileInfo objects
        data["input_files"] = [FileInfo(**f) for f in data.get("input_files", [])]
        data["output_files"] = [FileInfo(**f) for f in data.get("output_files", [])]
        return cls(**data)


class PipelineState:
    """Manages pipeline execution state for checkpoint/resume functionality."""

    STATE_FILE_NAME = ".variantcentrifuge_state.json"
    STATE_VERSION = "1.0"

    def __init__(self, output_dir: str, enable_checksum: bool = False):
        """Initialize pipeline state manager.

        Parameters
        ----------
        output_dir : str
            Directory where pipeline outputs and state file are stored
        enable_checksum : bool
            Whether to calculate file checksums (slower but more reliable)
        """
        self.output_dir = output_dir
        self.state_file = os.path.join(output_dir, self.STATE_FILE_NAME)
        self.enable_checksum = enable_checksum

        self.state = {
            "version": self.STATE_VERSION,
            "pipeline_version": None,  # Set by pipeline
            "start_time": None,
            "last_update": None,
            "configuration_hash": None,
            "steps": {},
            "metadata": {},
        }

        self.current_step: Optional[str] = None
        self._loaded_from_file = False

    def initialize(self, configuration: Dict[str, Any], pipeline_version: str) -> None:
        """Initialize a new pipeline run."""
        self.state["pipeline_version"] = pipeline_version
        self.state["start_time"] = time.time()
        self.state["last_update"] = time.time()
        self.state["configuration_hash"] = self._hash_configuration(configuration)
        self.state["metadata"]["command_line"] = configuration.get("command_line", "")
        self.save()

    def load(self) -> bool:
        """Load existing state from file.

        Returns
        -------
        bool
            True if state was loaded successfully, False otherwise
        """
        if not os.path.exists(self.state_file):
            return False

        try:
            with open(self.state_file, "r") as f:
                loaded_state = json.load(f)

            # Validate version compatibility
            if loaded_state.get("version") != self.STATE_VERSION:
                logger.warning(
                    f"State file version mismatch: "
                    f"{loaded_state.get('version')} != {self.STATE_VERSION}"
                )
                return False

            # Convert step dicts back to StepInfo objects
            loaded_state["steps"] = {
                name: StepInfo.from_dict(info)
                for name, info in loaded_state.get("steps", {}).items()
            }

            self.state = loaded_state
            self._loaded_from_file = True
            logger.info(f"Loaded pipeline state from {self.state_file}")
            return True

        except Exception as e:
            logger.error(f"Failed to load state file: {e}")
            return False

    def save(self) -> None:
        """Save current state to file."""
        self.state["last_update"] = time.time()

        # Convert to JSON-serializable format
        save_state = self.state.copy()
        save_state["steps"] = {name: step.to_dict() for name, step in self.state["steps"].items()}

        # Write atomically
        temp_file = self.state_file + ".tmp"
        with open(temp_file, "w") as f:
            json.dump(save_state, f, indent=2)
        os.replace(temp_file, self.state_file)

        logger.debug(f"Saved pipeline state to {self.state_file}")

    def can_resume(self, configuration: Dict[str, Any], pipeline_version: str) -> bool:
        """Check if pipeline can be resumed with given configuration.

        Parameters
        ----------
        configuration : Dict[str, Any]
            Current pipeline configuration
        pipeline_version : str
            Current pipeline version

        Returns
        -------
        bool
            True if resume is possible, False otherwise
        """
        if not self._loaded_from_file:
            return False

        # Check pipeline version compatibility
        if self.state.get("pipeline_version") != pipeline_version:
            logger.warning("Pipeline version mismatch, cannot resume")
            return False

        # Check configuration hash
        current_hash = self._hash_configuration(configuration)
        if self.state.get("configuration_hash") != current_hash:
            logger.warning("Configuration has changed, cannot resume")
            return False

        # Special case: If pipeline already completed, only validate final output
        if self.should_skip_step("final_output"):
            logger.info("Pipeline was already completed, only validating final output.")
            return True  # Let the pipeline.py handle checking for the final output file

        # Validate that output files still exist and match
        for step_name, step_info in self.state["steps"].items():
            if step_info.status == "completed":
                for file_info in step_info.output_files:
                    if not file_info.validate(self.enable_checksum):
                        logger.warning(
                            f"Output file validation failed for step '{step_name}': "
                            f"{file_info.path}"
                        )
                        return False

        return True

    def get_resume_point(self) -> Optional[str]:
        """Get the name of the last successfully completed step.

        Returns
        -------
        Optional[str]
            Name of last completed step, or None if no steps completed
        """
        completed_steps = []

        for name, info in self.state["steps"].items():
            if info.status == "completed":
                completed_steps.append((info.end_time or 0, name))

        if not completed_steps:
            return None

        # Return the most recently completed step
        completed_steps.sort(reverse=True)
        return completed_steps[0][1]

    def should_skip_step(self, step_name: str) -> bool:
        """Check if a step should be skipped (already completed).

        Parameters
        ----------
        step_name : str
            Name of the step to check

        Returns
        -------
        bool
            True if step should be skipped, False otherwise
        """
        if not self._loaded_from_file:
            return False

        step_info = self.state["steps"].get(step_name)
        return step_info is not None and step_info.status == "completed"

    def start_step(
        self,
        step_name: str,
        command_hash: Optional[str] = None,
        parameters: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Mark a step as started."""
        self.current_step = step_name

        step_info = StepInfo(
            name=step_name,
            status="running",
            start_time=time.time(),
            command_hash=command_hash,
            parameters=parameters or {},
        )

        self.state["steps"][step_name] = step_info
        self.save()

    def complete_step(
        self,
        step_name: str,
        input_files: Optional[List[str]] = None,
        output_files: Optional[List[str]] = None,
    ) -> None:
        """Mark a step as completed."""
        if step_name not in self.state["steps"]:
            logger.warning(f"Completing unstarted step: {step_name}")
            self.start_step(step_name)

        step_info = self.state["steps"][step_name]
        step_info.status = "completed"
        step_info.end_time = time.time()

        # Record input files
        if input_files:
            for filepath in input_files:
                if os.path.exists(filepath):
                    step_info.input_files.append(FileInfo.from_file(filepath, self.enable_checksum))

        # Record output files
        if output_files:
            for filepath in output_files:
                if os.path.exists(filepath):
                    step_info.output_files.append(
                        FileInfo.from_file(filepath, self.enable_checksum)
                    )

        self.save()
        logger.info(f"Completed step '{step_name}' in {step_info.duration:.2f}s")

    def fail_step(self, step_name: str, error: str) -> None:
        """Mark a step as failed."""
        if step_name not in self.state["steps"]:
            logger.warning(f"Failing unstarted step: {step_name}")
            self.start_step(step_name)

        step_info = self.state["steps"][step_name]
        step_info.status = "failed"
        step_info.end_time = time.time()
        step_info.error = error

        self.save()
        logger.error(f"Step '{step_name}' failed: {error}")

    def get_summary(self) -> str:
        """Get a human-readable summary of the pipeline state."""
        lines = ["Pipeline State Summary:"]
        lines.append(f"  State file: {self.state_file}")
        lines.append(f"  Pipeline version: {self.state.get('pipeline_version', 'unknown')}")

        if self.state.get("start_time"):
            start_time = datetime.fromtimestamp(self.state["start_time"])
            lines.append(f"  Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

        lines.append("\nSteps:")
        for name, info in self.state["steps"].items():
            status_symbol = {"completed": "✓", "failed": "✗", "running": "→", "pending": "·"}.get(
                info.status, "?"
            )

            duration_str = ""
            if info.duration:
                duration_str = f" ({info.duration:.1f}s)"

            lines.append(f"  {status_symbol} {name}{duration_str}")

            if info.status == "failed" and info.error:
                lines.append(f"    Error: {info.error}")

        return "\n".join(lines)

    def _hash_configuration(self, config: Dict[str, Any]) -> str:
        """Create a hash of the configuration for change detection."""
        # Select relevant configuration keys that affect pipeline behavior
        relevant_keys = [
            "gene_name",
            "gene_file",
            "vcf_file",
            "output_file",
            "preset",
            "filters",
            "fields",
            "threads",
            "no_replacement",
            "late_filtering",
            "chunk_size",
            "samples_file",
            "bcftools_prefilter",
            "final_filter",
            "scoring_config_path",
            "ped",
            "inheritance_mode",
            "annotate_bed",
            "annotate_gene_list",
            "annotate_json_genes",
        ]

        relevant_config = {k: v for k, v in config.items() if k in relevant_keys and v is not None}

        # Sort for consistent hashing
        config_str = json.dumps(relevant_config, sort_keys=True)
        return hashlib.sha256(config_str.encode()).hexdigest()


def checkpoint(
    step_name: str,
    input_files: Optional[Union[str, List[str]]] = None,
    output_files: Optional[Union[str, List[str]]] = None,
    parameters: Optional[Dict[str, Any]] = None,
):
    """Add checkpoint functionality to pipeline steps.

    Parameters
    ----------
    step_name : str
        Unique name for this pipeline step
    input_files : Optional[Union[str, List[str]]]
        Input file(s) to track
    output_files : Optional[Union[str, List[str]]]
        Output file(s) to track
    parameters : Optional[Dict[str, Any]]
        Additional parameters to record
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Get pipeline state from kwargs or create a dummy one
            pipeline_state: Optional[PipelineState] = kwargs.pop("_pipeline_state", None)

            if pipeline_state is None:
                # No checkpoint system active, run normally
                return func(*args, **kwargs)

            # Check if step should be skipped
            if pipeline_state.should_skip_step(step_name):
                logger.info(f"Skipping completed step: {step_name}")

                # Return expected output files if specified
                if output_files:
                    files = [output_files] if isinstance(output_files, str) else output_files
                    if len(files) == 1:
                        return files[0]
                    return files
                return None

            # Prepare file lists
            input_list = []
            if input_files:
                if callable(input_files):
                    # Dynamic input files
                    files = input_files(*args, **kwargs)
                else:
                    files = input_files
                input_list = [files] if isinstance(files, str) else files

            # Start the step
            command_hash = None
            if parameters:
                param_str = json.dumps(parameters, sort_keys=True)
                command_hash = hashlib.md5(param_str.encode()).hexdigest()

            pipeline_state.start_step(step_name, command_hash, parameters)

            try:
                # Execute the function
                result = func(*args, **kwargs)

                # Prepare output file list
                output_list = []
                if output_files:
                    if callable(output_files):
                        # Dynamic output files
                        files = output_files(*args, **kwargs, _result=result)
                    else:
                        files = output_files
                    output_list = [files] if isinstance(files, str) else files

                # Mark step as completed
                pipeline_state.complete_step(step_name, input_list, output_list)

                return result

            except Exception as e:
                # Mark step as failed
                pipeline_state.fail_step(step_name, str(e))
                raise

        return wrapper

    return decorator


class CheckpointContext:
    """Context manager for checkpoint operations within a code block."""

    def __init__(
        self,
        pipeline_state: PipelineState,
        step_name: str,
        command_hash: Optional[str] = None,
        parameters: Optional[Dict[str, Any]] = None,
    ):
        self.pipeline_state = pipeline_state
        self.step_name = step_name
        self.command_hash = command_hash
        self.parameters = parameters
        self.input_files = []
        self.output_files = []

    def __enter__(self):
        """Enter the checkpoint context."""
        if self.pipeline_state is None:
            # No pipeline state - just run normally
            self.skip = False
            return self

        if self.pipeline_state.should_skip_step(self.step_name):
            self.skip = True
            logger.info(f"Skipping completed step: {self.step_name}")
        else:
            self.skip = False
            self.pipeline_state.start_step(self.step_name, self.command_hash, self.parameters)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit the checkpoint context and save state if successful."""
        if self.skip or self.pipeline_state is None:
            return False

        if exc_type is None:
            # Success
            self.pipeline_state.complete_step(self.step_name, self.input_files, self.output_files)
        else:
            # Failure
            self.pipeline_state.fail_step(self.step_name, str(exc_val))

        return False  # Don't suppress exceptions

    def add_input_file(self, filepath: str) -> None:
        """Add an input file to track."""
        if filepath and os.path.exists(filepath):
            self.input_files.append(filepath)

    def add_output_file(self, filepath: str) -> None:
        """Add an output file to track."""
        if filepath and os.path.exists(filepath):
            self.output_files.append(filepath)
