"""Pipeline checkpoint and resume system for VariantCentrifuge.

This module provides a robust checkpoint system that tracks pipeline execution state,
allowing resumption from the last successful step in case of interruption or failure.
"""

import gzip
import hashlib
import json
import logging
import os
import threading
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from functools import wraps
from typing import Any, Callable, Dict, List, Optional, Union

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
        self._state_lock = threading.Lock()  # Thread-safe state operations

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
        """Save current state to file (thread-safe)."""
        with self._state_lock:
            self.state["last_update"] = time.time()

            # Convert to JSON-serializable format
            save_state = self.state.copy()
            save_state["steps"] = {
                name: step.to_dict() for name, step in self.state["steps"].items()
            }

            # Write atomically with unique temp file name to avoid conflicts
            import uuid

            temp_file = f"{self.state_file}.tmp.{uuid.uuid4().hex[:8]}"
            try:
                with open(temp_file, "w") as f:
                    json.dump(save_state, f, indent=2)
                os.replace(temp_file, self.state_file)
                logger.debug(f"Saved pipeline state to {self.state_file}")
            except Exception as e:
                # Clean up temp file if it exists
                if os.path.exists(temp_file):
                    try:
                        os.remove(temp_file)
                    except Exception:
                        pass
                logger.warning(f"Failed to save checkpoint state: {e}")
                raise

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
        if step_info is None:
            return False

        # If already completed, skip
        if step_info.status == "completed":
            return True

        # Handle stale "running" stages - these were started but pipeline was interrupted
        if step_info.status == "running":
            # Check if output files exist to determine if stage actually completed
            if step_info.output_files:
                # Check if all output files exist
                all_outputs_exist = all(
                    os.path.exists(file_info.path) for file_info in step_info.output_files
                )
                if all_outputs_exist:
                    # Stage actually completed, mark it as such
                    step_info.status = "completed"
                    step_info.end_time = step_info.start_time + 1.0  # Estimate completion time
                    logger.info(
                        f"Recovered completed stage '{step_name}' from stale 'running' state"
                    )
                    self.save()
                    return True

            # Stage has no output files - check if it might have actually completed
            # For stages that work in memory, we can't verify completion via files
            if not step_info.output_files:
                # If this is the very last stage that was running, assume it was interrupted
                # Otherwise, assume it completed since subsequent stages wouldn't have started
                all_step_times = [
                    (s.start_time or 0, name)
                    for name, s in self.state["steps"].items()
                    if s.start_time is not None
                ]
                all_step_times.sort()

                if all_step_times and all_step_times[-1][1] == step_name:
                    # This was the last stage started - likely interrupted
                    logger.warning(
                        f"Found stale running stage '{step_name}' (last started), will re-execute"
                    )
                    step_info.status = "failed"
                    step_info.error = "Pipeline was interrupted"
                    step_info.end_time = time.time()
                    self.save()
                else:
                    # Not the last stage - likely completed, just not marked properly
                    logger.info(
                        f"Recovered completed stage '{step_name}' from stale 'running' state "
                        f"(no output files)"
                    )
                    step_info.status = "completed"
                    step_info.end_time = step_info.start_time + 1.0
                    self.save()
                    return True
            else:
                # Stage has output files but they don't exist - definitely incomplete
                logger.warning(
                    f"Found stale running stage '{step_name}', output files missing, "
                    f"will re-execute"
                )
                step_info.status = "failed"
                step_info.error = "Pipeline was interrupted"
                step_info.end_time = time.time()
                self.save()

        return False

    def start_step(
        self,
        step_name: str,
        command_hash: Optional[str] = None,
        parameters: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Mark a step as started."""
        with self._state_lock:
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
        with self._state_lock:
            if step_name not in self.state["steps"]:
                logger.warning(f"Completing unstarted step: {step_name}")
                # Start step without saving yet
                step_info = StepInfo(
                    name=step_name,
                    status="running",
                    start_time=time.time(),
                    command_hash=None,
                    parameters={},
                )
                self.state["steps"][step_name] = step_info

            step_info = self.state["steps"][step_name]

            # If already completed, don't overwrite unless we have new file information
            if step_info.status == "completed":
                if output_files and not step_info.output_files:
                    # Add missing file information to completed step
                    for filepath in output_files:
                        if os.path.exists(filepath):
                            step_info.output_files.append(
                                FileInfo.from_file(filepath, self.enable_checksum)
                            )
                    logger.debug(
                        f"Added {len(output_files)} output files to completed step '{step_name}'"
                    )
                    self.save()  # Save the file updates
                return

            # Mark as completed
            step_info.status = "completed"
            step_info.end_time = time.time()

            # Record input files
            if input_files:
                for filepath in input_files:
                    if os.path.exists(filepath):
                        step_info.input_files.append(
                            FileInfo.from_file(filepath, self.enable_checksum)
                        )

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

    def clear_step_completion(self, step_name: str) -> None:
        """Clear completion status for a step to force re-execution.

        This method is used for restart functionality where you want to
        re-run a previously completed step.

        Parameters
        ----------
        step_name : str
            Name of the step to clear completion status for
        """
        if step_name not in self.state["steps"]:
            logger.debug(f"Step '{step_name}' not in checkpoint state, nothing to clear")
            return

        step_info = self.state["steps"][step_name]
        if step_info.status == "completed":
            # Reset to a state that will force re-execution
            step_info.status = "pending"
            step_info.end_time = None
            step_info.error = None
            # Keep output files for potential validation but clear the completion status
            logger.debug(f"Cleared completion status for step '{step_name}' (restart mode)")
            self.save()

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

    def get_completed_stages(self) -> List[tuple]:
        """Get all completed stages sorted by completion time.

        Returns
        -------
        List[tuple]
            List of (stage_name, StepInfo) tuples sorted by completion time
        """
        completed_stages = []

        for name, info in self.state["steps"].items():
            if info.status == "completed":
                completed_stages.append((name, info))

        # Sort by completion time (end_time)
        completed_stages.sort(key=lambda x: x[1].end_time or 0)
        return completed_stages

    def get_available_resume_points(self) -> List[str]:
        """Get stage names that can be used as resume points.

        Returns
        -------
        List[str]
            List of stage names that are safe resume points
        """
        completed_stages = self.get_completed_stages()
        return [stage_name for stage_name, _ in completed_stages]

    def validate_resume_from_stage(self, stage_name: str, available_stages: List[str]) -> tuple:
        """Validate that resuming from specific stage is possible.

        Parameters
        ----------
        stage_name : str
            Name of the stage to resume from
        available_stages : List[str]
            List of all available stage names in current configuration

        Returns
        -------
        tuple
            (is_valid: bool, error_message: str)
        """
        if not self._loaded_from_file:
            return False, "No checkpoint file loaded"

        # Check if stage name exists in available stages
        if stage_name not in available_stages:
            return False, f"Stage '{stage_name}' is not available in current configuration"

        # Check if we have any completed stages
        completed_stages = self.get_available_resume_points()
        if not completed_stages:
            return False, "No completed stages found in checkpoint"

        # If stage is already completed, it's a valid resume point
        if stage_name in completed_stages:
            return True, ""

        # If stage is not completed, check if it makes sense to resume from it
        # This would mean starting the pipeline from this stage without its dependencies
        return (
            False,
            f"Stage '{stage_name}' was not completed in previous run. "
            f"Cannot resume from an incomplete stage.",
        )

    def get_stages_to_execute(self, resume_from: str, all_stages: List[str]) -> List[str]:
        """Get ordered list of stages to execute when resuming from specific stage.

        Parameters
        ----------
        resume_from : str
            Name of the stage to resume from
        all_stages : List[str]
            List of all stage names in execution order

        Returns
        -------
        List[str]
            Ordered list of stages to execute starting from resume_from
        """
        if resume_from not in all_stages:
            return []

        # Find the index of the resume stage
        try:
            resume_index = all_stages.index(resume_from)
            # Return all stages from resume_from onwards
            return all_stages[resume_index:]
        except ValueError:
            return []

    def get_detailed_status(self) -> Dict[str, Any]:
        """Get detailed status information for enhanced display.

        Returns
        -------
        Dict[str, Any]
            Detailed status information including suggestions
        """
        if not self._loaded_from_file:
            return {"has_checkpoint": False, "message": "No checkpoint file found"}

        completed_stages = self.get_completed_stages()

        # Calculate total runtime
        total_runtime = sum(
            info.duration for _, info in completed_stages if info.duration is not None
        )

        # Get the most recent stage
        last_stage = completed_stages[-1] if completed_stages else None

        # Generate resume suggestions
        resume_suggestions = []
        if completed_stages:
            # Suggest resuming from the last stage for quick continuation
            last_stage_name = last_stage[0] if last_stage else None
            if last_stage_name:
                resume_suggestions.append(
                    {
                        "stage": last_stage_name,
                        "reason": "Continue from last completed stage",
                        "command": f"--resume-from {last_stage_name}",
                    }
                )

            # Suggest common resume points
            stage_names = [name for name, _ in completed_stages]
            if "variant_analysis" in stage_names:
                resume_suggestions.append(
                    {
                        "stage": "variant_analysis",
                        "reason": "Re-run analysis and all reports",
                        "command": "--resume-from variant_analysis",
                    }
                )
            if "tsv_output" in stage_names:
                resume_suggestions.append(
                    {
                        "stage": "tsv_output",
                        "reason": "Re-generate only reports (Excel, HTML, IGV)",
                        "command": "--resume-from tsv_output",
                    }
                )

        return {
            "has_checkpoint": True,
            "total_stages": len(self.state["steps"]),
            "completed_stages": len(completed_stages),
            "total_runtime": total_runtime,
            "last_stage": last_stage[0] if last_stage else None,
            "last_completed_time": last_stage[1].end_time if last_stage else None,
            "pipeline_version": self.state.get("pipeline_version"),
            "start_time": self.state.get("start_time"),
            "available_resume_points": [name for name, _ in completed_stages],
            "resume_suggestions": resume_suggestions,
            "stages": [
                {
                    "name": name,
                    "status": info.status,
                    "start_time": info.start_time,
                    "end_time": info.end_time,
                    "duration": info.duration,
                    "error": info.error,
                }
                for name, info in self.state["steps"].items()
            ],
        }

    def can_resume_from_stage(self, stage_name: str, available_stages: List[str]) -> bool:
        """Check if resuming from a specific stage is possible.

        Parameters
        ----------
        stage_name : str
            Name of the stage to resume from
        available_stages : List[str]
            List of all available stage names

        Returns
        -------
        bool
            True if resume is possible, False otherwise
        """
        is_valid, _ = self.validate_resume_from_stage(stage_name, available_stages)
        return is_valid

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

    def cleanup_stale_stages(self) -> None:
        """Clean up stages that were left in 'running' state from previous pipeline runs."""
        if not self._loaded_from_file:
            return

        stale_stages = []
        completed_stages = []

        for step_name, step_info in self.state["steps"].items():
            if step_info.status == "running":
                stale_stages.append(step_name)
            elif step_info.status == "completed":
                completed_stages.append(step_name)

        if stale_stages:
            logger.info(
                f"Found {len(stale_stages)} stale running stages: {', '.join(stale_stages)}"
            )
            for step_name in stale_stages:
                # Trigger the cleanup logic in should_skip_step
                self.should_skip_step(step_name)

            # If all stages were stale and none are actually completed, warn user
            if not completed_stages and len(stale_stages) > 5:
                logger.warning(
                    "Checkpoint appears corrupted (all stages were stale). "
                    "Consider removing checkpoint to start fresh."
                )

        completed_count = len([s for s in self.state["steps"].values() if s.status == "completed"])
        logger.info(f"After cleanup: {completed_count} completed stages")

    def __getstate__(self):
        """Exclude the lock from serialization."""
        state = self.__dict__.copy()
        # Remove the unpicklable lock
        del state["_state_lock"]
        return state

    def __setstate__(self, state):
        """Recreate the lock after unpickling."""
        self.__dict__.update(state)
        # Recreate the lock
        self._state_lock = threading.Lock()


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
