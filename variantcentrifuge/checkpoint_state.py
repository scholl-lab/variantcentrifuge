"""Pipeline state management and persistence for checkpoint/resume functionality."""

from __future__ import annotations

import contextlib
import json
import logging
import os
import threading
import time
import uuid
from typing import Any, ClassVar

from .checkpoint_helpers import (
    ANNOTATION_ALIAS_KEYS,
    CANONICAL_ANNOTATION_KEYS,
    LEGACY_RELEVANT_CONFIG_KEYS,
    RELEVANT_CONFIG_KEYS,
    format_state_summary,
    hash_configuration,
    hash_legacy_configuration,
    legacy_annotation_equivalent_configs,
    resume_configuration_hashes,
)
from .checkpoint_models import FileInfo, StepInfo

logger = logging.getLogger("variantcentrifuge")


class PipelineState:
    """Manages pipeline execution state for checkpoint/resume functionality."""

    STATE_FILE_NAME = ".variantcentrifuge_state.json"
    STATE_VERSION = "1.0"
    _RELEVANT_CONFIG_KEYS: ClassVar[tuple[str, ...]] = RELEVANT_CONFIG_KEYS
    _CANONICAL_ANNOTATION_KEYS: ClassVar[frozenset[str]] = CANONICAL_ANNOTATION_KEYS
    _LEGACY_RELEVANT_CONFIG_KEYS: ClassVar[tuple[str, ...]] = LEGACY_RELEVANT_CONFIG_KEYS
    _ANNOTATION_ALIAS_KEYS: ClassVar[frozenset[str]] = ANNOTATION_ALIAS_KEYS

    def __init__(self, output_dir: str, enable_checksum: bool = False):
        self.output_dir = output_dir
        self.state_file = os.path.join(output_dir, self.STATE_FILE_NAME)
        self.enable_checksum = enable_checksum
        self._state_lock = threading.Lock()

        self.state: dict[str, Any] = {
            "version": self.STATE_VERSION,
            "pipeline_version": None,
            "start_time": None,
            "last_update": None,
            "configuration_hash": None,
            "steps": {},
            "metadata": {},
        }

        self.current_step: str | None = None
        self._loaded_from_file = False

    def initialize(self, configuration: dict[str, Any], pipeline_version: str) -> None:
        """Initialize a new pipeline run."""
        self.state["pipeline_version"] = pipeline_version
        self.state["start_time"] = time.time()
        self.state["last_update"] = time.time()
        self.state["configuration_hash"] = self._hash_configuration(configuration)
        self.state["metadata"]["command_line"] = configuration.get("command_line", "")
        self.save()

    def load(self) -> bool:
        """Load existing state from file."""
        if not os.path.exists(self.state_file):
            return False

        try:
            with open(self.state_file) as f:
                loaded_state = json.load(f)

            if loaded_state.get("version") != self.STATE_VERSION:
                logger.warning(
                    f"State file version mismatch: "
                    f"{loaded_state.get('version')} != {self.STATE_VERSION}"
                )
                return False

            self.state = {
                "version": loaded_state.get("version", self.STATE_VERSION),
                "pipeline_version": loaded_state.get("pipeline_version"),
                "start_time": loaded_state.get("start_time"),
                "last_update": loaded_state.get("last_update"),
                "configuration_hash": loaded_state.get("configuration_hash"),
                "steps": {
                    name: StepInfo.from_dict(info)
                    for name, info in loaded_state.get("steps", {}).items()
                },
                "metadata": loaded_state.get("metadata", {}),
            }
            self._loaded_from_file = True
            logger.info(f"Loaded pipeline state from {self.state_file}")
            return True

        except Exception as e:
            logger.error(f"Failed to load state file: {e}")
            return False

    def save(self) -> None:
        """Save current state to file atomically (thread-safe)."""
        with self._state_lock:
            self.state["last_update"] = time.time()

            save_state = self.state.copy()
            save_state["steps"] = {
                name: step.to_dict() for name, step in self.state["steps"].items()
            }

            temp_file = f"{self.state_file}.tmp.{uuid.uuid4().hex[:8]}"
            try:
                with open(temp_file, "w") as f:
                    json.dump(save_state, f, indent=2)
                    f.flush()
                    os.fsync(f.fileno())
                os.replace(temp_file, self.state_file)
                logger.debug(f"Saved pipeline state to {self.state_file}")
            except Exception as e:
                if os.path.exists(temp_file):
                    with contextlib.suppress(Exception):
                        os.remove(temp_file)
                logger.warning(f"Failed to save checkpoint state: {e}")
                raise

    def can_resume(self, configuration: dict[str, Any], pipeline_version: str) -> bool:
        """Check if pipeline can be resumed with given configuration."""
        if not self._loaded_from_file:
            return False

        if self.state.get("pipeline_version") != pipeline_version:
            logger.warning("Pipeline version mismatch, cannot resume")
            return False

        current_hashes = self._resume_configuration_hashes(configuration)
        if self.state.get("configuration_hash") not in current_hashes:
            logger.warning("Configuration has changed, cannot resume")
            return False

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

    def get_resume_point(self) -> str | None:
        """Get the name of the last successfully completed step."""
        completed_steps = [
            (step_info.end_time or 0, name)
            for name, step_info in self.state["steps"].items()
            if step_info.status == "completed"
        ]
        if not completed_steps:
            return None

        completed_steps.sort(reverse=True)
        return str(completed_steps[0][1])

    def should_skip_step(self, step_name: str) -> bool:
        """Check if a step should be skipped (already completed with valid outputs)."""
        if not self._loaded_from_file:
            return False

        step_info = self.state["steps"].get(step_name)
        if step_info is None:
            return False

        if step_info.status == "completed":
            if step_info.output_files:
                for file_info in step_info.output_files:
                    if not file_info.validate(calculate_checksum=self.enable_checksum):
                        logger.warning(
                            f"Output file validation failed for '{file_info.path}' "
                            f"in step '{step_name}', marking for re-execution"
                        )
                        step_info.status = "pending"
                        step_info.end_time = None
                        self.save()
                        return False
            return True

        if step_info.status in ("running", "finalizing"):
            logger.warning(
                f"Found stale {step_info.status} stage '{step_name}', marking for re-execution"
            )
            if self.enable_checksum and step_info.output_files:
                logger.info(f"Attempting checksum-based recovery for stage '{step_name}'")
                all_valid = all(
                    f_info.validate(calculate_checksum=True) for f_info in step_info.output_files
                )
                if all_valid:
                    step_info.status = "completed"
                    base_time = step_info.start_time
                    step_info.end_time = base_time + 1.0 if base_time else time.time()
                    self.save()
                    return True

            step_info.status = "failed"
            step_info.error = "Pipeline was interrupted - stage will be re-executed for safety"
            step_info.end_time = time.time()
            self.save()

        return False

    def start_step(
        self,
        step_name: str,
        command_hash: str | None = None,
        parameters: dict[str, Any] | None = None,
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
        input_files: list[str] | None = None,
        output_files: list[str] | None = None,
    ) -> None:
        """Mark a step as completed."""
        with self._state_lock:
            if step_name not in self.state["steps"]:
                logger.warning(f"Completing unstarted step: {step_name}")
                step_info = StepInfo(
                    name=step_name,
                    status="running",
                    start_time=time.time(),
                )
                self.state["steps"][step_name] = step_info
            else:
                step_info = self.state["steps"][step_name]

            if input_files:
                for filepath in input_files:
                    if os.path.exists(filepath):
                        step_info.input_files.append(
                            FileInfo.from_file(filepath, self.enable_checksum)
                        )

            if output_files:
                for filepath in output_files:
                    if os.path.exists(filepath):
                        step_info.output_files.append(
                            FileInfo.from_file(filepath, self.enable_checksum)
                        )

            step_info.status = "completed"
            step_info.end_time = time.time()
            self.current_step = None

        self.save()
        logger.debug(f"Step '{step_name}' marked as completed")

    def finalize_step(self, step_name: str) -> None:
        """Mark a step as in finalizing phase."""
        with self._state_lock:
            if step_name in self.state["steps"]:
                self.state["steps"][step_name].status = "finalizing"
        self.save()

    def fail_step(self, step_name: str, error: str) -> None:
        """Mark a step as failed."""
        with self._state_lock:
            if step_name not in self.state["steps"]:
                step_info = StepInfo(
                    name=step_name,
                    status="failed",
                    start_time=time.time(),
                )
                self.state["steps"][step_name] = step_info
            else:
                step_info = self.state["steps"][step_name]

            step_info.status = "failed"
            step_info.end_time = time.time()
            step_info.error = error

        self.save()
        logger.error(f"Step '{step_name}' failed: {error}")

    def clear_step_completion(self, step_name: str, save: bool = True) -> None:
        """Clear completion status for a step to force re-execution."""
        if step_name not in self.state["steps"]:
            logger.debug(f"Step '{step_name}' not in checkpoint state, nothing to clear")
            return

        step_info = self.state["steps"][step_name]
        if step_info.status == "completed":
            step_info.status = "pending"
            step_info.end_time = None
            step_info.error = None
            logger.debug(f"Cleared completion status for step '{step_name}' (restart mode)")
            if save:
                self.save()

    def get_summary(self) -> str:
        """Get a human-readable summary of the pipeline state."""
        return format_state_summary(self.state, self.state_file)

    def get_completed_stages(self) -> list[tuple]:
        """Get all completed stages sorted by completion time."""
        completed_stages = [
            (name, info) for name, info in self.state["steps"].items() if info.status == "completed"
        ]
        completed_stages.sort(key=lambda x: x[1].end_time or 0)
        return completed_stages

    def get_available_resume_points(self) -> list[str]:
        """Get stage names that can be used as resume points."""
        return [stage_name for stage_name, _ in self.get_completed_stages()]

    def validate_resume_from_stage(self, stage_name: str, available_stages: list[str]) -> tuple:
        """Validate that resuming from specific stage is possible."""
        if not self._loaded_from_file:
            return False, "No checkpoint file loaded"
        if stage_name not in available_stages:
            return False, f"Stage '{stage_name}' is not available in current configuration"
        return True, ""

    def get_stages_to_execute(self, resume_from: str, all_stages: list[str]) -> list[str]:
        """Get ordered list of stages to execute when resuming from specific stage."""
        if resume_from not in all_stages:
            return []
        try:
            resume_index = all_stages.index(resume_from)
            return all_stages[resume_index:]
        except ValueError:
            return []

    def get_detailed_status(self) -> dict[str, Any]:
        """Get detailed status information for enhanced display."""
        if not self._loaded_from_file:
            return {"has_checkpoint": False, "message": "No checkpoint file found"}

        completed_stages = self.get_completed_stages()
        total_runtime = sum(
            info.duration for _, info in completed_stages if info.duration is not None
        )
        last_stage = completed_stages[-1] if completed_stages else None

        resume_suggestions = []
        if completed_stages and last_stage:
            resume_suggestions.append(
                {
                    "stage": last_stage[0],
                    "reason": "Continue from last completed stage",
                    "command": f"--resume-from {last_stage[0]}",
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

    def can_resume_from_stage(self, stage_name: str, available_stages: list[str]) -> bool:
        """Check if resuming from a specific stage is possible."""
        is_valid, _ = self.validate_resume_from_stage(stage_name, available_stages)
        return bool(is_valid)

    def _hash_configuration(self, config: dict[str, Any]) -> str:
        return hash_configuration(config)

    def _resume_configuration_hashes(self, config: dict[str, Any]) -> set[str]:
        return resume_configuration_hashes(config)

    def _hash_legacy_configuration(self, config: dict[str, Any]) -> str:
        return hash_legacy_configuration(config)

    def _legacy_annotation_equivalent_configs(self, config: dict[str, Any]) -> list[dict[str, Any]]:
        return legacy_annotation_equivalent_configs(config)

    def cleanup_stale_stages(self) -> None:
        """Clean up stages that were left in 'running' state from previous pipeline runs."""
        if not self._loaded_from_file:
            return

        for name, step_info in self.state["steps"].items():
            if step_info.status in ("running", "finalizing"):
                logger.warning(
                    f"Resetting stale {step_info.status} stage '{name}' to failed for re-execution"
                )
                step_info.status = "failed"
                step_info.error = "Pipeline interrupted during execution"
                step_info.end_time = time.time()
        self.save()

    def __getstate__(self) -> dict[str, Any]:
        state = self.__dict__.copy()
        state.pop("_state_lock", None)
        return state

    def __setstate__(self, state: dict[str, Any]) -> None:
        self.__dict__.update(state)
        self._state_lock = threading.Lock()
