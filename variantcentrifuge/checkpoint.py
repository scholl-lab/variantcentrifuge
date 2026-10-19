"""Pipeline checkpoint and resume system for VariantCentrifuge.

This module provides a robust checkpoint system that tracks pipeline execution state,
allowing resumption from the last successful step in case of interruption or failure.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
from collections.abc import Callable
from functools import wraps
from typing import Any

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
from .checkpoint_models import AtomicFileOperation, FileInfo, StepInfo
from .checkpoint_state import PipelineState

logger = logging.getLogger("variantcentrifuge")

__all__ = [
    "ANNOTATION_ALIAS_KEYS",
    "CANONICAL_ANNOTATION_KEYS",
    "LEGACY_RELEVANT_CONFIG_KEYS",
    "RELEVANT_CONFIG_KEYS",
    "AtomicFileOperation",
    "CheckpointContext",
    "FileInfo",
    "PipelineState",
    "StepInfo",
    "checkpoint",
    "format_state_summary",
    "hash_configuration",
    "hash_legacy_configuration",
    "legacy_annotation_equivalent_configs",
    "resume_configuration_hashes",
]


def checkpoint(
    step_name: str,
    input_files: str | list[str] | None = None,
    output_files: str | list[str] | None = None,
    parameters: dict[str, Any] | None = None,
) -> Callable[..., Any]:
    """Add checkpoint functionality to pipeline steps."""

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            pipeline_state: PipelineState | None = kwargs.pop("_pipeline_state", None)

            if pipeline_state is None:
                return func(*args, **kwargs)

            if pipeline_state.should_skip_step(step_name):
                logger.info(f"Skipping completed step: {step_name}")
                if output_files:
                    files = [output_files] if isinstance(output_files, str) else output_files
                    if len(files) == 1:
                        return files[0]
                    return files
                return None

            input_list: list[str] = []
            if input_files:
                in_files = input_files(*args, **kwargs) if callable(input_files) else input_files
                input_list = [in_files] if isinstance(in_files, str) else list(in_files)

            command_hash = None
            if parameters:
                param_str = json.dumps(parameters, sort_keys=True)
                command_hash = hashlib.md5(param_str.encode()).hexdigest()

            pipeline_state.start_step(step_name, command_hash, parameters)

            try:
                result = func(*args, **kwargs)
                output_list: list[str] = []
                if output_files:
                    if callable(output_files):
                        out_files = output_files(*args, **kwargs, _result=result)
                    else:
                        out_files = output_files
                    output_list = [out_files] if isinstance(out_files, str) else list(out_files)

                pipeline_state.complete_step(step_name, input_list, output_list)
                return result

            except Exception as e:
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
        command_hash: str | None = None,
        parameters: dict[str, Any] | None = None,
    ):
        self.pipeline_state = pipeline_state
        self.step_name = step_name
        self.command_hash = command_hash
        self.parameters = parameters
        self.input_files: list[str] = []
        self.output_files: list[str] = []
        self.skip: bool = False

    def __enter__(self) -> Any:
        if self.pipeline_state is None:
            self.skip = False
            return self

        if self.pipeline_state.should_skip_step(self.step_name):
            self.skip = True
            logger.info(f"Skipping completed step: {self.step_name}")
        else:
            self.skip = False
            self.pipeline_state.start_step(self.step_name, self.command_hash, self.parameters)
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        if self.skip or self.pipeline_state is None:
            return None

        if exc_type is None:
            self.pipeline_state.complete_step(self.step_name, self.input_files, self.output_files)
        else:
            self.pipeline_state.fail_step(self.step_name, str(exc_val))

        return None

    def add_input_file(self, filepath: str) -> None:
        if filepath and os.path.exists(filepath):
            self.input_files.append(filepath)

    def add_output_file(self, filepath: str) -> None:
        if filepath and os.path.exists(filepath):
            self.output_files.append(filepath)
