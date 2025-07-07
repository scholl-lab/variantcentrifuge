"""
PipelineContext - Single source of truth for pipeline state and data.

This module provides the PipelineContext dataclass that flows through all stages,
carrying configuration, state, and data artifacts throughout the pipeline execution.
"""

import argparse
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, TYPE_CHECKING
import threading

import pandas as pd

if TYPE_CHECKING:
    from .workspace import Workspace

logger = logging.getLogger(__name__)


class PicklableLock:
    """A lock that can be pickled by creating a new lock on unpickling."""

    def __init__(self):
        self._lock = threading.RLock()

    def __getstate__(self):
        # Return empty state - we don't need to save the lock state
        return {}

    def __setstate__(self, state):
        # Create a fresh lock on unpickling
        self._lock = threading.RLock()

    def __enter__(self):
        # Ensure lock exists (in case of edge cases)
        if not hasattr(self, "_lock"):
            self._lock = threading.RLock()
        return self._lock.__enter__()

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self._lock.__exit__(exc_type, exc_val, exc_tb)

    def acquire(self, blocking=True, timeout=-1):
        if not hasattr(self, "_lock"):
            self._lock = threading.RLock()
        return self._lock.acquire(blocking, timeout)

    def release(self):
        return self._lock.release()


@dataclass
class PipelineContext:
    """Container for all pipeline state and data - the single source of truth.

    This context object flows through all pipeline stages, maintaining state
    and providing access to configuration, workspace paths, and data artifacts.

    Attributes
    ----------
    args : argparse.Namespace
        Parsed command-line arguments
    config : Dict[str, Any]
        Merged configuration from file and CLI
    workspace : Workspace
        Manages all file paths for the pipeline run
    start_time : datetime
        Pipeline execution start time
    data : Any
        Primary data artifact that flows through pipeline (file path or DataFrame)
    completed_stages : Set[str]
        Names of stages that have completed successfully
    stage_results : Dict[str, Any]
        Optional results stored by stages
    vcf_samples : List[str]
        Sample names from VCF file
    pedigree_data : Optional[Dict]
        Loaded pedigree information
    phenotype_data : Optional[Dict]
        Loaded phenotype data
    scoring_config : Optional[Dict]
        Variant scoring configuration
    annotation_configs : Dict[str, Any]
        Custom annotation configurations
    gene_bed_file : Optional[Path]
        Generated BED file for gene regions
    extracted_vcf : Optional[Path]
        Path to extracted variants VCF
    filtered_vcf : Optional[Path]
        Path to filtered VCF
    extracted_tsv : Optional[Path]
        Path to extracted fields TSV
    current_dataframe : Optional[pd.DataFrame]
        Current DataFrame for analysis stages
    statistics : Dict[str, Any]
        Calculated statistics
    gene_burden_results : Optional[pd.DataFrame]
        Gene burden analysis results
    final_output_path : Optional[Path]
        Path to final output file
    report_paths : Dict[str, Path]
        Paths to generated reports
    checkpoint_state : Optional[Any]
        Checkpoint system state for resume capability
    """

    # --- Immutable Configuration ---
    args: argparse.Namespace
    config: Dict[str, Any]
    workspace: "Workspace"  # Forward reference, will be imported
    start_time: datetime = field(default_factory=datetime.now)

    # --- Mutable State & Data Artifacts ---
    # Core data that flows through pipeline
    data: Any = None  # Can be file path or DataFrame

    # Stage tracking
    completed_stages: Set[str] = field(default_factory=set)
    stage_results: Dict[str, Any] = field(default_factory=dict)

    # Loaded configurations
    vcf_samples: List[str] = field(default_factory=list)
    pedigree_data: Optional[Dict] = None
    phenotype_data: Optional[Dict] = None
    scoring_config: Optional[Dict] = None
    annotation_configs: Dict[str, Any] = field(default_factory=dict)

    # Processing artifacts
    gene_bed_file: Optional[Path] = None
    extracted_vcf: Optional[Path] = None
    filtered_vcf: Optional[Path] = None
    extracted_tsv: Optional[Path] = None

    # Analysis results
    current_dataframe: Optional[pd.DataFrame] = None
    statistics: Dict[str, Any] = field(default_factory=dict)
    gene_burden_results: Optional[pd.DataFrame] = None

    # Output paths
    final_output_path: Optional[Path] = None
    report_paths: Dict[str, Path] = field(default_factory=dict)

    # Checkpoint state (if enabled)
    checkpoint_state: Optional[Any] = None

    # Thread safety lock for parallel stages
    _lock: PicklableLock = field(default_factory=PicklableLock, init=False, repr=False)

    def mark_complete(self, stage_name: str, result: Any = None) -> None:
        """Mark a stage as complete with optional result storage.

        Parameters
        ----------
        stage_name : str
            Name of the stage to mark complete
        result : Any, optional
            Optional result to store for the stage
        """
        with self._lock:
            self.completed_stages.add(stage_name)
            if result is not None:
                self.stage_results[stage_name] = result
            logger.debug(f"Stage '{stage_name}' marked as complete")

    def is_complete(self, stage_name: str) -> bool:
        """Check if a stage has been completed.

        Parameters
        ----------
        stage_name : str
            Name of the stage to check

        Returns
        -------
        bool
            True if the stage has completed
        """
        with self._lock:
            return stage_name in self.completed_stages

    def get_result(self, stage_name: str) -> Optional[Any]:
        """Get the stored result for a completed stage.

        Parameters
        ----------
        stage_name : str
            Name of the stage

        Returns
        -------
        Any or None
            The stored result, or None if not found
        """
        with self._lock:
            return self.stage_results.get(stage_name)

    def update_data(self, new_data: Any) -> None:
        """Update the primary data artifact.

        Parameters
        ----------
        new_data : Any
            New data to set as the primary artifact
        """
        with self._lock:
            self.data = new_data

    def add_report_path(self, report_type: str, path: Path) -> None:
        """Add a generated report path.

        Parameters
        ----------
        report_type : str
            Type of report (e.g., 'html', 'excel', 'igv')
        path : Path
            Path to the generated report
        """
        with self._lock:
            self.report_paths[report_type] = path

    def get_execution_time(self) -> float:
        """Get the elapsed execution time in seconds.

        Returns
        -------
        float
            Elapsed time since pipeline start
        """
        return (datetime.now() - self.start_time).total_seconds()

    def __repr__(self) -> str:
        """Return string representation showing key state information."""
        return (
            f"PipelineContext("
            f"stages_completed={len(self.completed_stages)}, "
            f"has_data={self.data is not None}, "
            f"execution_time={self.get_execution_time():.1f}s)"
        )
