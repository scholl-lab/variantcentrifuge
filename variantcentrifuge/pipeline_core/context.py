"""
PipelineContext - Single source of truth for pipeline state and data.

This module provides the PipelineContext dataclass that flows through all stages,
carrying configuration, state, and data artifacts throughout the pipeline execution.
"""

import argparse
import logging
import threading
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

if TYPE_CHECKING:
    from .workspace import Workspace

logger = logging.getLogger(__name__)


class PicklableLock:
    """A lock that can be pickled by creating a new lock on unpickling."""

    def __init__(self):
        """Initialize a new picklable lock."""
        self._lock = threading.RLock()

    def __getstate__(self):
        """Return state for pickling, excluding the unpicklable RLock."""
        # Return empty state - we don't need to save the lock state
        return {}

    def __setstate__(self, state):
        """Restore state after unpickling by creating a fresh lock."""
        # Create a fresh lock on unpickling
        self._lock = threading.RLock()

    def __enter__(self):
        """Enter the context manager by acquiring the lock."""
        # Ensure lock exists (in case of edge cases)
        if not hasattr(self, "_lock"):
            self._lock = threading.RLock()
        return self._lock.__enter__()

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit the context manager by releasing the lock."""
        return self._lock.__exit__(exc_type, exc_val, exc_tb)

    def acquire(self, blocking=True, timeout=-1):
        """Acquire the lock with optional blocking and timeout."""
        if not hasattr(self, "_lock"):
            self._lock = threading.RLock()
        return self._lock.acquire(blocking, timeout)

    def release(self):
        """Release the lock."""
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
    genotype_replaced_tsv : Optional[Path]
        Path to TSV with genotypes replaced by sample IDs
    phenotypes_added_tsv : Optional[Path]
        Path to TSV with phenotype data integrated
    extra_columns_removed_tsv : Optional[Path]
        Path to TSV with extra columns removed
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
    config: dict[str, Any]
    workspace: "Workspace"  # Forward reference, will be imported
    start_time: datetime = field(default_factory=datetime.now)

    # --- Mutable State & Data Artifacts ---
    # Core data that flows through pipeline
    data: Any = None  # Can be file path or DataFrame

    # Stage tracking
    completed_stages: set[str] = field(default_factory=set)
    stage_results: dict[str, Any] = field(default_factory=dict)

    # Loaded configurations
    vcf_samples: list[str] = field(default_factory=list)
    pedigree_data: dict | None = None
    phenotype_data: dict | None = None
    scoring_config: dict | None = None
    annotation_configs: dict[str, Any] = field(default_factory=dict)

    # Processing artifacts
    gene_bed_file: Path | None = None
    extracted_vcf: Path | None = None
    filtered_vcf: Path | None = None
    extracted_tsv: Path | None = None
    genotype_replaced_tsv: Path | None = None
    phenotypes_added_tsv: Path | None = None
    extra_columns_removed_tsv: Path | None = None
    chunked_analysis_tsv: Path | None = None

    # Analysis results
    current_dataframe: pd.DataFrame | None = None
    variants_df: pd.DataFrame | None = None  # Optimized DataFrame for in-memory pass-through
    column_rename_map: dict[str, str] = field(default_factory=dict)  # Column rename mapping
    statistics: dict[str, Any] = field(default_factory=dict)
    gene_burden_results: pd.DataFrame | None = None
    association_results: pd.DataFrame | None = None

    # Output paths
    final_output_path: Path | None = None
    report_paths: dict[str, Path] = field(default_factory=dict)

    # Checkpoint state (if enabled)
    checkpoint_state: Any | None = None

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

            # Note: checkpoint completion is handled by the stage itself with file tracking

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

    def get_result(self, stage_name: str) -> Any | None:
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

    def merge_from(self, other: "PipelineContext") -> None:
        """Merge updates from another context (e.g., from parallel execution).

        This method is used to merge state updates from parallel stage execution
        back into the main context. Only mutable state is merged.

        Parameters
        ----------
        other : PipelineContext
            Context with updates to merge
        """
        with self._lock:
            # Merge completed stages
            self.completed_stages.update(other.completed_stages)

            # Merge stage results
            self.stage_results.update(other.stage_results)

            # Update file paths if they were set in the other context
            if other.gene_bed_file and not self.gene_bed_file:
                self.gene_bed_file = other.gene_bed_file
            if other.extracted_vcf and not self.extracted_vcf:
                self.extracted_vcf = other.extracted_vcf
            if other.filtered_vcf and not self.filtered_vcf:
                self.filtered_vcf = other.filtered_vcf
            if other.extracted_tsv and not self.extracted_tsv:
                self.extracted_tsv = other.extracted_tsv
            if other.genotype_replaced_tsv and not self.genotype_replaced_tsv:
                self.genotype_replaced_tsv = other.genotype_replaced_tsv
            if other.phenotypes_added_tsv and not self.phenotypes_added_tsv:
                self.phenotypes_added_tsv = other.phenotypes_added_tsv
            if other.extra_columns_removed_tsv and not self.extra_columns_removed_tsv:
                self.extra_columns_removed_tsv = other.extra_columns_removed_tsv

            # Update analysis results if present
            if other.current_dataframe is not None and self.current_dataframe is None:
                self.current_dataframe = other.current_dataframe
            if other.variants_df is not None and self.variants_df is None:
                self.variants_df = other.variants_df
            if other.column_rename_map:
                self.column_rename_map.update(other.column_rename_map)
            if other.statistics:
                self.statistics.update(other.statistics)
            if other.gene_burden_results is not None and self.gene_burden_results is None:
                self.gene_burden_results = other.gene_burden_results
            if other.association_results is not None and self.association_results is None:
                self.association_results = other.association_results

            # Merge config updates (important for ConfigurationLoadingStage)
            if other.config:
                # Log config merging for debugging
                logger.debug(f"Merging config updates: {len(other.config)} keys from other context")
                self.config.update(other.config)

            # Update configurations loaded by parallel stages
            if other.pedigree_data and not self.pedigree_data:
                self.pedigree_data = other.pedigree_data
            if other.phenotype_data and not self.phenotype_data:
                self.phenotype_data = other.phenotype_data
            if other.scoring_config and not self.scoring_config:
                self.scoring_config = other.scoring_config
            if other.annotation_configs:
                self.annotation_configs.update(other.annotation_configs)

            # Update output paths
            if other.final_output_path and not self.final_output_path:
                self.final_output_path = other.final_output_path
            if other.report_paths:
                self.report_paths.update(other.report_paths)

            # Update vcf_samples if populated
            if other.vcf_samples and not self.vcf_samples:
                self.vcf_samples = other.vcf_samples

            logger.debug(
                f"Merged context updates: "
                f"{len(other.completed_stages)} completed stages, "
                f"{len(other.stage_results)} stage results"
            )

    def __repr__(self) -> str:
        """Return string representation showing key state information."""
        return (
            f"PipelineContext("
            f"stages_completed={len(self.completed_stages)}, "
            f"has_data={self.data is not None}, "
            f"execution_time={self.get_execution_time():.1f}s)"
        )
