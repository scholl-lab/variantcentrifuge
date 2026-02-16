"""
Pipeline-wide Resource Manager for VariantCentrifuge

This module provides pipeline-wide resource management for any stage,
detecting available memory and CPUs, calculating optimal chunk sizes
and worker counts for parallel processing.
"""

import logging
import os
from pathlib import Path
from typing import Any

import psutil

logger = logging.getLogger(__name__)


class ResourceManager:
    """
    Pipeline-wide resource manager that:
    1. Detects actual system memory from CLI, SLURM, PBS, cgroups, or psutil
    2. Detects CPU cores
    3. Calculates optimal chunk sizes for dataset processing
    4. Calculates optimal worker counts for parallel processing
    5. Provides small-dataset threshold for parallelization decisions
    """

    def __init__(
        self,
        config: dict[str, Any] | None = None,
        overhead_factor: float = 3.0,
        memory_safety_factor: float = 0.80,
        min_items_for_parallel: int = 100,
    ):
        """
        Initialize the resource manager.

        Args:
            config: Configuration dictionary with memory settings
            overhead_factor: Memory overhead multiplier (default 3.0 for Pandas/analysis)
            memory_safety_factor: Fraction of memory to use (default 0.80 = 80%)
            min_items_for_parallel: Minimum items to enable parallelization (default 100)
        """
        self.config = config or {}
        self.overhead_factor = overhead_factor
        self.memory_safety_factor = memory_safety_factor
        self.min_items_for_parallel = min_items_for_parallel

        # Detect resources
        self.memory_gb = self._detect_memory()
        self.cpu_cores = self._detect_cpus()

        # Log detected resources
        source = self._get_memory_source()
        logger.info(
            f"ResourceManager: {self.cpu_cores} CPUs, {self.memory_gb:.1f}GB available ({source})"
        )

    def _get_memory_source(self) -> str:
        """Get description of memory detection source for logging."""
        if self.config.get("max_memory_gb"):
            return "CLI"
        if os.getenv("SLURM_MEM_PER_NODE"):
            return "SLURM"
        if os.getenv("PBS_RESC_MEM"):
            return "PBS"
        cgroup_limit = self._get_cgroup_memory_limit()
        if cgroup_limit:
            return "cgroup"
        return "psutil"

    def _detect_memory(self) -> float:
        """
        Detect available memory limit in GB.

        Priority:
        1. --max-memory-gb CLI parameter
        2. SLURM_MEM_PER_NODE environment variable
        3. PBS_RESC_MEM environment variable
        4. cgroup limits (v1 and v2)
        5. psutil available memory (fallback)

        Returns:
            Memory in GB
        """
        # 1. Check CLI parameter
        cli_memory = self.config.get("max_memory_gb")
        if cli_memory:
            logger.debug(f"Using CLI-specified memory limit: {cli_memory:.1f}GB")
            return float(cli_memory)

        # 2. Check SLURM
        slurm_mem = os.getenv("SLURM_MEM_PER_NODE")
        if slurm_mem:
            try:
                # SLURM memory is in MB
                slurm_gb = float(slurm_mem) / 1024
                logger.debug(f"Using SLURM allocated memory: {slurm_gb:.1f}GB")
                return slurm_gb
            except (ValueError, TypeError):
                logger.warning(f"Invalid SLURM_MEM_PER_NODE value: {slurm_mem}")

        # 3. Check PBS/Torque
        pbs_mem = os.getenv("PBS_RESC_MEM")
        if pbs_mem:
            try:
                # Parse PBS memory (could be like "16gb" or "16000mb")
                pbs_mem_lower = pbs_mem.lower()
                if pbs_mem_lower.endswith("gb"):
                    pbs_gb = float(pbs_mem_lower.replace("gb", ""))
                elif pbs_mem_lower.endswith("mb"):
                    pbs_gb = float(pbs_mem_lower.replace("mb", "")) / 1024
                else:
                    pbs_gb = float(pbs_mem) / 1024  # Assume MB
                logger.debug(f"Using PBS allocated memory: {pbs_gb:.1f}GB")
                return pbs_gb
            except (ValueError, TypeError):
                logger.warning(f"Invalid PBS_RESC_MEM value: {pbs_mem}")

        # 4. Check cgroup limits (containers/HPC)
        try:
            cgroup_limit = self._get_cgroup_memory_limit()
            if cgroup_limit:
                logger.debug(f"Using cgroup memory limit: {cgroup_limit:.1f}GB")
                return cgroup_limit
        except Exception as e:
            logger.debug(f"Could not read cgroup limits: {e}")

        # 5. Fallback to available system memory
        try:
            memory_info = psutil.virtual_memory()
            available_gb = memory_info.available / (1024**3)
            logger.debug(f"Using detected available memory: {available_gb:.1f}GB")
            return float(available_gb)
        except Exception as e:
            logger.warning(f"Could not detect memory: {e}. Using conservative 8GB")
            return 8.0

    def _get_cgroup_memory_limit(self) -> float | None:
        """
        Get memory limit from cgroup (containers/HPC).

        Returns:
            Memory limit in GB or None if not found
        """
        cgroup_paths = [
            "/sys/fs/cgroup/memory/memory.limit_in_bytes",  # cgroup v1
            "/sys/fs/cgroup/memory.max",  # cgroup v2
        ]

        for path in cgroup_paths:
            try:
                if Path(path).exists():
                    with open(path) as f:
                        limit_str = f.read().strip()
                        # cgroup v2 can have "max" string for unlimited
                        if limit_str == "max":
                            continue
                        limit_bytes = int(limit_str)
                        # Check if it's a real limit (not max value)
                        if limit_bytes < (1 << 62):  # Less than ~4 exabytes
                            return limit_bytes / (1024**3)
            except (OSError, ValueError) as e:
                logger.debug(f"Could not read {path}: {e}")

        return None

    def _detect_cpus(self) -> int:
        """
        Detect CPU core count.

        Returns:
            Number of physical CPU cores (fallback to 4 if detection fails)
        """
        try:
            # Prefer physical cores over logical (excludes hyperthreading)
            cores = psutil.cpu_count(logical=False)
            if cores:
                return cores
        except Exception as e:
            logger.debug(f"psutil.cpu_count(logical=False) failed: {e}")

        # Fallback to os.cpu_count
        try:
            cores = os.cpu_count()
            if cores:
                return cores
        except Exception as e:
            logger.debug(f"os.cpu_count() failed: {e}")

        # Conservative fallback
        logger.warning("Could not detect CPU count, using conservative 4 cores")
        return 4

    def auto_chunk_size(
        self,
        total_items: int,
        num_samples: int,
        bytes_per_item: int = 8,
    ) -> int:
        """
        Calculate optimal chunk size based on available memory.

        Args:
            total_items: Total number of items (e.g., variants)
            num_samples: Number of samples
            bytes_per_item: Bytes per item per sample (default 8 for float64)

        Returns:
            Optimal chunk size (number of items)
        """
        # Calculate safe memory budget
        safe_memory_bytes = self.memory_gb * self.memory_safety_factor * (1024**3)

        # Calculate memory per item (with overhead)
        memory_per_item = num_samples * bytes_per_item * self.overhead_factor

        # Calculate max items that fit in memory
        max_items = int(safe_memory_bytes / memory_per_item)

        # Apply reasonable bounds
        max_items = max(100, min(max_items, 1_000_000))  # 100 to 1M items

        # Return the smaller of total_items or max_items
        chunk_size = min(total_items, max_items)

        logger.debug(
            f"Auto chunk size: {chunk_size:,} items "
            f"(total={total_items:,}, samples={num_samples}, "
            f"memory={self.memory_gb:.1f}GB)"
        )

        return chunk_size

    def auto_workers(
        self,
        task_count: int,
        memory_per_task_gb: float,
    ) -> int:
        """
        Calculate optimal worker count for parallel processing.

        Args:
            task_count: Number of tasks to process
            memory_per_task_gb: Memory requirement per task in GB

        Returns:
            Optimal number of parallel workers
        """
        # Calculate safe memory budget
        safe_memory_gb = self.memory_gb * self.memory_safety_factor

        # How many tasks fit in memory?
        memory_constrained_workers = max(1, int(safe_memory_gb / memory_per_task_gb))

        # CPU-constrained workers
        cpu_constrained_workers = self.cpu_cores

        # Can't have more workers than tasks
        task_constrained_workers = task_count

        # Take minimum of all constraints
        workers = min(
            memory_constrained_workers,
            cpu_constrained_workers,
            task_constrained_workers,
        )

        logger.debug(
            f"Auto workers: {workers} "
            f"(memory_limit={memory_constrained_workers}, "
            f"cpu_limit={cpu_constrained_workers}, "
            f"task_limit={task_constrained_workers})"
        )

        return workers

    def should_parallelize(self, total_items: int) -> bool:
        """
        Determine if dataset is large enough to warrant parallelization.

        Args:
            total_items: Total number of items to process

        Returns:
            True if parallelization is recommended
        """
        return total_items >= self.min_items_for_parallel

    def estimate_memory(
        self,
        num_variants: int,
        num_samples: int,
        bytes_per_item: int = 8,
    ) -> float:
        """
        Estimate memory required for processing a dataset in GB.

        Args:
            num_variants: Number of variants
            num_samples: Number of samples
            bytes_per_item: Bytes per item per sample (default 8 for float64)

        Returns:
            Estimated memory in GB
        """
        memory_bytes = num_variants * num_samples * bytes_per_item * self.overhead_factor
        memory_gb = memory_bytes / (1024**3)
        return memory_gb

    def get_summary(self) -> dict[str, Any]:
        """
        Get summary of detected resources for logging.

        Returns:
            Dictionary with resource information
        """
        return {
            "memory_gb": self.memory_gb,
            "memory_source": self._get_memory_source(),
            "cpu_cores": self.cpu_cores,
            "memory_safety_factor": self.memory_safety_factor,
            "overhead_factor": self.overhead_factor,
            "min_items_for_parallel": self.min_items_for_parallel,
            "safe_memory_gb": self.memory_gb * self.memory_safety_factor,
        }
