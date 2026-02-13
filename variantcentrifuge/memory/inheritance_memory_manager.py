"""
Intelligent Memory Manager for Inheritance Analysis

This module provides intelligent memory management for inheritance analysis,
taking into account actual system memory and calculating optimal chunk sizes
and parallel processing limits.
"""

import logging
from pathlib import Path
from typing import Any

import psutil

logger = logging.getLogger(__name__)


class InheritanceMemoryManager:
    """
    Intelligent memory manager for inheritance analysis that:
    1. Detects actual system memory
    2. Calculates realistic memory requirements
    3. Determines optimal chunk sizes and parallel processing
    4. Adapts to different sample counts and system capabilities
    """

    def __init__(self, config: dict[str, Any] | None = None):
        """
        Initialize the memory manager.

        Args:
            config: Configuration dictionary with memory settings
        """
        self.config = config or {}

        # Get memory limit from config/CLI or detect available memory
        self._allocated_memory_gb = self._get_allocated_memory()
        self._current_available_gb = self._get_current_available_memory()

        # Memory safety factors - optimized for modern HPC/high-memory systems
        self.memory_safety_factor = self.config.get(
            "memory_safety_factor", 0.92
        )  # Use max 92% of allocated (optimized for dedicated systems)
        self.inheritance_memory_fraction = self.config.get(
            "inheritance_memory_fraction", 0.85
        )  # 85% of safe memory for inheritance (optimized for large cohorts)

        # Memory estimation constants (calibrated from actual usage)
        self.bytes_per_sample_column = 8  # float64 for genotype values
        self.pandas_overhead_factor = 2.5  # Pandas DataFrame overhead
        self.inheritance_analysis_factor = 3.0  # Additional memory during analysis

        logger.info(
            f"Memory Manager initialized: {self._allocated_memory_gb:.1f}GB allocated, "
            f"{self._current_available_gb:.1f}GB currently available"
        )

    def _get_allocated_memory(self) -> float:
        """
        Get allocated memory limit in GB.

        Priority:
        1. --max-memory-gb CLI parameter
        2. Environment variables (SLURM_MEM_PER_NODE, etc.)
        3. cgroup limits (containers/HPC)
        4. Available system memory (fallback)
        """
        # 1. Check CLI parameter
        cli_memory = self.config.get("max_memory_gb")
        if cli_memory:
            logger.info(f"Using CLI-specified memory limit: {cli_memory:.1f}GB")
            return float(cli_memory)

        # 2. Check environment variables (HPC schedulers)
        import os

        # SLURM
        slurm_mem = os.getenv("SLURM_MEM_PER_NODE")
        if slurm_mem:
            try:
                # SLURM memory is in MB
                slurm_gb = float(slurm_mem) / 1024
                logger.info(f"Using SLURM allocated memory: {slurm_gb:.1f}GB")
                return slurm_gb
            except (ValueError, TypeError):
                logger.warning(f"Invalid SLURM_MEM_PER_NODE value: {slurm_mem}")

        # PBS/Torque
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
                logger.info(f"Using PBS allocated memory: {pbs_gb:.1f}GB")
                return pbs_gb
            except (ValueError, TypeError):
                logger.warning(f"Invalid PBS_RESC_MEM value: {pbs_mem}")

        # 3. Check cgroup limits (containers)
        try:
            cgroup_limit = self._get_cgroup_memory_limit()
            if cgroup_limit:
                logger.info(f"Using cgroup memory limit: {cgroup_limit:.1f}GB")
                return cgroup_limit
        except Exception as e:
            logger.debug(f"Could not read cgroup limits: {e}")

        # 4. Fallback to available system memory
        try:
            memory_info = psutil.virtual_memory()
            available_gb = memory_info.available / (1024**3)
            logger.info(f"Using detected available memory: {available_gb:.1f}GB")
            return float(available_gb)
        except Exception as e:
            logger.warning(f"Could not detect memory: {e}. Using conservative 8GB")
            return 8.0

    def _get_cgroup_memory_limit(self) -> float | None:
        """Get memory limit from cgroup (containers/HPC)."""
        cgroup_paths = [
            "/sys/fs/cgroup/memory/memory.limit_in_bytes",  # cgroup v1
            "/sys/fs/cgroup/memory.max",  # cgroup v2
        ]

        for path in cgroup_paths:
            try:
                if Path(path).exists():
                    with open(path) as f:
                        limit_bytes = int(f.read().strip())
                        # Check if it's a real limit (not max value)
                        if limit_bytes < (1 << 62):  # Less than ~4 exabytes
                            return limit_bytes / (1024**3)
            except (OSError, ValueError) as e:
                logger.debug(f"Could not read {path}: {e}")

        return None

    def _get_current_available_memory(self) -> float:
        """Get currently available memory in GB (for monitoring)."""
        try:
            memory_info = psutil.virtual_memory()
            available_gb = memory_info.available / (1024**3)
            # Don't exceed our allocated limit
            return float(min(available_gb, self._allocated_memory_gb))
        except Exception as e:
            logger.warning(f"Could not get current memory: {e}. Using allocated limit")
            return self._allocated_memory_gb

    def estimate_chunk_memory_requirement(self, num_variants: int, num_samples: int) -> float:
        """
        Estimate memory requirement for a chunk in GB.

        Args:
            num_variants: Number of variants in chunk
            num_samples: Number of samples

        Returns:
            Estimated memory requirement in GB
        """
        # Base memory for sample columns (variants x samples x bytes_per_cell)
        sample_columns_memory = num_variants * num_samples * self.bytes_per_sample_column

        # Apply overhead factors
        with_pandas_overhead = sample_columns_memory * self.pandas_overhead_factor
        with_analysis_overhead = with_pandas_overhead * self.inheritance_analysis_factor

        # Convert to GB
        memory_gb = with_analysis_overhead / (1024**3)

        logger.debug(
            f"Memory estimate for {num_variants} variants x {num_samples} samples: "
            f"{memory_gb:.2f}GB"
        )

        return memory_gb

    def calculate_max_chunk_size(
        self, num_samples: int, target_memory_gb: float | None = None
    ) -> int:
        """
        Calculate maximum chunk size (variants) for given sample count.

        Args:
            num_samples: Number of samples
            target_memory_gb: Target memory per chunk (uses intelligent default if None)

        Returns:
            Maximum number of variants per chunk
        """
        if target_memory_gb is None:
            # Use fraction of allocated memory
            safe_memory = self._allocated_memory_gb * self.memory_safety_factor
            target_memory_gb = safe_memory * self.inheritance_memory_fraction

        # Calculate variants that fit in target memory
        memory_per_variant = (
            num_samples
            * self.bytes_per_sample_column
            * self.pandas_overhead_factor
            * self.inheritance_analysis_factor
        ) / (1024**3)

        max_variants = int(target_memory_gb / memory_per_variant)

        # Apply reasonable bounds
        max_variants = max(100, min(max_variants, 1_000_000))  # Between 100 and 1M variants

        logger.info(
            f"Calculated max chunk size for {num_samples} samples: "
            f"{max_variants:,} variants (target: {target_memory_gb:.1f}GB)"
        )

        return max_variants

    def calculate_optimal_parallelism(
        self, chunk_sizes: list, num_samples: int
    ) -> tuple[int, dict[str, Any]]:
        """
        Calculate optimal number of parallel workers based on memory constraints.

        Args:
            chunk_sizes: List of chunk sizes (number of variants per chunk)
            num_samples: Number of samples

        Returns:
            Tuple of (max_workers, memory_info_dict)
        """
        safe_memory_gb = self._allocated_memory_gb * self.memory_safety_factor
        inheritance_memory_gb = safe_memory_gb * self.inheritance_memory_fraction

        # Calculate memory requirement for largest chunk
        max_chunk_size = max(chunk_sizes) if chunk_sizes else 1000
        memory_per_chunk = self.estimate_chunk_memory_requirement(max_chunk_size, num_samples)

        # Calculate how many chunks can run in parallel
        max_parallel_chunks = max(1, int(inheritance_memory_gb / memory_per_chunk))

        # Consider system CPU cores as upper limit
        cpu_cores = psutil.cpu_count(logical=False) or 4
        max_workers = min(max_parallel_chunks, cpu_cores, len(chunk_sizes))

        memory_info = {
            "allocated_memory_gb": self._allocated_memory_gb,
            "current_available_gb": self._current_available_gb,
            "safe_memory_gb": safe_memory_gb,
            "inheritance_memory_gb": inheritance_memory_gb,
            "memory_per_chunk_gb": memory_per_chunk,
            "max_parallel_by_memory": max_parallel_chunks,
            "max_parallel_by_cpu": cpu_cores,
            "recommended_workers": max_workers,
        }

        logger.info(
            f"Optimal parallelism: {max_workers} workers "
            f"(memory limit: {max_parallel_chunks}, CPU cores: {cpu_cores})"
        )

        return max_workers, memory_info

    def should_process_chunk(
        self, num_variants: int, num_samples: int, force_processing: bool = False
    ) -> tuple[bool, str]:
        """
        Determine if a chunk should be processed based on memory constraints.

        Args:
            num_variants: Number of variants in chunk
            num_samples: Number of samples
            force_processing: If True, will try to process even if risky

        Returns:
            Tuple of (should_process, reason)
        """
        estimated_memory = self.estimate_chunk_memory_requirement(num_variants, num_samples)
        safe_memory = self._allocated_memory_gb * self.memory_safety_factor

        if estimated_memory <= safe_memory:
            return True, f"Safe: {estimated_memory:.2f}GB <= {safe_memory:.2f}GB"
        elif force_processing:
            if estimated_memory <= self._allocated_memory_gb:
                return (
                    True,
                    f"Risky but forced: {estimated_memory:.2f}GB <= "
                    f"{self._allocated_memory_gb:.2f}GB",
                )
            else:
                return (
                    False,
                    f"Cannot process: {estimated_memory:.2f}GB > "
                    f"{self._allocated_memory_gb:.2f}GB allocated",
                )
        else:
            return False, f"Unsafe: {estimated_memory:.2f}GB > {safe_memory:.2f}GB safe limit"

    def get_memory_strategy(self, total_variants: int, num_samples: int) -> dict[str, Any]:
        """
        Get comprehensive memory strategy for inheritance analysis.

        Args:
            total_variants: Total number of variants
            num_samples: Number of samples

        Returns:
            Dictionary with complete memory strategy
        """
        # Calculate if full dataset can be processed in memory
        full_memory = self.estimate_chunk_memory_requirement(total_variants, num_samples)
        safe_memory = self._allocated_memory_gb * self.memory_safety_factor

        can_process_full = full_memory <= safe_memory

        if can_process_full:
            strategy = {
                "approach": "full_dataset",
                "total_variants": total_variants,
                "estimated_memory_gb": full_memory,
                "chunks": 1,
                "max_workers": 1,
                "reason": f"Full dataset fits in memory "
                f"({full_memory:.2f}GB <= {safe_memory:.2f}GB)",
            }
        else:
            # Calculate chunking strategy
            max_chunk_size = self.calculate_max_chunk_size(num_samples)
            num_chunks = max(1, (total_variants + max_chunk_size - 1) // max_chunk_size)
            actual_chunk_sizes = []

            for i in range(num_chunks):
                start_idx = i * max_chunk_size
                end_idx = min((i + 1) * max_chunk_size, total_variants)
                chunk_size = end_idx - start_idx
                actual_chunk_sizes.append(chunk_size)

            max_workers, memory_info = self.calculate_optimal_parallelism(
                actual_chunk_sizes, num_samples
            )

            strategy = {
                "approach": "chunked",
                "total_variants": total_variants,
                "max_chunk_size": max_chunk_size,
                "chunks": num_chunks,
                "chunk_sizes": actual_chunk_sizes,
                "max_workers": max_workers,
                "estimated_memory_per_chunk_gb": self.estimate_chunk_memory_requirement(
                    max_chunk_size, num_samples
                ),
                "total_memory_budget_gb": safe_memory * self.inheritance_memory_fraction,
                "memory_info": memory_info,
                "reason": f"Dataset too large for memory "
                f"({full_memory:.2f}GB > {safe_memory:.2f}GB)",
            }

        logger.info(
            f"Memory strategy: {strategy['approach']} with {strategy['chunks']} chunk(s), "
            f"{strategy.get('max_workers', 1)} worker(s)"
        )

        return strategy

    def log_memory_status(self):
        """Log current memory status for debugging."""
        try:
            memory = psutil.virtual_memory()
            inheritance_budget = (
                self._allocated_memory_gb
                * self.memory_safety_factor
                * self.inheritance_memory_fraction
            )
            logger.info(
                f"Memory Status - System: {memory.total / (1024**3):.1f}GB, "
                f"System Available: {memory.available / (1024**3):.1f}GB, "
                f"Allocated: {self._allocated_memory_gb:.1f}GB, "
                f"Safe limit: {self._allocated_memory_gb * self.memory_safety_factor:.1f}GB, "
                f"Inheritance budget: {inheritance_budget:.1f}GB"
            )
        except Exception as e:
            logger.warning(f"Could not log memory status: {e}")
