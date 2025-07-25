"""
Enhanced parallel genotype replacement with intelligent chunking and streaming pipeline.

This module provides a high-performance alternative to the existing chunked vectorized approach,
focusing on optimal thread utilization and realistic memory management.
"""

import gzip
import logging
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from queue import Empty, PriorityQueue, Queue
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd
import psutil

logger = logging.getLogger(__name__)


def get_available_memory_gb() -> float:
    """Get available system memory in GB."""
    try:
        memory = psutil.virtual_memory()
        available_gb = memory.available / (1024**3)
        return available_gb
    except Exception as e:
        logger.warning(f"Could not detect available memory: {e}. Using default estimate of 8GB")
        return 8.0


def open_file(
    file_path: Union[str, Path], mode: str = "rt", compression: Optional[str] = None
) -> Any:
    """
    Open a file with optional compression support.

    Parameters
    ----------
    file_path : Union[str, Path]
        Path to the file
    mode : str, default 'rt'
        File open mode
    compression : str, optional
        Compression type ('gzip' for .gz files)

    Returns
    -------
    File handle
    """
    file_path = Path(file_path)

    if compression == "gzip" or str(file_path).endswith(".gz"):
        if "b" not in mode:
            # Convert text mode to binary equivalent for gzip
            if mode == "rt":
                return gzip.open(file_path, "rt", encoding="utf-8")
            elif mode == "wt":
                return gzip.open(file_path, "wt", encoding="utf-8")
            else:
                return gzip.open(file_path, mode, encoding="utf-8")
        else:
            return gzip.open(file_path, mode)
    else:
        return open(file_path, mode, encoding="utf-8" if "b" not in mode else None)


class AdaptiveMemoryManager:
    """
    Real-time memory profiling and adaptive chunk sizing for genotype replacement.

    This class learns from actual memory usage patterns to provide more accurate
    estimates than the conservative static calculations used previously.
    """

    def __init__(self):
        self.memory_samples = []
        self.processing_stats = {}
        self.baseline_memory = None

    def profile_chunk_memory(self, chunk_df: pd.DataFrame, sample_count: int) -> Dict[str, float]:
        """
        Profile actual memory usage during chunk processing.

        Parameters
        ----------
        chunk_df : pd.DataFrame
            Sample chunk to profile
        sample_count : int
            Number of samples being processed

        Returns
        -------
        Dict[str, float]
            Memory profiling statistics
        """
        process = psutil.Process()

        # Record baseline memory
        if self.baseline_memory is None:
            self.baseline_memory = process.memory_info().rss

        initial_memory = process.memory_info().rss

        # Create a copy to measure pandas overhead
        chunk_copy = chunk_df.copy()
        copy_memory = process.memory_info().rss

        # Simulate basic string operations
        if "GT" in chunk_copy.columns:
            _ = chunk_copy["GT"].str.split(",", expand=True)

        peak_memory = process.memory_info().rss

        # Calculate actual metrics
        copy_overhead = copy_memory - initial_memory
        processing_overhead = peak_memory - copy_memory
        total_memory = peak_memory - initial_memory

        stats = {
            "bytes_per_row_actual": total_memory / len(chunk_df) if len(chunk_df) > 0 else 0,
            "pandas_overhead_actual": total_memory / copy_overhead if copy_overhead > 0 else 3.0,
            "base_memory_per_row": copy_overhead / len(chunk_df) if len(chunk_df) > 0 else 0,
            "processing_memory_per_row": (
                processing_overhead / len(chunk_df) if len(chunk_df) > 0 else 0
            ),
            "sample_efficiency": sample_count / total_memory if total_memory > 0 else 1e-6,
            "total_memory_mb": total_memory / (1024 * 1024),
        }

        self.memory_samples.append(stats)

        # Keep only last 10 samples for adaptive learning
        while len(self.memory_samples) > 10:
            self.memory_samples.pop(0)

        return stats

    def get_realistic_memory_estimate(self, sample_count: int) -> Dict[str, float]:
        """
        Get realistic memory estimates based on profiling history.

        Parameters
        ----------
        sample_count : int
            Number of samples in the dataset

        Returns
        -------
        Dict[str, float]
            Memory estimation parameters
        """
        if not self.memory_samples:
            # Conservative fallback estimates
            return {
                "bytes_per_row": 500 + sample_count * 50,
                "pandas_overhead": 3.5,
                "confidence": 0.3,  # Low confidence without profiling
            }

        # Use median values from samples for stability
        bytes_per_row_samples = [s["bytes_per_row_actual"] for s in self.memory_samples]
        overhead_samples = [s["pandas_overhead_actual"] for s in self.memory_samples]

        return {
            "bytes_per_row": np.median(bytes_per_row_samples),
            "pandas_overhead": np.median(overhead_samples),
            "confidence": min(
                1.0, len(self.memory_samples) / 5.0
            ),  # Higher confidence with more samples
        }

    def calculate_safe_chunk_size(
        self,
        available_memory_gb: float,
        sample_count: int,
        target_chunks: int,
        file_rows_estimate: int,
    ) -> int:
        """
        Calculate optimal chunk size balancing memory safety and parallelization.

        Parameters
        ----------
        available_memory_gb : float
            Available system memory in GB
        sample_count : int
            Number of samples
        target_chunks : int
            Desired number of chunks for parallelization
        file_rows_estimate : int
            Estimated total rows in file

        Returns
        -------
        int
            Optimal chunk size
        """
        estimates = self.get_realistic_memory_estimate(sample_count)

        # Memory-based limit (use 60% of available memory for processing)
        safe_memory_bytes = available_memory_gb * 0.6 * (1024**3)
        memory_chunk_limit = int(safe_memory_bytes / estimates["bytes_per_row"])

        # Parallelization-based target
        parallel_chunk_target = max(1000, file_rows_estimate // target_chunks)

        # Performance-based limits - optimized for high-memory systems
        min_chunk_size = 1000  # Minimum for efficiency (increased from 500)
        max_chunk_size = (
            100000  # Maximum for responsiveness (increased from 25K for 256GB+ systems)
        )

        # Choose optimal size
        optimal_size = min(memory_chunk_limit, parallel_chunk_target)
        final_size = max(min_chunk_size, min(optimal_size, max_chunk_size))

        logger.info(
            f"Chunk sizing: memory_limit={memory_chunk_limit:,}, "
            f"parallel_target={parallel_chunk_target:,}, "
            f"final_size={final_size:,} (confidence={estimates['confidence']:.1f})"
        )

        return final_size


class IntelligentChunkCalculator:
    """Thread-aware chunking strategy that optimizes for CPU utilization."""

    @staticmethod
    def calculate_optimal_chunks(
        total_rows: int,
        threads: int,
        memory_gb: float,
        sample_count: int,
        memory_manager: Optional[AdaptiveMemoryManager] = None,
    ) -> Tuple[int, int]:
        """
        Calculate optimal chunk size and count for maximum parallelization.

        Parameters
        ----------
        total_rows : int
            Total number of rows to process
        threads : int
            Available thread count
        memory_gb : float
            Available memory in GB
        sample_count : int
            Number of samples
        memory_manager : AdaptiveMemoryManager, optional
            Memory profiling manager for adaptive sizing

        Returns
        -------
        Tuple[int, int]
            (chunk_size, expected_chunk_count)
        """
        # Target 3-4 chunks per thread for good load balancing
        target_chunks = threads * 3

        if memory_manager is not None:
            chunk_size = memory_manager.calculate_safe_chunk_size(
                memory_gb, sample_count, target_chunks, total_rows
            )
        else:
            # Fallback calculation without profiling
            estimated_bytes_per_row = 500 + sample_count * 60
            memory_bytes = memory_gb * 0.75 * (1024**3)  # Use 75% of memory (increased from 50%)
            memory_chunk_limit = int(memory_bytes / estimated_bytes_per_row)
            parallel_chunk_target = max(1000, total_rows // target_chunks)
            chunk_size = max(
                1000, min(memory_chunk_limit, parallel_chunk_target, 100000)
            )  # Increased limits

        expected_chunks = math.ceil(total_rows / chunk_size)

        logger.info(
            f"Calculated chunking: {total_rows:,} rows → {expected_chunks} chunks "
            f"of ~{chunk_size:,} rows each ({threads} threads available)"
        )

        return chunk_size, expected_chunks

    @staticmethod
    def estimate_file_rows(file_path: Union[str, Path]) -> int:
        """
        Fast estimation of file row count for chunking calculations.

        Parameters
        ----------
        file_path : Union[str, Path]
            Path to input file

        Returns
        -------
        int
            Estimated row count
        """
        try:
            # Sample-based estimation for large files
            sample_lines = 1000
            total_bytes = os.path.getsize(file_path)

            with open_file(file_path, "rt") as f:
                # Skip header
                header = f.readline()
                header_bytes = len(header.encode("utf-8"))

                # Sample lines to estimate average line length
                sample_bytes = 0
                lines_read = 0
                for line in f:
                    sample_bytes += len(line.encode("utf-8"))
                    lines_read += 1
                    if lines_read >= sample_lines:
                        break

            if lines_read == 0:
                return 0

            avg_line_bytes = sample_bytes / lines_read
            data_bytes = total_bytes - header_bytes
            estimated_rows = int(data_bytes / avg_line_bytes)

            logger.debug(f"File row estimation: {total_bytes} bytes → ~{estimated_rows:,} rows")
            return estimated_rows

        except Exception as e:
            logger.warning(f"Could not estimate file rows: {e}. Using conservative estimate.")
            return 100000  # Conservative fallback


class StreamingGenotypeProcessor:
    """
    High-performance streaming genotype processor with overlapped I/O and computation.

    This processor implements a producer-consumer pipeline where:
    - Producer threads handle file I/O
    - Worker processes handle genotype replacement
    - Consumer thread handles result writing
    """

    def __init__(self, threads: int = 1):
        self.threads = threads
        self.memory_manager = AdaptiveMemoryManager()

    def process_file_streaming(
        self,
        input_path: Union[str, Path],
        output_path: Union[str, Path],
        config: Dict[str, Any],
        chunk_size: Optional[int] = None,
    ) -> None:
        """
        Process genotype replacement with streaming pipeline architecture.

        Parameters
        ----------
        input_path : Union[str, Path]
            Input TSV file path
        output_path : Union[str, Path]
            Output file path
        config : Dict[str, Any]
            Processing configuration
        chunk_size : int, optional
            Override chunk size calculation
        """
        input_path = Path(input_path)
        output_path = Path(output_path)

        # Estimate file size and calculate optimal chunking
        if chunk_size is None:
            total_rows = IntelligentChunkCalculator.estimate_file_rows(input_path)
            sample_count = len(config.get("sample_list", "").split(","))
            available_memory = get_available_memory_gb()

            chunk_size, expected_chunks = IntelligentChunkCalculator.calculate_optimal_chunks(
                total_rows, self.threads, available_memory, sample_count, self.memory_manager
            )

        logger.info(f"Starting streaming genotype processing with {self.threads} threads")
        logger.info(f"Processing: {input_path} → {output_path}")

        # Thread allocation strategy
        io_threads = max(1, min(3, self.threads // 5))  # 15-20% for I/O
        cpu_workers = max(1, self.threads - io_threads - 1)  # Rest for processing, 1 for writing

        logger.info(f"Thread allocation: {io_threads} I/O, {cpu_workers} processing, 1 writing")

        # Process with streaming pipeline
        self._execute_streaming_pipeline(
            input_path, output_path, config, chunk_size, io_threads, cpu_workers
        )

    def _execute_streaming_pipeline(
        self,
        input_path: Path,
        output_path: Path,
        config: Dict[str, Any],
        chunk_size: int,
        io_threads: int,
        cpu_workers: int,
    ) -> None:
        """Execute the streaming pipeline with overlapped I/O and processing."""
        # Queues for pipeline stages
        chunk_queue = Queue(maxsize=cpu_workers * 2)  # Buffer chunks
        result_queue = PriorityQueue()  # Ordered results

        compression = "gzip" if str(input_path).endswith(".gz") else None
        output_compression = "gzip" if str(output_path).endswith(".gz") else None

        try:
            with ThreadPoolExecutor(io_threads + 1) as io_pool, ProcessPoolExecutor(
                cpu_workers
            ) as cpu_pool:

                # Start producer (chunk reader)
                producer_future = io_pool.submit(
                    self._chunk_producer, input_path, chunk_queue, chunk_size, compression
                )

                # Start consumers (chunk processors)
                consumer_futures = []
                for worker_id in range(cpu_workers):
                    future = cpu_pool.submit(
                        _process_chunks_worker, chunk_queue, result_queue, config, worker_id
                    )
                    consumer_futures.append(future)

                # Start result writer
                writer_future = io_pool.submit(
                    self._result_writer, result_queue, output_path, output_compression
                )

                # Monitor progress and handle completion
                self._monitor_pipeline_progress(
                    producer_future, consumer_futures, writer_future, result_queue
                )

        except Exception as e:
            logger.error(f"Streaming pipeline failed: {e}")
            raise

    def _chunk_producer(
        self, input_path: Path, chunk_queue: Queue, chunk_size: int, compression: Optional[str]
    ) -> None:
        """Read file and create chunks for processing."""
        try:
            chunk_id = 0
            profiling_done = False

            for chunk_df in pd.read_csv(
                input_path, sep="\t", dtype=str, compression=compression, chunksize=chunk_size
            ):
                # Profile first chunk for adaptive memory management
                if not profiling_done and len(chunk_df) > 0:
                    sample_count = len(chunk_df.columns) - 10  # Approximate sample count
                    self.memory_manager.profile_chunk_memory(chunk_df.head(100), sample_count)
                    profiling_done = True

                chunk_queue.put((chunk_id, chunk_df))
                chunk_id += 1

            # Signal completion
            for _ in range(chunk_queue.maxsize):  # Send stop signals to all workers
                chunk_queue.put((None, None))

            logger.info(f"Producer completed: {chunk_id} chunks queued")

        except Exception as e:
            logger.error(f"Producer failed: {e}")
            raise

    def _result_writer(
        self, result_queue: PriorityQueue, output_path: Path, compression: Optional[str]
    ) -> None:
        """Consumer thread: write results in order."""
        try:
            expected_chunk = 0
            pending_results = {}

            with open_file(output_path, "wt", compression=compression) as out_file:
                header_written = False

                while True:
                    try:
                        chunk_id, chunk_result = result_queue.get(timeout=5.0)

                        if chunk_result is None:  # Stop signal
                            break

                        pending_results[chunk_id] = chunk_result

                        # Write results in order
                        while expected_chunk in pending_results:
                            result_df = pending_results.pop(expected_chunk)

                            if not header_written:
                                result_df.to_csv(out_file, sep="\t", index=False, header=True)
                                header_written = True
                            else:
                                result_df.to_csv(out_file, sep="\t", index=False, header=False)

                            expected_chunk += 1

                    except Empty:
                        continue  # Check for more results

            logger.info(f"Writer completed: {expected_chunk} chunks written")

        except Exception as e:
            logger.error(f"Writer failed: {e}")
            raise

    def _monitor_pipeline_progress(
        self, producer_future, consumer_futures, writer_future, result_queue
    ):
        """Monitor pipeline execution and handle errors."""
        try:
            # Wait for producer to complete
            producer_future.result()
            logger.debug("Producer completed")

            # Wait for all consumers
            completed_consumers = 0
            for future in as_completed(consumer_futures):
                future.result()
                completed_consumers += 1
                logger.debug(f"Consumer {completed_consumers}/{len(consumer_futures)} completed")

            # Signal writer completion
            result_queue.put((None, None))  # Stop signal for writer

            # Wait for writer
            writer_future.result()
            logger.debug("Writer completed")

            logger.info("Streaming pipeline completed successfully")

        except Exception as e:
            logger.error(f"Pipeline monitoring failed: {e}")
            # Cancel remaining futures
            producer_future.cancel()
            for future in consumer_futures:
                future.cancel()
            writer_future.cancel()
            raise


def _process_chunks_worker(
    chunk_queue: Queue, result_queue: PriorityQueue, config: Dict[str, Any], worker_id: int
) -> None:
    """
    Worker process for chunk processing.

    This function runs in a separate process to handle CPU-intensive genotype replacement.
    """
    from .vectorized_replacer import VectorizedGenotypeReplacer

    logger.info(f"Worker {worker_id} started")
    replacer = VectorizedGenotypeReplacer(config)

    processed_chunks = 0

    try:
        while True:
            try:
                chunk_id, chunk_df = chunk_queue.get(timeout=2.0)

                if chunk_df is None:  # Stop signal
                    break

                # Process chunk
                start_time = time.time()
                processed_chunk = replacer.process_dataframe(chunk_df)
                process_time = time.time() - start_time

                # Send result
                result_queue.put((chunk_id, processed_chunk))
                processed_chunks += 1

                logger.debug(
                    f"Worker {worker_id}: processed chunk {chunk_id} "
                    f"({len(chunk_df)} rows in {process_time:.2f}s)"
                )

            except Empty:
                continue  # Keep checking for work

    except Exception as e:
        logger.error(f"Worker {worker_id} failed: {e}")
        raise
    finally:
        logger.info(f"Worker {worker_id} completed: {processed_chunks} chunks processed")


def process_streaming_parallel(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    config: Dict[str, Any],
    threads: int = 1,
) -> None:
    """
    High-level API for streaming parallel genotype replacement.

    This function provides a simple interface for the enhanced parallel processing
    that can be easily integrated into existing pipeline stages.

    Parameters
    ----------
    input_path : Union[str, Path]
        Path to input TSV file
    output_path : Union[str, Path]
        Path to output file
    config : Dict[str, Any]
        Processing configuration (compatible with vectorized_replacer)
    threads : int, default 1
        Number of threads to use for processing
    """
    processor = StreamingGenotypeProcessor(threads)
    processor.process_file_streaming(input_path, output_path, config)
