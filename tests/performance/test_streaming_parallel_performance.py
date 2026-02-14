"""
Performance tests for streaming parallel genotype replacement.

.. deprecated::
    Pre-GSD legacy performance tests. Superseded by pytest-benchmark-based
    benchmarks in benchmark_genotype_replacement.py. Retained temporarily
    to ensure no regression; scheduled for removal in v0.14.0.

These tests compare the performance of the new streaming parallel processor
against the existing chunked vectorized approach.
"""

import tempfile
import time
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.parallel_genotype_processor import process_streaming_parallel
from variantcentrifuge.vectorized_replacer import process_parallel_chunked_vectorized


@pytest.mark.performance
@pytest.mark.slow
class TestStreamingParallelPerformance:
    """Performance comparison tests."""

    def setup_method(self):
        """Create test data for performance testing."""
        self.test_dir = tempfile.mkdtemp()

        # Create increasingly large test datasets
        self.small_file = self._create_test_file("small", rows=1000, samples=50)
        self.medium_file = self._create_test_file("medium", rows=10000, samples=200)
        self.large_file = self._create_test_file("large", rows=50000, samples=500)

        self.config = {
            "sample_list": ",".join([f"sample_{i}" for i in range(500)]),
            "separator": ";",
            "extract_fields_separator": ",",
            "append_extra_sample_fields": False,
            "extra_sample_fields": [],
            "genotype_replacement_map": {},
        }

    def teardown_method(self):
        """Clean up test files."""
        import shutil

        shutil.rmtree(self.test_dir, ignore_errors=True)

    def _create_test_file(self, name: str, rows: int, samples: int) -> Path:
        """Create a test TSV file with specified dimensions."""
        file_path = Path(self.test_dir) / f"test_{name}.tsv"

        # Generate genotype data
        genotypes = []
        for _ in range(rows):
            gt_values = []
            for _ in range(samples):
                # Random genotypes
                import random

                gt = random.choice(["0/0", "0/1", "1/1", "./."])
                gt_values.append(gt)
            genotypes.append(",".join(gt_values))

        # Create DataFrame
        data = {
            "CHR": [str(i % 22 + 1) for i in range(rows)],
            "POS": list(range(100000, 100000 + rows)),
            "REF": ["A"] * rows,
            "ALT": ["T"] * rows,
            "GT": genotypes,
            "FILTER": ["PASS"] * rows,
            "QUAL": [30.0] * rows,
        }

        df = pd.DataFrame(data)
        df.to_csv(file_path, sep="\t", index=False)

        return file_path

    @pytest.mark.slow
    def test_streaming_vs_chunked_small_file(self):
        """Compare performance on small files."""
        output_streaming = Path(self.test_dir) / "small_streaming.tsv"
        output_chunked = Path(self.test_dir) / "small_chunked.tsv"

        # Test streaming parallel
        start_time = time.time()
        process_streaming_parallel(self.small_file, output_streaming, self.config, threads=4)
        streaming_time = time.time() - start_time

        # Test chunked vectorized (if available)
        start_time = time.time()
        try:
            process_parallel_chunked_vectorized(
                self.small_file,
                output_chunked,
                self.config,
                chunk_size=2000,
                max_workers=4,
                available_memory_gb=16.0,
            )
            chunked_time = time.time() - start_time

            print("Small file performance:")
            print(f"  Streaming parallel: {streaming_time:.2f}s")
            print(f"  Chunked vectorized: {chunked_time:.2f}s")
            print(f"  Speedup ratio: {chunked_time / streaming_time:.2f}x")

        except ImportError:
            print(f"Small file - Streaming parallel: {streaming_time:.2f}s")
            print("Chunked vectorized not available for comparison")

    @pytest.mark.slow
    def test_streaming_vs_chunked_medium_file(self):
        """Compare performance on medium files."""
        output_streaming = Path(self.test_dir) / "medium_streaming.tsv"
        output_chunked = Path(self.test_dir) / "medium_chunked.tsv"

        # Test streaming parallel
        start_time = time.time()
        process_streaming_parallel(self.medium_file, output_streaming, self.config, threads=8)
        streaming_time = time.time() - start_time

        # Test chunked vectorized (if available)
        start_time = time.time()
        try:
            process_parallel_chunked_vectorized(
                self.medium_file,
                output_chunked,
                self.config,
                chunk_size=5000,
                max_workers=8,
                available_memory_gb=32.0,
            )
            chunked_time = time.time() - start_time

            print("Medium file performance:")
            print(f"  Streaming parallel: {streaming_time:.2f}s")
            print(f"  Chunked vectorized: {chunked_time:.2f}s")
            print(f"  Speedup ratio: {chunked_time / streaming_time:.2f}x")

        except ImportError:
            print(f"Medium file - Streaming parallel: {streaming_time:.2f}s")
            print("Chunked vectorized not available for comparison")

    @pytest.mark.slow
    def test_streaming_thread_scaling(self):
        """Test how streaming parallel scales with thread count."""
        thread_counts = [1, 2, 4, 8, 16]
        times = []

        for threads in thread_counts:
            output_file = Path(self.test_dir) / f"scaling_{threads}t.tsv"

            start_time = time.time()
            process_streaming_parallel(self.medium_file, output_file, self.config, threads=threads)
            elapsed = time.time() - start_time
            times.append(elapsed)

            print(f"Threads: {threads:2d}, Time: {elapsed:.2f}s")

        # Calculate scaling efficiency
        base_time = times[0]  # Single thread time
        print("\nScaling efficiency:")
        for i, (threads, time_taken) in enumerate(zip(thread_counts, times, strict=False)):
            if i == 0:
                efficiency = 100.0
            else:
                ideal_speedup = threads
                actual_speedup = base_time / time_taken
                efficiency = (actual_speedup / ideal_speedup) * 100
            print(f"  {threads:2d} threads: {efficiency:.1f}% efficiency")

    def test_memory_profiling_overhead(self):
        """Test overhead of memory profiling."""
        output_with_profiling = Path(self.test_dir) / "with_profiling.tsv"

        # This would test the memory profiling overhead
        # For now, just verify it doesn't crash
        start_time = time.time()
        process_streaming_parallel(self.small_file, output_with_profiling, self.config, threads=2)
        time_with_profiling = time.time() - start_time

        # In a real test, we'd disable profiling and compare
        print(f"Processing with memory profiling: {time_with_profiling:.2f}s")
        assert time_with_profiling > 0

    @pytest.mark.slow
    def test_chunk_size_impact(self):
        """Test impact of different chunk sizes (if configurable)."""
        # This test would evaluate optimal chunk sizes
        # For now, just test that processing completes
        output_file = Path(self.test_dir) / "chunk_test.tsv"

        start_time = time.time()
        process_streaming_parallel(self.medium_file, output_file, self.config, threads=4)
        elapsed = time.time() - start_time

        print(f"Chunk size optimization test completed in {elapsed:.2f}s")
        assert output_file.exists()

    def test_resource_utilization_metrics(self):
        """Test resource utilization during processing."""
        import threading

        import psutil

        # Monitor resource usage during processing
        cpu_samples = []
        memory_samples = []
        monitoring = True

        def monitor_resources():
            while monitoring:
                cpu_samples.append(psutil.cpu_percent(interval=0.1))
                memory_samples.append(psutil.virtual_memory().percent)
                time.sleep(0.5)

        monitor_thread = threading.Thread(target=monitor_resources)
        monitor_thread.start()

        # Process file
        output_file = Path(self.test_dir) / "resource_test.tsv"
        start_time = time.time()
        process_streaming_parallel(self.small_file, output_file, self.config, threads=4)
        elapsed = time.time() - start_time

        # Stop monitoring
        monitoring = False
        monitor_thread.join()

        if cpu_samples and memory_samples:
            avg_cpu = sum(cpu_samples) / len(cpu_samples)
            max_memory = max(memory_samples)

            print("Resource utilization:")
            print(f"  Processing time: {elapsed:.2f}s")
            print(f"  Average CPU usage: {avg_cpu:.1f}%")
            print(f"  Peak memory usage: {max_memory:.1f}%")

            # Basic assertions
            assert avg_cpu > 0  # Should be using CPU
            assert max_memory < 90  # Shouldn't max out memory


@pytest.mark.performance
@pytest.mark.slow
class TestMemoryEfficiency:
    """Test memory efficiency of the new implementation."""

    def test_memory_growth_patterns(self):
        """Test memory usage patterns during processing."""
        # This would monitor memory growth over time
        # For now, just create a basic test structure
        import psutil

        process = psutil.Process()
        initial_memory = process.memory_info().rss

        # Create test data
        test_dir = tempfile.mkdtemp()
        test_file = Path(test_dir) / "memory_test.tsv"

        # Generate moderate-sized test data
        data = pd.DataFrame({"GT": ["0/1,1/1,0/0"] * 5000, "other": ["data"] * 5000})
        data.to_csv(test_file, sep="\t", index=False)

        config = {
            "sample_list": "sample1,sample2,sample3",
            "separator": ";",
            "extract_fields_separator": ",",
        }

        output_file = Path(test_dir) / "memory_output.tsv"
        process_streaming_parallel(test_file, output_file, config, threads=2)

        final_memory = process.memory_info().rss
        memory_growth = final_memory - initial_memory

        print(f"Memory growth during processing: {memory_growth / 1024 / 1024:.1f} MB")

        # Clean up
        import shutil

        shutil.rmtree(test_dir, ignore_errors=True)

        # Memory growth should be reasonable (less than 1GB for this small test)
        assert memory_growth < 1024 * 1024 * 1024  # 1GB
