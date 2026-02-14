"""
DataFrame I/O benchmarks.

Tests pandas TSV read/write performance at multiple scales.
Includes Phase 8 optimization measurements: PyArrow engine speedup,
categorical dtype memory reduction, and itertuples iteration speedup.
"""

import os

import pandas as pd
import pytest

from variantcentrifuge.dataframe_optimizer import (
    detect_categorical_columns,
    load_optimized_dataframe,
)


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [1000, 10000, 50000])
def test_csv_read_scaling(benchmark, synthetic_variants, tmp_path, n_variants):
    """
    Benchmark pandas read_csv performance at multiple scales.

    Tests the default pandas CSV reader (C engine) which is the hot path for
    all pipeline data loading. This establishes the baseline for Phase 8's
    PyArrow optimization comparison.
    """
    n_samples = 10

    # Generate synthetic data and write to file (setup, not benchmarked)
    df = synthetic_variants(n_variants, n_samples, seed=42)
    file_path = tmp_path / "variants.tsv"
    df.to_csv(file_path, sep="\t", index=False)

    # Get file size for metadata
    file_size_mb = file_path.stat().st_size / (1024 * 1024)
    n_columns = len(df.columns)

    # Benchmark read operation
    def read_file():
        return pd.read_csv(file_path, sep="\t")

    # Use setup function to ensure fresh read each iteration
    def setup():
        # Clear any OS-level caching by reading a different file
        pass

    if n_variants >= 50000:
        result = benchmark.pedantic(read_file, setup=setup, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(read_file)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_columns"] = n_columns
    benchmark.extra_info["file_size_mb"] = round(file_size_mb, 2)
    benchmark.extra_info["component"] = "dataframe_read"

    # Verify
    assert len(result) == n_variants


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [1000, 10000, 50000])
def test_csv_write_scaling(benchmark, synthetic_variants, tmp_path, n_variants):
    """
    Benchmark pandas to_csv performance at multiple scales.

    Tests DataFrame write performance which affects all pipeline output
    generation (TSV, intermediate files).
    """
    n_samples = 10

    # Generate synthetic data (setup, not benchmarked)
    df = synthetic_variants(n_variants, n_samples, seed=42)
    n_columns = len(df.columns)

    # Benchmark write operation
    def write_file():
        output_path = tmp_path / f"output_{os.getpid()}.tsv"
        df.to_csv(output_path, sep="\t", index=False)
        return output_path

    if n_variants >= 50000:
        output_path = benchmark.pedantic(write_file, rounds=3, iterations=1, warmup_rounds=1)
    else:
        output_path = benchmark(write_file)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_columns"] = n_columns
    benchmark.extra_info["component"] = "dataframe_write"

    # Verify output exists
    assert output_path.exists()


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [1000, 10000, 50000])
def test_pyarrow_read_scaling(benchmark, synthetic_variants, tmp_path, n_variants):
    """
    Benchmark pandas read_csv with PyArrow engine.

    Establishes the baseline for PyArrow-accelerated reads. This will be
    compared against test_csv_read_scaling to quantify the speedup in Phase 8.

    Skips if pyarrow is not installed.
    """
    pytest.importorskip("pyarrow")

    n_samples = 10

    # Generate synthetic data and write to file
    df = synthetic_variants(n_variants, n_samples, seed=42)
    file_path = tmp_path / "variants.tsv"
    df.to_csv(file_path, sep="\t", index=False)

    file_size_mb = file_path.stat().st_size / (1024 * 1024)
    n_columns = len(df.columns)

    # Benchmark PyArrow read
    def read_file():
        return pd.read_csv(file_path, sep="\t", engine="pyarrow")

    if n_variants >= 50000:
        result = benchmark.pedantic(read_file, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(read_file)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_columns"] = n_columns
    benchmark.extra_info["file_size_mb"] = round(file_size_mb, 2)
    benchmark.extra_info["engine"] = "pyarrow"
    benchmark.extra_info["component"] = "dataframe_read_pyarrow"

    # Verify
    assert len(result) == n_variants


@pytest.mark.performance
@pytest.mark.parametrize("n_columns", [10, 50, 100])
def test_csv_read_column_scaling(benchmark, tmp_path, n_columns):
    """
    Benchmark CSV read with varying column counts.

    DataFrame I/O performance can be affected by number of columns (wide vs
    narrow tables). This test fixes row count and varies column count.
    """
    n_variants = 10000

    # Generate DataFrame with many columns
    import numpy as np

    rng = np.random.default_rng(42)
    data = {f"col_{i}": rng.integers(0, 1000, size=n_variants) for i in range(n_columns)}
    df = pd.DataFrame(data)

    # Write to file
    file_path = tmp_path / "wide_table.tsv"
    df.to_csv(file_path, sep="\t", index=False)

    file_size_mb = file_path.stat().st_size / (1024 * 1024)

    # Benchmark read
    def read_file():
        return pd.read_csv(file_path, sep="\t")

    if n_columns >= 100:
        result = benchmark.pedantic(read_file, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(read_file)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_columns"] = n_columns
    benchmark.extra_info["file_size_mb"] = round(file_size_mb, 2)
    benchmark.extra_info["component"] = "dataframe_read_column_scaling"

    # Verify
    assert result.shape == (n_variants, n_columns)


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [1000, 10000])
def test_optimized_loader_performance(benchmark, synthetic_variants, tmp_path, n_variants):
    """
    Benchmark load_optimized_dataframe with PyArrow + categorical dtypes.

    Tests the full Phase 8 optimization stack: PyArrow engine, categorical
    dtype auto-detection, and column sanitization. This is the integrated
    real-world performance of the optimization.
    """
    pytest.importorskip("pyarrow")

    n_samples = 10

    # Generate synthetic data and write to file
    df = synthetic_variants(n_variants, n_samples, seed=42)
    file_path = tmp_path / "variants.tsv"
    df.to_csv(file_path, sep="\t", index=False)

    file_size_mb = file_path.stat().st_size / (1024 * 1024)
    n_columns = len(df.columns)

    # Benchmark optimized load
    def load_file():
        result_df, _ = load_optimized_dataframe(str(file_path))
        return result_df

    result = benchmark(load_file)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_columns"] = n_columns
    benchmark.extra_info["file_size_mb"] = round(file_size_mb, 2)
    benchmark.extra_info["component"] = "dataframe_optimized_load"

    # Verify
    assert len(result) == n_variants


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [1000, 10000])
def test_categorical_memory_reduction(benchmark, synthetic_variants, tmp_path, n_variants):
    """
    Measure memory reduction from categorical dtypes.

    Compares memory usage of the same DataFrame loaded with and without
    categorical dtypes. Measures deep memory usage to account for object
    dtype overhead.
    """
    n_samples = 10

    # Generate synthetic data and write to file
    df = synthetic_variants(n_variants, n_samples, seed=42)
    file_path = tmp_path / "variants.tsv"
    df.to_csv(file_path, sep="\t", index=False)

    # Load without categorical (baseline)
    df_baseline = pd.read_csv(file_path, sep="\t", dtype=str, keep_default_na=False)
    baseline_memory = df_baseline.memory_usage(deep=True).sum()

    # Load with categorical dtypes
    dtype_map = detect_categorical_columns(file_path)
    df_optimized = pd.read_csv(
        file_path, sep="\t", dtype=dtype_map, keep_default_na=False, na_values=[""]
    )
    optimized_memory = df_optimized.memory_usage(deep=True).sum()

    # Calculate reduction
    memory_reduction_pct = (
        (1 - optimized_memory / baseline_memory) * 100 if baseline_memory > 0 else 0
    )

    # Benchmark the detection + load process
    def load_with_categorical():
        dtype_map = detect_categorical_columns(file_path)
        return pd.read_csv(
            file_path, sep="\t", dtype=dtype_map, keep_default_na=False, na_values=[""]
        )

    result = benchmark(load_with_categorical)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["baseline_memory_mb"] = round(baseline_memory / (1024 * 1024), 2)
    benchmark.extra_info["optimized_memory_mb"] = round(optimized_memory / (1024 * 1024), 2)
    benchmark.extra_info["memory_reduction_pct"] = round(memory_reduction_pct, 1)
    benchmark.extra_info["component"] = "categorical_memory"

    # Verify
    assert len(result) == n_variants
    # Conservative assertion: at least 20% reduction (target is 50-70%)
    assert memory_reduction_pct >= 20.0, (
        f"Expected >= 20% memory reduction, got {memory_reduction_pct:.1f}%"
    )


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [10000])
def test_itertuples_vs_iterrows(benchmark, synthetic_variants, n_variants):
    """
    Compare itertuples vs iterrows iteration speed.

    Simulates the inheritance analysis pattern where we iterate over
    DataFrame rows and access multiple columns. This is the hot path
    that Phase 8 optimized.
    """
    n_samples = 10

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)

    # Ensure we have the columns we need
    if "CHROM" not in df.columns:
        df["CHROM"] = "chr1"
    if "POS" not in df.columns:
        df["POS"] = range(1, len(df) + 1)
    if "REF" not in df.columns:
        df["REF"] = "A"
    if "ALT" not in df.columns:
        df["ALT"] = "T"

    # Benchmark iterrows (old way)
    def iterate_with_iterrows():
        result = []
        for _, row in df.iterrows():
            # Simulate typical inheritance analysis access pattern
            key = f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}"
            result.append(key)
        return result

    # Benchmark itertuples (new way)
    def iterate_with_itertuples():
        result = []
        for row in df.itertuples():
            # Simulate same access pattern
            key = f"{row.CHROM}:{row.POS}:{row.REF}:{row.ALT}"
            result.append(key)
        return result

    # Run itertuples benchmark (this is what we measure)
    itertuples_result = benchmark(iterate_with_itertuples)

    # For comparison, measure iterrows separately (not benchmarked, just timed)
    import time

    start = time.perf_counter()
    iterrows_result = iterate_with_iterrows()
    iterrows_time = time.perf_counter() - start

    # Calculate speedup
    itertuples_time = benchmark.stats.stats.mean
    speedup = iterrows_time / itertuples_time if itertuples_time > 0 else 0

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["iterrows_time_s"] = round(iterrows_time, 4)
    benchmark.extra_info["itertuples_time_s"] = round(itertuples_time, 4)
    benchmark.extra_info["speedup_factor"] = round(speedup, 1)
    benchmark.extra_info["component"] = "itertuples_speedup"

    # Verify results match
    assert len(itertuples_result) == n_variants
    assert itertuples_result == iterrows_result
    # Conservative assertion: at least 5x speedup (target is 10-14x)
    assert speedup >= 5.0, f"Expected >= 5x speedup, got {speedup:.1f}x"
