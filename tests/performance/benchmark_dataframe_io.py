"""
DataFrame I/O benchmarks.

Tests pandas TSV read/write performance at multiple scales. Establishes baseline
for future PyArrow optimization (Phase 8).
"""

import os

import pandas as pd
import pytest


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
