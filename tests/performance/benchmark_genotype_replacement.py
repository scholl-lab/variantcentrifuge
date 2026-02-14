"""
Genotype replacement benchmarks.

Tests both vectorized and sequential genotype replacement implementations
at multiple scales with realistic sample counts.
"""

from io import StringIO
from pathlib import Path

import pytest

from variantcentrifuge.replacer import replace_genotypes
from variantcentrifuge.vectorized_replacer import replace_genotypes_vectorized


def _build_replacer_config(sample_list: list[str]) -> dict:
    """
    Build minimal config dict for genotype replacement functions.

    Parameters
    ----------
    sample_list : list[str]
        List of sample IDs

    Returns
    -------
    dict
        Config dict with required keys for replacer functions
    """
    return {
        "sample_list": ",".join(sample_list),
        "extract_fields_separator": ":",
        "separator": ";",
        "genotype_replacement_map": {r"[2-9]": "1"},
        "append_extra_sample_fields": False,
        "extra_sample_fields": [],
        "extra_sample_field_delimiter": ":",
    }


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [100, 1000, 10000])
def test_vectorized_replacement_scaling(benchmark, synthetic_variants, tmp_path, n_variants):
    """
    Benchmark vectorized genotype replacement at multiple scales.

    The vectorized replacer operates on files (reads TSV, writes TSV).
    This benchmark measures the full file I/O + replacement pipeline.
    """
    # Use realistic sample count (genotype replacement scales with sample count)
    n_samples = 100
    sample_list = [f"SAMPLE_{i:04d}" for i in range(1, n_samples + 1)]

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)

    # Transform GT column from comma-separated to colon-separated
    # The vectorized replacer expects: "0/1:0/0:1/1:..." (one genotype per sample)
    df["GT"] = df["GT"].str.replace(",", ":")

    # Write to temporary file
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    df.to_csv(input_path, sep="\t", index=False)

    # Build config
    cfg = _build_replacer_config(sample_list)

    # Benchmark vectorized replacement
    def run_replacement():
        replace_genotypes_vectorized(input_path, output_path, cfg)

    if n_variants >= 10000:
        benchmark.pedantic(run_replacement, rounds=3, iterations=1, warmup_rounds=1)
    else:
        benchmark(run_replacement)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["component"] = "genotype_replacement_vectorized"

    # Verify output exists
    assert output_path.exists()
    assert output_path.stat().st_size > 0


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [100, 1000])
def test_sequential_replacement_scaling(benchmark, synthetic_variants, n_variants):
    """
    Benchmark sequential (original) genotype replacement.

    The sequential replacer operates on an iterator of lines. Skip 10K variants
    as the sequential implementation is significantly slower.
    """
    # Use realistic sample count
    n_samples = 100
    sample_list = [f"SAMPLE_{i:04d}" for i in range(1, n_samples + 1)]

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)

    # Convert DataFrame to TSV lines iterator
    tsv_string = df.to_csv(sep="\t", index=False)
    lines = iter(tsv_string.splitlines(keepends=True))

    # Build config
    cfg = _build_replacer_config(sample_list)

    # Benchmark sequential replacement (must consume iterator)
    def run_replacement():
        # Reset iterator for each benchmark iteration
        lines_iter = iter(tsv_string.splitlines(keepends=True))
        result = list(replace_genotypes(lines_iter, cfg))
        return result

    if n_variants >= 1000:
        result = benchmark.pedantic(run_replacement, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(run_replacement)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["component"] = "genotype_replacement_sequential"

    # Verify output
    assert len(result) == n_variants + 1  # +1 for header


@pytest.mark.performance
@pytest.mark.parametrize("n_samples", [10, 100, 500])
def test_replacement_sample_scaling(benchmark, synthetic_variants, tmp_path, n_samples):
    """
    Benchmark genotype replacement with varying sample counts.

    Genotype replacement performance is sensitive to both variant count AND
    sample count. This test fixes variant count and varies sample count.
    """
    n_variants = 1000
    sample_list = [f"SAMPLE_{i:04d}" for i in range(1, n_samples + 1)]

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)

    # Transform GT column to colon-separated format
    df["GT"] = df["GT"].str.replace(",", ":")

    # Write to temporary file
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    df.to_csv(input_path, sep="\t", index=False)

    # Build config
    cfg = _build_replacer_config(sample_list)

    # Benchmark vectorized replacement
    def run_replacement():
        replace_genotypes_vectorized(input_path, output_path, cfg)

    if n_samples >= 500:
        benchmark.pedantic(run_replacement, rounds=3, iterations=1, warmup_rounds=1)
    else:
        benchmark(run_replacement)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["component"] = "genotype_replacement_sample_scaling"

    # Verify output
    assert output_path.exists()
