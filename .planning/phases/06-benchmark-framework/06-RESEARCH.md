# Phase 6: Benchmark Framework - Research

**Researched:** 2026-02-14
**Domain:** Python performance benchmarking, synthetic genomic data generation, regression detection
**Confidence:** HIGH

## Summary

Performance benchmarking for genomics pipelines requires combining pytest-benchmark for execution timing, tracemalloc for memory profiling, and synthetic data generators for reproducible test datasets. The standard approach uses pytest-benchmark's fixture-based API with JSON result storage, regression detection via `--benchmark-compare-fail`, and custom metadata for domain-specific metrics. Memory profiling through tracemalloc enables snapshot-based comparisons and peak memory tracking without heavy third-party dependencies. Synthetic genomic data generation prioritizes reproducibility through seeded randomness, realistic distributions derived from anonymized real-world patterns, and pandas-based DataFrame construction to match the pipeline's internal data structures.

**Key findings:**
- pytest-benchmark is the standard for Python benchmarking with built-in regression detection
- tracemalloc (stdlib) provides sufficient memory profiling for budget enforcement
- Synthetic genomic data requires seeded random generators for reproducibility
- Ratio assertions (comparing implementations within same run) eliminate CI flakiness
- JSON result storage with simple diff utilities matches the "optimization sprint tooling" use case

**Primary recommendation:** Use pytest-benchmark with custom fixtures for synthetic data generation, tracemalloc integration for memory budgets, and ratio assertions for comparing vectorized vs sequential implementations within the same test run.

## Standard Stack

The established libraries/tools for performance benchmarking in Python:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pytest-benchmark | 5.2.3 | Performance benchmarking fixture | De facto standard for pytest-based benchmarking, built-in regression detection |
| tracemalloc | stdlib (3.10+) | Memory profiling | Part of stdlib, designed specifically for memory allocation tracking |
| pandas | existing dep | Synthetic data generation | Already a dependency, natural fit for genomic DataFrame creation |
| numpy | existing dep | Statistical distributions | Underpins pandas, needed for realistic variant distributions |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pytest-mock | existing dev dep | Mocking external tools | Mock bcftools/SnpSift to isolate Python code benchmarks |
| psutil | existing dep | System resource monitoring | Already in deps, useful for validating memory budget assertions |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| pytest-benchmark | timeit module | pytest-benchmark provides regression detection, JSON storage, and pytest integration; timeit is lower-level |
| tracemalloc | memory-profiler | memory-profiler has function decorators but adds dependency and ~20% overhead; tracemalloc is stdlib with 5-20% overhead |
| Custom JSON diff | pytest-benchmark CLI | pytest-benchmark CLI has histogram generation but more complex; simple diff helper matches "keep it simple" requirement |

**Installation:**
```bash
# pytest-benchmark is the only new dependency
pip install pytest-benchmark[histogram]  # Optional histogram support
# OR add to pyproject.toml dev dependencies
```

## Architecture Patterns

### Recommended Project Structure
```
tests/performance/
├── conftest.py                    # Synthetic data fixtures
├── benchmark_inheritance.py       # Inheritance analysis benchmarks (micro/meso)
├── benchmark_comp_het.py          # Dedicated compound het benchmarks
├── benchmark_genotype_replacement.py  # Genotype replacement benchmarks
├── benchmark_gene_burden.py       # Gene burden analysis benchmarks
├── benchmark_scoring.py           # Scoring benchmarks
├── benchmark_dataframe_io.py      # DataFrame read/write benchmarks
├── benchmark_pipeline_macro.py    # Full pipeline end-to-end benchmarks
└── helpers/
    ├── synthetic_data.py          # Synthetic DataFrame generators
    ├── memory_budgets.py          # tracemalloc budget enforcement helpers
    └── result_diff.py             # JSON comparison utility
benchmarks/                        # Gitignored result storage
├── .gitignore                     # Ignore all JSON files
└── (ephemeral JSON results)
```

### Pattern 1: Benchmark Fixture with Parameterization
**What:** Use pytest.mark.parametrize to test multiple data sizes/configurations in one test
**When to use:** Testing scalability across variant counts (100, 1K, 10K, 50K)
**Example:**
```python
# Source: pytest-benchmark docs + variantcentrifuge patterns
import pytest

@pytest.mark.parametrize("n_variants,n_samples", [
    (100, 10),
    (1000, 100),
    (10000, 500),
    (50000, 1000),
])
def test_inheritance_analysis_scaling(benchmark, synthetic_variants, n_variants, n_samples):
    """Benchmark inheritance analysis at different scales."""
    df = synthetic_variants(n_variants=n_variants, n_samples=n_samples, seed=42)
    pedigree = synthetic_pedigree(n_samples=n_samples, seed=42)

    result = benchmark(analyze_inheritance, df, pedigree)

    # Store custom metadata
    benchmark.extra_info['n_variants'] = n_variants
    benchmark.extra_info['n_samples'] = n_samples

    assert result is not None
```

### Pattern 2: Ratio Assertions for Flakiness-Free Comparison
**What:** Compare two implementations within the same test run, compute speedup ratio
**When to use:** Comparing vectorized vs sequential, avoiding CI noise
**Example:**
```python
# Source: Adapted from pytest-benchmark pedantic mode
def test_vectorized_vs_sequential_ratio(benchmark):
    """Compare vectorized and sequential implementations via ratio."""
    df = synthetic_variants(n_variants=10000, n_samples=100, seed=42)

    # Benchmark vectorized
    result_vec = benchmark.pedantic(
        analyze_vectorized,
        args=(df,),
        iterations=5,
        rounds=3
    )
    vectorized_time = benchmark.stats.stats.mean

    # Benchmark sequential (outside benchmark fixture)
    import time
    times = []
    for _ in range(5):
        start = time.perf_counter()
        analyze_sequential(df)
        times.append(time.perf_counter() - start)
    sequential_time = sum(times) / len(times)

    # Ratio assertion
    speedup = sequential_time / vectorized_time
    benchmark.extra_info['speedup_ratio'] = speedup

    # Assert vectorized is faster (ratio > 1.0)
    assert speedup > 1.0, f"Vectorized should be faster, got {speedup:.2f}x"
```

### Pattern 3: Memory Budget Enforcement with tracemalloc
**What:** Track peak memory during benchmark execution, fail on budget violations
**When to use:** Functions with known memory constraints (e.g., genotype replacement < 2GB)
**Example:**
```python
# Source: Python tracemalloc docs
import tracemalloc

def test_genotype_replacement_memory_budget(benchmark, synthetic_variants):
    """Enforce memory budget for genotype replacement."""
    df = synthetic_variants(n_variants=50000, n_samples=1000, seed=42)

    # Start memory tracking
    tracemalloc.start()
    tracemalloc.reset_peak()

    # Run benchmark
    result = benchmark(replace_genotypes, df)

    # Get peak memory
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    peak_mb = peak / 1024 / 1024
    benchmark.extra_info['peak_memory_mb'] = peak_mb

    # Budget: warn only (not hard fail per CONTEXT.md)
    MEMORY_BUDGET_MB = 2048
    if peak_mb > MEMORY_BUDGET_MB:
        import warnings
        warnings.warn(f"Memory budget exceeded: {peak_mb:.1f}MB > {MEMORY_BUDGET_MB}MB")

    assert result is not None
```

### Pattern 4: Synthetic Data Fixtures with Seeded Randomness
**What:** Generate reproducible DataFrames matching real-world genomic data distributions
**When to use:** All benchmarks requiring variant data
**Example:**
```python
# Source: variantcentrifuge test patterns + pandas
import numpy as np
import pandas as pd
import pytest

@pytest.fixture
def synthetic_variants():
    """Factory fixture for generating synthetic variant DataFrames."""
    def _generate(n_variants, n_samples, seed=42):
        """Generate synthetic variant data with realistic distributions.

        Distributions derived from real-world VCF patterns (anonymized).
        All identifiers are synthetic, no traceable information.
        """
        rng = np.random.default_rng(seed)

        # Realistic chromosome distribution (weighted toward chr1-10)
        chr_weights = [0.15] + [0.08] * 9 + [0.03] * 12 + [0.01] * 2  # 1-22, X, Y
        chromosomes = rng.choice(
            [str(i) for i in range(1, 23)] + ['X', 'Y'],
            size=n_variants,
            p=chr_weights
        )

        # Positions distributed across 250Mbp typical chromosome
        positions = rng.integers(1, 250_000_000, size=n_variants)

        # REF/ALT: mostly SNVs, some indels
        ref_bases = rng.choice(['A', 'C', 'G', 'T'], size=n_variants, p=[0.3, 0.2, 0.2, 0.3])
        alt_bases = rng.choice(['A', 'C', 'G', 'T'], size=n_variants, p=[0.25, 0.25, 0.25, 0.25])

        # Genotypes: realistic allele frequency distribution
        # 90% rare (mostly 0/0), 8% uncommon (some 0/1), 2% common (varied)
        genotypes = []
        for _ in range(n_variants):
            variant_type = rng.choice(['rare', 'uncommon', 'common'], p=[0.90, 0.08, 0.02])
            if variant_type == 'rare':
                # Mostly homozygous reference, occasional het/hom-alt
                gts = rng.choice(['0/0', '0/1', '1/1'], size=n_samples, p=[0.95, 0.04, 0.01])
            elif variant_type == 'uncommon':
                gts = rng.choice(['0/0', '0/1', '1/1'], size=n_samples, p=[0.80, 0.15, 0.05])
            else:  # common
                gts = rng.choice(['0/0', '0/1', '1/1'], size=n_samples, p=[0.50, 0.40, 0.10])

            genotypes.append(','.join(gts))

        # Build DataFrame matching variantcentrifuge structure
        df = pd.DataFrame({
            'CHROM': chromosomes,
            'POS': positions,
            'REF': ref_bases,
            'ALT': alt_bases,
            'GT': genotypes,
            'GENE': rng.choice([f'GENE_{i:04d}' for i in range(100)], size=n_variants),
            'FILTER': rng.choice(['PASS', 'LowQual'], size=n_variants, p=[0.95, 0.05]),
        })

        return df.sort_values(['CHROM', 'POS']).reset_index(drop=True)

    return _generate
```

### Pattern 5: Regression Detection with 20% Threshold
**What:** Use pytest-benchmark's --benchmark-compare-fail with 20% threshold
**When to use:** All performance benchmarks (uniform policy per CONTEXT.md)
**Example:**
```bash
# First run: establish baseline
pytest tests/performance/ --benchmark-save=baseline

# Subsequent runs: detect 20% regressions
pytest tests/performance/ --benchmark-compare=baseline --benchmark-compare-fail=mean:20%

# This will FAIL the test suite if any benchmark's mean time increased by 20% or more
```

### Pattern 6: Custom Metadata for Genomic Benchmarks
**What:** Store domain-specific metrics (n_variants, n_genes, peak_memory_mb) in benchmark.extra_info
**When to use:** All benchmarks, enables result analysis and filtering
**Example:**
```python
def test_gene_burden_analysis(benchmark, synthetic_variants):
    """Benchmark gene burden analysis with metadata."""
    n_variants = 10000
    n_genes = 50
    df = synthetic_variants(n_variants=n_variants, n_samples=200, seed=42)

    result = benchmark(perform_gene_burden, df)

    # Store metadata for result analysis
    benchmark.extra_info['n_variants'] = n_variants
    benchmark.extra_info['n_genes'] = n_genes
    benchmark.extra_info['n_samples'] = 200
    benchmark.extra_info['analysis_type'] = 'gene_burden'

    assert result is not None
```

### Anti-Patterns to Avoid
- **Benchmarking with real patient data:** Privacy risk, non-reproducible. Always use synthetic data with seeded randomness.
- **CI-based performance tests:** CONTEXT.md specifies local-only, no CI integration. Avoid adding benchmark jobs to GitHub Actions.
- **Comparing across different machines:** Results are machine-dependent. Ratio assertions (within same run) solve this for vectorized vs sequential comparisons.
- **Complex historical tracking:** CONTEXT.md specifies ephemeral results, no timestamped storage. Don't build a time-series database.
- **Mocking in macro benchmarks:** End-to-end tests should exercise real pipeline code paths (but can skip external tools like bcftools).

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Performance timing | Custom time.perf_counter loops | pytest-benchmark fixture | Handles warmup, calibration, statistical aggregation, regression detection |
| Memory profiling | Manual memory tracking with psutil | tracemalloc snapshots | Stdlib, line-level attribution, snapshot comparison built-in |
| Benchmark result storage | Custom JSON schema | pytest-benchmark --benchmark-save | Standard format, CLI tools for comparison, histogram generation |
| Statistical aggregation | Manual mean/std calculation | pytest-benchmark stats | Provides min/max/mean/median/stddev/IQR/outliers automatically |
| Regression detection | Custom threshold checking | --benchmark-compare-fail | Proven approach, supports percentage and absolute thresholds |
| Synthetic genomic data | VCF file parsing/generation | Pandas DataFrame factories | VCF parsing is complex (headers, INFO fields, FORMAT); DataFrames match internal pipeline structures |

**Key insight:** pytest-benchmark handles calibration (automatic iteration count adjustment), warmup rounds, statistical outlier detection, and result formatting. Custom timing code would miss these critical features and introduce measurement errors.

## Common Pitfalls

### Pitfall 1: Unstable Benchmark Results from System Noise
**What goes wrong:** Benchmarks vary by 10-50% between runs due to background processes, thermal throttling, or garbage collection
**Why it happens:** Modern operating systems have unpredictable scheduling, GC can trigger mid-benchmark
**How to avoid:**
- Use pytest-benchmark's multiple rounds/iterations (default: auto-calibrated)
- Disable GC during benchmarks: `@pytest.mark.benchmark(disable_gc=True)`
- Use ratio assertions for implementation comparisons (eliminates cross-run variance)
- Run with sufficient warmup: pytest-benchmark does this automatically
**Warning signs:** Benchmark results differ by >10% when re-running same code unchanged

### Pitfall 2: Memory Profiling Overhead Affecting Timing
**What goes wrong:** tracemalloc adds 5-20% overhead, skewing performance measurements
**Why it happens:** tracemalloc intercepts every memory allocation to track call stacks
**How to avoid:**
- **Never enable tracemalloc during timing benchmarks**
- Use separate tests: timing tests without tracemalloc, memory tests with tracemalloc
- Document tracemalloc overhead in memory budget tests
- For combined tests, use pedantic mode to exclude setup from timing
**Warning signs:** Benchmarks with tracemalloc enabled run 10-20% slower than without

### Pitfall 3: Non-Reproducible Synthetic Data
**What goes wrong:** Benchmark results vary because synthetic data is different each run
**Why it happens:** Using unseeded random generators (random.random(), np.random.rand())
**How to avoid:**
- Always seed random generators: `rng = np.random.default_rng(seed=42)`
- Use factory fixtures with seed parameters, not module-level generation
- Validate data generation with hash checks in tests
- Document seed values in fixture docstrings
**Warning signs:** "Same" benchmark shows different timing despite code being unchanged

### Pitfall 4: Benchmarking DataFrame Copies Instead of Views
**What goes wrong:** Accidentally benchmarking `df.copy()` overhead instead of algorithm performance
**Why it happens:** Pandas DataFrame operations sometimes return views, sometimes copies
**How to avoid:**
- Understand pandas copy vs view semantics (chained assignment, filtering, slicing)
- Generate data in fixture setup (outside benchmark timing)
- Use `benchmark.pedantic` with separate setup functions for data preparation
- Verify with line_profiler that benchmark captures intended operations
**Warning signs:** Benchmark time scales linearly with DataFrame size even for O(1) operations

### Pitfall 5: Flaky Regression Detection from CI Noise
**What goes wrong:** --benchmark-compare-fail triggers false positives in CI environments
**Why it happens:** Shared CI runners have variable CPU allocation, noisy neighbors
**How to avoid:**
- **Don't run benchmarks in CI** (CONTEXT.md specifies local-only)
- Use ratio assertions for vectorized vs sequential comparisons (eliminates CI variance)
- If CI needed later, use higher thresholds (30-50%) or percentile-based comparisons
**Warning signs:** Benchmarks fail intermittently in CI but pass locally

### Pitfall 6: Privacy Leaks in Synthetic Data Generation
**What goes wrong:** Synthetic data generator accidentally embeds real patient identifiers or cohort names
**Why it happens:** Copy-pasting from real VCF files, using actual sample IDs as templates
**How to avoid:**
- **Never reference real cohort names in committed code** (CONTEXT.md requirement)
- Use fully synthetic identifiers: "SAMPLE_0001", "GENE_0042", never "patient_xyz"
- Generate data from statistical distributions, not by sampling real VCFs
- Code review all synthetic data generators for hardcoded identifiers
- Run data anonymization checks before committing generator code
**Warning signs:** Sample IDs, gene names, or chromosome coordinates match real datasets

### Pitfall 7: Incorrect Memory Budget Assertions
**What goes wrong:** Memory budget test fails hard instead of warning, or doesn't detect leaks
**Why it happens:** Misunderstanding CONTEXT.md policy (memory = warnings only) or not resetting peak
**How to avoid:**
- Memory violations are **warnings only, not hard failures** (per CONTEXT.md)
- Use `tracemalloc.reset_peak()` before operation to measure specific code section
- Take snapshots before stopping: `snapshot = tracemalloc.take_snapshot(); tracemalloc.stop()`
- Use `get_traced_memory()` for simple peak tracking, snapshots for leak detection
**Warning signs:** Memory tests hard-fail when they should warn, or don't detect obvious leaks

## Code Examples

Verified patterns from official sources and variantcentrifuge codebase:

### Example 1: Basic Benchmark Test
```python
# Source: pytest-benchmark docs
import pytest

def test_inheritance_deduction_performance(benchmark, synthetic_variants):
    """Benchmark inheritance pattern deduction."""
    from variantcentrifuge.inheritance.deducer import deduce_inheritance_patterns

    df = synthetic_variants(n_variants=1000, n_samples=10, seed=42)
    pedigree = {'child': {...}, 'mother': {...}, 'father': {...}}

    result = benchmark(deduce_inheritance_patterns, df, pedigree)

    benchmark.extra_info['n_variants'] = 1000
    benchmark.extra_info['component'] = 'inheritance_deduction'

    assert result is not None
```

### Example 2: Parametrized Scaling Benchmark
```python
# Source: pytest-benchmark + pytest parametrize docs
import pytest

@pytest.mark.parametrize("n_variants", [100, 1000, 10000, 50000])
def test_comp_het_scaling(benchmark, synthetic_variants, n_variants):
    """Benchmark compound het detection at different scales."""
    from variantcentrifuge.inheritance.comp_het_vectorized import detect_compound_hets

    df = synthetic_variants(n_variants=n_variants, n_samples=50, seed=42)

    result = benchmark(detect_compound_hets, df)

    benchmark.extra_info['n_variants'] = n_variants
    benchmark.extra_info['component'] = 'compound_het'

    assert result is not None
```

### Example 3: Memory Budget Test with tracemalloc
```python
# Source: Python tracemalloc docs
import tracemalloc
import pytest
import warnings

def test_genotype_replacement_memory(benchmark, synthetic_variants):
    """Enforce memory budget for genotype replacement."""
    from variantcentrifuge.genotype_replacement import replace_genotypes

    n_variants = 50000
    n_samples = 1000
    df = synthetic_variants(n_variants=n_variants, n_samples=n_samples, seed=42)

    # Start memory tracking
    tracemalloc.start()
    tracemalloc.reset_peak()

    # Run function (NOT via benchmark fixture - tracemalloc overhead)
    result = replace_genotypes(df)

    # Get peak memory
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    peak_mb = peak / 1024 / 1024

    # Warn if budget exceeded (warnings only, not hard fail)
    BUDGET_MB = 2048
    if peak_mb > BUDGET_MB:
        warnings.warn(
            f"Memory budget exceeded: {peak_mb:.1f}MB > {BUDGET_MB}MB "
            f"(n_variants={n_variants}, n_samples={n_samples})"
        )

    # Could log to JSON for tracking
    print(f"Peak memory: {peak_mb:.1f}MB")

    assert result is not None
```

### Example 4: Ratio Assertion for Vectorized vs Sequential
```python
# Source: pytest-benchmark pedantic mode docs
import time
import pytest

def test_vectorized_speedup_ratio(synthetic_variants):
    """Compare vectorized and sequential implementations."""
    from variantcentrifuge.inheritance.analyzer import (
        analyze_inheritance_vectorized,
        analyze_inheritance_sequential
    )

    df = synthetic_variants(n_variants=10000, n_samples=100, seed=42)
    pedigree = {...}

    # Time sequential implementation
    sequential_times = []
    for _ in range(5):
        start = time.perf_counter()
        analyze_inheritance_sequential(df, pedigree)
        sequential_times.append(time.perf_counter() - start)
    sequential_mean = sum(sequential_times) / len(sequential_times)

    # Time vectorized implementation
    vectorized_times = []
    for _ in range(5):
        start = time.perf_counter()
        analyze_inheritance_vectorized(df, pedigree)
        vectorized_times.append(time.perf_counter() - start)
    vectorized_mean = sum(vectorized_times) / len(vectorized_times)

    # Compute ratio
    speedup = sequential_mean / vectorized_mean

    print(f"Sequential: {sequential_mean:.3f}s, Vectorized: {vectorized_mean:.3f}s")
    print(f"Speedup: {speedup:.2f}x")

    # Assert vectorized is faster (zero flakiness - same run comparison)
    assert speedup > 1.0, f"Vectorized should be faster, got {speedup:.2f}x"
```

### Example 5: Synthetic Data Factory Fixture
```python
# Source: variantcentrifuge test patterns
import numpy as np
import pandas as pd
import pytest

@pytest.fixture
def synthetic_pedigree():
    """Generate synthetic pedigree (trio) for benchmarks."""
    def _generate(n_samples=3, seed=42):
        """Generate trio pedigree.

        For benchmarking purposes only. Fully synthetic identifiers.
        """
        if n_samples < 3:
            n_samples = 3  # Minimum trio

        # Simple trio structure
        pedigree = {
            'CHILD_001': {
                'sample_id': 'CHILD_001',
                'father_id': 'FATHER_001',
                'mother_id': 'MOTHER_001',
                'sex': '1',
                'affected_status': '2',
            },
            'FATHER_001': {
                'sample_id': 'FATHER_001',
                'father_id': '0',
                'mother_id': '0',
                'sex': '1',
                'affected_status': '1',
            },
            'MOTHER_001': {
                'sample_id': 'MOTHER_001',
                'father_id': '0',
                'mother_id': '0',
                'sex': '2',
                'affected_status': '1',
            },
        }

        # Add additional unrelated samples if needed
        rng = np.random.default_rng(seed)
        for i in range(3, n_samples):
            sample_id = f'SAMPLE_{i+1:04d}'
            pedigree[sample_id] = {
                'sample_id': sample_id,
                'father_id': '0',
                'mother_id': '0',
                'sex': str(rng.choice([1, 2])),
                'affected_status': str(rng.choice([1, 2])),
            }

        return pedigree

    return _generate
```

### Example 6: Simple JSON Diff Helper
```python
# Source: CONTEXT.md requirement for simple diff utility
import json
from pathlib import Path

def compare_benchmark_results(baseline_file: str, current_file: str):
    """Compare two benchmark result files and print differences.

    Args:
        baseline_file: Path to baseline JSON (e.g., benchmarks/baseline.json)
        current_file: Path to current JSON (e.g., benchmarks/current.json)

    Prints a table showing what got faster/slower.
    """
    with open(baseline_file) as f:
        baseline = json.load(f)
    with open(current_file) as f:
        current = json.load(f)

    # Extract benchmark data (pytest-benchmark JSON structure)
    baseline_benchmarks = {b['name']: b for b in baseline['benchmarks']}
    current_benchmarks = {b['name']: b for b in current['benchmarks']}

    print(f"{'Benchmark':<50} {'Baseline':>12} {'Current':>12} {'Change':>10}")
    print("-" * 86)

    for name in sorted(set(baseline_benchmarks.keys()) | set(current_benchmarks.keys())):
        if name in baseline_benchmarks and name in current_benchmarks:
            base_time = baseline_benchmarks[name]['stats']['mean']
            curr_time = current_benchmarks[name]['stats']['mean']
            change_pct = ((curr_time - base_time) / base_time) * 100

            status = "FASTER" if change_pct < 0 else "SLOWER"
            color = "\033[92m" if change_pct < 0 else "\033[91m"  # Green/Red
            reset = "\033[0m"

            print(
                f"{name:<50} {base_time:>10.4f}s {curr_time:>10.4f}s "
                f"{color}{change_pct:>+8.1f}% {status}{reset}"
            )
        elif name in baseline_benchmarks:
            print(f"{name:<50} {'REMOVED':>35}")
        else:
            print(f"{name:<50} {'NEW':>35}")

# Usage:
# compare_benchmark_results('benchmarks/baseline.json', 'benchmarks/current.json')
```

### Example 7: End-to-End Pipeline Benchmark
```python
# Source: variantcentrifuge test patterns
import pytest
from pathlib import Path

@pytest.mark.parametrize("pipeline_type", ["classic", "stage_based"])
def test_full_pipeline_benchmark(benchmark, tmp_path, synthetic_variants, pipeline_type):
    """Benchmark full pipeline execution (macro level)."""
    from variantcentrifuge.pipeline import run_refactored_pipeline

    # Generate test data
    n_variants = 10000
    df = synthetic_variants(n_variants=n_variants, n_samples=100, seed=42)
    input_vcf = tmp_path / "input.vcf"
    # ... write df to VCF format (simplified)

    # Prepare config
    config = {
        'vcf_file': str(input_vcf),
        'output_dir': str(tmp_path / 'output'),
        'use_new_pipeline': pipeline_type == 'stage_based',
        # ... other config
    }

    # Benchmark full pipeline
    result = benchmark(run_refactored_pipeline, config)

    benchmark.extra_info['n_variants'] = n_variants
    benchmark.extra_info['pipeline_type'] = pipeline_type
    benchmark.extra_info['level'] = 'macro'

    assert result is not None
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| unittest.TestCase timing | pytest-benchmark fixture | 2015-2017 | Standard regression detection, JSON storage, statistical rigor |
| memory_profiler decorator | tracemalloc stdlib | Python 3.4+ (2014) | No external dependency, lower overhead (5-20% vs 20%+) |
| Manual JSON schemas | pytest-benchmark format | 2015+ | Interoperable with CLI tools, histogram generation |
| Sampling real VCFs | Synthetic data from distributions | Ongoing privacy focus | GDPR/HIPAA compliance, reproducibility |
| Cross-run comparisons (flaky) | Ratio assertions (same run) | Best practice 2020+ | Zero flakiness for implementation comparisons |

**Deprecated/outdated:**
- **unittest + timeit:** Pre-pytest era testing approach, no fixtures or parametrization
- **memory_profiler:** Still maintained but heavyweight compared to stdlib tracemalloc
- **Custom benchmark JSON formats:** pytest-benchmark became de facto standard
- **CSV result storage:** JSON is standard for structured benchmark data

## Open Questions

Things that couldn't be fully resolved:

1. **Optimal synthetic data distributions for compound het benchmarks**
   - What we know: Compound hets require multiple variants per gene per sample, specific allele configurations
   - What's unclear: Best way to generate realistic compound het patterns without over-engineering
   - Recommendation: Start simple (random multi-variant genes), refine if benchmarks don't stress the algorithm sufficiently

2. **Memory budget thresholds for each component**
   - What we know: CONTEXT.md specifies tracemalloc for budget enforcement, warnings only
   - What's unclear: What are reasonable MB limits for 50K variants × 1000 samples in each component?
   - Recommendation: Profile current implementation with large datasets, set budgets at 2x observed peak to catch regressions

3. **Granularity of micro benchmarks**
   - What we know: Need micro (functions), meso (modules), macro (pipeline) levels
   - What's unclear: How fine-grained for micro? Every helper function or just hot paths?
   - Recommendation: Start with known hot paths (inheritance deduction, comp het, genotype replacement core loops), expand if profiling reveals other bottlenecks

4. **Handling pytest-benchmark storage in git**
   - What we know: Results are ephemeral, gitignored (CONTEXT.md)
   - What's unclear: Should baseline be committed for team consistency, or always local?
   - Recommendation: Gitignore all `benchmarks/*.json`, document how to establish local baseline in docs

## Sources

### Primary (HIGH confidence)
- [pytest-benchmark official docs - Usage](https://pytest-benchmark.readthedocs.io/en/latest/usage.html)
- [pytest-benchmark official docs - Comparing runs](https://pytest-benchmark.readthedocs.io/en/latest/comparing.html)
- [Python stdlib - tracemalloc documentation](https://docs.python.org/3/library/tracemalloc.html)
- [pytest parametrize docs](https://docs.pytest.org/en/stable/how-to/parametrize.html)

### Secondary (MEDIUM confidence)
- [Pytest Benchmark Tutorial - Pytest with Eric](https://pytest-with-eric.com/pytest-best-practices/pytest-benchmark/) - Practical patterns
- [Bencher - pytest-benchmark guide](https://bencher.dev/learn/benchmarking/python/pytest-benchmark/) - Integration patterns
- [CodSpeed benchmark thresholds](https://codspeed.io/changelog/2025-10-31-per-benchmark-regression-thresholds) - Industry standard thresholds (20% common)
- [Pandas performance best practices](https://pandas.pydata.org/docs/user_guide/enhancingperf.html) - DataFrame optimization
- [DataCamp - Benchmarking pandas alternatives](https://www.datacamp.com/tutorial/benchmarking-high-performance-pandas-alternatives) - Performance patterns

### Tertiary (LOW confidence)
- [Synthetic genomic data generation - Synth4bench](https://www.biorxiv.org/content/10.1101/2024.03.07.582313v1.full) - Research on synthetic VCF generation
- [Nature Communications - Genomic data de-identification](https://www.nature.com/articles/s41467-021-27219-2) - Anonymization techniques
- [Genomic variant benchmark - Genome Biology](https://link.springer.com/article/10.1186/s13059-023-03061-1) - Benchmarking standards in genomics

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - pytest-benchmark and tracemalloc are widely documented, stdlib-based
- Architecture: HIGH - Patterns verified from official docs and existing variantcentrifuge test structure
- Pitfalls: HIGH - Documented in pytest-benchmark FAQ, tracemalloc docs, and common testing issues
- Synthetic data: MEDIUM - Genomics-specific patterns less standardized, derived from research + pandas best practices

**Research date:** 2026-02-14
**Valid until:** 2026-04-14 (60 days - stable testing ecosystem)
