# Performance Benchmarking

This directory contains performance benchmarking tools for comparing the old and new pipeline implementations.

## Overview

The benchmarking framework measures:
- **Execution Time**: Total runtime for each configuration
- **Memory Usage**: Peak memory consumption during execution
- **CPU Utilization**: Average CPU usage percentage
- **Throughput**: Variants processed per second
- **Scalability**: Performance with different thread counts

## Quick Start

### Basic Benchmark
```bash
python tests/performance/benchmark_pipeline.py test_cohort.vcf.gz
```

### Standard Benchmark Suite
```bash
python tests/performance/benchmark_pipeline.py test_cohort.vcf.gz \
    --output-dir benchmark_results \
    --runs 3 \
    --config standard
```

### Quick Test (Single Gene Only)
```bash
python tests/performance/benchmark_pipeline.py test_cohort.vcf.gz \
    --config quick \
    --runs 1
```

## Benchmark Configurations

### Standard Tests

1. **Single Gene Tests**
   - `single_gene_1thread`: Baseline single-threaded performance
   - `single_gene_4thread`: Parallelization benefit for single gene

2. **Gene Set Tests**
   - `small_geneset_1thread`: 5 genes, single thread
   - `small_geneset_4thread`: 5 genes, parallel execution
   - `medium_geneset_1thread`: 50 genes, single thread
   - `medium_geneset_8thread`: 50 genes, high parallelization

3. **Feature Tests**
   - `complex_filter`: Multiple filter presets (rare + coding + pathogenic)
   - `with_scoring`: Variant scoring overhead

### Configuration Sets

- **quick**: Single gene tests only (2 tests)
- **standard**: All standard tests (8 tests)
- **extensive**: Standard + large gene sets (9+ tests)

## Output Files

The benchmark creates the following outputs:

```
benchmark_results/
├── benchmark_results.csv      # Raw timing data
├── performance_report.json    # Detailed analysis
├── plots/
│   ├── performance_comparison.png
│   ├── scalability_plot.png
│   └── speedup_comparison.png
├── results/                   # Pipeline output files
└── logs/                      # Pipeline log files
```

## Understanding Results

### Performance Report
The JSON report contains:
- System information (CPU, memory, platform)
- Summary statistics (mean, std, min, max times)
- Speedup calculations
- Memory usage comparisons

### Plots

1. **Performance Comparison**
   - Bar charts comparing execution time and memory usage
   - Error bars show standard deviation across runs

2. **Scalability Plot**
   - Shows performance scaling with thread count
   - Helps identify parallelization efficiency

3. **Speedup Comparison**
   - Visual representation of new vs old pipeline performance
   - Green bars: >20% improvement
   - Orange bars: Similar performance (±10%)
   - Red bars: Performance regression

## Example Results

```
BENCHMARK SUMMARY
==================
Overall speedup: 2.15x

Speedup by test:
  single_gene_1thread: 1.05x (4.8% improvement)
  single_gene_4thread: 2.85x (64.9% improvement)
  small_geneset_1thread: 1.12x (10.7% improvement)
  small_geneset_4thread: 3.42x (70.8% improvement)
  medium_geneset_8thread: 4.21x (76.2% improvement)
```

## Performance Tips

### For Optimal Benchmarking
1. Close other applications to reduce system noise
2. Run multiple iterations (--runs 3 or more)
3. Use a consistent test VCF file
4. Ensure system is not in power-saving mode

### Interpreting Results
- Small variations (±5%) between runs are normal
- Thread scaling depends on:
  - Number of genes being processed
  - Available CPU cores
  - I/O bottlenecks
- Memory usage should be relatively consistent

## Advanced Usage

### Custom Configurations
Edit `benchmark_pipeline.py` to add custom test configurations:

```python
BenchmarkConfig(
    name="custom_test",
    gene_names=["GENE1", "GENE2"],
    threads=16,
    preset="rare",
    extra_args=["--final-filter", "AF < 0.001"],
    description="Custom high-thread test",
)
```

### Automated Testing
Integrate into CI/CD:

```yaml
- name: Performance Benchmark
  run: |
    python tests/performance/benchmark_pipeline.py test.vcf.gz \
      --config quick \
      --runs 1 \
      --output-dir ci_benchmarks

    # Check for regressions
    python -c "
    import json
    with open('ci_benchmarks/performance_report.json') as f:
        report = json.load(f)
    speedup = report['summary']['overall_speedup']
    if speedup < 0.95:
        print(f'Performance regression detected: {speedup:.2f}x')
        exit(1)
    "
```

## Troubleshooting

### High Variability
If results show high standard deviation:
- Increase number of runs
- Check for background processes
- Verify consistent system load

### Memory Measurements
Memory tracking includes all child processes. Peaks may occur during:
- Initial VCF loading
- BED file creation
- Report generation

### Missing Dependencies
Ensure all required packages are installed:
```bash
pip install psutil matplotlib pandas numpy
```
