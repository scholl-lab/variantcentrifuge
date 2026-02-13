# Feature Landscape — Performance Optimization for Pandas-Heavy Bioinformatics Pipelines

**Domain:** Performance optimization of pandas-based genomic variant analysis pipelines
**Researched:** 2026-02-14
**Target pipeline:** variantcentrifuge (Python, pandas, 22 GB VCFs, 5,125 samples, 2,500 genes)

---

## Table Stakes

Features users expect from meaningful performance optimization. Missing = optimization feels incomplete or superficial.

| Feature | Why Expected | Complexity | Expected Impact | Notes |
|---------|--------------|------------|-----------------|-------|
| **Eliminate df.apply(axis=1)** | Standard pandas optimization practice; apply(axis=1) is 10-100x slower than vectorization | Medium | 40-60% of analysis time saved | Current: `row.to_dict()` creates 5,000+ key dicts per row. Solution: itertuples (13x faster) or vectorization (100-740x faster) |
| **Replace iterrows() loops** | iterrows is 13-50x slower than itertuples; universally considered anti-pattern | Low | 10-13x speedup on affected sections | Three instances in inheritance analysis. itertuples returns namedtuples vs full Series objects |
| **Add observed=True to groupby** | Prevents catastrophic 3,500x slowdown with categorical dtypes | Trivial | Prevents future regression | 12 instances across codebase. Currently no slowdown (GENE is string), but blocks dtype optimization |
| **Use typed dtypes (not dtype=str)** | dtype=str wastes 3-5x memory, prevents NumPy vectorization | Medium | 50-80% memory reduction | Convert CHROM→category (90% reduction), POS→int64, AF→float64. Enables numeric operations |
| **Enable PyArrow engine for CSV I/O** | PyArrow is 2-3x faster than default C engine; industry standard since 2024 | Low | 2-3x I/O speedup | Add `engine="pyarrow"` to 16 pd.read_csv calls. Multi-threaded decompression |
| **Remove dead code** | Dead code consumes CPU for no value; profiling should identify it | Low | 20-30% of gene burden stage | GT parsing loop (lines 220-249 in gene_burden.py) does nothing useful; results never used |

---

## Differentiators

Features that set optimization apart. Not expected everywhere, but deliver maximum impact or competitive advantage.

| Feature | Value Proposition | Complexity | Expected Impact | Notes |
|---------|-------------------|------------|-----------------|-------|
| **Vectorize inheritance deduction** | Inheritance analysis is 40-60% of total time; vectorization eliminates row-wise dict creation | High | 100-740x on hottest path | Replace lambda row: deduce_patterns_for_variant(row.to_dict(), ...) with NumPy matrix operations on genotype array |
| **Switch openpyxl→xlsxwriter for Excel** | openpyxl writes cells one-by-one; xlsxwriter streams data 2-5x faster | Low | 2-5x Excel generation speedup | Case study: 200K rows 9min→3min. Use xlsxwriter for initial write, openpyxl for finalization only |
| **Pre-parse GT column once** | GT string parsing happens per-gene; pre-parsing eliminates redundant work | Medium | Eliminates re-parsing overhead | Parse "Sample1(0/1:AD=20,30);Sample2(1/1)..." into structured dict ONCE at load time |
| **Categorical dtypes for low-cardinality** | CHROM (~25 values), IMPACT (4), FILTER (3) waste 90-99% memory as strings | Low | 50-90% memory reduction | Category dtype stores as int array + string mapping. 8x reduction demonstrated for country column |
| **Use Parquet/Feather for intermediate files** | CSV is 4-15x slower than columnar formats | Medium | 4-15x intermediate I/O speedup | Feather: 2.3s read vs 15s CSV (5M records). Parquet: 2x compression vs CSV. Replace write-then-read temp TSV |
| **Pipe fusion for external tools** | Sequential subprocess calls with temp files waste I/O | Medium | Eliminates 2 temp VCF files | Current: bcftools → temp → SnpSift → temp → extractFields. Proposed: bcftools \| SnpSift \| extractFields |
| **Pass DataFrames in memory** | Excel stage re-reads 2-10 GB TSV from disk that's already in memory | Trivial | Eliminates redundant disk read | Modify ExcelReportStage to accept df from PipelineContext instead of path |
| **Chunked DataFrame processing** | Enables processing datasets larger than memory | Medium | Unlimited dataset size | Use pd.read_csv(chunksize=...) for 10-100 MB chunks. Only 10K rows in memory at once |
| **Numba JIT compilation for hot loops** | Numba @jit compiles Python→native code for 30x speedup on numeric ops | Medium | 10-50x on JIT-compatible code | Best for NumPy array operations. First call has compilation overhead; subsequent calls cached. Requires ≥1M data points for benefit |
| **Use eval()/query() for complex expressions** | Numexpr evaluates expressions 50% faster without temp arrays | Low | 50% faster, 40% less memory | Only beneficial for DataFrames >100K rows. Use df.eval("new_col = col1 + col2") instead of df["new_col"] = df["col1"] + df["col2"] |

---

## Anti-Features

Features to explicitly NOT build. Common mistakes or premature optimizations in this domain.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Migrate to Polars** | 10-100x speedup potential but requires complete rewrite of 40+ stages; high risk | Polars is faster (5x CSV load, 11x sort, 4.6x filter), but pandas with optimizations achieves 3-4x target speedup without migration risk |
| **GPU acceleration (CUDA/cuDF)** | Genomic pipelines are I/O-bound and logic-heavy, not numeric-heavy; GPU overhead > benefit | VCF processing involves string parsing, conditional logic, pedigree lookups. GPU only helps for pure matrix operations (not present here) |
| **SharedMemory for parallel stages** | Adds complexity for minimal benefit when DataFrames are small (<100 MB per chunk) | Gene-level DataFrames are 1-20 variants each; pickle overhead is negligible. Only consider for DataFrames >1 GB passed between processes |
| **Cython/Rust for all Python code** | High development cost; only justifiable for proven bottlenecks after other optimizations | Profile first. Only rewrite in compiled language if: (1) vectorization impossible, (2) proven 80%+ time consumer, (3) performance still inadequate |
| **Custom C extensions for VCF parsing** | Reinventing bcftools/cyvcf2; bioinformatics community already solved this | Use cyvcf2 (wraps htslib in Python) for fast VCF parsing. Don't write custom parsers |
| **Async I/O for file operations** | Genomic files are large (GB scale); async doesn't help when operations are CPU/bandwidth-bound | ThreadPoolExecutor for parallel processing is already implemented. Async benefits network I/O, not local disk with large files |
| **In-memory databases (SQLite/DuckDB)** | Adds query layer overhead; pandas is already in-memory and optimized for DataFrame operations | Consider DuckDB only if queries involve complex joins across multiple large tables (not the case here) |
| **Modin (parallelized pandas)** | Drop-in replacement with mixed results; doesn't address root cause (inefficient operations) | Modin parallelizes pandas but doesn't fix apply(axis=1) or iterrows(). Fix the operations first; parallelism second |

---

## Feature Dependencies

```
Typed dtypes (table stakes)
    ↓
Categorical dtypes (differentiator) — Requires observed=True on groupby
    ↓
Memory reduction unlocks larger gene panels
    ↓
Vectorization (differentiator) — Requires typed dtypes for NumPy ops
    ↓
Numba JIT (differentiator) — Requires vectorized NumPy arrays

PyArrow engine (table stakes)
    ↓
Parquet/Feather intermediate files (differentiator)

Remove dead code (table stakes) — Prerequisite for accurate profiling
    ↓
Pre-parse GT column (differentiator) — Eliminates what dead code was attempting

Pipe fusion (differentiator) — Independent; modifies external tool orchestration
```

---

## MVP Recommendation

For **3-4x full pipeline speedup** with minimal risk, prioritize **in this order**:

### Phase 1: Quick Wins (1-2 hours, 30-40% speedup)
1. **Remove dead code** (gene_burden.py:220-249) — 20-30% of gene burden time
2. **Add observed=True** to all 12 groupby calls — Prevents future 3,500x slowdown
3. **Enable PyArrow engine** for CSV reads — 2-3x I/O speedup
4. **Fix temp file leak** (gene_bed.py:153) — Correctness, not performance

### Phase 2: Core Optimizations (2-3 days, 2-3x additional speedup)
5. **Replace iterrows() with itertuples()** — 10-13x on passes 2-3 of inheritance analysis
6. **Replace apply(axis=1) with vectorized operations** — 100-740x on pass 1 (40-60% of total time)
7. **Categorical dtypes** for CHROM, IMPACT, FILTER — 50-90% memory reduction
8. **Pass DataFrame directly to Excel stage** — Eliminates 2-10 GB disk read

### Phase 3: Strategic Enhancements (1-2 weeks, pushes to 4x target)
9. **Pre-parse GT column once** — Eliminates redundant parsing across stages
10. **xlsxwriter for Excel generation** — 2-5x Excel write speed
11. **Pipe fusion** for bcftools→SnpSift→extractFields — Eliminates 2 intermediate VCF files
12. **Parquet intermediate files** — Replace TSV write-then-read with 4-15x faster format

---

## Defer to Post-MVP

| Feature | Reason to Defer | Consider When |
|---------|----------------|---------------|
| **Numba JIT compilation** | Requires extensive vectorization first; only beneficial for >1M rows | After Phase 2 vectorization complete; if profiling shows numeric bottlenecks remain |
| **Chunked processing** | Only needed if DataFrames exceed available memory | Dataset exceeds 50% of system RAM; currently 22 GB VCF → 2-10 GB DataFrame fits in 32 GB RAM |
| **eval()/query() expressions** | Only 50% speedup; smaller impact than other optimizations | After Phase 3 complete; if complex filter expressions identified in profiling |
| **Polars migration** | Major rewrite; pandas optimizations achieve target speedup | Only if optimized pandas still 3x too slow; require compelling 10x+ improvement case |
| **Cython/Rust extensions** | High development cost; only for proven bottlenecks after other fixes | Genotype replacement kernel if still >30% of pipeline time after Phase 3 |

---

## Complexity Ratings Explained

| Rating | Definition | Examples |
|--------|------------|----------|
| **Trivial** | Single-line change, search-and-replace | Add observed=True, pass df in memory |
| **Low** | 1-4 hours, localized changes, minimal testing | iterrows→itertuples, add PyArrow engine, xlsxwriter |
| **Medium** | 1-3 days, requires logic refactoring, moderate testing | Vectorize simple operations, typed dtypes, pre-parse GT column, Parquet files |
| **High** | 1-2 weeks, complex logic transformation, extensive testing | Vectorize inheritance deduction (100+ lines, nested conditionals), Numba integration |

---

## Expected Speedup: Full Pipeline Estimate

**Baseline:** 90 minutes for GCKD (22 GB VCF, 5,125 samples, 2,500 genes)

| Phase | Optimizations | Estimated Time | Cumulative Speedup |
|-------|---------------|----------------|-------------------|
| Baseline | Current implementation | 90 min | 1.0x |
| + Phase 1 | Dead code removal, observed=True, PyArrow | 60 min | 1.5x |
| + Phase 2 | itertuples, vectorization, categoricals, in-memory df | 25 min | 3.6x |
| + Phase 3 | Pre-parse GT, xlsxwriter, pipe fusion, Parquet | 18 min | 5.0x |

**Target achieved:** 3-4x speedup with Phase 1 + Phase 2 (2-4 days effort)

**Stretch goal:** 5x speedup with Phase 3 (2-3 weeks total)

---

## Feature Prioritization Matrix

| Feature | Impact | Effort | Priority | Rationale |
|---------|--------|--------|----------|-----------|
| Remove dead code | High | Trivial | **P0** | Zero-risk, immediate 20-30% gain |
| observed=True | Critical (future) | Trivial | **P0** | Blocks categorical dtype optimization |
| PyArrow engine | Medium | Low | **P0** | 2-3x I/O speedup, low risk |
| iterrows→itertuples | High | Low | **P1** | 10-13x speedup, proven technique |
| Categorical dtypes | High | Low-Med | **P1** | 50-90% memory reduction |
| Vectorize apply(axis=1) | Critical | High | **P1** | 100-740x on hottest path; highest impact |
| Pass df in memory | Medium | Trivial | **P2** | Easy win, eliminates redundant I/O |
| Pre-parse GT column | Medium-High | Medium | **P2** | Eliminates cross-stage redundancy |
| xlsxwriter | Medium | Low | **P2** | 2-5x Excel speedup; 10% of pipeline |
| Pipe fusion | Low-Medium | Medium | **P3** | Saves 2 temp files; one-time per run |
| Parquet files | Medium | Medium | **P3** | 4-15x intermediate I/O; smaller impact |
| Numba JIT | Medium | High | **P4** | Requires vectorization first |
| Polars migration | High | Very High | **P5** | High risk; defer until pandas optimizations exhausted |

---

## Performance Validation Strategy

### Regression Testing
```python
# tests/performance/test_benchmarks.py
import pytest
from pytest_benchmark.fixture import BenchmarkFixture

@pytest.fixture
def sample_dataframe():
    """5,125 samples × 50K variants (simulates GCKD scale)."""
    return pd.DataFrame(...)

def test_inheritance_deduction_performance(benchmark: BenchmarkFixture, sample_dataframe):
    """Baseline: 30 min → Target: 3 min (10x speedup)."""
    result = benchmark.pedantic(
        analyze_inheritance,
        args=(sample_dataframe, pedigree_data, sample_list),
        rounds=3,
        warmup_rounds=1,
    )
    assert benchmark.stats.median < 180  # 3 minutes

def test_gene_burden_performance(benchmark: BenchmarkFixture, burden_df):
    """Baseline: 15 min → Target: 5 min (3x speedup)."""
    result = benchmark(compute_gene_burden, burden_df, stats_config)
    assert benchmark.stats.median < 300  # 5 minutes
```

### Profiling Checkpoints
```bash
# Profile before each optimization phase
scalene --html --outfile phase1_before.html variantcentrifuge/cli.py -- <args>
# ... apply Phase 1 optimizations ...
scalene --html --outfile phase1_after.html variantcentrifuge/cli.py -- <args>
# Compare: diff_scalene.py phase1_before.html phase1_after.html
```

### Memory Profiling
```bash
# Track memory reduction from dtype optimizations
python -c "
import tracemalloc
tracemalloc.start()
df = pd.read_csv('large_file.tsv', sep='\t', dtype=str)  # Baseline
snapshot_before = tracemalloc.take_snapshot()

df = pd.read_csv('large_file.tsv', sep='\t', dtype=typed_dtypes)  # Optimized
snapshot_after = tracemalloc.take_snapshot()

for stat in snapshot_after.compare_to(snapshot_before, 'lineno')[:10]:
    print(stat)
"
```

---

## Sources

### Pandas Performance Optimization
- [Enhancing performance — pandas 3.0.0 documentation](https://pandas.pydata.org/docs/user_guide/enhancingperf.html)
- [Pandas Performance Best Practices | Compile N Run](https://www.compilenrun.com/docs/library/pandas/pandas-performance/pandas-performance-best-practices/)
- [Performance Optimization in Pandas: Techniques & Best Practices | Medium](https://medium.com/python-for-engineers/performance-optimization-in-pandas-techniques-best-practices-part-1-690d8eddae21)
- [Advanced Pandas Techniques for Data Processing and Performance | Towards Data Science](https://towardsdatascience.com/advanced-pandas-techniques-for-data-processing-and-performance-23f9026d21f5/)

### Vectorization vs apply()
- [Pandas Apply vs Vectorization: The Benchmark Story | Medium](https://medium.com/@hadiyolworld007/pandas-apply-vs-vectorization-the-benchmark-story-bcd981ffabb0)
- [Pandas Apply vs. Vectorization: The Performance Showdown | Medium](https://medium.com/@ThinkingLoop/pandas-apply-vs-vectorization-the-performance-showdown-730043f1d920)
- [Pandas vectorization: faster code, slower code, bloated memory](https://pythonspeed.com/articles/pandas-vectorization/)
- [Tutorial: 700x speed improvement on Pandas with Vectorization vs. Iterrows/Apply](https://www.linkedin.com/pulse/tutorial-basic-vectorization-pandas-iterrows-apply-duc-lai-trung-minh-75d4c)

### iterrows vs itertuples
- [Why Pandas itertuples() Is Faster Than iterrows() | Medium](https://medium.com/swlh/why-pandas-itertuples-is-faster-than-iterrows-and-how-to-make-it-even-faster-bc50c0edd30d)
- [For the Love of God, Stop Using iterrows() – r y x, r](https://ryxcommar.com/2020/01/15/for-the-love-of-god-stop-using-iterrows/)
- [Looping over Pandas data | Architecture & Performance](https://aetperf.github.io/2018/07/03/Looping-over-Pandas-data.html)

### Dtype Optimization & Categorical
- [Mastering Memory Optimization for Pandas DataFrames](https://thinhdanggroup.github.io/pandas-memory-optimization/)
- [Reducing Pandas memory usage #1: lossless compression](https://pythonspeed.com/articles/pandas-load-less-data/)
- [Using The Pandas Category Data Type - Practical Business Python](https://pbpython.com/pandas_dtypes_cat.html)
- [Seven Killer Memory Optimization Techniques Every Pandas User Should Know | Towards Data Science](https://towardsdatascience.com/seven-killer-memory-optimization-techniques-every-pandas-user-should-know-64707348ab20/)

### groupby observed=True
- [pandas.DataFrame.groupby — pandas 3.0.0 documentation](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.groupby.html)
- [What's new in 3.0.0 (January 21, 2026) — pandas documentation](https://pandas.pydata.org/docs/whatsnew/v3.0.0.html)
- [Be Careful When Using Pandas Groupby with Categorical Data Type | Towards Data Science](https://towardsdatascience.com/be-careful-when-using-pandas-groupby-with-categorical-data-type-a1d31f66b162/)
- [Pandas GroupBy at Speed: Pitfalls and Power Moves | Medium](https://medium.com/@2nick2patel2/pandas-groupby-at-speed-pitfalls-and-power-moves-8b27ca7ccc5a)

### File Formats (Parquet, Feather, CSV)
- [Saving Pandas DataFrames Efficiently — Parquet vs Feather vs ORC vs CSV | Towards Data Science](https://towardsdatascience.com/saving-pandas-dataframes-efficiently-and-quickly-parquet-vs-feather-vs-orc-vs-csv-26051cc98f2e/)
- [Speed Comparison - CSV, Feather, Pickle, HDF5, and Parquet](https://kimmonzon.com/csv-feather-pickle-hdf5-and-parquet/)
- [Feather vs Parquet vs CSV vs Jay | Medium](https://bawaji94.medium.com/feather-vs-parquet-vs-csv-vs-jay-55206a9a09b0)

### Excel Performance
- [Performance — openpyxl 3.1.4 documentation](https://openpyxl.readthedocs.io/en/3.1/performance.html)
- [Openpyxl vs XlsxWriter: The Ultimate Showdown for Excel Automation](https://hive.blog/python/@geekgirl/openpyxl-vs-xlsxwriter-the-ultimate-showdown-for-excel-automation)
- [Optimizing Excel Report Generation: From OpenPyXL to XLSMerge | Medium](https://mass-software-solutions.medium.com/optimizing-excel-report-generation-from-openpyxl-to-xlsmerge-processing-52-columns-200k-rows-5b5a03ecbcd4)
- [Benchmark of several Python Excel writing modules · GitHub](https://gist.github.com/jmcnamara/ba25c2bf4ba0777065eb)

### Numba & JIT Compilation
- [Enhancing performance — pandas 3.0.0 documentation](https://pandas.pydata.org/docs/user_guide/enhancingperf.html)
- [A ~5 minute guide to Numba](https://numba.readthedocs.io/en/stable/user/5minguide.html)
- [This Decorator will Make Python 30 Times Faster | Towards Data Science](https://towardsdatascience.com/this-decorator-will-make-python-30-times-faster-715ca5a66d5f/)

### Polars Comparison
- [Polars vs pandas: What's the Difference? – Real Python](https://realpython.com/polars-vs-pandas/)
- [Pandas vs Polars: A Comprehensive Performance Comparison | Medium](https://medium.com/@neurogenou/pandas-vs-polars-a-comprehensive-performance-comparison-31a9296e4cd4)
- [Polars vs. Pandas — An Independent Speed Comparison | Towards Data Science](https://towardsdatascience.com/polars-vs-pandas-an-independent-speed-comparison/)

### Chunked Processing
- [Scaling to large datasets — pandas 3.0.0 documentation](https://pandas.pydata.org/docs/user_guide/scale.html)
- [Scalable Data Processing with Pandas: Handling Large CSV Files in Chunks | Medium](https://medium.com/itversity/scalable-data-processing-with-pandas-handling-large-csv-files-in-chunks-b15c5a79a3e3)
- [Reducing Pandas memory usage #3: Reading in chunks](https://pythonspeed.com/articles/chunking-pandas/)

### eval() and query()
- [High-Performance Pandas: eval() and query() | Python Data Science Handbook](https://jakevdp.github.io/PythonDataScienceHandbook/03.12-performance-eval-and-query.html)
- [Get More Efficient in your Data Analysis with Pandas query() and eval() Methods | Towards Data Science](https://towardsdatascience.com/get-more-efficient-in-your-data-analysis-with-pandas-query-and-eval-methods-3646317e591f/)

### Bioinformatics-Specific
- [How to Spot (and Fix) 5 Common Performance Bottlenecks in pandas Workflows | NVIDIA Technical Blog](https://developer.nvidia.com/blog/how-to-spot-and-fix-5-common-performance-bottlenecks-in-pandas-workflows/)
- [Analysis-ready VCF at Biobank scale using Zarr | GigaScience](https://doi.org/10.1093/gigascience/giaf049)
- [cyvcf2: fast, flexible variant analysis with Python | Bioinformatics](https://academic.oup.com/bioinformatics/article/33/12/1867/2971439)

### Dead Code Detection & Profiling
- [vulture · PyPI](https://pypi.org/project/vulture/)
- [GitHub - albertas/deadcode: Find and fix unused Python code](https://github.com/albertas/deadcode)
- [Top 7 Python Profiling Tools for Performance](https://daily.dev/blog/top-7-python-profiling-tools-for-performance)
