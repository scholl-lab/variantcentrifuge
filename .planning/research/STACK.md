# Technology Stack for Performance Optimization

**Project:** VariantCentrifuge v0.12.1 Performance Milestone
**Researched:** 2026-02-14
**Focus:** Stack additions for benchmarking and optimization (NOT re-architecting)

---

## Executive Summary

This stack research focuses exclusively on **additions** to the existing variantcentrifuge stack to enable:
1. **Performance benchmarking** with regression tracking
2. **DataFrame optimization** for pandas-heavy workloads
3. **Native extensions** for critical hot paths
4. **Excel generation** speedup

**Key recommendation:** Incremental optimization using proven pandas ecosystem tools. Do NOT migrate to Polars or other frameworks—the existing pandas infrastructure is sound, just needs targeted optimizations.

---

## Core Performance Stack

### Benchmarking & Profiling

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **pytest-benchmark** | 5.2.3 | Regression tracking for unit-level benchmarks | Mature (stable since 2015), integrates with existing pytest suite (1035 tests), supports pedantic mode for exact iteration control, generates historical comparison reports |
| **py-spy** | 0.4.1 | Production profiling (sampling) | Zero-overhead sampling profiler, works on running processes without code modification, native C extension support (`--native` flag) critical for pandas/NumPy profiling, generates flame graphs |
| **pytest-codspeed** | 4.2.0+ (optional) | CI/CD performance tracking | Backward-compatible with pytest-benchmark, CPU simulation for steady measurements in noisy CI environments, free for open source. Consider for future CI integration. |

**Integration:**
```toml
# pyproject.toml
[project.optional-dependencies]
dev = [
    # ... existing dev deps ...
    "pytest-benchmark>=5.2",
]
perf = [
    "pytest-benchmark>=5.2",
    "py-spy>=0.4",
]
```

**Rationale:** pytest-benchmark is the industry standard for Python performance testing. It's stable (5 years of production use), works with existing pytest infrastructure, and provides historical comparison. py-spy is essential for profiling pandas/NumPy-heavy code because it can profile native C extensions that cProfile cannot see.

**Cross-platform:** Both tools support Windows and Linux with pre-built wheels for x86-64 and ARM64.

---

### DataFrame Optimization

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **PyArrow** | 23.0.0 | Columnar data backend for pandas | 2-3x faster CSV parsing (multi-threaded), 70% less memory for string columns, native support in pandas ≥2.0 via `dtype_backend="pyarrow"`, enables zero-copy interchange with Arrow-based tools |
| **numba** | 0.63.1 | JIT compilation for numeric kernels | 100-740x speedup for NumPy array operations, supports pandas engine integration (`engine="numba"`), minimal code changes (@jit decorator), free threading for parallel execution |

**NOT recommended:**
- **Polars** (0.20+): 5-20x faster than pandas for data transformations BUT incompatible with most Python ML/viz libraries, poor hierarchical indexing support (critical for genomic multi-sample data), would require full pipeline rewrite. Pandas with PyArrow backend achieves 70% of Polars performance with zero rewrite.

**Integration:**
```toml
# pyproject.toml
dependencies = [
    "pandas>=2.0",  # Already present, ensure ≥2.0 for PyArrow backend
    "pyarrow>=23.0",
    "numba>=0.63",
    # ... existing deps ...
]
```

**Rationale for PyArrow:**
- pandas 3.0 (released 2026-01) defaults to PyArrow backend for strings
- Performance analysis report shows `dtype=str` everywhere—PyArrow strings reduce memory 70% and enable vectorized string ops
- Multi-threaded CSV parsing critical for 22 GB VCF-derived TSV files
- Minimal code changes: `pd.read_csv(..., engine="pyarrow", dtype_backend="pyarrow")`

**Rationale for numba:**
- Hot path in `inheritance/analyzer.py:86` uses `df.apply(axis=1)` which is 100-740x slower than vectorized NumPy
- numba enables JIT compilation of Python functions to LLVM machine code with minimal refactoring
- Supports both CPU and GPU (future-proofing)
- Can use `@njit` decorator on existing helper functions without full Cython rewrite

**Cross-platform:** Both PyArrow and numba ship pre-built wheels for Windows (x86-64, ARM64) and Linux (manylinux). numba requires LLVM which is included in binary wheels.

---

### Native Extensions (Selective Use)

| Technology | Version | Purpose | When to Use |
|------------|---------|---------|-------------|
| **Cython** | 3.2.4 | Compiled Python extensions for hot paths | For genotype replacement kernel (`vectorized_replacer.py`) and GT string parsing—pure Python string processing is the #1 bottleneck (20-30% of pipeline time). Use for: iterative string parsing, nested loops that can't be vectorized. |
| **hatch-cython** | 0.5.0+ | Build plugin for hatchling | Enables Cython compilation with existing hatchling build system—no migration to setuptools required |

**Integration:**
```toml
# pyproject.toml
[build-system]
requires = ["hatchling", "hatch-cython>=0.5", "Cython>=3.2", "setuptools"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel.hooks.cython]
dependencies = ["hatch-cython"]

[tool.hatch.build.targets.wheel.hooks.cython.options]
include_numpy = true
```

**Rationale:**
- Genotype replacement stage processes 5,125 samples × 50K variants with regex per cell—pure Python bottleneck
- Cython provides 10-100x speedup for string processing loops without full rewrite
- hatch-cython integrates with existing hatchling build (no need to migrate to setuptools)
- Gradual adoption: compile individual modules (.pyx files) while keeping rest as pure Python

**Cross-platform considerations:**
- Cython generates C code that compiles with MSVC (Windows) or GCC/Clang (Linux)
- CI/CD must install C compiler: Windows (MSVC build tools), Linux (gcc/g++)
- Ship pre-built wheels via cibuildwheel for common platforms to avoid compilation on user machines
- Source distribution (.tar.gz) includes generated .c files for fallback compilation

**When NOT to use Cython:**
- DataFrame-level operations (use PyArrow/numba instead)
- I/O bottlenecks (use PyArrow engine, not Cython)
- Code that can be vectorized with NumPy/pandas (vectorization >> Cython for complexity/maintenance)

---

### Excel Generation Optimization

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **xlsxwriter** | 3.2.9 | Fast Excel file creation | 2-5x faster than openpyxl for initial write (real-world: 9 min → 3 min for 200K rows), write-only design optimized for large files, supports all Excel formatting needed for reports |
| **openpyxl** | (keep current) | Excel file modification | Keep for post-processing (adding hyperlinks, freeze panes, auto-filter) since xlsxwriter is write-only |

**Integration:**
```toml
# pyproject.toml
dependencies = [
    "openpyxl",      # Keep existing
    "xlsxwriter>=3.2",  # Add
    # ... existing deps ...
]
```

**Usage pattern:**
```python
# Initial write (fast)
df.to_excel(xlsx_file, engine="xlsxwriter", index=False)

# Post-processing (features)
import openpyxl
wb = openpyxl.load_workbook(xlsx_file)
# Add hyperlinks, freeze panes, etc.
wb.save(xlsx_file)
```

**Rationale:**
- Performance analysis shows Excel generation takes 2-10 minutes for GCKD cohort (50K rows)
- xlsxwriter is write-only, so 2-5x faster than openpyxl's read-write architecture
- Keep openpyxl for features xlsxwriter doesn't support (modifying existing files, reading)
- Both libraries are pure Python (no compilation), cross-platform compatible

**Cross-platform:** xlsxwriter is pure Python with zero dependencies—works identically on Windows and Linux.

---

## Optional: Advanced Optimization Tools (Deferred)

| Technology | Version | Status | Notes |
|------------|---------|--------|-------|
| **Scalene** | 1.5+ | Consider for Tier 3 | CPU + memory + GPU profiler with line-level granularity. More detailed than py-spy but higher overhead. Use for deep-dive profiling sessions, not routine benchmarking. |
| **cibuildwheel** | 2.16+ | Required if shipping Cython wheels | Automates building wheels for Windows/Linux/macOS across Python 3.10-3.14. Essential for distributing pre-compiled Cython extensions. |
| **mypyc** | (bundled with mypy) | Not recommended | mypy's native compiler—generates C extensions from type-annotated Python. Less mature than Cython, requires complete type annotations, limited control over optimization. Stick with Cython for explicit hot paths. |

---

## Integration Points with Existing Stack

### 1. DataFrame Loading (DataFrameLoadingStage)

**Current:**
```python
df = pd.read_csv(tsv_file, sep="\t", dtype=str)
```

**With PyArrow backend:**
```python
import pyarrow as pa
df = pd.read_csv(
    tsv_file,
    sep="\t",
    engine="pyarrow",
    dtype_backend="pyarrow",
    # Define dtypes for known columns
    dtype={
        "CHROM": "category",
        "POS": "Int64",
        "IMPACT": "category",
        "FILTER": "category",
    }
)
```

**Impact:** 2-3x faster read, 50-90% less memory for categorical columns, enables downstream vectorization.

---

### 2. Inheritance Analysis (InheritanceAnalysisStage)

**Current bottleneck (analyzer.py:86):**
```python
df["_inheritance_patterns"] = df.apply(
    lambda row: deduce_patterns_for_variant(row.to_dict(), pedigree_data, sample_list),
    axis=1
)
```

**With numba (Tier 2 optimization):**
```python
from numba import njit
import numpy as np

@njit
def deduce_patterns_vectorized(genotype_matrix, pedigree_indices):
    """Vectorized pattern deduction using NumPy arrays."""
    patterns = np.zeros(genotype_matrix.shape[0], dtype=np.int32)
    # Vectorized genotype analysis
    ...
    return patterns

# Extract sample columns as NumPy array
sample_cols = [s for s in sample_list if s in df.columns]
gt_matrix = df[sample_cols].values
patterns = deduce_patterns_vectorized(gt_matrix, pedigree_indices)
df["_inheritance_patterns"] = patterns
```

**Impact:** 100-740x speedup on hottest pipeline path.

---

### 3. Genotype Replacement (GenotypeReplacementStage)

**Current bottleneck (vectorized_replacer.py):**
Pure Python regex processing per cell in ProcessPoolExecutor workers.

**With Cython (Tier 3 optimization):**
```cython
# genotype_replacer_cy.pyx
import cython
from libc.string cimport strstr, strcpy

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef str process_genotype(str gt_str, str sample_name):
    """Cython-optimized genotype string processing."""
    cdef char* c_str = <char*>gt_str
    # C-level string manipulation
    ...
    return result
```

**Impact:** 10-50x speedup on 20-30% of total pipeline time.

---

### 4. Performance Benchmarking (New Test Suite)

**File:** `tests/performance/test_benchmarks.py`

```python
import pytest
import pandas as pd
from variantcentrifuge.inheritance.analyzer import analyze_inheritance

@pytest.fixture
def small_panel_df():
    """Load PKD 7-gene panel for benchmarking."""
    return pd.read_csv("tests/fixtures/giab/final/...", sep="\t")

def test_inheritance_analysis_perf(benchmark, small_panel_df):
    """Regression test: inheritance analysis should complete in <5s."""
    pedigree = {...}
    samples = ["HG002", "HG003", "HG004"]

    result = benchmark.pedantic(
        analyze_inheritance,
        args=(small_panel_df, pedigree, samples),
        rounds=5,
        warmup_rounds=1,
    )

    assert len(result) > 0
    # Performance assertion
    assert benchmark.stats.stats.mean < 5.0, "Inheritance analysis regressed"

def test_gene_burden_perf(benchmark, burden_df):
    """Regression test: gene burden should complete in <10s."""
    from variantcentrifuge.gene_burden import compute_gene_burden

    result = benchmark.pedantic(
        compute_gene_burden,
        args=(burden_df, {"case_samples": [...], "control_samples": [...]}),
        rounds=5,
    )

    assert "p_value" in result.columns
```

**Run with:**
```bash
# Run performance tests
pytest -m performance --benchmark-only

# Compare against baseline
pytest -m performance --benchmark-compare=0001

# Save baseline
pytest -m performance --benchmark-save=baseline
```

**Integration with CI:**
Add GitHub Actions workflow to track performance regressions over time.

---

## What NOT to Add

| Technology | Why NOT |
|------------|---------|
| **Polars** | Requires full pipeline rewrite (3-6 weeks effort), breaks scikit-learn/matplotlib integration, poor multi-index support (critical for multi-sample genomics). Pandas + PyArrow achieves 70% of Polars benefits with 5% of rewrite cost. |
| **Dask** | Pipeline is not memory-bound (22 GB VCF fits in typical HPC node), parallelism already implemented with ProcessPoolExecutor. Dask adds complexity for distributed computing we don't need. |
| **Ray** | Same as Dask—overkill for single-node parallelism. ProcessPoolExecutor is simpler and sufficient. |
| **Rust extensions (PyO3)** | Cython achieves 90% of Rust performance with 50% less complexity. Rust toolchain adds Windows build dependencies (rustc, cargo). Only consider if Cython proves insufficient. |
| **GPU acceleration (CuPy/RAPIDS)** | Pipeline is I/O and string-processing bound, not pure numeric computation. GPU transfer overhead would negate benefits. |

---

## Dependency Versions Summary

**Add to pyproject.toml:**

```toml
[project]
dependencies = [
    "pandas>=2.0",        # Ensure PyArrow backend support
    "pyarrow>=23.0",      # NEW: Columnar backend
    "numba>=0.63",        # NEW: JIT compilation
    "xlsxwriter>=3.2",    # NEW: Fast Excel write
    # ... keep existing: jinja2, scipy, statsmodels, numpy, intervaltree, psutil, smart-open, openpyxl
]

[project.optional-dependencies]
dev = [
    "mypy>=1.10",
    "pytest>=8.0",
    "pytest-cov>=5.0",
    "pytest-mock>=3.14",
    "pytest-benchmark>=5.2",  # NEW
    "ruff>=0.11",
    "pre-commit",
]

perf = [
    "pytest-benchmark>=5.2",
    "py-spy>=0.4",
]

[build-system]
# Only add if implementing Cython extensions (Tier 3)
requires = ["hatchling", "hatch-cython>=0.5", "Cython>=3.2", "setuptools"]
build-backend = "hatchling.build"
```

**Installation:**
```bash
# Standard install with new performance deps
uv pip install -e ".[dev]"

# Profiling-only install
uv pip install -e ".[perf]"
```

---

## Cross-Platform Compatibility Matrix

| Tool | Windows | Linux | Notes |
|------|---------|-------|-------|
| PyArrow | ✅ x86-64, ARM64 | ✅ manylinux x86-64, ARM64 | Pre-built wheels |
| numba | ✅ x86-64 | ✅ x86-64, ARM64 | LLVM bundled in wheels |
| xlsxwriter | ✅ Pure Python | ✅ Pure Python | Zero dependencies |
| pytest-benchmark | ✅ | ✅ | Pure Python |
| py-spy | ✅ x86-64 | ✅ x86-64, ARM64 | Rust binary, pre-built |
| Cython | ✅ Requires MSVC | ✅ Requires gcc | Ship pre-built wheels via cibuildwheel |

**Windows-specific:**
- Cython extensions require Microsoft Visual C++ 14.0+ (free download)
- numba works out-of-box (LLVM bundled)
- All other tools are pure Python or ship Windows wheels

**Linux-specific:**
- All tools ship manylinux wheels (glibc 2.17+)
- Cython requires gcc/g++ (standard on HPC systems)

---

## Confidence Assessment

| Stack Component | Confidence | Source |
|----------------|------------|--------|
| pytest-benchmark | HIGH | Official PyPI (5.2.3), pandas docs, 5+ years production use |
| PyArrow | HIGH | Official PyPI (23.0.0), pandas 3.0 default backend, Apache Arrow official docs |
| numba | HIGH | Official PyPI (0.63.1), pandas official docs recommend for performance |
| xlsxwriter | HIGH | Official PyPI (3.2.9), multiple real-world benchmarks (9min→3min) |
| Cython | HIGH | Official PyPI (3.2.4), pandas/NumPy/SciPy use internally |
| hatch-cython | MEDIUM | Community plugin (active development), modern hatchling integration |
| py-spy | HIGH | Official PyPI (0.4.1), widely used in production Python profiling |

**All versions verified:** 2026-02-14 via PyPI official sources.

---

## Sources

### Benchmarking & Profiling
- [pytest-benchmark documentation](https://pytest-benchmark.readthedocs.io/)
- [pytest-benchmark PyPI](https://pypi.org/project/pytest-benchmark/)
- [How To Measure And Improve Code Efficiency with Pytest Benchmark](https://pytest-with-eric.com/pytest-best-practices/pytest-benchmark/)
- [py-spy GitHub](https://github.com/benfred/py-spy)
- [Why Is My Code So Slow? A Guide to Py-Spy Python Profiling](https://towardsdatascience.com/why-is-my-code-so-slow-a-guide-to-py-spy-python-profiling/)
- [pytest-codspeed PyPI](https://pypi.org/project/pytest-codspeed/)

### DataFrame Optimization
- [PyArrow Functionality — pandas 3.0.0 documentation](https://pandas.pydata.org/docs/user_guide/pyarrow.html)
- [Python Pandas Ditches NumPy for Speedier PyArrow](https://thenewstack.io/python-pandas-ditches-numpy-for-speedier-pyarrow/)
- [Utilizing PyArrow to improve pandas and Dask workflows](https://towardsdatascience.com/utilizing-pyarrow-to-improve-pandas-and-dask-workflows-2891d3d96d2b/)
- [Enhancing performance — pandas 3.0.0 documentation](https://pandas.pydata.org/docs/user_guide/enhancingperf.html)
- [Polars vs pandas: What's the Difference?](https://realpython.com/polars-vs-pandas/)
- [Pandas vs. Polars: A Comparative Analysis](https://www.kdnuggets.com/pandas-vs-polars-a-comparative-analysis-of-python-dataframe-libraries)

### Categorical Dtypes
- [Categorical data — pandas 3.0.0 documentation](https://pandas.pydata.org/docs/user_guide/categorical.html)
- [Mastering Memory Optimization for Pandas DataFrames](https://thinhdanggroup.github.io/pandas-memory-optimization/)
- [A Powerful Tool for Saving Memory Space in Pandas: Categorical Types](https://www.oreateai.com/blog/a-powerful-tool-for-saving-memory-space-in-pandas-a-detailed-explanation-of-categorical-types/f9f0bec66f58a7997336de0e70e49058)

### Excel Generation
- [Openpyxl vs XlsxWriter: The Ultimate Showdown](https://hive.blog/python/@geekgirl/openpyxl-vs-xlsxwriter-the-ultimate-showdown-for-excel-automation)
- [Optimizing Excel Report Generation: From OpenPyXL to XLSXWriter](https://mass-software-solutions.medium.com/optimizing-excel-report-generation-from-openpyxl-to-xlsmerge-processing-52-columns-200k-rows-5b5a03ecbcd4)
- [xlsxwriter PyPI](https://pypi.org/project/xlsxwriter/)

### Native Extensions
- [Cython official site](https://cython.org/)
- [Building Cython code — Cython 3.3.0 documentation](https://cython.readthedocs.io/en/latest/src/quickstart/build.html)
- [hatch-cython PyPI](https://pypi.org/project/hatch-cython/)
- [CythonExtensionsOnWindows Wiki](https://github.com/cython/cython/wiki/CythonExtensionsOnWindows)

### Numba
- [Fast vectorization with Numba](https://tutorial.xarray.dev/advanced/apply_ufunc/numba-vectorization.html)
- [Making Python extremely fast with Numba](https://medium.com/@mflova/making-python-extremely-fast-with-numba-advanced-deep-dive-2-3-f809b43f8300)
- [numba PyPI](https://pypi.org/project/numba/)
