# Research Summary: Performance Optimization for variantcentrifuge v0.13.0

**Project:** variantcentrifuge (Clinical genomic variant analysis pipeline)
**Milestone:** v0.13.0 Performance Optimization
**Target:** 3-4x full pipeline speedup on large cohorts (5,125 samples, 22 GB VCFs)
**Research Date:** 2026-02-14
**Confidence:** HIGH

---

## Key Findings

This research covers performance optimization of a pandas-heavy clinical genomics pipeline. The codebase currently processes 22 GB VCFs with 5,125 samples across 2,500 genes, taking 90 minutes. The 3-4x speedup target is achievable through **incremental optimization** without architectural rewrites. Critical insight: the tool uses `dtype=str` everywhere as a defensive measure, and optimizations must preserve byte-identical output for clinical correctness.

Research identified 27+ specific optimizations across 7 integration points, organized into 3 tiers. Tier 1 optimizations (remove dead code, add `observed=True`, PyArrow engine, categorical dtypes) deliver 1.5-2x speedup in 1-2 weeks with low risk. Tier 2 (vectorize iterrows/apply, xlsxwriter, DataFrame pass-through) pushes to 3.6x in 2-4 weeks. Tier 3 (pipe fusion, Cython, full vectorization) reaches 5x but has higher complexity. The phased approach minimizes risk by keeping old implementations as feature-flagged fallbacks.

The performance analysis report revealed inheritance analysis consumes 40-60% of total time (dominated by `df.apply(axis=1)` with 5,000-key dicts per row), gene burden has 20-30% dead code (GT parsing loop results never used), and Excel generation re-reads multi-GB TSV from disk. These bottlenecks are addressable with proven pandas techniques rather than migration to Polars or GPU frameworks.

---

## Stack Additions

Add these dependencies to `pyproject.toml` for benchmarking and optimization:

**Core Performance Stack (required):**
- `pyarrow>=23.0` — Columnar backend for pandas; 2-3x faster CSV I/O, 70% memory reduction for strings
- `numba>=0.63` — JIT compilation for numeric kernels; 100-740x speedup on vectorized array operations
- `xlsxwriter>=3.2` — Fast Excel generation; 2-5x faster than openpyxl for write-only operations

**Benchmarking and Profiling (dev dependencies):**
- `pytest-benchmark>=5.2` — Regression tracking with pedantic mode for exact iteration control
- `py-spy>=0.4` — Production profiling with native C extension support for pandas/NumPy profiling

**Build System (optional, only if Cython extensions needed in Tier 3):**
- `hatch-cython>=0.5` — Integrates Cython with existing hatchling build system
- `Cython>=3.2` — Compiled extensions for hot paths (genotype string processing)

**Version Requirements:**
- Ensure `pandas>=2.0` (already present) to support PyArrow backend
- All tools ship pre-built wheels for Windows and Linux (x86-64, ARM64)

**NOT Recommended:**
- Polars (requires 3-6 week full rewrite, breaks scikit-learn/matplotlib integration)
- Dask/Ray (pipeline not memory-bound, ProcessPoolExecutor already handles parallelism)
- GPU acceleration (I/O and string-processing bound, not numeric-heavy)
- Rust extensions (Cython achieves 90% of Rust performance with 50% less complexity)

---

## Feature Categories

### Table Stakes (Expected, Zero-Risk)

| Feature | Impact | Effort | Rationale |
|---------|--------|--------|-----------|
| **Eliminate df.apply(axis=1)** | 40-60% time saved | Medium | Standard pandas optimization; apply(axis=1) is 10-100x slower than vectorization |
| **Replace iterrows() loops** | 10-13x speedup | Low | Three instances in inheritance analysis; itertuples is drop-in replacement |
| **Add observed=True to groupby** | Prevents 3500x slowdown | Trivial | 12 call sites; critical when categorical dtypes added |
| **Typed dtypes (not dtype=str)** | 50-80% memory reduction | Medium | Convert CHROM→category, POS→int64, AF→float64; enables NumPy vectorization |
| **PyArrow engine for CSV I/O** | 2-3x I/O speedup | Low | Industry standard since 2024; 16 pd.read_csv calls need `engine="pyarrow"` |
| **Remove dead code** | 20-30% of gene burden | Low | Lines 220-249 in gene_burden.py parse GT but results never used |

### Differentiators (Maximum Impact)

| Feature | Impact | Effort | Rationale |
|---------|--------|--------|-----------|
| **Vectorize inheritance deduction** | 100-740x on hottest path | High | Replace lambda row: deduce_patterns_for_variant(row.to_dict(), ...) with NumPy matrix ops |
| **xlsxwriter for Excel** | 2-5x Excel speedup | Low | openpyxl writes cells one-by-one; xlsxwriter streams 2-5x faster |
| **Pre-parse GT column once** | Eliminates re-parsing overhead | Medium | GT parsing happens per-gene; do once at load time |
| **Categorical dtypes** | 50-90% memory reduction | Low | CHROM (~25 values), IMPACT (4), FILTER (3) waste 90-99% memory as strings |
| **Pass DataFrame in memory** | Eliminates 2-10 GB disk read | Trivial | Excel stage re-reads TSV that's already in context.current_dataframe |
| **Pipe fusion for external tools** | Eliminates 2 temp VCFs | Medium | Current: bcftools → temp → SnpSift → temp → extractFields; use pipes instead |
| **Parquet/Feather intermediates** | 4-15x I/O speedup | Medium | Feather: 2.3s vs 15s CSV for 5M records; replace write-then-read temp TSV |
| **Numba JIT for hot loops** | 10-50x on numeric ops | Medium | Best for NumPy array operations; first-call compilation overhead, subsequent cached |

### Anti-Features (Explicitly Avoid)

- **Migrate to Polars** — 10-100x speedup but requires complete rewrite; pandas + optimizations hits target
- **GPU acceleration** — I/O-bound and logic-heavy pipeline, not numeric-heavy; transfer overhead negates benefit
- **SharedMemory for parallel stages** — Complexity not justified for small gene-level DataFrames (1-20 variants each)
- **Cython/Rust for all code** — Only justified for proven bottlenecks after other optimizations exhausted
- **Custom C extensions for VCF** — Use cyvcf2 (wraps htslib); don't reinvent bioinformatics wheel
- **Async I/O** — Genomic files are GB-scale; async helps network I/O, not local disk bandwidth-bound
- **In-memory databases (SQLite/DuckDB)** — Query layer overhead; pandas already in-memory and optimized
- **Modin** — Drop-in parallelized pandas but doesn't fix root cause (inefficient apply/iterrows)

---

## Architecture Integration

Optimizations integrate through 7 integration points across the existing stage-based pipeline:

**Integration Points:**

1. **DataFrameLoadingStage** (analysis_stages.py)
   - Add `load_optimized_dataframe()` utility
   - Apply PyArrow engine, categorical dtypes, pre-parsed GT column
   - Risk: MEDIUM (PyArrow edge cases, categorical requires observed=True everywhere)

2. **InheritanceAnalysisStage** (analysis_stages.py)
   - Replace `apply(axis=1)` with itertuples (10x) or vectorization (100-740x)
   - Replace Pass 3 iterrows with itertuples
   - Use `observed=True` on `groupby("GENE")`
   - Risk: HIGH (complex logic, must preserve clinical correctness)

3. **StatisticsGenerationStage** (analysis_stages.py)
   - Add `observed=True, sort=False` to all groupby calls
   - Risk: LOW (trivial change, no logic modifications)

4. **GeneBurdenAnalysisStage** (analysis_stages.py)
   - Remove dead code (lines 220-249)
   - Add `observed=True` to groupby
   - Use pre-parsed GT column if available
   - Risk: VERY LOW (dead code removal has zero functional impact)

5. **ExcelReportStage** (output_stages.py)
   - Pass DataFrame from context (no re-read)
   - Use xlsxwriter for initial write, openpyxl for finalization
   - Risk: LOW (DataFrame already in memory, xlsxwriter well-tested)

6. **Subprocess Pipe Fusion** (processing_stages.py)
   - Pipe SnpSift → bgzip (eliminates temp VCF)
   - Check all return codes (both processes)
   - Risk: MEDIUM (pipe errors can fail silently)

7. **Shared Utilities Module** (NEW: optimization/)
   - `dataframe_loader.py`, `dtype_optimizer.py`, `gt_parser.py`, `groupby_helper.py`, `pipe_fusion.py`
   - Isolates all optimizations from clinical logic
   - Risk: LOW (new code, doesn't modify existing stages directly)

**Build Order (Phased):**

1. **Phase 1 (Week 1): Quick Wins** — Remove dead code, add observed=True, gc.collect()
   - Impact: 30-40% speedup
   - Risk: VERY LOW

2. **Phase 2 (Week 2): DataFrame Loading** — PyArrow engine, categorical dtypes, numeric downcasting
   - Impact: 50-70% memory reduction, 2-3x I/O speedup
   - Risk: LOW (feature-flagged)

3. **Phase 3 (Week 3): Inheritance Itertuples** — Replace apply/iterrows with itertuples
   - Impact: 5-10x inheritance speedup
   - Risk: MEDIUM (keep old implementation as fallback)

4. **Phase 4 (Week 4): Excel Pass-through** — DataFrame to Excel, xlsxwriter
   - Impact: 2-5x Excel speedup, eliminates disk re-read
   - Risk: LOW (feature-flagged)

5. **Phase 5 (Week 5, Optional): Pipe Fusion** — SnpSift → bgzip piping
   - Impact: 10-30s faster, saves 1-2 GB temp file
   - Risk: MEDIUM (subprocess error handling)

6. **Phase 6 (Week 6+, Advanced): Vectorized Inheritance** — Full vectorization of Pass 1
   - Impact: 100-740x on hottest path, 2-3x overall
   - Risk: HIGH (extensive validation required, opt-in only)

**Expected Cumulative Impact:**
- Baseline: 90 min
- After Phase 1: 60 min (1.5x)
- After Phase 2: 25 min (3.6x)
- After Phase 4: 18 min (5.0x)

Target achieved at Phase 2-3 (2-4 weeks effort).

---

## Critical Pitfalls (Top 5)

### 1. Pandas 3.0 String Dtype Breaking Changes with PyArrow (CRITICAL)

**Risk:** PyArrow engine in pandas 3.0 defaults to `string[pyarrow]` instead of object dtype. Missing value semantics change from `np.nan` to `pd.NA`, breaking string comparisons (`pd.NA == "value"` returns `pd.NA`, not `False`). Genotype matching code like `if genotype == "0/1"` fails if genotype is `pd.NA`.

**Prevention:**
- Keep `dtype=str` when using PyArrow engine: `pd.read_csv(path, sep="\t", dtype=str, engine="pyarrow")`
- Test with pandas 3.0+ in CI
- Use `pd.isna(value) or value == ""` instead of `if value == ""`
- Regression test all string comparisons

**Detection:** Tests fail with `TypeError: boolean value of NA is ambiguous`, different variant counts between pandas 2.x and 3.0

**Affected:** #4 (PyArrow engine), #5 (categorical dtypes)

### 2. Categorical Groupby 3500x Slowdown Without observed=True (CRITICAL)

**Risk:** Converting columns to categorical dtype without adding `observed=True` to groupby creates groups for every possible category (e.g., 2,500 genes) even if only 3 genes have variants. Each iteration processes empty DataFrames, causing catastrophic slowdown.

**Prevention:**
- ALWAYS use `df.groupby("GENE", observed=True, sort=False)`
- Search-and-replace across entire codebase (12 call sites)
- Add lint rule to detect missing `observed=True`
- Regression test: assert elapsed < 5.0s for categorical gene burden

**Detection:** Tests hang or timeout, memory spikes to GB for small DataFrames, FutureWarning about observed=False

**Affected:** #3 (add observed=True), #5 (categorical dtypes)

### 3. Vectorization Changing Row Processing Order (CRITICAL)

**Risk:** Replacing `iterrows()` with vectorized operations can break order-dependent logic. Compound heterozygous detection expects sequential gene analysis; vectorization processes all rows simultaneously without sequential state.

**Prevention:**
- Audit for side effects before vectorizing (mutation of external state is red flag)
- Preserve order-dependent logic: `for gene, gene_df in df.groupby("GENE", sort=False):`
- Test determinism: run 5 times, assert identical output
- Add row checksums to output

**Detection:** Different output on repeated runs, compound het pairs missing in vectorized version

**Affected:** #1 (vectorize inheritance), #6 (replace iterrows), #11 (full vectorization)

### 4. Dead Code Removal Breaking Validation Logic (HIGH)

**Risk:** Lines 220-249 in gene_burden.py parse GT strings but results are never used. However, the parsing loop may serve as validation that GT format is correct. Removing it could allow malformed data to pass silently.

**Prevention:**
- Confirm dead code is truly unused (instrument with flag before removing)
- Extract validation if needed: `validate_gt_format(gt_string)` applied once, fail fast
- Add regression test for malformed GT input
- Check git blame and PR history to understand why code was written

**Detection:** Tests pass after removal but production crashes on edge cases, no ValueError on malformed GT

**Affected:** #2 (remove dead code)

### 5. xlsxwriter Incompatibility with openpyxl Post-Processing (HIGH)

**Risk:** Pipeline uses openpyxl to finalize Excel files (add hyperlinks, freeze panes, auto-filters). xlsxwriter is write-only and cannot open existing files for modification. Switching breaks finalization step.

**Prevention:**
- Do all formatting during initial write via xlsxwriter API (freeze_panes, autofilter, write_url)
- OR use hybrid: write with xlsxwriter, reopen with openpyxl for finalization
- Benchmark both engines before committing
- Test Excel output structure: assert hyperlinks, freeze panes, auto-filters present

**Detection:** Missing hyperlinks/freeze panes/auto-filters, InvalidFileException from openpyxl

**Affected:** #8 (xlsxwriter for Excel)

---

## Recommended Phase Structure

### Phase 1: Quick Wins (1-2 hours, 1.5x speedup)

**What:** Remove dead code, add `observed=True`, enable PyArrow engine

**Why First:** Zero risk, immediate 30-40% speedup, unblocks categorical dtypes

**Deliverables:**
- Remove gene_burden.py:220-249
- Add `observed=True, sort=False` to 12 groupby calls
- Add `gc.collect()` between stages in runner.py
- Fix temp file leak (gene_bed.py:153)

**Validation:** All tests pass, output byte-identical, gene burden 20-30% faster

**Research Flags:** None (standard patterns)

### Phase 2: Core Optimizations (2-3 days, 3.6x cumulative speedup)

**What:** DataFrame loading optimizations, itertuples, categorical dtypes, in-memory pass-through

**Why Second:** Low-medium risk, delivers target speedup, enables Phase 3+

**Deliverables:**
- Create `optimization/` module with dataframe_loader.py, dtype_optimizer.py
- Add feature flags: `--enable-pyarrow`, `--enable-categorical-dtypes`
- Replace iterrows with itertuples in inheritance analyzer
- Pass DataFrame directly to Excel stage

**Validation:** Output byte-identical, memory reduced 50-70%, inheritance 5-10x faster

**Research Flags:** Needs research for PyArrow edge cases with malformed CSV

### Phase 3: Strategic Enhancements (1-2 weeks, 5x cumulative speedup)

**What:** Pre-parse GT, xlsxwriter, pipe fusion, Parquet intermediates

**Why Third:** Pushes beyond target, larger refactoring, optional for 3-4x goal

**Deliverables:**
- Pre-parse GT column at load time
- Hybrid xlsxwriter + openpyxl for Excel
- Pipe SnpSift → bgzip → extractFields
- Replace TSV write-then-read with Parquet

**Validation:** Output functionally identical (allow 1e-9 float tolerance), Excel has hyperlinks, subprocess errors caught

**Research Flags:** Needs research for subprocess piping error handling on Windows vs Linux

### Defer to Post-MVP

**Numba JIT compilation** — Requires extensive vectorization first; only beneficial after Phase 2 complete

**Chunked processing** — Only needed if DataFrames exceed 50% system RAM; currently 22 GB → 2-10 GB DataFrame fits in 32 GB

**Polars migration** — Only if optimized pandas still 3x too slow; require compelling 10x+ improvement case

**Cython/Rust extensions** — Only for genotype replacement kernel if still >30% of pipeline time after Phase 3

---

## Confidence Assessment

| Area | Confidence | Source Quality | Notes |
|------|------------|----------------|-------|
| Stack (PyArrow, numba, xlsxwriter) | HIGH | Official PyPI, pandas docs, 5+ years production use | All versions verified 2026-02-14 |
| Features (table stakes vs differentiators) | HIGH | Performance analysis report, pandas official docs, real-world benchmarks | 27+ optimizations from actual bottleneck analysis |
| Architecture (7 integration points) | HIGH | Existing codebase analysis, stage-based pipeline structure | Clear separation via `optimization/` module |
| Pitfalls (pandas 3.0, categorical groupby) | HIGH | Official pandas docs, GitHub issues with reproducible examples | Critical pitfalls well-documented |
| Phase structure (1→2→3 ordering) | HIGH | Dependency analysis, risk assessment | Phase 1-2 achieves target (3.6x), Phase 3 is stretch goal (5x) |

**Gaps to Address:**
- PyArrow engine behavior with malformed CSV requires testing on real data
- Subprocess piping error handling needs platform-specific validation (Windows vs Linux)
- Benchmark stability on GitHub Actions (2.66% variance) requires ratio assertions instead of absolute thresholds

**Overall Confidence:** HIGH — Research is specific to planned optimizations with actionable prevention strategies for all critical pitfalls.

---

## Sources

**Technology Stack:**
- [PyArrow Functionality — pandas 3.0.0 documentation](https://pandas.pydata.org/docs/user_guide/pyarrow.html)
- [pytest-benchmark documentation](https://pytest-benchmark.readthedocs.io/)
- [py-spy GitHub](https://github.com/benfred/py-spy)
- [xlsxwriter PyPI](https://pypi.org/project/xlsxwriter/)
- [Cython official site](https://cython.org/)
- [hatch-cython PyPI](https://pypi.org/project/hatch-cython/)

**Performance Optimization:**
- [Enhancing performance — pandas 3.0.0 documentation](https://pandas.pydata.org/docs/user_guide/enhancingperf.html)
- [Pandas vectorization: faster code, slower code, bloated memory](https://pythonspeed.com/articles/pandas-vectorization/)
- [Tutorial: 700x speed improvement on Pandas with Vectorization](https://www.linkedin.com/pulse/tutorial-basic-vectorization-pandas-iterrows-apply-duc-lai-trung-minh-75d4c)
- [For the Love of God, Stop Using iterrows()](https://ryxcommar.com/2020/01/15/for-the-love-of-god-stop-using-iterrows/)

**Categorical Dtypes:**
- [Categorical data — pandas 3.0.0 documentation](https://pandas.pydata.org/docs/user_guide/categorical.html)
- [Be Careful When Using Pandas Groupby with Categorical Data Type](https://towardsdatascience.com/be-careful-when-using-pandas-groupby-with-categorical-data-type-a1d31f66b162/)
- [Pandas GroupBy at Speed: Pitfalls and Power Moves](https://medium.com/@2nick2patel2/pandas-groupby-at-speed-pitfalls-and-power-moves-8b27ca7ccc5a)

**Excel Generation:**
- [Openpyxl vs XlsxWriter: The Ultimate Showdown](https://hive.blog/python/@geekgirl/openpyxl-vs-xlsxwriter-the-ultimate-showdown-for-excel-automation)
- [Optimizing Excel Report Generation: From OpenPyXL to XLSXWriter](https://mass-software-solutions.medium.com/optimizing-excel-report-generation-from-openpyxl-to-xlsmerge-processing-52-columns-200k-rows-5b5a03ecbcd4)

**Pitfalls:**
- [pandas 3.0.0 whatsnew](https://pandas.pydata.org/docs/whatsnew/v3.0.0.html)
- [String dtype breaking changes](https://github.com/pandas-dev/pandas/issues/59328)
- [Categorical groupby slowdown](https://github.com/pandas-dev/pandas/issues/32976)
- [AMP/CAP NGS bioinformatics validation guidelines](https://www.sciencedirect.com/science/article/pii/S1525157817303732)

**Benchmarking:**
- [GitHub Actions performance stability](https://aakinshin.net/posts/github-actions-perf-stability/)
- [Practical performance tests](https://solidean.com/blog/2025/practical-performance-tests/)

---

## Ready for Roadmap Creation

This summary synthesizes findings from 4 parallel research streams (STACK, FEATURES, ARCHITECTURE, PITFALLS) into actionable guidance. The roadmapper agent can now structure phases with clear dependencies, risk levels, and validation criteria. All critical pitfalls have prevention strategies, and the phased approach enables incremental delivery with rollback safety.

**Key Recommendations for Roadmapper:**
1. Start with Phase 1 (trivial wins) to build momentum
2. Phase 2 delivers target speedup; Phase 3 is optional stretch goal
3. Every optimization behind feature flag with old implementation as fallback
4. Byte-identical output validation required at every phase
5. Cross-platform testing (Windows + Linux) mandatory
6. Test with pandas 3.0+ to catch string dtype issues early
