# Architecture Integration for Performance Optimizations

**Project:** variantcentrifuge
**Researched:** 2026-02-14
**Context:** Subsequent milestone — adding performance optimizations to existing pipeline
**Confidence:** HIGH (based on existing codebase analysis and performance report)

---

## Executive Summary

Performance optimizations must integrate seamlessly with the existing **stage-based pipeline architecture** while preserving **clinical correctness** — the tool's primary constraint. This document provides the architectural blueprint for safely integrating 27+ identified optimizations across 7 integration points, with a phased build order that minimizes risk and maintains test coverage throughout.

**Key architectural constraints:**
- Stage-based pipeline passes `PipelineContext` between 40+ stages
- DataFrame flows through `context.current_dataframe` — most optimizations touch this
- Inheritance analyzer has 3 sequential passes that must remain deterministic
- All changes must produce byte-identical output (clinical requirement)
- Existing test suite (1035 tests) must pass throughout optimization work

**Recommended approach:** Build optimizations as **isolated utility modules** that stages call into, rather than modifying stage logic directly. This enables:
- Independent testing of each optimization
- Easy rollback if issues arise
- A/B comparison of old vs new implementations
- Gradual migration with feature flags

---

## Current Architecture Overview

### Pipeline Flow

```
CLI Entry (variantcentrifuge.cli:main)
    ↓
create_parser() → argparse.Namespace
    ↓
Config mapping → dict[str, Any]
    ↓
run_refactored_pipeline()
    ↓
PipelineRunner.run(stages, context)
    ↓
┌─────────────────────────────────────────────────────────────┐
│ Stage Execution (dependency-ordered, parallel where safe)   │
├─────────────────────────────────────────────────────────────┤
│ Phase 1: Setup (7 stages)                                   │
│   - ConfigurationLoadingStage                               │
│   - PedigreeLoadingStage                                    │
│   - PhenotypeLoadingStage                                   │
│   - ScoringConfigLoadingStage                               │
│   - GeneBedCreationStage                                    │
│                                                             │
│ Phase 2: VCF Processing (3-4 stages)                        │
│   - VariantExtractionStage (bcftools)                       │
│   - SnpSiftFilterStage (SnpSift)                            │
│   - FieldExtractionStage (SnpSift extractFields)            │
│   - GenotypeReplacementStage (vectorized_replacer.py)       │
│                                                             │
│ Phase 3: DataFrame Analysis (6 stages)                      │
│   - DataFrameLoadingStage ← OPTIMIZATION POINT 1            │
│   - InheritanceAnalysisStage ← OPTIMIZATION POINT 2         │
│   - VariantScoringStage                                     │
│   - StatisticsGenerationStage ← OPTIMIZATION POINT 3        │
│   - GeneBurdenAnalysisStage ← OPTIMIZATION POINT 4          │
│   - ClinVarPM5Stage                                         │
│                                                             │
│ Phase 4: Output (5 stages)                                  │
│   - TSVOutputStage                                          │
│   - ExcelReportStage ← OPTIMIZATION POINT 5                 │
│   - HTMLReportStage                                         │
│   - IGVReportStage                                          │
│   - ArchiveStage                                            │
└─────────────────────────────────────────────────────────────┘
    ↓
PipelineContext (final state)
```

### PipelineContext Data Flow

```
PipelineContext carries data through stages:

Stage 1-7 (Setup):
  context.config = {...}
  context.pedigree_data = {...}
  context.vcf_samples = ["HG002", "HG003", ...]
  context.gene_bed_file = Path("genes.bed")

Stage 8-11 (Processing):
  context.extracted_vcf = Path("extracted.vcf.gz")
  context.filtered_vcf = Path("filtered.vcf.gz")
  context.extracted_tsv = Path("fields.tsv")
  context.genotype_replaced_tsv = Path("replaced.tsv")

Stage 12 (DataFrame Loading):
  context.data = Path("replaced.tsv")  # file path
  ↓
  df = pd.read_csv(context.data, sep="\t", dtype=str)  ← SLOW
  ↓
  context.current_dataframe = df  # DataFrame in memory

Stage 13-18 (Analysis):
  df = context.current_dataframe
  # Each stage modifies df in place or creates new df
  df = analyze_inheritance(df, context.pedigree_data, ...)
  df = apply_scoring(df, context.scoring_config, ...)
  context.current_dataframe = df  # Updated DataFrame

Stage 19-23 (Output):
  df = context.current_dataframe
  df.to_csv(context.final_output_path)  # TSV output

  # Excel stage RE-READS from disk (INEFFICIENT)
  df_excel = pd.read_csv(context.final_output_path, ...)  ← REDUNDANT
  df_excel.to_excel(xlsx_path)
```

**Critical insight:** Most optimizations target the DataFrame after it enters `context.current_dataframe` at Stage 12. This is the integration point.

---

## Integration Points

### 1. DataFrame Loading (DataFrameLoadingStage)

**File:** `variantcentrifuge/stages/analysis_stages.py:DataFrameLoadingStage._process()`

**Current implementation:**
```python
df = pd.read_csv(tsv_file, sep="\t", dtype=str)
context.current_dataframe = df
```

**Optimizations to integrate:**
- PyArrow CSV engine (10-40x faster I/O)
- Categorical dtypes for low-cardinality columns (70-97% memory reduction)
- Numeric downcasting (50-87% memory on numeric columns)
- Pre-parsed GT column (eliminates re-parsing in gene burden)

**New architecture:**

```
DataFrameLoadingStage._process()
    ↓
df = load_optimized_dataframe(tsv_file, config)  ← NEW UTILITY
    ↓
    ├─ pd.read_csv(engine="pyarrow") if available, else "c"
    ├─ apply_genomic_dtypes(df)  ← NEW: categorical + numeric
    ├─ preparse_gt_column(df) if gene_burden_enabled  ← NEW
    └─ return df
    ↓
context.current_dataframe = df
```

**New component required:**
- `variantcentrifuge/optimization/dataframe_loader.py`
  - `load_optimized_dataframe(path, config) -> pd.DataFrame`
  - `apply_genomic_dtypes(df) -> pd.DataFrame`
  - `preparse_gt_column(df) -> pd.DataFrame`

**Risk:** MEDIUM
- PyArrow engine has edge cases with malformed CSV
- Categorical dtypes require `observed=True` on all groupby (12 call sites)
- Pre-parsed GT changes gene burden input format

**Mitigation:**
- Feature flag: `--enable-pyarrow` (default: False initially)
- Feature flag: `--enable-categorical-dtypes` (default: False)
- Comprehensive tests for each optimization independently
- Validate byte-identical output with/without optimizations

**Testing strategy:**
```python
# tests/unit/test_dataframe_loader.py
def test_load_optimized_vs_standard_identical():
    """Verify optimized loader produces identical DataFrame."""
    df_standard = pd.read_csv(tsv_file, sep="\t", dtype=str)
    df_optimized = load_optimized_dataframe(tsv_file, config={})
    # Compare values, allowing dtype differences
    pd.testing.assert_frame_equal(df_standard, df_optimized, check_dtype=False)

def test_categorical_dtypes_memory_reduction():
    """Measure memory savings from categorical dtypes."""
    with MemoryBudget(max_mb=100) as mem:
        df = load_optimized_dataframe(large_tsv, {"enable_categorical": True})
    assert mem.peak_mb < 50  # Expect ~50% reduction
```

---

### 2. Inheritance Analysis (InheritanceAnalysisStage)

**File:** `variantcentrifuge/stages/analysis_stages.py:InheritanceAnalysisStage._process()`

**Current implementation:**
```python
df = context.current_dataframe
df = analyze_inheritance(df, pedigree_data, vcf_samples)
context.current_dataframe = df
```

**Calls into:** `variantcentrifuge/inheritance/analyzer.py:analyze_inheritance()`

**3-pass structure (MUST PRESERVE):**
```python
# Pass 1: Pattern deduction (CRITICAL BOTTLENECK)
df["_inheritance_patterns"] = df.apply(
    lambda row: deduce_patterns_for_variant(row.to_dict(), ...),
    axis=1  # ← 50K Python function calls, each with 5K-key dict
)

# Pass 2: Compound heterozygous detection
for gene, gene_df in df.groupby("GENE"):
    comp_het_results = analyze_gene_for_compound_het(gene_df, ...)

# Pass 3: Apply comp_het results + prioritization
for idx, row in df.iterrows():  # ← 2 sequential iterrows passes
    # Apply Pass 2 results
    # Prioritize patterns
```

**Optimizations to integrate:**
- Replace `apply(axis=1)` with `itertuples()` (10x faster) or full vectorization (100-740x)
- Replace `iterrows()` with `itertuples()` in Pass 3 (10x faster)
- Use vectorized comp_het (already exists: `comp_het_vectorized.py`, make default)
- Add `observed=True` to `groupby("GENE")` (prevents 3500x slowdown)

**New architecture:**

```
InheritanceAnalysisStage._process()
    ↓
df = context.current_dataframe
    ↓
df = analyze_inheritance_optimized(df, pedigree_data, vcf_samples, config)
    ↓
    Pass 1: Pattern deduction
    ├─ if config.get("use_vectorized_inheritance"):
    │     patterns = deduce_patterns_vectorized(df, pedigree, samples)  ← NEW
    └─ else:
          patterns = [deduce_patterns_for_variant(row._asdict(), ...)
                      for row in df.itertuples(index=False)]  ← OPTIMIZED

    Pass 2: Compound het
    ├─ Use comp_het_vectorized.py (already exists, make default)
    └─ Add observed=True, sort=False to groupby

    Pass 3: Apply + prioritize
    └─ Replace iterrows with vectorized map/apply operations  ← NEW
    ↓
context.current_dataframe = df
```

**New components required:**
- `variantcentrifuge/inheritance/deducer_vectorized.py`
  - `deduce_patterns_vectorized(df, pedigree, samples) -> list[dict]`
  - Operates on NumPy arrays extracted from sample columns
- `variantcentrifuge/inheritance/analyzer_optimized.py`
  - `analyze_inheritance_optimized(df, ...) -> pd.DataFrame`
  - Wraps old analyzer, adds feature flags for new implementations

**Risk:** HIGH
- Inheritance logic is complex with many edge cases
- Must produce identical results (clinical correctness)
- Vectorization requires reimplementing deduction logic

**Mitigation:**
- Build `deducer_vectorized.py` independently, test exhaustively
- Implement as opt-in: `--use-vectorized-inheritance` flag
- Add determinism tests comparing old vs new on real data
- Keep old implementation as fallback indefinitely

**Testing strategy:**
```python
# tests/unit/test_inheritance_optimized.py
@pytest.mark.parametrize("variant_count", [100, 1000, 10000])
def test_inheritance_itertuples_vs_apply(variant_count, benchmark):
    """Verify itertuples optimization produces identical results."""
    df = make_variant_df(variant_count)

    # Old implementation
    df1 = df.copy()
    df1 = analyze_inheritance(df1, pedigree, samples)

    # New implementation
    df2 = df.copy()
    df2 = analyze_inheritance_optimized(df2, pedigree, samples,
                                         {"use_itertuples": True})

    # Results must be byte-identical
    pd.testing.assert_frame_equal(df1, df2)

def test_vectorized_inheritance_correctness():
    """Verify vectorized implementation matches reference."""
    # Use known test cases with expected outputs
    df = load_test_variants("giab_trio_brca1.tsv")

    df_ref = analyze_inheritance(df.copy(), pedigree, samples)
    df_vec = analyze_inheritance_optimized(df.copy(), pedigree, samples,
                                            {"use_vectorized_inheritance": True})

    # Compare inheritance pattern columns
    pd.testing.assert_series_equal(df_ref["Inheritance_Pattern"],
                                    df_vec["Inheritance_Pattern"])
```

---

### 3. Statistics Generation (StatisticsGenerationStage)

**File:** `variantcentrifuge/stages/analysis_stages.py:StatisticsGenerationStage._process()`

**Current implementation:**
```python
from ..stats import (
    compute_gene_level_stats,
    compute_impact_summary,
    compute_variant_type_summary,
)

stats = {}
stats["gene_stats"] = compute_gene_level_stats(df)
stats["impact_summary"] = compute_impact_summary(df)
stats["variant_type_summary"] = compute_variant_type_summary(df)
context.statistics = stats
```

**Calls into:** `variantcentrifuge/stats.py`

**Current bottleneck:**
```python
# stats.py:117, 151, 178
df.groupby("GENE")  # Missing observed=True, sort=False
df.groupby(["GENE", "IMPACT"])
df.groupby(["GENE", "EFFECT"])
```

**Optimization to integrate:**
- Add `observed=True, sort=False` to all 12 groupby calls

**New architecture:**

```
No new components needed — direct fix in stats.py

stats.py:compute_gene_level_stats()
    ↓
gene_stats = df.groupby("GENE", observed=True, sort=False).agg({...})
                                 ↑↑↑↑↑↑↑↑↑↑↑↑  ↑↑↑↑↑↑↑↑↑↑↑
                                 ADD THIS      AND THIS
```

**Risk:** LOW
- Trivial code change (search-and-replace)
- No logic changes
- `observed=True` only matters when GENE is categorical (Integration Point 1)

**Mitigation:**
- Add regression test verifying stats are identical with/without categoricals
- Test with both string and categorical GENE column

**Testing strategy:**
```python
# tests/unit/test_stats.py
def test_stats_with_categorical_gene():
    """Verify stats work correctly with categorical GENE column."""
    df = make_variant_df(1000)
    df["GENE"] = df["GENE"].astype("category")

    stats = compute_gene_level_stats(df)
    assert len(stats) > 0
    assert "BRCA1" in stats.index

def test_stats_identical_string_vs_categorical():
    """Verify stats are identical regardless of GENE dtype."""
    df = make_variant_df(1000)

    df_str = df.copy()
    df_str["GENE"] = df_str["GENE"].astype(str)
    stats_str = compute_gene_level_stats(df_str)

    df_cat = df.copy()
    df_cat["GENE"] = df_cat["GENE"].astype("category")
    stats_cat = compute_gene_level_stats(df_cat)

    pd.testing.assert_frame_equal(stats_str, stats_cat)
```

---

### 4. Gene Burden Analysis (GeneBurdenAnalysisStage)

**File:** `variantcentrifuge/stages/analysis_stages.py:GeneBurdenAnalysisStage._process()`

**Current implementation:**
```python
from ..gene_burden import perform_gene_burden_analysis

burden_results = perform_gene_burden_analysis(df, context)
context.gene_burden_results = burden_results
```

**Calls into:** `variantcentrifuge/gene_burden.py:compute_gene_burden()`

**Current bottleneck (CRITICAL):**
```python
# gene_burden.py:207
for gene, gene_df in df.groupby("GENE"):  # Missing observed=True
    # gene_burden.py:220-249 — DEAD CODE
    for _, row in gene_df.iterrows():  # Parse GT column
        gt_value = str(row.get("GT", ""))
        for sample_entry in gt_value.split(";"):
            # Parse "Sample1(0/1);Sample2(0/0);..."
            # BUT RESULTS ARE NEVER USED (lines 251-256 use pre-computed counts)
```

**Performance report finding:** Lines 220-249 are provably dead code. Entire nested loop produces values that are immediately discarded.

**Optimizations to integrate:**
- Remove dead code (lines 220-249)
- Add `observed=True, sort=False` to groupby
- Use pre-parsed GT column from Integration Point 1 (if available)

**New architecture:**

```
GeneBurdenAnalysisStage._process()
    ↓
df = context.current_dataframe
    ↓
burden_results = perform_gene_burden_analysis(df, context)
    ↓
    gene_burden.py:compute_gene_burden()
        ↓
        for gene, gene_df in df.groupby("GENE", observed=True, sort=False):
            # REMOVE lines 220-249 entirely (dead code)
            # Use existing aggregated counts (lines 251-256)
            case_with_var = gene_df[case_mask].shape[0]
            control_with_var = gene_df[control_mask].shape[0]
            # Fisher's exact test, FDR correction, etc.
    ↓
context.gene_burden_results = burden_results
```

**Risk:** VERY LOW
- Removing dead code has zero functional impact
- Adding `observed=True` is safe (only matters with categoricals)
- No logic changes

**Mitigation:**
- Add regression test verifying burden results identical before/after
- Test with both string and categorical GENE column

**Testing strategy:**
```python
# tests/unit/test_gene_burden.py
def test_gene_burden_dead_code_removed():
    """Verify gene burden works after removing dead GT parsing loop."""
    df = load_test_variants("gene_burden_test.tsv")

    burden_old = compute_gene_burden_with_dead_code(df, config)
    burden_new = compute_gene_burden(df, config)

    pd.testing.assert_frame_equal(burden_old, burden_new)

def test_gene_burden_with_preparsed_gt():
    """Verify gene burden can use pre-parsed GT column."""
    df = load_test_variants("gene_burden_test.tsv")
    df = preparse_gt_column(df)  # From Integration Point 1

    burden = compute_gene_burden(df, config)
    assert "p_value" in burden.columns
    assert all(burden["p_value"] >= 0) and all(burden["p_value"] <= 1)
```

---

### 5. Excel Output (ExcelReportStage)

**File:** `variantcentrifuge/stages/output_stages.py:ExcelReportStage._process()`

**Current implementation:**
```python
from ..converter import convert_to_excel, finalize_excel_file

# RE-READS TSV FROM DISK (inefficient)
xlsx_path = convert_to_excel(str(tsv_output), context.config)
finalize_excel_file(xlsx_path, context.config)
```

**Calls into:** `variantcentrifuge/converter.py:tsv_to_xlsx()`

**Current bottleneck:**
```python
# converter.py:55
def tsv_to_xlsx(tsv_file, ...):
    df = pd.read_csv(tsv_file, sep="\t", low_memory=False)  ← RE-READ
    df.to_excel(xlsx_file, index=False)  # openpyxl is slow
```

**Optimizations to integrate:**
- Pass DataFrame directly from context (eliminate disk read)
- Use xlsxwriter for initial write (2-5x faster than openpyxl)
- Keep openpyxl only for finalization (hyperlinks, freeze panes)

**New architecture:**

```
ExcelReportStage._process()
    ↓
df = context.current_dataframe  ← USE IN-MEMORY DATAFRAME
    ↓
xlsx_path = convert_dataframe_to_excel(df, context.config, output_path)  ← NEW
    ↓
    converter.py:convert_dataframe_to_excel()
        ↓
        # Use xlsxwriter for fast initial write
        df.to_excel(xlsx_path, engine="xlsxwriter", index=False)
        ↓
        # Use openpyxl for finalization (hyperlinks, etc.)
        finalize_excel_file(xlsx_path, config)
    ↓
context.add_report_path("excel", xlsx_path)
```

**New component required:**
- Modify `variantcentrifuge/converter.py`
  - Add `convert_dataframe_to_excel(df, config, output_path) -> Path`
  - Keep `tsv_to_xlsx()` for backward compatibility

**Risk:** LOW
- DataFrame is already in memory
- xlsxwriter is well-tested for basic writes
- openpyxl finalization unchanged

**Mitigation:**
- Feature flag: `--excel-use-dataframe` (default: False initially)
- Regression test comparing TSV-based vs DataFrame-based Excel output
- Verify hyperlinks, freeze panes, auto-filter still work

**Testing strategy:**
```python
# tests/unit/test_excel_output.py
def test_excel_from_dataframe_vs_tsv():
    """Verify DataFrame-based Excel identical to TSV-based."""
    df = make_variant_df(100)
    tsv_path = tmp_path / "test.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)

    # Old approach: read from TSV
    xlsx_old = convert_to_excel(str(tsv_path), config)

    # New approach: use DataFrame
    xlsx_new = convert_dataframe_to_excel(df, config, tmp_path / "new.xlsx")

    # Compare Excel files (openpyxl)
    wb_old = load_workbook(xlsx_old)
    wb_new = load_workbook(xlsx_new)

    # Verify same data
    assert wb_old["Results"].max_row == wb_new["Results"].max_row
    assert wb_old["Results"]["A1"].value == wb_new["Results"]["A1"].value

def test_xlsxwriter_performance(benchmark):
    """Benchmark xlsxwriter vs openpyxl."""
    df = make_variant_df(10000)

    # Measure openpyxl
    time_openpyxl = benchmark.pedantic(
        lambda: df.to_excel("test_openpyxl.xlsx", engine="openpyxl"),
        rounds=3, iterations=1
    )

    # Measure xlsxwriter
    time_xlsxwriter = benchmark.pedantic(
        lambda: df.to_excel("test_xlsxwriter.xlsx", engine="xlsxwriter"),
        rounds=3, iterations=1
    )

    # Expect 2-5x speedup
    assert time_openpyxl / time_xlsxwriter >= 2.0
```

---

### 6. Subprocess Pipe Fusion (Processing Stages)

**Files:**
- `variantcentrifuge/stages/processing_stages.py:SnpSiftFilterStage`
- `variantcentrifuge/filters.py:apply_snpsift_filter()`

**Current implementation:**
```python
# SnpSift filter writes temp file
subprocess.run(["SnpSift", "filter", expr, input_vcf], stdout=temp_vcf)
# Then separately compress
subprocess.run(["bgzip", temp_vcf])
# Then index
subprocess.run(["bcftools", "index", temp_vcf_gz])
```

**Optimization to integrate:**
- Pipe SnpSift output directly to bgzip (eliminates temp file)

**New architecture:**

```
SnpSiftFilterStage._process()
    ↓
filtered_vcf = apply_snpsift_filter_piped(variant_file, expr, config, output_file)
    ↓
    filters.py:apply_snpsift_filter_piped()
        ↓
        # Start SnpSift process
        snpsift = subprocess.Popen(["SnpSift", "filter", expr, input_vcf],
                                   stdout=subprocess.PIPE)
        # Pipe to bgzip
        with open(output_vcf_gz, "wb") as out:
            bgzip = subprocess.Popen(["bgzip", "-c"],
                                     stdin=snpsift.stdout,
                                     stdout=out)
        snpsift.stdout.close()  # Allow SIGPIPE
        bgzip.wait()
        snpsift.wait()
        # Index result
        subprocess.run(["bcftools", "index", output_vcf_gz])
    ↓
context.filtered_vcf = filtered_vcf
```

**Risk:** MEDIUM
- Subprocess piping can fail silently
- Need to check return codes for both processes
- SIGPIPE handling on Windows vs Linux

**Mitigation:**
- Feature flag: `--enable-pipe-fusion` (default: False)
- Comprehensive error checking (both return codes)
- Test on both Windows and Linux
- Keep old implementation as fallback

**Testing strategy:**
```python
# tests/unit/test_subprocess_piping.py
def test_pipe_fusion_vs_temp_file():
    """Verify piped output identical to temp-file approach."""
    vcf_in = "tests/fixtures/test.vcf.gz"
    filter_expr = "( AF < 0.01 )"

    # Old approach: temp file
    vcf_old = apply_snpsift_filter(vcf_in, filter_expr, config)

    # New approach: piped
    vcf_new = apply_snpsift_filter_piped(vcf_in, filter_expr, config)

    # Compare outputs (should be byte-identical)
    assert vcf_old.read_bytes() == vcf_new.read_bytes()

def test_pipe_fusion_error_handling():
    """Verify piped approach handles errors correctly."""
    vcf_in = "tests/fixtures/test.vcf.gz"
    filter_expr = "( INVALID SYNTAX )"

    with pytest.raises(subprocess.CalledProcessError):
        apply_snpsift_filter_piped(vcf_in, filter_expr, config)
```

---

### 7. Shared Utilities (New Module)

**New file:** `variantcentrifuge/optimization/__init__.py`

**Purpose:** Central location for all optimization utilities, keeping them isolated from core logic.

**Structure:**
```
variantcentrifuge/
    optimization/
        __init__.py
        dataframe_loader.py       # Integration Point 1
        dtype_optimizer.py        # Categorical + numeric dtypes
        gt_parser.py              # Pre-parse GT column
        groupby_helper.py         # Wrapper adding observed=True
        pipe_fusion.py            # Subprocess piping utilities
```

**Benefits:**
- All optimizations in one place
- Easy to enable/disable entire module
- Independent testing
- Clear separation from clinical logic

**Example usage:**
```python
# In stages/analysis_stages.py
from ..optimization import load_optimized_dataframe, apply_genomic_dtypes

class DataFrameLoadingStage(Stage):
    def _process(self, context):
        if context.config.get("enable_optimizations", False):
            df = load_optimized_dataframe(tsv_file, context.config)
        else:
            df = pd.read_csv(tsv_file, sep="\t", dtype=str)

        context.current_dataframe = df
        return context
```

---

## Build Order (Phased Approach)

### Phase 1: Quick Wins (Tier 1) — Week 1

**Goal:** Immediate impact, zero risk, no architectural changes.

| # | Change | Integration Point | Files Modified | Tests Added |
|---|--------|-------------------|----------------|-------------|
| 1 | Remove gene burden dead code | IP 4 | `gene_burden.py` | `test_gene_burden.py` |
| 2 | Add `observed=True` to all groupby | IP 3, 4 | `stats.py`, `gene_burden.py`, `analyzer.py` | `test_stats.py` |
| 3 | Add gc.collect() between stages | — | `runner.py` | None (memory test) |

**Validation:**
```bash
# Run full test suite
pytest -m "not slow and not integration"

# Run on GIAB trio
variantcentrifuge --vcf giab_trio.vcf.gz --genes BRCA1 BRCA2

# Compare output byte-for-byte
diff old_output.tsv new_output.tsv  # Should be identical
```

**Expected impact:**
- Gene burden: 20-30% faster (dead code removed)
- No memory/speed regression
- Zero functional changes

---

### Phase 2: DataFrame Loading Optimizations (Tier 1-2) — Week 2

**Goal:** Optimize I/O and memory without touching analysis logic.

| # | Change | Integration Point | Files Modified | Tests Added |
|---|--------|-------------------|----------------|-------------|
| 4 | PyArrow CSV engine | IP 1 | `optimization/dataframe_loader.py`, `analysis_stages.py` | `test_dataframe_loader.py` |
| 5 | Categorical dtypes | IP 1 | `optimization/dtype_optimizer.py` | `test_dtype_optimizer.py` |
| 6 | Numeric downcasting | IP 1 | `optimization/dtype_optimizer.py` | `test_dtype_optimizer.py` |

**New files:**
- `variantcentrifuge/optimization/__init__.py`
- `variantcentrifuge/optimization/dataframe_loader.py`
- `variantcentrifuge/optimization/dtype_optimizer.py`
- `tests/unit/test_dataframe_loader.py`
- `tests/unit/test_dtype_optimizer.py`

**Feature flags (in config or CLI):**
```python
config = {
    "enable_pyarrow": False,  # Default off until validated
    "enable_categorical_dtypes": False,
    "enable_numeric_downcasting": False,
}
```

**Validation:**
```bash
# Test with each optimization independently
variantcentrifuge --enable-pyarrow --vcf test.vcf.gz
variantcentrifuge --enable-categorical-dtypes --vcf test.vcf.gz
variantcentrifuge --enable-numeric-downcasting --vcf test.vcf.gz

# Test with all optimizations
variantcentrifuge --enable-pyarrow --enable-categorical-dtypes \
                  --enable-numeric-downcasting --vcf test.vcf.gz

# Memory benchmark
pytest tests/performance/test_bench_dataframe_ops.py --benchmark-only
```

**Expected impact:**
- CSV read: 10-40x faster (PyArrow)
- Memory: 50-70% reduction (categoricals + downcasting)
- Zero functional changes (output identical)

---

### Phase 3: Inheritance Itertuples (Tier 2) — Week 3

**Goal:** Low-risk inheritance speedup without full vectorization.

| # | Change | Integration Point | Files Modified | Tests Added |
|---|--------|-------------------|----------------|-------------|
| 7 | Replace `apply(axis=1)` with itertuples | IP 2 | `inheritance/analyzer.py` | `test_inheritance_optimized.py` |
| 8 | Replace iterrows with itertuples (Pass 3) | IP 2 | `inheritance/analyzer.py` | `test_inheritance_optimized.py` |
| 9 | Make comp_het_vectorized default | IP 2 | `inheritance/analyzer.py` | Existing tests |

**Modified files:**
- `variantcentrifuge/inheritance/analyzer.py` (update Pass 1 and Pass 3)
- `tests/unit/test_inheritance_optimized.py` (new)

**Validation:**
```bash
# Determinism test (output must be identical)
for i in {1..10}; do
    variantcentrifuge --vcf test.vcf.gz --seed $i > output_$i.tsv
done
md5sum output_*.tsv  # All should match

# Compare old vs new implementation
pytest tests/unit/test_inheritance_optimized.py -v

# Performance benchmark
pytest tests/performance/test_bench_inheritance.py --benchmark-only
```

**Expected impact:**
- Inheritance Pass 1: 10x faster (apply → itertuples)
- Inheritance Pass 3: 10x faster (iterrows → itertuples)
- Combined: 5-10x inheritance speedup
- Output: byte-identical

---

### Phase 4: Excel DataFrame Pass-through (Tier 2) — Week 4

**Goal:** Eliminate redundant disk I/O in Excel generation.

| # | Change | Integration Point | Files Modified | Tests Added |
|---|--------|-------------------|----------------|-------------|
| 10 | Pass DataFrame to Excel stage | IP 5 | `converter.py`, `output_stages.py` | `test_excel_output.py` |
| 11 | Use xlsxwriter for initial write | IP 5 | `converter.py` | `test_excel_output.py` |

**Modified files:**
- `variantcentrifuge/converter.py` (add `convert_dataframe_to_excel()`)
- `variantcentrifuge/stages/output_stages.py` (use DataFrame directly)
- `tests/unit/test_excel_output.py` (new)

**Feature flag:**
```python
config = {"excel_use_dataframe": False}  # Default off initially
```

**Validation:**
```bash
# Compare Excel outputs
variantcentrifuge --excel-use-dataframe --vcf test.vcf.gz
openpyxl_compare old_output.xlsx new_output.xlsx

# Verify hyperlinks still work
pytest tests/unit/test_excel_output.py::test_hyperlinks -v

# Performance benchmark
pytest tests/performance/test_bench_excel.py --benchmark-only
```

**Expected impact:**
- Excel generation: 2-5x faster (xlsxwriter + no re-read)
- Eliminates redundant 1-2 GB disk read
- Output: functionally identical (hyperlinks, formatting preserved)

---

### Phase 5: Subprocess Pipe Fusion (Tier 3) — Week 5 (Optional)

**Goal:** Eliminate intermediate VCF temp files in processing.

| # | Change | Integration Point | Files Modified | Tests Added |
|---|--------|-------------------|----------------|-------------|
| 12 | Pipe SnpSift → bgzip | IP 6 | `filters.py`, `processing_stages.py` | `test_subprocess_piping.py` |

**Modified files:**
- `variantcentrifuge/filters.py` (add `apply_snpsift_filter_piped()`)
- `variantcentrifuge/stages/processing_stages.py` (use piped version)
- `tests/unit/test_subprocess_piping.py` (new)

**Feature flag:**
```python
config = {"enable_pipe_fusion": False}
```

**Validation:**
```bash
# Compare VCF outputs (byte-identical)
variantcentrifuge --enable-pipe-fusion --vcf test.vcf.gz
bcftools isec -c none old.vcf.gz new.vcf.gz  # Should be identical

# Test error handling
pytest tests/unit/test_subprocess_piping.py::test_error_handling -v
```

**Expected impact:**
- Filtering: 10-30s faster (eliminates temp I/O)
- Disk space: saves 1-2 GB intermediate file
- Output: byte-identical

---

### Phase 6: Vectorized Inheritance (Tier 3) — Week 6+ (Advanced)

**Goal:** Full vectorization of inheritance deduction (100-740x speedup).

**WARNING:** High complexity, high risk. Only pursue if Phases 1-5 are insufficient.

| # | Change | Integration Point | Files Modified | Tests Added |
|---|--------|-------------------|----------------|-------------|
| 13 | Vectorize Pass 1 deduction | IP 2 | `inheritance/deducer_vectorized.py` (new) | Extensive |

**New files:**
- `variantcentrifuge/inheritance/deducer_vectorized.py`
- `tests/unit/test_deducer_vectorized.py`
- `tests/integration/test_inheritance_vectorized_e2e.py`

**Feature flag:**
```python
config = {"use_vectorized_inheritance": False}  # Opt-in only
```

**Validation:**
```bash
# Extensive testing with known truth sets
pytest tests/unit/test_deducer_vectorized.py -v --tb=short

# Compare with reference implementation on real data
pytest tests/integration/test_inheritance_vectorized_e2e.py -v

# Benchmark on large dataset
pytest tests/performance/test_bench_inheritance.py \
    --benchmark-compare --benchmark-autosave
```

**Expected impact:**
- Inheritance Pass 1: 100-740x faster (vectorized NumPy)
- Overall pipeline: 2-3x faster (if inheritance is bottleneck)
- Output: must be byte-identical (extensive validation required)

**Risk mitigation:**
- Keep old implementation indefinitely as reference
- Extensive property-based testing (hypothesis)
- Test on all inheritance patterns (de novo, AR, AD, X-linked, mitochondrial)
- Validate on real clinical data before production use

---

## Data Flow Changes Summary

### Before Optimizations

```
VCF → bcftools → SnpSift → extractFields → TSV
                                            ↓
                                    pd.read_csv(dtype=str)
                                            ↓
                                    DataFrame (all string objects)
                                            ↓
                        ┌───────────────────┴───────────────────┐
                        ↓                                       ↓
                analyze_inheritance                    gene_burden (re-parse GT)
                (apply + iterrows)                     (iterrows nested loop)
                        ↓                                       ↓
                    scoring                                statistics
                        ↓                                       ↓
                        └───────────────────┬───────────────────┘
                                            ↓
                                    df.to_csv(TSV)
                                            ↓
                                    pd.read_csv(TSV)  ← RE-READ
                                            ↓
                                    df.to_excel(openpyxl)
```

### After Optimizations (All Phases)

```
VCF → bcftools → SnpSift | bgzip → extractFields → TSV
                    ↑                                 ↓
                    └─ PIPED (no temp file)   pd.read_csv(engine="pyarrow")
                                                      ↓
                                    apply_genomic_dtypes() (categorical + numeric)
                                                      ↓
                                    preparse_gt_column() (if gene burden enabled)
                                                      ↓
                                DataFrame (optimized dtypes, pre-parsed GT)
                                                      ↓
                        ┌───────────────────┴───────────────────┐
                        ↓                                       ↓
                analyze_inheritance                    gene_burden (use pre-parsed GT)
                (itertuples/vectorized)                (no re-parsing, observed=True)
                        ↓                                       ↓
                    scoring                                statistics
                        ↓                                (observed=True on groupby)
                        └───────────────────┬───────────────────┘
                                            ↓
                                    df.to_csv(TSV)
                                            ↓
                                    (in-memory DataFrame)  ← NO RE-READ
                                            ↓
                                    df.to_excel(xlsxwriter + openpyxl)
```

**Key changes:**
1. Pipe fusion eliminates temp VCF
2. PyArrow + optimized dtypes reduce memory 50-70%
3. Pre-parsed GT eliminates redundant string parsing
4. itertuples/vectorized speeds up iteration 10-740x
5. `observed=True` prevents categorical slowdowns
6. DataFrame pass-through eliminates redundant disk I/O

---

## Risk Assessment by Phase

| Phase | Risk Level | Rollback Strategy |
|-------|-----------|------------------|
| 1 (Quick wins) | VERY LOW | Trivial (dead code removal, flag additions) |
| 2 (DataFrame loading) | LOW | Feature flags disable all optimizations independently |
| 3 (Inheritance itertuples) | MEDIUM | Keep old implementation, feature flag to revert |
| 4 (Excel pass-through) | LOW | Feature flag reverts to TSV-based approach |
| 5 (Pipe fusion) | MEDIUM | Feature flag reverts to temp files |
| 6 (Vectorized inheritance) | HIGH | Keep old implementation, opt-in only, extensive validation |

**Overall risk mitigation:**
- Every optimization is feature-flagged
- Old implementations retained indefinitely
- Comprehensive test coverage for each optimization
- Byte-identical output validation
- Phased rollout (one phase at a time)

---

## Testing Strategy

### Unit Tests (per optimization)

```
tests/unit/
    test_dataframe_loader.py        # Phase 2
    test_dtype_optimizer.py         # Phase 2
    test_inheritance_optimized.py   # Phase 3
    test_excel_output.py            # Phase 4
    test_subprocess_piping.py       # Phase 5
    test_deducer_vectorized.py      # Phase 6
```

### Integration Tests (E2E validation)

```
tests/integration/
    test_pipeline_with_optimizations.py
        - Run full pipeline with all optimizations enabled
        - Compare output to baseline (byte-identical)
        - Verify all stages complete successfully
```

### Performance Benchmarks

```
tests/performance/
    test_bench_dataframe_ops.py     # Phase 2
    test_bench_inheritance.py       # Phase 3, 6
    test_bench_excel.py             # Phase 4
    test_bench_gene_burden.py       # Phase 1
```

### Regression Tests (output validation)

```python
# tests/regression/test_output_identity.py
def test_optimizations_preserve_output():
    """Verify optimized pipeline produces identical output."""

    # Baseline: no optimizations
    result_baseline = run_pipeline(vcf_file, config={})

    # With all optimizations
    result_optimized = run_pipeline(vcf_file, config={
        "enable_pyarrow": True,
        "enable_categorical_dtypes": True,
        "enable_numeric_downcasting": True,
        "excel_use_dataframe": True,
        "enable_pipe_fusion": True,
    })

    # Compare TSV outputs (must be byte-identical)
    assert result_baseline.read_text() == result_optimized.read_text()
```

---

## Success Criteria

**Phase completion criteria:**
- [ ] All new unit tests pass
- [ ] All existing tests continue to pass
- [ ] Output validation: byte-identical results with/without optimizations
- [ ] Performance benchmarks show expected improvement
- [ ] Memory profiling shows expected reduction (if applicable)
- [ ] Code review completed
- [ ] Documentation updated

**Overall milestone completion:**
- [ ] All 6 phases complete (or Phases 1-5 if Phase 6 deferred)
- [ ] Full test suite passes (1035+ tests)
- [ ] E2E integration tests pass with all optimizations enabled
- [ ] Performance benchmarks show 2-3x overall speedup
- [ ] Memory usage reduced 50-70%
- [ ] Clinical correctness preserved (output byte-identical)

---

## Dependencies Between Optimizations

```
Phase 1 (Quick wins)
    ↓ (independent)
Phase 2 (DataFrame loading)
    ↓ (categorical dtypes enable observed=True benefits)
Phase 3 (Inheritance itertuples)
    ↓ (independent)
Phase 4 (Excel pass-through)
    ↓ (independent)
Phase 5 (Pipe fusion)
    ↓ (independent)
Phase 6 (Vectorized inheritance) — optional, can be deferred
```

**Critical path:** Phase 2 → Phase 3 (categorical dtypes make `observed=True` critical)

**Parallelizable:** Phases 4, 5, 6 can be developed independently after Phase 1-2 complete

---

## Architecture Patterns

### 1. Optimization Wrapper Pattern

```python
# Instead of modifying core functions directly:
def analyze_inheritance_optimized(df, pedigree, samples, config):
    """Wrapper that selects implementation based on config."""
    if config.get("use_vectorized_inheritance"):
        return analyze_inheritance_vectorized(df, pedigree, samples)
    elif config.get("use_itertuples"):
        return analyze_inheritance_itertuples(df, pedigree, samples)
    else:
        return analyze_inheritance(df, pedigree, samples)  # Original
```

**Benefits:**
- Old implementation always available
- Easy A/B testing
- Feature flags control behavior
- Zero risk to existing code

### 2. Utility Module Pattern

```python
# Isolated optimization utilities in variantcentrifuge/optimization/

from ..optimization import (
    load_optimized_dataframe,
    apply_genomic_dtypes,
    preparse_gt_column,
)

# Stages call into utilities:
df = load_optimized_dataframe(tsv_file, config)
```

**Benefits:**
- Optimizations isolated from clinical logic
- Easy to disable entire module
- Independent testing
- Clear code organization

### 3. Feature Flag Pattern

```python
# All optimizations behind config flags:
config = {
    "enable_pyarrow": False,
    "enable_categorical_dtypes": False,
    "enable_numeric_downcasting": False,
    "excel_use_dataframe": False,
    "enable_pipe_fusion": False,
    "use_vectorized_inheritance": False,
}

# Gradual rollout: enable one at a time, test, then enable next
```

**Benefits:**
- Safe rollout
- Easy rollback
- Per-optimization A/B testing
- Production-ready from day 1

---

## Conclusion

Performance optimizations integrate cleanly into the existing stage-based architecture through **7 integration points** across **6 build phases**. The phased approach ensures:

1. **Safety:** Every optimization is feature-flagged and independently testable
2. **Correctness:** Byte-identical output validation at every phase
3. **Rollback:** Old implementations retained indefinitely
4. **Gradual adoption:** One phase at a time, with full test coverage
5. **Minimal risk:** Phases 1-4 are low risk, Phase 5 is medium, Phase 6 is opt-in only

**Recommended execution order:** Phases 1 → 2 → 3 → 4, then evaluate. Phases 5-6 only if needed.

**Expected outcome after Phases 1-4 (4-5 weeks):**
- 2-3x faster overall pipeline
- 50-70% less memory usage
- Zero functional changes
- All 1035+ tests passing
- Production-ready optimizations

**Phase 6 (vectorized inheritance) is high-risk/high-reward and should only be pursued if:**
- Phases 1-5 are insufficient for performance needs
- Team has 2+ weeks for extensive validation
- Clinical validation resources are available
- Old implementation will be retained indefinitely as reference
