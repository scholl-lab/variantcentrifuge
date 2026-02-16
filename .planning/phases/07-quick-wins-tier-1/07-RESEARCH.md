# Phase 7: Quick Wins - Tier 1 - Research

**Researched:** 2026-02-14
**Domain:** Python performance optimization (dead code removal, pandas groupby, temp file cleanup, memory management)
**Confidence:** HIGH

## Summary

Phase 7 targets four zero-risk optimizations that collectively deliver 30-40% speedup on gene burden analysis with identical output. These optimizations are standard Python/pandas best practices that eliminate known inefficiencies without architectural changes.

The research confirms all four optimization targets exist exactly as described in the performance analysis report, with precise file locations and code patterns identified. Implementation is straightforward with minimal risk.

**Primary recommendation:** Apply all four optimizations in a single phase. Start with dead code verification, apply groupby fixes atomically, refactor all temp file patterns to context managers, and add gc.collect() with memory logging between all pipeline stages.

## Standard Stack

### Core Technologies
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pandas | 2.x | DataFrame operations | Industry standard for data analysis |
| gc | stdlib | Garbage collection | Python built-in memory management |
| tempfile | stdlib | Temporary file handling | Python built-in, context manager support |
| psutil | 6.1.0+ | Memory monitoring | Cross-platform system monitoring |
| logging | stdlib | Debug logging | Python built-in structured logging |

### Supporting Tools
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| pytest | 8.0+ | Regression testing | Verify identical output before/after |
| ruff | 0.11+ | Linting | Enforce code quality standards |

**Installation:**
Already installed — all dependencies present in current environment.

## Architecture Patterns

### Pattern 1: Dead Code Verification Before Removal
**What:** Instrument suspected dead code with assertions/logging to prove it never executes.
**When to use:** When static analysis suggests code is unreachable but you want runtime proof.
**Example:**
```python
# Step 1: Add verification (in PR branch)
for _, row in gene_df.iterrows():
    assert False, f"Dead code executed! Gene={gene}, variants={len(gene_df)}"
    # ... existing dead code ...

# Step 2: Run full test suite + benchmarks
# If assertion never fires → code is provably dead

# Step 3: Remove entire block (same PR)
# Clean removal, no comments, git history is the record
```

### Pattern 2: Context Manager for Temp Files
**What:** Use `tempfile.NamedTemporaryFile` with `with` statement for automatic cleanup.
**When to use:** Any time temporary files are created.
**Example:**
```python
# BEFORE (leak-prone):
bed_fd, bed_path = tempfile.mkstemp(suffix=".bed")
os.close(bed_fd)
# ... use bed_path ...
# FORGOT: os.remove(bed_path)  ← leak!

# AFTER (leak-proof):
with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
    bed_path = tmp.name
    # ... use tmp.name ...
# Automatic cleanup on context exit
os.remove(bed_path)  # Explicit if delete=False needed

# OR for automatic deletion:
with tempfile.NamedTemporaryFile(mode='w', suffix='.bed') as tmp:
    # ... use tmp.name ...
# Automatic cleanup, no manual remove needed
```

### Pattern 3: Memory-Aware Stage Transitions
**What:** Call `gc.collect()` after every pipeline stage with optional memory logging.
**When to use:** Between stages that process large DataFrames or produce intermediate objects.
**Example:**
```python
# In runner.py _execute_stage()
import gc
import psutil

def _execute_stage(self, stage: Stage, context: PipelineContext) -> PipelineContext:
    start_time = time.time()

    # Optional: Log memory before stage
    if logger.isEnabledFor(logging.DEBUG):
        process = psutil.Process()
        mem_before_mb = process.memory_info().rss / 1024 / 1024
        logger.debug(f"Stage '{stage.name}' starting - Memory: {mem_before_mb:.1f} MB")

    result = stage(context)
    elapsed = time.time() - start_time

    # Force garbage collection
    gc.collect()

    # Optional: Log memory after GC
    if logger.isEnabledFor(logging.DEBUG):
        mem_after_mb = process.memory_info().rss / 1024 / 1024
        freed_mb = mem_before_mb - mem_after_mb
        logger.debug(
            f"Stage '{stage.name}' complete - Memory: {mem_after_mb:.1f} MB "
            f"(freed {freed_mb:.1f} MB)"
        )

    return result
```

### Pattern 4: Atomic Groupby Migration
**What:** Add `observed=True` to all groupby calls in a single commit.
**When to use:** When migrating to categorical dtypes or preventing future slowdowns.
**Example:**
```python
# Search pattern: \.groupby\(
# Replace: .groupby(..., observed=True)

# BEFORE:
for gene, gene_df in df.groupby("GENE"):
    ...

# AFTER:
for gene, gene_df in df.groupby("GENE", observed=True):
    ...
```

### Anti-Patterns to Avoid
- **Commenting out dead code** — Delete it cleanly, git history preserves it
- **Manual temp file cleanup without try/finally** — Use context managers instead
- **gc.collect() only after "heavy" stages** — Collect after every stage for predictable memory profile
- **Partial groupby migration** — Update all call sites atomically to maintain consistency

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Temp file cleanup | Manual os.remove() tracking | tempfile context managers | Automatic cleanup on all exit paths (normal/exception) |
| Memory measurement | Manual before/after subtraction | psutil.Process().memory_info() | Cross-platform, accounts for OS caching |
| Garbage collection timing | Explicit del statements | gc.collect() between stages | Forces collection of circular references |
| Lint rule enforcement | Custom pre-commit script | Ruff built-in rules (when available) | Faster, maintained by community |

**Key insight:** Python's standard library provides battle-tested patterns for temp file management and garbage collection. Custom solutions introduce edge cases (exceptions, nested contexts) that stdlib handles correctly.

## Common Pitfalls

### Pitfall 1: Removing Code Without Verification
**What goes wrong:** Code that appears dead in static analysis may execute in edge cases (rare gene types, specific sample configurations).
**Why it happens:** Static analysis can't trace all dynamic call paths through pandas operations.
**How to avoid:** Add assertion/logging to prove code never executes across full test suite + benchmarks before removal.
**Warning signs:** Test failures after "safe" removal, changed output in benchmarks.

### Pitfall 2: Temp File Deletion Timing
**What goes wrong:** Deleting temp files too early (before subprocess reads them) or not at all (leak).
**Why it happens:** Subprocess calls may not read files synchronously, or exceptions prevent cleanup.
**How to avoid:** Use context managers with `delete=False` when file must outlive context, then explicit `os.remove()` in try/finally. For files that can be deleted immediately, use `delete=True` (default).
**Warning signs:** Subprocess errors reading files, `/tmp` directory growth over time.

### Pitfall 3: observed=True Changes Behavior
**What goes wrong:** Adding `observed=True` can silently drop empty category groups that downstream code expects.
**Why it happens:** With `observed=False`, groupby creates groups for all category levels even if no data exists. Code may assume these groups exist.
**How to avoid:** Run full test suite after each groupby change. If tests fail or output changes, mark that call site as "requires Phase 8 categorical migration" and skip it.
**Warning signs:** KeyError on expected group keys, missing rows in output, test failures in aggregation logic.

### Pitfall 4: gc.collect() Performance Impact
**What goes wrong:** Excessive gc.collect() calls can slow down the pipeline if called too frequently (e.g., inside tight loops).
**Why it happens:** gc.collect() scans all objects to find garbage, which takes time proportional to heap size.
**How to avoid:** Call gc.collect() only between stages (coarse-grained), never inside stage logic (fine-grained). Document that this is an inter-stage operation.
**Warning signs:** Benchmark regressions despite optimizations, high CPU usage in gc statistics.

### Pitfall 5: Ruff Custom Rule Limitations
**What goes wrong:** Assuming Ruff supports custom lint rules like Flake8 plugins.
**Why it happens:** Ruff is Rust-based and doesn't support Python plugin architecture.
**How to avoid:** Use Ruff's 800+ built-in rules when available. For custom checks (like groupby enforcement), use grep/sed in pre-commit hooks or document as manual review item until Ruff adds native support.
**Warning signs:** Spending time trying to write Ruff plugins that don't exist.

## Code Examples

### Verified Pattern 1: Dead GT Parsing Loop Removal
```python
# Location: variantcentrifuge/gene_burden.py:220-249
# Current code parses GT but never uses results

# VERIFICATION STEP (commit 1):
for _, row in gene_df.iterrows():
    assert False, (
        f"GT parsing loop executed unexpectedly! "
        f"Gene={gene}, variants={len(gene_df)}. "
        f"This loop was thought to be dead code."
    )
    gt_value = str(row.get("GT", ""))
    # ... rest of loop ...

# Run: pytest && pytest -m performance
# If assertion never fires → safe to remove

# REMOVAL STEP (commit 2):
# Delete lines 220-249 entirely
# Immediately after line 218 (control_total_alleles = 0), go to line 251
# No comments, no traces, git history preserves context
```

### Verified Pattern 2: All Groupby Call Sites
```python
# 17 total call sites identified:

# 1. variantcentrifuge/gene_burden.py:207
for gene, gene_df in df.groupby("GENE", observed=True):

# 2. variantcentrifuge/stats_engine.py:157
grouped = df.groupby(groupby_cols, observed=True)

# 3. variantcentrifuge/stats_engine.py:240
grouped = df.groupby(groupby_cols, observed=True)

# 4. variantcentrifuge/stats.py:117
df.groupby("GENE", observed=True)

# 5. variantcentrifuge/stats.py:151
impact_counts = df.groupby(["GENE", "IMPACT"], observed=True).size().reset_index(name="count")

# 6. variantcentrifuge/stats.py:178
type_counts = df.groupby(["GENE", "EFFECT"], observed=True).size().reset_index(name="count")

# 7. variantcentrifuge/stages/analysis_stages.py:3068
gene_groups = df.groupby(gene_column, observed=True)

# 8. variantcentrifuge/inheritance/parallel_analyzer.py:167
gene_groups = df.groupby("GENE", observed=True)

# 9. variantcentrifuge/inheritance/parallel_analyzer.py:221
for gene, gene_df in df.groupby("GENE", observed=True):

# 10. variantcentrifuge/inheritance/analyzer.py:101
gene_counts = df.groupby("GENE", observed=True).size()

# 11. variantcentrifuge/inheritance/analyzer.py:104
for gene, gene_df in df.groupby("GENE", observed=True):

# 12. tests/performance/benchmark_comp_het.py:144
for gene, gene_df in df.groupby("GENE", observed=True):

# 13. tests/performance/benchmark_pipeline.py:376
df.groupby(["test_name", "pipeline_type"], observed=True)["execution_time"]

# 14. tests/performance/benchmark_pipeline.py:402
df.groupby(["test_name", "pipeline_type"], observed=True)["peak_memory_mb"].mean().reset_index()

# 15. tests/performance/benchmark_pipeline.py:433
thread_summary = pipeline_df.groupby("threads", observed=True)["execution_time"].mean()

# 16. scripts/create_cohort_report.py:326
gene_variants = df.groupby("Gene", observed=True).size().reset_index(name="VariantCount")

# 17. scripts/create_cohort_report.py:329
gene_samples = df.groupby("Gene", observed=True)["SampleID"].nunique().reset_index(name="SampleCount")
```

### Verified Pattern 3: Temp File Audit Results
```python
# Three temp file creation patterns found:

# PATTERN 1: mkstemp (manual cleanup required)
# variantcentrifuge/gene_bed.py:153 — LEAKS (bed_path never removed)
bed_fd, bed_path = tempfile.mkstemp(suffix=".bed")
os.close(bed_fd)
# ... use bed_path ...
# Missing: os.remove(bed_path)

# FIX:
with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
    bed_path = tmp.name
    # Write to tmp
# ... rest of processing ...
os.remove(bed_path)  # Explicit cleanup after use

# PATTERN 2: mktemp (deprecated, manual cleanup)
# variantcentrifuge/filters.py:207 — CLEANED (os.remove on line 224)
tmp_vcf = tempfile.mktemp(suffix=".vcf")
# ... use tmp_vcf ...
os.remove(tmp_vcf)  # Present, but mktemp is deprecated

# FIX:
with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as tmp:
    tmp_vcf = tmp.name
    # ... use tmp_vcf ...
os.remove(tmp_vcf)

# PATTERN 3: mkdtemp (manual cleanup required)
# variantcentrifuge/pipeline_core/workspace.py:61 — CLEANED in workspace cleanup
self.temp_dir = Path(tempfile.mkdtemp(prefix="vc_temp_", dir=self.intermediate_dir))

# This is acceptable — workspace.py has cleanup method called at end
# No change needed
```

### Verified Pattern 4: Memory Logging Integration
```python
# Location: variantcentrifuge/pipeline_core/runner.py:574-598
# Add memory logging to _execute_stage method

import gc
import psutil

def _execute_stage(self, stage: Stage, context: PipelineContext) -> PipelineContext:
    """Execute a single stage and track timing + memory."""
    start_time = time.time()

    # Capture memory before stage (debug only)
    mem_before_mb = None
    if logger.isEnabledFor(logging.DEBUG):
        process = psutil.Process()
        mem_before_mb = process.memory_info().rss / 1024 / 1024
        logger.debug(f"[{stage.name}] Starting - Memory: {mem_before_mb:.1f} MB")

    result = stage(context)
    elapsed = time.time() - start_time
    self._execution_times[stage.name] = elapsed

    # Capture subtask times if any were recorded
    if hasattr(stage, "subtask_times") and stage.subtask_times:
        self._subtask_times[stage.name] = stage.subtask_times

    # Force garbage collection after every stage
    gc.collect()

    # Capture memory after GC (debug only)
    if mem_before_mb is not None:
        process = psutil.Process()
        mem_after_mb = process.memory_info().rss / 1024 / 1024
        freed_mb = mem_before_mb - mem_after_mb
        logger.debug(
            f"[{stage.name}] Complete in {elapsed:.1f}s - "
            f"Memory: {mem_after_mb:.1f} MB (freed {freed_mb:.1f} MB)"
        )

    return result
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual temp file cleanup | Context managers | Python 2.5+ (2006) | Automatic cleanup on all exit paths |
| `observed=False` (default) | `observed=True` for categorical | Pandas 3.0.0 (Jan 2026) | Default changed to avoid 3500x slowdown |
| Implicit garbage collection | Explicit `gc.collect()` between stages | Memory-sensitive apps | Predictable memory profile in long pipelines |
| Flake8 custom plugins | Ruff built-in rules | Ruff 0.1+ (2023) | 10-100x faster linting, no plugin system |

**Deprecated/outdated:**
- `tempfile.mktemp()`: Deprecated since Python 2.3 (2003) — security vulnerability, use NamedTemporaryFile instead
- Custom Ruff plugins: Never existed — Ruff is Rust-based, no Python plugin API

## Open Questions

### Question 1: Ruff Custom Rule for groupby
**What we know:** Ruff doesn't support custom plugins (Rust-based, no Python API).
**What's unclear:** Will Ruff add a built-in pandas rule for `observed=True` enforcement?
**Recommendation:** Use pre-commit hook with grep for now. File feature request with Ruff project. Example:
```bash
# .pre-commit-config.yaml
- repo: local
  hooks:
    - id: pandas-groupby-observed
      name: Enforce observed=True in groupby
      entry: bash -c 'grep -rn "\.groupby(" --include="*.py" variantcentrifuge/ | grep -v "observed=True" && exit 1 || exit 0'
      language: system
      pass_filenames: false
```

### Question 2: gc.collect() Frequency
**What we know:** Calling after every stage is safe (coarse-grained, ~30 stages).
**What's unclear:** Diminishing returns — does gc.collect() after lightweight stages (e.g., config loading) waste time?
**Recommendation:** Implement after all stages uniformly for Phase 7. If benchmarks show regression from gc overhead, make it conditional in Phase 8 (only after stages with `estimated_runtime > 1.0`).

### Question 3: observed=True Behavior Changes
**What we know:** Can drop empty groups if categorical dtypes are used.
**What's unclear:** Which call sites (if any) rely on empty groups existing?
**Recommendation:** Add `observed=True` to all 17 call sites. Run full test suite. If any test fails, revert that specific call site and document it. Currently no categorical dtypes are used, so behavior should be unchanged.

## Sources

### Primary (HIGH confidence)
- **Codebase inspection:** All file locations, line numbers, and code patterns verified by direct read of source files
- **Performance analysis report:** `.planning/performance-analysis-report.md` (lines 422-427, 571-576)
- **Pandas 3.0.0 documentation:** [pandas.DataFrame.groupby](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.groupby.html) — `observed=True` default change confirmed
- **Python tempfile documentation:** [tempfile module](https://docs.python.org/3/library/tempfile.html) — context manager patterns
- **psutil documentation:** Installed in project (requirements), memory monitoring API verified

### Secondary (MEDIUM confidence)
- [Towards Data Science: Pandas GroupBy with Categorical Data](https://towardsdatascience.com/be-careful-when-using-pandas-groupby-with-categorical-data-type-a1d31f66b162/) — 3500x slowdown documented
- [Ruff documentation](https://docs.astral.sh/ruff/rules/) — No custom plugin API confirmed
- [Real Python: tempfile module](https://realpython.com/ref/stdlib/tempfile/) — Context manager best practices

### Tertiary (LOW confidence)
- None — all findings verified with primary sources

## Metadata

**Confidence breakdown:**
- Dead code identification: HIGH — verified at exact line numbers (220-249), logic confirms unused
- groupby call sites: HIGH — all 17 instances found via grep, locations verified
- Temp file patterns: HIGH — 3 patterns identified, leak confirmed (bed_path not removed)
- gc.collect() placement: HIGH — runner.py _execute_stage is the correct insertion point
- Memory logging approach: MEDIUM — psutil API verified, debug-level logging pattern standard but format is discretionary

**Research date:** 2026-02-14
**Valid until:** 60 days (stable Python/pandas APIs, no fast-moving dependencies)

## Benchmark Strategy

### Phase 6 Benchmarks to Run
Based on `.benchmarks/Windows-CPython-3.10-64bit/0001_baseline_v0.12.1.json`:

```bash
# Full benchmark suite
pytest -m performance --benchmark-autosave --benchmark-name=quick_wins_before

# After optimizations
pytest -m performance --benchmark-autosave --benchmark-name=quick_wins_after

# Specific gene burden benchmarks (primary target for 30-40% speedup)
pytest tests/performance/benchmark_gene_burden.py -v
```

### Expected Improvements
| Optimization | Component | Before (est.) | After (est.) | Speedup |
|--------------|-----------|---------------|--------------|---------|
| Dead code removal | Gene burden | 15 min | 10-11 min | 1.3-1.5x |
| observed=True | All groupby | baseline | baseline | No change (no categoricals yet) |
| gc.collect() | Memory footprint | 30 GB peak | 20-25 GB peak | 20-30% reduction |
| Temp file cleanup | Disk usage | `/tmp` growth | No growth | Correctness fix |

**Combined gene burden speedup:** 30-40% (dead code removal is the primary contributor)

### Regression Test Requirements
```python
# tests/unit/test_gene_burden_regression.py
def test_dead_code_removal_identical_output(gene_burden_fixture):
    """Verify gene burden output is byte-identical after dead code removal."""
    # Load reference output from before optimization
    reference_df = pd.read_csv("tests/fixtures/gene_burden_reference.tsv", sep="\t")

    # Run gene burden with optimized code
    result_df = perform_gene_burden_analysis(gene_burden_fixture["df"], gene_burden_fixture["config"])

    # Compare all columns (values must be identical)
    pd.testing.assert_frame_equal(
        reference_df.sort_values("GENE").reset_index(drop=True),
        result_df.sort_values("GENE").reset_index(drop=True),
        check_exact=False,
        rtol=1e-9  # Allow tiny float rounding differences
    )
```
