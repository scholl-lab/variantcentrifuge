# Phase 37: Association Resource Management & Memory Streaming - Research

**Researched:** 2026-02-25
**Domain:** Python pipeline resource management, DataFrame lifecycle, memory streaming
**Confidence:** HIGH

## Summary

Phase 37 completes the association performance hardening by implementing three focused
structural improvements to the pipeline: a shared ResourceManager in PipelineContext (Fix 4),
elimination of the GT column drop/recover antipattern (Fix 5), and streaming genotype matrix
construction (Fix 6). Fixes 1-3 (vectorized GT parsing, groupby, auto-workers) were applied
in Phase 34.

The research is entirely internal — all code under modification is in this repository, all
APIs are well-understood, and no external library changes are needed. The primary research
questions are "how exactly does the current code work" and "what are the precise change
boundaries?" Both are answered by direct code inspection.

The standard approach for all three fixes is already articulated in the performance
investigation document. This research verifies the current state of the codebase to confirm
what is already done (Fix 3 via `_resolve_association_workers`) and what remains (Fix 4/5/6),
documents the precise change boundaries, and identifies edge cases.

**Primary recommendation:** The three fixes are independent and can be planned as three
separate tasks, each touching well-isolated code paths. All changes are refactoring with no
new external dependencies.

---

## Standard Stack

No new libraries needed. All changes use existing imports already present in the codebase.

### Core (already present)
| Component | Location | Purpose | Status |
|-----------|----------|---------|--------|
| `ResourceManager` | `variantcentrifuge/memory/resource_manager.py` | Memory/CPU detection and worker allocation | READY — well-designed, 378 lines, fully tested |
| `PipelineContext` | `variantcentrifuge/pipeline_core/context.py` | Single source of truth for pipeline state | NEEDS: `resource_manager` field added |
| `_find_per_sample_gt_columns` | `variantcentrifuge/stages/output_stages.py:42` | Detect per-sample GT columns | READY — already used by multiple stages |
| `reconstruct_gt_column` | `variantcentrifuge/stages/output_stages.py:61` | Pack per-sample GT → packed GT string | READY — move call site to output stages only |
| `build_genotype_matrix` | `variantcentrifuge/association/genotype_matrix.py:163` | Build per-gene genotype matrix | READY — already vectorized (Fix 1 done) |

### Supporting (already in use)
| Component | Location | Purpose |
|-----------|----------|---------|
| `numpy` | stdlib | Matrix operations in genotype matrix |
| `pandas groupby` | stdlib | Per-gene grouping (Fix 2 already applied) |
| `dataclasses.field` | stdlib | PipelineContext field additions |

**Installation:** No new packages needed.

---

## Architecture Patterns

### Current State (What Exists Before Phase 37)

**Fix 3 is already implemented** in `_resolve_association_workers()` at
`analysis_stages.py:2083`:

```python
# ALREADY DONE (Fix 3):
def _resolve_association_workers(cfg: dict, _get) -> int:
    raw = _get("association_workers", default=0, nullable=False)
    if raw == 0:
        threads = int(cfg.get("threads", 1))  # auto from --threads
        return threads
    if raw == -1:
        cpus = os.cpu_count() or 1  # all cores
        return cpus
    return int(raw)
```

This already correctly derives `association_workers` from `--threads` when not explicitly set.
The issue investigation's Fix 3 is complete. **Phase 37 only covers Fix 4/5/6.**

### Fix 4: Shared ResourceManager in PipelineContext

**Current anti-pattern:** Each stage that needs ResourceManager instantiates its own:

```python
# IN analysis_stages.py:928-930 (InheritanceAnalysisStage):
from ..memory import ResourceManager
rm = ResourceManager(config=context.config)

# IN analysis_stages.py:3019-3031 (ChunkedAnalysisStage):
from ..memory import ResourceManager
rm = ResourceManager(config=context.config)

# IN analysis_stages.py:3486-3491 (another location):
from ..memory import ResourceManager
rm = ResourceManager(config=context.config)

# IN analysis_stages.py:3758-3761 (yet another):
from ..memory import ResourceManager
rm = ResourceManager(config=context.config)
```

Each instantiation calls `psutil.virtual_memory()` and reads environment variables. With 4+
instantiations, this multiplies resource detection overhead and can produce inconsistent
results if environment changes between calls (unlikely but possible in containers).

**Target pattern:**

```python
# pipeline_core/context.py — add field:
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ..memory.resource_manager import ResourceManager

@dataclass
class PipelineContext:
    ...
    resource_manager: "ResourceManager | None" = None  # Add this field
```

```python
# pipeline.py — initialize once after context creation:
from .memory import ResourceManager
context = PipelineContext(args=args, config=initial_config, workspace=workspace)
context.resource_manager = ResourceManager(config=initial_config)
```

```python
# Each stage — replace local instantiation with context access:
rm = context.resource_manager or ResourceManager(config=context.config)  # fallback
```

**Pickling concern:** PipelineContext uses `PicklableLock` to handle pickling for parallel
execution. ResourceManager is a pure computation object with no locks — it is safely picklable
as long as it holds no file handles. The `_detect_memory()` method only reads env vars and
files during `__init__`, so a pickled+unpickled ResourceManager retains correct values.
This is SAFE.

**merge_from concern:** `PipelineContext.merge_from()` at line 258 merges contexts from
parallel stage execution. The ResourceManager should NOT be merged from parallel contexts
(we want the single shared instance). The merge_from method should simply not touch
`resource_manager` — this is the correct behavior since the field will be initialized on
the parent context and parallel contexts should inherit it.

### Fix 5: Eliminate GT Column Drop/Recover Cycle

**Current flow (the antipattern):**

```
DataFrameLoadingStage
  → context.current_dataframe has 5,125 per-sample GEN_N__GT columns
  → context.variants_df = df (copy stored for fallback — 16 GB duplicate)

VariantAnalysisStage (analysis_stages.py:1464-1471):
  → if "GT" not in df.columns:
      gt_cols = _find_per_sample_gt_columns(df)
      df = reconstruct_gt_column(df.copy(), context.vcf_samples)
      # DROPS per-sample GEN_N__GT columns, replaces with packed "Sample1(0/1);..."

GeneBurdenAnalysisStage (analysis_stages.py:1827-1835):
  → if "GT" not in df.columns:
      df = reconstruct_gt_column(df.copy(), context.vcf_samples)
      context.current_dataframe = df
      # SAME — but at this point GT may already exist from VariantAnalysis

AssociationAnalysisStage (analysis_stages.py:2404-2446):
  → needs_regression = "coast" in tests or "skat" in tests or ...
  → if "GT" not in df.columns:
      if needs_regression:
          df_with_per_sample_gt = df   # Save BEFORE dropping
      df = reconstruct_gt_column(df.copy(), ...)
  → elif "GT" in df.columns and needs_regression:
      # RECOVERY PATH: try context.variants_df fallback (16 GB duplicate)
      fallback_df = context.variants_df   # ← THE 16 GB DUPLICATE
      gt_cols_fb = _find_per_sample_gt_columns(fallback_df)
      df_with_per_sample_gt = fallback_df  # Recover per-sample cols
```

**Problem:** `context.variants_df` exists solely as a safety net to recover per-sample GT
columns that VariantAnalysisStage already destroyed. This is a 16 GB memory waste.

**Target flow:**

```
DataFrameLoadingStage
  → context.current_dataframe has per-sample GEN_N__GT columns (unchanged)
  → context.variants_df = None (no longer needed as fallback)

VariantAnalysisStage:
  → DO reconstruct GT for analyze_variants (it needs packed GT for TSV writing)
  → BUT: reconstruct on a separate copy just for the analysis function call
  → DO NOT store reconstructed df back to context.current_dataframe
  → context.current_dataframe KEEPS per-sample columns

GeneBurdenAnalysisStage:
  → Works with context.current_dataframe (has per-sample GT columns)
  → Calls reconstruct_gt_column only on local temp copy

AssociationAnalysisStage:
  → context.current_dataframe still has per-sample GT columns
  → No recovery needed — df_with_per_sample_gt = df directly
  → No fallback to context.variants_df needed

TSVOutputStage (output_stages.py:574-578):
  → reconstruct_gt_column called here (already does this!)
  → Per-sample columns are dropped AT OUTPUT TIME only

ExcelReportStage (output_stages.py:717-720):
  → reconstruct_gt_column called here (already does this!)
```

The key insight: `TSVOutputStage` and `ExcelReportStage` already call `reconstruct_gt_column`.
The fix is to STOP calling it in analysis stages and trust the output stages to handle it.

**Critical detail for VariantAnalysisStage:** The stage writes a temp TSV for `analyze_variants`
at line 1490. The `analyze_variants` function expects a packed GT column. The fix here is to:
1. Build the temp TSV with a reconstructed GT column (local operation, don't save to context)
2. After getting back analysis results, merge them onto the ORIGINAL per-sample df (not the
   reconstructed one)
3. context.current_dataframe retains per-sample columns

**Impact on context.variants_df:** After Fix 5, `context.variants_df` is no longer needed
as a fallback. It can be removed from the recovery paths. However, `DataFrameLoadingStage`
sets it for "in-memory pass-through" (line 588). This is a different use case (checkpoint
skip). The field should be kept but the recovery fallback in AssociationStage removed.

### Fix 6: Stream Genotype Matrix Construction

**Current flow:**

```python
# analysis_stages.py:2593-2648 — ALL matrices built before engine.run_all()
for gene_data in gene_burden_data:  # 5,084 iterations
    gene_df = gene_groups.get_group(gene_name)
    geno, mafs, sample_mask, gt_warnings = build_genotype_matrix(...)  # Each ~50MB
    gene_data["genotype_matrix"] = geno  # STORED IN DICT — never freed
    gene_data["variant_mafs"] = mafs

# After all 5K matrices built (16 GB peak):
result_df = engine.run_all(gene_burden_data)  # All matrices still in memory
```

**Problem:** All 5,084 matrices (each ~50MB for 5K samples x ~10 variants) are built and
stored before the engine processes any. Peak memory: 5,084 x 50MB ≈ 254 GB (actual observed:
16 GB due to smaller average variant count per gene, but still excessive).

**Target flow — lazy callback pattern:**

```python
# Build a lazy callback instead of pre-building matrices:
def _make_matrix_builder(gene_df, vcf_samples_list, gt_columns_for_matrix,
                          is_binary, assoc_config, phenotype_vector):
    """Returns a callable that builds the genotype matrix on demand."""
    def build():
        return build_genotype_matrix(
            gene_df,
            vcf_samples_list,
            gt_columns_for_matrix,
            is_binary=is_binary,
            missing_site_threshold=assoc_config.missing_site_threshold,
            missing_sample_threshold=assoc_config.missing_sample_threshold,
            phenotype_vector=phenotype_vector,
        )
    return build

# Store callback, not the matrix:
gene_data["_genotype_matrix_builder"] = _make_matrix_builder(...)
# No "genotype_matrix" key until engine calls the builder
```

Then in `engine.py`, the engine calls the builder just before running tests:

```python
# engine.py — in the per-gene loop:
if "_genotype_matrix_builder" in gene_data and "genotype_matrix" not in gene_data:
    builder = gene_data.pop("_genotype_matrix_builder")
    geno, mafs, sample_mask, warnings = builder()
    # Apply sample mask, MAC check...
    gene_data["genotype_matrix"] = geno
    gene_data["variant_mafs"] = mafs
    # After test.run(), the matrix can be discarded:
    gene_data.pop("genotype_matrix", None)
    gene_data.pop("variant_mafs", None)
```

**Alternative approach:** Move the matrix building logic entirely into the engine's gene
loop by passing `gt_source_df` and `gene_groups` to the engine. This is cleaner but
requires engine API changes.

**Recommended approach:** The lazy callback pattern is cleanest because:
1. It keeps AssociationAnalysisStage responsible for the matrix building logic
2. The engine remains test-focused (receives data, not data sources)
3. MAC check and sample masking stay in AssociationAnalysisStage (where they belong)
4. The callback captures all needed references (gene_df, configs) in its closure

**Memory impact:** Peak memory drops from (n_genes x matrix_size) to (1-2 x matrix_size)
= ~50-100MB instead of 16+ GB.

**Parallel execution note:** The engine's parallel path uses `ProcessPoolExecutor` and
pickles `gene_data` dicts. Lambdas cannot be pickled. The callback must be a regular
(named) function or a picklable callable class. Using `functools.partial` or a named
closure function works. A dataclass with `__call__` is safest for pickling.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Memory detection | Custom psutil wrapper | `ResourceManager` (already exists) | Has 5-priority chain: CLI > SLURM > PBS > cgroup > psutil |
| Worker count | Manual min/max math | `ResourceManager.auto_workers()` | Already accounts for memory + CPU constraints |
| GT column detection | Custom regex | `_find_per_sample_gt_columns()` (output_stages.py:42) | Handles both `GEN_0__GT` and `GEN[0].GT` formats |
| Per-gene DataFrame filter | `df[df.GENE == gene]` loop | `groupby("GENE").get_group()` | Fix 2 already applied — O(1) per lookup |

---

## Common Pitfalls

### Pitfall 1: Pickling Lambdas in Parallel Engine

**What goes wrong:** If the lazy callback is implemented as a lambda, pickling for the
ProcessPoolExecutor will fail with `AttributeError: Can't pickle local object`.

**Why it happens:** Python's pickle module cannot serialize lambdas or locally-defined
functions by name.

**How to avoid:** Use `functools.partial` with a module-level function, or a named class
with `__call__`. Example:

```python
# WRONG — lambda not picklable:
gene_data["_builder"] = lambda: build_genotype_matrix(gene_df, ...)

# RIGHT — functools.partial is picklable:
from functools import partial
gene_data["_builder"] = partial(build_genotype_matrix, gene_df, vcf_samples_list, ...)

# ALSO RIGHT — named class with __call__:
@dataclass
class GenotypeMatrixBuilder:
    gene_df: pd.DataFrame
    vcf_samples: list[str]
    # ... other fields
    def __call__(self) -> tuple[np.ndarray, np.ndarray, list[bool], list[str]]:
        return build_genotype_matrix(self.gene_df, self.vcf_samples, ...)
```

**Warning signs:** `AttributeError: Can't pickle local object` during parallel association.

### Pitfall 2: merge_from Overwriting ResourceManager

**What goes wrong:** `PipelineContext.merge_from()` updates fields from parallel contexts.
If `resource_manager` is added without excluding it from merge, parallel contexts (which
also have the shared instance) will overwrite the parent's instance unnecessarily.

**Why it happens:** The current `merge_from()` method has explicit field-by-field merging.
Adding `resource_manager` to the dataclass without adding a merge exclusion means it would
be silently left unmerged (no entry in `merge_from`), which is actually correct behavior.

**How to avoid:** Simply do not add a `resource_manager` merge case in `merge_from()`.
The field will remain at whatever the parent context set it to, which is correct.

**Warning signs:** Inconsistent memory/CPU detection between stages.

### Pitfall 3: VariantAnalysisStage Temp TSV Needs Packed GT

**What goes wrong:** `analyze_variants()` is called with a temp TSV. The function expects
a packed `GT` column (e.g., `Sample1(0/1);Sample2(1/1)`). If per-sample columns are kept
in `context.current_dataframe` and not reconstructed for the temp TSV, `analyze_variants`
will fail.

**Why it happens:** The `analyze_variants` function is legacy code that expects the packed
GT format. Fix 5 preserves per-sample columns in the context DataFrame, but the temp TSV
written to disk still needs packed GT.

**How to avoid:** In VariantAnalysisStage, reconstruct GT on a LOCAL COPY for temp TSV
writing only. Do NOT assign the reconstructed DataFrame back to `context.current_dataframe`.

```python
# CORRECT pattern:
temp_df = reconstruct_gt_column(df.copy(), context.vcf_samples)  # local only
temp_df.to_csv(temp_tsv, ...)

# After analyze_variants, merge results onto original df (with per-sample cols):
# context.current_dataframe stays as df (with per-sample columns)
```

**Warning signs:** `analyze_variants` returning no results or missing GT-dependent columns
in the final output.

### Pitfall 4: GeneBurdenAnalysisStage GT Requirement

**What goes wrong:** `GeneBurdenAnalysisStage` at line 1837 checks `if "GT" not in
df.columns` and returns early with an error. After Fix 5, per-sample columns are preserved
but there is no packed `GT` column.

**Why it happens:** Gene burden uses the `GT` column via `assign_case_control_counts` which
iterates GT strings. This function expects packed format.

**How to avoid:** GeneBurdenAnalysisStage should reconstruct GT for its own local use
(same pattern as VariantAnalysisStage after Fix 5). The reconstruction should be on a
local copy, not stored back to context.

### Pitfall 5: ResourceManager Memory Init Cost

**What goes wrong:** Each ResourceManager `__init__` calls `psutil.virtual_memory()` which
is a syscall. This is fast (~1ms) but if ResourceManager is instantiated N times in tests,
it adds up.

**Why it happens:** Currently 4+ instantiations in analysis_stages.py.

**How to avoid:** Fix 4 addresses this — one shared instance. For tests, mock
`ResourceManager` at the class level or use the config-based override
(`{"max_memory_gb": 8.0, "SLURM_CPUS_PER_TASK": None}`).

### Pitfall 6: context.variants_df Still Used by Other Code

**What goes wrong:** After Fix 5, removing `context.variants_df` usage from AssociationStage
is safe, but `context.variants_df` may still be used by DataFrameLoadingStage (set) and
potentially other code.

**Why it happens:** `variants_df` serves two purposes: (1) in-memory pass-through for
checkpoint skip (DataFrameLoadingStage), (2) fallback for per-sample GT recovery
(AssociationStage). Fix 5 eliminates purpose (2) but purpose (1) must be preserved.

**How to avoid:** Only remove the RECOVERY path (`fallback_df = context.variants_df` in
AssociationStage lines 2430-2437). Do not remove the field from PipelineContext or the
setting in DataFrameLoadingStage.

---

## Code Examples

### Shared ResourceManager in PipelineContext

```python
# Source: variantcentrifuge/pipeline_core/context.py — current state (no resource_manager)
# Target addition:

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ..memory.resource_manager import ResourceManager  # avoids circular import

@dataclass
class PipelineContext:
    # ... existing fields ...
    resource_manager: "ResourceManager | None" = None  # new field
```

```python
# Source: variantcentrifuge/pipeline.py — initialization after context creation
# After line 527 (context = PipelineContext(...)):
from .memory import ResourceManager
context.resource_manager = ResourceManager(config=initial_config)
```

```python
# Source: variantcentrifuge/stages/analysis_stages.py — stage usage pattern
# Replace all local ResourceManager instantiations with:
rm = context.resource_manager
if rm is None:
    from ..memory import ResourceManager
    rm = ResourceManager(config=context.config)  # fallback for tests/edge cases
```

### GT Column Preservation Pattern (Fix 5)

```python
# Source: variantcentrifuge/stages/analysis_stages.py — VariantAnalysisStage
# BEFORE (drops per-sample cols from context):
if "GT" not in df.columns and context.vcf_samples:
    gt_cols = _find_per_sample_gt_columns(df)
    if gt_cols:
        df = reconstruct_gt_column(df.copy(), context.vcf_samples)
        # df now has packed GT, per-sample cols dropped
# writes df to temp_tsv, runs analyze_variants

# AFTER (reconstructs locally for temp TSV only):
df_for_analysis = df
if "GT" not in df.columns and context.vcf_samples:
    gt_cols = _find_per_sample_gt_columns(df)
    if gt_cols:
        df_for_analysis = reconstruct_gt_column(df.copy(), context.vcf_samples)
        # df_for_analysis used for temp TSV; df (with per-sample cols) unchanged
# writes df_for_analysis to temp_tsv, runs analyze_variants
# context.current_dataframe = df  (per-sample cols preserved)
```

```python
# Source: variantcentrifuge/stages/analysis_stages.py — AssociationAnalysisStage
# BEFORE (complex recovery from variants_df):
df_with_per_sample_gt: pd.DataFrame | None = None
if "GT" not in df.columns and context.vcf_samples:
    gt_cols = _find_per_sample_gt_columns(df)
    if gt_cols:
        if needs_regression:
            df_with_per_sample_gt = df  # save before drop
        df = reconstruct_gt_column(df.copy(), context.vcf_samples)
elif "GT" in df.columns and needs_regression:
    fallback_df = context.variants_df  # 16 GB duplicate!
    ...

# AFTER (per-sample cols always present, no recovery needed):
# context.current_dataframe has per-sample GT columns (Fix 5 preserved them)
gt_cols = _find_per_sample_gt_columns(df)
if needs_regression and gt_cols:
    df_with_per_sample_gt = df  # direct reference, no copy needed
if "GT" not in df.columns:
    df = reconstruct_gt_column(df.copy(), context.vcf_samples)
```

### Lazy Genotype Matrix Builder (Fix 6)

```python
# Source: variantcentrifuge/stages/analysis_stages.py — picklable callable
from dataclasses import dataclass
import pandas as pd
import numpy as np

@dataclass
class _GenotypeMatrixBuilder:
    """Picklable lazy builder for genotype matrices (Fix 6)."""
    gene_df: pd.DataFrame
    vcf_samples: list[str]
    gt_columns: list[str]
    is_binary: bool
    missing_site_threshold: float
    missing_sample_threshold: float
    phenotype_vector: "np.ndarray | None"

    def __call__(self):
        from ..association.genotype_matrix import build_genotype_matrix
        return build_genotype_matrix(
            self.gene_df,
            self.vcf_samples,
            self.gt_columns,
            is_binary=self.is_binary,
            missing_site_threshold=self.missing_site_threshold,
            missing_sample_threshold=self.missing_sample_threshold,
            phenotype_vector=self.phenotype_vector,
        )
```

```python
# Source: variantcentrifuge/association/engine.py — consume in gene loop
# In run_all(), before test.run():
if "_genotype_matrix_builder" in gene_data:
    builder = gene_data.pop("_genotype_matrix_builder")
    geno, mafs, sample_mask, gt_warnings = builder()
    # Apply sample mask and MAC check (moved from analysis_stages.py):
    ...
    gene_data["genotype_matrix"] = geno
    gene_data["variant_mafs"] = mafs

# After all tests for this gene:
gene_data.pop("genotype_matrix", None)
gene_data.pop("variant_mafs", None)
```

**Note on MAC check and sample mask:** Currently in AssociationAnalysisStage
(lines 2626-2648). If the builder is called inside the engine, these checks should
move into the engine's gene loop OR remain in the builder's `__call__`. The cleanest
split: builder builds the raw matrix, engine applies MAC check. But this requires
engine to know about sample masking, which is business logic. Alternatively, the builder
`__call__` returns the fully processed matrix (post-filter, post-mask) and also stores
the processed `phenotype_vector` and `covariate_matrix` in the dict.

---

## Current State of Fixes (Verified)

| Fix | Status | Evidence |
|-----|--------|---------|
| Fix 1: Vectorize GT parsing | DONE (Phase 34) | `genotype_matrix.py:37-98` has `_vectorized_parse_gt_columns()` |
| Fix 2: GroupBy for per-gene filter | DONE (Phase 34) | `analysis_stages.py:2600` has `gt_source_df.groupby("GENE")` |
| Fix 3: auto-workers from --threads | DONE (Phase 34) | `analysis_stages.py:2083-2101` has `_resolve_association_workers()` |
| Fix 4: Shared ResourceManager | NOT DONE | `PipelineContext` has no `resource_manager` field; stages create local instances |
| Fix 5: GT lifecycle cleanup | NOT DONE | `analysis_stages.py:2430` still has `context.variants_df` recovery fallback |
| Fix 6: Stream matrix construction | NOT DONE | `analysis_stages.py:2601-2648` builds all matrices before engine call |

### Verified File Locations for Each Fix

**Fix 4 changes:**
- `variantcentrifuge/pipeline_core/context.py` — add `resource_manager` field (dataclass)
- `variantcentrifuge/pipeline.py` — initialize `ResourceManager` after line 527
- `variantcentrifuge/stages/analysis_stages.py` — replace 4 local instantiations:
  - Line 930 (InheritanceAnalysisStage)
  - Line 3031 (ChunkedAnalysisStage `_calculate_chunk_size`)
  - Lines 3489-3491 (another `_calculate_chunk_size` variant)
  - Lines 3761 (ChunkedAnalysisStage `_run_chunked_analysis`)

**Fix 5 changes:**
- `variantcentrifuge/stages/analysis_stages.py`:
  - `VariantAnalysisStage._process()` lines 1464-1471 — reconstruct only for temp TSV
  - `GeneBurdenAnalysisStage._process()` lines 1827-1835 — reconstruct only for local use
  - `AssociationAnalysisStage._process()` lines 2404-2446 — remove recovery fallback
- No changes needed to `output_stages.py` (TSVOutputStage and ExcelReportStage already
  call `reconstruct_gt_column` at output time — lines 574-578 and 717-720)

**Fix 6 changes:**
- `variantcentrifuge/stages/analysis_stages.py` lines 2593-2648 — replace matrix building
  with builder creation
- `variantcentrifuge/association/engine.py` — consume builder in gene loop(s):
  - Sequential path: line 417 (`for gene_data in sorted_data:`)
  - Parallel path: `_run_gene_worker()` function (line 86)

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| iterrows() GT parsing | Vectorized column-wise | Phase 34 | 50-100x speedup |
| Per-gene linear scan | groupby().get_group() | Phase 34 | ~3000x speedup |
| Hardcoded association_workers=1 | Auto from --threads | Phase 34 | N-core parallelism |
| All matrices upfront (16 GB) | Still upfront (Phase 37 will fix) | Phase 37 target | 16x memory reduction |
| Multiple ResourceManager instances | Still multiple (Phase 37 will fix) | Phase 37 target | Consistency + minor perf |
| GT drop/recover cycle | Still present (Phase 37 will fix) | Phase 37 target | Eliminate 16 GB duplicate |

---

## Open Questions

1. **MAC check and sample mask location after Fix 6**
   - What we know: MAC check (lines 2636-2643) and sample mask application (2626-2634)
     currently live in AssociationAnalysisStage AFTER matrix building
   - What's unclear: Should these move into the engine's gene loop or into the builder?
   - Recommendation: Move into the builder's `__call__` method. The builder becomes a
     "processed matrix provider" — it builds, filters, masks, and returns the final matrix
     along with the masked phenotype/covariate vectors. This keeps engine test-focused.

2. **context.variants_df retention after Fix 5**
   - What we know: `variants_df` is set in DataFrameLoadingStage for "in-memory pass-through"
     (checkpoint skip functionality), and also used as recovery fallback in AssociationStage
   - What's unclear: After removing the recovery fallback, does `variants_df` still serve
     any purpose?
   - Recommendation: Keep the field and the DataFrameLoadingStage assignment (checkpoint
     feature). Only remove the recovery USAGE in AssociationStage lines 2430-2437.

3. **Engine parallel path and per-gene matrix discard**
   - What we know: The parallel path in `engine.run_all()` pickles `gene_data` dicts and
     sends them to worker processes via `ProcessPoolExecutor`
   - What's unclear: If `gene_data["_genotype_matrix_builder"]` is called inside the worker,
     the `gene_df` captured by the builder closure will be pickled. For 5K samples x all-gene
     variants, this may be large.
   - Recommendation: In the parallel path, the builder should be called BEFORE pickling
     (i.e., in AssociationAnalysisStage before `engine.run_all()`), or the engine's parallel
     path should build matrices before `ProcessPoolExecutor.map()`. The sequential path
     (default for most tests) is straightforward.

---

## Sources

### Primary (HIGH confidence)
All findings are from direct code inspection of this repository.

- `variantcentrifuge/memory/resource_manager.py` — ResourceManager API, 378 lines, all methods
- `variantcentrifuge/pipeline_core/context.py` — PipelineContext dataclass, all fields
- `variantcentrifuge/pipeline.py` — context creation at line 527
- `variantcentrifuge/stages/analysis_stages.py` — all 4 ResourceManager instantiations,
  GT lifecycle (lines 1464-1471, 1827-1835, 2404-2446), matrix building loop (2593-2648)
- `variantcentrifuge/stages/output_stages.py` — `reconstruct_gt_column` (line 61),
  `_find_per_sample_gt_columns` (line 42), TSVOutputStage GT handling (574-578),
  ExcelReportStage GT handling (717-720)
- `variantcentrifuge/association/engine.py` — sequential loop (417), parallel path (370-414)
- `variantcentrifuge/association/genotype_matrix.py` — `build_genotype_matrix` (163),
  `_vectorized_parse_gt_columns` (37)
- `.planning/association-performance-investigation.md` — Fix 4/5/6 specification

### Secondary (MEDIUM confidence)
- Python `pickle` documentation on lambda/closure pickling limitations — well-known constraint
- `dataclasses` module behavior for `field(default=None)` — standard Python

---

## Metadata

**Confidence breakdown:**
- Fix 4 (shared ResourceManager): HIGH — code is trivial, risks are minimal and identified
- Fix 5 (GT lifecycle): HIGH — all call sites verified by inspection; output stages already
  handle reconstruction at output time
- Fix 6 (streaming matrices): HIGH for approach; MEDIUM for parallel path interaction
  (pickling of large DataFrames in builder closures needs care)

**Research date:** 2026-02-25
**Valid until:** 2026-03-25 (stable internal codebase, no external dependencies)
