---
phase: 37-association-resource-management-memory-streaming
plan: 02
subsystem: association
tags: [pandas, genotype, memory, gt-columns, gene-burden, association]

# Dependency graph
requires:
  - phase: 37-01
    provides: PipelineContext resource manager scaffolding for PERF-05 context
provides:
  - GT column lifecycle fix (PERF-05): per-sample GEN_N__GT columns preserved through all analysis stages
  - VariantAnalysisStage: key-column merge to re-attach GT cols after analyze_variants
  - GeneBurdenAnalysisStage: local-copy reconstruction, no context.current_dataframe writeback
  - AssociationAnalysisStage: direct GT column access, removed context.variants_df recovery fallback
  - Consolidated _find_gt_columns into single implementation (_find_per_sample_gt_columns from output_stages)
  - GT lifecycle test suite (10 tests)
affects:
  - 37-03 (context.variants_df memory reduction - fallback removal unlocks refactoring)
  - AssociationAnalysisStage tests
  - Any future stage that uses GT columns

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Local-copy reconstruction: reconstruct_gt_column(df.copy(), ...) — never reconstruct in-place and write back to context"
    - "Key-column merge for GT re-attachment: use CHROM/POS/REF/ALT/GENE as stable variant identity keys"
    - "Top-level module imports over repeated local function imports"

key-files:
  created:
    - tests/unit/test_gt_lifecycle.py
  modified:
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/gene_burden.py

key-decisions:
  - "VariantAnalysisStage re-attaches GT via key-column merge (not positional), so row reordering/filtering by analyze_variants is handled safely"
  - "_find_gt_columns in gene_burden.py is now an import alias for _find_per_sample_gt_columns from output_stages — single canonical implementation"
  - "AssociationAnalysisStage removed context.variants_df recovery fallback entirely (Fix 5 guarantees per-sample cols are always present)"
  - "reconstruct_gt_column/find_per_sample_gt_columns promoted to top-level imports in analysis_stages.py to avoid redefinition lint errors"

patterns-established:
  - "GT Lifecycle Invariant: per-sample GEN_N__GT columns must survive in context.current_dataframe through ALL analysis stages"
  - "Reconstruction-on-copy: always use df.copy() before reconstruct_gt_column, never overwrite df in context"

# Metrics
duration: 25min
completed: 2026-02-25
---

# Phase 37 Plan 02: GT Column Lifecycle Fix Summary

**GT drop/recover antipattern eliminated: per-sample GEN_N__GT columns now survive VariantAnalysis, GeneBurden, and Association stages via local-copy reconstruction and key-column merge re-attachment**

## Performance

- **Duration:** 25 min
- **Started:** 2026-02-25T09:52:24Z
- **Completed:** 2026-02-25T10:17:40Z
- **Tasks:** 2
- **Files modified:** 3 (analysis_stages.py, gene_burden.py, + 1 created)

## Accomplishments
- VariantAnalysisStage saves per-sample GT columns before reconstruction, re-attaches via CHROM/POS/REF/ALT/GENE merge after analyze_variants completes
- GeneBurdenAnalysisStage reconstructs GT on `df_for_burden` (local copy), never writes back to `context.current_dataframe`
- AssociationAnalysisStage accesses per-sample GT columns directly from `context.current_dataframe`; removed 30-line `context.variants_df` recovery fallback and removed `context.current_dataframe = df` at old line 2424
- `_find_gt_columns` in `gene_burden.py` replaced with import alias for `_find_per_sample_gt_columns` from `output_stages.py` — single canonical implementation
- Added 10-test GT lifecycle suite covering finder function, local-copy pattern, and consolidation

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix VariantAnalysisStage, GeneBurdenAnalysisStage, and AssociationAnalysisStage to preserve per-sample GT columns** - `82728af` (fix)
2. **Task 2: Consolidate duplicate _find_gt_columns and add GT lifecycle tests** - `8caa3c3` (refactor)

## Files Created/Modified
- `variantcentrifuge/stages/analysis_stages.py` - GT lifecycle cleanup across 3 stages; top-level import of _find_per_sample_gt_columns/reconstruct_gt_column
- `variantcentrifuge/gene_burden.py` - Removed duplicate _find_gt_columns function; added import alias from output_stages
- `tests/unit/test_gt_lifecycle.py` - 10 unit tests for GT column lifecycle invariant (created)

## Decisions Made
- Used key-column merge (CHROM/POS/REF/ALT/GENE) to re-attach GT cols in VariantAnalysisStage — safer than positional alignment because analyze_variants may filter or reorder rows
- `_find_per_sample_gt_columns` and `reconstruct_gt_column` promoted to top-level imports in analysis_stages.py to satisfy ruff's "no redefinition" lint rule
- AssociationAnalysisStage recovery fallback removed entirely rather than left as dead code — Fix 5 makes it unnecessary

## Deviations from Plan

None — plan executed exactly as written. One minor structural deviation: the plan showed `from ..stages.output_stages import` for local imports, but since analysis_stages.py is already inside `stages/`, the correct relative path is `from .output_stages import`. Updated all local imports and promoted to top-level to eliminate lint warnings.

## Issues Encountered
- Ruff lint error: local re-imports of `_find_per_sample_gt_columns` inside methods were flagged as redefinitions of the top-level import. Resolved by removing all local imports and using the top-level import throughout.

## Next Phase Readiness
- Per-sample GT columns now guaranteed in context.current_dataframe through all analysis stages
- `context.variants_df` recovery fallback removed from AssociationAnalysisStage — this unblocks 37-03 (context.variants_df memory reduction)
- 37-03 can now safely reduce/eliminate context.variants_df without breaking association GT column access

---
*Phase: 37-association-resource-management-memory-streaming*
*Completed: 2026-02-25*
