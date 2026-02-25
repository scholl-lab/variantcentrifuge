---
phase: 37-association-resource-management-memory-streaming
plan: 03
subsystem: association
tags: [memory-optimization, streaming, genotype-matrix, lazy-evaluation, dataclass, PERF-06]

# Dependency graph
requires:
  - phase: 37-02
    provides: GT column lifecycle fix ensuring per-sample GT columns survive through all stages
  - phase: 37-01
    provides: ResourceManager initialization fix for parallel association execution
provides:
  - _GenotypeMatrixBuilder picklable dataclass in analysis_stages.py (module-level)
  - Lazy per-gene matrix building invoked by engine just before tests
  - Matrix discard after each gene's tests complete in sequential path
  - Eager batch build + discard in parallel path to avoid double-pickling
  - 7 unit tests for builder pickling, empty/MAC-filtered genes, engine lifecycle
affects: [phase-35, phase-36, future-association-work]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "PERF-06: Lazy builder callback pattern — O(1 gene) peak memory for genotype matrix construction"
    - "Dataclass-as-callable pattern for picklable lazy computation"
    - "Pop-after-use pattern for matrix discard in engine gene loop"

key-files:
  created:
    - tests/unit/test_streaming_matrix.py
  modified:
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/association/engine.py

key-decisions:
  - "[37-03] _GenotypeMatrixBuilder is a module-level dataclass (not nested) — required for picklability with ProcessPoolExecutor"
  - "[37-03] Sequential path: build per-gene, discard after tests — O(1) peak memory"
  - "[37-03] Parallel path: build eagerly for all remaining genes before dispatch — avoids pickling both builder (with gene_df) and matrix to workers"
  - "[37-03] CADD/REVEL annotation extraction stays in the loop (cheap, needed early); site filter mask computed independently without pre-built matrix"
  - "[37-03] phenotype_vector and covariate_matrix retained in gene_data after discard — shared references, not per-gene copies"

patterns-established:
  - "Matrix builder stores all inputs as dataclass fields, __call__ builds and returns result dict"
  - "Engine checks for _genotype_matrix_builder key and genotype_matrix key before invoking builder"
  - "gene_data.pop() removes builder and matrix after use to enable GC"

# Metrics
duration: 74min
completed: 2026-02-25
---

# Phase 37 Plan 03: Streaming Genotype Matrix Construction Summary

**_GenotypeMatrixBuilder picklable dataclass replaces bulk pre-build loop: peak matrix memory drops from O(all_genes) to O(1 gene) in sequential association path**

## Performance

- **Duration:** 74 min
- **Started:** 2026-02-25T10:42:50Z
- **Completed:** 2026-02-25T11:57:00Z
- **Tasks:** 2
- **Files modified:** 3 (plus 1 created)

## Accomplishments

- Implemented `_GenotypeMatrixBuilder` dataclass at module level in `analysis_stages.py`, encapsulating all matrix build logic (build, sample-mask, MAC check) in a picklable callable
- Replaced bulk pre-build loop (built all 5K+ matrices before engine.run_all) with per-gene builder creation loop; engine invokes builder just before running each gene's tests
- Engine sequential path: build per-gene, discard matrix + mafs after tests → peak memory O(1 gene)
- Engine parallel path: build eagerly for remaining batch before dispatch (avoids pickling both builder + DataFrame to workers), discard after workers return
- 7 unit tests covering pickling, empty gene, MAC filter, result key completeness, matrix discard verification, and Fisher-only path (no builder)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create _GenotypeMatrixBuilder and replace pre-build loop** - `45c5702` (feat)
2. **Task 2: Engine consumes builder in sequential and parallel paths with matrix discard** - `0531cc3` (feat)

## Files Created/Modified

- `variantcentrifuge/stages/analysis_stages.py` - Added `_GenotypeMatrixBuilder` dataclass at module level (lines 44-113); replaced bulk matrix pre-build loop with lazy builder creation loop; added `numpy` and `dataclasses` imports
- `variantcentrifuge/association/engine.py` - Builder invocation in sequential gene loop (per-gene, with discard), builder invocation in parallel path (eager batch build before dispatch, discard after workers return), first-gene sequential run also invokes/discards builder
- `tests/unit/test_streaming_matrix.py` - 7 unit tests: `TestGenotypeMatrixBuilder` (5 tests) and `TestEngineBuilderConsumption` (2 tests)

## Decisions Made

- **Builder at module level, not nested**: Python's pickle requires classes to be importable by dotted path. A nested class inside `AssociationAnalysisStage._process()` would not be picklable, breaking the `ProcessPoolExecutor` parallel path.

- **Sequential vs parallel builder strategy differs**: Sequential path gets full O(1) benefit — build one gene, run tests, discard, then build next. Parallel path cannot do this because worker processes need the matrix already materialized when they receive the `gene_data` dict (can't pickle a builder-with-DataFrame alongside a matrix without double cost). Solution: build all remaining-batch matrices before dispatch, discard after workers return.

- **CADD/REVEL annotation extraction stays in the loop**: These annotations are cheap to extract and needed independently of the matrix. The site filter mask is now computed independently using the same vectorized logic as `build_genotype_matrix` (instead of reading it from `gene_data["variant_mafs"]` which no longer exists at loop time).

- **phenotype_vector and covariate_matrix not discarded**: After engine invokes builder and stores pv/cm in `gene_data`, these are kept after tests complete. They are shared references (same arrays for all genes), not per-gene allocations, so no memory benefit from discarding them.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] CADD/REVEL site filter mask computed without pre-built variant_mafs**

- **Found during:** Task 1 (Create _GenotypeMatrixBuilder and replace pre-build loop)
- **Issue:** Existing CADD/REVEL annotation alignment code used `len(gene_data["variant_mafs"])` to get `_n_kept` (post-filter variant count). With lazy builder, `gene_data["variant_mafs"]` doesn't exist at loop time.
- **Fix:** Removed the `_n_kept < _n_df` conditional shortcut; always compute the site filter mask directly from the DataFrame using the same vectorized logic. Added `_keep_mask_ann = None` when mask passes all variants (equivalent to no filtering).
- **Files modified:** `variantcentrifuge/stages/analysis_stages.py`
- **Verification:** Lint and all 2092 unit tests pass; logic is equivalent to before.
- **Committed in:** `45c5702` (Task 1 commit)

**2. [Rule 3 - Blocking] Fixed import of AssociationConfig in test file**

- **Found during:** Task 2 (creating tests/unit/test_streaming_matrix.py)
- **Issue:** Plan pseudocode used `from variantcentrifuge.association.config import AssociationConfig` but the module is `variantcentrifuge.association.base`.
- **Fix:** Changed import to `from variantcentrifuge.association.base import AssociationConfig`.
- **Files modified:** `tests/unit/test_streaming_matrix.py`
- **Verification:** All 7 new tests pass.
- **Committed in:** `0531cc3` (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (1 bug, 1 blocking)
**Impact on plan:** Both auto-fixes necessary for correctness. No scope creep.

## Issues Encountered

- Ruff lint/format violations required two fix rounds (line-too-long on engine.py condition, format-check on test file set literal formatting). Both resolved by reformatting.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 37 is now complete (all 3 plans done): ResourceManager fix (37-01), GT lifecycle fix (37-02), streaming matrix builder (37-03)
- The "genotype matrix is streamed per-gene" architecture invariant is now formally implemented — not just a stated constraint
- Peak memory for association analysis on 5K+ gene panels is O(1 gene matrix) in sequential mode (the default)
- Future phases 35 and 36 can assume the lazy builder pattern exists; if they need matrix access, they should follow the engine's pop-after-use pattern

---
*Phase: 37-association-resource-management-memory-streaming*
*Completed: 2026-02-25*
