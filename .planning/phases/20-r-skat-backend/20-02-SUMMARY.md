---
phase: 20-r-skat-backend
plan: 02
subsystem: association
tags: [rpy2, skat, skat-binary, skat-o, r, gene-burden, association, thread-safety, memory-management]

# Dependency graph
requires:
  - phase: 20-01
    provides: RSKATBackend skeleton, NullModelResult dataclass, SKATBackend ABC, RSKATTest stub
  - phase: 19
    provides: genotype matrix builder, AssociationConfig with trait_type/skat_method/variant_weights

provides:
  - RSKATBackend.fit_null_model(): SKAT_Null_Model with Adjustment=TRUE, global-env R object pattern
  - RSKATBackend.test_gene(): SKATBinary (binary) or SKAT (continuous), SKAT-O rho extraction
  - RSKATBackend memory management: _run_r_gc(), _check_r_heap(), cleanup() with timing summary
  - RSKATTest.run(): lazy null model, test_gene per gene, GC every 100 genes, progress logging
  - RSKATTest.prepare()/finalize(): lifecycle hooks for panel warning, timing, backend cleanup
  - AssociationTest ABC: default no-op prepare() and finalize() methods
  - engine.run_all(): calls prepare() before and finalize() after gene loop
  - AssociationAnalysisStage.parallel_safe = False (SKAT-08 requirement)
affects: ["20-03", "21-pure-python-skat", "22-acat-o-diagnostics"]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "R global-env pattern: assign to ._vc_* globals, call SKAT, rm(list=ls(pattern='\\._vc_'))"
    - "withCallingHandlers in R code string for warning capture without Python overhead"
    - "Lazy null model caching: fit once on first gene, reuse for entire cohort"
    - "Lifecycle hooks: prepare(n_genes)/finalize() called by engine around gene loop"
    - "GC_INTERVAL=100 constant for periodic R gc() + heap monitoring"
    - "NA_Real guard: from rpy2.rinterface import NA_Real; p = None if p_val_r is NA_Real"

key-files:
  created: []
  modified:
    - variantcentrifuge/association/backends/r_backend.py
    - variantcentrifuge/association/tests/skat_r.py
    - variantcentrifuge/association/engine.py
    - variantcentrifuge/association/base.py
    - variantcentrifuge/stages/analysis_stages.py

key-decisions:
  - "withCallingHandlers embedded in R code string (not rpy2 Python-side callbacks) — cleaner for string-based R execution"
  - "GC triggered in RSKATTest.run() every GC_INTERVAL genes (not only in finalize) to prevent R heap accumulation mid-run"
  - "_genes_processed tracked independently in both RSKATBackend and RSKATTest — backend tracks for cleanup() summary"
  - "prepare()/finalize() as ABC no-ops in AssociationTest — Fisher/burden tests are unaffected; only RSKATTest overrides"
  - "parallel_safe=False on AssociationAnalysisStage is unconditional (not gated on skat_backend config key)"

patterns-established:
  - "R backend always cleans up ._vc_* globals after each call (fit_null_model and test_gene)"
  - "SKAT-O rho extracted only when method='SKATO'; None otherwise (no key error)"
  - "Progress logging: max(10, min(50, n//10)) interval — 10 genes log every gene, 500 genes log every 50"
  - "Large panel warning at >2000 genes emitted in prepare() before any R calls"

# Metrics
duration: 5min
completed: 2026-02-20
---

# Phase 20 Plan 02: RSKATBackend Full Implementation Summary

**SKAT_Null_Model(Adjustment=TRUE) + SKATBinary/SKAT dispatch + SKAT-O rho + R heap monitoring via rpy2 global-env pattern**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-20T15:27:34Z
- **Completed:** 2026-02-20T15:32:50Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- RSKATBackend fully operational: fit_null_model uses SKAT_Null_Model(Adjustment=TRUE) with global-env R object pattern;
  test_gene dispatches SKATBinary (binary) or SKAT (continuous); SKAT-O extracts optimal rho; statistical NA mapped to
  p_value=None via NA_Real guard; R warnings captured via withCallingHandlers in R code string
- Memory management complete: _run_r_gc() calls R gc() + Python gc.collect() every 100 genes;
  _check_r_heap() monitors R Vcells heap and emits WARNING if > 4 GB; cleanup() logs timing summary
- RSKATTest lifecycle wired to engine: prepare(n_genes) emits large panel warning + start message;
  finalize() logs aggregate timing + calls backend.cleanup(); engine.run_all() calls both hooks;
  AssociationTest ABC gets default no-op prepare()/finalize() so Fisher/burden tests are unaffected
- AssociationAnalysisStage.parallel_safe = False explicitly (SKAT-08 requirement — rpy2 is not thread-safe)

## Task Commits

Each task was committed atomically:

1. **Task 1: RSKATBackend fit_null_model + test_gene + memory management** - `b5c6f7b` (feat)
2. **Task 2: RSKATTest gene loop + stage parallel_safe + lifecycle hooks** - `d79a7b7` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified
- `variantcentrifuge/association/backends/r_backend.py` - Full RSKATBackend: fit_null_model, test_gene,
  _run_r_gc, _check_r_heap, cleanup; GC_INTERVAL=100, R_HEAP_WARNING_GB=4.0, LARGE_PANEL_THRESHOLD=2000
- `variantcentrifuge/association/tests/skat_r.py` - RSKATTest.run() with GC/progress; prepare()/finalize() hooks;
  _parse_weights_beta updated with "uniform" -> (1.0, 1.0)
- `variantcentrifuge/association/engine.py` - run_all() calls test.prepare(n) before loop, test.finalize() after
- `variantcentrifuge/association/base.py` - AssociationTest ABC gets default no-op prepare() and finalize()
- `variantcentrifuge/stages/analysis_stages.py` - AssociationAnalysisStage.parallel_safe property = False

## Decisions Made

- **withCallingHandlers in R string**: R warning capture is implemented via an R code string passed to ro.r() that
  includes `withCallingHandlers`. This is cleaner than rpy2 Python-side warning handlers and matches the research
  pattern from 20-RESEARCH.md.
- **Dual gene counter**: RSKATBackend and RSKATTest each track `_genes_processed` independently. The backend's counter
  feeds cleanup() timing summary; RSKATTest's counter feeds progress logging and GC triggering.
- **GC in run() not finalize()**: R GC is triggered every 100 genes inside RSKATTest.run() (not just at finalize) to
  prevent R heap accumulation during long runs. This is the correct place since test_gene() is the source of R objects.
- **parallel_safe=False unconditional**: The property always returns False regardless of which tests are active.
  Even if only Fisher is running alongside SKAT, the stage is correctly marked non-parallel-safe.

## Deviations from Plan

None — plan executed exactly as written. The run() method in skat_r.py was refactored for clarity (separated
existing stub code) but the logic matches the plan specification. `_parse_weights_beta` received the "uniform"
case as specified in the plan.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- RSKATBackend is functionally complete and ready for Plan 20-03 (mocked unit tests)
- The R warning capture (withCallingHandlers) and NA_Real guard will be tested with mocked rpy2 in Plan 20-03
- The null model fit-once invariant, SKATBinary vs SKAT dispatch, and SKAT-O rho extraction are all ready for
  mocked validation in Plan 20-03
- No blockers for Plan 20-03

---
*Phase: 20-r-skat-backend*
*Completed: 2026-02-20*
