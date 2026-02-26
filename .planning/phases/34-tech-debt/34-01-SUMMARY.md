---
phase: 34-tech-debt
plan: "01"
subsystem: pipeline
tags: [pipeline, association, gene-burden, diagnostics, lambda-gc, config]

# Dependency graph
requires:
  - phase: 33-gene-level-fdr-weighting
    provides: AssociationAnalysisStage and engine.run_all() as primary association entrypoints
provides:
  - create_stages_from_config() correctly activates AssociationAnalysisStage and GeneBurdenAnalysisStage from config
  - Integration test suite for config key -> stage mapping correctness
  - compute_lambda_gc() documented with Fisher exemption and over-correction risk guidance
affects:
  - 35-diagnostics-wiring (depends on association stage activation being correct)
  - any caller using create_stages_from_config() API

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "create_stages_from_config() docstring maps config keys to stages for maintainability"
    - "Integration tests guard config key -> stage activation regressions"

key-files:
  created:
    - tests/integration/test_create_stages_from_config.py
  modified:
    - variantcentrifuge/pipeline.py
    - variantcentrifuge/association/diagnostics.py

key-decisions:
  - "igv_report vs igv naming mismatch fixed: create_stages_from_config now sets args.igv (matching build_pipeline_stages check) not args.igv_report"
  - "Fisher lambda_GC computed for diagnostic output only; correction must not be applied to Fisher p-values"

patterns-established:
  - "Integration tests for config -> stage mapping guard against silent activation failures"

# Metrics
duration: 8min
completed: 2026-02-24
---

# Phase 34 Plan 01: Tech Debt Summary

**Silent config bugs fixed: create_stages_from_config() now correctly maps perform_association and perform_gene_burden keys, plus Fisher-exemption documented in compute_lambda_gc()**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-24T20:34:28Z
- **Completed:** 2026-02-24T20:42:51Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments

- Fixed TD-02: `create_stages_from_config({'perform_association': True})` now returns a stage list containing `AssociationAnalysisStage` (was silently ignored before)
- Fixed TD-02: `create_stages_from_config({'perform_gene_burden': True})` now returns `GeneBurdenAnalysisStage` (was silently ignored before)
- Fixed `igv_report` vs `igv` naming mismatch (create_stages_from_config set wrong attribute name)
- Added missing `no_stats` mapping
- Added docstring table in `create_stages_from_config()` mapping config keys to activated stages
- Added 8 integration tests covering all fixed config keys
- Documented TD-05: Fisher's exact test exemption from lambda_GC correction, over-correction risk, and inflation threshold guidance in `compute_lambda_gc()` docstring

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix create_stages_from_config mappings and add integration test (TD-02)** - `6db72ac` (fix)
2. **Task 2: Document lambda_GC Fisher exemption and over-correction risk (TD-05)** - `5202eff` (docs)

**Plan metadata:** (committed below)

## Files Created/Modified

- `variantcentrifuge/pipeline.py` - Added missing config key mappings (perform_association, perform_gene_burden, no_stats), fixed igv_report->igv mismatch, added docstring table
- `variantcentrifuge/association/diagnostics.py` - Extended compute_lambda_gc() docstring with Fisher exemption, applicable tests, inflation thresholds, over-correction risk; inline comment in write_diagnostics() loop
- `tests/integration/test_create_stages_from_config.py` - 8 integration tests for config key -> stage activation

## Decisions Made

- `igv_report` renamed to `igv` in create_stages_from_config: build_pipeline_stages() checks `getattr(args, "igv", False)` — the old `igv_report` attribute never matched. Fixed by aligning attribute name.
- Fisher lambda_GC is written to diagnostics output for all tests including Fisher — this is intentional (diagnostic metric), but callers must NOT use it for p-value correction. Documented via docstring Notes section and inline comment.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed igv_report -> igv naming mismatch**

- **Found during:** Task 1 (auditing create_stages_from_config against build_pipeline_stages)
- **Issue:** create_stages_from_config set `args.igv_report` but build_pipeline_stages checks `getattr(args, "igv", False)` — IGV report stage was never activated from config
- **Fix:** Changed attribute to `args.igv = config.get("igv", False)` to match check in build_pipeline_stages; also updated the docstring table to reflect `igv` key name
- **Files modified:** variantcentrifuge/pipeline.py
- **Verification:** test_igv_true_activates_igv_stage integration test passes
- **Committed in:** 6db72ac (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug)
**Impact on plan:** Necessary fix discovered during mandatory audit step in task. No scope creep.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 34 plan 01 complete; config key -> stage mapping is now correct and tested
- Phase 35 (diagnostics wiring) can proceed with confidence that association stages activate correctly from config
- No blockers

---
*Phase: 34-tech-debt*
*Completed: 2026-02-24*
