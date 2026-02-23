---
phase: 29-classic-pipeline-deprecation-and-removal
plan: 01
subsystem: pipeline
tags: [pipeline, refactoring, naming, cleanup]

# Dependency graph
requires:
  - phase: 27-performance-optimizations
    provides: Stage-based pipeline complete and production-ready; ProcessPoolExecutor parallelization

provides:
  - run_pipeline() as the canonical stage-based pipeline entry point (renamed from run_refactored_pipeline)
  - pipeline.py with clean module docstring and no vestigial demo code
  - All test and production imports updated to use run_pipeline

affects:
  - 29-02: Remove --use-new-pipeline flag (follows from this naming cleanup)
  - 29-03: Documentation cleanup

# Tech tracking
tech-stack:
  added: []
  patterns: []

key-files:
  created: []
  modified:
    - variantcentrifuge/pipeline.py
    - variantcentrifuge/cli.py
    - tests/integration/test_pipeline_with_mocked_tools.py
    - tests/test_gene_list_integration.py
    - tests/unit/test_cli_debug_logging.py

key-decisions:
  - "Retain local variable refactored_args in cli.py — user-facing impact is nil, rename would touch many lines for no benefit"
  - "Remove main() demo entirely — vestigial, not called from production code, duplicates CLI"

patterns-established:
  - "pipeline.py = THE stage-based pipeline; no legacy duality in naming"

# Metrics
duration: 7min
completed: 2026-02-23
---

# Phase 29 Plan 01: Pipeline Naming Cleanup Summary

**Renamed run_refactored_pipeline to run_pipeline everywhere (pipeline.py, cli.py, 3 test files) and removed the vestigial main() demo function from pipeline.py**

## Performance

- **Duration:** ~7 min
- **Started:** 2026-02-23T06:32:00Z
- **Completed:** 2026-02-23T06:46:30Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- Updated module docstring in pipeline.py from "Refactored pipeline using new Stage architecture" to clean "Stage-based pipeline orchestration" description
- Renamed `run_refactored_pipeline()` to `run_pipeline()` in pipeline.py with updated docstring
- Removed vestigial `main()` demo function and `if __name__ == "__main__":` block from pipeline.py (21 lines removed)
- Updated cli.py import and call site (line 12 and call at ~1665)
- Updated 3 test files: test_pipeline_with_mocked_tools.py (7 occurrences), test_gene_list_integration.py (mock function + monkeypatch target), test_cli_debug_logging.py (2 @patch decorators)
- 2001 tests pass, lint clean

## Task Commits

Each task was committed atomically:

1. **Task 1: Rename pipeline.py functions and clean up docstrings** - `44ab13b` (refactor)
2. **Task 2: Update test imports referencing run_refactored_pipeline** - already incorporated (tests were clean after Task 1 commit was applied; prior work in branch)

**Plan metadata:** TBD (docs commit)

## Files Created/Modified
- `variantcentrifuge/pipeline.py` - Updated docstring, renamed function, removed main() demo
- `variantcentrifuge/cli.py` - Updated import and call site
- `tests/integration/test_pipeline_with_mocked_tools.py` - Updated import and 6 call sites
- `tests/test_gene_list_integration.py` - Updated mock function name, docstring, and monkeypatch target
- `tests/unit/test_cli_debug_logging.py` - Updated 2 @patch decorators

## Decisions Made
- Retain local variable name `refactored_args` in cli.py — it is not user-facing, renaming would touch many lines for no observable benefit and introduces bug risk with no reward.
- Remove `main()` demo entirely rather than deprecating it — it was never called from production code, purely a demonstration/dev-time artifact.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- The test files (Task 2) appeared already clean after the Task 1 commit was staged — the prior session had committed test changes in `bb83400`. No rework was needed; verification confirmed all references were gone.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- pipeline.py naming is now clean; ready for Plan 02 (remove --use-new-pipeline flag)
- No blockers or concerns

---
*Phase: 29-classic-pipeline-deprecation-and-removal*
*Completed: 2026-02-23*
