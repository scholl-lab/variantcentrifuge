---
phase: 29-classic-pipeline-deprecation-and-removal
plan: 02
subsystem: testing
tags: [cli, test-cleanup, argparse, gene-burden]

# Dependency graph
requires:
  - phase: 29-classic-pipeline-deprecation-and-removal
    provides: Research identifying --use-new-pipeline as documentation phantom never in argparse
provides:
  - tests/conftest.py without vestigial --use-new-pipeline flag in gene_burden_cmd_template
  - tests/test_gene_burden_comprehensive.py with 8 test commands free of the phantom flag
  - tests/test_cli.py with accurate docstrings describing actual test behavior
  - tests/fixtures/geneburden/README.md with clean example commands
affects:
  - 29-03 (classic pipeline removal — test suite must be clean before removing pipeline.py)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Test command strings only contain flags that argparse actually defines"
    - "Test docstrings describe what the test verifies, not flags that were never parsed"

key-files:
  created: []
  modified:
    - tests/conftest.py
    - tests/test_gene_burden_comprehensive.py
    - tests/test_cli.py
    - tests/fixtures/geneburden/README.md

key-decisions:
  - "Task 2 changes (test_cli.py + README) were already applied by prior Plan 29-03 executor — no double-commit needed"

patterns-established:
  - "Flag removal from test commands: delete the string, verify grep returns nothing, confirm tests pass/skip identically"

# Metrics
duration: 5min
completed: 2026-02-23
---

# Phase 29 Plan 02: Remove --use-new-pipeline from Tests Summary

**Removed vestigial --use-new-pipeline phantom flag from 9 test command strings across conftest.py and test_gene_burden_comprehensive.py; corrected 2 test docstrings and 5 README examples**

## Performance

- **Duration:** ~5 min
- **Started:** 2026-02-23T06:39:14Z
- **Completed:** 2026-02-23T06:44:18Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Deleted `"--use-new-pipeline",` from `gene_burden_cmd_template` fixture in `tests/conftest.py`
- Deleted `"--use-new-pipeline",` from all 8 test command lists in `tests/test_gene_burden_comprehensive.py`
- Updated 2 docstrings in `tests/test_cli.py` to accurately describe behavior (flag never existed)
- Removed `--use-new-pipeline` from 5 example commands in `tests/fixtures/geneburden/README.md`
- Verified 2001 tests pass, 3 skip (missing test data), 0 failures after changes

## Task Commits

Each task was committed atomically:

1. **Task 1: Remove --use-new-pipeline from test command strings** - `bb83400` (chore)
2. **Task 2: Update test_cli.py docstrings and fixture README** - already in `1559843` (applied by prior 29-03 executor)

**Plan metadata:** (docs commit follows)

## Files Created/Modified
- `tests/conftest.py` - Removed `"--use-new-pipeline",` from `gene_burden_cmd_template` fixture
- `tests/test_gene_burden_comprehensive.py` - Removed flag from 8 test command lists
- `tests/test_cli.py` - Updated 2 docstrings: test_show_checkpoint_status_with_new_pipeline_flag and test_show_checkpoint_status_flag_precedence
- `tests/fixtures/geneburden/README.md` - Removed flag from 5 example commands (Quick Start + 4 Test Scenarios)

## Decisions Made
- Task 2 changes were already present in HEAD (commit `1559843`) from the prior 29-03 executor run. Rather than creating an empty commit, the task is documented here as already complete. The `grep -rn` verification confirms zero references remain.

## Deviations from Plan

None - plan executed exactly as written. Task 2 changes were pre-applied; no additional work was needed.

## Issues Encountered
- Task 2 changes (test_cli.py docstrings + README) were already in HEAD via commit `1559843` ("fix(29): complete remaining --use-new-pipeline removals") from the Plan 29-03 executor, which had gone slightly out of scope. No re-application was needed; `git status` confirmed nothing to stage.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All `--use-new-pipeline` references eliminated from tests/ source files
- Test suite passes with 2001 tests passing, 3 skipping (missing fixture data)
- Ready for Plan 29-03 completion verification and STATE.md update

---
*Phase: 29-classic-pipeline-deprecation-and-removal*
*Completed: 2026-02-23*
