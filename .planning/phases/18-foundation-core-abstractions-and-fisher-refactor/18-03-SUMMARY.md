---
phase: 18-foundation-core-abstractions-and-fisher-refactor
plan: "03"
subsystem: cli, output
tags: [cli, argparse, excel, association, fisher, stage-pipeline]

# Dependency graph
requires:
  - phase: 18-02
    provides: AssociationAnalysisStage wired into pipeline_core with association_output config key
provides:
  - CLI flags --perform-association, --association-tests, --skat-backend accepted and validated
  - context.config keys perform_association, association_tests, skat_backend set from CLI
  - ExcelReportStage generates "Association" sheet when perform_association + association_output active
affects:
  - 18-04 (parity tests will exercise CLI parsing and output path)
  - 19-onwards (all downstream phases use perform_association CLI flag)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Soft dependency pattern: ExcelReportStage.soft_dependencies includes both gene_burden_analysis and association_analysis"
    - "Association sheet mirrors Gene Burden sheet pattern verbatim (guards: config flag + file exists + size > 0)"
    - "Validation pattern: --association-tests without --perform-association triggers parser.error()"

key-files:
  created: []
  modified:
    - variantcentrifuge/cli.py
    - variantcentrifuge/stages/output_stages.py

key-decisions:
  - "Used parser.error() (not sys.exit) for --association-tests without --perform-association, consistent with argparse conventions"
  - "skat_backend uses getattr(args, 'skat_backend', 'auto') defensively in config mapping"
  - "Association sheet block is exact mirror of Gene Burden sheet block for maintainability"

patterns-established:
  - "New analysis flags: add to stats_group in create_parser(), map to cfg[] after gene burden block in main()"
  - "New Excel sheet: add to soft_dependencies + add sheet block in _add_additional_sheets() following Gene Burden pattern"

# Metrics
duration: 17min
completed: 2026-02-19
---

# Phase 18 Plan 03: CLI + Excel Association Interface Summary

**--perform-association CLI flag with validation + Association sheet in ExcelReportStage mirroring Gene Burden pattern**

## Performance

- **Duration:** 17 min
- **Started:** 2026-02-19T08:03:27Z
- **Completed:** 2026-02-19T08:20:30Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Users can invoke `--perform-association --association-tests fisher` from CLI
- Invalid combination `--association-tests` without `--perform-association` raises a clear parser error
- ExcelReportStage gains "Association" sheet generation with same guard logic as "Gene Burden" sheet
- Both Gene Burden and Association sheets can coexist in the same Excel output when both flags are active
- 858 unit tests pass with zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Add CLI arguments for association framework** - `31a3c3d` (feat)
2. **Task 2: Add Association sheet to ExcelReportStage** - `8cea44e` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified

- `variantcentrifuge/cli.py` - Added --perform-association, --association-tests, --skat-backend args; config mapping; validation
- `variantcentrifuge/stages/output_stages.py` - association_analysis in soft_dependencies; Association sheet block in _add_additional_sheets()

## Decisions Made

- Used `parser.error()` (not `sys.exit`) for the `--association-tests` without `--perform-association` validation — argparse convention produces correctly formatted usage message.
- `getattr(args, 'skat_backend', 'auto')` defensive access in cfg mapping to future-proof against namespace differences.
- Association sheet block is a verbatim mirror of Gene Burden sheet block — explicit duplication preferred over abstraction for parallel maintainability.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- CLI interface for association framework is complete: --perform-association, --association-tests, --skat-backend all wired
- ExcelReportStage ready to produce Association sheet when AssociationAnalysisStage produces output
- Plan 18-04 (parity tests: FisherExactTest vs gene_burden.py cross-validation) can proceed immediately
- All Phase 18 must-haves for CLI and Excel output are satisfied

---
*Phase: 18-foundation-core-abstractions-and-fisher-refactor*
*Completed: 2026-02-19*
