---
phase: 32-region-restriction-and-pca-wiring
plan: 02
subsystem: association
tags: [pca, akt, pipeline, cli, argparse, subprocess, covariate-adjustment]

# Dependency graph
requires:
  - phase: 32-01
    provides: RegionRestrictionStage and --regions-bed CLI flag
  - phase: 23-pca-functional-weights-allelic-series-json-config
    provides: pca.py module with load_pca_file/merge_pca_covariates and AssociationConfig.pca_file

provides:
  - PCAComputationStage class in processing_stages.py with AKT subprocess and file pass-through
  - Unified --pca CLI flag replacing --pca-file and --pca-tool with backward compat aliases
  - PCAComputationStage wired into build_pipeline_stages before AssociationAnalysisStage
  - AssociationAnalysisStage reads pca_file from pca_computation stage result
  - "pca" added to VALID_ASSOCIATION_KEYS for JSON config support
  - 9 unit tests covering all PCA stage code paths

affects:
  - phase: 35-weighted-coast-and-diagnostics
  - phase: 36-sparse-matrix-optimization

# Tech tracking
tech-stack:
  added: []
  patterns:
    - PCA computation as a pipeline stage with cache reuse (skip AKT subprocess if eigenvec file exists)
    - Unified flag with hidden deprecated aliases using argparse dest= sharing
    - Soft dependency: AssociationAnalysisStage uses pca_computation result via setdefault (non-destructive)

key-files:
  created:
    - .planning/phases/32-region-restriction-and-pca-wiring/32-02-SUMMARY.md
  modified:
    - variantcentrifuge/stages/processing_stages.py
    - variantcentrifuge/cli.py
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/pipeline.py
    - tests/unit/stages/test_processing_stages.py

key-decisions:
  - "--pca unified flag uses dest='pca' for --pca-file and --pca-tool aliases; no argparse conflict"
  - "PCAComputationStage sets cfg['pca_file'] for backward compat with _build_assoc_config_from_context"
  - "AssociationAnalysisStage uses setdefault so explicit pca_file= in config.json still overrides stage result"
  - "Cache reuse: skip AKT subprocess if output eigenvec file already exists and is non-empty"
  - "pca_computation added as soft_dependency not hard dependency (PCA optional feature)"

patterns-established:
  - "Stage output handoff: stage sets context.config['pca_file'] AND marks complete with result dict"
  - "Downstream stage reads result via context.get_result() and uses setdefault to not override explicit config"

# Metrics
duration: 22min
completed: 2026-02-23
---

# Phase 32 Plan 02: PCA Wiring Summary

**PCAComputationStage wired into pipeline — unified --pca flag auto-detects file path vs 'akt' tool, with AKT subprocess caching and AssociationAnalysisStage soft dependency**

## Performance

- **Duration:** 22 min
- **Started:** 2026-02-23T20:32:06Z
- **Completed:** 2026-02-23T20:53:47Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments

- PCAComputationStage added to processing_stages.py with `_run_akt` subprocess invocation, pre-computed file pass-through, and output caching
- Unified `--pca` CLI flag replaces `--pca-file` and `--pca-tool`; deprecated aliases use `argparse dest='pca'` for zero-friction backward compatibility
- PCAComputationStage wired into `build_pipeline_stages()` before AssociationAnalysisStage; added to `create_stages_from_config()` for JSON config support
- AssociationAnalysisStage picks up `pca_file` from `pca_computation` stage result via soft dependency; uses `setdefault` to not override explicit config
- 9 unit tests cover all PCA stage paths (file, akt, cache, invalid, errors, properties)

## Task Commits

1. **Task 1: Create PCAComputationStage and unify --pca CLI flag** - `29e411d` (feat)
2. **Task 2: Wire PCAComputationStage into pipeline and AssociationAnalysisStage** - `0136034` (feat)

**Plan metadata:** (committed with summary)

## Files Created/Modified

- `variantcentrifuge/stages/processing_stages.py` - Added PCAComputationStage class at end of file
- `variantcentrifuge/cli.py` - Replaced --pca-file/--pca-tool with unified --pca; hidden compat aliases; cfg['pca'] mapping; single validation check
- `variantcentrifuge/stages/analysis_stages.py` - Added 'pca_computation' to soft_dependencies; added pca_file pickup in _process; added 'pca' to VALID_ASSOCIATION_KEYS
- `variantcentrifuge/pipeline.py` - Added PCAComputationStage import and wiring in build_pipeline_stages; args.pca in create_stages_from_config
- `tests/unit/stages/test_processing_stages.py` - Added TestPCAComputationStage with 9 tests

## Decisions Made

- **Unified --pca flag**: Single argument auto-detects file path (`os.path.isfile`) vs tool name string; 'akt' is the only recognized tool. Invalid values raise `ValueError` with clear message.
- **Backward compat via argparse dest sharing**: `--pca-file` and `--pca-tool` both use `dest='pca'` and `help=argparse.SUPPRESS`. No code changes required downstream since args.pca_file/args.pca_tool never existed in cfg after Phase 32-01 — cfg['pca'] is the canonical key now.
- **Stage sets cfg['pca_file'] for downstream compat**: `_build_assoc_config_from_context` reads `pca_file` key, so PCAComputationStage writes to that key. AssociationAnalysisStage also reads stage result via `get_result()` as belt-and-suspenders.
- **setdefault not assignment**: AssociationAnalysisStage uses `context.config.setdefault('pca_file', ...)` so explicit `pca_file:` in JSON association config still takes precedence over the computed stage result.
- **AKT cache**: If `{base_name}.pca.eigenvec` exists and is non-empty, skip subprocess call — important for pipeline resume/checkpoint use cases.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Dead code in processing_stages.py after DataSortingStage._process**

- **Found during:** Task 1 (adding PCAComputationStage to processing_stages.py)
- **Issue:** Initial Edit of `processing_stages.py` matched a non-unique pattern and appended PCAComputationStage in the wrong location, causing `DataSortingStage._process` to lose its `context.extracted_tsv = sorted_tsv` update lines and `return context`. The dead code also contained unreachable `sorted_tsv` references.
- **Fix:** Used `git checkout HEAD -- processing_stages.py` to restore the original, then correctly appended PCAComputationStage after `DataSortingStage.get_output_files` at the actual end of the file.
- **Files modified:** variantcentrifuge/stages/processing_stages.py
- **Verification:** `grep -n "class.*Stage"` confirmed correct class/method structure; lint and tests passed
- **Committed in:** 29e411d (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug: Edit tool target collision)
**Impact on plan:** Required git restore and re-apply of change. No scope creep.

## Issues Encountered

Edit tool pattern collision: the Edit tool matched ambiguous context lines inside `DataSortingStage._process` rather than at the end of the file. Resolution: used `git checkout` to restore the file, then appended PCAComputationStage at the definitive end-of-file position (after `DataSortingStage.get_output_files`) to avoid any ambiguity.

## User Setup Required

None - no external service configuration required. AKT tool must be in PATH when `--pca akt` is used, but this is validated at runtime with a clear error message.

## Next Phase Readiness

- PCA wiring complete: `--pca akt` auto-computes eigenvectors; `--pca /path/file` uses pre-computed file
- Phase 32 plans complete (01: region restriction, 02: PCA wiring)
- Ready for Phase 33+ (FDR weighting, case-confidence)
- No blockers introduced

---
*Phase: 32-region-restriction-and-pca-wiring*
*Completed: 2026-02-23*
