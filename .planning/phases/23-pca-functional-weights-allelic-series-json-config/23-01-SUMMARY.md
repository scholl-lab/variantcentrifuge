---
phase: 23-pca-functional-weights-allelic-series-json-config
plan: 01
subsystem: association
tags: [pca, population-stratification, akt, plink, eigenvec, covariates, numpy, subprocess]

# Dependency graph
requires:
  - phase: 19-covariate-system-burden-tests
    provides: load_covariates() pattern, AssociationConfig, covariate_matrix handling in AssociationAnalysisStage
  - phase: 22-acat-o-diagnostics
    provides: AssociationConfig with diagnostics fields, base for PCA field extension
provides:
  - PCA file loading with auto-detection of PLINK .eigenvec, AKT stdout, and generic TSV formats
  - merge_pca_covariates() for combining PCA + existing covariate matrices
  - PCAComputationStage wrapping AKT subprocess (registered in processing stage registry)
  - CLI args: --pca-file, --pca-tool, --pca-components
  - AssociationAnalysisStage PCA integration (inline merge into covariate_matrix)
  - 35 unit tests covering all formats, alignment, merge, and stage behavior
affects:
  - 23-02 (functional weights), 23-03 (allelic series) — AssociationConfig PCA fields available

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "PCA file auto-detection via _detect_pca_format() using first-line heuristics"
    - "Sample alignment: reindex to vcf_samples order; ValueError on missing samples"
    - "PCA merged inline into local covariate_matrix (NOT stored in PipelineContext)"
    - "ToolNotFoundError (hard error) when AKT binary missing — no silent skip"

key-files:
  created:
    - variantcentrifuge/association/pca.py
    - tests/unit/test_pca.py
  modified:
    - variantcentrifuge/association/base.py
    - variantcentrifuge/stages/processing_stages.py
    - variantcentrifuge/stages/stage_registry.py
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/cli.py

key-decisions:
  - "PCA matrix merged inline into local covariate_matrix — never stored in context.config or PipelineContext"
  - "Missing AKT binary is a hard ToolNotFoundError — not a skip/warning"
  - "PLINK eigenvec: always use IID (column 2) not FID (column 1) for sample ID"
  - "Format detection via first-line heuristics: #FID/FID header vs two-column vs single-column"

patterns-established:
  - "Format auto-detection: read first 5 lines, check column count and numeric-ness of cells"
  - "Sample alignment pattern mirrors covariates.py: reindex, assert no NaN, ValueError on missing"

# Metrics
duration: 13min
completed: 2026-02-21
---

# Phase 23 Plan 01: PCA Integration Summary

**PCA file loading (PLINK .eigenvec, AKT output, generic TSV) with auto-format detection, sample alignment validation, and AKT subprocess stage; merged inline into covariate matrix for association tests**

## Performance

- **Duration:** 13 min
- **Started:** 2026-02-21T22:09:03Z
- **Completed:** 2026-02-21T22:22:21Z
- **Tasks:** 2/2
- **Files modified:** 7

## Accomplishments

- Created `variantcentrifuge/association/pca.py` with `load_pca_file()` handling 3 file formats, `merge_pca_covariates()`, and private `_detect_pca_format()` heuristic
- Extended `AssociationConfig` with `pca_file`, `pca_tool`, `pca_components` fields (Phase 23 section)
- Added `PCAComputationStage` wrapping AKT subprocess, registered in processing stage registry as `pca_computation`
- Wired PCA loading inline into `AssociationAnalysisStage._process()` after covariate loading; merged directly into local `covariate_matrix` variable
- Added 35 unit tests covering all formats, alignment, merge behavior, `>20` PC warning, and `ToolNotFoundError`
- 1432 unit tests pass (up from 1249 baseline; no regressions)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create pca.py module, extend AssociationConfig, add CLI args** - `03d7a9e` (feat)
2. **Task 2: PCAComputationStage, stage integration, and unit tests** - `707861b` (feat)

**Plan metadata:** (to be added with docs commit)

## Files Created/Modified

- `variantcentrifuge/association/pca.py` (created) — `load_pca_file()`, `merge_pca_covariates()`, `_detect_pca_format()`
- `variantcentrifuge/association/base.py` (modified) — Added `pca_file`, `pca_tool`, `pca_components` to `AssociationConfig`
- `variantcentrifuge/stages/processing_stages.py` (modified) — Added `PCAComputationStage` at end of file
- `variantcentrifuge/stages/stage_registry.py` (modified) — Registered `PCAComputationStage` in `_register_processing_stages()`
- `variantcentrifuge/stages/analysis_stages.py` (modified) — PCA integration block after covariate loading; PCA fields in `assoc_config`; renamed `cov_names` to `covariate_col_names`
- `variantcentrifuge/cli.py` (modified) — `--pca-file`, `--pca-tool`, `--pca-components` args; validation; cfg wiring
- `tests/unit/test_pca.py` (created) — 35 unit tests in 7 test classes

## Decisions Made

- **PCA matrix inline merge:** PCA columns merged directly into the local `covariate_matrix` variable in `_process()` — they feed into `gene_data["covariate_matrix"]` naturally without extra storage in `context.config`. Avoids the "genotype matrix never in PipelineContext" memory architecture invariant.
- **ToolNotFoundError is hard error:** Per CONTEXT.md decision, missing `akt` binary raises `ToolNotFoundError` immediately, not a logged skip. This prevents silent failures when users expect AKT to run.
- **IID not FID for PLINK:** PLINK .eigenvec uses FID as family ID (not sample-level unique); IID is the VCF-level sample identifier. Always use column index 1 (IID).
- **Format detection heuristic:** Two non-numeric leading columns → `plink_nohdr`; one non-numeric leading column → `akt_or_generic`; header starting with `#FID`/`FID` → `plink_header`. Default: `generic`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed test content construction bug in test_more_than_20_components_warns**
- **Found during:** Task 2 (unit test execution)
- **Issue:** Test content rows used `f"0.{s}{i:02d}"` where `s` was a character (e.g. 'A'), producing non-numeric `0.A01` values which were parsed as NaN by pandas, triggering an internal assertion
- **Fix:** Changed to `f"0.{j}{i:02d}"` where `j` is the enumerate index (integer), producing valid numeric values
- **Files modified:** `tests/unit/test_pca.py`
- **Verification:** 35/35 tests pass
- **Committed in:** `707861b`

**2. [Rule 3 - Blocking] Fixed 4 ruff lint errors in test file (SIM117, RUF005)**
- **Found during:** Task 2 (lint check)
- **Issue:** Nested `with` statements flagged by SIM117; list concatenation `[header] + rows` flagged by RUF005
- **Fix:** Combined nested `with` statements; replaced `+` list concatenation with unpacking `[header, *rows]`
- **Files modified:** `tests/unit/test_pca.py`
- **Verification:** `make lint && make format` clean
- **Committed in:** `707861b`

---

**Total deviations:** 2 auto-fixed (1 bug in test data, 1 lint compliance)
**Impact on plan:** Both minor; no scope change. Plan executed as specified.

## Issues Encountered

None beyond the two auto-fixed deviations above.

## User Setup Required

None - no external service configuration required. AKT must be installed by the user if `--pca-tool akt` is used.

## Next Phase Readiness

- PCA infrastructure complete: `load_pca_file()`, `merge_pca_covariates()`, `PCAComputationStage`, CLI args all ready
- `AssociationConfig` has PCA fields populated from context.config in `AssociationAnalysisStage`
- Phase 23 Plan 02 (functional weights: CADD/REVEL) can proceed immediately
- Architecture invariant maintained: PCA matrix never stored in PipelineContext

---
*Phase: 23-pca-functional-weights-allelic-series-json-config*
*Completed: 2026-02-21*
