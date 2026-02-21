---
phase: 22-acat-o-and-diagnostics
plan: 02
subsystem: association-analysis
tags: [diagnostics, lambda-gc, qq-plot, genomic-inflation, sample-size-warnings, cli]

# Dependency graph
requires:
  - phase: 22-01
    provides: AssociationConfig with min_cases/max_case_control_ratio/min_case_carriers/diagnostics_output fields; engine.run_all() returning results_df with acat_o_p_value
provides:
  - diagnostics.py module: compute_lambda_gc(), compute_qq_data(), emit_sample_size_warnings(), compute_per_gene_warnings(), write_diagnostics()
  - --diagnostics-output CLI argument with validation (requires --perform-association)
  - Per-gene warnings column in association TSV (semicolon-separated flag strings)
  - lambda_gc.tsv, qq_data.tsv, summary.txt diagnostics files
affects:
  - Phase 23: PCA + functional weights (can use diagnostics module for QC output)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Lazy import for diagnostics module in AssociationAnalysisStage (same pattern as covariates)"
    - "Per-gene warnings dict built from gene_burden_data before engine.run_all()"
    - "Hardcoded _EXPECTED_CHI2_MEDIAN constant avoids import-time scipy call"

key-files:
  created:
    - variantcentrifuge/association/diagnostics.py
  modified:
    - variantcentrifuge/cli.py
    - variantcentrifuge/stages/analysis_stages.py

key-decisions:
  - "gene_burden_data dicts use 'GENE' (uppercase) key; results_df uses 'gene' (lowercase)"
  - "Hazen quantile formula i/(n+1) for expected QQ quantiles"
  - "QQ data sorted ascending by expected_neg_log10_p (non-significant end first)"
  - "Cohort warnings emitted by emit_sample_size_warnings() and also written to summary.txt"
  - "Per-gene warnings column added to results_df before TSV write so it appears in output"

patterns-established:
  - "Diagnostics: three-file output pattern (lambda_gc.tsv, qq_data.tsv, summary.txt)"
  - "Warning flags: semicolon-separated uppercase strings in 'warnings' column"

# Metrics
duration: 9min
completed: 2026-02-21
---

# Phase 22 Plan 02: Diagnostics Module Summary

**lambda_GC genomic inflation factor, QQ TSV data, sample size warnings, and per-gene carrier warnings integrated into association pipeline via --diagnostics-output CLI arg**

## Performance

- **Duration:** 9 min
- **Started:** 2026-02-21T17:54:30Z
- **Completed:** 2026-02-21T18:04:15Z
- **Tasks:** 2
- **Files modified:** 3 (1 created, 2 modified)

## Accomplishments

- Created `variantcentrifuge/association/diagnostics.py` with five functions covering all DIAG requirements
- Added `--diagnostics-output` CLI argument with proper validation (requires `--perform-association`)
- Integrated per-gene warnings column into association TSV output
- Integrated diagnostics file writing into `AssociationAnalysisStage._process()` after TSV write
- All 1622 tests pass, CI checks green

## Task Commits

Each task was committed atomically:

1. **Task 1: Diagnostics module** - `e868a9f` (feat)
2. **Task 2: CLI arg + stage integration** - `9ec9b66` (feat)

**Plan metadata:** (forthcoming in docs commit)

## Files Created/Modified

- `variantcentrifuge/association/diagnostics.py` - Full diagnostics module: lambda_GC computation, QQ data generation, sample size warnings, per-gene warnings, and write_diagnostics() orchestrator
- `variantcentrifuge/cli.py` - Added `--diagnostics-output` argument and config mapping + validation
- `variantcentrifuge/stages/analysis_stages.py` - AssociationConfig extended with Phase 22 thresholds; per-gene warnings dict; diagnostics write after TSV; lazy import of diagnostics functions

## Decisions Made

- **gene_burden_data key is "GENE" (uppercase):** The dict returned by all three aggregation paths (`_aggregate_gene_burden_from_columns`, `_aggregate_gene_burden_from_gt`, `_aggregate_gene_burden_legacy`) uses `"GENE"` uppercase. The results DataFrame column is `"gene"` (lowercase, set by engine). Used `d.get("GENE", d.get("gene", ""))` for robustness.
- **QQ sort ascending:** Sorted by `expected_neg_log10_p` ascending so the non-significant end (small values) is first — matches typical sequential rendering order for QQ plots.
- **Hardcoded chi2 median:** `_EXPECTED_CHI2_MEDIAN = 0.45493642311957174` avoids calling `chi2.ppf()` at import time, which would add ~20ms cold import overhead.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed KeyError on "gene" key in gene_burden_data lookup**
- **Found during:** Task 2 (stage integration), discovered by test failure
- **Issue:** Plan specified `{d["gene"]: d for d in gene_burden_data}` but gene_burden_data dicts use uppercase `"GENE"` key (all three aggregation paths set this). `KeyError: 'gene'` broke test `test_sets_association_results_on_context`.
- **Fix:** Changed to `{d.get("GENE", d.get("gene", "")): d for d in gene_burden_data}` for robustness across both casings
- **Files modified:** variantcentrifuge/stages/analysis_stages.py
- **Verification:** `pytest tests/unit/test_association_stage.py -x -q` — all 21 pass
- **Committed in:** 9ec9b66 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug)
**Impact on plan:** Necessary correctness fix. No scope creep.

## Issues Encountered

- Ruff lint flagged three issues in diagnostics.py: two `UP037` (quoted type annotations under `from __future__ import annotations` are redundant) and one `F541` (f-string without placeholders). Fixed with `ruff check --fix` and `ruff format`.

## User Setup Required

None - no external service configuration required. All diagnostics are written to a local directory specified by `--diagnostics-output`.

## Next Phase Readiness

- Phase 22 complete: ACAT-O + all diagnostics (lambda_GC, QQ data, sample size warnings, per-gene warnings, summary.txt)
- Phase 23 (PCA + Functional Weights + Allelic Series + JSON Config) can proceed
- The diagnostics module can be extended in Phase 23 for matplotlib QQ plots (DIAG-04)
- No blockers

---
*Phase: 22-acat-o-and-diagnostics*
*Completed: 2026-02-21*
