---
phase: 31-coast-fix
plan: 02
subsystem: association
tags: [coast, allelic-series, scoring, formula-engine, classification, cli, snpeff]

# Dependency graph
requires:
  - phase: 31-coast-fix/31-01
    provides: PurePythonCOASTTest and COASTTest with partial-category fix
  - phase: 24-pure-python-coast-backend
    provides: PythonCOASTBackend used by classify_variants formula engine path
provides:
  - COAST-02: auto-field injection adds required annotation fields when COAST selected
  - COAST-04: _resolve_effect() resolves '&'-concatenated SnpEff multi-transcript effect strings
  - COAST-05: --coast-classification CLI option selects classification model by name
  - COAST-06: scoring/coast_classification/ with sift_polyphen, cadd, and custom_template models
  - COAST-07: _resolve_effect implements PTV > missense > first priority resolution
  - Formula engine path in classify_variants: delegates to read_scoring_config + apply_scoring
  - coast_classification field on AssociationConfig propagated from cfg
  - 15 unit tests for effect resolution, formula engine, CLI option, and field injection
affects: [35-weighted-coast, 36-sparse-matrices]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "COAST classification via formula engine: scoring/coast_classification/{model}/ configs delegate to apply_scoring()"
    - "Column normalization before formula engine: COAST_EFFECT/COAST_IMPACT/COAST_SIFT/COAST_POLYPHEN/COAST_CADD"
    - "Effect string resolution: _resolve_effect() applied in both hardcoded and formula engine paths"
    - "CLI auto-injection: COAST fields injected before cfg['fields_to_extract'] is set"

key-files:
  created:
    - scoring/coast_classification/sift_polyphen/variable_assignment_config.json
    - scoring/coast_classification/sift_polyphen/formula_config.json
    - scoring/coast_classification/cadd/variable_assignment_config.json
    - scoring/coast_classification/cadd/formula_config.json
    - scoring/coast_classification/custom_template/variable_assignment_config.json
    - scoring/coast_classification/custom_template/formula_config.json
    - tests/unit/test_coast_classification.py
  modified:
    - variantcentrifuge/association/tests/allelic_series.py
    - variantcentrifuge/association/tests/allelic_series_python.py
    - variantcentrifuge/association/base.py
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/cli.py

key-decisions:
  - "Variable assignment configs use COAST_* normalized column names (not raw VCF names); classify_variants normalizes before calling apply_scoring"
  - "coast_classification in AssociationConfig stores absolute path (None = use hardcoded logic); cli.py resolves model name to path"
  - "Auto-injection filters out COAST_* internal names from vcf_fields since they map to normalized columns, not raw VCF INFO fields"
  - "_classify_via_formula_engine is a private helper; classify_variants() is the single public API"
  - "CADD_COLUMN_CANDIDATES added to allelic_series.py for auto-detection of CADD column variants"

patterns-established:
  - "Column normalization before formula engine: build work DataFrame with COAST_* prefixed standard names"
  - "Model resolution: model name -> absolute path in cli.py; path stored in cfg['coast_classification']"
  - "Formula engine fallback safety: if model_dir load fails, return zeros (log error, don't crash pipeline)"

# Metrics
duration: 17min
completed: 2026-02-23
---

# Phase 31 Plan 02: COAST Classification Config and Effect Resolution Summary

**Configurable COAST variant classification via formula engine (sift_polyphen/cadd/custom), with '&'-concatenated SnpEff effect resolution and auto-injection of required annotation fields**

## Performance

- **Duration:** ~17 min
- **Started:** 2026-02-23T18:07:03Z
- **Completed:** 2026-02-23T18:24:14Z
- **Tasks:** 2
- **Files modified:** 11 (7 created + 4 modified source + 1 new test)

## Accomplishments

- COAST-04/07: `_resolve_effect()` resolves `stop_gained&splice_region_variant` to `stop_gained` (PTV > missense > first); applied in both hardcoded and formula engine paths
- COAST-06: Three classification configs in `scoring/coast_classification/` — `sift_polyphen` (replicates prior hardcoded logic), `cadd` (CADD>=15 → DMV, 0<CADD<15 → BMV), `custom_template` (annotated placeholder)
- COAST-02/05: `--coast-classification` CLI option auto-injects required annotation fields into field extraction when COAST test is selected; model name resolved to absolute path stored in `AssociationConfig.coast_classification`
- Both `COASTTest.run()` and `PurePythonCOASTTest.run()` pass `model_dir` from `config.coast_classification` to `classify_variants()`
- All 142 COAST tests pass, 1584 unit tests pass, CI green

## Task Commits

1. **Task 1: Classification configs, _resolve_effect, formula engine path** - `bf7d828` (feat)
2. **Task 2: CLI option, auto-field injection, AssociationConfig, and tests** - `593121f` (feat)

**Plan metadata:** (see docs commit below)

## Files Created/Modified

- `scoring/coast_classification/sift_polyphen/variable_assignment_config.json` - Maps COAST_EFFECT/IMPACT/SIFT/POLYPHEN to formula variables
- `scoring/coast_classification/sift_polyphen/formula_config.json` - Replicates hardcoded SIFT+PolyPhen logic; outputs coast_category
- `scoring/coast_classification/cadd/variable_assignment_config.json` - Maps COAST_EFFECT/IMPACT/CADD to formula variables
- `scoring/coast_classification/cadd/formula_config.json` - CADD threshold-based BMV/DMV/PTV classification
- `scoring/coast_classification/custom_template/variable_assignment_config.json` - Annotated template for custom models
- `scoring/coast_classification/custom_template/formula_config.json` - Template with IMPACT-based placeholder and comments
- `variantcentrifuge/association/tests/allelic_series.py` - Added `_resolve_effect()`, `_classify_via_formula_engine()`, updated `classify_variants()` signature with model_dir/diagnostics_rows, added CADD_COLUMN_CANDIDATES
- `variantcentrifuge/association/tests/allelic_series_python.py` - Updated classify_variants() call to pass model_dir from config
- `variantcentrifuge/association/base.py` - Added `coast_classification: str | None` field to AssociationConfig
- `variantcentrifuge/stages/analysis_stages.py` - Added `coast_classification` to `_build_assoc_config_from_context()`
- `variantcentrifuge/cli.py` - Added `--coast-classification` argument, auto-field injection block, cfg propagation
- `tests/unit/test_coast_classification.py` - 15 unit tests covering all COAST-02/04/05/06/07 features

## Decisions Made

- **Column normalization approach:** The scoring configs use `COAST_*` prefixed normalized column names instead of raw VCF field names (e.g., `ANN[0].EFFECT`). `classify_variants()` builds a normalized working copy before calling `apply_scoring()`. This avoids column sanitization issues (sanitized names like `ANN_0__EFFECT` vs raw `ANN[0].EFFECT`) and decouples the formula engine from VCF field naming conventions.
- **coast_classification stores absolute path:** cli.py resolves the model name to an absolute path and stores it in `cfg["coast_classification"]`. AssociationConfig stores the path (or None for default). This avoids path resolution at classification time inside test.run() methods.
- **Auto-injection filters COAST_* names:** The variable_assignment_config keys are COAST_* internal names (not raw VCF fields), so the auto-injection block filters them out — no raw VCF fields are injected for the built-in models. This is intentional; users with custom models using raw VCF field names would get those injected.
- **FutureWarning fix for CADD column:** Used `.where(notna(), other=0)` instead of `.fillna(0)` to avoid pandas FutureWarning about downcasting object dtype arrays.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed FutureWarning in CADD column numeric conversion**

- **Found during:** Task 1 (testing CADD model)
- **Issue:** `pd.to_numeric(series.fillna(0), errors='coerce')` produced FutureWarning about object dtype downcasting
- **Fix:** Changed to `.where(series.notna(), other=0)` to avoid the intermediate `.fillna(0)` on object dtype
- **Files modified:** `variantcentrifuge/association/tests/allelic_series.py`
- **Verification:** No FutureWarning emitted when tests run with `warnings.filterwarnings('error')`
- **Committed in:** `bf7d828` (Task 1 commit)

**2. [Rule 1 - Bug] Variable assignment configs redesigned to use COAST_* normalized names**

- **Found during:** Task 1 (designing variable_assignment_config)
- **Issue:** Original plan used `ANN[0].EFFECT` as config keys, but at classify_variants() call time, column names are sanitized to `ANN_0__EFFECT` or `EFFECT`. The scoring engine looks for exact key matches, so config keys would not match any column.
- **Fix:** Redesigned configs to use `COAST_*` internal names; `_classify_via_formula_engine()` builds a normalized working copy with these standard names before calling `apply_scoring()`
- **Files modified:** All 6 `variable_assignment_config.json` files, `allelic_series.py`
- **Verification:** Both models produce correct codes on test DataFrames; sift_polyphen matches hardcoded logic exactly
- **Committed in:** `bf7d828` (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (Rule 1 - Bug: FutureWarning fix; Rule 1 - Bug: column normalization redesign)
**Impact on plan:** Both fixes were necessary for correctness. The column normalization redesign is architecturally cleaner than the original plan's approach.

## Issues Encountered

- `os` not imported at top of cli.py (only locally imported via `import os as _os` inside the COAST injection block) — resolved by using `_os` alias
- Ruff I001 (unsorted imports) in test file — resolved with `ruff check --fix`
- Unused imports (`MagicMock`, `patch`) left from test scaffolding — removed

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- COAST classification is now fully configurable via `--coast-classification sift_polyphen|cadd|<custom>`
- `_resolve_effect()` handles real-world multi-transcript SnpEff annotations
- Phase 35 (weighted COAST) can build on this without modification (coast_classification is independent of coast_weights)
- Phase 35 may want to add diagnostics output (the `diagnostics_rows` parameter of `classify_variants()` is ready but not yet wired to file output)
- Custom model support: users can drop a directory under `scoring/coast_classification/` and pass its name via `--coast-classification`

---
*Phase: 31-coast-fix*
*Completed: 2026-02-23*
