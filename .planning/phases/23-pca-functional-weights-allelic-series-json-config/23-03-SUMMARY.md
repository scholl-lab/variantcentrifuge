---
phase: 23-pca-functional-weights-allelic-series-json-config
plan: 03
subsystem: association
tags: [coast, allelic-series, rpy2, R, AllelicSeries, variant-classification, BMV, DMV, PTV, association-testing]

# Dependency graph
requires:
  - phase: 23-02
    provides: CADD/REVEL functional weights, variant_weight_params, AssociationConfig fields
  - phase: 22
    provides: ACAT-O omnibus framework that COAST feeds into
  - phase: 20
    provides: rpy2 pattern (RSKATTest) that COASTTest mirrors

provides:
  - COASTTest: R AllelicSeries::COAST() wrapper registered in engine as "coast"
  - classify_variants(): BMV/DMV/PTV classification from EFFECT/IMPACT/SIFT/PolyPhen columns
  - coast_weights: configurable via --coast-weights CLI and JSON config
  - gene_df: per-variant annotation DataFrame now stored in gene_data for COAST access
  - 43 unit tests covering classification, mocked R invocation, engine integration

affects:
  - 23-04 (JSON config phase: coast_weights in JSON config schema)
  - Future diagnostics (coast_burden_p_value, coast_skat_p_value diagnostic output)

# Tech tracking
tech-stack:
  added: [AllelicSeries (R package, runtime dependency)]
  patterns:
    - "COASTTest mirrors RSKATTest pattern: parallel_safe=False, periodic GC, prepare/finalize lifecycle"
    - "classify_variants auto-detects SIFT/PolyPhen column names from 6 candidates each"
    - "gene_df stored in gene_data dict for per-variant annotation access by tests beyond Fisher/SKAT"

key-files:
  created:
    - variantcentrifuge/association/tests/allelic_series.py
    - tests/unit/test_coast.py
  modified:
    - variantcentrifuge/association/tests/__init__.py
    - variantcentrifuge/association/base.py
    - variantcentrifuge/association/engine.py
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/cli.py
    - pyproject.toml

key-decisions:
  - "classify_variants returns include_mask: missense without SIFT/PolyPhen gets code 0, excluded from COAST only"
  - "COASTTest.run() returns p_value=None when gene missing any category (BMV/DMV/PTV), not spurious p"
  - "gene_df stored in gene_data dict (reset_index) for COAST annotation column access"
  - "Annotation/genotype mismatch detected by shape check: len(anno_codes) != geno.shape[1]"
  - "coast marker added to pyproject.toml for selective test running"

patterns-established:
  - "Annotation tests (beyond Fisher/SKAT) access per-variant columns via gene_data['gene_df']"
  - "classify_variants multi-value cell splitting on ';' or ',' for dbNSFP fields"

# Metrics
duration: 11min
completed: 2026-02-21
---

# Phase 23 Plan 03: COAST Allelic Series Test Summary

**COAST allelic series test with BMV/DMV/PTV classification via SnpEff/dbNSFP annotations, R AllelicSeries wrapper, ACAT-O feeding, configurable weights, and 43 unit tests (1819 total passing)**

## Performance

- **Duration:** 11 min
- **Started:** 2026-02-21T22:41:16Z
- **Completed:** 2026-02-21T22:52:30Z
- **Tasks:** 2
- **Files modified:** 8

## Accomplishments

- Created `classify_variants()` with BMV/DMV/PTV classification from EFFECT/IMPACT/SIFT/PolyPhen: 6 candidate column names for SIFT, 6 for PolyPhen, multi-value cell splitting, unclassified variants excluded from COAST only (include_mask=False)
- Created `COASTTest` class wrapping R AllelicSeries::COAST() via rpy2: parallel_safe=False, all 3 categories required, missing categories return p_value=None, omnibus p feeds ACAT-O, burden/SKAT components in extra dict
- Registered "coast" in engine._build_registry(), added "coast" to needs_regression tuple, added gene_df to gene_data dict, added coast_weights to AssociationConfig construction in stage
- CLI: --coast-weights (comma-separated floats, validated to exactly 3 values), updated --association-tests help text, coast pytest marker in pyproject.toml

## Task Commits

1. **Task 1: COASTTest with classify_variants, rpy2 wrapper, CLI** - `fc750e4` (feat)
2. **Task 2: Engine registration, ACAT-O integration, unit tests** - `910cfbd` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified

- `variantcentrifuge/association/tests/allelic_series.py` - COASTTest + classify_variants (730 lines)
- `tests/unit/test_coast.py` - 43 unit tests for COAST classification and engine integration
- `variantcentrifuge/association/tests/__init__.py` - Added COASTTest lazy loader and __all__
- `variantcentrifuge/association/base.py` - Added coast_weights field to AssociationConfig
- `variantcentrifuge/association/engine.py` - Registered COASTTest as "coast" in _build_registry()
- `variantcentrifuge/stages/analysis_stages.py` - Added "coast" to needs_regression, coast_weights to config construction, gene_df to gene_data
- `variantcentrifuge/cli.py` - Added --coast-weights CLI arg with validation, updated --association-tests help
- `pyproject.toml` - Added "coast" pytest marker

## Decisions Made

- **IMPL-47: classify_variants returns include_mask** — missense variants without SIFT/PolyPhen predictions get code 0 and include_mask=False; they remain in SKAT/burden genotype matrix but are excluded from COAST classification. This preserves SKAT/burden power while enforcing COAST's ordered-alternative assumption.
- **IMPL-48: gene_df stored in gene_data dict** — COASTTest needs per-variant annotation columns (EFFECT, IMPACT, SIFT, PolyPhen) that are in the gene DataFrame but not summarized in the standard gene_burden aggregation keys. Stored as gene_df (reset_index) in gene_data; COASTTest detects annotation/genotype mismatch via shape check.
- **IMPL-49: Annotation/genotype mismatch skips gene** — If len(anno_codes) != geno.shape[1] (can happen when site-filter removes variants that build_genotype_matrix filtered but gene_df still has), COAST returns p_value=None with ANNOTATION_GENOTYPE_MISMATCH skip reason rather than misaligning annotations.
- **IMPL-50: p_value=None when any COAST category missing** — Missing BMV, DMV, or PTV category means the ordered-alternative test cannot be applied. Returns None (skip) not a spurious p-value. Skip reason includes which categories are missing and their counts.

## Deviations from Plan

None — plan executed exactly as written, with one implementation detail (annotation/genotype mismatch detection) added as a correctness guard.

## Issues Encountered

- Minor: 9 ruff lint errors in test_coast.py (7 unused `mask` variables renamed to `_mask`, 2 nested `with` statements combined). Fixed before commit.

## User Setup Required

None — no external service configuration required. R AllelicSeries package is a runtime dependency (already documented in check_dependencies() error messages with install instructions).

## Next Phase Readiness

- COAST test fully implemented and registered; Phase 23 Plan 04 (JSON config) can add coast_weights to the JSON config schema
- All 1819 tests passing; no regressions
- Remaining Phase 23 plans: 04 (JSON config for association parameters)

---
*Phase: 23-pca-functional-weights-allelic-series-json-config*
*Completed: 2026-02-21*
