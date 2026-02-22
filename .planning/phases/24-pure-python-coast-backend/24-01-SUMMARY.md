---
phase: 24
plan: 01
subsystem: association-testing
tags: [coast, allelic-series, python-backend, thread-safe, skat, burden, cauchy]
requires:
  - 21-pure-python-skat-backend
  - 22-acat-o-and-diagnostics
  - 23-pca-functional-weights-allelic-series-json-config
provides:
  - PythonCOASTBackend with test_gene() computing 7 COAST components
  - PurePythonCOASTTest wrapping PythonCOASTBackend with parallel_safe=True
decisions:
  - IMPL-51: Uniform weights [1,1,1] for baseline 3-df burden test; coast_weights for sum/max
  - IMPL-52: SKAT variance uses aaf*(1-aaf) NOT 2*aaf*(1-aaf) (matches AllelicSeries R source)
  - IMPL-53: coast_burden_p_value in extra is Cauchy of 6 burden components (not the 7-way omnibus)
  - IMPL-54: Lazy null model fits via PythonSKATBackend.fit_null_model() on first gene call
affects:
  - Future phase: engine registry swap (PurePythonCOASTTest as default coast implementation)
tech-stack:
  added: []
  patterns:
    - Pure Python COAST backend pattern (mirrors PurePythonSKATTest/PythonSKATBackend pattern)
    - Annotation-aware SKAT weights using classify_variants() codes
    - 3-layer COAST decomposition (6 burden + 1 SKAT + Cauchy omnibus)
key-files:
  created:
    - variantcentrifuge/association/backends/coast_python.py
    - variantcentrifuge/association/tests/allelic_series_python.py
  modified: []
metrics:
  duration: 4 minutes
  completed: 2026-02-22
---

# Phase 24 Plan 01: Pure Python COAST Backend Summary

**One-liner:** Pure Python COAST backend (6 burden + 1 SKAT + Cauchy with [1,1,1,1,1,1,6] weights) wrapped in PurePythonCOASTTest with parallel_safe=True, no R dependency.

## What Was Built

Two new files implementing the full COAST algorithm (McCaw et al. AJHG 2023) without R or rpy2:

### `variantcentrifuge/association/backends/coast_python.py`
**PythonCOASTBackend** — pure math layer:
- `_aggregate_by_category()`: sums genotype dosages per COAST category (BMV/DMV/PTV), supports indicator encoding and three aggregation methods (none=baseline, sum, max)
- `_run_burden_test()`: OLS (quantitative) or Logit LRT (binary) with 1-df or 3-df options
- `_compute_allelic_skat_weights()`: annotation-aware SKAT weights via `sqrt(coast_weight / aaf*(1-aaf))`
- `_run_allelic_skat()`: reuses `PythonSKATBackend._compute_eigenvalues_filtered()` + `compute_pvalue()` from davies.py
- `test_gene()`: orchestrates all 7 components and calls `cauchy_combination()` with weights `[1,1,1,1,1,1,6]`

### `variantcentrifuge/association/tests/allelic_series_python.py`
**PurePythonCOASTTest** — AssociationTest wrapper:
- `parallel_safe = True` (thread-safe, no rpy2)
- `name = "coast"` (same as COASTTest for registry swap)
- Reuses `classify_variants()` from `allelic_series.py`
- Lazy null model via `PythonSKATBackend.fit_null_model()`
- All 7 skip-condition guards from COASTTest replicated with identical `coast_skip_reason` values
- `extra` dict matches COASTTest output keys: `coast_burden_p_value`, `coast_skat_p_value`, `coast_n_bmv/dmv/ptv`

## Decisions Made

| ID | Decision | Rationale |
|----|----------|-----------|
| IMPL-51 | Baseline 3-df burden uses uniform weights [1,1,1], NOT coast_weights | Matches R AllelicSeries: baseline tests per-category effect equality without ordering assumption |
| IMPL-52 | SKAT variance = aaf*(1-aaf) NOT 2*aaf*(1-aaf) | Matches AllelicSeries R source; the factor of 2 is absorbed by Q=score'score/2 convention |
| IMPL-53 | coast_burden_p_value = Cauchy of 6 burden components separately | Preserves COASTTest's output contract; downstream consumers expect standalone burden sub-summary |
| IMPL-54 | Null model fits lazily on first gene call | Same pattern as PurePythonSKATTest; cohort-level singleton avoids repeated fitting |

## Verification Results

- Both files import without errors
- `ruff check` and `ruff format --check` pass on both files
- `PurePythonCOASTTest.parallel_safe is True`
- `PurePythonCOASTTest.name == "coast"`
- All 43 existing `tests/unit/test_coast.py` tests pass (no regressions)

## Deviations from Plan

None — plan executed exactly as written.

## Next Phase Readiness

- Phase 24 plan 02 can now add engine/registry integration (swap `"coast"` to PurePythonCOASTTest by default or on `--coast-backend python`)
- PurePythonCOASTTest is ready for end-to-end testing with synthetic data against R AllelicSeries
- No blockers; PythonCOASTBackend is fully self-contained

## Files Created

| File | Lines | Purpose |
|------|-------|---------|
| `variantcentrifuge/association/backends/coast_python.py` | 472 | PythonCOASTBackend with all 3 COAST layers |
| `variantcentrifuge/association/tests/allelic_series_python.py` | 480 | PurePythonCOASTTest AssociationTest wrapper |

## Commits

| Hash | Message |
|------|---------|
| f2c7354 | feat(24-01): implement PythonCOASTBackend pure math layer |
| 8a54817 | feat(24-01): implement PurePythonCOASTTest AssociationTest wrapper |
