---
phase: 23-pca-functional-weights-allelic-series-json-config
plan: 02
subsystem: association
tags: [cadd, revel, functional-weights, variant-annotation, burden-test, numpy, scipy]

# Dependency graph
requires:
  - phase: 23-01
    provides: PCA integration, extended AssociationConfig, genotype matrix loop in stage
  - phase: 19
    provides: get_weights(), beta_maf_weights(), uniform_weights(), AssociationConfig.variant_weights

provides:
  - cadd_weights(): Beta(MAF) x min(CADD_phred/cap, 1.0) with NaN-safe score parsing
  - revel_weights(): Beta(MAF) x REVEL_score (no normalization, REVEL is [0,1])
  - combined_weights(): CADD preferred, fallback to REVEL or Beta(MAF)-only
  - get_weights() extended with keyword-only cadd_scores/revel_scores/variant_effects/weight_params
  - LOF_EFFECTS / MISSENSE_EFFECTS frozenset constants for type-aware fallback logging
  - AssociationConfig.variant_weight_params field
  - CLI --variant-weight-params JSON arg + --variant-weights updated help text
  - Stage annotation extraction: CADD/REVEL/EFFECT per-gene in genotype matrix loop
  - logistic_burden.py and linear_burden.py wired to pass functional scores to get_weights()
  - 45 unit tests in test_functional_weights.py

affects:
  - 23-03 (allelic series — COAST test may need variant effect annotations)
  - 23-04 (JSON config — weight_params configurable via JSON)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "NaN-safe score parsing via _parse_scores_to_float() handles None/'.'/NA/nan"
    - "Keyword-only extension of get_weights() preserves backward compatibility for all existing callers"
    - "Type-aware fallback logging groups missing scores by LoF/missense/other"
    - "Site-filter mask replication in stage for annotation alignment with build_genotype_matrix output"

key-files:
  created:
    - tests/unit/test_functional_weights.py
  modified:
    - variantcentrifuge/association/weights.py
    - variantcentrifuge/association/base.py
    - variantcentrifuge/cli.py
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/association/tests/logistic_burden.py
    - variantcentrifuge/association/tests/linear_burden.py

key-decisions:
  - "Keyword-only kwargs on get_weights() for backward compatibility: all existing callers with 2 positional args unchanged"
  - "LOF_EFFECTS/MISSENSE_EFFECTS as module-level frozensets for O(1) membership testing"
  - "CADD cap default=40.0 matches dbNSFP convention; configurable via weight_params dict"
  - "combined_weights prefers CADD over REVEL when both provided"
  - "Stage annotation extraction replicates build_genotype_matrix site-filter logic inline to align n_kept_variants"
  - "Annotation extraction only runs when variant_weights in ('cadd', 'revel', 'combined') — zero overhead for other specs"

patterns-established:
  - "Functional weight functions follow same structure: maf_w = beta_maf_weights(); nan_mask from parsed scores; functional = np.where(nan_mask, 1.0, score_normalized)"
  - "Missing scores always default to functional=1.0 (conservative: no penalty, no up-weighting)"

# Metrics
duration: 10min
completed: 2026-02-21
---

# Phase 23 Plan 02: Functional Variant Weights Summary

**CADD/REVEL annotation-based variant weights with type-aware missing-score fallback, wired end-to-end from CLI through stage to burden tests**

## Performance

- **Duration:** ~10 min
- **Started:** 2026-02-21T22:26:48Z
- **Completed:** 2026-02-21T22:36:41Z
- **Tasks:** 2/2
- **Files modified:** 6 files modified, 1 created

## Accomplishments

- Extended `weights.py` with `cadd_weights`, `revel_weights`, `combined_weights` functions and updated `get_weights()` with keyword-only args for full backward compatibility
- Wired annotation extraction (CADD/REVEL/EFFECT columns) per-gene in the stage's genotype matrix loop, with correct alignment to `variant_mafs` via replicated site-filter mask
- Updated `logistic_burden.py` and `linear_burden.py` to pass functional scores from `contingency_data` to `get_weights()`
- Added 45 unit tests covering all weight schemes, type-aware fallback logging, numerical difference from uniform, dispatch, and regression safety

## Task Commits

Each task was committed atomically:

1. **Task 1: Add functional weight functions to weights.py and extend CLI** - `6a9ebeb` (feat)
2. **Task 2: Stage integration (annotation extraction + weight wiring) and unit tests** - `f8bf773` (feat)

## Files Created/Modified

- `variantcentrifuge/association/weights.py` — Added cadd_weights, revel_weights, combined_weights, _parse_scores_to_float, _log_missing_score_counts, LOF_EFFECTS, MISSENSE_EFFECTS; extended get_weights() with keyword-only optional args
- `variantcentrifuge/association/base.py` — Added variant_weight_params: dict | None to AssociationConfig
- `variantcentrifuge/cli.py` — Updated --variant-weights help text; added --variant-weight-params JSON arg with validation; wired variant_weight_params to context.config
- `variantcentrifuge/stages/analysis_stages.py` — Wired variant_weight_params into AssociationConfig constructor; added CADD/REVEL/EFFECT annotation extraction with site-filter alignment in genotype matrix loop
- `variantcentrifuge/association/tests/logistic_burden.py` — Updated get_weights() call to pass cadd_scores, revel_scores, variant_effects, weight_params from contingency_data
- `variantcentrifuge/association/tests/linear_burden.py` — Same update as logistic_burden.py
- `tests/unit/test_functional_weights.py` — 45 unit tests for functional weight schemes (created)

## Decisions Made

- **IMPL-43: Keyword-only kwargs on get_weights()** — Added `cadd_scores`, `revel_scores`, `variant_effects`, `weight_params` as keyword-only args with None defaults. All existing callers `get_weights(mafs, spec)` continue to work unchanged; new kwargs are ignored for "beta:*" and "uniform" specs.
- **IMPL-44: Site-filter mask replicated in stage for annotation alignment** — `build_genotype_matrix` internally applies a `keep_variants_mask` that filters out high-missing variants. To align annotation arrays with the returned `mafs` (n_kept_variants), the stage replicates the same logic using `parse_gt_to_dosage` when `len(mafs) < len(gene_df)`. In the common case (no filtering), this is skipped entirely.
- **IMPL-45: combined_weights prefers CADD over REVEL** — When both scores are available, CADD is used (higher CADD phred scores = more damaging). REVEL is used only when CADD is absent.
- **IMPL-46: Annotation extraction conditional on weight spec** — Extraction block only executes when `variant_weights in ("cadd", "revel", "combined")`; zero overhead for other specs.

## Deviations from Plan

None - plan executed exactly as written.

The site-filter alignment issue (annotation array length mismatch with variant_mafs) was anticipated: the plan notes "annotation arrays have length = n_variants (same as variant_mafs)". The implementation correctly replicates `build_genotype_matrix`'s `keep_variants_mask` logic to achieve this.

## Issues Encountered

None — implementation proceeded without unexpected problems.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- Functional weight schemes (CADD, REVEL, combined) are ready to use with `--variant-weights cadd|revel|combined` and `--variant-weight-params '{"cadd_cap": 30}'`
- CADD/REVEL annotation columns extracted per-gene and passed into burden test runs
- 1776 tests passing (45 added from this plan, no regressions)
- Phase 23 Plan 03 (allelic series — COAST test) can proceed; variant_effects arrays are now extracted and available in gene_data

---
*Phase: 23-pca-functional-weights-allelic-series-json-config*
*Completed: 2026-02-21*
