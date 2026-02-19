---
phase: 18-foundation-core-abstractions-and-fisher-refactor
plan: "02"
subsystem: analysis
tags: [association, pipeline, stage, gene-burden, fisher, correction]

# Dependency graph
requires:
  - phase: 18-01
    provides: association/ package with AssociationEngine, AssociationConfig, FisherExactTest, apply_correction
provides:
  - AssociationAnalysisStage wired into stage-based pipeline (analysis_stages.py)
  - PipelineContext.association_results field with merge_from() support
  - AssociationAnalysisStage registered in stage_registry (aliases: association_analysis, association)
  - pipeline.py conditionally includes AssociationAnalysisStage when --perform-association is set
  - gene_burden.py correction rewired to use association/correction.py (with smm fallback)
affects:
  - Phase 19 (covariate system + burden tests — adds tests to AssociationEngine registry)
  - Phase 20 (R SKAT backend — registers "skat" test in AssociationEngine)
  - Phase 21 (pure Python SKAT — registers "skat_python" test)
  - Phase 22 (ACAT-O omnibus — extends AssociationEngine with "acat_o")
  - Phase 23 (PCA + JSON config — configures association via JSON)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "AssociationAnalysisStage mirrors GeneBurdenAnalysisStage: identical aggregation strategy selection, vcf_samples pass-through, checkpoint/output patterns"
    - "Stage guards on config key (perform_association) independent from perform_gene_burden"
    - "Correction import with graceful fallback: try .association.correction, fall back to inline smm"

key-files:
  created:
    - variantcentrifuge/.planning/phases/18-foundation-core-abstractions-and-fisher-refactor/18-02-SUMMARY.md
  modified:
    - variantcentrifuge/pipeline_core/context.py
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/stages/stage_registry.py
    - variantcentrifuge/pipeline.py
    - variantcentrifuge/gene_burden.py

key-decisions:
  - "AssociationAnalysisStage reuses _aggregate_gene_burden_from_columns/_from_gt/_legacy from gene_burden.py (not duplicated) ensuring bit-identical contingency data"
  - "gene_burden.py correction uses _apply_correction from association.correction when available, falls back to inline smm — zero behavioral change"
  - "Both GeneBurdenAnalysisStage and AssociationAnalysisStage are fully independent; both can run in same pipeline invocation without interference"

patterns-established:
  - "New association tests added in later phases only need: (1) implement AssociationTest subclass, (2) register in engine._build_registry() — no stage changes needed"
  - "Association output stored in context.association_results (DataFrame) and context.config['association_output'] (TSV path)"

# Metrics
duration: 22min
completed: 2026-02-19
---

# Phase 18 Plan 02: Pipeline Wiring for AssociationAnalysisStage Summary

**AssociationAnalysisStage wired into stage-based pipeline with PipelineContext.association_results field, gene_burden.py correction rewired to association/correction.py with safe fallback**

## Performance

- **Duration:** 22 min
- **Started:** 2026-02-19T07:38:59Z
- **Completed:** 2026-02-19T08:01:12Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments

- PipelineContext extended with `association_results: pd.DataFrame | None` field and `merge_from()` support
- AssociationAnalysisStage implemented in analysis_stages.py — mirrors GeneBurdenAnalysisStage exactly (same aggregation strategy, vcf_samples pass-through, checkpoint/output patterns), guards on `perform_association` config key
- Stage registered in stage_registry with aliases `["association_analysis", "association"]`, priority 30.0
- pipeline.py conditionally includes AssociationAnalysisStage when `args.perform_association` is set
- gene_burden.py correction path rewired to use `association.correction.apply_correction` with safe fallback to inline smm — zero behavioral change

## Task Commits

Each task was committed atomically:

1. **Task 1: Add association_results to PipelineContext and create AssociationAnalysisStage** - `cc84518` (feat)
2. **Task 2: Register stage, add to pipeline.py, and rewire gene_burden.py correction** - `8cf0e14` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified

- `variantcentrifuge/pipeline_core/context.py` - Added `association_results` field; updated `merge_from()` to merge association_results
- `variantcentrifuge/stages/analysis_stages.py` - Added AssociationAnalysisStage class (170 lines); added imports from association/ package and gene_burden private functions
- `variantcentrifuge/stages/stage_registry.py` - Added AssociationAnalysisStage to _register_analysis_stages()
- `variantcentrifuge/pipeline.py` - Added AssociationAnalysisStage import and conditional stage inclusion
- `variantcentrifuge/gene_burden.py` - Added _apply_correction import from association.correction; rewired correction block with fallback

## Decisions Made

- **Reuse aggregation functions from gene_burden.py** — AssociationAnalysisStage imports `_aggregate_gene_burden_from_columns`, `_aggregate_gene_burden_from_gt`, `_aggregate_gene_burden_legacy`, and `_find_gt_columns` directly. This ensures bit-identical contingency data between --perform-gene-burden and --perform-association paths.
- **Correction rewiring is additive, not replacing** — gene_burden.py keeps the existing smm import as fallback. If association package is unavailable (ImportError), the original code path still works. This is zero-risk.
- **Both stages fully independent** — GeneBurdenAnalysisStage is completely unmodified. Both stages can run in the same pipeline invocation without interference.

## Deviations from Plan

None - plan executed exactly as written.

One minor auto-fix during execution: ruff linted an if/else block in the significance counter as preferring a ternary expression (SIM108). Fixed inline during lint pass.

## Issues Encountered

None - all imports resolved cleanly, all 858 unit tests passed without modification.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Association pipeline wiring is complete; ready for Phase 19 (covariate system + burden tests)
- Phase 19 only needs to: (1) implement burden regression AssociationTest subclass, (2) add "burden" to `_build_registry()` in engine.py — no stage changes needed
- AssociationAnalysisStage handles all test orchestration via AssociationEngine.from_names()
- Context field `association_results` is ready for downstream stages (e.g., future Excel output)

---
*Phase: 18-foundation-core-abstractions-and-fisher-refactor*
*Completed: 2026-02-19*
