---
phase: 30
plan: "01"
subsystem: stages
tags: [dead-code, refactor, cleanup, stages, naming]

dependency-graph:
  requires: []
  provides:
    - clean-stage-source-files
    - clean-stage-registry
    - normalized-pipeline-naming
  affects:
    - phase-31-coast-fix
    - phase-32-region-restriction
    - all-subsequent-v0.16.0-phases

tech-stack:
  added: []
  patterns:
    - removed-dead-stage-classes

key-files:
  created: []
  modified:
    - variantcentrifuge/stages/processing_stages.py
    - variantcentrifuge/stages/analysis_stages.py
    - variantcentrifuge/stages/output_stages.py
    - variantcentrifuge/stages/stage_registry.py
    - variantcentrifuge/stages/__init__.py
    - variantcentrifuge/pipeline.py
    - variantcentrifuge/pipeline_core/runner.py
    - variantcentrifuge/pipeline_core/error_handling.py
    - variantcentrifuge/analyze_variants.py
    - tests/unit/stages/test_processing_stages_critical.py
    - tests/unit/stages/test_processing_stages.py
    - tests/unit/test_parallel_stage_resume.py
    - tests/unit/test_pca.py
    - tests/unit/stages/test_analysis_stages.py
    - tests/unit/stages/test_output_stages_simple.py
    - tests/unit/stages/test_parallel_complete_processing.py

decisions:
  - id: keep-pca-assoc-tests
    description: >
      test_pca.py tests both association/pca.py (live module) and PCAComputationStage
      (dead stage). Only the TestPcaComputationStage class was removed; the rest of the
      file was kept intact.
  - id: fix-test-assertions
    description: >
      Two tests in test_parallel_stage_resume.py and test_parallel_complete_processing.py
      asserted that context.mark_complete("parallel_variant_extraction") was called. After
      removing that call from ParallelCompleteProcessingStage._process and
      _handle_checkpoint_skip, these test assertions were updated to match current behavior.

metrics:
  duration: "~19 minutes"
  completed: "2026-02-23"
---

# Phase 30 Plan 01: Dead Stage Class Removal Summary

**One-liner:** Removed 8 dead stage classes plus all "refactored_pipeline"/"old pipeline"/"Refactored" naming artifacts.

## What Was Done

### Task 1: Remove 8 Dead Stage Classes

Removed 5 dead classes from `processing_stages.py`:
- `ParallelVariantExtractionStage` (parallel VCF extraction, never used)
- `BCFToolsPrefilterStage` (standalone bcftools filter, superseded by VariantExtractionStage)
- `TranscriptFilterStage` (transcript-level filtering, never wired in pipeline)
- `StreamingDataProcessingStage` (streaming processing alternative, never enabled)
- `PCAComputationStage` (AKT PCA wrapper, moved to Phase 32 scope)

Removed 2 dead classes from `analysis_stages.py`:
- `GenotypeFilterStage` (genotype-mode filtering, never called by pipeline)
- `ParallelAnalysisOrchestrator` (gene-parallel analysis, never enabled)

Removed 1 dead class from `output_stages.py`:
- `ParallelReportGenerationStage` (parallel report generation, never used in pipeline.py)

**Cross-reference cleanup:**
- Removed `"parallel_variant_extraction"` from `SnpSiftFilterStage.soft_dependencies`
- Removed `"transcript_filtering"` from `FieldExtractionStage.soft_dependencies`
- Removed `"genotype_filtering"` from `VariantAnalysisStage.soft_dependencies`
- Removed two `context.mark_complete("parallel_variant_extraction")` calls from
  `ParallelCompleteProcessingStage`

**Registry cleanup (`stage_registry.py`):**
- Removed all 8 dead stage imports and `register_stage()` calls from
  `_register_processing_stages()`, `_register_analysis_stages()`, `_register_output_stages()`

**`__init__.py` cleanup:**
- Removed imports: `ParallelAnalysisOrchestrator`, `ParallelReportGenerationStage`,
  `BCFToolsPrefilterStage`, `ParallelVariantExtractionStage`, `StreamingDataProcessingStage`
- Removed corresponding `__all__` entries

**Test cleanup:**
- `test_processing_stages_critical.py`: removed `TestBCFToolsPrefilterStage`,
  `TestStreamingDataProcessingStage`
- `test_processing_stages.py`: removed `TestParallelVariantExtractionStage`
- `test_parallel_stage_resume.py`: removed `TestParallelVariantExtractionStage`;
  fixed `TestParallelCompleteProcessingStage` assertion
- `test_pca.py`: removed `TestPcaComputationStage` class only (rest of file kept)
- `test_analysis_stages.py`: removed `TestParallelAnalysisOrchestrator`
- `test_output_stages_simple.py`: removed `TestParallelReportGenerationStage`
- `test_parallel_complete_processing.py`: fixed assertion for removed `mark_complete` call

### Task 2: Normalize Naming Artifacts

**CLEAN-04/TD-06 — "refactored_pipeline" defaults (2 occurrences):**
- `pipeline.py`: `"refactored_pipeline"` → `"pipeline"` in `initial_config.get(...)`
- `pipeline_core/runner.py`: `"refactored_pipeline"` → `"pipeline"` in `context.config.get(...)`

**CLEAN-06 — "Refactored" docstrings (2 occurrences):**
- `pipeline_core/error_handling.py`: `"Enhanced error handling utilities for the refactored pipeline."` → `"Enhanced error handling utilities for the pipeline."`
- `analyze_variants.py`: `"Refactored to be more modular:"` → `"Modular design:"`

**CLEAN-05 — "old pipeline" comments (9 occurrences in stage files):**
- `output_stages.py` (7): removed "like the old pipeline" / "to match old pipeline" phrasing
- `analysis_stages.py` (1): removed "like the old pipeline does" from stats path comment
- `processing_stages.py` (2): reworded docstring and sort comment

## Verification Results

All requirements met:

1. `grep -rn "ParallelVariantExtractionStage|BCFToolsPrefilterStage|..."` — 0 matches
2. `grep -rn "refactored_pipeline|old pipeline"` — 0 matches
3. `grep -rn "Refactored" error_handling.py analyze_variants.py` — 0 matches
4. `grep -n "refactored_args" cli.py` — still present (local variable, preserved)
5. `make test-fast` — 1967 passed, 0 failures
6. `make ci-check` — ALL CI CHECKS PASSED

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed test assertions for parallel_variant_extraction mark_complete**

- **Found during:** Task 1 — running `make test-fast` after Task 1
- **Issue:** Two tests in `test_parallel_stage_resume.py` and
  `test_parallel_complete_processing.py` asserted that
  `context.mark_complete("parallel_variant_extraction")` was called.
  After removing that call from `ParallelCompleteProcessingStage`, these assertions
  failed.
- **Fix:** Updated both test files to remove the `parallel_variant_extraction` assertion.
  The behavior is now correct — `ParallelCompleteProcessingStage` marks `variant_extraction`,
  `snpsift_filtering`, and `field_extraction` complete, not the removed dead stage.
- **Files modified:** `tests/unit/test_parallel_stage_resume.py`,
  `tests/unit/stages/test_parallel_complete_processing.py`

**2. [Rule 2 - Missing Critical] Removed unused `filter_final_tsv_by_genotype` import**

- **Found during:** Task 1 — ruff detected unused import after `GenotypeFilterStage` removal
- **Issue:** `from ..filters import filter_final_tsv_by_genotype` was only used by
  `GenotypeFilterStage`. After class removal, the import became unused.
- **Fix:** Removed the import. ruff handled this automatically via `make format`.
- **Files modified:** `variantcentrifuge/stages/analysis_stages.py`

## Next Phase Readiness

Phase 30 continues with plans 02-04 (further cleanup tasks as defined in 30-RESEARCH.md).
All subsequent v0.16.0 phases benefit from this cleanup — the codebase is now free of
dead stage noise.
