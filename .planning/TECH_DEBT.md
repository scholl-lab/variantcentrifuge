# Tech Debt Backlog

Updated 2026-02-26 (full codebase audit via 6 parallel agents).

## Resolved in v0.16.0

- ~~TD-01: PCAComputationStage not wired~~ → Phase 32
- ~~TD-02: create_stages_from_config() missing association mapping~~ → Phase 34
- ~~TD-03: COAST R golden values~~ → Phase 34
- ~~TD-04: skat_o_p_value naming~~ → Phase 34
- ~~TD-06: "refactored_pipeline" strings~~ → Phase 30

## Open — Dead Code

### stage_info.py — entire dead module (509 lines)
`variantcentrifuge/stage_info.py` — never imported. Delete.

### --coast-backend r — no implementation
`cli.py:441` offers `"r"` choice but `coast_r.py` doesn't exist. Runtime crash. Remove from choices.

### ~9 dead inheritance functions
`prioritizer.py`: adjust_pattern_score, filter_compatible_patterns, get_pattern_description, group_patterns_by_category, resolve_conflicting_patterns
`analyzer.py`: export_inheritance_report, filter_by_inheritance_pattern, get_inheritance_summary, process_inheritance_output

### Redundant --gzip-intermediates flag
`cli.py:938` — store_true with default=True is no-op. Only --no-gzip-intermediates has effect.

## Open — Stale Docs

### faq.md references 4 removed CLI flags
`docs/source/faq.md`: --bcftools-filter (→ --bcftools-prefilter), --chunks (removed), --interval-expansion (doesn't exist), --no-filtering (doesn't exist)

### association_testing.md column name
Lines 334, 417, 439: skat_pvalue → skat_o_pvalue

### changelog.md classic pipeline reference
Line 51: "classic and stage-based pipeline modes" — classic removed.

## Open — Minor

### 3 "refactored pipeline" docstrings in integration tests
test_basic_pipeline.py:1, test_pipeline_with_mocked_tools.py:1, test_inheritance_analysis_integration.py:19

### 3 "original pipeline" comments in source
pipeline.py:423, setup_stages.py:130, processing_stages.py:1485

### Missing __all__ exports
stages/__init__.py: VariantAnalysisStage, AssociationAnalysisStage, ParallelCompleteProcessingStage

### 2 dead integration test methods
test_pipeline_with_mocked_tools.py: TestBCFToolsPrefilter (500-763), test_parallel_variant_extraction (490-497)

### TODO(12-02) legacy chunking check
analysis_stages.py:760

### TODO intelligent batching
runner.py:725

### TD-05: Fisher lambda_GC doc clarification
