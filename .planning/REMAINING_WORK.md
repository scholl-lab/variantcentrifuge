# Remaining Work — Prioritized (v0.16.0 → next)

**Date:** 2026-02-26 (verified via 6 parallel codebase audits)
**Current version:** 0.16.0
**Tests passing:** 2228
**Open issues:** 3 (#58, #60, #85)

---

## P1 — Performance

### Inheritance Pass 2: ThreadPoolExecutor still GIL-bound
`parallel_analyzer.py:235` uses `ThreadPoolExecutor`. Matrix reuse and vectorized het detection are done, but pair enumeration + parent checking remain GIL-bound.
**Options:** ProcessPoolExecutor + shared memory, or Numba @njit(parallel=True)
**Impact:** ~8-16x on large cohorts (1,452+ samples)

---

## P2 — Features (Deferred)

### #85: Case-confidence weighted status
Deferred from v0.16.0 Phase 35. Needs HPO infrastructure.

---

## P3 — Dead Code Removal

### stage_info.py — entire module dead
`variantcentrifuge/stage_info.py` (509 lines). Never imported anywhere. Stage info display utilities that were never integrated. **Delete.**

### --coast-backend r — no implementation exists
`cli.py:441` offers `"r"` as a choice for `--coast-backend`, but no `coast_r.py` file exists. Will crash at runtime if selected. **Remove `"r"` from choices.**

### Dead inheritance functions (~9 functions)
`inheritance/prioritizer.py`: `adjust_pattern_score()`, `filter_compatible_patterns()`, `get_pattern_description()`, `group_patterns_by_category()`, `resolve_conflicting_patterns()`
`inheritance/analyzer.py`: `export_inheritance_report()`, `filter_by_inheritance_pattern()`, `get_inheritance_summary()`, `process_inheritance_output()`
Never called. Verify via git blame then **delete**.

### Redundant --gzip-intermediates flag
`cli.py:938` — `store_true` with `default=True` is a no-op. Only `--no-gzip-intermediates` has effect. **Remove the redundant positive flag** or fix the pattern.

---

## P4 — Stale Documentation

### faq.md — 4 removed CLI flags referenced
`docs/source/faq.md`:
- Line 146: `--bcftools-filter` (actual flag: `--bcftools-prefilter`)
- Lines 151, 168: `--chunks` (removed in Phase 12)
- Lines 161, 179: `--interval-expansion` (doesn't exist)
- Line 184: `--no-filtering` (doesn't exist)

### association_testing.md — skat_pvalue column name
`docs/source/guides/association_testing.md` lines 334, 417, 439: says `skat_pvalue` but default output is `skat_o_pvalue`.

### changelog.md — classic pipeline reference
`docs/source/changelog.md:51`: "classic and stage-based pipeline modes" — classic pipeline is removed.

---

## P5 — Minor Cleanup

### 3 stale "refactored pipeline" docstrings
- `tests/integration/test_basic_pipeline.py:1`
- `tests/integration/test_pipeline_with_mocked_tools.py:1`
- `tests/integration/test_inheritance_analysis_integration.py:19`

### 3 stale "original pipeline" comments
- `variantcentrifuge/pipeline.py:423`
- `variantcentrifuge/stages/setup_stages.py:130`
- `variantcentrifuge/stages/processing_stages.py:1485`

### Missing __all__ exports in stages/__init__.py
`VariantAnalysisStage`, `AssociationAnalysisStage`, `ParallelCompleteProcessingStage` — used in pipeline.py but not in `__all__`.

### 2 dead integration test methods (#88 remnant)
`tests/integration/test_pipeline_with_mocked_tools.py`:
- `TestBCFToolsPrefilter` class (lines 500-763)
- `test_parallel_variant_extraction` (xfail'd)

### TODO(12-02): Legacy chunking check
`analysis_stages.py:760` — `config.get("chunks")` guard for removed CLI flag.

### TODO: Intelligent batching
`pipeline_core/runner.py:725` — "Implement intelligent batching based on estimated runtime"

### TD-05: Fisher lambda_GC doc clarification
ROADMAP criterion assumes score-based test behavior; Fisher is inherently conservative.

---

## P6 — Testing Polish

### #58: Performance CI thresholds
13 benchmarks exist, no automated regression enforcement.

### #60: Test dataset expansion
GIAB trio + GCKD golden + synthetic exist. Could expand.

---

## Deferred (Future Milestones)

| Item | Reason |
|------|--------|
| Sparse genotype matrices (Phase 36) | Streaming solves OOM; centering destroys sparsity |
| Remove R backends for SKAT/COAST | v0.17.0 |
| Single eigendecomposition for SKAT-O | ~5x speedup, medium risk |
| Kinship matrix / mixed models | Requires SAIGE-style sparse GRM |
| IHW (FDR-04) | No Python implementation exists |
