---
phase: 37-association-resource-management-memory-streaming
verified: 2026-02-25T12:05:46Z
status: passed
score: 4/4 must-haves verified
gaps: []
---

# Phase 37: Association Resource Management & Memory Streaming Verification Report

**Phase Goal:** The pipeline has a single shared ResourceManager in PipelineContext used by all stages, the GT column drop/recover antipattern is eliminated so genotype data flows cleanly through analysis stages, and genotype matrices are streamed per-gene (build-test-discard) to prevent OOM on large panels.
**Verified:** 2026-02-25T12:05:46Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | All stages use `context.resource_manager` with fallback — no standalone ResourceManager instantiations | VERIFIED | 4 locations in analysis_stages.py all use `rm = context.resource_manager; if rm is None: rm = ResourceManager(...)` pattern |
| 2 | Per-sample GT columns preserved through analysis stages — local-copy reconstruction only | VERIFIED | VariantAnalysisStage saves GT backup and re-attaches via key-column merge (lines 1545-1756); GeneBurdenAnalysisStage uses `df_for_burden` local copy (lines 1958-1963); AssociationAnalysisStage removed `context.current_dataframe = df` at old line 2424 |
| 3 | Genotype matrices built per-gene via `_GenotypeMatrixBuilder`, discarded after tests | VERIFIED | `_GenotypeMatrixBuilder` dataclass at module level (lines 43-114 of analysis_stages.py); engine invokes builder per-gene in sequential path and discards after use (engine.py lines 464-485); 25 tests pass |
| 4 | `--association-workers` defaults to `--threads` value via `_resolve_association_workers` | VERIFIED | `_resolve_association_workers` at line 2211 returns `threads` when `raw == 0`; engine uses `assoc_config.association_workers` directly from resolved config |

**Score:** 4/4 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `variantcentrifuge/pipeline_core/context.py` | `resource_manager` field on PipelineContext | VERIFIED | Line 172: `resource_manager: "ResourceManager | None" = None`; absent from `merge_from()` |
| `variantcentrifuge/pipeline.py` | One-time ResourceManager initialization | VERIFIED | Lines 529-532: `context.resource_manager = ResourceManager(config=initial_config)` after context creation, before checkpoint init |
| `variantcentrifuge/stages/analysis_stages.py` | Fallback pattern in 4 locations; `_GenotypeMatrixBuilder` at module level; local-copy GT reconstruction | VERIFIED | 4 fallback patterns confirmed; `@dataclass _GenotypeMatrixBuilder` at lines 43-114; `df_for_burden`/`_gt_backup` local-copy patterns confirmed |
| `variantcentrifuge/association/engine.py` | Builder invocation in sequential and parallel paths; matrix discard | VERIFIED | Sequential path lines 464-485; parallel path lines 422-457; `gene_data.pop("genotype_matrix")` in both |
| `variantcentrifuge/gene_burden.py` | `_find_gt_columns` is import alias, no duplicate implementation | VERIFIED | Line 31: `from .stages.output_stages import _find_per_sample_gt_columns as _find_gt_columns`; no `def _find_gt_columns` function body |
| `tests/unit/test_pipeline_context_resource_manager.py` | 8 unit tests for shared ResourceManager pattern | VERIFIED | 149 lines, 8 tests, all passing |
| `tests/unit/test_gt_lifecycle.py` | GT lifecycle tests | VERIFIED | 113 lines, 9 tests, all passing |
| `tests/unit/test_streaming_matrix.py` | Builder pickling, MAC filter, engine lifecycle tests | VERIFIED | 230 lines, 8 tests, all passing |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `pipeline.py` | `pipeline_core/context.py` | `context.resource_manager = ResourceManager(config=...)` | WIRED | Line 532 confirmed |
| `analysis_stages.py` (4 locations) | `pipeline_core/context.py` | `rm = context.resource_manager` with fallback | WIRED | Lines 1005, 3110, 3572, 3844 all use fallback pattern |
| `analysis_stages.py` (AssociationAnalysisStage) | `association/engine.py` | `gene_data["_genotype_matrix_builder"] = builder` | WIRED | Line 2726 stores builder; engine.py lines 464, 381, 427 invoke it |
| `analysis_stages.py` (_GenotypeMatrixBuilder.__call__) | `association/genotype_matrix.py` | `build_genotype_matrix(...)` | WIRED | Line 80 in `_GenotypeMatrixBuilder.__call__` |
| `analysis_stages.py` (_resolve_association_workers) | `engine.py` (n_workers) | `assoc_config.association_workers` = `threads` when unset | WIRED | `_resolve_association_workers` line 2218-2222 returns `threads`; engine line 356 uses it |

### Requirements Coverage

All four success criteria from the phase goal are satisfied:

| Requirement | Status | Notes |
|-------------|--------|-------|
| Single shared ResourceManager via `context.resource_manager` | SATISFIED | 4 fallback patterns confirmed; field in PipelineContext; initialized once in pipeline.py |
| GT columns preserved (no drop/recover cycle in context) | SATISFIED | Local-copy reconstruction in all 3 stages; key-column merge re-attachment in VariantAnalysisStage |
| Genotype matrices streamed per-gene | SATISFIED | `_GenotypeMatrixBuilder` dataclass; pop-after-use pattern in engine; 25 tests pass |
| `--threads` governs workers when `--association-workers` absent | SATISFIED | `_resolve_association_workers` correctly falls back to `threads` |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `analysis_stages.py` | 760 | `# TODO(12-02): Remove this check...` | Info | Pre-existing, unrelated to phase 37 work |

No blockers or warnings introduced by phase 37 changes.

### Human Verification Required

None — all success criteria are structurally verifiable:

- ResourceManager wiring is verified by code inspection
- GT column lifecycle is verified by pattern inspection and unit tests
- Matrix streaming is verified by `_GenotypeMatrixBuilder` structure + engine pop-after-use + passing tests
- Worker count resolution is verified by `_resolve_association_workers` implementation

### Test Results

```
25 tests in phase-37 specific test files — all passed (15.2s)
2092 fast unit tests — all passed (178s)
0 regressions
```

### Gaps Summary

No gaps. All four success criteria are fully implemented and verified against the actual codebase:

1. **ResourceManager sharing**: `PipelineContext.resource_manager` field exists, initialized once in `pipeline.py`, accessed via fallback pattern in all 4 stage locations, absent from `merge_from()`.

2. **GT column lifecycle**: `VariantAnalysisStage` saves GT columns before reconstruction and re-attaches via key-column merge after `analyze_variants`. `GeneBurdenAnalysisStage` uses `df_for_burden` local copy. `AssociationAnalysisStage` removed the old `context.variants_df` recovery fallback and the `context.current_dataframe = df` overwrite that dropped per-sample columns.

3. **Per-gene streaming**: `_GenotypeMatrixBuilder` picklable dataclass at module level encapsulates build/mask/MAC logic. Engine invokes and discards per-gene in sequential path (O(1) peak memory), builds eagerly before worker dispatch in parallel path.

4. **Worker count from `--threads`**: `_resolve_association_workers` returns `threads` value when `--association-workers` is not specified (value 0), gating all association parallelism decisions through the shared config.

---

_Verified: 2026-02-25T12:05:46Z_
_Verifier: Claude (gsd-verifier)_
