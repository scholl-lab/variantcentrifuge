# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-26)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.17.0 — Tech Debt Cleanup & Compound Het Parallelization

## Current Position

Phase: 39 — Compound Het Parallelization — In Progress
Plan: 1/2 complete
Status: Plan 39-01 complete
Last activity: 2026-02-26 — Completed 39-01-PLAN.md (benchmark + baseline)

Progress: ██████████░░░░░░ 56% (1.5/2 phases complete, 4/5 plans done)

Next: Execute Phase 39 Plan 02 (optimization)

## Performance Metrics

**Milestone velocity (v0.16.0):**
- 6 active phases, 12 plans, 3 days (2026-02-23 → 2026-02-26)
- 62 files changed, 6,372 insertions, 2,251 deletions
- 227 new tests (2,001 → 2,228)

## Accumulated Context

### Decisions

| Plan  | Decision | Rationale |
|-------|----------|-----------|
| 38-01 | create_inheritance_details renamed to _create_inheritance_details | Internal helper with no external callers; underscore prefix signals private |
| 38-01 | PATTERN_CATEGORIES removed with get_pattern_category | Constant had only one consumer (the deleted function) |
| 38-03 | stages/__init__.py __all__ sorted alphabetically | ruff RUF022 requires isort-style sort |
| 38-03 | TODO intelligent batching replaced with deferred-work NOTE | Current fixed-size batching sufficient; revisit at ~100 stages |
| 39-01 | Benchmark measures total wall-time (all 3 passes), not isolated Pass 2 | Pass 1+3 fast relative to Pass 2 for large datasets; total time is reasonable proxy |
| 39-01 | Baseline confirms GIL contention: parallel 0.21x-0.92x vs sequential pre-optimization | Empirical evidence supports CONTEXT.md hypothesis; optimization target: parallel >= sequential |

### Architecture Invariants

- ResourceManager: initialized once in pipeline.py, shared via context.resource_manager; fallback pattern in stages for test compatibility
- R backend: parallel_safe=False; rpy2 calls only from main thread
- Binary traits: always SKATBinary — never continuous-trait SKAT on binary phenotypes
- FDR strategy: single pass on ACAT-O p-values across all genes (not per-test)
- Genotype matrix: streamed per-gene via _GenotypeMatrixBuilder (build-test-discard); peak memory O(1 gene)
- P-value computation: always through compute_pvalue() — never call Liu/Kuonen directly
- Weighted BH: weights renormalized to mean=1.0 at correction time using testable genes subset
- COAST classification: model_dir=None → hardcoded SIFT/PolyPhen; model_dir=path → formula engine
- GT Lifecycle Invariant: per-sample GEN_N__GT columns survive through all analysis stages; reconstruct_gt_column only called on local copies

### Blockers/Concerns

(None)

## Session Continuity

Last session: 2026-02-26 22:02Z
Stopped at: Completed 39-01-PLAN.md (benchmark created, baseline captured)
Resume file: None
Next: Execute Phase 39 Plan 02 (GIL contention optimization)
