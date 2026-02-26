# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-26)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.17.0 — Tech Debt Cleanup & Compound Het Parallelization

## Current Position

Phase: 38 — Codebase Cleanup
Plan: 01 complete (re-executed); 02, 03 already complete
Status: In progress — Plan 38-01 dead code removal finalized
Last activity: 2026-02-26 — Completed 38-01-PLAN.md (dead code removal — stage_info.py, 11 dead functions, CLI cleanup)

Progress: ░░░░░░░░░░░░░░░░ (phase in progress)

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
| 38-01 | Public alias create_inheritance_details removed | After test cleanup no callers remain; alias was added by 38-02 as temp bridge |
| 38-01 | PATTERN_CATEGORIES/_PATTERN_CATEGORIES removed with get_pattern_category | Constant had only one consumer (the deleted function) |
| 38-02 | filter_by_inheritance_pattern extended with min_confidence kwarg | Test expected confidence threshold filtering; added as optional param |
| 38-03 | stages/__init__.py __all__ sorted alphabetically (no section comments) | ruff RUF022 requires isort-style sort; section comments break it; imports at top of file already group by module |
| 38-03 | TODO intelligent batching replaced with deferred-work NOTE | Current fixed-size batching is sufficient; explicit note sets expectation to revisit at ~100 stages |

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

(None — dead functions fully removed; all 2066 tests passing; make ci-check clean)

## Session Continuity

Last session: 2026-02-26T17:34:24Z
Stopped at: Completed 38-01-PLAN.md (dead code removal — source + tests fully cleaned)
Resume file: None
Next: Continue Phase 38 remaining plans (38-02 and 38-03 already complete by prior agents)
