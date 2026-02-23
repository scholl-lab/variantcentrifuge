# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-23)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.16.0 — Association Hardening & Multi-Cohort Features (Phase 30)

## Current Position

Phase: 30 of 36 (Dead Code Cleanup)
Plan: 0 of 10 (across v0.16.0)
Status: Ready to plan
Last activity: 2026-02-23 — Roadmap created for v0.16.0

Progress: ░░░░░░░░░░ 0%

## Performance Metrics

**Velocity:**
- Total plans completed (v0.16.0): 0
- Prior milestone (v0.15.0): 35 plans, 12 phases, 5 days

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| v0.16.0 — not started | - | - | - |

*Updated after each plan completion*

## Accumulated Context

### Decisions

- [Roadmap] TD-06 maps to Phase 30 (covered by CLEAN-04; same change, not double-counted)
- [Roadmap] Phase 32 groups Region Restriction + PCA Wiring (both standalone, both small)
- [Roadmap] Phase 35 depends on Phase 31 (COAST path must work before weighted COAST runs are validated)
- [Roadmap] Phase 36 last — opt-in sparse matrices have lowest urgency at 5K-sample GCKD scale

### Architecture Invariants

- R backend: parallel_safe=False; rpy2 calls only from main thread
- Binary traits: always SKATBinary — never continuous-trait SKAT on binary phenotypes
- FDR strategy: single pass on ACAT-O p-values across all genes (not per-test)
- Genotype matrix: never stored in PipelineContext (memory constraint)
- P-value computation: always through compute_pvalue() — never call Liu/Kuonen directly
- Weighted BH: weights MUST be renormalized to mean=1.0 at load time (Genovese 2006 FDR guarantee)
- Case-confidence weights: must be applied to null model (var_weights in GLM), not to residuals post-hoc

### Blockers/Concerns

- [Phase 35] Weighted SKAT phi adjustment is MEDIUM confidence — validate with permutation before shipping
- [Phase 36] SKAT-O "single eigendecomposition" intent is ambiguous — clarify before implementing (may already be done)
- [Phase 36] Sparse matrix breakeven threshold (500K cells, <20% density) is an estimate; profile on GCKD first

## Session Continuity

Last session: 2026-02-23 — Roadmap created
Stopped at: ROADMAP.md and STATE.md written; REQUIREMENTS.md traceability updated
Resume file: None
Next: /gsd:plan-phase 30
