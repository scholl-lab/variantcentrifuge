# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-23)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.16.0 — Association Hardening & Multi-Cohort Features (Phase 30)

## Current Position

Phase: 30 of 36 (Dead Code Cleanup)
Plan: 1 of 10 (across v0.16.0)
Status: In progress
Last activity: 2026-02-23 — Completed 30-01-PLAN.md (dead stage class removal)

Progress: █░░░░░░░░░ 10%

## Performance Metrics

**Velocity:**
- Total plans completed (v0.16.0): 1
- Prior milestone (v0.15.0): 35 plans, 12 phases, 5 days

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 30 — Dead Code Cleanup | 1 done | 4 plans est. | ~19 min |

*Updated after each plan completion*

## Accumulated Context

### Decisions

- [Roadmap] TD-06 maps to Phase 30 (covered by CLEAN-04; same change, not double-counted)
- [Roadmap] Phase 32 groups Region Restriction + PCA Wiring (both standalone, both small)
- [Roadmap] Phase 35 depends on Phase 31 (COAST path must work before weighted COAST runs are validated)
- [Roadmap] Phase 36 last — opt-in sparse matrices have lowest urgency at 5K-sample GCKD scale
- [30-01] PCAComputationStage removed from processing_stages.py; PCA wiring belongs in Phase 32
- [30-01] test_pca.py kept (tests association/pca.py module); only TestPcaComputationStage class removed
- [30-01] ParallelCompleteProcessingStage.mark_complete no longer marks "parallel_variant_extraction" — tests updated

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

Last session: 2026-02-23
Stopped at: Completed 30-01-PLAN.md (dead stage class removal + naming normalization)
Resume file: None
Next: Continue Phase 30 plans (30-02 through 30-04)
