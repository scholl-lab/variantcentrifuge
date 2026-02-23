# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-23)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.16.0 — Association Hardening & Multi-Cohort Features (Phase 31 COMPLETE)

## Current Position

Phase: 31 of 36 (COAST Fix — COMPLETE)
Plan: 3 of 10 (across v0.16.0)
Status: Phase 31 complete — ready for Phase 32
Last activity: 2026-02-23 — Completed 31-02-PLAN.md (COAST-02, COAST-04, COAST-05, COAST-06, COAST-07)

Progress: ███░░░░░░░ 30%

## Performance Metrics

**Velocity:**
- Total plans completed (v0.16.0): 3
- Prior milestone (v0.15.0): 35 plans, 12 phases, 5 days

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 30 — Dead Code Cleanup | 1/1 done | ~19 min | ~19 min |
| 31 — COAST Fix | 2/2 done | ~26 min | ~13 min |

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
- [31-01] COAST-03: skip threshold changed from "any category missing" to "ALL categories missing" — matches R reference (insitro/AllelicSeries drop_empty=TRUE)
- [31-01] coast_status ('complete'/'partial'/'skipped') added to TestResult.extra for diagnostic transparency
- [31-01] GT matrix fallback order: variants_df first, then current df (variants_df preferred as it preserves pre-reconstruction state)
- [31-02] COAST classification configs use COAST_* normalized column names; classify_variants() normalizes before apply_scoring()
- [31-02] coast_classification in AssociationConfig stores absolute path (None = hardcoded logic); cli.py resolves model name to path
- [31-02] Auto-injection filters out COAST_* internal names from vcf_fields (built-in models use normalized columns, not raw VCF fields)
- [31-02] diagnostics_rows parameter added to classify_variants() — ready for Phase 35 diagnostics wiring but not yet connected to file output

### Architecture Invariants

- R backend: parallel_safe=False; rpy2 calls only from main thread
- Binary traits: always SKATBinary — never continuous-trait SKAT on binary phenotypes
- FDR strategy: single pass on ACAT-O p-values across all genes (not per-test)
- Genotype matrix: never stored in PipelineContext (memory constraint)
- P-value computation: always through compute_pvalue() — never call Liu/Kuonen directly
- Weighted BH: weights MUST be renormalized to mean=1.0 at load time (Genovese 2006 FDR guarantee)
- Case-confidence weights: must be applied to null model (var_weights in GLM), not to residuals post-hoc
- COAST classification: model_dir=None → hardcoded SIFT/PolyPhen; model_dir=path → formula engine

### Blockers/Concerns

- [Phase 35] Weighted SKAT phi adjustment is MEDIUM confidence — validate with permutation before shipping
- [Phase 35] diagnostics_rows in classify_variants() is ready but not wired to file output — Phase 35 should complete this
- [Phase 36] SKAT-O "single eigendecomposition" intent is ambiguous — clarify before implementing (may already be done)
- [Phase 36] Sparse matrix breakeven threshold (500K cells, <20% density) is an estimate; profile on GCKD first

## Session Continuity

Last session: 2026-02-23T18:24:14Z
Stopped at: Completed 31-02-PLAN.md — COAST-02/04/05/06/07 (classification configs, effect resolution, CLI option)
Resume file: None
Next: Phase 32 — Region Restriction + PCA Wiring
