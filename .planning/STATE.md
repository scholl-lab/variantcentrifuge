# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-23)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.16.0 — Association Hardening & Multi-Cohort Features (Phase 33 in progress)

## Current Position

Phase: 33 of 36 (Gene-Level FDR Weighting — Plan 01 COMPLETE)
Plan: 6 of 10 (across v0.16.0)
Status: Phase 33 plan 01 complete — ready for next plan or phase
Last activity: 2026-02-24 — Completed 33-01-PLAN.md (weighted BH, --gene-prior-weights CLI, fdr_weight column)

Progress: ██████░░░░ 60%

## Performance Metrics

**Velocity:**
- Total plans completed (v0.16.0): 6
- Prior milestone (v0.15.0): 35 plans, 12 phases, 5 days

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 30 — Dead Code Cleanup | 1/1 done | ~19 min | ~19 min |
| 31 — COAST Fix | 2/2 done | ~26 min | ~13 min |
| 32 — Region Restriction + PCA Wiring | 2/2 done | ~44 min | ~22 min |
| 33 — Gene-Level FDR Weighting | 1/1 done | ~12 min | ~12 min |

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
- [32-01] --regions-bed CLI flag added; RegionRestrictionStage wired into GeneBedCreationStage via _intersect_with_restriction_bed
- [32-02] --pca unified flag auto-detects file path vs 'akt' tool name; --pca-file/--pca-tool are hidden deprecated aliases via dest='pca'
- [32-02] PCAComputationStage sets cfg['pca_file'] for compat with _build_assoc_config_from_context; AssociationAnalysisStage uses setdefault to not override explicit config
- [32-02] AKT cache: skip subprocess if {base_name}.pca.eigenvec already exists and is non-empty
- [33-01] Weighted BH renormalization happens at correction time (inside apply_weighted_correction), not at load time — allows reuse across different testable gene subsets
- [33-01] fdr_weight column shows NORMALIZED weights (post-renormalization), not raw file weights; column absent when --gene-prior-weights not used (backward compatible)
- [33-01] Weight loading and weighted correction stay in engine.run_all(), not AssociationAnalysisStage — engine owns all correction logic
- [33-01] write_fdr_weight_diagnostics() called from engine when both gene_prior_weights AND diagnostics_output are set
- [33-01] IHW not implemented — no flag, no stub, no error message (Python-first policy; IHW deferred to backlog)

### Architecture Invariants

- R backend: parallel_safe=False; rpy2 calls only from main thread
- Binary traits: always SKATBinary — never continuous-trait SKAT on binary phenotypes
- FDR strategy: single pass on ACAT-O p-values across all genes (not per-test)
- Genotype matrix: never stored in PipelineContext (memory constraint)
- P-value computation: always through compute_pvalue() — never call Liu/Kuonen directly
- Weighted BH: weights renormalized to mean=1.0 at correction time using testable genes subset (Genovese 2006 FDR guarantee)
- Case-confidence weights: must be applied to null model (var_weights in GLM), not to residuals post-hoc
- COAST classification: model_dir=None → hardcoded SIFT/PolyPhen; model_dir=path → formula engine
- PCA stage handoff: PCAComputationStage sets context.config['pca_file'] AND marks complete with result dict; AssociationAnalysisStage reads both

### Blockers/Concerns

- [Phase 35] Weighted SKAT phi adjustment is MEDIUM confidence — validate with permutation before shipping
- [Phase 35] diagnostics_rows in classify_variants() is ready but not wired to file output — Phase 35 should complete this
- [Phase 36] SKAT-O "single eigendecomposition" intent is ambiguous — clarify before implementing (may already be done)
- [Phase 36] Sparse matrix breakeven threshold (500K cells, <20% density) is an estimate; profile on GCKD first

## Session Continuity

Last session: 2026-02-24T08:00:50Z
Stopped at: Completed 33-01-PLAN.md — weighted BH correction, --gene-prior-weights CLI, engine wiring, fdr_weight column
Resume file: None
Next: Phase 34+ — Case Confidence / remaining v0.16.0 phases

Note: Phase 33 has only 1 plan (33-01), which is now complete.
