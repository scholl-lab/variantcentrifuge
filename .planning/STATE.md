# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-27)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** Planning next milestone

## Current Position

Phase: None — v0.17.0 complete
Plan: N/A
Status: Ready for next milestone
Last activity: 2026-02-27 — v0.17.0 milestone archived

Next: New milestone via /gsd:new-milestone

## Performance Metrics

**Milestone velocity (v0.17.0):**
- 2 phases, 5 plans, 2 days (2026-02-26 → 2026-02-27)
- 45 files changed, 3,214 insertions, 1,580 deletions
- 14/16 requirements shipped (2 intentionally abandoned)

## Accumulated Context

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

Last session: 2026-02-27
Stopped at: v0.17.0 milestone archived
Resume file: None
Next: New milestone via /gsd:new-milestone
