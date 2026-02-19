# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-19)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.15.0 Modular Rare Variant Association Framework

## Current Position

Phase: 18 — Foundation: Core Abstractions and Fisher Refactor
Plan: —
Status: Ready to plan Phase 18
Last activity: 2026-02-19 — Roadmap created for v0.15.0

Progress: ░░░░░░░░░░░░░░░░░░░░░ 0%

## Milestone Overview

**v0.15.0 — Modular Rare Variant Association Framework** (Phases 18-23)

| Phase | Goal | Requirements | Status |
|-------|------|--------------|--------|
| 18. Foundation | Core abstractions + Fisher refactor; bit-identical output validation | CORE-01..08 (8) | Pending |
| 19. Covariate System + Burden Tests | Logistic/linear burden tests with covariate adjustment and genotype matrix builder | COV-01..04, BURDEN-01..03, WEIGHT-01..02 (9) | Pending |
| 20. R SKAT Backend | R SKAT via rpy2 as gold standard oracle; SKATBinary + moment adjustment | SKAT-01..04, SKAT-08..09 (6) | Pending |
| 21. Pure Python SKAT Backend | Davies ctypes + saddlepoint + Liu fallback; validated against R within 10% | SKAT-05..07, SKAT-10 (4) | Pending |
| 22. ACAT-O + Diagnostics | ACAT-O omnibus; single FDR; lambda_GC; QQ TSV; sample size warnings | OMNI-01..03, DIAG-01..03, DIAG-05..06 (8) | Pending |
| 23. PCA + Functional Weights + Allelic Series + JSON Config | PCA file loading + AKT stage; CADD/REVEL weights; COAST test; JSON config; matplotlib plots | DIAG-04, PCA-01..04, SERIES-01..02, CONFIG-01..02, WEIGHT-03..05 (12) | Pending |

**Total requirements:** 47 mapped across 6 phases

## Accumulated Context

### Decisions

| ID | Decision | Plan | Rationale |
|----|----------|------|-----------|
| MILE-01 | Dual SKAT backend (R via rpy2 + pure Python) | Milestone | R is gold standard for validation; Python for portability |
| MILE-02 | Compiled Davies method via ctypes (qfc.c) | Milestone | Exact p-values without R dependency; Liu fallback for safety |
| MILE-03 | Support both AKT and PLINK for PCA | Milestone | AKT for VCF-native workflows, PLINK for pre-existing BED files |
| MILE-04 | Full scope (Steps 1-7 from design doc) | Milestone | Include allelic series, functional weights, quantitative traits, JSON config |
| ARCH-01 | R backend declares parallel_safe=False | Phase 20 | rpy2 is not thread-safe; ThreadPoolExecutor causes segfaults with no traceback |
| ARCH-02 | Binary trait gating locked at null model construction (Phase 18) | Phase 18 | SKATBinary vs SKAT decision cannot be retrofitted later; affects all downstream phases |
| ARCH-03 | Single FDR on ACAT-O across genes | Phase 22 | Applying FDR separately per test is statistically incorrect |
| ARCH-04 | DIAG-04 (matplotlib plots) in Phase 23, not Phase 22 | Roadmap | Optional and architecturally independent; all other diagnostics must work first |

### Architecture Invariants (from research)

- R backend: parallel_safe=False; rpy2 calls only from main thread (segfault risk otherwise)
- Binary traits: always SKATBinary — never continuous-trait SKAT on binary phenotypes
- Davies defaults: acc=1e-9, lim=10^6 (corrected from Liu et al. 2016)
- Covariate alignment: always reindex to vcf_samples order; assert no NaN after reindex
- FDR strategy: single pass on ACAT-O p-values across all genes (not per-test)
- Genotype matrix: never stored in PipelineContext (5K samples x 50K variants = 1.6 GB)
- Eigenvalue stability: scipy.linalg.eigh; threshold max(eigenvalues, 0); skip if matrix_rank < 2
- Python version: recommend bumping requires-python to >=3.11 (scipy 1.16 dropped 3.10)

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See archived REQUIREMENTS.md Future Requirements.
- **RESEARCH-01** (before Phase 20): Validate whether parallel_safe=False on the stage is sufficient for rpy2 thread safety, or whether a dedicated R worker process with queue is required.
- **RESEARCH-02** (before Phase 21): Saddlepoint approximation algorithm for middle tier of Davies fallback chain — not in scipy stdlib; need reference from fastSKAT (PMC4375394) before implementation.

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-19
Stopped at: Roadmap created — ready to plan Phase 18
Resume file: None
Next: `/gsd:plan-phase 18`
