# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-19)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.15.0 Modular Rare Variant Association Framework

## Current Position

Phase: 19 — Covariate System + Burden Tests
Plan: 1/2 complete
Status: In progress — Plan 19-01 complete
Last activity: 2026-02-20 — Completed 19-01-PLAN.md (Covariate System and Data Infrastructure)

Progress: █████░░░░░░░░░░░░░░░░ ~20% (Phase 18 complete, Phase 19 plan 1/2 complete)

## Milestone Overview

**v0.15.0 — Modular Rare Variant Association Framework** (Phases 18-23)

| Phase | Goal | Requirements | Status |
|-------|------|--------------|--------|
| 18. Foundation | Core abstractions + Fisher refactor; bit-identical output validation | CORE-01..08 (8) | Complete |
| 19. Covariate System + Burden Tests | Logistic/linear burden tests with covariate adjustment and genotype matrix builder | COV-01..04, BURDEN-01..03, WEIGHT-01..02 (9) | Pending |
| 20. R SKAT Backend | R SKAT via rpy2 as gold standard oracle; SKATBinary + moment adjustment | SKAT-01..04, SKAT-08..09 (6) | Pending |
| 21. Pure Python SKAT Backend | Davies ctypes + saddlepoint + Liu fallback; validated against R within 10% | SKAT-05..07, SKAT-10 (4) | Pending |
| 22. ACAT-O + Diagnostics | ACAT-O omnibus; single FDR; lambda_GC; QQ TSV; sample size warnings | OMNI-01..03, DIAG-01..03, DIAG-05..06 (8) | Pending |
| 23. PCA + Functional Weights + Allelic Series + JSON Config | PCA file loading + AKT stage; CADD/REVEL weights; COAST test; JSON config; matplotlib plots | DIAG-04, PCA-01..04, SERIES-01..02, CONFIG-01..02, WEIGHT-03..05 (12) | Pending |

**Total requirements:** 47 mapped across 6 phases (8 complete, 39 pending)

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
| IMPL-01 | Clean reimplementation (not delegation) for FisherExactTest | 18-01 | fisher.py doesn't import gene_burden.py; correct coupling direction for future deprecation |
| IMPL-02 | p_value=None for zero-variant genes (not 1.0) | 18-01 | Semantically distinct: skip vs tested-with-no-signal; zero-variant genes excluded from output |
| IMPL-03 | Lazy test registry via _build_registry() | 18-01 | Defers FisherExactTest import to avoid circular imports at package load time |
| IMPL-04 | AssociationAnalysisStage reuses gene_burden.py aggregation functions | 18-02 | Ensures bit-identical contingency data between --perform-gene-burden and --perform-association paths |
| IMPL-05 | gene_burden.py correction rewired to association/correction.py with smm fallback | 18-02 | Zero behavioral change; if association package unavailable, original code path still works |
| IMPL-06 | GeneBurdenAnalysisStage and AssociationAnalysisStage are fully independent | 18-02 | Both guard on separate config keys; both can run in same pipeline invocation without interference |
| IMPL-07 | parser.error() for --association-tests without --perform-association | 18-03 | argparse convention; produces correctly formatted usage message |
| IMPL-08 | Association sheet mirrors Gene Burden sheet pattern verbatim | 18-03 | Explicit duplication preferred over abstraction for parallel maintainability |
| TEST-01 | Bit-identity uses == (exact equality) not pytest.approx for Fisher p-values | 18-04 | Same scipy call chain guarantees floating-point reproducibility; tolerance would hide regressions |
| TEST-02 | CORE-05 verified via source inspection (inspect.getsource) | 18-04 | Structural proof that GeneBurdenAnalysisStage._process() never references perform_association key |
| IMPL-09 | parse_gt_to_dosage returns (int\|None, bool) not int\|None | 19-01 | Multi-allelic flag needed to emit 'run bcftools norm' warning without second parse pass |
| IMPL-10 | load_covariates returns (np.ndarray, list[str]) tuple | 19-01 | Column names returned alongside matrix for diagnostics; callers can ignore second element |
| IMPL-11 | build_genotype_matrix: sample_mask is list[bool], all samples remain in geno | 19-01 | Callers (logistic burden test) decide whether to exclude high-missing samples |

### Architecture Invariants (from research)

- R backend: parallel_safe=False; rpy2 calls only from main thread (segfault risk otherwise)
- Binary traits: always SKATBinary — never continuous-trait SKAT on binary phenotypes
- Davies defaults: acc=1e-9, lim=10^6 (corrected from Liu et al. 2016)
- Covariate alignment: always reindex to vcf_samples order; assert no NaN after reindex
- FDR strategy: single pass on ACAT-O p-values across all genes (not per-test)
- Genotype matrix: never stored in PipelineContext (5K samples x 50K variants = 1.6 GB)
- Eigenvalue stability: scipy.linalg.eigh; threshold max(eigenvalues, 0); skip if matrix_rank < 2
- Python version: recommend bumping requires-python to >=3.11 (scipy 1.16 dropped 3.10)
- apply_correction([]) returns empty array (statsmodels multipletests raises ZeroDivisionError on empty; guarded in correction.py)

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See archived REQUIREMENTS.md Future Requirements.
- **RESEARCH-01** (before Phase 20): Validate whether parallel_safe=False on the stage is sufficient for rpy2 thread safety, or whether a dedicated R worker process with queue is required.
- **RESEARCH-02** (before Phase 21): Saddlepoint approximation algorithm for middle tier of Davies fallback chain — not in scipy stdlib; need reference from fastSKAT (PMC4375394) before implementation.

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-20
Stopped at: Completed 19-01-PLAN.md (Covariate System and Data Infrastructure)
Resume file: None
Next: Execute Plan 19-02 (LogisticBurdenTest + LinearBurdenTest)
