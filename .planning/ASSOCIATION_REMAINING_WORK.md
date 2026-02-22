# Association Framework — Remaining Work & Polish Plan

**Date:** 2026-02-22 (updated)
**Scope:** Post-Phase 24 work to finalize v0.15.0 milestone
**Consolidates:** ASSOCIATION-FRAMEWORK-DESIGN.md (open items), SKAT_OPTIMIZATION_PLAN.md (all 7 opts)

---

## Current State

Phases 18-26 are complete. Phase 27 is in progress (2/3 plans done).

The association framework is **functionally complete and documented**:
- 6 tests: Fisher, logistic burden, linear burden, SKAT, SKAT-O, COAST
- ACAT-O omnibus combination + ACAT-V per-variant score test
- Pure Python backends default for both SKAT and COAST (R deprecated)
- Saddlepoint-before-Liu fallback for extreme tail accuracy
- Covariates, PCA file loading, Beta/CADD/REVEL weights
- Diagnostics: lambda_GC, QQ data, QQ plot (matplotlib)
- FDR/Bonferroni correction
- JSON config with `"association"` section and 27 validated keys
- Comprehensive user guide, API stubs, FAQ, changelog, cross-references
- 2000+ tests passing, CI clean

---

## Part 1: Make Python Default, R Optional — COMPLETE ✅

Done in Phase 25.

- ✅ 1.1 Swap default backend to Python (`skat_backend="python"`, `coast_backend="python"`)
- ✅ 1.2 Mark R backends deprecated (DeprecationWarning on use)
- 🔮 1.3 Remove R backends — deferred to v0.17.0

---

## Part 2: Documentation — COMPLETE ✅

Done in Phase 26. Fact-checked and corrected in post-phase review.

- ✅ 2.1 Association testing guide (`docs/source/guides/association_testing.md`, ~850 lines)
- ✅ 2.2 Updated existing docs (usage.md, cohort_analysis.md, faq.md, index.md, README.md)
- ✅ 2.3 API reference stubs (`docs/source/api/association.md`)
- ✅ 2.4 Changelog v0.15.0 entry
- ✅ Post-review: 10 factual errors corrected, 2 stale API stubs removed (replacer.md, phenotype_filter.md)

---

## Part 3: Performance Optimizations — IN PROGRESS (Phase 27)

### 3.1 Saddlepoint Before Liu Fallback — COMPLETE ✅

Done in Phase 25.

### 3.2 Gene-Level Parallelization — IN PROGRESS

Phase 27 plans 01 + 02 complete:
- ✅ 27-01: Gauss-Legendre quadrature for SKAT-O integration + `parallel_safe=True` on all Python-backend tests
- ✅ 27-02: `--association-workers` CLI arg, `AssociationConfig.association_workers` field, stage builder plumbing, unit tests

Remaining:
- ⬜ 27-03: Actual `ProcessPoolExecutor` implementation in `engine.py`, BLAS thread pinning, worker initializer, integration tests

### 3.3 Cache/Interpolate Davies in Omnibus Integration — SUPERSEDED

Replaced by Gauss-Legendre quadrature (27-01), which eliminates adaptive overhead entirely.

### 3.4 Single Eigendecomposition for SKAT-O — NOT STARTED

**Impact:** ~5x speedup on SKAT-O eigenvalue step (7 eigh calls → 1)
**Effort:** ~80 lines in `python_backend.py`
**Risk:** Medium — requires careful linear algebra
**Status:** Not in current Phase 27 plans. Candidate for post-v0.15.0.

### 3.5 Sparse Genotype Matrices — NOT STARTED

**Impact:** Memory and speed for n > 10K samples with rare variants (MAF < 1%)
**Status:** Not in current Phase 27 plans. Candidate for post-v0.15.0.

### 3.6 ACAT-V Per-Variant Score Test — COMPLETE ✅

Done in Phase 25.

---

## Part 4: Minor Feature Gaps — MOSTLY RESOLVED

### 4.1 JSON Config for Association — COMPLETE ✅

Done in Phase 23. `"association"` section in config.json with 27 validated keys,
CLI override precedence, type/enum validation.

### 4.2 AKT/PLINK Pipeline Stage Wrappers — DEFERRED

Users run AKT/PLINK externally and provide eigenvec file. Defer to post-v0.15.0.

### 4.3 Kinship Matrix / Mixed Models — DEFERRED

Not in scope for v0.15.0. Would require SAIGE-style sparse GRM support.

---

## Remaining Work for v0.15.0

| # | Task | Status | Phase |
|---|------|--------|-------|
| 27-03 | ProcessPoolExecutor in engine.py | ⬜ Not started | 27 |

**After 27-03 completes:** Phase 27 verification → milestone audit → v0.15.0 release.

---

## Deferred (post-v0.15.0)

- 1.3 Remove R backends (v0.17.0)
- 3.4 Single eigendecomposition for SKAT-O
- 3.5 Sparse genotype matrices
- 4.2 AKT/PLINK wrappers
- 4.3 Kinship / mixed models

---

## References

- [regenie (Mbatchou 2021)](https://www.nature.com/articles/s41588-021-00870-7) — gene-level parallelization pattern
- [SAIGE-GENE+ (Zhou 2022)](https://www.nature.com/articles/s41588-022-01178-w) — SPA for extreme tails, sparse genotypes
- [FastSKAT (Lumley 2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6129408/) — partial eigenvalue decomposition
- [Das 2024](https://arxiv.org/abs/2404.05062) — new methods for generalized chi-square distribution
- [MCMC-CE (2025)](https://www.biorxiv.org/content/10.1101/2025.03.16.643492v1) — extreme tail p-values
- [Lumley 2024 blog](https://notstatschat.rbind.io/2024/09/13/two-approaches-to-approximating-sums-of-chisquareds/) — SPA vs Liu guidance
- [Liu & Xie 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407498/) — ACAT/ACAT-V
- [STAARpipeline](https://www.nature.com/articles/s41592-022-01641-w) — ACAT-O as recommended omnibus
