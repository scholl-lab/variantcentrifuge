# Roadmap: VariantCentrifuge

## Milestones

- [x] **v0.12.1 Baseline** - Phases 1-5 (shipped 2026-02-14)
- [x] **v0.13.0 Performance Optimization** - Phases 6-12 (shipped 2026-02-16)
- [x] **v0.14.0 Report UX Overhaul** - Phases 13-17 (shipped 2026-02-19)
- [x] **v0.15.0 Modular Rare Variant Association Framework** - Phases 18-29 (shipped 2026-02-23)
- [ ] **v0.16.0 Association Hardening & Multi-Cohort Features** - Phases 30-37 (in progress)

## Phases

<details>
<summary>Phases 1-29 (v0.12.1 through v0.15.0) — SHIPPED</summary>

Phases 1-5: v0.12.1 Baseline (pre-GSD tracking)
Phases 6-12: v0.13.0 Performance Optimization
Phases 13-17: v0.14.0 Report UX Overhaul
Phases 18-29: v0.15.0 Modular Rare Variant Association Framework

See MILESTONES.md for details.

</details>

---

### v0.16.0 Association Hardening & Multi-Cohort Features (In Progress)

**Milestone Goal:** Harden the v0.15.0 association framework for real-world multi-cohort use — fix COAST p=None bugs that block scientific validation, add gene-level FDR priors, case-confidence sample weighting, BED region restriction, and sparse matrix support. Every new feature is backward-compatible; existing CLI and output formats unchanged.

- [x] **Phase 30: Dead Code Cleanup** — Remove 8 vestigial stage classes and normalize legacy naming artifacts
- [x] **Phase 31: COAST Fix** — Fix genotype matrix reconstruction, partial-category fallback, and configurable classification scoring
- [x] **Phase 32: Region Restriction and PCA Wiring** — Add BED-based region prefilter and wire PCAComputationStage into the pipeline
- [x] **Phase 33: Gene-Level FDR Weighting** — Implement weighted Benjamini-Hochberg with per-gene biological priors
- [x] **Phase 34: Tech Debt** — Fix association config mapping, generate COAST golden values, fix column naming, clarify diagnostics
- [ ] **Phase 35: Case-Confidence Weights** — Add per-sample confidence weights to burden and SKAT null models
- [ ] **Phase 36: Performance — Sparse Genotype Matrices** — Add opt-in sparse matrix path for large rare-variant cohorts
- [x] **Phase 37: Association Resource Management & Memory Streaming** — Shared ResourceManager in PipelineContext, eliminate GT column drop/recover antipattern, stream genotype matrices to prevent OOM

---

## Phase Details

### Phase 30: Dead Code Cleanup

**Goal:** The codebase is free of the 8 dead stage classes and all vestigial "refactored_pipeline" / "old pipeline" naming artifacts that accumulated before v0.15.0.

**Depends on:** Nothing (standalone cleanup; reduces noise for subsequent phases)

**Requirements:** CLEAN-01, CLEAN-02, CLEAN-03, CLEAN-04, CLEAN-05, CLEAN-06, TD-06

**Success Criteria** (what must be TRUE):
1. Running `make test-fast` passes with 0 failures after each stage class is removed (no silent dependencies were missed)
2. No string "refactored_pipeline" appears anywhere in pipeline.py, runner.py, or checkpoint files written by the pipeline
3. No "old pipeline" or "Refactored" comments or docstrings remain in error_handling.py, analyze_variants.py, or pipeline stage files
4. Stage registry imports and `__all__` entries are consistent — no reference to removed stage class names exists anywhere in the codebase

**Plans:** 1 plan

Plans:
- [x] 30-01-PLAN.md — Remove dead stage classes, synchronize registry/imports, normalize naming artifacts

---

### Phase 31: COAST Fix

**Goal:** COAST produces valid p-values for real-world cohorts where most genes lack all three variant categories (BMV/DMV/PTV), genotype matrices are reliably available to the test, and multi-transcript SnpEff effect strings are handled correctly.

**Depends on:** Nothing (bug fix; establishes working baseline needed before Phase 35)

**Requirements:** COAST-01, COAST-02, COAST-03, COAST-04, COAST-05, COAST-06, COAST-07

**Success Criteria** (what must be TRUE):
1. COAST returns a numeric p-value (not None) for genes where only 1 or 2 of the 3 variant categories are present, with a warning logged naming the missing categories
2. COAST classification correctly handles SnpEff "&"-concatenated effect strings (e.g., `stop_gained&splice_region_variant` is recognized as PTV)
3. When `--association-tests coast` is specified, SIFT and PolyPhen fields are present in the extracted data without requiring the user to add them to field lists manually
4. User can specify `--coast-classification cadd` and the pipeline uses CADD-score thresholds for BMV/DMV/PTV classification without error
5. A `scoring/coast_classification/` config with BMV/DMV/PTV thresholds exists and is the authoritative source for classification logic

**Plans:** 2 plans

Plans:
- [x] 31-01-PLAN.md — Fix genotype matrix fallback (COAST-01) and partial-category logic (COAST-03)
- [x] 31-02-PLAN.md — Fix &-concatenated effect strings, configurable classification scoring, auto-field injection (COAST-02, COAST-04, COAST-05, COAST-06, COAST-07)

---

### Phase 32: Region Restriction and PCA Wiring

**Goal:** Users can restrict variant extraction to callable regions via a BED file, and PCA covariates computed by the pipeline are automatically available to association tests without manual file passing.

**Depends on:** Nothing (both features are standalone)

**Requirements:** REGION-01, REGION-02, REGION-03, REGION-04, PCA-01, PCA-02, PCA-03

**Success Criteria** (what must be TRUE):
1. Passing `--regions-bed capture_kit.bed` restricts bcftools variant extraction to those regions; variants outside the BED do not appear in any output
2. If chromosome naming mismatches between the restriction BED and the VCF are detected, the pipeline exits with a clear error message (not silent 0-variant output)
3. When `--pca akt` is set, PCA eigenvectors are computed automatically during the pipeline run and appear as covariates in association test output without any additional user steps
4. The restriction BED intersection happens once, upstream of chunk splitting, so all gene chunks respect the same region mask

**Plans:** 2 plans

Plans:
- [x] 32-01-PLAN.md — Add --regions-bed CLI flag and BED intersection in GeneBedCreationStage with chr mismatch detection
- [x] 32-02-PLAN.md — Add PCAComputationStage, unified --pca CLI flag, pipeline wiring, AssociationAnalysisStage integration

---

### Phase 33: Gene-Level FDR Weighting

**Goal:** Users with biological prior knowledge (pLI scores, GWAS signals) can increase statistical power by providing per-gene weights for weighted Benjamini-Hochberg FDR correction, with correct weight normalization and transparent diagnostics.

**Depends on:** Nothing (standalone; pure correction.py addition)

**Requirements:** FDR-01, FDR-02, FDR-03, FDR-05, FDR-06 (FDR-04 IHW excluded per context discussion)

**Success Criteria** (what must be TRUE):
1. Passing `--gene-prior-weights weights.tsv` applies weighted BH with weights renormalized to mean=1.0 internally; the q-values differ from unweighted BH in proportion to weight variation
2. Genes absent from the weight file receive weight=1.0 (neutral); if more than 50% of tested genes are missing from the file, a warning is emitted
3. Effective number of tests (`sum(w)^2 / sum(w^2)`) appears in diagnostics output
4. Running without `--gene-prior-weights` produces identical results to the current plain BH (backward compatible)
5. IHW is NOT implemented in this phase — no CLI flag, no stub, no error message (per CONTEXT.md decision)

**Plans:** 1 plan

Plans:
- [x] 33-01-PLAN.md — Implement weighted BH in correction.py, CLI flag, weight file loader, and diagnostics

---

### Phase 34: Tech Debt

**Goal:** The association subsystem config mapping is correct, COAST has R golden values for regression testing, SKAT-O column naming is consistent, and Fisher lambda_GC scope is documented.

**Depends on:** Phase 31 (TD-03 COAST golden values require a working COAST implementation)

**Requirements:** TD-02, TD-03, TD-04, TD-05

**Success Criteria** (what must be TRUE):
1. `create_stages_from_config()` correctly activates association and gene burden stages when the corresponding config dict keys are set (confirmed by integration test)
2. COAST R golden p-values are hardcoded in the test file and a CI test asserts Python output is within 10% relative difference
3. SKAT-O results output uses `skat_o_pvalue` as the column name (not a variant of `skat_o_p`)
4. The Fisher lambda_GC behavior is documented in code comments: excluded from Fisher p-values, applied only to score-based tests (SKAT, burden)

**Plans:** 3 plans

Plans:
- [x] 34-01-PLAN.md — Fix create_stages_from_config mapping (TD-02) and document lambda_GC Fisher exemption (TD-05)
- [x] 34-02-PLAN.md — Standardize column naming: _pvalue/_qvalue convention across engine, diagnostics, tests, and docs (TD-04)
- [x] 34-03-PLAN.md — Generate COAST R golden values from GCKD data and add 10% tolerance comparison test (TD-03)

---

### Phase 35: Case-Confidence Weights

**Goal:** Users can provide per-sample confidence weights to improve statistical power in cohorts with uncertain phenotyping (EHR-derived diagnoses), with correct null model integration validated by permutation testing.

**Depends on:** Phase 31 (COAST path must work to validate weighted COAST runs end-to-end)

**Requirements:** CONF-01, CONF-02, CONF-03, CONF-04, CONF-05, CONF-06, CONF-07

**Success Criteria** (what must be TRUE):
1. Passing `--case-confidence file --case-confidence-column confidence_score` loads per-sample weights from the phenotype file and applies them to logistic burden via `var_weights` in the GLM null model
2. Passing `--case-confidence count` auto-computes confidence weights from HPO term overlap counts
3. Linear burden test uses WLS (weighted least squares) with sample weights; logistic burden uses weighted GLM
4. Fisher's exact test ignores sample weights and logs an info message explaining why (integer counts only)
5. Effective sample size (`sum(w)^2 / sum(w^2)`) appears in diagnostics output
6. `--case-confidence equal` (default) produces uniform weights and identical results to the current unweighted run (backward compatible)

**Plans:** TBD

Plans:
- [ ] 35-01: Implement sample_weights.py, CLI flags, and weighted GLM null model integration
- [ ] 35-02: Permutation validation — confirm lambda_GC ~1.0 on permuted phenotype with weights held fixed

---

### Phase 36: Performance — Sparse Genotype Matrices

**Goal:** Users working with large rare-variant cohorts can opt into sparse genotype matrix representation that reduces memory usage without changing p-value results, validated by density profiling and numerical equivalence tests.

**Depends on:** Nothing (opt-in feature; no correctness dependency)

**Requirements:** PERF-01, PERF-02, PERF-03

**Success Criteria** (what must be TRUE):
1. Passing `--sparse-genotypes` stores genotype matrices as `scipy.sparse.csr_array` when enabled; memory consumption is measurably lower than dense representation at >2 standard GCKD gene windows
2. SKAT kernel computation produces numerically identical p-values (within 1e-10 relative) whether the input genotype matrix is sparse or dense
3. SKAT-O eigendecomposition behavior is profiled and documented: either confirmed that the existing precomputed `A` matrix is already the single-decomposition optimization (ticket closed as done), or a validated alternative is implemented

**Plans:** TBD

Plans:
- [ ] 36-01: Profile SKAT-O eigendecomposition and clarify intent; implement sparse genotype matrix opt-in path

---

### Phase 37: Association Resource Management & Memory Streaming

**Goal:** The pipeline has a single shared ResourceManager in PipelineContext used by all stages, the GT column drop/recover antipattern is eliminated so genotype data flows cleanly through analysis stages, and genotype matrices are streamed per-gene (build-test-discard) to prevent OOM on large panels.

**Depends on:** Phase 34 (builds on Fix 1-3 already applied to genotype_matrix.py, analysis_stages.py, cli.py)

**Requirements:** PERF-04, PERF-05, PERF-06

**Success Criteria** (what must be TRUE):
1. All stages use `context.resource_manager.auto_workers()` for thread/worker allocation — no stage creates its own ResourceManager instance
2. Per-sample GT columns are preserved through analysis stages and only dropped at output time (TSV/Excel); no reconstruct/recover cycle exists
3. Genotype matrices are built per-gene inside the association engine loop and released after test execution — peak memory for 5K genes x 5K samples stays under 8 GB (vs 21+ GB before)
4. Running without `--association-workers` uses the `--threads` value (already implemented); ResourceManager gates all parallelism decisions

**Plans:** 3 plans

Plans:
- [x] 37-01-PLAN.md — Shared ResourceManager in PipelineContext (PERF-04)
- [x] 37-02-PLAN.md — GT lifecycle cleanup: eliminate drop/recover antipattern (PERF-05)
- [x] 37-03-PLAN.md — Streaming genotype matrix construction per-gene (PERF-06)

**Details:**
See `.planning/association-performance-investigation.md` for full analysis.
Fixes 1-3 (vectorize GT parsing, groupby, auto-workers) already applied.
This phase covers Fix 4 (shared ResourceManager), Fix 5 (GT lifecycle), Fix 6 (streaming matrices).

---

## Progress

**Execution Order:** 30 -> 31 -> 32 -> 33 -> 34 -> 35 -> 36 -> 37

Note: Phases 32 and 33 are independent of Phase 31 and of each other; they may be developed in parallel if needed.

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 30. Dead Code Cleanup | v0.16.0 | 1/1 | Complete | 2026-02-23 |
| 31. COAST Fix | v0.16.0 | 2/2 | Complete | 2026-02-23 |
| 32. Region Restriction and PCA Wiring | v0.16.0 | 2/2 | Complete | 2026-02-23 |
| 33. Gene-Level FDR Weighting | v0.16.0 | 1/1 | Complete | 2026-02-24 |
| 34. Tech Debt | v0.16.0 | 3/3 | Complete | 2026-02-24 |
| 35. Case-Confidence Weights | v0.16.0 | 0/2 | Not started | - |
| 36. Performance — Sparse Genotype Matrices | v0.16.0 | 0/1 | Not started | - |
| 37. Association Resource Management & Memory Streaming | v0.16.0 | 3/3 | Complete | 2026-02-25 |
