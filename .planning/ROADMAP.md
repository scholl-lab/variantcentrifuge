# Roadmap: VariantCentrifuge

## Milestones

- SHIPPED **v0.12.1 Baseline** — Phases 1-5 (shipped 2026-02-14, pre-GSD)
- SHIPPED **v0.13.0 Performance Optimization** — Phases 6-12 (shipped 2026-02-16) — [archive](milestones/v0.13.0-ROADMAP.md)
- SHIPPED **v0.14.0 Report UX Overhaul** — Phases 13-17 (shipped 2026-02-19)
- ACTIVE **v0.15.0 Modular Rare Variant Association Framework** — Phases 18-23

## Phases

<details>
<summary>v0.12.1 Baseline (Phases 1-5) — SHIPPED 2026-02-14</summary>

**Milestone Goal:** Full-featured variant analysis pipeline with two architectures, 40+ stages, inheritance analysis, gene burden, scoring, multi-format output. All 30 historical issues resolved. 1035 tests passing. CI/CD with Docker, docs, and multi-platform testing.

**Note:** Phases 1-5 represent completed work before GSD tracking. Details available in git history.

</details>

<details>
<summary>v0.13.0 Performance Optimization (Phases 6-12) — SHIPPED 2026-02-16</summary>

- [x] Phase 6: Benchmark Framework (4/4 plans) — completed 2026-02-14
- [x] Phase 7: Quick Wins - Tier 1 (3/3 plans) — completed 2026-02-14
- [x] Phase 8: DataFrame Optimization (4/4 plans) — completed 2026-02-14
- [x] Phase 9: Inheritance Analysis Optimization (5/5 plans) — completed 2026-02-14
- [x] Phase 10: Output Optimization (3/3 plans) — completed 2026-02-15
- [x] Phase 11: Pipeline I/O Elimination (3/3 plans) — completed 2026-02-15
- [x] Phase 12: Parallelization & Chunking (4/4 plans) — completed 2026-02-16

Full details: [milestones/v0.13.0-ROADMAP.md](milestones/v0.13.0-ROADMAP.md)

</details>

<details>
<summary>v0.14.0 Report UX Overhaul (Phases 13-17) — SHIPPED 2026-02-19</summary>

- [x] Phase 13: JS Stack Modernization (3/3 plans) — completed 2026-02-16
- [x] Phase 14: Information Hierarchy and Semantic Color Coding (3/3 plans) — completed 2026-02-16
- [x] Phase 15: Table Redesign (3/3 plans) — completed 2026-02-17
- [x] Phase 16: Column-Level Filtering and Visualization (3/3 plans) — completed 2026-02-18
- [x] Phase 17: Accessibility and Print/PDF (3/3 plans) — completed 2026-02-17

</details>

### v0.15.0 Modular Rare Variant Association Framework (Phases 18-23)

Upgrade the existing Fisher's exact gene burden pipeline to a multi-test association engine supporting SKAT, SKAT-O, ACAT-O, logistic/linear burden tests, covariate adjustment, PCA integration, and variant weighting — all backward compatible with existing `--perform-gene-burden` behavior.

---

#### Phase 18: Foundation — Core Abstractions and Fisher Refactor

**Goal:** Users can run `--perform-association --association-tests fisher` and receive output bit-identical to `--perform-gene-burden`, establishing the association package skeleton and pipeline integration that all subsequent phases build on.

**Dependencies:** None (foundation phase)

**Requirements:** CORE-01, CORE-02, CORE-03, CORE-04, CORE-05, CORE-06, CORE-07, CORE-08

**Plans:** 4 plans
Plans:
- [x] 18-01-PLAN.md — association/ package skeleton (base.py, engine.py, correction.py, fisher.py)
- [x] 18-02-PLAN.md — Pipeline integration (AssociationAnalysisStage, PipelineContext, stage registry, gene_burden.py correction rewire)
- [x] 18-03-PLAN.md — CLI args (--perform-association, --association-tests, --skat-backend) and ExcelReportStage Association sheet
- [x] 18-04-PLAN.md — Unit tests: bit-identity Fisher parity, edge cases, correction, stage behavior

**Success Criteria:**

1. Running `--perform-association --association-tests fisher` produces output byte-for-byte identical to `--perform-gene-burden` across all existing integration test fixtures
2. Running `--perform-gene-burden` without `--perform-association` produces zero behavioral change — GeneBurdenAnalysisStage executes normally and its existing tests all pass
3. The `association/` package is importable without R, rpy2, or matplotlib installed — all optional dependencies are imported lazily inside the functions that need them
4. `ExcelReportStage` generates an "Association" sheet when `--perform-association` is active
5. FDR and Bonferroni correction functions in `association/correction.py` are re-exported from `gene_burden.py` with no change to existing callers

---

#### Phase 19: Covariate System and Burden Tests

**Goal:** Users can run logistic or linear regression burden tests with covariate adjustment, and the genotype matrix builder handles all real-world genotype encodings correctly.

**Dependencies:** Phase 18 (AssociationEngine, AssociationTest ABC, pipeline integration)

**Requirements:** COV-01, COV-02, COV-03, COV-04, BURDEN-01, BURDEN-02, BURDEN-03, WEIGHT-01, WEIGHT-02

**Plans:** 3 plans estimated
- 19-01: covariates.py — covariate file loading, sample ID alignment with VCF order, one-hot encoding, multicollinearity check; genotype matrix builder with multi-allelic, missing, and phased genotype handling
- 19-02: logistic_burden.py and linear_burden.py tests; Beta(MAF;1,25) and uniform weight modules; CLI args --covariate-file, --covariates, --trait-type, --variant-weights
- 19-03: Unit tests for covariate alignment with shuffled sample order, genotype edge cases (1/2, ./., 0|1), logistic/linear burden test outputs validated against manual statsmodels calls

**Success Criteria:**

1. Providing a covariate file with rows in a different order than the VCF sample list produces identical results to a covariate file in VCF sample order — the reindex assertion catches any mismatch
2. The logistic burden test reports OR + 95% CI for a binary trait with covariates, and results match manually running `statsmodels.Logit` on the same genotype matrix and covariate matrix
3. The linear burden test reports beta + SE for a quantitative trait, matching `statsmodels.OLS` on the same inputs
4. Beta(MAF; 1, 25) weights up-weight rare variants relative to common ones; uniform weights produce the same result as setting all weights to 1.0
5. Genotype strings `1/2`, `./.`, and `0|1` are parsed correctly — multi-allelic hets counted as 1 dosage, missing values imputed to 2*MAF, phased genotypes summed to 0/1/2

---

#### Phase 20: R SKAT Backend

**Goal:** Users with R and the SKAT package installed can run SKAT and SKAT-O via rpy2, with SKATBinary used automatically for binary traits, moment adjustment enabled by default for small samples, and R memory managed explicitly to prevent heap exhaustion across thousands of genes.

**Dependencies:** Phase 19 (genotype matrix builder, covariate matrix, weight vector — all inputs to SKAT)

**Requirements:** SKAT-01, SKAT-02, SKAT-03, SKAT-04, SKAT-08, SKAT-09

**Plans:** 3 plans estimated
- 20-01: backends/base.py SKATBackend ABC and NullModel container; backends/__init__.py get_skat_backend() factory (lazy, never at module import); r_backend.py RSKATBackend with rpy2 import guard, R/SKAT detection, graceful fallback
- 20-02: RSKATBackend.fit_null_model() using SKAT_Null_Model_MomentAdjust for binary traits; RSKATBackend.test_gene() dispatching to SKATBinary vs SKAT by trait type; SKAT-O with method="optimal.adj"; parallel_safe=False on AssociationAnalysisStage when R backend active; explicit del + gc() every 100 genes
- 20-03: Synthetic test fixtures (100 cases/100 controls, 50 genes, 5 genes with injected burden signal); validate R backend detects signal genes at p < 0.05; validate R memory usage stays bounded across 500-gene run

**Success Criteria:**

1. On a system with R and the SKAT package, `--skat-backend r` runs SKAT and SKAT-O and reports p-values; on a system without R, the same command falls back gracefully to the Python backend with an informative log message
2. A binary trait phenotype always uses SKATBinary — the continuous-trait SKAT formulation is never called for binary outcomes, verified by inspecting which R function was invoked
3. SKAT-O reports the optimal rho value alongside the p-value, and uses `method="optimal.adj"` correction (not the uncorrected minimum p)
4. Running the R backend across 500 synthetic genes does not cause R heap exhaustion — R objects are deleted after each gene and `gc()` is called every 100 genes
5. Calling the R backend from a ThreadPoolExecutor worker thread raises an explicit error rather than causing a segfault or silent crash

---

#### Phase 21: Pure Python SKAT Backend

**Goal:** Users without R can run SKAT and SKAT-O via a pure Python implementation that matches R output within 10% relative difference on log10(p), using Davies ctypes for exact p-values with a Liu moment-matching fallback when compilation is unavailable.

**Dependencies:** Phase 20 (R backend as correctness oracle; must run both backends simultaneously for validation)

**Requirements:** SKAT-05, SKAT-06, SKAT-07, SKAT-10

**Plans:** 3 plans estimated
- 21-01: backends/davies.py — Liu moment-matching fallback first (pure scipy, always available); then ctypes lazy compilation of bundled qfc.c with acc=1e-9, lim=10^6; fallback chain: Davies -> saddlepoint -> Liu; p_method metadata in every result; data/qfc.c bundled source; pyproject.toml artifacts field
- 21-02: backends/python_backend.py PythonSKATBackend — score test (SKAT linear/logistic formulation from null model residuals); SKAT eigenvalue computation with scipy.linalg.eigh, threshold max(eigenvalues,0), skip if matrix_rank < 2; SKAT-O rho grid search only after SKAT validates against R
- 21-03: Validation test suite comparing Python vs R p-values across 50+ genes on same synthetic fixtures — assert 10% relative difference on log10(p) for p > 0.001; document expected larger divergence for p < 0.001 in test comments

**Success Criteria:**

1. Python SKAT p-values are within 10% relative difference of R SKAT p-values on log10(p) scale for p-values greater than 0.001, validated across at least 50 synthetic genes
2. Davies ctypes compilation succeeds on both Linux and Windows; when gcc is absent, Liu moment-matching activates automatically and the `p_method` column records "liu"
3. The `p_method` output column records "davies", "saddlepoint", or "liu" for every gene, giving users visibility into which numerical method was used
4. SKAT skips genes where the kernel matrix rank is less than 2 rather than returning a spurious p-value, and reports `p_value=NA` with a diagnostic warning for those genes
5. Running `--skat-backend python` without R installed completes successfully and produces association results — the entire pipeline is functional without any R dependency

---

#### Phase 22: ACAT-O and Diagnostics

**Goal:** Users receive omnibus ACAT-O p-values combining burden and SKAT results per gene, a single FDR-corrected set of ACAT-O p-values across all genes, and a diagnostics directory containing lambda_GC per test and QQ plot data TSV.

**Dependencies:** Phase 19 (burden p-values), Phase 20 or 21 (SKAT p-values); ACAT-O requires both test types to exist

**Requirements:** OMNI-01, OMNI-02, OMNI-03, DIAG-01, DIAG-02, DIAG-03, DIAG-05, DIAG-06

**Plans:** 3 plans estimated
- 22-01: association/tests/acat.py — ACAT-V Cauchy combination of marginal score test p-values per gene; ACAT-O combining burden + SKAT + ACAT-V p-values; single FDR applied to ACAT-O across all genes (not per-test)
- 22-02: association/diagnostics.py — lambda_GC per test (median chi2 / 0.4549); QQ data TSV (observed vs expected -log10(p)); per-gene TSV with all standard output columns; CLI arg --diagnostics-output; sample size warnings when n_cases < 200, case:control > 1:20, or case_carriers < 10 per gene
- 22-03: Tests validating ACAT-O Cauchy formula against published values; lambda_GC within [0.95, 1.05] on permuted null phenotype; output TSV column completeness check

**Success Criteria:**

1. The per-gene association TSV contains all standard columns: fisher_p, burden_p, skat_p, skat_o_p, acat_o_p, effect sizes, confidence intervals, and variant/carrier counts
2. ACAT-O p-values are computed using a single Cauchy combination per gene and a single FDR correction pass across all genes — not separate FDR corrections per test type
3. The `--diagnostics-output` directory contains `lambda_gc.txt` with one genomic inflation factor per active test and `qq_data.tsv` with observed vs expected -log10(p) columns
4. The pipeline logs a warning for any gene where `case_carriers < 10` and emits a summary warning when `n_cases < 200` or `case:control ratio > 1:20`, flagging those genes in the output TSV
5. Lambda_GC computed on a permuted null phenotype (random case/control assignment) falls within [0.95, 1.05] for Fisher, burden, and SKAT tests — indicating no systematic inflation in the implementation

---

#### Phase 23: PCA Integration, Functional Weights, Allelic Series, and JSON Config

**Goal:** Users can supply a pre-computed PCA file or trigger AKT-based PCA computation as a pipeline stage, apply CADD/REVEL functional variant weights, run the COAST allelic series test, and specify all association options via a JSON config file — with optional matplotlib QQ plot generation for users who have it installed.

**Dependencies:** Phase 19 (covariate system for PCA merging; weights architecture for functional weights), Phase 22 (diagnostics architecture for QQ plot)

**Requirements:** DIAG-04, PCA-01, PCA-02, PCA-03, PCA-04, SERIES-01, SERIES-02, CONFIG-01, CONFIG-02, WEIGHT-03, WEIGHT-04, WEIGHT-05

**Plans:** 4 plans estimated
- 23-01: association/pca.py — PLINK .eigenvec parser, AKT output parser, generic TSV parser; sample ID alignment validation; PCAComputationStage in processing_stages.py wrapping akt subprocess; CLI args --pca-file, --pca-tool, --pca-components; covariates.py updated to merge PCA eigenvectors; warn if >20 components requested
- 23-02: Functional weights — weights/functional_weights.py with CADD-normalized and REVEL-based schemes; weights/combined.py for Beta(MAF) x functional score; CLI arg --variant-weights cadd/revel/combined with --variant-weight-params
- 23-03: COAST allelic series test — tests/allelic_series.py classifying variants as BMV/DMV/PTV from PolyPhen/SIFT columns; ordered alternative test; configurable category weights (default w=1,2,3)
- 23-04: JSON config mode — association/config.py JSON loader mapping all CLI options to AssociationConfig fields; CLI arg --association-config; optional matplotlib QQ plot PNG/SVG in diagnostics.py with lazy import and matplotlib.use("Agg") for headless HPC

**Success Criteria:**

1. Providing `--pca-file pca.eigenvec` merges the requested number of PCs as additional covariate columns, with sample alignment verified against VCF sample order; requesting >20 PCs triggers a logged warning
2. Running `--pca-tool akt` invokes AKT as a pipeline stage, stores eigenvectors in context, and uses them in subsequent association tests — when AKT is not in PATH, the stage skips gracefully with an informative warning
3. Setting `--variant-weights cadd` applies CADD-normalized weights to variants, and `--variant-weights combined` applies Beta(MAF) x CADD weights; both produce different p-values than uniform weights on the same data, confirming the weights are applied
4. Running `--association-config config.json` with a valid JSON file produces the same results as passing equivalent CLI flags directly — the JSON config is a complete alternative to CLI args for all association options
5. When matplotlib is installed, `--diagnostics-output` produces a QQ plot PNG in addition to the data TSV; when matplotlib is absent, the pipeline completes without error and logs a single INFO message noting the plot was skipped

---

## Progress

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1-5. Baseline | v0.12.1 | N/A | Complete | 2026-02-14 |
| 6. Benchmark Framework | v0.13.0 | 4/4 | Complete | 2026-02-14 |
| 7. Quick Wins - Tier 1 | v0.13.0 | 3/3 | Complete | 2026-02-14 |
| 8. DataFrame Optimization | v0.13.0 | 4/4 | Complete | 2026-02-14 |
| 9. Inheritance Analysis Optimization | v0.13.0 | 5/5 | Complete | 2026-02-14 |
| 10. Output Optimization | v0.13.0 | 3/3 | Complete | 2026-02-15 |
| 11. Pipeline I/O Elimination | v0.13.0 | 3/3 | Complete | 2026-02-15 |
| 12. Parallelization & Chunking | v0.13.0 | 4/4 | Complete | 2026-02-16 |
| 13. JS Stack Modernization | v0.14.0 | 3/3 | Complete | 2026-02-16 |
| 14. Information Hierarchy and Semantic Color Coding | v0.14.0 | 3/3 | Complete | 2026-02-16 |
| 15. Table Redesign | v0.14.0 | 3/3 | Complete | 2026-02-17 |
| 16. Column-Level Filtering and Visualization | v0.14.0 | 3/3 | Complete | 2026-02-18 |
| 17. Accessibility and Print/PDF | v0.14.0 | 3/3 | Complete | 2026-02-17 |
| 18. Foundation — Core Abstractions and Fisher Refactor | v0.15.0 | 4/4 | Complete | 2026-02-19 |
| 19. Covariate System and Burden Tests | v0.15.0 | 0/3 | Pending | — |
| 20. R SKAT Backend | v0.15.0 | 0/3 | Pending | — |
| 21. Pure Python SKAT Backend | v0.15.0 | 0/3 | Pending | — |
| 22. ACAT-O and Diagnostics | v0.15.0 | 0/3 | Pending | — |
| 23. PCA Integration, Functional Weights, Allelic Series, and JSON Config | v0.15.0 | 0/4 | Pending | — |
