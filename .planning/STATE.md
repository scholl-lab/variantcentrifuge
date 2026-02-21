# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-19)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.15.0 Modular Rare Variant Association Framework

## Current Position

Phase: 23 — PCA + Functional Weights + Allelic Series + JSON Config
Plan: 2/4 complete
Status: In progress — Phase 23 Plan 02 complete: CADD/REVEL functional weights, CLI args, stage integration, 45 tests, 1776 passing
Last activity: 2026-02-21 — Completed 23-02-PLAN.md (functional weights)

Progress: █████████████████░░░░ ~85% (Phases 18-22 complete, Phase 23 Plan 2/4 done)

## Milestone Overview

**v0.15.0 — Modular Rare Variant Association Framework** (Phases 18-23)

| Phase | Goal | Requirements | Status |
|-------|------|--------------|--------|
| 18. Foundation | Core abstractions + Fisher refactor; bit-identical output validation | CORE-01..08 (8) | Complete |
| 19. Covariate System + Burden Tests | Logistic/linear burden tests with covariate adjustment and genotype matrix builder | COV-01..04, BURDEN-01..03, WEIGHT-01..02 (9) | Complete ✓ |
| 20. R SKAT Backend | R SKAT via rpy2 as gold standard oracle; SKATBinary + moment adjustment | SKAT-01..04, SKAT-08..09 (6) | Complete ✓ |
| 21. Pure Python SKAT Backend | Davies ctypes + saddlepoint + Liu fallback; validated against R within 10% | SKAT-05..07, SKAT-10 (4) | Complete ✓ |
| 22. ACAT-O + Diagnostics | ACAT-O omnibus; single FDR; lambda_GC; QQ TSV; sample size warnings | OMNI-01..03, DIAG-01..03, DIAG-05..06 (8) | Complete ✓ |
| 23. PCA + Functional Weights + Allelic Series + JSON Config | PCA file loading + AKT stage; CADD/REVEL weights; COAST test; JSON config; matplotlib plots | DIAG-04, PCA-01..04, SERIES-01..02, CONFIG-01..02, WEIGHT-03..05 (12) | Pending |

**Total requirements:** 47 mapped across 6 phases (35 complete, 12 pending)
<!-- Note: 8 CORE + 9 Phase 19 + 6 Phase 20 + 4 Phase 21 + 8 Phase 22 = 35 complete -->

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
| TEST-03 | CI validation uses sm.OLS conf_int() reference, not true-beta coverage | 19-03 | Finite-sample bias causes true beta to fall outside 95% CI with n=100; statsmodels reference is correct check |
| TEST-04 | FIRTH_CONVERGE_FAIL test replaced with always-returns-TestResult invariant | 19-03 | All-carrier gene exposes sm.add_constant removing intercept edge case; robustness invariant is more useful |
| IMPL-09 | parse_gt_to_dosage returns (int\|None, bool) not int\|None | 19-01 | Multi-allelic flag needed to emit 'run bcftools norm' warning without second parse pass |
| IMPL-10 | load_covariates returns (np.ndarray, list[str]) tuple | 19-01 | Column names returned alongside matrix for diagnostics; callers can ignore second element |
| IMPL-11 | build_genotype_matrix: sample_mask is list[bool], all samples remain in geno | 19-01 | Callers (logistic burden test) decide whether to exclude high-missing samples |
| IMPL-12 | LogisticBurdenTest builds design matrix inline (not in genotype_matrix.py) | 19-02 | Firth + separation checks are logistic-specific; keeps coupling clean |
| IMPL-13 | Tiered sample size check at n_cases<10 aborts with logger.error() not exception | 19-02 | Stage returns context cleanly; exception would crash pipeline; logger.error signals severity |
| IMPL-14 | ~~linear_burden effect_size=beta; engine column named *_or contains beta~~ RESOLVED | 19-02 | Originally deferred to Phase 22; resolved early via effect_column_names() polymorphism (commit 48a6e68) |
| IMPL-15 | Burden tests report beta+SE (not OR); Fisher keeps OR columns | 19 (post-verify) | Per-unit burden OR misleadingly close to 1.0 due to Beta(MAF;1,25) weights; beta+SE matches SKAT/SAIGE-GENE convention |
| IMPL-16 | AssociationAnalysisStage recovers per-sample GT from context.variants_df | 19 (post-verify) | gene_burden_analysis at same level drops per-sample GT columns before association can use them |
| IMPL-17 | RSKATTest.check_dependencies() hardcodes backend='r'; no auto-detect at test level | 20-01 | RSKATTest IS the R SKAT test; auto-selection is at the factory level. Separate PurePythonSKATTest for Phase 21. |
| IMPL-18 | Extra columns written with bare key names (skat_o_rho, not skat_skat_o_rho) | 20-01 | Keys in TestResult.extra are already namespaced by test; double-prefixing would produce unreadable names |
| IMPL-19 | effect_column_names() return type is dict[str, str \| None] not dict[str, str] | 20-01 | SKAT has no effect size; all four slots are None. Type broadened to accommodate without breaking existing tests. |
| IMPL-20 | withCallingHandlers embedded in R code string (not rpy2 Python callbacks) | 20-02 | Cleaner for string-based R execution; avoids rpy2 Python-side callback wiring complexity |
| IMPL-21 | GC triggered in RSKATTest.run() every 100 genes (not only in finalize) | 20-02 | Prevents R heap accumulation during long runs; finalize only handles terminal cleanup |
| IMPL-22 | prepare()/finalize() as no-ops in AssociationTest ABC | 20-02 | Fisher/burden tests unaffected; only RSKATTest overrides for R-specific lifecycle management |
| IMPL-23 | parallel_safe=False on AssociationAnalysisStage is unconditional | 20-02 | Not gated on skat_backend config; rpy2 safety applies regardless of test mix in same invocation |
| TEST-05 | rpy2 mock hierarchy requires parent attribute linking: mock_rpy2.robjects = mock_ro | 20-03 | bare sys.modules injection fails for nested submodule imports inside method bodies; parent mock attribute must point to child mock |
| TEST-06 | NA_Real sentinel: create unique object() per test; inject via sys.modules rpy2.rinterface.NA_Real | 20-03 | identity check `p_val_r is NA_Real` requires same object; fresh sentinel per test prevents cross-test contamination |
| FIX-01 | rpy2 3.6.x removed `rpy2.__version__`; use `importlib.metadata.version("rpy2")` | 20 (live) | rpy2 3.6.4 raises AttributeError on `rpy2.__version__`; importlib.metadata is stdlib since 3.8 |
| FIX-02 | R cleanup pattern must be `'^\\\\._vc_'` not `'\\._vc_'` | 20 (live) | Unescaped dot in regex matches any char; caret anchors to variable name start |
| FIX-03 | SKAT-O rho: `param$rho_est` is optimal rho; `param$rho` is the search grid | 20 (live) | `param$rho[0]` always returns 0.0 (first grid value); `param$rho_est` is the actual estimate |
| FIX-04 | Remove `.tolist()` before `FloatVector()` — numpy arrays accepted directly | 20 (live) | Unnecessary copy; rpy2 FloatVector accepts numpy arrays natively |
| IMPL-24 | CFFI set_source header must use extern "C" brackets for C++/C linkage bridging | 21-01 | Without extern "C" in the CFFI wrapper's forward declaration, C++ name mangling makes qfc() unresolvable at link time |
| IMPL-25 | qfc.cpp R headers replaced with standard C++ headers (math identical) | 21-01 | <R.h> and "Rmath.h" were included but unused; standalone compilation requires standard headers only |
| IMPL-26 | compute_pvalue() uses proactive saddlepoint at p<=1e-5 even when Davies ifault=0 | 21-01 | GMMAT pattern: Davies can produce false convergence near integration singularity for extreme p-values |
| IMPL-27 | ~~SKAT-O uses minimum-p approach~~ RESOLVED: Full Lee et al. (2012) SKAT-O implemented | 21-02, post-verify | Analytical R.M^{1/2} eigenvalue computation + omnibus chi2(1) integration; matches R SKAT exactly on GCKD cohort |
| IMPL-28 | _parse_weights_beta moved to shared tests/_utils.py to prevent rpy2 transitive import | 21-02 | If it stayed in skat_r.py, importing it from skat_python.py would transitively import rpy2 |
| IMPL-29 | Backend-aware swap in from_names() runs BEFORE unknown-name check | 21-02 | Critical ordering: 'skat' must resolve to correct class before validation; swap changes registry[skat] target |
| IMPL-30 | Davies C ext returns CDF not SF for some eigenvalue/Q combinations | 21-03 | qfc() returns ~1.0 for large Q where true p is small; fallback chain handles via saddlepoint/Liu; tests use _liu_pvalue directly for chi2 ground-truth validation |
| FIX-05 | Zero-variant guard added before matrix_rank in _test_skat/_test_skato | 21-03 | np.linalg.matrix_rank raises ValueError on (n,0) matrices; guard returns p_value=None,skip_reason=rank_deficient |
| FIX-06 | SKAT projection: project Z_tilde (not G_w) through hat matrix; eigenvalues /2 | 21 post-verify | R SKAT projects diag(phi)@G_w, not G_w; K <- W/2 halves eigenvalues to match Q=score'score/2 |
| FIX-07 | Davies compute_pvalue matches R Get_PValue.Lambda: acc=1e-6, lim=10000, keep non-converged | 21 post-verify | R SKAT uses acc=1e-6 and lim=10000 (not 1e-9/1M); non-converged Davies kept if 0<p<=1 |
| IMPL-31 | SKAT-O eigenvalues via analytical R.M^{1/2} (not Cholesky) | 21 post-verify | R.M = (1-rho)*I + rho*J has known eigenvalues; sqrt computed analytically avoiding Cholesky instability at high rho; handles rho=1.0 |
| IMPL-32 | R's rho >= 0.999 capping applied in SKAT-O eigenvalue loop | 21 post-verify | R SKAT caps rho at 0.999 to avoid rank-deficient correlation matrix |
| IMPL-33 | cauchy_combination() uses 1/(p*pi) approximation for p < 1e-16 | 22-01 | tan((0.5-p)*pi) overflows to ±inf for tiny p; approximation is numerically equivalent (Liu & Xie 2020 Section 2.2) |
| IMPL-34 | ACAT-O NOT in _TEST_REGISTRY — post-loop meta-test only | 22-01 | ACAT-O has no genotype input; it combines primary test results; adding to registry would allow nonsensical from_names(['acat_o']) calls |
| IMPL-35 | Single valid p-value returns as pass-through in cauchy_combination() | 22-01 | k=1 case per CONTEXT.md decision; one test p-value is informative enough to surface without modification |
| IMPL-36 | gene_burden_data dicts use "GENE" (uppercase) key; results_df column is "gene" (lowercase) | 22-02 | All three aggregation paths set uppercase key; engine sets lowercase column name; lookup uses .get("GENE", .get("gene", "")) for robustness |
| IMPL-37 | QQ data sorted ascending by expected_neg_log10_p (non-significant end first) | 22-02 | Smallest expected values (bottom-left of QQ plot) come first; matches sequential rendering convention |
| IMPL-38 | _EXPECTED_CHI2_MEDIAN hardcoded (not computed at import) | 22-02 | chi2.ppf(0.5, df=1) = 0.45493642311957174; avoids scipy cold import overhead at module load |
| IMPL-39 | PCA matrix merged inline into local covariate_matrix — not stored in context | 23-01 | Follows genotype matrix memory invariant (never in PipelineContext); feeds gene_data["covariate_matrix"] naturally |
| IMPL-40 | Missing AKT binary raises ToolNotFoundError (hard error, not skip) | 23-01 | Silent skip would let users believe AKT ran when it didn't; hard error surfaces misconfiguration immediately |
| IMPL-41 | PLINK eigenvec always uses IID (column 2) not FID (column 1) | 23-01 | FID is family ID, not sample-level identifier; IID matches VCF sample names |
| IMPL-42 | PCA format detection via first-line heuristics (column count + numeric checks) | 23-01 | Three-way detection: #FID/FID header, two non-numeric columns (PLINK nohdr), one non-numeric column (AKT/generic) |
| IMPL-43 | Keyword-only kwargs on get_weights() for backward compatibility | 23-02 | cadd_scores/revel_scores/variant_effects/weight_params added as keyword-only with None defaults; all existing callers `get_weights(mafs, spec)` unchanged |
| IMPL-44 | Site-filter mask replicated in stage for annotation alignment | 23-02 | build_genotype_matrix applies keep_variants_mask internally; stage replicates same logic via parse_gt_to_dosage when len(mafs) < len(gene_df) to align annotation arrays |
| IMPL-45 | combined_weights prefers CADD over REVEL when both provided | 23-02 | Higher CADD phred = more damaging; CADD is primary functional score, REVEL secondary |
| IMPL-46 | Annotation extraction conditional on weight spec in ("cadd","revel","combined") | 23-02 | Zero overhead for beta:*/uniform specs; extraction only runs when functional annotations are actually needed |

### Architecture Invariants (from research)

- R backend: parallel_safe=False; rpy2 calls only from main thread (segfault risk otherwise)
- Binary traits: always SKATBinary — never continuous-trait SKAT on binary phenotypes
- Davies defaults: davies_pvalue() uses acc=1e-9, lim=1_000_000; compute_pvalue() (SKAT) uses acc=1e-6, lim=10_000 matching R SKAT Get_PValue.Lambda
- Covariate alignment: always reindex to vcf_samples order; assert no NaN after reindex
- FDR strategy: single pass on ACAT-O p-values across all genes (not per-test)
- Genotype matrix: never stored in PipelineContext (5K samples x 50K variants = 1.6 GB)
- Eigenvalue stability: scipy.linalg.eigh; threshold max(eigenvalues, 0); skip if matrix_rank < 2
- Python version: recommend bumping requires-python to >=3.11 (scipy 1.16 dropped 3.10)
- apply_correction([]) returns empty array (statsmodels multipletests raises ZeroDivisionError on empty; guarded in correction.py)
- Firth NR fallback: self-contained 130-line Newton-Raphson; no external package; step-halving on penalized log-likelihood
- Separation detection: both mle_retvals['converged']==False AND bse.max()>100 needed (statsmodels may not raise exception)
- P-value computation: always through compute_pvalue() — never call Liu/Kuonen directly in SKAT backend code

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See archived REQUIREMENTS.md Future Requirements.
- ~~**RESEARCH-01** (before Phase 20): Validate whether parallel_safe=False on the stage is sufficient for rpy2 thread safety~~ — RESOLVED: parallel_safe=False + _assert_main_thread() guard confirmed sufficient (unit tests verify RuntimeError from worker threads)
- ~~**RESEARCH-02** (before Phase 21): Saddlepoint approximation algorithm for middle tier of Davies fallback chain~~ — RESOLVED: Kuonen Lugannani-Rice saddlepoint implemented in 21-01 from GENESIS variantSetTests.R reference
- ~~**COLUMN-01** (Phase 22): Rename `linear_burden_or` column to `linear_burden_beta`~~ — DONE (resolved via effect_column_names() in commit 48a6e68)

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-21
Stopped at: Completed 23-02-PLAN.md (functional weights — cadd_weights, revel_weights, combined_weights, stage integration, 45 tests)
Resume file: None
Next: Execute Phase 23 Plan 03 (allelic series — COAST test)
