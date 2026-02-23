# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-19)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** v0.15.0 Modular Rare Variant Association Framework

## Current Position

Phase: 28 — Tech Debt Cleanup
Plan: 1/2 complete
Status: In progress
Last activity: 2026-02-23 — Completed 28-01-PLAN.md (RSKATTest parallel_safe + lambda_gc gap)

Progress: ████████████████████░░░ 85% (Phases 18-27 complete; 28-01 complete; 28-02 and 29 pending)

## Milestone Overview

**v0.15.0 — Modular Rare Variant Association Framework** (Phases 18-24)

| Phase | Goal | Requirements | Status |
|-------|------|--------------|--------|
| 18. Foundation | Core abstractions + Fisher refactor; bit-identical output validation | CORE-01..08 (8) | Complete |
| 19. Covariate System + Burden Tests | Logistic/linear burden tests with covariate adjustment and genotype matrix builder | COV-01..04, BURDEN-01..03, WEIGHT-01..02 (9) | Complete ✓ |
| 20. R SKAT Backend | R SKAT via rpy2 as gold standard oracle; SKATBinary + moment adjustment | SKAT-01..04, SKAT-08..09 (6) | Complete ✓ |
| 21. Pure Python SKAT Backend | Davies ctypes + saddlepoint + Liu fallback; validated against R within 10% | SKAT-05..07, SKAT-10 (4) | Complete ✓ |
| 22. ACAT-O + Diagnostics | ACAT-O omnibus; single FDR; lambda_GC; QQ TSV; sample size warnings | OMNI-01..03, DIAG-01..03, DIAG-05..06 (8) | Complete ✓ |
| 23. PCA + Functional Weights + Allelic Series + JSON Config | PCA file loading + AKT stage; CADD/REVEL weights; COAST test; JSON config; matplotlib plots | DIAG-04, PCA-01..04, SERIES-01..02, CONFIG-01..02, WEIGHT-03..05 (12) | Complete ✓ |
| 24. Pure Python COAST Backend | Pure Python COAST matching R AllelicSeries; parallel_safe=True; no R dependency | COAST-PY-01..05 (5) | Complete ✓ |
| 25. Python Default + Quick Wins | Python default backends, R deprecated, saddlepoint fallback, ACAT-V | — | Complete ✓ |
| 26. Documentation | Association testing guide, update existing docs, API stubs, changelog | TBD | Complete ✓ |
| 27. Performance Optimizations | GL quadrature (46x SKAT-O), ProcessPoolExecutor gene parallelization | — | Complete ✓ |
| 28. Tech Debt Cleanup | RSKATTest attribute, CLI args for skat_method + diagnostic thresholds | Audit gaps | Pending |
| 29. Classic Pipeline Deprecation | Remove pipeline.py, make stage-based default, remove --use-new-pipeline | DEPR-01 | Pending |

**Total requirements:** 52 mapped across 7 phases (52 complete); phases 25-27 requirements TBD

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
| IMPL-47 | classify_variants returns include_mask: missense without SIFT/PolyPhen gets code 0, excluded from COAST only | 23-03 | Preserves SKAT/burden power while enforcing COAST ordered-alternative assumption |
| IMPL-48 | gene_df stored in gene_data dict (reset_index) for COASTTest annotation column access | 23-03 | COASTTest needs per-variant EFFECT/IMPACT/SIFT/PolyPhen columns not in standard aggregation keys |
| IMPL-49 | Annotation/genotype mismatch skips gene with p_value=None | 23-03 | Site-filter removes variants from build_genotype_matrix but gene_df has all variants; shape check detects misalignment |
| IMPL-50 | p_value=None when any COAST category (BMV/DMV/PTV) missing | 23-03 | Ordered-alternative test requires all 3 categories; missing any produces spurious or undefined results |
| IMPL-51 | VALID_ASSOCIATION_KEYS uses unprefixed field names; CLI prefixed keys mapped via _get(cli_key, json_key) | 23-04 | Separates CLI and JSON key naming conventions; e.g. CLI "association_min_cases" maps to JSON "min_cases" |
| IMPL-52 | nullable=True/False in _get() closure: None in context.config means "not set" for nullable; key presence for non-nullable | 23-04 | Correctly handles distinction between CLI writing None (flag absent) vs CLI writing non-None (flag set) |
| IMPL-53 | matplotlib mock uses sys.modules["matplotlib"]=None; not builtins.__import__ patch | 23-04 | builtins.__import__ patch causes RecursionError via pandas internal imports during write_diagnostics |
| IMPL-54 | association_tests not in AssociationConfig; read with JSON fallback separately in _process() | 23-04 | Test names feed AssociationEngine not AssociationConfig; separate resolution with or-chain fallback |
| IMPL-55 | Baseline 3-df burden uses uniform weights [1,1,1]; sum/max use coast_weights | 24-01 | Matches R AllelicSeries: baseline tests per-category effect equality without ordering assumption |
| IMPL-56 | Allelic SKAT variance = aaf*(1-aaf) NOT 2*aaf*(1-aaf) | 24-01 | Matches AllelicSeries R source; factor of 2 absorbed by Q=score'score/2 convention |
| IMPL-57 | coast_burden_p_value in extra = Cauchy of 6 burden components (not 7-way omnibus) | 24-01 | Preserves COASTTest output contract; downstream expects standalone burden sub-summary |
| IMPL-58 | PurePythonCOASTTest fits null model lazily via PythonSKATBackend.fit_null_model() | 24-01 | Same pattern as PurePythonSKATTest; cohort-level singleton avoids repeated fitting |
| IMPL-59 | COAST backend swap placement: after SKAT swap, before unknown-name check | 24-02 | Follows IMPL-29; "coast" must resolve to correct class before validation |
| IMPL-60 | COAST auto mode probes both rpy2 AND AllelicSeries importr() | 24-02 | rpy2 presence alone doesn't guarantee AllelicSeries is installed; SKAT auto only checks rpy2 |
| IMPL-61 | coast_backend uses nullable=False in _get() (same as skat_backend) | 24-02 | CLI always writes the key with its default; non-nullable is correct for typed string fields |
| TEST-07 | Predictor/phenotype must use different seeds in _run_burden_test tests | 24-03 | Identical seeds produce identical arrays via numpy.random.default_rng; perfect correlation gives p=0.0 which is valid but misleading; different seeds ensure independent predictor/phenotype |
| TEST-08 | Engine _tests is dict[str, AssociationTest] not list; access via ["coast"] not [0] | 24-03 | AssociationEngine.__init__ builds {t.name: t} dict; list-style access causes KeyError |
| IMPL-62 | compute_acat_v uses single-pass loop: vectorized score/variance, per-variant loop for valid pairs | 25-02 | Clear filtering logic; avoids intermediate array allocation for edge case handling |
| IMPL-63 | mu_hat: np.ndarray \| None; all uses guarded by binary+not-None check inside compute_acat_v | 25-02 | Quantitative path always uses sigma2*I regardless of mu_hat; None is valid for quantitative |
| IMPL-64 | auto backend selector maps directly to PythonSKATBackend (no R probe) | 25-01 | R probe was expensive and error-prone; Python is always available |
| IMPL-65 | Unified in("python","auto") pattern in engine.py for backend swap | 25-01 | Single code path for auto and python avoids divergence if default changes later |
| IMPL-66 | DeprecationWarning only on RSKATTest/COASTTest __init__, not RSKATBackend | 25-01 | RSKATTest.check_dependencies() constructs RSKATBackend — would double-warn if both had warnings |
| IMPL-67 | Out-of-range saddlepoint fallback only in compute_pvalue (no proactive threshold) | 25-01 | IMPL-26 proactive threshold deferred; Phase 25 scope = out-of-range fallback only |
| IMPL-68 | acat_v_p=None in ALL early-return extra dicts in PurePythonSKATTest.run() | 25-02 | Key always present regardless of code path; prevents KeyError when engine reads res.extra |
| IMPL-69 | ACAT-V block in engine: AFTER test_pvals collection loop, BEFORE compute_acat_o() call | 25-02 | Critical ordering — if inserted after compute_acat_o(), ACAT-V is silently dropped from omnibus |
| IMPL-70 | association_workers uses nullable=False in _get() (same as coast_backend, skat_backend) | 27-02 | CLI always writes the key with its default (1); non-nullable is correct for typed int fields |
| PERF-01 | 128-node GL quadrature for SKAT-O omnibus integration; bounds [0,40] matching R upper=40 | 27-01 | 46x speedup (379ms -> 8ms per gene); 128 nodes sufficient for smooth SKAT-O integrands |
| PERF-02 | chi2(1) singularity at x=0 precludes direct GL accuracy validation; use exp(-x/20) instead | 27-01 | chi2(1) pdf diverges at x=0; GL achieves 1e-10 on smooth functions; SKAT-O integrand cancels singularity via (1-cdf) factor |
| ARCH-05 | parallel_safe=True on all Python-backend test classes (Fisher, LogisticBurden, LinearBurden, PurePythonSKAT) | 27-01 | Prerequisite for Plan 02 ProcessPoolExecutor dispatch; R-backend classes retain parallel_safe=False |
| PERF-03 | First gene runs sequentially before parallel dispatch to trigger lazy null model fitting | 27-03 | SKAT/COAST fit null models lazily on first run(); pre-fitting before pickle avoids redundant fitting in each worker |
| PERF-04 | Worker initializer sets OPENBLAS/MKL/OMP NUM_THREADS=1 | 27-03 | Prevents N_workers * BLAS_threads oversubscription causing CPU thrashing on multi-core machines |
| PERF-05 | use_parallel requires: n_workers != 1 AND all_parallel_safe AND len(sorted_data) > 1 | 27-03 | Single-gene panels skip parallel overhead; R-backend tests fall back to sequential with warning |
| PERF-06 | Worker under-provisioning guard: actual_workers = max(1, len(remaining)//2) when remaining < actual_workers * 2 | 27-03 | Prevents spawning more workers than useful for small remaining panels |
| TECH-01 | `parallel_safe: bool = False` as class attribute (not instance) in RSKATTest | 28-01 | Matches pattern of all other subclasses; engine uses getattr() so class-level is semantically correct |

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

### Roadmap Evolution

- Phase 25 added: Python Default Backends and Quick Wins (swap defaults, R deprecation, saddlepoint fallback, ACAT-V)
- Phase 26 added: Association Testing Documentation (guide, existing docs update, API stubs, changelog)
- Phase 27 added: Association Performance Optimizations (gene parallelization, Davies caching, single eigendecomposition)

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See archived REQUIREMENTS.md Future Requirements.
- ~~**RESEARCH-01** (before Phase 20): Validate whether parallel_safe=False on the stage is sufficient for rpy2 thread safety~~ — RESOLVED: parallel_safe=False + _assert_main_thread() guard confirmed sufficient (unit tests verify RuntimeError from worker threads)
- ~~**RESEARCH-02** (before Phase 21): Saddlepoint approximation algorithm for middle tier of Davies fallback chain~~ — RESOLVED: Kuonen Lugannani-Rice saddlepoint implemented in 21-01 from GENESIS variantSetTests.R reference
- ~~**COLUMN-01** (Phase 22): Rename `linear_burden_or` column to `linear_burden_beta`~~ — DONE (resolved via effect_column_names() in commit 48a6e68)

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-02-23T06:33:51Z – 2026-02-23T06:38:40Z
Stopped at: Completed 28-01-PLAN.md — RSKATTest parallel_safe attribute + lambda_gc verification gap
Resume file: None
Next: 28-02-PLAN.md (CLI args for skat_method and diagnostic thresholds)
Next: Plan Phase 28 (Tech Debt Cleanup). Then Phase 29 (Classic Pipeline Deprecation). After both complete, re-audit and complete milestone.
