# Research Summary: v0.15.0 Modular Rare Variant Association Framework

**Project:** variantcentrifuge — v0.15.0 milestone
**Domain:** Rare variant association testing for clinical genomics (nephrology, GCKD cohort)
**Researched:** 2026-02-19
**Confidence:** HIGH (primary literature + direct codebase reading)

---

## Executive Summary

This milestone adds a modular rare variant association framework to variantcentrifuge, upgrading
the existing Fisher's exact test in `gene_burden.py` to a multi-test engine supporting SKAT,
SKAT-O, ACAT-O, logistic/linear burden tests, covariate adjustment, and PCA integration. The GCKD
cohort (5,125 samples) is the primary target; all existing Mendelian analysis capabilities
(inheritance patterns, compound het, scoring) remain unchanged. The framework must be backward-
compatible: `--perform-gene-burden` users see zero behavioral change. The entire new `association/`
package is declared as an optional extra in pyproject.toml; the main package remains importable
without R or matplotlib.

The recommended approach is a six-phase build that starts with foundation abstractions and a
Fisher-equivalent refactor (bit-identical to current output), then adds covariate-adjusted burden
tests (using existing statsmodels), the R SKAT backend via rpy2 as a validated oracle, the pure-
Python SKAT fallback validated against R, ACAT-O combination and diagnostics, and finally PCA
computation as a pipeline stage. The dual-backend strategy (R via rpy2 + pure Python) is the key
architectural decision: it gives users without R a validated fallback while using the canonical R
SKAT package as the ground truth during development. Approximately 80% of the new statistical
functionality uses already-required dependencies (scipy, statsmodels, numpy).

The three most dangerous failure modes are: (1) Davies method silently returning wrong p-values
for extreme significance levels without raising any error; (2) using continuous-trait SKAT on
binary phenotypes without small-sample correction, which causes 5-fold type I error inflation with
imbalanced case/control ratios; and (3) rpy2 called from ThreadPoolExecutor worker threads, which
crashes the pipeline via segfault with no Python traceback. All three must be addressed as
architecture decisions in Phase 1, before any statistical code is written.

---

## Key Findings

### Stack Additions

The existing stack (pandas, numpy, scipy, statsmodels) handles approximately 80% of the new
functionality without new required dependencies. Only two optional dependencies are added:
`rpy2>=3.6.4` for the R SKAT backend and `matplotlib>=3.10` for diagnostic plots. Both are
declared under a new `[association]` extra — not required dependencies. The Davies method C
extension (`qfc.c` from the CompQuadForm R package) is bundled as source and compiled at runtime
via ctypes stdlib; if gcc is unavailable, Liu moment-matching (implemented in 20 lines of
existing scipy primitives) serves as fallback.

**Core additions:**

- `rpy2>=3.6.4` (optional, `[association]` extra) — R bridge for canonical SKAT package; 3.6.0 fixed R 4.3+ `Rcomplex` C-API bug; broad `except Exception` required at import, not just `ImportError` (R initialization failures raise `OSError` and `RRuntimeError`, not only `ImportError`)
- `ctypes` (stdlib) — wraps `qfc.c` for Davies method; no new dependency; runtime lazy compile with Liu fallback if gcc absent
- `scipy.stats.ncx2 / norm` (existing) — Liu 2009 moment-matching fallback; 20-line implementation from published algorithm; no new scipy requirement
- `scipy.linalg.eigh` (existing) — SKAT kernel eigenvalues; use `eigh` not `eig` for symmetric matrices (guaranteed real eigenvalues, faster)
- `statsmodels.api.Logit / OLS` (existing, 0.14.6) — logistic/linear burden tests; use `wald_test(scalar=True)` explicitly for forward compatibility with 0.15.x API change
- `matplotlib>=3.10` (optional) — QQ/Manhattan plots; lazy import only; `matplotlib.use("Agg")` must precede pyplot import for headless HPC nodes

**Critical Python version decision (must resolve before Phase 1):**

scipy 1.16 dropped Python 3.10. The current `requires-python = ">=3.10"` with no upper scipy
bound means Python 3.10 users get scipy 1.15.x (maintenance-only) while Python 3.11+ gets
scipy 1.17.x (active development, January 2026). Recommendation: bump `requires-python = ">=3.11"`.
HPC clusters have standardized on Python 3.11+ by 2025 and scipy 1.17.0 includes linalg
improvements relevant to SKAT kernel computation.

**Do NOT add:** momentchi2 (abandoned since 2021; Liu via scipy is equivalent), cffi (no
advantage over ctypes for one stable function), hatch-cython (build-time gcc requirement defeats
portability goal), pystatgen/limix (heavy BLAS/MKL deps, poor Windows support).

**pyproject.toml delta:**

```toml
[project.optional-dependencies]
association = ["rpy2>=3.6.4", "matplotlib>=3.10"]

[tool.hatch.build.targets.wheel]
artifacts = ["variantcentrifuge/data/*.c"]  # Bundle qfc.c source
```

### Feature Categories

**Table stakes — publishable credibility requires all of these:**

- Fisher's exact test (refactored, bit-identical output to current `--perform-gene-burden`)
- Logistic regression burden test with covariate adjustment (OR + 95% CI)
- SKAT and SKAT-O via dual backend (R via rpy2 + pure Python fallback)
- Beta(MAF; 1, 25) variant weights as default (the dominant standard across SKAT/SAIGE/REGENIE)
- FDR (BH) and Bonferroni multiple testing correction (already in gene_burden.py; carry over)
- Per-gene TSV output with p-values, effect sizes, variant counts, confidence intervals
- Lambda_GC genomic inflation factor diagnostic (expected in any analysis submission)
- QQ plot data TSV (observed vs expected -log10(p); plot optional, TSV mandatory)
- Binary trait support (case/control is the primary GCKD design)

**Standard output columns (seqMeta/SAIGE-GENE+ convention):**

```
GENE  n_variants  n_cases  n_controls  case_carriers  control_carriers
fisher_p  fisher_or  fisher_ci_lower  fisher_ci_upper
burden_p  burden_beta  burden_se  burden_or  burden_or_ci_lower  burden_or_ci_upper
skat_p  skat_stat
skat_o_p  skat_o_rho_opt  skat_o_stat
acat_o_p
corrected_fisher_p  corrected_skat_o_p  corrected_acat_o_p
```

**Differentiators — post-MVP, high scientific value:**

- ACAT-O omnibus test (analytically simple once SKAT/burden p-values exist; combine with Cauchy combination then single FDR)
- ACAT-V per-variant Cauchy combination
- Quantitative trait support (linear burden + linear SKAT for continuous phenotypes like eGFR)
- Functional variant weights (CADD-normalized, REVEL-based — scores already in pipeline columns from SnpSift extraction)
- PCA file integration (PLINK .eigenvec, AKT, generic TSV; 10 PCs default)
- JSON config mode for reproducible analysis SOP
- `p_method` column in output indicating Davies/saddlepoint/Liu (transparency on numerical precision)
- `SKAT_Null_Model_MomentAdjust()` for small-sample binary trait correction (R backend default)

**Defer to v2+ (explicit anti-features — do not build):**

- Meta-analysis (RAREMETAL, Meta-SAIGE) — different problem, out of scope for single-cohort tool
- Mixed model / GRM (SAIGE-GENE approach) — biobank-scale; PCs + kinship exclusion sufficient for GCKD (<10K samples)
- Phased haplotype tests — existing comp_het.py already handles compound heterozygous detection
- DeepRVAT neural network weights — requires massive training data, no nephrology model exists
- Conditional analysis — for known-signal follow-up; not needed for discovery
- Adaptive permutation p-values — high complexity; flag low-MAC genes in output as mitigation
- Allelic series test (COAST) — requires BMV/DMV/PTV classification; PolyPhen/SIFT available but v2+ scope
- REGENIE-style whole-genome split — designed for >100K samples; not needed for lab cohorts

**Small-sample handling (GCKD-specific):** With N_cases potentially 50-500, use
`SKAT_Null_Model_MomentAdjust()` (R backend) and `method="optimal.adj"` (SKAT-O) as defaults.
Flag genes where `case_carriers < 10` in output for cautious interpretation.

### Architecture Integration

The association framework slots into the existing stage-based pipeline as two new stages
(`AssociationAnalysisStage` and `PCAComputationStage`) without modifying any existing stage.
`GeneBurdenAnalysisStage` becomes a conditional shim that skips execution when
`--perform-association` is active. The `gene_burden.py` aggregation utilities (`_find_gt_columns`,
`_gt_to_dosage`, `_aggregate_gene_burden_from_columns`) are reused directly by the new engine
rather than duplicated. The R backend is isolated via a factory pattern (`get_skat_backend()`) so
the entire `association/` package remains importable without R.

**New package layout:**

```
variantcentrifuge/
  association/           # NEW — statistical engine
    engine.py            # AssociationEngine orchestrator (per-gene loop)
    config.py            # AssociationConfig dataclass
    results.py           # TestResult, AssociationResults container
    correction.py        # FDR/Bonferroni (refactored from gene_burden.py)
    covariates.py        # Covariate file loading + PCA merging + alignment validation
    pca.py               # PCA file parser + AKT/PLINK wrappers
    diagnostics.py       # Lambda_GC, QQ plot data
    tests/               # One module per statistical test
      base.py            # Abstract AssociationTest + TestResult
      fisher.py          # Thin wrapper around gene_burden logic
      logistic_burden.py # statsmodels.Logit; OR + 95% CI
      linear_burden.py   # statsmodels.OLS; beta + SE
      skat.py            # SKAT/SKAT-O dispatch to backend factory
      acat.py            # ACAT-V, ACAT-O Cauchy combination
    weights/             # Variant weighting schemes
      beta_weights.py    # Beta(MAF; 1, 25) via scipy.stats.beta.pdf
      uniform.py         # Uniform weights
    backends/            # SKAT-specific computation only
      r_backend.py       # rpy2 bridge to R SKAT (module-level rpy2 import OK here)
      python_backend.py  # Pure Python SKAT
      davies.py          # Davies ctypes(qfc.c) + saddlepoint + Liu fallback

  stages/
    analysis_stages.py   # MODIFY: add AssociationAnalysisStage; GeneBurdenAnalysisStage becomes shim
    processing_stages.py # MODIFY: add PCAComputationStage

  pipeline_core/
    context.py           # MODIFY: add association_results field

  data/
    qfc.c                # NEW: bundled Davies method C source
```

**Key integration points:**

- `AssociationAnalysisStage` hard-depends on `{"dataframe_loading", "sample_config_loading"}` and soft-depends on `{"custom_annotation", "pca_computation"}`; declared `parallel_safe = False` when R backend active
- `PCAComputationStage` lives in `processing_stages.py`, hard-depends on `{"variant_extraction"}`, stores eigenvectors in `context.stage_results["pca_computation"]`
- `PipelineContext` gains one new field: `association_results: Any | None = None`
- `get_skat_backend()` factory is called inside `AssociationEngine`, never at module import time
- `ExcelReportStage` gets a new "Association" sheet and optional "QQ Data" sheet
- Output: `{base}.association.tsv` + optional `{base}.diagnostics/` directory
- Genotype matrix is NEVER stored in PipelineContext — extracted per-gene inside the engine loop and discarded immediately (memory safety: 5K samples × 50K variants = 1.6 GB if stored)
- `gene_burden.py` is not removed or restructured during this milestone; `perform_gene_burden_analysis()` remains the public API for backward compat

**Backward compatibility guarantee:** `--perform-gene-burden` alone produces zero behavioral change. `GeneBurdenAnalysisStage` shim exits early only when `perform_association=True`.

### Critical Pitfalls (Top 5)

**1. R thread safety — crashes parallel pipeline (CRITICAL)**

R is not thread-safe. Calling rpy2 from `ThreadPoolExecutor` worker threads (the existing parallel
execution path in `pipeline_core/runner.py`) causes segfaults or `ValueError: signal only works
in main thread`. No Python traceback is produced — the crash is silent from the pipeline's
perspective. Prevention: declare any stage using the R backend as `parallel_safe = False`. R
backend calls must originate from the main thread only. Use `ProcessPoolExecutor` (separate
interpreter per process) only for the pure Python backend. This constraint must be locked in as
an architecture invariant before any SKAT code is written.

**2. Binary trait type I error inflation (CRITICAL)**

Standard SKAT (continuous-trait formulation) on binary phenotypes with case/control imbalance
produces up to 5-fold inflated type I error. The asymptotic null distribution is inaccurate for
sparse genotype data with imbalanced samples. Prevention: always use `SKATBinary` (R package) for
binary traits, never `SKAT`. Implement moment-adjustment (`SKAT_Null_Model_MomentAdjust`) as the
default for small samples. Default to `method="optimal.adj"` for SKAT-O. Warn when `n_cases < 200`
or `case:control ratio > 1:20`. Gate on trait_type at null model construction time in Phase 1
— this decision cannot be retrofitted later.

**3. Davies method silent numerical failure (CRITICAL)**

The Davies method returns plausible-looking wrong p-values (not an error, no exception) when the
true p-value falls below 10^-6 with default accuracy settings (acc=1e-6, lim=10^4). It may also
return exactly `0.0` due to floating-point underflow, which is uninformative. Prevention: use
corrected defaults from Liu et al. 2016 (acc=1e-9, lim=10^6). Implement hybrid fallback chain:
Davies (primary) → saddlepoint approximation (secondary) → Liu (last resort only). Never pass
`p=0.0` downstream; log as `"p < precision_floor"`. Record computation method in `p_method` output
column. Liu is anti-conservative for p < 10^-4 — it must not be treated as a routine fallback.

**4. rpy2 memory leaks in per-gene loops (HIGH)**

R objects created during each per-gene SKAT call accumulate in R's heap because Python's GC cannot
see R-managed memory. For 5K samples × thousands of genes, peak memory eventually exceeds HPC node
limits. Prevention: convert R results to Python primitives immediately (`float(result[0])`, not
`result[0]`), `del` R object references after each gene, call `rpy2.robjects.r('gc()')` every
100 genes. Profile memory on a 100-gene run before declaring the R backend production-ready.

**5. Sample/covariate mismatch — silent analysis corruption (MEDIUM but extremely dangerous)**

Genotype data row order is determined by the VCF header. Covariate file row order is arbitrary.
If alignment is incorrect, each sample's covariates are assigned to the wrong genotype row — the
analysis runs without errors and produces plausible but completely wrong results. Prevention:
always reindex the covariate DataFrame to `context.vcf_samples` order explicitly. Assert no NaN
after reindex (any VCF sample absent from the covariate file is an error, not a warning). Test
with shuffled covariate file row order and assert identical results.

**Additional high-severity pitfalls (tracked in PITFALLS.md):**

- Eigenvalue instability from near-singular kernel matrices: use `scipy.linalg.eigh`, threshold `eigenvalues = np.maximum(eigenvalues, 0)`, skip SKAT for genes with `matrix_rank < 2`
- Genotype matrix construction: existing `_gt_to_dosage` must be audited for `1/2` multi-allelic hets and `./. ` missing genotypes (mean-impute to 2×MAF, not 0)
- SKAT-O rho grid correctness: the covariance of rho-specific statistics is non-trivial; implement SKAT before SKAT-O; validate Python SKAT-O via extensive null simulation before deploying
- Multiple testing strategy: apply ACAT-O Cauchy combination within genes, then single FDR correction across genes — applying FDR separately per test is statistically wrong
- PC overcorrection: >10 PCs can inflate type I error; default to 10; warn at >20; validate lambda_GC under permuted null

---

## Implications for Roadmap

The dependency graph from FEATURES.md and ARCHITECTURE.md dictates a clear 6-phase build order.
Each phase is independently testable and delivers scientific value before the next begins.

### Phase 1: Foundation — Core Abstractions + Fisher Refactor

**Rationale:** All subsequent phases depend on `AssociationAnalysisStage` being registered in the
pipeline and the `AssociationEngine` orchestration pattern being established. The Fisher refactor
provides a validatable baseline with no new math — the validation gate is bit-identical output.
Critical architecture decisions (R thread safety, binary trait gating, correction strategy) must
be locked in here before any statistical code is written.

**Delivers:**
- `association/` package skeleton (`__init__.py`, `config.py`, `results.py`, `tests/base.py`)
- `association/tests/fisher.py` wrapping existing `gene_burden.perform_gene_burden_analysis`
- `AssociationEngine` with single-test execution loop
- `AssociationAnalysisStage` (Fisher-only) registered in `stage_registry.py`
- `GeneBurdenAnalysisStage` shim check (`perform_association` guard)
- `PipelineContext.association_results` field added to `context.py`
- `ExcelReportStage` "Association" sheet extension
- CLI args: `--perform-association`, `--association-tests`
- `association/correction.py` (FDR/Bonferroni, refactored from gene_burden.py, re-exported there)

**Validation gate:** `--perform-association --association-tests fisher` produces bit-identical
output to `--perform-gene-burden` across all existing integration tests. All existing
`--perform-gene-burden` integration tests pass unchanged.

**Addresses pitfalls:** R thread safety architecture decision, multiple testing strategy design,
sample mismatch validation pattern established in covariates.py foundation.

**Research flag:** Standard patterns — no additional research needed. Stage integration and
factory patterns are well-established in the existing codebase.

### Phase 2: Covariate System + Logistic/Linear Burden

**Rationale:** Uses only existing dependencies (statsmodels is already required). Provides
immediate scientific value (covariate-adjusted burden tests) before the complex SKAT math.
Establishes the covariate matrix assembly pattern that all later phases reuse. The genotype matrix
builder audit for multi-allelic/missing genotypes belongs here (it affects burden tests as well
as SKAT).

**Delivers:**
- `association/covariates.py` — covariate file loading, sample ID alignment to `vcf_samples`, missingness validation, condition-number multicollinearity check
- `association/weights/` — `beta_weights.py` (Beta(MAF; 1, 25) via `scipy.stats.beta.pdf(maf, 1, 25)`), `uniform.py`
- `association/tests/logistic_burden.py` — `statsmodels.Logit` with Wald and LRT tests, OR + 95% CI; `wald_test(scalar=True)` for API forward-compat
- `association/tests/linear_burden.py` — `statsmodels.OLS` with beta + SE for quantitative traits
- Genotype matrix builder audited and tested for `1/2`, `.|.`, `0|1`, `./.` edge cases; mean imputation for missing
- CLI args: `--covariate-file`, `--covariates`, `--variant-weights`, `--trait-type`

**Avoids:** Covariate multicollinearity (condition number check), sample mismatch (explicit
reindex with assertion), PC overcorrection (document 10-PC default before PCA stage exists).

**Research flag:** Standard patterns — statsmodels Logit/OLS API well-documented; wald_test
scalar note already captured in STACK.md.

### Phase 3: R SKAT Backend (Oracle)

**Rationale:** The R SKAT package (leelabsg/SKAT on CRAN) is the peer-reviewed gold standard
with 10+ years of validation. Build and validate it before writing the Python backend — this
order allows the R backend to serve as the correctness oracle during Phase 4 development.

**Delivers:**
- `association/backends/base.py` — `SKATBackend` ABC, `NullModel` container
- `association/backends/__init__.py` — `get_skat_backend()` factory (lazy, never at module import time)
- `association/backends/r_backend.py` — `RSKATBackend` via rpy2; `SKATBinary` for binary traits by default; `SKAT_Null_Model_MomentAdjust` for small-sample binary adjustment; explicit R GC every 100 genes; `parallel_safe = False` on the stage
- `association/tests/skat.py` — dispatch to backend for both SKAT and SKAT-O
- Synthetic cohort test fixtures (100 cases/controls, 50 genes, 5 signal genes with injected burden)
- CLI arg: `--skat-backend auto|r|python`

**Avoids:** R thread safety (enforce `parallel_safe = False`, main-thread-only R calls; no
ThreadPoolExecutor for R-using stages), binary trait type I error (SKATBinary + moment adjustment
by default), R memory leaks (explicit del + periodic gc()), R_HOME detection failure (informative
error message + automatic fallback to Python backend).

**Research flag:** Needs phase research for rpy2 thread safety solution before implementation
commits. Specifically: whether `parallel_safe=False` on the stage is sufficient given the existing
runner architecture, or whether a dedicated R worker process with queue-based serialization is
required.

### Phase 4: Pure Python SKAT Backend

**Rationale:** Requires the R backend as a correctness reference. The iterative development loop
(implement → compare against R on same inputs → fix discrepancies) needs both backends running
simultaneously. SKAT before SKAT-O — validate the simpler test first, then add the rho grid.

**Delivers:**
- `association/backends/davies.py` — Liu fallback first (pure scipy, always available), then ctypes qfc.c lazy compile; fallback order: Davies → saddlepoint → Liu; `p_method` metadata in every result
- `association/backends/python_backend.py` — `PythonSKATBackend`; SKAT implementation first; SKAT-O only after SKAT validates against R across 50+ genes
- `variantcentrifuge/data/qfc.c` — bundled qfc.c source from CompQuadForm; `pyproject.toml` artifacts field updated
- Validation test suite: Python vs R p-values within 10% relative difference on log10(p) scale for p > 0.001; larger differences expected and documented for p < 0.001
- Eigenvalue stability: `scipy.linalg.eigh`, threshold `np.maximum(eigenvalues, 0)`, skip if `matrix_rank < 2`

**Avoids:** Davies silent failure (acc=1e-9, lim=10^6; saddlepoint fallback; p=0.0 logged as
floor, not passed downstream), SKAT-O rho correctness (SKAT first; null simulation before SKAT-O
deployment), Liu anti-conservation (Liu is last resort with `p_method="liu"` flag), ctypes
compilation failures (Liu always available; Davies is enhancement).

**Research flag:** Needs phase research for saddlepoint approximation algorithm (the middle tier
of the Davies fallback chain). Not in scipy stdlib; requires derivation from Davies 1980 or
reference implementation from fastSKAT literature before Phase 4 starts.

### Phase 5: ACAT-O + Diagnostics

**Rationale:** Analytically simple once SKAT and burden p-values exist (ACAT-O is a Cauchy
combination formula, a few lines of scipy). Diagnostics require all p-values to be collected
across all genes, making this naturally the last computation step.

**Delivers:**
- `association/tests/acat.py` — ACAT-V (per-variant Cauchy via `scipy.stats.cauchy.sf`) and ACAT-O (combine SKAT + burden + ACAT-V p-values; one Cauchy combination per gene)
- `association/diagnostics.py` — lambda_GC per test (`median(chi2) / 0.4549`), QQ data TSV, optional matplotlib QQ plot PNG with lazy headless import
- `ExcelReportStage` "QQ Data" sheet addition
- Single FDR correction applied to ACAT-O p-values across genes (not separately per test)
- Per-test raw p-values reported alongside combined ACAT-O for full transparency
- CLI arg: `--diagnostics-output` (path for diagnostics/ directory)

**Avoids:** Multiple testing correction on wrong set (ACAT-O combination then single FDR),
ACAT-O with non-uniform inputs (validate all inputs are analytical p-values, not permutation-
based), Liu anti-conservation (ACAT-O inputs come from Davies/saddlepoint path).

**Research flag:** Standard patterns — ACAT methodology is well-published (Liu et al. 2019,
PMC6407498). Cauchy combination formula is straightforward.

### Phase 6: PCA Integration

**Rationale:** Optional and architecturally independent. Pre-computed PCA files already work from
Phase 2 onward via `--covariate-file`. The `PCAComputationStage` is purely additive. Deferring it
allows Phases 1-5 to validate the statistical core without the external tool dependency.

**Delivers:**
- `association/pca.py` — PLINK `.eigenvec` parser, AKT output parser, generic TSV parser; sample ID alignment validation
- `PCAComputationStage` in `processing_stages.py` — wraps akt and plink2 via subprocess; stores eigenvectors in `context.stage_results["pca_computation"]`; registered in stage_registry under "processing" category
- `covariates.py` updated to merge eigenvectors from `stage_results["pca_computation"]` alongside covariate file columns
- Lambda_GC validation with and without PCA to detect overcorrection
- CLI args: `--pca-file`, `--pca-tool akt|plink`, `--pca-components` (default 10; warn if >20)
- Documentation: PCA must be computed on common variants (MAF > 1%), not on the rare variants being tested

**Avoids:** PC overcorrection (default 10; warn at >20; lambda diagnostic), PCA computed on rare
variants (explicit documentation and optional validation check).

**Research flag:** Standard patterns if loading pre-computed files. Needs phase research if
automating AKT/PLINK2 invocations — tool versions, exact command-line flags, and output format
consistency across versions are not fully verified.

### Phase Ordering Rationale

- Phases 1-2 use zero new dependencies (statsmodels/scipy already required) — validate the new architecture before adding optional extras
- Phase 3 before Phase 4 — R backend is the oracle; correctness of the Python backend cannot be established without R comparison
- Phase 5 after Phases 3 and 4 — ACAT-O requires p-values from both SKAT and burden tests to exist
- Phase 6 last — PCA is optional and architecturally independent; pre-computed file path available from Phase 2 via `--covariate-file`
- Binary trait gating (SKATBinary vs SKAT) is a Phase 1 architecture decision that cannot be retrofitted — it affects the null model interface used by all later phases

### Research Flags

**Needs deeper research before implementation:**

- **Phase 3:** rpy2 thread safety solution in the existing PipelineRunner — whether `parallel_safe=False` on the stage is sufficient or whether a dedicated R worker process is required
- **Phase 4:** Saddlepoint approximation implementation — the mathematical derivation for the intermediate Davies fallback tier is not in any stdlib and needs a verified reference implementation

**Standard patterns (skip research-phase):**

- **Phase 1:** Stage integration, factory pattern, abstract base classes — well-established patterns in the existing codebase (see ARCHITECTURE.md for precise integration points)
- **Phase 2:** statsmodels Logit/OLS — well-documented; API notes already captured
- **Phase 5:** ACAT-O Cauchy combination — simple formula from published paper (PMC6407498)
- **Phase 6 (file loading path):** PLINK .eigenvec format is fixed and well-documented

---

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack additions | HIGH | All versions verified on PyPI 2026-02-19; ctypes vs cffi decision technically clear; scipy version constraint documented from release notes |
| Feature requirements (table stakes) | HIGH | From original SKAT/SKAT-O papers (PMC3135811, PMC3415556) and major tool documentation (SAIGE, REGENIE) |
| Feature requirements (differentiators) | HIGH | ACAT from PMC6407498; COAST from PMC10432147; Beta(1,25) confirmed in multiple tool docs |
| Architecture integration | HIGH | Based on direct codebase reading of context.py, analysis_stages.py, gene_burden.py, stage_registry.py |
| Critical pitfalls (Davies, binary traits) | HIGH | Primary literature (Liu et al. 2016, Lee et al. 2012) and SKAT package official docs |
| rpy2 thread safety | HIGH | Multiple independent sources; documented in rpy2 code and confirmed by external reports |
| SKAT-O rho correctness implementation details | MEDIUM | SKAT-O paper is definitive for the math; specific Python re-implementation bugs are speculative until implemented |
| PCA tooling (akt/plink exact flags) | MEDIUM | Tool versions and output format consistency not verified against current releases |
| Madsen-Browning current adoption | LOW | Single WebSearch source; deferred appropriately to anti-features |

**Overall confidence:** HIGH for architecture decisions and critical statistical pitfalls; MEDIUM
for implementation details that only surface during Phase 4 (SKAT-O rho, saddlepoint derivation).

### Gaps to Address During Implementation

- **Saddlepoint approximation algorithm:** The middle tier of the Davies fallback chain. Not in scipy stdlib. Derivation from Davies 1980 or reference from fastSKAT (PMC4375394) needed before Phase 4 starts.
- **rpy2 thread safety solution:** Queue-based serialization is the documented solution but adds complexity. Validate whether `parallel_safe=False` on the stage is sufficient in the existing PipelineRunner, or whether a dedicated R worker process is required.
- **GCKD cryptic relatedness:** The GCKD cohort (kidney disease study) may include related individuals. Standard SKAT assumes independence. If IBD analysis reveals kinship > 0.05 for any pairs, standard SKAT is invalid and the analysis plan requires mixed-model SKAT (GENESIS package) or sample pruning. This is a study design question that should be assessed before Phase 3.
- **Python version bump validation:** The `requires-python = ">=3.11"` change must be confirmed acceptable before Phase 1. Active users on Python 3.10 need a migration communication plan.
- **qfc.c compilation on HPC job nodes:** Runtime compilation via gcc works on development machines; HPC job nodes often block subprocess execution. If the GCKD cluster does not permit gcc during job execution, qfc.c must be pre-compiled into the package or the saddlepoint fallback must be treated as the production Davies path (not an exception case).

---

## Sources

### Primary (HIGH confidence)

- [SKAT original paper — Wu et al. 2011 (PMC3135811)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3135811/) — SKAT formulation, Beta(1,25) weights, kernel definition
- [SKAT-O paper — Lee et al. 2012 (PMC3415556)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/) — SKAT-O rho grid, small-sample adjustment for binary traits, SKATBinary requirement
- [Liu et al. 2016 — efficient SKAT p-value (PMC4761292)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4761292/) — hybrid Davies/saddlepoint approach, corrected accuracy defaults (acc=1e-9, lim=10^6)
- [ACAT paper — Liu et al. 2019 (PMC6407498)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407498/) — Cauchy combination formula, ACAT-O, multiple testing strategy
- [COAST allelic series — PMC10432147](https://pmc.ncbi.nlm.nih.gov/articles/PMC10432147/) — BMV/DMV/PTV classification, COAST test; correctly deferred to v2+
- [PC stratification in rare variant tests — PMC6283567](https://pmc.ncbi.nlm.nih.gov/articles/PMC6283567/) — empirical evidence for >10 PC overcorrection
- [Controlling stratification — PMC8463695](https://pmc.ncbi.nlm.nih.gov/articles/PMC8463695/) — PC selection guidance
- [rpy2 official docs — version 3.6.4](https://rpy2.github.io/doc/latest/html/introduction.html) — importr pattern, numpy2ri, memory management
- [rpy2 changelog (Rcomplex fix in 3.6.0)](https://rpy2.github.io/doc/latest/html/changes.html) — R 4.3+ compatibility requirement
- [SKAT R package documentation (RDocumentation)](https://www.rdocumentation.org/packages/SKAT/versions/2.0.1) — SKATBinary, SKAT_Null_Model_MomentAdjust, method="optimal.adj"
- [SAIGE-GENE+ documentation](https://saigegit.github.io/SAIGE-doc/) — output format conventions
- [REGENIE gene-based test options](https://rgcgithub.github.io/regenie/options/) — supported tests, ACATO as peer to SKATO
- [CompQuadForm CRAN source — qfc.c function signature](https://rdrr.io/cran/CompQuadForm/src/R/davies.R) — verified 2026-02-19

### Secondary (MEDIUM confidence)

- [RAVA-FIRST CADD weighting — PMC9518893](https://pmc.ncbi.nlm.nih.gov/articles/PMC9518893/) — functional weight adoption
- [DeepRVAT — Nature Genetics 2024](https://www.nature.com/articles/s41588-024-01919-z) — neural network weights, correctly excluded as anti-feature
- [fastSKAT eigenvalue approximation — PMC4375394](https://pmc.ncbi.nlm.nih.gov/articles/PMC4375394/) — saddlepoint reference for Phase 4
- [rpy2 R 4.4 CentOS 7 issue (GitHub #1107)](https://github.com/rpy2/rpy2/issues/1107) — platform-specific build issue scope (CentOS 7 only)
- [rpy2 thread safety confirmed (Streamlit issue)](https://github.com/streamlit/streamlit/issues/2618) — external corroboration
- [rpy2 memory management docs](https://rpy2.github.io/doc/v3.2.x/html/rinterface-memorymanagement.html) — gc() cadence guidance
- [scipy PyPI (1.17.0, January 2026)](https://pypi.org/project/SciPy/) — version verified
- [statsmodels PyPI (0.14.6, December 2025)](https://pypi.org/project/statsmodels/) — version and wald_test API note verified
- [matplotlib PyPI (3.10.8, December 2025)](https://pypi.org/project/matplotlib/) — version verified
- [seqMeta output conventions](https://github.com/genepi-freiburg/seqmeta) — column naming standard
- [Multi-allelic variant association (PMC7288273)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7288273/) — genotype matrix construction guidance
- [Population stratification — PMC6283567](https://pmc.ncbi.nlm.nih.gov/articles/PMC6283567/) — PC count empirical evidence

### Tertiary (LOW confidence)

- Madsen-Browning / WSS current adoption rate — single WebSearch source; deferred to anti-features
- ACAT-O replacing SKAT-O as primary — no evidence found; both are peer options in REGENIE; ACAT-O is additional test, not a replacement

---

*Research completed: 2026-02-19*
*Synthesized from: STACK.md, FEATURES.md, ARCHITECTURE.md, PITFALLS.md (4 research streams)*
*Ready for roadmap: yes*
