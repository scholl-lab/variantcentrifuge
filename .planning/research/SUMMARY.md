# Research Summary: v0.16.0 Association Hardening & Multi-Cohort Features

**Project:** variantcentrifuge — v0.16.0 milestone
**Domain:** Rare variant association framework hardening and multi-cohort features
**Researched:** 2026-02-23
**Confidence:** HIGH

---

## Executive Summary

v0.16.0 is a hardening and capability-expansion milestone for the rare variant association framework introduced in v0.15.0. The v0.15.0 foundation is solid — SKAT, SKAT-O, ACAT-O, COAST, covariate adjustment, dual R/Python backends, and FDR correction all exist and work at the GCKD cohort scale (5,125 samples). v0.16.0 adds seven targeted improvements that address real bugs (COAST p=None for most genes due to a genotype matrix reconstruction failure and a hard 3-category gate), missing infrastructure (PCAComputationStage exists but is not wired into the pipeline), and methodological gaps needed for multi-cohort use (weighted FDR for gene-level priors, BED-based region restriction, sparse genotype matrix support, case-confidence sample weighting, and SKAT-O eigendecomposition optimization). The dominant theme is: the code is already there, the wiring or correctness is the gap.

The critical finding from architecture research is that ALL eight v0.16.0 features have well-defined, narrow integration points. No feature requires new Python package dependencies — the installed stack (scipy 1.14.1, statsmodels 0.14.4, numpy 2.2.6) supports every proposed API change. SKAT-O eigendecomposition optimization may be partially a non-issue: the existing per-rho eigenvalue structure is already mathematically sound, and the claimed "single decomposition" opportunity requires clarification before any implementation begins. Three features are clearly high-priority bug fixes (COAST genotype matrix reconstruction, COAST 3-category gate, PCAComputationStage wiring) that block scientific credibility. Four features are genuine capability additions (weighted BH, case-confidence weights, sparse matrices, region restriction).

The key risks are statistical correctness risks, not implementation complexity. Weighted BH without weight renormalization silently violates FDR control with no error raised. Case-confidence weights applied to residuals rather than the null model inflate type I error. BED region files with chromosome prefix mismatches produce 0-variant output with no warning. These are "looks right, gives wrong science" failure modes that require permutation testing and input validation guards — not just code review. The recommended phase order follows the dependency structure: COAST fix first (unblocks validation of other features), then parallel standalone phases (FDR weighting, region restriction, PCA wiring), then dependent features (case-confidence weights, sparse matrices), then deferred optimizations.

---

## Key Findings

### Stack Additions

Zero new Python package dependencies are required for any v0.16.0 feature. The installed stack is already sufficient:

**Core technologies (unchanged from v0.15.0):**
- `scipy 1.14.1` — `scipy.sparse.csr_array` available since scipy 1.7; CSR sparse format recommended for genotype matrices (24x memory reduction at 2% fill, expected for rare variants with MAF < 1%)
- `statsmodels 0.14.4` — `sm.GLM(var_weights=...)` supports case-confidence weights (documented and verified in 0.15.0 dev docs)
- `numpy 2.2.6` — weighted Benjamini-Hochberg implementable in ~12 lines; neither statsmodels `multipletests` nor `scipy.stats.false_discovery_control` has a `weights` parameter
- `bcftools` (external, already required) — `--regions-file` flag handles BED restriction natively
- `akt` / `plink2` (external, optional) — output already parsed by `pca.py`; new `PCAComputationStage` wiring only

**Explicitly rejected additions:**
- IHW Python package — no mature PyPI implementation exists as of February 2026; defer to rpy2 bridge on demand
- `scikit-allel` — `scipy.sparse.csr_array` already covers the sparse genotype matrix need
- `polars` — Pandas deeply embedded; mixing frames is a maintenance burden
- Any version bumps — no new feature requires a newer version than currently installed

**Critical version constraint to preserve:** Do NOT bump scipy to 1.16+ — it requires Python 3.11+, which conflicts with `requires-python = ">=3.10"`.

See `.planning/research/STACK.md` for detailed API verification, breakeven analysis, and version compatibility notes.

### Feature Categories

**Must ship — table stakes (fix real bugs, enable primary use cases):**
- COAST genotype matrix reconstruction fix — Path B fallback (`analysis_stages.py` lines 2473-2486) reaches into `context.variants_df` which may have dropped per-sample GT columns; add Path C that reconstructs from packed GT string via `create_sample_columns_from_gt_vectorized()`
- COAST partial-category fallback — hard 3-category gate at lines 551-581 of `allelic_series.py` returns p=None for any gene missing one category; this describes most real-world genes
- COAST configurable classification — support CADD/REVEL thresholds when SIFT/PolyPhen unavailable; configurable via a new `scoring/coast_classification/` JSON scoring model
- Weighted BH (static, Genovese 2006) — ~12 lines of numpy; add `apply_weighted_correction()` to `correction.py`
- Per-gene FDR weight file loader — new `--fdr-weight-file` CLI flag; TSV format; normalize to mean=1 internally
- `--restrict-regions` BED prefilter — bcftools-native; extend `GeneBedCreationStage` with bedtools intersect

**Defer to post-MVP — differentiators (harder, or need validation first):**
- IHW integration (R/rpy2 bridge, cross-validation complexity, thread-safety constraints)
- Weighted SKAT per-sample phenotype weights (research-grade; no published reference implementation; needs permutation validation)
- Sparse genotype matrices (memory nicety at GCKD 5K scale; becomes essential at 50K+ scale)
- PCAComputationStage wiring (documentation gap, not a bug; users can run PLINK2 manually)
- SKAT-O single eigendecomposition (clarify feature intent first; may already be implemented)

**Anti-features (explicitly do not build):**
- Mixed model / GRM (SAIGE-GENE approach) — overkill for GCKD cohort; PCA + kinship exclusion is sufficient
- Meta-analysis — different problem domain (cross-cohort combining)
- Mandatory biological priors for FDR — default must remain flat BH to avoid blocking basic usage
- Deep learning variant weights (DeepRVAT) — requires massive phenotype-specific training data

See `.planning/research/FEATURES.md` for the complete feature dependency graph, v0.16.0 MVP recommendation, and anti-feature rationale.

### Architecture Integration

All seven v0.16.0 features integrate into the existing stage-based pipeline at narrow, well-defined points. The association subsystem uses a clean layering: `AssociationConfig` dataclass (all parameters) → `AssociationEngine` (per-gene loop, FDR, ACAT-O) → `AssociationTest` ABC subclasses → `PythonSKATBackend` (SKAT/SKAT-O/burden). New features fit into this hierarchy without structural changes.

**Key integration points by feature:**
1. COAST fix — `stages/analysis_stages.py` only: add Path C genotype matrix fallback (reconstruct from packed GT string), split EFFECT column on `&` before `isin()` matching
2. Weighted BH — `association/correction.py`: new `apply_weighted_correction()` leaf function; `association/engine.py`: routing when `config.fdr_weights` set
3. Region restriction — `stages/processing_stages.py` `GeneBedCreationStage._process()`: add `bedtools intersect` call when `restrict_regions_bed` configured; `cli.py`: `--restrict-regions` flag
4. PCA wiring — `pipeline.py`: 4-line import + conditional append; `stages/processing_stages.py` `PCAComputationStage._process()`: use prefiltered VCF fallback chain
5. Case-confidence weights — new file `association/sample_weights.py`; `association/backends/python_backend.py`: `fit_null_model()` gains `sample_weights` param, `_compute_eigenvalues_filtered()` adjusts phi
6. Sparse matrices — `association/genotype_matrix.py`: optional `return_sparse: bool`; `python_backend.py`: sparse-aware matmuls; `association/base.py`: `use_sparse_geno` flag
7. COAST scoring model — new JSON config files only: `scoring/coast_classification/formula_config.json` and `variable_assignment_config.json`; no Python changes

**Architecture patterns to follow (established in v0.15.0):**
- Config field expansion via `AssociationConfig` dataclass with backward-compatible defaults
- Gene data dict as feature bag: store new data in augmentation loop (lines 2634-2751), read in `run()` with `.get()` and sensible fallback
- Skip guards returning `p_value=None`: every test method begins with skip checks; engine handles `None` gracefully
- Leaf module isolation: `correction.py`, `covariates.py`, `sample_weights.py` import only stdlib/numpy/pandas/statsmodels
- Lazy null model caching: fit once per cohort in `self._null_model`, reuse per gene

**New files required:** `association/sample_weights.py`, `scoring/coast_classification/formula_config.json`, `scoring/coast_classification/variable_assignment_config.json`.

See `.planning/research/ARCHITECTURE.md` for exact line references, code samples, data flow diagrams, and the build order analysis.

### Critical Pitfalls (Top 5)

1. **Weighted BH weights not renormalized to mean=1** — The Genovese (2006) FDR guarantee holds only when `w.mean() == 1`. User weight files contain raw scores (pLI, GWAS odds ratios), not normalized weights. This is a silent violation: q-values look plausible, no error is raised, but the true FDR exceeds nominal alpha in proportion to how far the mean weight deviates from 1. Prevention: normalize at load time (`w = w / w.mean()`), assert `abs(w.mean() - 1.0) < 1e-9`, add permutation test validating null FDR.

2. **Case-confidence SKAT weights applied to residuals, not null model** — Multiplying `residuals * sample_weights` post-hoc changes the score statistic but not the eigenvalue null distribution, inflating type I error in proportion to weight variance. The null model must be refit with `sm.GLM.fit(var_weights=sample_weights)` AND phi must be adjusted to `sqrt(mu*(1-mu) * sample_weights)` in `_compute_eigenvalues_filtered()`. Validate with permutation: lambda_GC on permuted phenotype with weights held fixed must be ~1.0.

3. **SKAT-O "single eigendecomposition" is either wrong or already implemented** — Per-rho eigenvalues ARE required; the rho-dependent kernel K_rho has genuinely different eigenspectra for each rho value. If "single decomposition" means reusing rho=0 eigenvalues for all rho, it is mathematically incorrect and will bias rho selection toward rho=0. The current code already precomputes the base `A = z1_half.T @ z1_half` once and derives per-rho k_sym analytically from A. Clarify intent before any implementation begins.

4. **BED region restriction silent chromosome name failure** — `bcftools view -R <bed>` requires exact chr name matching (chr1 vs 1). A mismatch produces 0-variant output with no error or warning. Prevention: validate chromosome names by comparing the first BED record to VCF header contigs before the bcftools call; add post-filter guard emitting WARNING when 0 variants pass and `--restrict-regions` was specified.

5. **COAST multi-transcript SnpEff `&`-concatenated effect strings** — `classify_variants()` uses exact `isin()` matching on the EFFECT column. SnpEff produces compound consequences like `"stop_gained&splice_region_variant"` — this string is NOT in `PTV_EFFECTS`, so a PTV variant is classified as code 0 (excluded from COAST). Prevention: split effect strings on `&` before matching; add unit tests with compound SnpEff effect strings.

See `.planning/research/PITFALLS.md` for 12 total pitfalls including: missing genes in prior file excluded from FDR (Pitfall 2), PCAComputationStage registry/pipeline coupling (Pitfall 7), dead code registry/import coupling (Pitfall 8), sparse matrix eigendecomposition requiring dense input (Pitfall 9), COAST weight vector length mismatch (Pitfall 10), IHW R thread restriction (Pitfall 11), and ultra-rare variant p=1.0 instead of skip (Pitfall 12).

---

## Implications for Roadmap

Eight v0.16.0 features map to a 6-phase structure ordered by dependencies and risk. Phases 1 and 5 address critical bugs. Phases 2-4 are standalone capability additions that can be developed in parallel. Phase 6 addresses the highest-risk implementation. Deferred optimizations (single eigendecomposition, dead code cleanup) are explicitly out of scope.

### Phase 1: COAST Fix — Genotype Matrix Reconstruction and Partial-Category Fallback

**Rationale:** The COAST p=None bug blocks scientific use of the most differentiated feature in the framework. It also blocks validation of any feature that touches the COAST path (case-confidence weights in Phase 5, COAST scoring model validation). Fix this first to establish a working baseline. The root cause is fully confirmed from direct code reading; the fix is contained to `analysis_stages.py`.

**Delivers:** Working COAST p-values for real-world cohorts where not all genes have all three variant categories; Path C genotype matrix fallback using `create_sample_columns_from_gt_vectorized()`; SnpEff `&`-concatenated effect string handling; configurable BMV/DMV/PTV classification thresholds via JSON config.

**Addresses:** FEATURES.md COAST configurable classification (table stakes); FEATURES.md partial-category fallback (table stakes); PITFALLS 5 (multi-transcript effect strings) and 10 (weight vector length mismatch).

**Avoids:** ARCHITECTURE anti-pattern 1 (reaching into `variants_df` as genotype source).

**Research flag:** Standard patterns — root cause confirmed by direct code reading at exact line numbers. No deeper research needed.

### Phase 2: Gene-Level FDR Weighting (Weighted BH)

**Rationale:** Standalone; no dependency on Phase 1. Short implementation (~20 lines in `correction.py`, new CLI flag, weight file loader). High scientific value: up to 2x more discoveries when good biological priors are available (Huang et al. 2022). Pure numpy — no new dependencies.

**Delivers:** `apply_weighted_correction()` in `correction.py`; `--fdr-weight-file` CLI flag; weight normalization with `mean=1` enforcement and assertion; missing-gene imputation to weight=1.0 with renormalization; permutation validation test.

**Addresses:** FEATURES.md v0.16.0 Feature 1 weighted BH (table stakes).

**Avoids:** PITFALL 1 (weights not renormalized to mean=1), PITFALL 2 (missing genes in prior file silently excluded), PITFALL 11 (IHW thread safety — IHW stays as future rpy2 bridge, not v0.16.0 scope).

**Research flag:** Standard patterns — Genovese (2006) algorithm is well-documented (PMC3458887); statsmodels confirmed to lack a `weights` parameter (weighted BH in numpy is the correct approach).

### Phase 3: Region Restriction (`--restrict-regions`)

**Rationale:** Standalone; no dependency on other phases. Essential for cross-cohort analysis where different exome capture kits define different callable regions. bcftools-native via `bedtools intersect` + `GeneBedCreationStage` extension. Critical for the multi-cohort hardening theme of this milestone.

**Delivers:** `--restrict-regions` CLI flag; bedtools intersect call in `GeneBedCreationStage._process()`; chromosome name format validation before bcftools call; post-filter zero-variant warning when restriction is active; `.bed` file extension enforcement to avoid coordinate system ambiguity.

**Addresses:** FEATURES.md v0.16.0 Feature 4 region restriction (table stakes).

**Avoids:** PITFALL 6 (silent chr/no-chr prefix mismatch); ARCHITECTURE recommendation to use Option A (extend `GeneBedCreationStage`) rather than Option B (new stage).

**Research flag:** Standard patterns — bcftools `--regions-file` and `bedtools intersect` are well-documented.

### Phase 4: PCAComputationStage Wiring

**Rationale:** Standalone; 4-line change in `pipeline.py` plus fixing the VCF source fallback chain in `PCAComputationStage`. The stage is already implemented and registered in `stage_registry.py` (line 481); it simply is not connected to the pipeline execution graph. Medium priority: users can manually run PLINK2 and pass `--pca-file`, but automatic wiring removes a documentation-only workaround.

**Delivers:** `PCAComputationStage` active in pipeline when `--pca-tool akt` or `--pca-tool plink2` is specified; uses prefiltered VCF fallback chain (`prefiltered_vcf` → `extracted_vcf` → `vcf_file`); integration test that asserts PCA covariates appear in association results (this test must fail before the fix, pass after).

**Addresses:** FEATURES.md v0.16.0 Feature 7 PCAComputationStage wiring (post-MVP differentiator).

**Avoids:** PITFALL 7 (unwired stage producing no covariates despite appearing to run); PITFALL 8 (registry/import coupling — write the failing integration test first, then fix).

**Research flag:** Standard patterns — code reading confirms exact 4-line wiring gap and fix.

### Phase 5: Case-Confidence Sample Weights

**Rationale:** Depends on Phase 1 (COAST path must work to validate weighted COAST runs). Medium complexity due to the null model and eigenvalue computation changes. High scientific value for nephrology cohorts with uncertain phenotyping (EHR-derived diagnoses). Flag SKAT path as experimental pending permutation validation.

**Delivers:** `variantcentrifuge/association/sample_weights.py` (new file; mirrors `covariates.py` pattern); `--case-confidence-file` CLI flag; weighted GLM null model via `var_weights`; adjusted `phi = sqrt(mu*(1-mu) * sample_weights)` in `_compute_eigenvalues_filtered()`; effective sample size (ESS) diagnostic; permutation-based type I error validation confirming lambda_GC ~1.0 on permuted phenotype.

**Addresses:** FEATURES.md v0.16.0 Feature 2 case-confidence weights (table stakes for burden path; differentiator/experimental for SKAT path).

**Avoids:** PITFALL 3 (weights applied to residuals post-hoc rather than null model); ARCHITECTURE lazy null model caching pattern must be preserved (cohort-constant `sample_weights` passed consistently to all gene calls).

**Research flag:** Needs validation — weighted SKAT phi derivation (`phi = sqrt(mu*(1-mu) * sample_weights)`) is MEDIUM confidence (mathematical derivation from literature; no published reference Python implementation found for any genomics tool). Do not ship weighted SKAT without permutation validation.

### Phase 6: Sparse Genotype Matrices (Opt-In)

**Rationale:** Last because it has the most call-site audit risk (boolean indexing behavior differs between dense and sparse arrays; all sites in the augmentation loop must be audited). Memory benefit at GCKD scale (5K samples) is modest; essential only at UK Biobank scale (50K+ samples). Implement as opt-in via `--sparse-genotype-matrix` CLI flag or automatic threshold; do not make sparse the default for v0.16.0.

**Delivers:** `return_sparse: bool` parameter in `build_genotype_matrix()`; sparse-aware score vector computation (`score_vec = geno_weighted.T @ residuals`); dense conversion only at the eigendecomposition step (preserves the existing dimension-switching optimization); automatic threshold (enable when `n_samples * n_variants > 500,000` AND density < 20%); dense/sparse numerical equivalence tests.

**Addresses:** FEATURES.md v0.16.0 Feature 6 sparse genotype matrices (post-MVP, deferred).

**Avoids:** PITFALL 9 (scipy.linalg.eigh rejects sparse input — densify only at the eigendecomposition step after dimension-switching, not before score_vec computation); ARCHITECTURE anti-pattern 5 (never store full multi-gene genotype tensor in context).

**Research flag:** Needs profiling — the 500,000-cell breakeven threshold is an estimate from CSR literature (MEDIUM confidence). Profile `build_genotype_matrix()` on representative GCKD genes before setting the automatic threshold default.

### Deferred Features

**Single eigendecomposition optimization:** Clarify feature intent with the issue author before scoping. If it means reusing rho=0 eigenvalues for all rho, do not implement (mathematically incorrect; produces biased rho selection). If it means confirming the existing A = z1_half.T @ z1_half precomputation is already optimal, close the ticket as already-implemented. From STACK.md performance analysis: eigendecomposition on typical gene windows (p=5-50 variants) takes ~10 microseconds; 7 rho values × 20,000 genes = 1.4 seconds total wall time — not a measurable bottleneck given the 7-hour genotype_replacement stage.

**Dead code cleanup:** Dedicated sub-phase, not bundled with feature work. Run `vulture` or `deadcode` first; remove one stage at a time with full `make test-fast` between each removal; search for stage class names as strings (registry refs, logging, test assertions) before deleting.

### Phase Ordering Rationale

- Phase 1 (COAST fix) is a prerequisite for validating Phase 5 (case-confidence weights integrate with COAST); must come first
- Phases 2, 3, 4 are fully independent of each other and of Phase 1; can be developed in parallel
- Phase 5 (case-confidence weights) depends on Phase 1 being complete so that the COAST integration path can be validated end-to-end
- Phase 6 (sparse matrices) comes last because it has the most call-site audit scope and the lowest urgency at current cohort scale
- The COAST classification scoring model (JSON files only) can be developed in parallel with any phase since it requires no Python changes

### Research Flags

**Phases needing validation during implementation:**
- **Phase 5 (case-confidence SKAT weights):** Weighted phi adjustment (`phi = sqrt(mu*(1-mu) * sample_weights)`) is MEDIUM confidence — mathematical derivation only; no reference Python implementation found in any published genomics tool. Must validate via permutation: empirical lambda_GC on permuted phenotype with weights held fixed must be ~1.0. Do not ship without this validation.
- **Phase 6 (sparse matrices):** Breakeven threshold (500,000 cells, density < 20%) is a literature estimate. Profile on actual GCKD data before committing to automatic threshold behavior.

**Phases with standard patterns (skip research-phase):**
- **Phase 1 (COAST fix):** Root cause confirmed by direct code reading at exact line numbers (`analysis_stages.py` 2462-2487, 2630; `allelic_series.py` 551-581). Fix is mechanical.
- **Phase 2 (weighted BH):** Genovese (2006) algorithm is fully specified in PMC3458887; weight normalization requirement confirmed; statsmodels `multipletests` lacking `weights` param verified against docs.
- **Phase 3 (region restriction):** bcftools `--regions-file` is documented; bedtools intersect is standard; chromosome validation logic is straightforward.
- **Phase 4 (PCA wiring):** Code reading confirms the exact 4-line `pipeline.py` change. No algorithmic complexity.

---

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack — zero new dependencies | HIGH | All APIs verified against official docs (statsmodels 0.15.0, scipy 1.17.0); `scipy.sparse.csr_array` confirmed in 1.14.1 |
| Features — COAST bug root cause | HIGH | Direct code reading: `analysis_stages.py` lines 2462-2487 and 2630; failure mechanism fully traced |
| Features — Weighted BH design | HIGH | Genovese (2006) PMC3458887 verified; statsmodels lacking `weights` param confirmed against 0.15.0 dev docs |
| Features — COAST classification reference | HIGH | PMC10432147 full text; CRAN AllelicSeries vignette confirms annotations are user-supplied |
| Features — Region restriction (bcftools) | HIGH | Official bcftools manual; known duplicate-region issue #221 verified |
| Architecture — Integration points | HIGH | All modifications identified from direct codebase reading with line numbers confirmed |
| Architecture — PCAComputationStage gap | HIGH | Confirmed absent from `pipeline.py` imports; registered in `stage_registry.py` line 481 |
| Architecture — Sample weight null model math | MEDIUM | Mathematical derivation verified; no reference Python implementation found in any genomics tool |
| Pitfalls — BH renormalization | HIGH | FDR budget constraint verified against Genovese (2006); permutation test design is standard practice |
| Pitfalls — SKAT-O eigendecomp intent | HIGH | Per-rho eigenvalue requirement confirmed in SKAT-O paper (PMC3415556) and verified in code |
| Pitfalls — Sparse matrix breakeven | MEDIUM | Estimate from CSR literature and spVCF paper; not profiled on this specific codebase |
| Pitfalls — Case-confidence SKAT type I error | MEDIUM | No published reference implementation; type I error concern is theoretically sound but empirically unverified |

**Overall confidence:** HIGH for what to build, in what order, and at which exact code locations. MEDIUM for the statistical correctness of case-confidence weights in SKAT (requires permutation validation before shipping).

### Gaps to Address During Implementation

- **SKAT-O "single eigendecomposition" feature intent:** Clarify with the issue author before Phase 6 planning. The existing code already computes `A = z1_half.T @ z1_half` once and derives per-rho k_sym analytically from A. If this IS the intended optimization, close the ticket as already-done. If the intent is something different (e.g., rank-1 eigenvalue perturbation), assess mathematical correctness first.

- **Permutation validation infrastructure for case-confidence weights:** Phase 5 requires running SKAT with permuted phenotypes (weights held fixed) to validate lambda_GC ~1.0. This test harness does not exist in the current `tests/integration/` or `tests/performance/` directories and must be created as part of Phase 5 acceptance criteria.

- **Sparse matrix profiling baseline:** Before finalizing the Phase 6 automatic threshold, profile `build_genotype_matrix()` memory and time on a representative GCKD gene (typical size: 5125 samples × 5-50 rare variants). The 500,000-cell threshold is an estimate; actual GCKD characteristics may differ.

- **IHW future path:** IHW via rpy2 is explicitly deferred from v0.16.0. If users request it after release, the implementation path is clear: wrap the Bioconductor `IHW` R package via existing `rpy2` bridge, guard with `parallel_safe = False` (same constraint as R SKAT), add `--fdr-method ihw --ihw-covariate <col>` CLI flags. Document this extension point in the weighted BH implementation's docstring.

- **Missing-gene weight imputation default:** When a user provides a prior weight file that covers only a subset of tested genes, the unrepresented genes must receive a default weight. Weight=1.0 (neutral prior) is the recommended default per PITFALLS.md. Add `--missing-gene-weight` CLI flag with default 1.0 for transparency.

---

## Sources

### Primary (HIGH confidence)

- Genovese, Roeder & Wasserman (2006) weighted BH — budget constraint and normalization: [PMC3458887](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458887/)
- Wu et al. SKAT (2011) — score statistic formula, Beta(1,25) weights: [PMC3135811](https://pmc.ncbi.nlm.nih.gov/articles/PMC3135811/)
- Lee et al. SKAT-O (2012) — rho grid, per-rho eigenvalue requirement: [PMC3415556](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/)
- Liu et al. ACAT (2019) — Cauchy combination, ACAT-O: [PMC6407498](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407498/)
- McCaw et al. COAST (AJHG 2023) — BMV/DMV/PTV classification, user-supplied annotations: [PMC10432147](https://pmc.ncbi.nlm.nih.gov/articles/PMC10432147/)
- AllelicSeries CRAN vignette — annotations are user-supplied, not dictated by COAST method: [cran.r-project.org](https://cran.r-project.org/web/packages/AllelicSeries/vignettes/coast.html)
- Ignatiadis et al. IHW (Nature Methods 2016): [PMC4930141](https://pmc.ncbi.nlm.nih.gov/articles/PMC4930141/)
- SAIGE-GENE+ sparse matrix performance (Nature Genetics 2022): [Nature Genetics](https://www.nature.com/articles/s41588-022-01178-w)
- statsmodels GLM `var_weights` docs (0.15.0 dev) — verified parameter available: [statsmodels.org](https://www.statsmodels.org/dev/generated/statsmodels.genmod.generalized_linear_model.GLM.html)
- statsmodels `multipletests` docs (0.15.0) — confirmed no `weights` parameter: [statsmodels.org](https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html)
- scipy `false_discovery_control` (1.17.0) — confirmed no weights parameter: [scipy.org](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html)
- scipy sparse docs — CSR recommendation, new array API: [scipy.org](https://docs.scipy.org/doc/scipy/reference/sparse.html)
- bcftools manual — chromosome exact-match requirement, BED 0-based handling, indexed VCF requirement: [samtools.github.io](https://samtools.github.io/bcftools/bcftools.html)
- Codebase `analysis_stages.py` lines 2457-2487, 2630 — COAST genotype matrix failure paths (HIGH, direct read)
- Codebase `python_backend.py` lines 447-466 — SKAT-O rho loop eigenvalue structure (HIGH, direct read)
- Codebase `pipeline.py` — `PCAComputationStage` absent from imports (HIGH, direct read)
- Codebase `stage_registry.py` line 481 — `PCAComputationStage` registered but not wired (HIGH, direct read)
- Codebase `allelic_series.py` lines 551-581 — hard 3-category gate returning None (HIGH, direct read)
- Codebase `association/tests/allelic_series.py` — `&`-concatenated effect string not handled in `isin()` (HIGH, direct read)

### Secondary (MEDIUM confidence)

- Huang et al. gene-level IHW for rare variants (Genes 2022) — IHW covariate approach: [PMC8872452](https://pmc.ncbi.nlm.nih.gov/articles/PMC8872452/)
- spVCF sparse genotype matrix paper (Bioinformatics 2021) — density figures for rare variants: [Oxford Academic](https://academic.oup.com/bioinformatics/article/36/22-23/5537/6029516)
- Liability threshold model GWAS (Nature Genetics 2025) — case-confidence workaround: [Nature Genetics](https://www.nature.com/articles/s41588-025-02370-4)
- Phenotype misclassification SuperControl design (medRxiv 2025): [medRxiv](https://www.medrxiv.org/content/10.64898/2025.12.14.25342213v1.full)
- Population stratification in rare variant analysis — 10 PC default: [PMC6283567](https://pmc.ncbi.nlm.nih.gov/articles/PMC6283567/)
- scipy 1.16.0 release notes — Python 3.11+ constraint: [scipy.org](https://docs.scipy.org/doc/scipy/release/1.16.0-notes.html)
- REGENIE overview and binary trait handling: [rgcgithub.github.io/regenie](https://rgcgithub.github.io/regenie/overview/)

### Tertiary (LOW confidence)

- Case-confidence weighted SKAT phi derivation — theoretical extension; no published reference implementation in any genomics tool as of February 2026
- Madsen-Browning / WSS current adoption rate — single WebSearch source; correctly deferred as anti-feature

---

*Research completed: 2026-02-23*
*Synthesized from: STACK.md, FEATURES.md, ARCHITECTURE.md, PITFALLS.md (4 research streams)*
*Ready for roadmap: yes*
