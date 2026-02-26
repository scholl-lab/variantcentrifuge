# Pitfalls Research: v0.16.0 Association Hardening & Multi-Cohort Features

**Domain:** Adding hardening features to an existing rare variant association framework
**Focus:** v0.16.0 specific — weighted FDR, case-confidence weights, sparse matrices, single eigendecomposition, COAST classification fix, --restrict-regions, PCAComputationStage wiring, dead code cleanup
**Researched:** 2026-02-23
**Context:** VariantCentrifuge, 2001 tests passing, GCKD cohort (5,125 samples, 22 GB VCF). Clinical tool — correctness is paramount. All features added to an existing working system.

---

## Executive Summary

All eight features in v0.16.0 are additions or fixes to an already-working association framework. The danger is not "can we build this?" but "can we add this without breaking what exists?" The most dangerous pitfalls are:

1. **Weighted BH without renormalization** — weights that don't average to 1 silently violate the FDR guarantee. No error is raised; q-values look plausible.
2. **Case-confidence weights applied post-hoc to residuals** — the null model must be refit with sample weights; shortcuts inflate type I error.
3. **SKAT-O "single eigendecomposition" optimization** — the existing per-rho eigenvalue computation is already correct and fast; attempting to optimize it to a single decomposition is mathematically wrong.
4. **bcftools region filtering silent failures** — chromosome naming mismatches produce 0 variants without any error.
5. **Dead code removal registry coupling** — 8 stage classes are registered in two places; removing from one but not both causes startup failures.

---

## Critical Pitfalls

### Pitfall 1: Weighted BH — Weights That Do Not Average to 1

**What goes wrong:**
The Genovese (2006) weighted Benjamini-Hochberg procedure controls FDR at level alpha *only* when weights satisfy the budget constraint: the mean of all weights equals 1 (equivalently, the sum equals the number of hypotheses m). If gene-level prior weights are loaded from a file and applied without renormalization, FDR control is lost. When the mean weight exceeds 1, the procedure is anti-conservative (true FDR > alpha). When it is below 1, power is wasted.

The current `apply_correction()` in `variantcentrifuge/association/correction.py` takes raw p-values and applies standard (unweighted) BH via statsmodels. Extending it to weighted BH means computing `p_weighted_i = p_i / w_i` before the BH ranking step — but only after ensuring `w.mean() == 1`.

**Why it happens:**
Prior files typically store raw functional scores or odds ratios, not normalized weights. A researcher provides `weights.tsv` with values like `[0.1, 0.5, 2.0, 3.5]`. The code loads these and divides each p-value by the corresponding weight. The sum of weights is 6.1, not the number of genes, so the budget constraint is violated.

**Consequences:**
- Mean weight = 2.0 → effective FDR threshold doubles → 2x more false positives than nominal alpha.
- The bug is invisible: reported q-values look plausible and are in [0, 1].
- No error or warning is raised. Detection requires simulation.

**Prevention:**
Normalize weights at load time: `w = w / w.mean()`. Assert `abs(w.mean() - 1.0) < 1e-9` immediately after normalization. Add this assertion as a unit test. Document the normalization in the function docstring.

**Warning signs:**
- More genes pass FDR 5% in the weighted run than the unweighted run even for genes with unfavorable weights.
- Under permuted phenotypes (null), the weighted procedure produces more than 5% false positives.

**Detection test:** Under pure null (1000 permutations of phenotype), run weighted BH with the proposed weight scheme. Empirical FDR must be <= nominal alpha.

**Which phase:** Phase implementing gene-level prior weighting for FDR.

---

### Pitfall 2: Weighted BH — Missing Genes in the Prior File

**What goes wrong:**
Prior files (e.g., pLI scores, OMIM disease genes) cover a subset of tested genes. If the code drops genes not in the prior file, it silently removes genes from the analysis. If it assigns weight = 0, those genes can never pass any threshold. If it assigns weight = 1 without renormalizing the full vector, the budget constraint is violated because the represented genes collectively push the mean above or below 1.

**Why it happens:**
The implementation does a left-join of p-values onto the weight file, filling missing entries with `NaN` or dropping them. Neither behavior maintains FDR control.

**Consequences:**
- Systematic exclusion of novel (unrepresented) genes, which are the most likely candidates for new discoveries.
- FDR control lost if weights are not renormalized after imputation.

**Prevention:**
After loading the prior file, assign weight = 1.0 (or a configurable `--missing-gene-weight` default) to all genes not present in the file — representing "no prior information." Then renormalize the *full* weight vector (represented + unrepresented) so the mean is 1. Log a warning with the count of imputed genes: "Assigned weight=1.0 to 7,211 of 15,000 genes not in prior file."

**Warning signs:**
- Gene count in weighted results differs from gene count in unweighted results.
- Genes with missing weights produce `NaN` or `inf` corrected p-values.
- The prior file is described as covering 8,000 genes, but the analysis reports results for exactly 8,000 genes (missing the other 7,000).

**Which phase:** Phase implementing gene-level prior weighting for FDR.

---

### Pitfall 3: SKAT Sample-Level Case-Confidence Weights — Null Model Misspecification

**What goes wrong:**
Adding per-sample "case confidence" weights (e.g., phenotyping certainty scores) to SKAT requires refitting the null model with those weights — not multiplying residuals by the weights post-hoc. The current `PythonSKATBackend.fit_null_model()` in `variantcentrifuge/association/backends/python_backend.py` fits a `statsmodels GLM` without sample weights and stores `resid_response`. This cached residual vector is reused for all genes.

If sample weights are applied by multiplying `residuals * sample_weights` before computing `score_vec = geno_weighted.T @ residuals`, the score statistic changes but the eigenvalue distribution of the kernel matrix was computed assuming the unweighted residuals. The null distribution of Q is then wrong — it no longer matches the mixture of chi-squared distributions used to compute p-values.

**Why it happens:**
The shortcut is tempting because the null model is fit once and cached. Passing weights through to the score computation looks cheap and "feels right." The R SKAT `Get_Logistic_Weights()` function is designed for a different purpose (case-control imbalance adjustment for SKATBinary), not arbitrary per-sample confidence weights.

**Consequences:**
- Type I error inflation. Under permuted phenotypes with weights active, SKAT produces more rejections at alpha=0.05 than expected.
- Inflation scales with the variance of the weight distribution. Extreme weights (e.g., 0.1 to 10) produce large inflation.
- The bug is only detectable via permutation-based null simulation.

**Prevention:**
If sample weights are added, they must be passed to `statsmodels GLM.fit(freq_weights=sample_weights)` during null model fitting. The `_compute_eigenvalues_filtered()` method must also incorporate the weights into the projection (the hat matrix computation). This is a significant API change to the null model lifecycle — it requires version-bumping the `NullModelResult` dataclass and updating all call sites.

Validate with permutation: empirical type I error must match nominal alpha under permuted phenotype with weights held fixed.

**Warning signs:**
- QQ-plot inflation (lambda > 1.05) when running on permuted phenotypes with sample weights active.
- Weighted results on null data show systematically lower p-values than unweighted.

**Which phase:** Phase implementing sample-level case-confidence weights for weighted logistic/SKAT.

---

### Pitfall 4: SKAT-O "Single Eigendecomposition" — Mathematically Incorrect Optimization

**Background:**
The existing `_skato_get_pvalue()` in `variantcentrifuge/association/backends/python_backend.py` already computes *separate* eigenvalues for each rho via the analytical R.M^{1/2} formula. It does NOT do a single eigendecomposition. The code correctly computes the symmetric kernel `k_sym` for each rho and calls `_get_lambda(k_sym)` per rho (lines ~447-466 of `python_backend.py`).

**What would go wrong if optimized to a single decomposition:**
The kernel K_rho = G W R_rho W G' has rho-dependent eigenvalues through R_rho = (1-rho)I + rho*J. Only at rho=0 does K_rho reduce to the standard SKAT kernel. For rho > 0, the correlation structure R_rho non-trivially modifies the eigenspectrum. A "single eigendecomposition" would use the rho=0 eigenvalues for all rho values, producing incorrect per-rho p-values and biasing rho selection toward rho=0 (pure SKAT).

**Why this pitfall is listed despite the correct implementation:**
The milestone mentions "Single eigendecomposition for SKAT-O" as a feature. If this is interpreted as "compute eigenvalues once and reuse," it is incorrect. If it refers to something else (e.g., a single eigendecomposition of the base kernel A = Z1_half.T @ Z1_half, from which per-rho kernels are derived analytically), then the existing code already does exactly that — the base matrix A is computed once and k_sym is derived analytically per rho.

**Prevention:**
Clarify the feature intent before implementing. If "single eigendecomposition" means reusing rho=0 eigenvalues, do not implement it. If it means using the base A matrix (which is what the code already does), document that the feature is already implemented and close the ticket.

**Warning signs:**
- A PR proposes caching `lambda_all` outside the rho loop in `_skato_get_pvalue()`.
- Benchmark results show suspiciously identical p-values for rho=0.01 vs rho=0.04.
- SKAT-O results diverge from R SKAT after the optimization.

**Which phase:** SKAT-O optimization planning — needs clarification before the phase is scoped.

---

### Pitfall 5: COAST Classification — Multi-Transcript Annotation Conflicts

**What goes wrong:**
The `classify_variants()` function in `variantcentrifuge/association/tests/allelic_series.py` (line ~103) classifies variants using the EFFECT and IMPACT columns. The code does:
```python
is_ptv_effect = effect_series.str.strip().isin(PTV_EFFECTS)
```
where `PTV_EFFECTS = frozenset({"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"})`.

SnpEff can produce concatenated effects separated by `&`, e.g., `"splice_region_variant&synonymous_variant"` or `"stop_gained&splice_region_variant"`. An exact `isin()` match on the concatenated string fails: `"stop_gained&splice_region_variant"` is not in `PTV_EFFECTS`. The variant is then misclassified as code 0 (excluded from COAST) rather than PTV.

**Why it happens:**
SnpEff's annotation format produces compound consequence strings in certain annotation configurations. The exact format depends on the SnpEff version and the VCF annotation pipeline settings. This was likely not present in the initial test cases but appears in production data.

**Consequences:**
- PTV variants classified as code 0 and excluded from COAST.
- COAST produces `MISSING_CATEGORIES:PTV` skip reason for genes that clinically have known PTVs.
- N_PTV counts in COAST output are systematically lower than expected.
- A SnpEff version upgrade silently changes annotation format and causes a drop in classifiable COAST genes.

**Prevention:**
Change the PTV effect matching to split on `&` before matching:
```python
def _any_effect_matches(effect_str: str, effect_set: frozenset) -> bool:
    parts = effect_str.replace(",", "&").split("&")
    return any(p.strip() in effect_set for p in parts)
```
Apply the same split logic to the EFFECT column for all classification rules (is_ptv_effect, is_missense). Add unit tests with compound SnpEff effect strings to prevent regression.

**Warning signs:**
- `coast_skip_reason = "MISSING_CATEGORIES:PTV"` appearing for a large fraction of genes.
- N_PTV << expected based on gene burden table variant counts.
- Classification results change after a SnpEff version upgrade.

**Which phase:** COAST classification fix phase.

---

### Pitfall 6: --restrict-regions BED Intersection — Silent Chromosome Naming Failure

**What goes wrong:**
`bcftools view -R <bed_file>` requires exact chromosome name matching between the BED file and the VCF. The bcftools manual states: "Note that sequence names must match exactly, 'chr20' is not the same as '20'." If the VCF uses `chr1, chr2, ...` (UCSC style) and the BED file uses `1, 2, ...` (Ensembl style), bcftools silently skips all regions. The output VCF is empty or contains only variants from chromosomes with names that happen to match. No error is raised.

**Second issue — coordinate system:** BED files are 0-based, half-open intervals. bcftools detects BED format by file extension (`.bed` or `.bed.gz`, case-insensitive). If the user provides a regions file without the `.bed` extension, bcftools treats it as 1-based tab-delimited format. Variants at the exact start position of a region would be included in one mode and excluded in the other — a 1-position off-by-one error.

**Third issue — indexed VCF requirement:** `bcftools view -R` requires the input VCF/BCF to be indexed. If the prefilter intermediate file is unindexed, bcftools fails with a cryptic error about missing index rather than a clear "please index your VCF" message.

**Why it happens:**
Users provide BED files from various sources (UCSC Genome Browser, Ensembl, custom gene lists). GCKD VCFs use `chr` prefixes. An Ensembl-format BED without `chr` prefixes produces silent empty output.

**Consequences:**
- 0 qualifying variants in every gene after `--restrict-regions` is used. The pipeline reports no associations, which is indistinguishable from a legitimate null result.
- Off-by-one exclusion at region boundaries for improperly typed region files.

**Prevention:**
1. Before the bcftools call, validate chromosome names: read the first few lines of the BED file and the VCF header, compare chromosome name formats. If mismatch is detected, either normalize automatically (add/remove `chr` prefix) or raise a clear error.
2. Enforce `.bed` file extension or add explicit `--bed-format` vs `--regions-format` flag.
3. Post-filter sanity check: if `--restrict-regions` was specified and 0 variants pass prefiltering, emit a `WARNING: 0 variants passed region filter — check chromosome name format in BED file` before proceeding.
4. Verify the prefilter VCF is indexed before attempting region-filtered query.

**Warning signs:**
- 0 qualifying variants in all genes after `--restrict-regions` is added.
- The prefiltered VCF file exists but is empty or very small.
- Results with `--restrict-regions` differ from without by more than the expected regional subset.

**Which phase:** `--restrict-regions` BED intersection phase.

---

## Moderate Pitfalls

### Pitfall 7: PCAComputationStage — Unwired Stage Producing No Covariates

**What goes wrong:**
`PCAComputationStage` exists in `variantcentrifuge/stages/processing_stages.py` (line 2075) and is registered in `stage_registry.py` (line 481). However, it is NOT imported in `variantcentrifuge/pipeline.py` and is therefore never added to the pipeline's stage execution list. The stage runs AKT PCA and writes an eigenvec file to disk, but nothing downstream reads that file and passes it to `load_pca_file()` / `merge_pca_covariates()` in `variantcentrifuge/association/pca.py`.

The current PCA loading path uses `load_pca_file()` directly inside `AssociationAnalysisStage` when `--pca-file` is explicitly specified by the user. `PCAComputationStage` is intended to compute PCA on-the-fly and feed the result into that path automatically — but the context artifact (`context.pca_output_path` or similar) that downstream stages would read from is never defined or populated.

**Why it happens:**
The stage was created in anticipation of a full wiring but never connected to the pipeline graph. The stage registry registration gives a false sense that the stage is active. The `pipeline.py` import list is the authoritative source of what actually runs.

**Consequences:**
- Users who see `PCAComputationStage` in logs, documentation, or `--help` output assume PCA is computed automatically when `--pca-tool` is set.
- Actual behavior: the stage has an early-return guard when `pca_tool` is not set; when it is set, the stage runs and writes the file but PCA covariates are never loaded into the association analysis.
- Wiring requires: (1) importing in `pipeline.py`, (2) defining the dependency relationship with `BCFToolsPrefilterStage` or similar, (3) defining the context attribute that carries the output path, (4) reading that attribute in `AssociationAnalysisStage`.

**Prevention:**
Do not treat stage registry presence as equivalent to "wired." Verify wiring by checking `pipeline.py` imports and the stage graph. Before implementing the wiring, write an integration test that:
1. Runs the pipeline with `--pca-tool akt`.
2. Asserts that PCA covariates appear in the association results.
This test will fail on the current code and pass only when wiring is complete.

**Warning signs:**
- `PCAComputationStage` completes in logs but PCA covariates are absent from SKAT results.
- `context.pca_output_path` does not exist or is `None` after the stage runs.
- The stage's output file exists on disk but is not referenced by any other stage.

**Which phase:** PCAComputationStage wiring phase.

---

### Pitfall 8: Dead Code Cleanup — Registry and Import Coupling for 8 Stage Classes

**What goes wrong:**
The pipeline has two registration layers: `pipeline.py` (explicit import + instantiation in the stage list) and `stage_registry.py` (the global registry). Removing a stage class requires updating both layers in coordination. Missing one produces:

- **If removed from module but not from `pipeline.py` import:** `ImportError` at startup, breaking the entire pipeline.
- **If removed from module but registry still calls `register_stage(RemovedClass, ...)`:** `NameError` at registry initialization time, also breaking startup.
- **If removed from both but tests reference the class name:** test failures in `tests/unit/stages/` or integration tests.
- **If removed but checkpoint files reference the stage name:** existing checkpoints may fail to resume cleanly.

The stage files are large (analysis_stages.py is 4162 lines, processing_stages.py is 2138 lines). Dead stage classes are interspersed with live ones. Some dead stages (like `PCAComputationStage`, `StreamingDataProcessingStage`, `ParallelAnalysisOrchestrator`, `GenotypeFilterStage`) are registered in the registry but NOT imported in `pipeline.py`. Others may be imported in `pipeline.py` but never instantiated in the stage list.

**Why it happens:**
Dead code in a pipeline framework is hard to identify mechanically. The registry's `register_stage()` calls at the bottom of `stage_registry.py` look like a flat list — easy to miss a reference when deleting a class. String references (in logging, in test assertions, in stage names stored in checkpoint files) are not caught by import analysis.

**Consequences:**
- `ImportError` or `NameError` at startup if any reference is missed.
- Integration tests fail if they reference removed class names.
- Existing user checkpoints fail to resume if stage names change.

**Prevention:**
1. Run `vulture` or `deadcode` on the codebase before deciding which stages are "dead." These tools find unreferenced symbols.
2. Remove stages one at a time with a full `make test-fast` run between each removal, not as a batch.
3. For each stage being removed, search for its name as a string (for registry refs, logging, and test assertions): `grep -rn "ClassName\|class_name" variantcentrifuge/ tests/`.
4. Check `tests/unit/stages/` and integration tests for explicit stage class references before removing.
5. Do not rename any stage that has appeared in checkpoint files — add a compatibility alias if needed.

**Warning signs:**
- `NameError` at import time after removing a class from a module.
- `KeyError` in registry lookup at startup.
- Integration tests fail because they instantiate or reference the removed stage class.
- CI passes on the refactored code but existing user checkpoints break on resume.

**Which phase:** Dead code cleanup phase. The cleanup should be a dedicated phase, not bundled with other feature work.

---

### Pitfall 9: Sparse Genotype Matrix — Eigendecomposition Requires Dense Input

**What goes wrong:**
For a sparse CSR genotype matrix G (n_samples x n_variants), the expression `G @ G.T` is mathematically correct for computing the kernel. However, the eigenvalue step in `_compute_eigenvalues_filtered()` calls `scipy.linalg.eigh(eig_mat / 2.0, ...)`. `scipy.linalg.eigh` requires a dense array. Passing a sparse matrix either fails with a `TypeError` or silently calls `.toarray()` internally, potentially triggering an out-of-memory allocation for large n.

Additionally, the existing code already has the optimization for dimension choice:
```python
eig_mat = z_adj.T @ z_adj if n_variants <= geno_weighted.shape[0] else z_adj @ z_adj.T
```
This ensures the smaller matrix is decomposed (p x p vs n x n). If sparse matrices are added naively, this dimension check must be preserved. For the GCKD cohort (n=5125), decomposing a 5125 x 5125 dense matrix for every gene defeats the purpose of sparsity.

**Consequences:**
- `TypeError` when `scipy.linalg.eigh` receives a sparse matrix.
- Memory explosion if `.toarray()` is called on a n x n sparse matrix for large n.
- Incorrect results if `scipy.sparse.linalg.eigsh` is used instead (it returns only `k` eigenvalues; missing small eigenvalues change the null distribution tail).

**Prevention:**
Keep sparse representation only for the score vector computation (`score_vec = geno_weighted.T @ residuals`), where sparsity genuinely reduces computation. At the eigenvalue computation step, explicitly call `.toarray()` after computing the smaller dimension matrix (p x p for most rare variant genes with few variants). Preserve the n_variants <= n_samples branching. For a gene with 10 variants and 5125 samples, the dense p x p matrix is only 10x10 — the toarray() cost is trivial.

**Warning signs:**
- Memory usage does not decrease when sparse matrices are enabled for large genes.
- `TypeError: input must be a square 2-D array` from scipy.linalg.eigh.
- P-values differ from dense mode on the same data after sparse conversion.

**Which phase:** Sparse genotype matrix implementation phase.

---

### Pitfall 10: COAST Configurable Scoring Model — Weight Vector Length Mismatch

**What goes wrong:**
The COAST burden tests use `coast_weights = [1.0, 2.0, 3.0]` (BMV, DMV, PTV) with `_CATEGORY_ORDER = [1, 2, 3]` defined as a module-level constant in `variantcentrifuge/association/backends/coast_python.py`. The function `_aggregate_by_category()` uses `zip(weights, ...)` (strict=True) which raises `ValueError` if the weight vector length does not match `len(_CATEGORY_ORDER) = 3`.

If the configurable scoring model allows different numbers of categories (e.g., splitting DMV into "possibly damaging" and "probably damaging" for a 4-category scheme), the weight vector must be updated in lockstep with `_CATEGORY_ORDER`, `_compute_allelic_skat_weights()`, and the Cauchy combination weight vector `_COAST_CAUCHY_WEIGHTS`.

**Prevention:**
Define a single `N_COAST_CATEGORIES` constant. Validate `len(coast_weights) == N_COAST_CATEGORIES` at the entry point of `test_gene()` with a clear `ValueError`. If the number of categories changes, update all four places atomically and add a unit test that checks all vectors are the same length.

**Warning signs:**
- `ValueError: zip() arguments have different lengths` during COAST testing.
- COAST returns different numbers of burden components across runs.

**Which phase:** COAST configurable scoring model phase.

---

## Minor Pitfalls

### Pitfall 11: IHW Optional Dependency — R Heap and Thread Safety

**What goes wrong:**
IHW (Independent Hypothesis Weighting) is a Bioconductor R package. Integrating it via `rpy2` would introduce the same single-threaded R constraint that currently forces `COASTTest` and `RSKATTest` to declare `parallel_safe=False`. If IHW is added without isolating it behind the existing R-backend conditional path, it may silently break the parallel execution model.

**Prevention:**
Keep IHW in the R backend path, parallel-guarded the same way as other R tests. Provide Python-native weighted BH (simple normalization approach) as the default, making IHW an opt-in advanced mode. Users on HPC who need high-performance IHW can use the R backend; Python-backend users get the simpler but correct weighted BH.

**Which phase:** Weighted BH / IHW phase.

---

### Pitfall 12: Ultra-Rare Variants — Low MAC Produces p=1.0 Instead of Skip

**What goes wrong:**
For genes where all qualifying variants are carried by only 1-2 samples out of 5125, the SKAT score vector is near-zero for all variants. The Q statistic at rho=0 approaches zero. The omnibus integration in `_skato_integrate_davies()` may underflow and return 0.0 for the integral, producing `pvalue = 1.0 - 0.0 = 1.0`. After the Bonferroni guard (`pmin * 7 < pvalue`), the final clipped result is `1.0`.

This is statistically correct (a gene with only 1-2 carriers has essentially no power) but misleading: the user sees `p_value = 1.0` rather than `p_value = None` (skip). The gene is included in multiple testing correction with a trivially non-significant p-value instead of being excluded as untestable.

**Prevention:**
Add a minimum allele count (MAC) guard before running SKAT: if `geno.sum() < 5` (total alternate allele count across all variants in the gene), return `skip_reason = "low_mac"` with `p_value = None`. This matches standard practice in rare variant analysis (most tools require MAC >= 5 for reliable p-values) and matches how the existing rank check works. Log at DEBUG level with the MAC value.

**Which phase:** Any phase that changes data access patterns — specifically relevant if sparse matrix support changes how monomorphic variants are handled.

---

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation |
|-------------|---------------|------------|
| Weighted BH (gene prior) | Weights not normalized to mean=1 | Normalize at load, assert mean=1, permutation test |
| Weighted BH (missing genes) | Genes absent from prior file excluded | Impute weight=1, renormalize full vector, log count |
| Case-confidence SKAT | Weights applied to residuals, not null model | Pass weights to `GLM.fit(freq_weights=...)`, permutation validate |
| SKAT-O single eigendecomposition | Reusing rho=0 eigenvalues for all rho | DO NOT do this — per-rho eigenvalues are required; clarify feature scope |
| COAST classification | SnpEff "&" concatenated effects excluded | Split effect strings on "&" before `isin()` matching |
| COAST classification | ANN[0] ordering changes with SnpEff config | Validate HIGH-impact variants appear in PTV category |
| COAST configurable scoring | Weight vector length mismatch with categories | Validate len(weights) == N_COAST_CATEGORIES at entry point |
| --restrict-regions | chr vs no-chr prefix mismatch (silent!) | Validate chromosome names before bcftools call; warn on 0 output |
| --restrict-regions | File extension drives coordinate interpretation | Enforce .bed extension; document in CLI help |
| Sparse genotype matrices | scipy.linalg.eigh rejects sparse input | Densify after dimension-switching; keep sparse only for score_vec |
| PCAComputationStage wiring | Stage registered but not wired in pipeline.py | Define context contract, write failing integration test first |
| Dead code cleanup (8 stages) | Registry + import coupling | One stage at a time, full test run between each |
| Dead code cleanup | Checkpoint files reference removed stage names | Search for stage names as strings, not just imports |
| Ultra-rare variants | p=1.0 instead of skip for MAC < 5 | Add MAC guard before SKAT; return skip_reason="low_mac" |
| IHW optional dependency | R heap thread restriction | Keep in R backend path only; Python default uses simple weighted BH |

---

## Sources

- [False discovery control with p-value weighting (Genovese 2006)](https://academic.oup.com/biomet/article-abstract/93/3/509/380687) — budget constraint: weights must average to 1
- [Weighted False Discovery Rate Controlling Procedures — PMC6474384](https://pmc.ncbi.nlm.nih.gov/articles/PMC6474384/) — extreme weight caution; pre-specification requirement
- [SKAT Package — CRAN PDF](https://cran.r-project.org/web/packages/SKAT/SKAT.pdf) — sample weights in null model; type I error control for binary traits
- [Optimal Unified Approach for Rare-Variant Association (SKAT-O) — PMC3415556](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/) — rho grid and eigenvalue structure; confirms per-rho eigenvalues are required
- [IHW — Bioconductor introduction](https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html) — covariate independence requirements
- [Data-driven hypothesis weighting (IHW) — PMC4930141](https://pmc.ncbi.nlm.nih.gov/articles/PMC4930141/) — IHW pitfalls with correlated covariates in genomics
- [bcftools manual](https://samtools.github.io/bcftools/bcftools.html) — chromosome name exact-match requirement; BED 0-based half-open interval handling; indexed VCF requirement
- [Annotating variants with VEP — HUMU 2022](https://onlinelibrary.wiley.com/doi/10.1002/humu.24298) — multi-transcript consequence conflicts
- [COAST allelic series test (McCaw 2023) — PMC10432147](https://pmc.ncbi.nlm.nih.gov/articles/PMC10432147/) — BMV/DMV/PTV classification requirements
- Codebase: `variantcentrifuge/association/backends/python_backend.py` — SKAT-O rho grid implementation (verified: per-rho eigenvalues already computed correctly)
- Codebase: `variantcentrifuge/association/correction.py` — current FDR implementation (unweighted; normalization not present)
- Codebase: `variantcentrifuge/association/tests/allelic_series.py` — COAST classification logic (verified: `&` concatenation not handled in `isin()` call)
- Codebase: `variantcentrifuge/pipeline.py` — stage wiring (verified: `PCAComputationStage` absent from imports and stage list)
- Codebase: `variantcentrifuge/stages/stage_registry.py` — registry (verified: `PCAComputationStage` registered at line 481 but not wired)
