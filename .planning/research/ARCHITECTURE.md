# Architecture Research: v0.16.0

**Domain:** Association hardening and multi-cohort features added to existing v0.15.0 framework
**Researched:** 2026-02-23
**Sources:** Direct codebase read — all findings are HIGH confidence unless noted

---

## Existing Architecture Baseline (v0.15.0 State)

### Pipeline Execution Flow (Relevant Portion)

```
CLI args -> build_pipeline_stages() -> PipelineRunner -> Stage DAG

  ConfigurationLoadingStage
    -> SampleConfigLoadingStage
      -> bcftools_prefilter (VariantExtractionStage)
        -> [PCAComputationStage — defined in processing_stages.py, NOT WIRED in pipeline.py]
        -> SnpSiftFilterStage -> FieldExtractionStage -> GenotypeReplacementStage
          -> DataFrameLoadingStage
            -> AssociationAnalysisStage (dependencies: {"dataframe_loading", "sample_config_loading"})
              -> write .association.tsv[.gz]
```

### Association Subsystem Components

| Module | Location | Role |
|--------|----------|------|
| `AssociationEngine` | `association/engine.py` | Per-gene loop, parallel dispatch, ACAT-O, FDR |
| `AssociationConfig` | `association/base.py` | Dataclass — all test parameters; all fields default to backward-compatible values |
| `AssociationTest` (ABC) | `association/base.py` | One subclass per test; `run(gene, contingency_data, config) -> TestResult` |
| `TestResult` | `association/base.py` | Per-gene result container with `extra` dict for test-specific data |
| `apply_correction` | `association/correction.py` | Leaf module: FDR/Bonferroni via statsmodels, no vc imports |
| `load_covariates` | `association/covariates.py` | Covariate file -> aligned float64 matrix |
| `build_genotype_matrix` | `association/genotype_matrix.py` | Per-sample GT columns -> (n_samples, n_variants) dense ndarray |
| `PythonSKATBackend` | `association/backends/python_backend.py` | SKAT/Burden/SKAT-O numpy/scipy/statsmodels |
| `PurePythonCOASTTest` | `association/tests/allelic_series_python.py` | COAST test; `parallel_safe=True` |
| `AssociationAnalysisStage` | `stages/analysis_stages.py` (line 2302) | Pipeline stage wiring all of the above |
| `PCAComputationStage` | `stages/processing_stages.py` (line 2075) | AKT PCA runner — EXISTS but not wired |

### Key Data Structures in AssociationAnalysisStage

The stage passes a `gene_data` dict (one per gene) to `engine.run_all()`:

```python
gene_data = {
    "GENE": str,
    "proband_count": int,
    "control_count": int,
    "proband_carrier_count": int,
    "control_carrier_count": int,
    "proband_allele_count": int,
    "control_allele_count": int,
    "n_qualifying_variants": int,
    # Added for regression tests:
    "genotype_matrix": np.ndarray,     # (n_samples, n_variants), float64
    "variant_mafs": np.ndarray,        # (n_variants,)
    "phenotype_vector": np.ndarray,    # (n_samples,), 0.0/1.0
    "covariate_matrix": np.ndarray | None,  # (n_samples, k)
    "gene_df": pd.DataFrame,           # per-variant annotations for COAST
    "vcf_samples": list[str],          # ordered sample IDs
    # Optional functional annotations for weight schemes:
    "cadd_scores": np.ndarray | None,
    "revel_scores": np.ndarray | None,
    "variant_effects": np.ndarray | None,
}
```

---

## Integration Point Analysis

### 1. COAST Fix: Genotype Matrix Availability

**Root cause confirmed by reading analysis_stages.py lines 2457-2487:**

The stage has two code paths for obtaining per-sample GT columns:

**Path A (line 2462-2472):** GT not yet in df AND vcf_samples present.
- Calls `reconstruct_gt_column(df.copy(), vcf_samples)` which drops per-sample columns.
- Saves original df as `df_with_per_sample_gt` first. This path works correctly.

**Path B (line 2473-2486):** GT already in df AND needs_regression.
- Tries `context.variants_df` as a fallback for per-sample GT columns.
- This path is fragile: `variants_df` may not contain per-sample columns if the
  DataFrame has been processed by another stage that dropped them.

**The failure mode:** If `GeneBurdenAnalysisStage` runs before `AssociationAnalysisStage`
(both can be in the pipeline), or if the DataFrame loading path does not preserve
per-sample GT columns, `df_with_per_sample_gt` stays `None`, `gt_columns_for_matrix`
is empty, and the genotype matrix augmentation block at line 2630 is skipped entirely:

```python
# analysis_stages.py line 2630
if needs_regression and gt_columns_for_matrix and vcf_samples_list:
    # This block is entirely skipped when gt_columns_for_matrix is empty
```

Each `gene_data` dict then has no `genotype_matrix` key. COAST's `run()` method
(allelic_series_python.py line 213) returns `TestResult(extra={"coast_skip_reason": "NO_GENOTYPE_MATRIX"})` for every gene.

**The fix:**

A third fallback should reconstruct per-sample columns from the packed GT string
using `create_sample_columns_from_gt_vectorized()` which is already defined in
`analysis_stages.py` (lines 66-150) and handles the `"Sample(0/1);Sample2(0/0)"`
format.

```python
# After Path A and Path B checks, add Path C:
elif "GT" in df.columns and needs_regression and vcf_samples_list:
    from ..stages.output_stages import _find_per_sample_gt_columns
    if not _find_gt_columns(df):  # no per-sample columns in current df
        logger.info("Reconstructing per-sample GT columns from packed GT for genotype matrix")
        df_gt_expanded = create_sample_columns_from_gt_vectorized(
            df.copy(), vcf_samples_list
        )
        df_with_per_sample_gt = df_gt_expanded
```

**Files to modify:** `variantcentrifuge/stages/analysis_stages.py` only.

---

### 2. Gene-Level FDR Weighting

**Current FDR path (engine.py lines 443-449):**

```python
raw_pvals = [acat_o_results[g].p_value for g in testable_genes]
corrected = apply_correction(raw_pvals, self._config.correction_method)
```

`apply_correction()` in `correction.py` is a clean leaf module that calls
`statsmodels.stats.multitest.multipletests()`. It takes a flat list with no weights.

**Integration point:** Weighted BH requires dividing each p-value by its weight
before applying standard BH. `correction.py` is the right module (it already has
the `smm` import and the docstring pattern). Add a second function:

```python
def apply_weighted_correction(
    pvals: list[float] | np.ndarray,
    weights: list[float] | np.ndarray,
    method: str = "fdr",
) -> np.ndarray:
    """
    Apply IHW-style weighted Benjamini-Hochberg correction.

    Adjusts each p-value by dividing by its weight (p_adj = p / w),
    then applies standard BH to the adjusted values. Weights should
    be normalized to have mean 1.0 (so the effective number of tests
    is unchanged).
    """
    pvals_arr = np.asarray(pvals, dtype=float)
    weights_arr = np.asarray(weights, dtype=float)
    # Normalize weights so mean = 1
    weights_arr = weights_arr / weights_arr.mean()
    adjusted = pvals_arr / np.maximum(weights_arr, 1e-10)
    return apply_correction(np.minimum(adjusted, 1.0), method)
```

`AssociationEngine.run_all()` checks `self._config.fdr_weights` and dispatches to
`apply_weighted_correction()` when weights are provided.

**Files to modify:**
- `variantcentrifuge/association/correction.py` — new `apply_weighted_correction()` function
- `variantcentrifuge/association/base.py` — new `AssociationConfig` field:
  `fdr_weights: str | None = None` (values: `"n_variants"`, `"uniform"`, path to gene-weight file)
- `variantcentrifuge/association/engine.py` — dispatch to weighted function when configured

**No new files needed.**

---

### 3. Case-Confidence Weights (Sample-Level Weights)

**Current state:** `phenotype_vector` is a hard 0/1 array (analysis_stages.py line 2571).
No per-sample weight support anywhere in the engine, backend, or test classes.

**Where weights enter the computation:**

Statsmodels GLM supports `var_weights` in `.fit()`. For probabilistic phenotype
(confidence that a sample is truly a case), `var_weights` is the appropriate
mechanism — it scales each observation's contribution to the log-likelihood.

The score statistic must also be adjusted: `Q = ||G_w^T (r * sqrt(var_w))||^2 / 2`,
where `r` are the weighted residuals from GLM and the eigenvalue computation uses
`phi = sqrt(mu*(1-mu)*var_w)` instead of plain `sqrt(mu*(1-mu))`.

**Data flow:**

```
AssociationConfig.case_confidence_file: str | None
    |
    v
AssociationAnalysisStage._process()
    load_sample_weights(filepath, vcf_samples_list)
    -> weight_vector: np.ndarray(n_samples,)    <- float in [0, 1]
    |
    v [in gene augmentation loop, line ~2743]
    gene_data["sample_weights"] = weight_vector  (or masked version)
    |
    v
engine.run_all() -> test.run(gene, gene_data, config)
    |
    v [in SKAT/burden tests only]
    null_model = backend.fit_null_model(
        phenotype, covariates, trait_type,
        sample_weights=gene_data.get("sample_weights")   <- NEW parameter
    )
    |
    v [inside fit_null_model — python_backend.py line ~619]
    glm = sm.GLM(phenotype, design_matrix, family=family)
    fit_result = glm.fit(var_weights=sample_weights)
    |
    v [inside _compute_eigenvalues_filtered — line ~742]
    phi = np.sqrt(mu_hat * (1.0 - mu_hat) * sample_weights)  <- weight-adjusted
```

**Files to modify:**
- `variantcentrifuge/association/base.py` — `AssociationConfig` + 2 new fields:
  `case_confidence_file: str | None = None`, `case_confidence_column: str | None = None`
- `variantcentrifuge/association/backends/python_backend.py` — `fit_null_model()`,
  `_compute_eigenvalues_filtered()`, `_test_skato()` (Z1 projection must use weighted phi)
- `variantcentrifuge/stages/analysis_stages.py` — load weight file, pass to gene_data

**Suggested new file:** `variantcentrifuge/association/sample_weights.py`
(same pattern as `covariates.py`: read file, align to vcf_samples order, return ndarray).
This keeps the I/O concern separated from the covariate matrix concern.

Fisher test ignores sample weights (count-based, no regression). All burden/SKAT
tests benefit. No changes to `correction.py`, `engine.py` routing, or output schema.

---

### 4. Sparse Genotype Matrices

**Current construction (genotype_matrix.py):**

`build_genotype_matrix()` returns a dense `np.ndarray(n_samples, n_variants, float64)`.
For n_samples=5000 and n_variants=100, this is 4 MB per gene — manageable.
For n_variants=1000+ (pan-exome runs), this is 40 MB per gene times parallel workers.

**Impact on python_backend.py kernel operations:**

The key operations in `_test_skat()` (lines 863-900):

```python
score_vec = geno_weighted.T @ residuals  # (p,) — matmul
q_stat = score_vec @ score_vec / 2.0    # scalar
lambdas = self._compute_eigenvalues_filtered(geno_weighted, null_model)  # eigh
```

For SKAT-O (lines 1052-1094):

```python
a_mat = z1_half.T @ z1_half             # (p, p) — matmul
```

scipy sparse CSC matrix supports `@` and `T@` but the result of `sparse.T @ sparse`
is dense for small p. The benefit is in the matmul with residuals: for a matrix with
5% fill (rare variants in a large cohort), `sparse.T @ vector` saves ~20x memory and
~10x time.

**Practical threshold:** Sparse format pays off when n_variants/n_samples * fill_fraction
produces enough zeros. For rare variant tests, fill fraction is typically 1-5%
(most samples are wild-type). Sparse format is worth implementing for n_samples > 1000
or n_variants > 200.

**Implementation path:**

1. `build_genotype_matrix()` gains a `return_sparse: bool` parameter.
2. When `return_sparse=True`, return `scipy.sparse.csc_array` instead of ndarray.
3. `PythonSKATBackend.test_gene()` accepts `Union[np.ndarray, csc_array]` and converts
   with `np.asarray(geno)` only at the eigendecomposition step (eigh requires dense).
4. The `geno_weighted = geno * weights` step must use `geno.multiply(weights)` for
   sparse input.

**Files to modify:**
- `variantcentrifuge/association/genotype_matrix.py` — optional sparse output
- `variantcentrifuge/association/backends/python_backend.py` — sparse-aware matmuls
- `variantcentrifuge/association/base.py` — `use_sparse_geno: bool = False`
- `variantcentrifuge/stages/analysis_stages.py` — pass flag to `build_genotype_matrix()`

**Risk:** The sample mask application (`geno[mask_arr]`) uses boolean indexing that works
differently on sparse arrays (`geno[mask_arr, :]` required). All sites in the
augmentation loop at lines 2657-2678 that reshape `geno` must be audited.

---

### 5. Single Eigendecomposition Across SKAT-O Rho Grid

**Current path confirmed in python_backend.py lines 447-466:**

`_skato_get_pvalue()` loops over 7 rho values in `_SKATO_RHO_GRID`. For each rho, it
constructs `k_sym` (p x p matrix) and calls `_get_lambda(k_sym)` which calls
`scipy.linalg.eigh(k_sym)`. This is 7 separate O(p^3) eigendecompositions.

**The mathematical structure:**

The base kernel `A = z1_half.T @ z1_half` is the same for all rho (computed at
line 426, `a_mat = z1_half.T @ z1_half`). Each `k_sym(rho)` is:

```
k_sym(rho) = s^2 * A  +  s*delta * (J@A + A@J)  +  delta^2 * J@A@J
```

where `s = sqrt(1-rho)`, `delta = (sqrt(1-rho+p*rho) - sqrt(1-rho)) / p`.

If `A = V @ diag(d) @ V.T` (precomputed once), then the eigenvalues of `k_sym(rho)`
can be found from the eigenvalues of a rank-1 update of `diag(d)` in the V basis.
Specifically, letting `u = V.T @ z_mean_normalized`, the eigenvalues of `k_sym(rho)`
are the eigenvalues of `s^2 * diag(d) + (correction terms)` which reduces to a
rank-1 update solvable in O(p^2) via the secular equation.

**Implementation location:** `_skato_get_pvalue()` in `python_backend.py`.

The eigendecomposition of A can be performed once before the rho loop:

```python
# Before the rho loop (line 447):
eig_vals_a, eig_vecs_a = scipy.linalg.eigh(a_mat, driver="evr")

# Inside the loop: compute k_sym eigenvalues from A's eigenvectors
# using rank-1 update instead of full eigh
```

This is a pure internal optimization: same API, same outputs, no schema changes.
Measured speedup is O(7x) at the eigendecomposition step, significant for large p.

**Files to modify:** `variantcentrifuge/association/backends/python_backend.py` only.
No API changes, no new files.

---

### 6. `--restrict-regions` Integration

**Current bcftools invocation (filters.py lines 157-184):**

```python
cmd = [
    "bcftools", "view",
    "--threads", threads,
    "-W",           # auto-index
    "-R", bed_file, # genomic regions from gene BED
]
if cfg.get("bcftools_prefilter"):
    cmd.extend(["-i", cfg["bcftools_prefilter"]])
cmd.extend(["-Oz", "-o", output_file, vcf_file])
```

The `-R` flag restricts to gene BED intervals. A `--restrict-regions` feature would
add a second BED mask (e.g. exome capture intervals, callable regions) that intersects
with the gene BED before passing to bcftools.

**Why bcftools cannot take two `-R` flags:** It does not support multiple region files.
The solution is pre-intersection via `bedtools intersect`.

**Integration options and recommendation:**

Option A — Extend `GeneBedCreationStage`: After creating the gene BED, if
`restrict_regions_bed` is configured, call `bedtools intersect -a gene.bed -b restrict.bed`
and store the result as `context.config["bed_file"]`. The existing `VariantExtractionStage`
already uses `context.config["bed_file"]` unchanged.

Option B — New `RegionRestrictionStage`: A dedicated stage between `GeneBedCreationStage`
and `VariantExtractionStage` with dependency on `"gene_bed_creation"`.

**Recommendation: Option A** — minimal diff, no new stage registration, uses existing
bedtools dependency. `GeneBedCreationStage` already writes to a path and sets
`context.config["bed_file"]`; the restriction can be applied in-place.

**Files to modify:**
- `variantcentrifuge/stages/processing_stages.py` — `GeneBedCreationStage._process()`:
  add bedtools intersect call when `restrict_regions_bed` is set
- `variantcentrifuge/cli.py` — add `--restrict-regions` argument
- `variantcentrifuge/association/base.py` — OR general config; this is a pre-processing
  concern not test-specific; keeping it in general `context.config` is cleaner

`filters.py` needs no changes.

---

### 7. PCAComputationStage Wiring

**Current state confirmed by code inspection:**

`PCAComputationStage` (processing_stages.py line 2075):
- `dependencies = {"bcftools_prefilter"}` — runs after bcftools prefilter VCF is ready
- `parallel_safe = True` — subprocess only
- Writes eigenvec output to `context.config["pca_file"]`
- Uses `context.config.get("vcf_file")` which is the RAW input VCF, not the prefiltered VCF

`pipeline.py` `build_pipeline_stages()`:
- Does NOT import `PCAComputationStage`
- Does NOT reference `PCAComputationStage` anywhere

`AssociationAnalysisStage._process()` at lines 2544-2561:
- Reads `assoc_config.pca_file` and merges via `load_pca_file()` / `merge_pca_covariates()`
- This already works for pre-computed PCA files
- `PCAComputationStage` would populate `context.config["pca_file"]` for this code path

**Wiring requires two changes:**

Change 1 — pipeline.py (4 lines):

```python
# Add to imports at top of pipeline.py:
from .stages.processing_stages import (
    ...,
    PCAComputationStage,
)

# In build_pipeline_stages(), after VariantExtractionStage:
pca_tool = getattr(args, "pca_tool", None) or config.get("pca_tool")
if pca_tool:
    stages.append(PCAComputationStage())
```

Change 2 — PCAComputationStage._process() (3 lines):

The stage currently uses `context.config.get("vcf_file")` (raw input VCF). For
population stratification, the prefiltered VCF (which contains the same sample set
as the association analysis) is more appropriate. Change to:

```python
vcf_file = (
    context.config.get("prefiltered_vcf")
    or context.config.get("extracted_vcf")
    or context.config.get("vcf_file")
)
```

**Files to modify:**
- `variantcentrifuge/pipeline.py` — import + conditional append
- `variantcentrifuge/stages/processing_stages.py` — `PCAComputationStage._process()`:
  use prefiltered VCF, add fallback chain

---

### 8. COAST Classification as a Scoring Model

**Scoring system pattern (read from `scoring/nephro_candidate_score/formula_config.json`
and `variantcentrifuge/scoring.py`):**

`apply_scoring(df, scoring_config)` evaluates pandas-eval formulas row-by-row on the
variant-level DataFrame. Input and output are both row-level (one row = one variant).

**COAST classification output taxonomy:**

- Variant-level: BMV (1), DMV (2), PTV (3), unclassified (0). These come from
  `classify_variants()` in `allelic_series.py` using EFFECT, IMPACT, SIFT, PolyPhen
  columns — all present in the per-variant DataFrame.
- Gene-level: COAST p-value, n_bmv, n_dmv, n_ptv. These come from the association
  engine after all samples are processed. They are NOT per-variant.

**Can COAST classification be a scoring model?**

Yes, for the variant-level classification. No, for the gene-level test result.

The `classify_variants()` logic maps each variant row to a category code using
EFFECT/IMPACT columns. This is expressible as a scoring model formula:

```json
// scoring/coast_classification/formula_config.json
{
  "output_scores": ["coast_variant_class"],
  "formulas": [
    {
      "is_ptv": "(impact == 'HIGH') & (effect.str.contains('stop_gained|frameshift_variant|splice_donor|splice_acceptor|start_lost', na=False))"
    },
    {
      "is_dmv": "(effect == 'missense_variant') & ((sift_pred == 'D') | (polyphen_pred.str.startswith('probably', na=False)) | (polyphen_pred.str.startswith('possibly', na=False)))"
    },
    {
      "is_bmv": "(effect == 'missense_variant') & (~is_dmv)"
    },
    {
      "coast_variant_class": "is_ptv * 3 + is_dmv * 2 + is_bmv * 1"
    }
  ]
}
```

This would appear in the main variant TSV and Excel output as a `coast_variant_class`
column, giving users per-variant functional classification without running association.

**Recommended architecture:** Two-layer design:

1. **Variant-level scoring model** (`scoring/coast_classification/`): New config files,
   no Python changes. Integrates with `VariantScoringStage` like all scoring models.
   Available to all users regardless of association analysis.

2. **Gene-level COAST test** remains in `AssociationAnalysisStage` as a registered
   `AssociationTest`. COAST p-values appear only in `.association.tsv`.

**Files to create:**
- `scoring/coast_classification/formula_config.json` (above pattern)
- `scoring/coast_classification/variable_assignment_config.json` (maps EFFECT/IMPACT/SIFT/PolyPhen column names)

**No Python changes needed for the scoring model.**

---

## Data Flow: Sample Weights Through the Engine

```
[--case-confidence-file samples_confidence.tsv]
    |
    v
AssociationAnalysisStage._process()
    sample_weights.load_sample_weights(filepath, vcf_samples_list)
    -> weight_vector: np.ndarray(n_samples,)   <- float in [0.0, 1.0]
    |
    v [in gene augmentation loop, before engine.run_all()]
    # Apply sample mask same as phenotype_vector masking (lines 2657-2678):
    if not all(sample_mask):
        sw = weight_vector[mask_arr]
    else:
        sw = weight_vector
    gene_data["sample_weights"] = sw
    |
    v
engine.run_all(gene_burden_data)
    -> for each test: test.run(gene, gene_data, config)
    |
    v [inside SKAT/burden test.run() methods]
    if self._null_model is None:
        self._null_model = backend.fit_null_model(
            phenotype=phenotype,
            covariates=covariate_matrix,
            trait_type=config.trait_type,
            sample_weights=gene_data.get("sample_weights"),   <- NEW
        )
    |
    v [inside PythonSKATBackend.fit_null_model()]
    glm = sm.GLM(phenotype, design_matrix, family=family)
    fit_result = glm.fit(var_weights=sample_weights)
    # residuals from weighted fit incorporate sample uncertainty
    |
    v [inside _compute_eigenvalues_filtered()]
    # Binary: phi = sqrt(mu_hat * (1-mu_hat) * sample_weights)
    # instead of phi = sqrt(mu_hat * (1-mu_hat))
```

Fisher test: ignores `sample_weights` (count-based, no regression framework).
Burden/SKAT/COAST tests: use weighted residuals from GLM fit.

**COAST lazy null model:** `PurePythonCOASTTest` fits its null model lazily on the
first gene call (allelic_series_python.py lines 383-391). The
`backend._skat_backend.fit_null_model()` call at line 388 would receive
`sample_weights` from the first gene's `contingency_data`. Since the null model is
cohort-level (not gene-level), `sample_weights` must be a cohort-constant array
passed through `contingency_data` consistently.

---

## New vs Modified Components

### New Files

| File | Purpose |
|------|---------|
| `variantcentrifuge/association/sample_weights.py` | Load + align sample confidence weights (mirrors covariates.py pattern) |
| `scoring/coast_classification/formula_config.json` | Variant-level BMV/DMV/PTV classification scoring model |
| `scoring/coast_classification/variable_assignment_config.json` | Column mappings for EFFECT, IMPACT, SIFT, PolyPhen |

### Modified Files

| File | What Changes | Scope |
|------|-------------|-------|
| `stages/analysis_stages.py` | COAST genotype matrix fix (Path C fallback reconstruct); sample weight loading and gene_data injection | Medium |
| `pipeline.py` | Import + wire `PCAComputationStage` when `pca_tool` configured | Small (4 lines) |
| `stages/processing_stages.py` | `PCAComputationStage`: use prefiltered VCF via fallback chain; optional `RegionRestrictionStage` | Small |
| `association/base.py` | `AssociationConfig` new fields: `case_confidence_file`, `fdr_weights`, `use_sparse_geno` | Small |
| `association/correction.py` | Add `apply_weighted_correction()` function | Small |
| `association/engine.py` | Route to weighted correction when `fdr_weights` configured; pass gene weights | Small |
| `association/backends/python_backend.py` | `fit_null_model()` sample_weights param; weighted phi in `_compute_eigenvalues_filtered()`; single eigendecomp in `_skato_get_pvalue()`; sparse input handling | Medium-Large |
| `association/genotype_matrix.py` | Optional sparse output path when `return_sparse=True` | Medium |
| `cli.py` | `--restrict-regions`, `--case-confidence-file`, `--pca-tool` (if not already wired) | Small |
| `filters.py` | No changes needed | None |

---

## Architecture Patterns to Follow

### Pattern 1: Config Field Expansion

All new feature flags go into `AssociationConfig` as optional fields with defaults
that preserve backward compatibility. The `_build_assoc_config_from_context()` helper
in `analysis_stages.py` reads from `context.config` — no other call sites need
changing when a new field is added to the dataclass.

### Pattern 2: Gene Data Dict as Feature Bag

The `gene_data` dict passed to each test's `run()` method carries all per-gene data.
New data (sample weights, sparse flags, annotation arrays) follows the same pattern:
store in the dict during the augmentation loop (lines 2634-2751 in analysis_stages.py),
read in `run()` with `.get()` and a sensible fallback (usually `None`).

### Pattern 3: Skip Guards Returning p_value=None

Every test method begins with a chain of skip guards:

```python
if "genotype_matrix" not in contingency_data:
    return TestResult(..., extra={"xxx_skip_reason": "REASON"})
```

New failure modes must add a guard here rather than raise. The engine handles
`p_value=None` gracefully (excludes from ACAT-O and FDR pass).

### Pattern 4: Leaf Module Isolation

`correction.py`, `covariates.py`, and `sample_weights.py` must import only
stdlib/numpy/pandas/statsmodels. No imports from other `variantcentrifuge` modules.
This keeps these modules testable in isolation.

### Pattern 5: Lazy Null Model Caching

`PurePythonCOASTTest` (and `PurePythonSKATTest`) fit the null model lazily on the
first gene and cache it in `self._null_model`. This is the established pattern for
cohort-level fit + per-gene test. New backend parameters (sample_weights) must be
consistent across all gene calls — the null model is fit once from the first gene's
`contingency_data` and reused.

---

## Anti-Patterns to Avoid

### Anti-Pattern 1: Reaching Into variants_df as Genotype Source

The current Path B fallback (analysis_stages.py line 2473-2486) reaches into
`context.variants_df` looking for per-sample GT columns. This creates an implicit
ordering dependency between stages. The correct fix (Path C) reconstructs per-sample
columns from the packed GT string directly — no dependency on `variants_df` state.

### Anti-Pattern 2: Gene-Level Association Results in Per-Variant Scoring

COAST p-values and n_bmv/n_dmv/n_ptv counts are gene-level outputs from the
association engine. Do not propagate them back into the per-variant DataFrame via the
scoring system. The separation between `.association.tsv` (gene-level) and the main
variant TSV (variant-level) is intentional and correct.

### Anti-Pattern 3: Null Model Re-Fit Per Gene

`fit_null_model()` is an O(n_samples^2) operation (GLM fitting). It must be called
once per cohort and cached. The lazy caching pattern in `PurePythonCOASTTest` and
`PurePythonSKATTest` is correct. New tests that call the backend must replicate this
pattern, not call `fit_null_model()` inside the per-gene loop.

### Anti-Pattern 4: Parallel Workers With Rpy2

`AssociationEngine.run_all()` already guards: `all_parallel_safe` check before
`ProcessPoolExecutor`. New tests using only numpy/scipy must declare
`parallel_safe = True`. Any test using rpy2 must declare `parallel_safe = False`
(or it will cause segfaults in worker processes).

### Anti-Pattern 5: Storing Full Genotype Matrix in Context

`gene_data["genotype_matrix"]` is per-gene and discarded after each test call. The
full multi-gene genotype tensor is never stored. For a 5000-sample cohort with 20000
genes, storing all genotype matrices would require ~80 GB. The current design correctly
builds the matrix per-gene and lets it be garbage collected.

---

## Build Order (Suggested Phase Sequence)

### Phase 1 — COAST Fix (prerequisite for all COAST validation)

- Fix genotype matrix reconstruction fallback (Path C) in `AssociationAnalysisStage`
- Integration test: COAST runs when `GeneBurdenAnalysisStage` runs first, confirming
  the `genotype_matrix` key is present in `gene_data`
- No new files; one function change in `analysis_stages.py`

**Blocks:** Any COAST-dependent feature validation (case-confidence weights with
COAST, COAST scoring model testing with the association path).

### Phase 2 — PCAComputationStage Wiring (infrastructure prerequisite)

- Wire `PCAComputationStage` into `pipeline.py` (4-line change)
- Fix `PCAComputationStage` to use prefiltered VCF
- Add `--pca-tool akt` CLI argument if not present
- Test: `--pca-tool akt --perform-association` runs without error

**Standalone:** No dependency on Phase 1.

### Phase 3 — Gene-Level FDR Weighting

- `apply_weighted_correction()` in `correction.py`
- `AssociationConfig.fdr_weights` field
- Engine routing to weighted function
- Tests: weighted correction produces expected ranking change vs unweighted

**Standalone:** No dependency on Phase 1 or 2.

### Phase 4 — Case-Confidence Weights

- `sample_weights.py` (new file: load + align)
- `AssociationConfig.case_confidence_file` field
- `fit_null_model()` `sample_weights` parameter
- `_compute_eigenvalues_filtered()` and `_test_skato()` weighted phi adjustments
- Tests: weighted null model with uniform weights produces identical results to
  unweighted (numerical regression test)

**Depends on:** Phase 1 (to validate COAST integration with weights).

### Phase 5 — Region Restriction

- `--restrict-regions` CLI arg
- Extend `GeneBedCreationStage` with bedtools intersect call
- Test: intersected BED reduces variant count as expected

**Standalone:** No dependency on other phases.

### Phase 6 — Sparse Genotype Matrices

- `genotype_matrix.py` sparse construction option
- `python_backend.py` sparse-aware matmuls (score_vec computation)
- `AssociationConfig.use_sparse_geno` flag
- All sample mask application sites audited for sparse indexing correctness
- Tests: sparse and dense paths produce identical p-values on same input

**High risk due to many call site changes. Do last.**

### Phase 7 — Single Eigendecomposition Optimization

- Pure internal refactor of `_skato_get_pvalue()` in `python_backend.py`
- No API changes, no new files, no config changes
- Tests: p-values before and after refactor are within 1e-10 relative tolerance

**Can be done independently at any point.**

### Phase 8 — COAST Classification Scoring Model

- New `scoring/coast_classification/` config files
- No Python changes required
- Tests: `--scoring-config scoring/coast_classification/` adds `coast_variant_class`
  column to variant output

**Can be developed in parallel with any phase.**

---

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| COAST genotype matrix bug root cause | HIGH | Read exact code paths; failure mechanism confirmed in analysis_stages.py lines 2462-2486 and 2630 |
| PCAComputationStage not wired | HIGH | Confirmed absent from pipeline.py imports and build_pipeline_stages() |
| Sample weight statsmodels integration | HIGH | Read fit_null_model(); sm.GLM.fit(var_weights=...) is documented statsmodels API |
| Sparse matrix impact on backend ops | HIGH | Read all kernel computations; assessed by inspection of matmul patterns |
| Single eigendecomp optimization | HIGH | Mathematical structure confirmed by reading _skato_get_pvalue() loop at lines 447-466 |
| FDR weighting approach | HIGH | apply_correction() signature confirmed; IHW weighted BH is established method |
| restrict-regions integration | HIGH | extract_variants() bcftools call confirmed; -R flag and bedtools intersect approach verified |
| Weighted phi in SKAT eigenvalue computation | MEDIUM | Mathematical derivation requires validation against R SKAT reference with sample weights |
| Sparse array indexing correctness at all call sites | MEDIUM | Must audit all numpy boolean indexing in augmentation loop |
