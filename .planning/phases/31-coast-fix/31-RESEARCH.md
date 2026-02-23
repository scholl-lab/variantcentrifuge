# Phase 31: COAST Fix - Research

**Researched:** 2026-02-23
**Domain:** COAST allelic series association test — bug fixes and configurability
**Confidence:** HIGH

## Summary

Phase 31 fixes five distinct problems in the existing COAST implementation and adds two new features (auto-field injection and configurable classification models). All work is internal to the association testing subsystem: `variantcentrifuge/association/tests/allelic_series*.py`, `variantcentrifuge/association/backends/coast_python.py`, and `variantcentrifuge/stages/analysis_stages.py`.

The codebase already has the entire scaffold in place. The bugs are small and localized. The new features (COAST-02, COAST-05, COAST-06, COAST-07) require the most new code but reuse the existing scoring formula engine (`scoring.py`) and CLI pattern from `--coast-backend`.

**Primary recommendation:** Fix the five bugs first (COAST-01 through COAST-04), then add the three new features (COAST-02, COAST-05/COAST-06/COAST-07) as independent units.

## Bug Inventory

### COAST-01: Genotype Matrix Not Available to COAST

**Root cause — located in `analysis_stages.py` lines 2364-2388 and 2530-2532:**

```python
# Line 2364: branch 1 — per-sample GT columns exist, GT not yet reconstructed
if "GT" not in df.columns and context.vcf_samples:
    gt_cols = _find_per_sample_gt_columns(df)
    if gt_cols:
        if needs_regression:
            df_with_per_sample_gt = df  # <-- saves df before reconstruction
        df = reconstruct_gt_column(df.copy(), context.vcf_samples)

# Line 2375: branch 2 — GT already reconstructed (e.g. by prior gene_burden stage)
elif "GT" in df.columns and needs_regression and context.vcf_samples:
    fallback_df = context.variants_df
    if fallback_df is not None:
        gt_cols_fb = _find_per_sample_gt_columns(fallback_df)
        if gt_cols_fb:
            df_with_per_sample_gt = fallback_df  # <-- saves variants_df

# Line 2530: genotype matrix construction
gt_source_df = df_with_per_sample_gt if df_with_per_sample_gt is not None else df
gt_columns_for_matrix = _find_gt_columns(gt_source_df)
if needs_regression and gt_columns_for_matrix and vcf_samples_list:
    # build genotype matrix...
```

**The failure mode:** When `df` already has GT reconstructed (branch 2 applies) but `context.variants_df` is None or has no per-sample GT columns, `df_with_per_sample_gt` stays None. Then `gt_source_df = df` — but `df` has had its per-sample GT columns replaced by the packed GT string column by `reconstruct_gt_column`. So `_find_gt_columns(df)` finds zero columns, the `if needs_regression and gt_columns_for_matrix` guard fires False, and no genotype matrix is ever added to `gene_data`. COAST then returns `NO_GENOTYPE_MATRIX` for every gene.

**Fix strategy:** The fallback in branch 2 must also check `df` itself before giving up. When `context.variants_df` has no GT columns but `df` still has per-sample GT columns (before any reconstruction), save `df` as the source. The correct guard should be:

```python
elif "GT" in df.columns and needs_regression and context.vcf_samples:
    # Check variants_df first, then fall back to df itself
    fallback_df = context.variants_df
    gt_cols_fb = _find_per_sample_gt_columns(fallback_df) if fallback_df is not None else []
    if gt_cols_fb:
        df_with_per_sample_gt = fallback_df
    else:
        # df itself may still have per-sample columns if GT reconstruction hasn't run yet
        gt_cols_df = _find_per_sample_gt_columns(df)
        if gt_cols_df:
            df_with_per_sample_gt = df
```

Additionally, the comment in the per-gene loop (line 2523) still says "logistic_burden, linear_burden, skat, skat_python" — it should include "coast".

### COAST-03: Partial-Category Fallback

**Current behavior (both `allelic_series.py` and `allelic_series_python.py`):**

```python
# Lines 551-581 (R backend) / 342-373 (Python backend):
if n_bmv == 0 or n_dmv == 0 or n_ptv == 0:
    # ... return TestResult(p_value=None, extra={"coast_skip_reason": "MISSING_CATEGORIES:..."})
```

**R reference behavior (from `insitro/AllelicSeries` source):**

The `Aggregator()` function uses `drop_empty=TRUE` (default) to remove zero-sum columns before regression. This means:
- 1 or 2 missing categories: drop those columns, run COAST on remaining, return valid p-value
- All categories empty: return NA (which maps to our `p_value=None` / `coast_skip_reason`)

The impact is enormous: in real cohorts most genes lack all 3 categories. Currently they all get `p_value=None`. After the fix, genes with 1 or 2 categories will produce valid (reduced-power) p-values.

**Fix in `coast_python.py` — `_aggregate_by_category`:** After building `cat_matrix`, the function already fills zero columns for absent categories. The downstream burden test functions need to handle rank-deficient matrices (already tested). No change needed in `_aggregate_by_category` itself.

**Fix in `allelic_series_python.py` (and `allelic_series.py`) — `run()` method:**

Replace the "all 3 required" guard with a "proceed with present categories; skip only if all empty" guard:

```python
# Count variants per category
n_bmv = int(np.sum(anno_filtered == 1))
n_dmv = int(np.sum(anno_filtered == 2))
n_ptv = int(np.sum(anno_filtered == 3))

missing = [cat for cat, n in [("BMV", n_bmv), ("DMV", n_dmv), ("PTV", n_ptv)] if n == 0]
present = [cat for cat, n in [("BMV", n_bmv), ("DMV", n_dmv), ("PTV", n_ptv)] if n > 0]

if len(present) == 0:
    # All categories empty — skip (existing behavior)
    return TestResult(p_value=None, extra={"coast_skip_reason": "NO_CLASSIFIABLE_VARIANTS", ...})

# Run COAST with whatever categories are present (drop_empty equivalent)
coast_status = "complete" if len(missing) == 0 else "partial"
if missing:
    logger.debug(f"COAST [{gene}]: partial — missing {', '.join(missing)} ...")
```

The `_aggregate_by_category` already handles missing categories by leaving those columns as zeros. The R equivalent of "drop_empty" maps to the existing zero-column behavior in burden tests — rank-deficient matrices already fall through gracefully via the `try/except` in `_run_burden_test`.

**Extra dict additions:**
```python
extra = {
    "coast_n_bmv": n_bmv,
    "coast_n_dmv": n_dmv,
    "coast_n_ptv": n_ptv,
    "coast_status": coast_status,  # "complete" or "partial"
}
if missing:
    extra["coast_missing_categories"] = ",".join(missing)
```

**INFO-level summary:** The `finalize()` method needs counters for complete/partial/skipped:

```python
logger.info(
    f"COAST: {n_complete + n_partial} genes tested "
    f"({n_complete} complete, {n_partial} partial, {n_skipped} skipped)"
)
```

### COAST-04: SnpEff "&"-Concatenated Effect Strings

**Current bug in `classify_variants()` (`allelic_series.py` lines 193-195):**

```python
is_ptv_effect = effect_series.str.strip().isin(PTV_EFFECTS)
is_missense = effect_series.str.strip() == MISSENSE_EFFECT
```

SnpEff emits multi-transcript effects as `"stop_gained&splice_region_variant"` when a variant affects multiple transcripts. The `isin()` and `==` checks fail on these concatenated strings.

**Fix strategy:** Split on `"&"` and check if ANY part matches. The decision (per CONTEXT.md) is to use split-and-prioritize: split on `"&"`, take the highest-priority effect.

Priority order (highest first):
1. PTV effects: stop_gained, frameshift_variant, splice_acceptor_variant, splice_donor_variant
2. missense_variant
3. Everything else

```python
def _resolve_effect(effect_str: str) -> str:
    """Resolve '&'-concatenated SnpEff effect string to single highest-priority effect."""
    parts = [p.strip() for p in effect_str.split("&")]
    # PTV takes priority
    for part in parts:
        if part in PTV_EFFECTS:
            return part
    # Missense next
    if MISSENSE_EFFECT in parts:
        return MISSENSE_EFFECT
    # Return first part as fallback
    return parts[0] if parts else effect_str
```

Then in `classify_variants()`:
```python
effect_series = gene_df[effect_col].astype(object).fillna("").astype(str).map(_resolve_effect)
```

This ensures `"stop_gained&splice_region_variant"` resolves to `"stop_gained"` (PTV), while `"missense_variant&splice_region_variant"` resolves to `"missense_variant"` (missense, then classified by SIFT/PolyPhen).

Logging at DEBUG per the decisions.

## New Features

### COAST-02: Auto-Field Injection

**Where injection happens:** In `analysis_stages.py` `GeneBurdenAnalysisStage.run()` **before** `FieldExtractionStage` runs — the injection must happen at the config level so SnpSift/bcftools extracts the fields.

**The challenge:** The fields are injected into `context.config["fields_to_extract"]` which has already been set. The injection must happen before `FieldExtractionStage` executes, i.e., during pipeline setup or in an early stage.

**Recommended approach:** Inject in a new early-pipeline step or in `GeneBurdenAnalysisStage` as a pre-check that modifies `context.config["fields_to_extract"]` if COAST is in the test list. This requires that `GeneBurdenAnalysisStage` runs BEFORE `FieldExtractionStage`. Looking at the pipeline order, the analysis stage runs after extraction.

**Correct injection point:** In `cli.py`'s `main()` function, after `fields` is resolved but before the pipeline runs. Alternatively, as a new `CoastFieldInjectionStage` that runs before `FieldExtractionStage`.

**Fields needed by the default (`sift_polyphen`) model:**
- `ANN[0].EFFECT` (EFFECT)
- `ANN[0].IMPACT` (IMPACT)
- `dbNSFP_SIFT_pred`
- `dbNSFP_Polyphen2_HDIV_pred`

**Fields needed by a CADD-based model:**
- `ANN[0].EFFECT` (EFFECT)
- `ANN[0].IMPACT` (IMPACT)
- `dbNSFP_CADD_phred`

The `variable_assignment_config.json` maps VCF field names (left side) to internal variable names. The keys on the left are the field names to inject.

**Injection pattern** (reusing `ensure_fields_in_extract`):

```python
# In cli.py or a setup stage:
if "coast" in test_names:
    classification_model = getattr(args, "coast_classification", "sift_polyphen")
    model_dir = scoring_base_dir / "coast_classification" / classification_model
    var_config = read_variable_assignment_config(model_dir)
    required_fields = list(var_config.get("variables", {}).keys())
    # Check VCF header for missing fields
    vcf_fields = get_vcf_info_fields(vcf_file)  # existing utility
    missing_vcf_fields = [f for f in required_fields if f not in vcf_fields]
    if missing_vcf_fields:
        raise SystemExit(f"COAST classification model '{classification_model}' requires "
                         f"fields not in VCF: {missing_vcf_fields}")
    fields = ensure_fields_in_extract(fields, required_fields)
    logger.info(f"COAST auto-injected fields: {required_fields}")
```

**Where VCF field checking exists:** `variantcentrifuge/vcf_utils.py` or similar — check for `parse_vcf_header` usage in CLI.

### COAST-05/COAST-06: Classification Model Selection and Config

**New CLI option:** `--coast-classification <model-name>` (default: `sift_polyphen`)

Pattern: identical to `--coast-backend python/r/auto`. Add to `argparse` in `cli.py`, propagate to `cfg["coast_classification"]`, propagate to `AssociationConfig`.

**New config directory structure:**

```
scoring/coast_classification/
├── sift_polyphen/          # Default model (migrated from hardcoded classify_variants)
│   ├── variable_assignment_config.json
│   └── formula_config.json
├── cadd/                   # CADD-threshold model
│   ├── variable_assignment_config.json
│   └── formula_config.json
└── custom_template/        # Empty template with documentation comments
    ├── variable_assignment_config.json
    └── formula_config.json
```

**Formula engine reuse:** The existing `read_scoring_config()` + `apply_scoring()` in `scoring.py` works on per-row DataFrames. For COAST classification, the formula outputs a single integer column (0/1/2/3) per variant row.

**sift_polyphen formula_config.json example:**

```json
{
  "output_scores": ["coast_category"],
  "formulas": [
    {
      "is_ptv": "(impact == 'HIGH') & effect.isin(['stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant'])"
    },
    {
      "is_missense": "effect == 'missense_variant'"
    },
    {
      "sift_dam": "sift_pred.str.contains('deleterious|^D$', na=False, regex=True)"
    },
    {
      "polyphen_dam": "polyphen_pred.str.contains('probably_damaging|possibly_damaging|^D$|^P$', na=False, regex=True)"
    },
    {
      "sift_ben": "sift_pred.str.contains('tolerated|^T$', na=False, regex=True)"
    },
    {
      "polyphen_ben": "polyphen_pred.str.contains('benign|^B$', na=False, regex=True)"
    },
    {
      "is_dmv": "is_missense & (sift_dam | polyphen_dam)"
    },
    {
      "is_bmv": "is_missense & sift_ben & polyphen_ben"
    },
    {
      "coast_category": "is_ptv * 3 + is_dmv * 2 + is_bmv * 1"
    }
  ]
}
```

Note: The formula engine uses `pandas.eval` with `engine="python"` and allows series method calls.

**sift_polyphen variable_assignment_config.json:**

```json
{
  "variables": {
    "ANN[0].EFFECT": "effect|default:''",
    "ANN[0].IMPACT": "impact|default:''",
    "dbNSFP_SIFT_pred": "sift_pred|default:''",
    "dbNSFP_Polyphen2_HDIV_pred": "polyphen_pred|default:''"
  }
}
```

**cadd variable_assignment_config.json:**

```json
{
  "variables": {
    "ANN[0].EFFECT": "effect|default:''",
    "ANN[0].IMPACT": "impact|default:''",
    "dbNSFP_CADD_phred": "cadd_phred|default:0.0"
  }
}
```

**cadd formula_config.json (CADD thresholds from literature):**
- BMV: missense AND 0 < CADD_phred < 15 (below damaging threshold)
- DMV: missense AND CADD_phred >= 15 (Kircher 2014 threshold)
- PTV: HIGH impact + PTV effects (same as default)

```json
{
  "output_scores": ["coast_category"],
  "formulas": [
    {
      "is_ptv": "(impact == 'HIGH') & effect.isin(['stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant'])"
    },
    {
      "is_missense": "effect == 'missense_variant'"
    },
    {
      "is_dmv": "is_missense & (cadd_phred >= 15)"
    },
    {
      "is_bmv": "is_missense & (cadd_phred > 0) & (cadd_phred < 15)"
    },
    {
      "coast_category": "is_ptv * 3 + is_dmv * 2 + is_bmv * 1"
    }
  ]
}
```

**Integration in `classify_variants()`:** The function needs to accept an optional `model_dir` parameter. When provided, it calls `read_scoring_config(model_dir)` + `apply_scoring(gene_df, config)` and reads the `coast_category` output column. When not provided, it uses the hardcoded logic (backward compat) or defaults to the `sift_polyphen` model.

**Formula engine limitation:** `apply_scoring()` renames columns for evaluation and renames them back. It expects all input columns to be in the DataFrame. Since we're classifying per variant (gene_df already has all annotation columns), this works naturally. The output column `coast_category` will be an integer 0/1/2/3.

### COAST-07: Multi-Transcript Resolution Strategy

Per CONTEXT.md, this is handled by the `_resolve_effect()` function in COAST-04 above. The model config can declare a `resolution_strategy` key (default: `majority_vote` or `priority_order`). For simplicity and correctness, priority-order (PTV > missense > other) is recommended since it matches how SnpEff annotation priority works.

The formula configs handle this implicitly: by using `isin()` on the pre-resolved effect string, if the classification model config itself calls `_resolve_effect()` before classification. Since formula eval uses pandas operations on strings, the `&`-split resolution must happen before passing to the formula engine.

### COAST Diagnostics: `coast_classification.tsv`

When `--diagnostics-output` is set, write a per-variant TSV audit file.

**Column order:** `gene`, `variant_id`, `original_effect`, `resolved_effect`, `assigned_category`, `score_used`

Where `variant_id` = `CHROM:POS:REF:ALT` assembled from columns or a `VARIANT_ID` column if present.
`score_used` = name of classification model used (e.g., `"sift_polyphen"`).

## Architecture Patterns

### Existing Pattern: Scoring Formula Engine

The scoring engine (`scoring.py`) runs `pandas.eval()` with `engine="python"` on a per-row DataFrame. It:
1. Reads `variable_assignment_config.json` (column mapping + defaults)
2. Renames columns to variable names
3. Evaluates each formula in sequence (results become new columns)
4. Renames columns back
5. Returns final output column(s)

COAST classification reuses this with one output column: `coast_category` (int 0/1/2/3).

**Constraint:** The formula engine uses `scored_df.eval()`, which means string operations need `.str.` method calls, not direct operations. The existing formulas in `nephro_candidate_score` use `str.contains()` which is valid.

**Constraint:** The `output_scores` list in `formula_config.json` controls which intermediate columns are kept. For COAST classification, only `coast_category` should be in `output_scores` — all intermediate boolean columns (`is_ptv`, `is_missense`, etc.) will be dropped.

### Existing Pattern: CLI Option Propagation

New CLI option → `cfg["key"]` in `cli.py` → `AssociationConfig.field` → accessible in test via `getattr(config, "field", default)`.

Existing examples:
- `--coast-backend` → `cfg["coast_backend"]` → `AssociationConfig.coast_backend`
- `--coast-weights` → `cfg["coast_weights"]` → `AssociationConfig.coast_weights`

New option follows same pattern:
- `--coast-classification` → `cfg["coast_classification"]` → `AssociationConfig.coast_classification`

### Existing Pattern: TestResult Extra Dict

All COAST-specific output goes into `extra` dict on `TestResult`. The engine writes all extra keys directly to the output DataFrame (engine.py lines 498-499):

```python
for extra_key, extra_val in res.extra.items():
    row[extra_key] = extra_val
```

New extra keys (`coast_status`, `coast_missing_categories`) will automatically appear in association output.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Field presence check in fields_to_extract | Custom string parser | `ensure_fields_in_extract()` in `utils.py` | Already handles deduplication |
| Scoring formula evaluation | Custom eval engine | `read_scoring_config()` + `apply_scoring()` in `scoring.py` | Already handles column rename, defaults, error handling |
| Argparse CLI option | Manual sys.argv | Follow existing `--coast-backend` pattern in `cli.py` | Consistent with all existing association options |
| VCF header field list | Shell subprocess | `parse_vcf_header()` in `variantcentrifuge/cli.py` or equivalent | Already used for `--show-vcf-annotations` |
| Building genotype matrix | Custom numpy loop | `build_genotype_matrix()` in `association/genotype_matrix.py` | Handles missing, imputation, multi-allelic |

## Common Pitfalls

### Pitfall 1: Formula Engine Categorical Columns

**What goes wrong:** When `gene_df` has EFFECT as a pandas Categorical column (memory optimization), `pandas.eval()` with string operations like `.isin()` can fail or behave unexpectedly.

**How to avoid:** The existing `classify_variants()` already handles this: `effect_series = gene_df[effect_col].astype(object).fillna("").astype(str)`. The formula configs should document that inputs may be Categorical.

**Alternative:** In the scoring integration, pre-cast string columns to `object` dtype before calling `apply_scoring()`.

### Pitfall 2: Partial-Category Burden Tests with Zero-Variance Predictors

**What goes wrong:** When only PTV category is present, the baseline 3-df test gets a 3-column predictor where columns 0 and 1 (BMV, DMV) are all zeros. Logistic regression will fail to converge or produce degenerate F-tests.

**How to avoid:** The `_run_burden_test()` function wraps fits in `try/except` and returns `None` on failure. The Cauchy combination in `cauchy_combination()` skips `None` values. So degenerate columns produce `None` p-values which are silently excluded from the omnibus — this is correct behavior. **No additional code needed.**

**Warning sign:** Very high omnibus p-values when only 1 category present — this is expected (reduced power).

### Pitfall 3: Genotype Matrix / gene_df Dimension Mismatch After Partial Classification

**What goes wrong:** When `build_genotype_matrix()` applies the missing-site filter and drops some variants, `geno.shape[1]` < `len(gene_df)`. The current `ANNOTATION_GENOTYPE_MISMATCH` check skips the gene in this case.

**Note:** This existing check is correct and should remain. The partial-category fix does NOT touch this check. The dimension mismatch is a separate issue from missing categories.

### Pitfall 4: Formula Config with Missing Columns Causes Silent Default

**What goes wrong:** When a user requests the CADD classification model but their VCF has no CADD annotations, `apply_scoring()` creates the column with `default: 0.0`. All variants get CADD=0, which means all missense variants become BMV (CADD < 15). This silently produces wrong results.

**How to avoid:** The auto-field injection step (COAST-02) must check VCF header fields BEFORE pipeline runs. If CADD fields are missing from the VCF, fail early with a clear error: `"COAST classification model 'cadd' requires 'dbNSFP_CADD_phred' which is not in the VCF. Add CADD annotations or use --coast-classification sift_polyphen"`.

**Implementation:** Use `parse_vcf_header()` (already in CLI for `--show-vcf-annotations`) to get available fields, cross-reference against `variable_assignment_config.json` keys.

### Pitfall 5: apply_scoring() Output Column Dtype

**What goes wrong:** The formula engine wraps results: `pd.to_numeric(scored_df[score_name], errors="coerce")`. The `coast_category` column becomes float64 (0.0, 1.0, 2.0, 3.0). The existing code uses `anno_filtered == 1` (integer comparison), which works with float equality but is fragile.

**How to avoid:** After calling `apply_scoring()`, cast the output column: `anno_codes = gene_df["coast_category"].fillna(0).astype(int).to_numpy()`.

## Code Examples

### Current Partial-Category Skip (to be replaced):

```python
# Source: variantcentrifuge/association/tests/allelic_series_python.py lines 342-373
if n_bmv == 0 or n_dmv == 0 or n_ptv == 0:
    missing = []
    if n_bmv == 0:
        missing.append("BMV")
    if n_dmv == 0:
        missing.append("DMV")
    if n_ptv == 0:
        missing.append("PTV")
    logger.debug(...)
    return TestResult(
        gene=gene, p_value=None,
        extra={"coast_skip_reason": f"MISSING_CATEGORIES:{','.join(missing)}", ...}
    )
```

### Replacement Partial-Category Logic:

```python
# After COAST-03 fix
missing = [cat for cat, n in [("BMV", n_bmv), ("DMV", n_dmv), ("PTV", n_ptv)] if n == 0]
present = [cat for cat, n in [("BMV", n_bmv), ("DMV", n_dmv), ("PTV", n_ptv)] if n > 0]

if len(present) == 0:
    return TestResult(p_value=None, extra={"coast_skip_reason": "ALL_CATEGORIES_EMPTY", ...})

coast_status = "complete" if not missing else "partial"
if missing:
    logger.debug(f"COAST [{gene}]: partial — missing {', '.join(missing)} "
                 f"(BMV={n_bmv}, DMV={n_dmv}, PTV={n_ptv})")

# Proceed to run COAST with present categories (drop_empty handled implicitly)
# extra dict additions:
extra_base = {
    "coast_n_bmv": n_bmv,
    "coast_n_dmv": n_dmv,
    "coast_n_ptv": n_ptv,
    "coast_status": coast_status,
}
if missing:
    extra_base["coast_missing_categories"] = ",".join(missing)
```

### SnpEff "&" Resolution:

```python
# Source: new function in allelic_series.py
def _resolve_effect(effect_str: str) -> str:
    """Resolve '&'-concatenated SnpEff effect string to single highest-priority effect."""
    if "&" not in effect_str:
        return effect_str  # fast path for common case
    parts = [p.strip() for p in effect_str.split("&")]
    for part in parts:
        if part in PTV_EFFECTS:
            return part
    if MISSENSE_EFFECT in parts:
        return MISSENSE_EFFECT
    return parts[0] if parts else effect_str
```

### Finalize Summary Counter Pattern:

```python
# In prepare():
self._n_complete = 0
self._n_partial = 0
self._n_skipped = 0

# In run() — after successful test:
if coast_status == "complete":
    self._n_complete += 1
else:
    self._n_partial += 1

# In run() — on skip:
self._n_skipped += 1

# In finalize():
n_tested = self._n_complete + self._n_partial
logger.info(
    f"COAST: {n_tested} genes tested "
    f"({self._n_complete} complete, {self._n_partial} partial, "
    f"{self._n_skipped} skipped)"
)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| R AllelicSeries via rpy2 | Pure Python COAST backend (default) | Phase 24 | Thread-safe, no R dependency |
| Hardcoded classify_variants() | Formula engine config (Phase 31) | This phase | Configurable classification |
| Skip on missing categories | Partial-category fallback (Phase 31) | This phase | Valid p-values for real cohorts |

**Deprecated (in this phase):** The hardcoded `classify_variants()` logic in `allelic_series.py` is migrated to formula config. The function itself is kept for backward compatibility but delegates to the formula engine when a model config is available.

## Open Questions

1. **FieldExtractionStage timing for auto-injection**
   - What we know: `fields_to_extract` is set in `cli.py` and passed to the extraction stage
   - What's unclear: Whether injecting in `cli.py` vs. a pre-extraction stage is better
   - Recommendation: Inject in `cli.py` `main()` after `fields` is resolved (before `validate_mandatory_parameters`). This is the simplest path — no new stage needed.

2. **Formula engine &-string handling in CADD model**
   - What we know: The formula engine evaluates on the raw `effect` column, which may contain `"stop_gained&splice_region_variant"`
   - What's unclear: Whether formulas using `.isin()` or `==` will fail on &-strings
   - Recommendation: Pre-resolve `effect` column with `_resolve_effect()` before passing gene_df to `apply_scoring()`, so formula configs don't need to handle concatenated strings.

3. **diagnostics_output directory vs. file path**
   - What we know: `diagnostics_output` in `AssociationConfig` is currently a path to a directory
   - What's unclear: Whether to write `coast_classification.tsv` inside that directory or to a fixed path
   - Recommendation: Write to `{diagnostics_output}/coast_classification.tsv` (consistent with existing diagnostics pattern).

## Sources

### Primary (HIGH confidence)
- Direct codebase inspection: `variantcentrifuge/association/tests/allelic_series.py` — classify_variants(), COASTTest.run()
- Direct codebase inspection: `variantcentrifuge/association/tests/allelic_series_python.py` — PurePythonCOASTTest.run()
- Direct codebase inspection: `variantcentrifuge/association/backends/coast_python.py` — PythonCOASTBackend, _aggregate_by_category
- Direct codebase inspection: `variantcentrifuge/stages/analysis_stages.py` lines 2354-2653 — genotype matrix construction
- Direct codebase inspection: `variantcentrifuge/scoring.py` — formula engine
- Direct codebase inspection: `variantcentrifuge/association/base.py` — AssociationConfig, TestResult

### Secondary (MEDIUM confidence)
- WebFetch of insitro/AllelicSeries R source (`allelic_series.R`): Confirmed `Aggregator(drop_empty=TRUE)` behavior — empty categories dropped before regression, COAST returns NA only when all categories empty
- CRAN AllelicSeries documentation: Confirmed partial category handling intent

### Tertiary (LOW confidence)
- WebSearch for AllelicSeries behavior: confirmed partial category support exists but full source code not directly read

## Metadata

**Confidence breakdown:**
- Bug identification (COAST-01, COAST-03, COAST-04): HIGH — traced through source code
- Fix implementations: HIGH — patterns follow existing code exactly
- Formula engine config design (COAST-06): HIGH — reuses existing scoring engine
- Auto-injection implementation (COAST-02): MEDIUM — field injection point needs validation against full pipeline execution path
- R reference partial-category behavior: MEDIUM — confirmed via WebFetch but not verified against live R execution

**Research date:** 2026-02-23
**Valid until:** 2026-03-23 (30 days — stable codebase)
