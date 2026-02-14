---
phase: 09
plan: 02
subsystem: inheritance-analysis-vectorization
tags: [vectorization, numpy, performance, pattern-deduction]
requires:
  - "09-01: Golden file validation infrastructure"
  - "Phase 8: DataFrame optimization (column sanitization, encode_genotypes)"
provides:
  - "Vectorized Pass 1 pattern deduction (10-100x speedup expected)"
  - "NumPy boolean mask operations for inheritance patterns"
  - "Clinically equivalent output validated by golden files"
affects:
  - "09-03: Compound het vectorization (will complete Pass 2)"
  - "09-04: Full vectorization validation (will validate entire pipeline)"
tech-stack:
  added: []
  patterns:
    - "NumPy boolean mask operations for genetic pattern detection"
    - "Per-sample vectorized logic with genotype matrix encoding"
    - "Fallback pattern tracking per sample (carrier vs unknown)"
key-files:
  created:
    - variantcentrifuge/inheritance/vectorized_deducer.py
  modified:
    - variantcentrifuge/inheritance/analyzer.py
decisions:
  - decision: "Build genotype matrix once before all pattern checks"
    rationale: "Encode all sample genotypes (n_variants x n_samples) once using int8, reuse for all checks"
    alternatives: ["Encode per-sample", "Encode per-pattern-check"]
    impact: "Memory efficient (int8), eliminates redundant encoding, enables pure NumPy operations"
  - decision: "Per-sample pattern tracking with patterns_before_sample"
    rationale: "Original applies fallback per-sample (if sample added no patterns, add carrier/unknown)"
    alternatives: ["Global fallback after all samples", "No fallback tracking"]
    impact: "Matches original behavior exactly, adds carrier/unknown for samples with variant but no specific patterns"
  - decision: "Use isinstance(dtype, pd.CategoricalDtype) for categorical check"
    rationale: "pd.api.types.is_categorical_dtype() is deprecated in pandas 2.x"
    alternatives: ["Keep deprecated function", "Use hasattr check"]
    impact: "Future-proof, no deprecation warnings, follows pandas best practices"
  - decision: "Reuse encode_genotypes from comp_het_vectorized.py"
    rationale: "Already tested, handles missing/phased genotypes, returns int8 array"
    alternatives: ["Create new encoding function", "Use string operations"]
    impact: "Code reuse, consistent encoding across modules, proven correctness"
  - decision: "Clear variable names and genetic logic comments"
    rationale: "Per 09-CONTEXT.md readability requirement: de_novo_mask not m1, father_gts not fg"
    alternatives: ["Compact variable names", "Minimal comments"]
    impact: "Maintainable code, understandable genetic logic for future developers"
metrics:
  duration: "15 minutes"
  completed: "2026-02-14"
---

# Phase 09 Plan 02: Deducer Vectorization

**One-liner:** Replaced df.apply(deduce_patterns_for_variant) with vectorized_deduce_patterns using NumPy boolean masks for 10-100x speedup while preserving clinical correctness

## Objective

Vectorize Pass 1 (inheritance pattern deduction) by replacing per-row `df.apply()` calls with NumPy boolean mask operations. This is the highest-impact optimization in Phase 9, as Pass 1 is the primary bottleneck in inheritance analysis.

## What Was Built

### 1. Vectorized Deducer Module (`vectorized_deducer.py`)

**802-line module implementing vectorized pattern deduction:**

**Core function:**
```python
def vectorized_deduce_patterns(
    df: pd.DataFrame,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str]
) -> list[list[str]]
```

Returns list of pattern lists (one per variant), matching the exact output format of the original `deduce_patterns_for_variant`.

**Architecture:**

1. **Genotype Matrix Encoding** (once, upfront)
   - Build (n_variants x n_samples) int8 matrix
   - Encode ALL sample genotypes using `encode_genotypes()` from `comp_het_vectorized.py`
   - -1 = missing, 0 = ref (0/0), 1 = het (0/1), 2 = hom_alt (1/1)
   - Handle samples not in DataFrame (fill with -1)

2. **Single Sample / No Pedigree Handler**
   - If `not pedigree_data or len(sample_list) == 1`
   - Apply simple rules: ref → "reference", hom_alt → "homozygous", het → "unknown"
   - Return early (no genetic logic needed)

3. **Per-Sample Vectorized Pattern Deduction**
   - For each sample in sample_list:
     - Track patterns_before_sample (for fallback logic)
     - Run all pattern checks using NumPy boolean masks
     - Apply fallback if sample added no patterns

**Pattern Detection Functions:**

| Function | Pattern | Logic |
|----------|---------|-------|
| `_check_de_novo_vectorized` | De novo | `has_variant & (father==0) & (mother==0)` |
| `_check_dominant_vectorized` | Autosomal dominant | `has_variant & affected & (parent_variant & parent_affected)` |
| `_check_recessive_vectorized` | Autosomal recessive | `hom_alt & affected & (father_variant & mother_variant)` |
| `_check_xlr_male_vectorized` | X-linked recessive (male) | `x_chrom & has_variant & affected & mother_carrier & ~father_has` |
| `_check_xlr_female_vectorized` | X-linked recessive (female) | `x_chrom & hom_alt & affected & father_has & mother_has` |
| `_check_xld_vectorized` | X-linked dominant | `x_chrom & has_variant & affected & parent_variant & parent_affected` |
| `_check_mitochondrial_vectorized` | Mitochondrial | `mt_chrom & has_variant & (mother_has \| no_mother_data)` |
| `_apply_fallback_vectorized` | Unknown / Carrier | If sample added no patterns: affected → "unknown", else → "carrier" |

**Critical correctness features:**
- Use bitwise operators (`&`, `|`, `~`) on NumPy arrays, not logical (`and`, `or`, `not`)
- Use `np.isin()` for categorical-safe chromosome comparison
- Handle missing genotypes (-1) correctly (never treat as ref or variant)
- Per-sample fallback tracking (matches original behavior)
- Deduplication per variant (remove duplicates while preserving order)

### 2. Analyzer Integration

**Updated `analyzer.py` to use vectorized deducer:**

**Before (lines 86-88):**
```python
df["_inheritance_patterns"] = df.apply(
    lambda row: deduce_patterns_for_variant(row.to_dict(), pedigree_data, sample_list),
    axis=1
)
```

**After (line 86):**
```python
df["_inheritance_patterns"] = vectorized_deduce_patterns(df, pedigree_data, sample_list)
```

- Removed unused import of `deduce_patterns_for_variant` from `analyzer.py`
- Kept import in `deducer.py` (still used by `parallel_analyzer.py` — will be updated in Plan 04)

## Verification Results

All verification criteria met:

✅ `python scripts/validate_inheritance.py compare` — Exit code 0, all 10 scenarios pass
✅ `pytest tests/test_inheritance/ -v` — All 140 tests pass (5 skipped, no failures)
✅ `pytest -m unit` — All 599 unit tests pass (no regressions)
✅ `make ci-check` — All 1071 tests pass, linting clean, typecheck clean
✅ `grep "df.apply.*deduce_patterns" analyzer.py` — Returns no matches (old code removed)
✅ Module imports cleanly, passes linting, no deprecation warnings

## Technical Implementation

### Genotype Matrix Building

```python
def _build_genotype_matrix(df, sample_list):
    n_variants = len(df)
    n_samples = len(sample_list)
    gt_matrix = np.full((n_variants, n_samples), -1, dtype=np.int8)

    for i, sample_id in enumerate(sample_list):
        if sample_id in df.columns:
            gt_matrix[:, i] = encode_genotypes(df[sample_id])

    return gt_matrix, sample_to_idx
```

**Memory efficiency:** int8 encoding uses 1 byte per genotype (vs 3-7 bytes for strings)

### Per-Sample Fallback Logic

**Key insight:** Original applies fallback PER SAMPLE, not globally.

**Implementation:**
```python
# Before this sample's analysis
patterns_before_sample = [len(p) for p in patterns_per_variant]

# ... run all pattern checks for this sample ...

# After: check if THIS SAMPLE added any patterns
for idx in variant_indices_with_variant:
    patterns_added_by_sample = len(patterns_per_variant[idx]) - patterns_before_sample[idx]
    if patterns_added_by_sample == 0:
        # This sample has variant but added no patterns
        if affected:
            patterns_per_variant[idx].append("unknown")
        else:
            patterns_per_variant[idx].append("carrier")
```

This matches the original logic:
- Father with variant adds "autosomal_dominant" → no fallback
- Father with variant adds no patterns → fallback to "unknown" (if affected) or "carrier"

### X-linked Pattern Detection

**Male X-linked recessive:**
```python
male_xlr_base_mask = x_mask & has_variant_mask & affected
mother_has_variant_mask = mother_gts > 0
father_violates_mask = father_gts > 0  # Males can't pass X to sons

classic_xlr_mask = male_xlr_base_mask & mother_has_variant_mask & ~father_violates_mask
```

**Female X-linked recessive:**
```python
female_xlr_base_mask = x_mask & (sample_gts == 2) & affected  # hom_alt required
classic_xlr_mask = female_xlr_base_mask & (father_gts > 0) & (mother_gts > 0)
```

**X chromosome detection:**
```python
x_mask = np.isin(chrom_array, ["X", "CHRX", "23", "x", "chrX", "chrx"])
```

Handles all common X chromosome representations.

## Golden File Validation

All 10 scenarios produce clinically equivalent output:

| Scenario | Variants | Key Pattern | Result |
|----------|----------|-------------|--------|
| `trio_denovo` | 1 | De novo (child het, parents ref) | ✅ PASS |
| `trio_dominant` | 1 | Autosomal dominant (father affected) | ✅ PASS |
| `trio_recessive` | 1 | Autosomal recessive (child hom_alt) | ✅ PASS |
| `trio_denovo_candidate` | 1 | De novo candidate (missing GT) | ✅ PASS |
| `single_sample` | 2 | No pedigree (het → "unknown") | ✅ PASS |
| `extended_family` | 1 | Multi-generation dominant | ✅ PASS |
| `x_linked` | 1 | X-linked recessive (male) | ✅ PASS |
| `mitochondrial` | 1 | Mitochondrial (MT chromosome) | ✅ PASS |
| `compound_het` | 2 | Compound heterozygous (trans) | ✅ PASS |
| `edge_cases` | 2 | Missing GTs, single variants | ✅ PASS |

**Total:** 13 variants across 10 scenarios — all pass

## Deviations from Plan

**Auto-fixed issues (Deviation Rule 1 & 2):**

1. **Missing genotype handling in fallback**
   - **Found during:** Task 1 implementation
   - **Issue:** Original checks `not patterns` for this sample's patterns, not global patterns
   - **Fix:** Track `patterns_before_sample` per sample, compare before/after to determine fallback
   - **Files modified:** `vectorized_deducer.py`
   - **Commit:** feat(09-02): create vectorized Pass 1 pattern deduction

2. **Categorical dtype deprecation warning**
   - **Found during:** Task 2 testing (pytest showed 301 deprecation warnings)
   - **Issue:** `pd.api.types.is_categorical_dtype()` is deprecated in pandas 2.x
   - **Fix:** Changed to `isinstance(chrom_array.dtype, pd.CategoricalDtype)`
   - **Files modified:** `vectorized_deducer.py`
   - **Commit:** feat(09-02): wire vectorized Pass 1 into analyzer.py

3. **Code formatting (line length)**
   - **Found during:** make ci-check
   - **Issue:** Some NumPy boolean expressions exceeded 100 char line length
   - **Fix:** Applied `make format` (ruff auto-formatting)
   - **Files modified:** `vectorized_deducer.py`
   - **Commit:** Automatic via make format

## Challenges Encountered

1. **Per-sample fallback logic**
   - **Challenge:** Initial implementation applied fallback globally (if variant has NO patterns from ANY sample)
   - **Root cause:** Misread original logic — fallback is per-sample, not per-variant
   - **Solution:** Track `patterns_before_sample` and compare to detect if THIS SAMPLE added patterns
   - **Impact:** 4 golden file scenarios failed before fix, all pass after

2. **Chromosome comparison with categorical dtype**
   - **Challenge:** `chrom_array == "X"` fails on categorical dtype
   - **Root cause:** Phase 8 introduced categorical dtypes for low-cardinality columns
   - **Solution:** Use `np.isin(chrom_array, ["X", "CHRX", ...])` which handles categorical correctly
   - **Impact:** X-linked and mitochondrial tests would fail without this

## Performance Impact

**Expected speedup:** 10-100x for Pass 1 (deduction)

**Benchmarking deferred to Plan 05** (Benchmark Verification), but preliminary analysis:

- **Eliminated:** `df.apply()` with per-row lambda (Python interpreter overhead)
- **Replaced with:** NumPy vectorized operations (compiled C code)
- **Genotype encoding:** Once instead of per-variant-per-sample
- **Memory:** int8 matrix (1 byte per genotype) vs string Series

**Baseline (from STATE.md):**
- 10K variants: 3.2 seconds total inheritance analysis
- Pass 1 (deduction) is ~50-70% of this time (~1.6-2.2 seconds)

**Expected after vectorization:**
- 10K variants Pass 1: ~0.02-0.2 seconds (10-100x speedup)
- Total inheritance analysis: ~1.2-1.6 seconds (50-70% overall speedup)

Will be confirmed in Plan 05 benchmarks.

## Code Readability

Per 09-CONTEXT.md requirements, code uses:

**Clear variable names:**
- `de_novo_mask` not `m1`
- `father_gts` not `fg`
- `has_variant_mask` not `hv`
- `mother_has_variant_mask` not `mh`

**Genetic logic comments:**
```python
# Classic de novo: child has variant, both parents ref
de_novo_mask = has_variant_mask & (father_gts == 0) & (mother_gts == 0)

# Father check: males can't pass X to sons, so father must NOT have variant
father_violates_mask = father_gts > 0
```

**Step-by-step operations:**
```python
# Step 1: Build genotype matrix
gt_matrix, sample_to_idx = _build_genotype_matrix(df, sample_list)

# Step 2: Per-sample vectorized pattern deduction
for sample_id in sample_list:
    # a. De novo check
    # b. Dominant check
    # c. Recessive check
    # ... etc
```

## Test Coverage

- **Existing tests:** All 140 inheritance tests pass unchanged
- **New tests:** None (vectorization is implementation detail, interface unchanged)
- **Test markers:** Tests marked as `unit` and `inheritance`
- **Golden file validation:** 10 scenarios, all pass

## Next Phase Readiness

**Phase 9 Plan 03 (Compound Het Vectorization):**
- ✅ Golden files validated with vectorized Pass 1
- ✅ Pass 1 vectorization complete and tested
- ✅ Ready for Pass 2 (compound het) vectorization
- ✅ Can reuse golden files for validation

**Phase 9 Plan 04 (Full Vectorization Validation):**
- ✅ Vectorized deducer ready for integration
- ✅ Will update `parallel_analyzer.py` to use vectorized deducer
- ✅ Will run final golden file validation

**Phase 9 Plan 05 (Benchmark Verification):**
- ✅ Vectorized Pass 1 ready for performance measurement
- ✅ Expected 10-100x speedup on Pass 1
- ✅ Expected 50-70% overall inheritance analysis speedup

## Files Changed

**Created:**
- `variantcentrifuge/inheritance/vectorized_deducer.py` (802 lines)

**Modified:**
- `variantcentrifuge/inheritance/analyzer.py` (removed df.apply, added vectorized_deduce_patterns import and call)

## Commits

1. **feat(09-02): create vectorized Pass 1 pattern deduction with NumPy boolean masks**
   - Implement vectorized_deduce_patterns() using NumPy arrays
   - Encode all genotypes once using int8 matrix
   - Vectorized pattern detection: de novo, AD, AR, X-linked, mitochondrial
   - Handle edge cases: single sample, missing genotypes, incomplete pedigree
   - Commit: a04d261

2. **feat(09-02): wire vectorized Pass 1 into analyzer.py and validate**
   - Replace df.apply(deduce_patterns_for_variant) with vectorized_deduce_patterns()
   - Remove unused import
   - Fix per-sample fallback logic
   - Fix categorical dtype deprecation warning
   - All golden file scenarios pass
   - Commit: 4584b07

## Conclusion

Successfully vectorized Pass 1 (inheritance pattern deduction) using NumPy boolean masks, achieving:

1. **Clinical correctness:** All 10 golden file scenarios pass (clinically equivalent output)
2. **Code quality:** Clean linting, type checking, readable variable names, genetic logic comments
3. **Test coverage:** All 140 inheritance tests + 599 unit tests pass with no regressions
4. **Performance potential:** Expected 10-100x speedup (to be confirmed in Plan 05 benchmarks)

**Ready to proceed with Phase 9 Plan 03 (Compound Het Vectorization).**
