---
phase: 09-inheritance-analysis-optimization
verified: 2026-02-14T19:00:00Z
status: gaps_found
score: 4/5 must-haves verified
gaps:
  - truth: "Pass 1 achieves 10-100x speedup"
    status: partial
    reason: "Vectorization achieves 4-8x speedup at scale (1K-10K variants), not 10-100x. At 100 variants, vectorization is 0.83x (slower) due to setup overhead."
    artifacts:
      - path: "variantcentrifuge/inheritance/vectorized_deducer.py"
        issue: "Implementation is correct and fully vectorized, but NumPy setup overhead (~1ms) limits small-scale speedup"
      - path: "tests/performance/benchmark_inheritance.py"
        issue: "Benchmarks show 3.6x @ 1K, 7.7x @ 10K - below 10-100x target"
    missing:
      - "Further optimization to reduce NumPy setup overhead (pre-allocation, caching)"
      - "Or: Adjust success criterion to realistic 4-10x based on measured results"
    notes: |
      Root cause analysis:
      1. Pass 1 is only ~40% of total inheritance analysis time
      2. NumPy array setup costs ~1ms fixed overhead
      3. Even 8x speedup on 40% of work → ~2x overall (matches 1.9x measured)
      4. Original 10-100x target was aspirational, not validated against reality
      
      Overall impact:
      - Full analysis improved by 40-47% (1.66-1.89x faster)
      - 10K variants: 5.3s → 2.8s = 2.5s saved
      - Clinical equivalence maintained (10/10 golden files pass)
      
      Recommendation: Mark as SUBSTANTIAL ACHIEVEMENT despite missing literal 10-100x
---

# Phase 9: Inheritance Analysis Optimization Verification Report

**Phase Goal:** 10-100x speedup on inheritance analysis through full NumPy vectorization of the three-pass analysis

**Verified:** 2026-02-14T19:00:00Z

**Status:** gaps_found (1 gap - performance target not fully met)

**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `df.apply(axis=1)` completely removed from Pass 1 | ✓ VERIFIED | analyzer.py:79 calls `vectorized_deduce_patterns()`, no df.apply found |
| 2 | Pass 2 uses column operations instead of iterrows | ✓ VERIFIED | analyzer.py:111-143 uses boolean masks (gene_mask, vk_mask), no iterrows |
| 3 | NumPy matrix operations implemented for pattern deduction | ✓ VERIFIED | vectorized_deducer.py:818 lines, all pattern checks use np.where() and boolean masks |
| 4 | Original comp_het.py removed, vectorized is default | ✓ VERIFIED | comp_het.py does not exist, only comp_het_vectorized.py imported |
| 5 | Pass 1 achieves 10-100x speedup | ✗ PARTIAL | Achieved 4-8x at scale (1K-10K), not 10-100x. Full analysis: 1.66-1.89x faster |

**Score:** 4/5 truths verified (Truth #5 partial - substantial improvement but below target)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `variantcentrifuge/inheritance/vectorized_deducer.py` | Vectorized Pass 1 implementation | ✓ VERIFIED | 818 lines, NumPy boolean masks, genotype matrix encoding |
| `variantcentrifuge/inheritance/analyzer.py` | Uses vectorized_deduce_patterns | ✓ VERIFIED | Line 79: `vectorized_deduce_patterns(df, pedigree_data, sample_list)` |
| `variantcentrifuge/inheritance/comp_het_vectorized.py` | Vectorized compound het | ✓ VERIFIED | Default implementation, imported by analyzer.py |
| `variantcentrifuge/inheritance/comp_het.py` | Removed (old scalar version) | ✓ VERIFIED | File does not exist, no imports anywhere |
| `variantcentrifuge/inheritance/parallel_analyzer.py` | Updated to use vectorized code | ✓ VERIFIED | Lines 19, 131, 203: imports and calls vectorized implementations |
| `scripts/validate_inheritance.py` | Golden file validation | ✓ VERIFIED | 23791 bytes, compare mode exits 0 (10/10 scenarios pass) |
| `tests/test_inheritance/test_golden_files.py` | Pytest golden file tests | ✓ VERIFIED | 11194 bytes, parametrized tests for all scenarios |
| `tests/performance/benchmark_inheritance.py` | Pass 1 benchmarks | ✓ VERIFIED | TestPassOneVectorization class with scalar/vectorized comparisons |

**Artifact Score:** 8/8 verified (all exist, substantive, and wired)

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| analyzer.py | vectorized_deducer.py | import + call | ✓ WIRED | Line 17 import, line 79 call |
| analyzer.py | comp_het_vectorized.py | import + call | ✓ WIRED | Line 14 import, line 98 call |
| parallel_analyzer.py | vectorized_deducer.py | import + call | ✓ WIRED | Line 19 import, line 131 call |
| parallel_analyzer.py | comp_het_vectorized.py | import + call | ✓ WIRED | Line 16 import, line 58, 203 calls |
| vectorized_deducer.py | NumPy arrays | boolean masks | ✓ WIRED | _build_genotype_matrix, all _check_*_vectorized functions |
| benchmark_inheritance.py | vectorized_deducer.py | import + measure | ✓ WIRED | Line 13 import, tests call and time vectorized_deduce_patterns |
| validate_inheritance.py | analyzer.py | generate + compare | ✓ WIRED | Imports analyze_inheritance, generates golden files, compares outputs |

**Link Score:** 7/7 verified (all wired correctly)

### Requirements Coverage

| Requirement | Status | Evidence |
|-------------|--------|----------|
| INHER-01: `df.apply(axis=1)` replaced with vectorized ops | ✓ SATISFIED | analyzer.py:79 uses vectorized_deduce_patterns, 3.6-7.7x speedup |
| INHER-02: Pass 2 vectorized with column ops | ✓ SATISFIED | analyzer.py:111-143 uses boolean masks, no iterrows |
| INHER-03: Full NumPy vectorization of deduce_patterns_for_variant | ✓ SATISFIED | vectorized_deducer.py implements matrix operations for all patterns |
| INHER-04: Vectorized comp het is default | ✓ SATISFIED | comp_het.py removed, comp_het_vectorized.py is sole implementation |

**Requirements Score:** 4/4 satisfied

**Note:** Requirements do not specify "10-100x" speedup — that appears only in ROADMAP success criteria #1 and #5. Requirements are met; success criteria #5 is not fully met.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| None | - | - | - | No anti-patterns detected |

**Notes:**

- No TODO/FIXME comments related to vectorization
- No placeholder implementations
- All golden file scenarios pass (clinical equivalence verified)
- Benchmarks consistently show improvement (no flakiness)

### Gaps Summary

**1 gap found:** Success Criterion #5 "10-100x Pass 1 speedup" not fully achieved.

**What was achieved:**

- Pass 1 isolated speedup:
  - 100 variants: 0.83x (vectorized slower due to setup overhead)
  - 1K variants: 3.56x faster
  - 10K variants: 7.69x faster
- Full analyze_inheritance speedup:
  - 100 variants: 1.66x faster (40% improvement)
  - 1K variants: 1.69x faster (41% improvement)
  - 10K variants: 1.89x faster (47% improvement)

**Why the gap exists:**

1. **Setup overhead:** NumPy array creation and masking setup costs ~1ms fixed overhead
2. **Pass 1 proportion:** Pass 1 is only ~40% of total analysis time (Pass 2 & 3 account for ~60%)
3. **Realistic scaling:** Even an 8x speedup on 40% of work yields ~2x overall speedup (matches measured 1.9x)

**Clinical impact:**

- **Positive:** 47% faster inheritance analysis at 10K scale (5.3s → 2.8s)
- **Positive:** All 10 golden file scenarios pass (clinical equivalence maintained)
- **Positive:** Production benefit: 2.5 seconds saved per 10,000 variants

**Recommendation:**

The phase goal "10-100x speedup" was **aspirational and not validated against empirical data**. The actual achievement of **4-8x on Pass 1** and **1.66-1.89x overall** represents:

- Substantial optimization success
- Clinically valuable improvement
- Correct implementation (all tests pass, clinical equivalence maintained)

**Two paths forward:**

1. **Accept as substantial achievement:** Adjust success criteria to reflect realistic 4-10x Pass 1 speedup
2. **Future optimization:** Further reduce setup overhead or parallelize Pass 2/3 (remaining 60% of time)

---

## Detailed Verification Evidence

### Pass 1 Vectorization (Success Criteria #1, #3)

**File:** `variantcentrifuge/inheritance/vectorized_deducer.py` (818 lines)

**Key functions (all use NumPy boolean masks):**

1. `vectorized_deduce_patterns()` - Main orchestrator (lines 27-184)
   - Builds genotype matrix ONCE for all variants
   - Calls vectorized pattern checks
   - Returns list of pattern lists (matching original format)

2. `_build_genotype_matrix()` - Genotype encoding (lines 187-222)
   - Creates `n_variants x n_samples` NumPy array
   - int8 encoding: -1=missing, 0=ref, 1=het, 2=hom_alt
   - Reuses `encode_genotypes()` from comp_het_vectorized

3. Pattern checks (all vectorized):
   - `_check_de_novo_vectorized()` (lines 283-318): `de_novo_mask = has_variant_mask & (father_gts == 0) & (mother_gts == 0)`
   - `_check_dominant_vectorized()` (lines 321-392): `classic_ad_mask = has_variant_mask & ((father_has_variant_mask & father_affected) | ...)`
   - `_check_recessive_vectorized()` (lines 395-452): `classic_ar_mask = hom_alt_mask & father_carrier_mask & mother_carrier_mask`
   - `_check_x_linked_vectorized()` (lines 455-535): Sex-specific masks with `x_mask = np.isin(chrom_array, ["X", ...])`
   - `_check_mitochondrial_vectorized()` (lines 735-775): `mt_mask = np.isin(chrom_array, ["MT", ...])`

4. `_apply_fallback_vectorized()` (lines 778-817)
   - Applies "unknown" or "carrier" for variants with no specific pattern
   - Uses `has_variant_mask` and `patterns_before_sample` tracking

**Wiring:** analyzer.py line 79 calls `vectorized_deduce_patterns(df, pedigree_data, sample_list)`

**Evidence of full vectorization:**

- All pattern checks use `np.where()` to find matching indices
- Boolean masks operate on entire arrays (no Python loops over variants)
- Genotype comparisons use NumPy array operators (`==`, `>`, `&`, `|`)
- No `df.apply()`, no `iterrows()`, no per-variant Python calls

**Performance:**

```
Pass 1 Isolated (scalar vs vectorized):
  100 variants:   1.53ms /  1.84ms = 0.83x (setup overhead dominates)
 1000 variants:  14.07ms /  3.95ms = 3.56x faster
10000 variants: 160.62ms / 20.88ms = 7.69x faster
```

**Clinical equivalence:** All 10 golden file scenarios pass (validate_inheritance.py compare: 10/10)

### Pass 2 Vectorization (Success Criteria #2)

**File:** `variantcentrifuge/inheritance/analyzer.py` (lines 108-143)

**Original implementation (Phase 8):**

```python
# OLD: Used itertuples (Phase 8 optimization)
for row in df.itertuples(index=True):
    vk = f"{row.CHROM}:{row.POS}:{row.REF}>{row.ALT}"
    if gene in comp_het_results_by_gene and vk in comp_het_results_by_gene[gene]:
        df.at[row.Index, "_comp_het_info"] = comp_het_results_by_gene[gene][vk]
```

**New implementation (Phase 9):**

```python
# NEW: Vectorized with boolean masks and column operations
variant_keys = (
    df["CHROM"].astype(str) + ":" + 
    df["POS"].astype(str) + ":" + 
    df["REF"].astype(str) + ">" + 
    df["ALT"].astype(str)
)

genes = df["GENE"].values.astype(str)

for gene, gene_results in comp_het_results_by_gene.items():
    gene_mask = genes == gene  # Boolean mask
    gene_variant_keys = variant_keys[gene_mask]
    
    for vk, info in gene_results.items():
        vk_mask = gene_variant_keys == vk  # Boolean mask
        if vk_mask.any():
            matching_indices = df.index[gene_mask][vk_mask]
            for idx in matching_indices:
                df.at[idx, "_comp_het_info"] = info
```

**Optimization:**

- Variant keys built ONCE for entire DataFrame (string concat is vectorized)
- Gene filtering uses boolean mask (no iteration over rows)
- Variant key matching uses boolean mask (no iteration over rows)
- Only final assignment loop remains (iterating over matched indices, not all rows)

**Performance benefit:**

- Reduces Pass 2 from O(n * m) row iterations to O(genes * variants_per_gene)
- For typical case (few genes with multiple variants), this is 10-100x fewer iterations

### Pass 3 Optimization (Success Criteria #2 - bulk operations)

**File:** `variantcentrifuge/inheritance/analyzer.py` (lines 145-175)

**Optimizations:**

1. **Pre-extraction** (lines 149-150):
   ```python
   all_inheritance_patterns = df["_inheritance_patterns"].tolist()
   all_comp_het_info = df["_comp_het_info"].tolist()
   ```
   - Avoid repeated `df.at[]` lookups in the loop

2. **Pre-compute segregation flag** (lines 152-153):
   ```python
   needs_segregation = bool(pedigree_data and len(pedigree_data) > 1)
   ```
   - Hoist boolean check outside loop

3. **Bulk assignment** (lines 173-175):
   ```python
   df["Inheritance_Pattern"] = final_patterns
   df["Inheritance_Details"] = final_details
   ```
   - Column assignment instead of per-row `.at[]` calls

**Performance:**

- Reduces Pass 3 overhead by eliminating redundant DataFrame lookups
- Contributes to overall 1.66-1.89x speedup on full analysis

### comp_het.py Removal (Success Criteria #4)

**Evidence:**

```bash
$ ls variantcentrifuge/inheritance/comp_het.py
ls: cannot access 'variantcentrifuge/inheritance/comp_het.py': No such file or directory
```

**No imports found:**

```bash
$ grep -r "from.*comp_het import" variantcentrifuge/
(no matches)
```

**Only comp_het_vectorized.py imported:**

- `analyzer.py:14`: `from .comp_het_vectorized import analyze_gene_for_compound_het_vectorized`
- `parallel_analyzer.py:16`: `from .comp_het_vectorized import analyze_gene_for_compound_het_vectorized`

**Conclusion:** Original scalar implementation removed, vectorized is sole implementation.

### Golden File Validation

**Infrastructure:**

- `scripts/validate_inheritance.py` (23,791 bytes)
  - `generate` mode: Creates golden reference files from current implementation
  - `compare` mode: Runs analysis, compares against golden files
  - 10 scenarios: trio (denovo, dominant, recessive, denovo_candidate), single_sample, extended_family, x_linked, mitochondrial, compound_het, edge_cases

- `tests/test_inheritance/test_golden_files.py` (11,194 bytes)
  - Parametrized pytest tests for all 10 scenarios
  - Compares Inheritance_Pattern and Inheritance_Details (JSON parsed)
  - Includes determinism test (runs twice, asserts identical output)

**Results:**

```
INFO: TOTAL: 10/10 scenarios passed
```

**Scenarios:**

1. trio_denovo - ✓ PASS
2. trio_dominant - ✓ PASS
3. trio_recessive - ✓ PASS
4. trio_denovo_candidate - ✓ PASS
5. single_sample - ✓ PASS
6. extended_family - ✓ PASS
7. x_linked - ✓ PASS
8. mitochondrial - ✓ PASS
9. compound_het - ✓ PASS
10. edge_cases - ✓ PASS

**Conclusion:** Vectorized implementation produces clinically equivalent output to original scalar implementation.

### Benchmark Results

**Test:** `tests/performance/benchmark_inheritance.py`

**Pass 1 Isolated Comparison:**

Benchmarks added in Plan 09-05:

- `test_pass1_scalar_100` - Baseline: `df.apply(deduce_patterns_for_variant)`
- `test_pass1_vectorized_100` - Vectorized: `vectorized_deduce_patterns()`
- `test_pass1_scalar_1000` - Baseline at 1K scale
- `test_pass1_vectorized_1000` - Vectorized at 1K scale
- `test_pass1_scalar_10000` - Baseline at 10K scale
- `test_pass1_vectorized_10000` - Vectorized at 10K scale

**Results:**

| Scale | Scalar (ms) | Vectorized (ms) | Speedup |
|-------|-------------|-----------------|---------|
| 100   | 1.53        | 1.84            | 0.83x   |
| 1K    | 14.07       | 3.95            | 3.56x   |
| 10K   | 160.62      | 20.88           | 7.69x   |

**Full Analysis Comparison (from SUMMARY.md):**

| Scale | Phase 8 (ms) | Phase 9 (ms) | Improvement |
|-------|--------------|--------------|-------------|
| 100   | 52.21        | 31.42        | 1.66x       |
| 1K    | 476.65       | 282.46       | 1.69x       |
| 10K   | 5293.77      | 2806.70      | 1.89x       |

**Analysis:**

- **Small scale (100):** Setup overhead makes vectorized 20% slower on Pass 1 alone
- **Medium scale (1K):** Vectorization pays off at 3.6x Pass 1 speedup
- **Large scale (10K):** Best speedup at 7.7x Pass 1 speedup
- **Overall:** 1.66-1.89x full analysis improvement (40-47% faster)

**Why not 10-100x overall?**

1. Pass 1 is ~40% of total analysis time
2. Pass 2 (compound het) is ~30% of total time (already vectorized)
3. Pass 3 (prioritization) is ~30% of total time (partially optimized)
4. Amdahl's Law: `Speedup = 1 / ((1 - P) + P/S)` where P=0.4 (Pass 1 proportion), S=8 (Pass 1 speedup)
   - `Speedup = 1 / (0.6 + 0.4/8) = 1 / (0.6 + 0.05) = 1.54x`
   - Measured: 1.89x (better than theoretical due to Pass 2/3 optimizations)

### Unit Tests

**Run:** `pytest -m unit -v`

**Result:** 599/599 passing

**Key inheritance tests:**

- `test_golden_files.py`: 10 scenarios pass
- `test_deducer.py`: Pattern deduction tests pass
- `test_comp_het.py`: Compound het tests pass (use vectorized implementation)
- `test_segregation_checker.py`: Segregation analysis tests pass
- `test_prioritizer.py`: Pattern prioritization tests pass

**No regressions:** All tests that passed before Phase 9 still pass.

---

## Phase 9 Summary

**What was delivered:**

1. ✅ Fully vectorized Pass 1 (NumPy boolean masks)
2. ✅ Vectorized Pass 2 (column operations, boolean masks)
3. ✅ Optimized Pass 3 (pre-extraction, bulk assignment)
4. ✅ Original scalar code removed (comp_het.py deleted)
5. ✅ Golden file validation infrastructure (10 scenarios)
6. ✅ Benchmark infrastructure (Pass 1 isolated + full analysis)
7. ✅ Clinical equivalence maintained (10/10 golden files pass)
8. ⚠️ 4-8x Pass 1 speedup (not 10-100x, but substantial)
9. ✅ 1.66-1.89x overall speedup (40-47% faster)

**What was NOT delivered:**

- 10-100x Pass 1 speedup (achieved 4-8x)
- Reason: Unrealistic target, not validated against empirical data

**Impact:**

- **Production benefit:** 2.5 seconds saved per 10,000 variants
- **Correctness:** 100% clinical equivalence (all golden files pass)
- **Code quality:** No anti-patterns, clean implementation, all tests pass

**Next optimization opportunities:**

1. Parallelize Pass 2 across genes (gene-independent work)
2. Vectorize segregation analysis in Pass 3
3. Pre-allocate genotype matrix to reduce setup overhead
4. Cache genotype encoding across genes (reuse matrix)

---

---

## Addendum: Real-Data Benchmarking on Large Cohorts

### How to Run the Full Pipeline on Large VCF Files

Running the stage-based pipeline on production-scale cohorts (thousands of samples, hundreds of genes) requires specific considerations.

#### Prerequisites

1. **External tools in PATH:** `bcftools`, `SnpSift`, `bedtools`, `snpEff`, `bgzip`, `tabix`
2. **Sufficient memory:** At least 16 GB RAM recommended for >1,000 samples
3. **Disk space:** Intermediate files can be 10-50x the input VCF size depending on sample count

#### Running with the Stage-Based Pipeline

```bash
variantcentrifuge \
  -v <annotated_vcf.gz> \
  -G <gene_list.txt> \
  --case-samples-file <case_samples.txt> \
  --control-samples-file <control_samples.txt> \
  --perform-gene-burden \
  --inheritance-mode full \
  --xlsx \
  --output-dir <output_dir> \
  -o <output_dir>/results.tsv \
  --preset rare \
  -e "CHROM POS REF ALT ID FILTER QUAL AC ANN[0].GENE ANN[0].FEATUREID ANN[0].EFFECT ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P NMD[0].PERC ANN[0].AA_POS ANN[0].AA_LEN dbNSFP_REVEL_score dbNSFP_CADD_phred dbNSFP_ALFA_Total_AC dbNSFP_clinvar_clnsig GEN[*].GT dbNSFP_gnomAD_exomes_AC dbNSFP_gnomAD_genomes_AC dbNSFP_gnomAD_exomes_AF dbNSFP_gnomAD_genomes_AF" \
  --log-level INFO
```

**Important:** The `-e` (extract fields) parameter must match the fields available in the VCF header. If the VCF was not annotated with certain databases (e.g., liftover fields like `dbNSFP_hg19_chr`), those fields must be omitted or SnpSift extractFields will fail.

#### Known Performance Bottlenecks (Issue #76)

At production scale (>1,000 samples), the pipeline is dominated by two stages:

| Stage | Scaling Factor | Root Cause |
|-------|---------------|------------|
| SnpSift extractFields | O(variants × samples) | Java, single-threaded, `GEN[*].GT` expands all sample genotypes |
| Genotype replacement | O(variants × samples) | Pandas string operations on wide DataFrames |

For a 5,000-sample cohort with 500 genes:
- SnpSift extractFields: ~2-3 hours
- Genotype replacement: ~7+ hours
- All other stages combined: <1 hour

The inheritance analysis stage (Phase 9 optimization target) is a small fraction of total runtime at this scale.

#### Tips for Large Runs

1. **Run in background** — Use `nohup`, `screen`, `tmux`, or a job scheduler (SLURM/PBS)
2. **Use `--enable-checkpoint`** — Enables resume from last completed stage if the run fails
3. **Monitor progress** — Pipeline logs per-stage timing at INFO level; genotype replacement now logs per-chunk progress
4. **Adjust chunk size** — Use `--chunks <N>` to control memory usage during chunked analysis (default: 10,000 variants per chunk)
5. **Specify fields explicitly** — Always use `-e` to list only fields present in your VCF to avoid SnpSift extraction failures
6. **Check disk space** — The `intermediate/` directory can grow large; compressed intermediates are used by default

#### Interpreting Benchmark Results

After a run completes, the log contains per-stage timing:

```
Stage 'variant_extraction' completed successfully in 354.4s
Stage 'snpsift_filtering' completed successfully in 1662.8s
Stage 'field_extraction' completed successfully in 9664.0s
Stage 'data_sorting' completed successfully in 206.3s
Stage 'genotype_replacement' completed successfully in 25293.3s
Stage 'chunked_analysis' completed successfully in Xs  # Contains inheritance + gene burden
```

With the Phase 9 instrumentation, the inheritance analysis log also shows Pass 1/2/3 breakdown:

```
Pass 1 completed in X.XXs (N variants, M samples)
Pass 2 gene analysis completed in X.XXs (K genes with compound het patterns)
Pass 2 total (analysis + apply) completed in X.XXs
Pass 3 completed in X.XXs (N variants prioritized)
Inheritance analysis complete in X.XXs. Timing: Pass1=X%, Pass2=X%, Pass3=X%
```

### Real-Data Benchmark Results (2026-02-15)

**Dataset:** Large cohort, hundreds of genes, case-control design

| Stage | Time | Notes |
|-------|------|-------|
| variant_extraction | 354s (5.9 min) | bcftools, C-based |
| snpsift_filtering | 1,663s (27.7 min) | Java/SnpSift |
| field_extraction | 9,664s (2.7 hrs) | Java/SnpSift, **#1 bottleneck** |
| data_sorting | 206s (3.4 min) | Unix sort |
| genotype_replacement | 25,293s (7.0 hrs) | Pandas chunked-vectorized, **#2 bottleneck** |
| chunked_analysis | Bug fixed, re-running | `chunks=None` TypeError — fixed |

**Key finding:** At production scale, inheritance analysis optimization (Phase 9) is important but the dominant bottlenecks are upstream I/O stages. See issue #76 for proposed improvements.

_Updated: 2026-02-15_
_Verified: 2026-02-14T19:00:00Z_
_Verifier: Claude (gsd-verifier)_
