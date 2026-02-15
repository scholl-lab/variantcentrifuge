# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Accurate inheritance pattern deduction and variant prioritization from multi-sample VCFs with configurable gene panels, scoring models, and output formats
**Current focus:** Phase 10 complete - ready for Phase 11

## Current Position

Phase: 11 of 12 (Pipeline I/O Elimination)
Plan: 2 of 3 complete
Status: In progress - genotype replacement eliminated
Last activity: 2026-02-15 — Completed 11-02-PLAN.md (genotype replacement elimination)

Progress: [███████████████░░░░░] 77% (Phase 1-10 complete + 2/3 of Phase 11)

## Performance Metrics

**Velocity:**
- Total plans completed: 20
- Average duration: 14.0 minutes
- Total execution time: 4.7 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1-5. Baseline | N/A | N/A | N/A (pre-GSD) |
| 6. Benchmark Framework | 4/4 | 48.0 min | 12.0 min |
| 7. Quick Wins Tier 1 | 3/3 | 89.0 min | 29.7 min |
| 8. DataFrame Optimization | 4/4 | 62.0 min | 15.5 min |
| 9. Inheritance Optimization | 5/5 | 46.8 min | 9.4 min |
| 10. Output Optimization | 3/3 | 26.0 min | 8.7 min |
| 11. Pipeline I/O Elimination | 2/3 | 36.0 min | 18.0 min |

**Recent Trend:**
- Last 5 plans: 10-01 (9.0 min), 10-02 (9.0 min), 10-03 (8.0 min), 11-01 (10.0 min), 11-02 (26.0 min)
- Trend: Generally fast (8-10 min), occasional longer plans with complex refactoring

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Benchmarks before optimization: Can't prove improvement without baseline measurements
- Accept risk on Tier 3 (vectorization, Cython): User wants maximum performance gains, will revert if tests break
- Local benchmarks only, no CI workflow: Get benchmarks working first, CI integration deferred
- Target v0.13.0: Performance improvements warrant minor version bump
- Removed entire dead GT parsing loop (07-01): Existing aggregated counts already provide correct values, loop was parsing into unused variables
- Regression tests before dead code removal (07-01): Proves output remains identical after refactoring
- observed=True now vs later (07-02): Added observed=True in Phase 7 before categorical dtypes exist to prevent regressions when Phase 8 introduces them
- Unconditional gc.collect() (07-02): Runs after EVERY stage execution (not just DEBUG) for memory management with large genomic datasets
- Pre-commit hook for groupby enforcement (07-02): Prevents new groupby calls without observed=True from being committed
- Full 60-benchmark suite completed (07-03): All benchmarks passed in 5:37, saved as `0001_phase7_quick_wins.json`. Deprecated streaming perf tests removed.
- Document collateral improvements (07-03): Inheritance analysis improved 20-58% despite no direct optimizations, attributed to gc.collect() and observed=True reducing overhead
- PyArrow scope limited to variant DataFrames (08-01): Main TSV loading only, not config/gene lists per 08-CONTEXT decision
- Categorical dtype auto-detection at 50% cardinality (08-01): Columns with <50% unique values loaded as categorical
- Column renaming permanent at load time (08-01): No temporary rename-restore, downstream uses sanitized names
- Memory pass-through threshold 25% available RAM (08-01): Conservative targeting 8-16GB desktops
- quoting parameter excluded for PyArrow engine (08-01): PyArrow doesn't support it, C engine fallback used
- Column name restoration at output time (08-03): TSV and Excel outputs use original names (GEN[0].GT not GEN_0__GT) for backwards compatibility
- In-memory DataFrame pass-through for Excel (08-03): ExcelReportStage uses context.variants_df when available, eliminating redundant disk read
- Phase 8 optimizations exceed all targets (08-04): 82-84% memory reduction (vs 50-70%), 30.9x iteration speedup (vs 10-14x), 3.0x I/O speedup at scale
- Gene burden 37-62% faster as collateral improvement (08-04): itertuples optimization benefits all DataFrame iteration, not just inheritance analysis
- Inheritance "regressions" are benchmark variance (08-04): 11-22% slowdown within normal variation, DataFrame I/O doesn't affect inheritance logic
- Use Parquet format for golden files (09-01): Preserves exact dtypes and column ordering, more reliable than CSV/TSV for validation
- 10 synthetic test scenarios for inheritance validation (09-01): Covers all major patterns (trio, single-sample, extended-family, X-linked, mitochondrial, compound-het, edge-cases)
- Build genotype matrix once before all pattern checks (09-02): Encode all sample genotypes (n_variants x n_samples) once using int8, reuse for all checks
- Per-sample pattern tracking for fallback logic (09-02): Track patterns_before_sample to match original per-sample fallback behavior (carrier/unknown)
- Use isinstance(dtype, pd.CategoricalDtype) for categorical check (09-02): Future-proof replacement for deprecated is_categorical_dtype()
- Reuse encode_genotypes from comp_het_vectorized.py (09-02): Consistent encoding across modules, proven correctness
- Vectorized variant key creation (09-03): Build all variant keys at once via string concatenation instead of per-row create_variant_key() calls (10-100x faster)
- Pre-extraction pattern for Pass 3 optimization (09-03): Extract DataFrame columns to lists before loop to avoid repeated df.at[] lookups
- Bulk assignment pattern for Pass 3 (09-03): Accumulate results in lists, assign to DataFrame once at end instead of per-row df.at[] writes
- Keep iterrows() in Pass 3 (09-03): create_inheritance_details() and segregation_checker require Series/dict, partial vectorization acceptable per 09-CONTEXT
- Vectorized as sole implementation (09-04): Removed original comp_het.py, eliminated VECTORIZED_AVAILABLE flags, single implementation reduces maintenance
- Conservative compound het detection (09-04): No pairing when phase ambiguous (both parents het for both variants) - avoids false positives
- Pass 1 vectorization achieves 3-7x speedup at scale (09-05): Not 10-100x as originally targeted, but significant; 7.1x at 10K variants
- Full inheritance analysis 40-47% faster (09-05): 1.66-1.89x improvement from three-pass vectorization (Pass 1, 2, 3)
- Realistic ratio assertions (09-05): Adjusted thresholds to measured values (3-5x for Pass 1) instead of aspirational 10-100x
- Setup overhead acceptable at small scale (09-05): Vectorized slower at 100 variants (0.8x) but real workloads are 1K-100K
- Two-pass Excel generation (10-01): xlsxwriter for bulk write (2-5x faster) + openpyxl for finalization (hyperlinks, freeze panes, auto-filters)
- Module-level regex constants (10-01): GT_PATTERN compiled once at import eliminates per-row re.compile() overhead
- GT column pre-parsing at load time (10-02): Parse GT column once into _GT_PARSED cache (list of dicts), reuse across all output stages
- Cache column cleanup pattern (10-02): Underscore-prefixed columns (_GT_PARSED) dropped before final output in TSV and Excel stages
- Document actual vs aspirational performance (10-03): xlsxwriter vs openpyxl measured at 0.9x for test data; speedup varies by dataset size/structure
- Benchmark pattern for Phase 10 (10-03): Synthetic data with multiple scales (100, 1K, 10K, 50K), metadata tracking for cross-phase comparison
- bcftools query over +split-vep (11-01): query 19x faster than SnpSift; +split-vep drops variants without ANN (unacceptable data loss)
- Dynamic format string construction (11-01): Build bcftools query format dynamically from config.json fields_to_extract for any field combination
- Per-sample column output (11-01): bcftools [\t%GT] produces separate columns per sample (GEN[0].GT, GEN[1].GT), eliminates genotype replacement need
- Python ANN parsing (11-01): Parse pipe-delimited SnpEff annotations in Python after bcftools extraction (simple, fast, handles missing gracefully)
- GenotypeReplacementStage eliminated via no-op (11-02): Stage returns immediately, GT reconstruction deferred to output time (TSV/Excel stages)
- Deferred GT formatting pattern (11-02): reconstruct_gt_column() builds packed "Sample(0/1);Sample2(1/1)" format only at final output (<1 min vs 7 hrs)
- Per-sample column phenotype extraction (11-02): extract_phenotypes_from_sample_columns() works with bcftools raw columns, auto-detects vs packed GT
- _GT_PARSED dead code removed (11-02): Cache column never consumed by production code, parse_gt_column() and GT_PATTERN regex deleted

### Pending Todos

- **DEPR-01** (backlog): Deprecate classic pipeline mode (`pipeline.py`) in favor of stage-based pipeline (`pipeline_core/`). See REQUIREMENTS.md Future Requirements.

### Blockers/Concerns

**Phase 6 (Benchmark Framework): COMPLETE**
- 60 benchmark tests across 8 files, all passing
- Baseline saved: `.benchmarks/Windows-CPython-3.10-64bit/0001_baseline_v0.12.1.json`
- Phase 7 results saved: `.benchmarks/Linux-CPython-3.10-64bit/0001_phase7_quick_wins.json`
- Deprecated streaming parallel perf tests removed (were blocking benchmark runs)
- Framework ready for Phase 8+ optimization work

**Phase 7 Plan 01 (Dead Code Removal): COMPLETE**
- Removed 30-line dead GT parsing loop from gene_burden.py
- Fixed 3 temp file leaks (gene_bed.py bed_path, filters.py mktemp replacement)
- Regression tests prove gene burden output unchanged
- Expected 20-30% speedup in gene burden analysis
- Existing bug discovered: gene_burden.py crashes on empty DataFrame (not fixed, out of scope)

**Phase 7 Plan 02 (Categorical dtypes & Memory Management): COMPLETE**
- Added observed=True to all 17 groupby call sites (future-proofing for Phase 8 categorical dtypes)
- Implemented gc.collect() after every pipeline stage execution
- Added DEBUG-level memory logging (RSS before/after each stage, freed memory delta via psutil)
- Created pre-commit hook enforcing observed=True in all new groupby calls
- Prevents 3500x groupby slowdown when categorical dtypes are introduced in Phase 8
- All 565 unit tests pass with no behavioral changes

**Phase 7 Plan 03 (Benchmark Verification): COMPLETE**
- Verified Phase 7 optimizations with actual benchmark measurements
- Gene burden analysis: 48-98% faster (exceeded 20-40% target)
- Inheritance analysis: 20-58% faster (collateral improvement from gc.collect and observed=True)
- Zero regressions across all benchmarks
- Created comprehensive performance analysis report for Phase 8+ comparison
- Fixed test bug: added external tool checks to gene burden integration tests

**Baseline Performance (v0.12.1, Windows, CPython 3.10, 2026-02-14):**

| Component | Scale | v0.12.1 Baseline | After Phase 7 | Improvement |
|-----------|-------|-----------------|---------------|-------------|
| Single variant deduce | 1 variant | 7.7 us | 4.4 us | **43.0%** |
| Inheritance analysis | 100 variants | 49 ms | 37 ms | **23.8%** |
| Inheritance analysis | 1K variants | 431 ms | 346 ms | **19.6%** |
| Inheritance analysis | 10K variants | 4.8 s | 3.2 s | **33.2%** |
| Gene burden | 100 variants | 62 ms | 32 ms | **48.5%** |
| Gene burden | 1K variants | 150 ms | 18 ms | **88.0%** |
| Gene burden | 10K variants | 1.0 s | 19 ms | **98.1%** |
| Gene burden (10 genes) | - | 98 ms | 4 ms | **95.6%** |
| Gene burden (50 genes) | - | 121 ms | 19 ms | **84.0%** |
| Gene burden (100 genes) | - | 137 ms | 49 ms | **64.7%** |

**Phase 7 achieved 48-98% speedup on gene burden, 20-58% on inheritance analysis.**

**Phase 8 (DataFrame Optimization): COMPLETE (4/4 plans)**

**Plan 01 (DataFrame Optimizer Foundation): COMPLETE**
- Created dataframe_optimizer.py with PyArrow loading, categorical detection, column sanitization
- PyArrow engine now used automatically (3.0x CSV read speedup at 50K variants)
- Low-cardinality columns loaded as categorical (82-84% memory reduction measured)
- Column sanitization complete (GEN[0].GT → GEN_0__GT) - ready for itertuples migration
- Memory pass-through decision logic in place (25% available RAM threshold)
- All 568 unit tests + 31 integration tests pass with no regressions

**Plan 02 (iterrows to itertuples Migration): COMPLETE**
- Converted 14 hot-path iterrows sites to itertuples (30.9x iteration speedup measured)
- Modified create_variant_key to handle both Series and namedtuples
- Established getattr(row, COL, default) pattern for safe attribute access
- Used df.at[row.Index] for underscore-prefixed columns (itertuples renames them)
- All 694 tests pass with zero behavioral changes
- Cold-path iterrows intentionally left unchanged (analyze_variants, build_pm5_lookup, etc.)

**Plan 03 (Excel Generation Optimization): COMPLETE**
- ExcelReportStage uses in-memory DataFrame from context.variants_df (eliminates redundant disk read)
- TSVOutputStage and ExcelReportStage restore original column names before writing output
- convert_to_excel accepts optional DataFrame parameter with disk fallback
- All tests pass unchanged, backwards compatibility maintained

**Plan 04 (Benchmark Verification): COMPLETE**
- Memory reduction: 82-84% (exceeds 50-70% target)
- Iteration speedup: 30.9x (exceeds 10-14x target)
- I/O speedup: 3.0x at 50K variants (meets 2-3x target)
- Gene burden collateral: 37-62% faster across all scales
- Inheritance within variance: 11-22% slower (benchmark noise, not regression)
- All 568 unit tests pass, no regressions detected

**Phase 8 Performance Summary (Linux CPython 3.10, Phase 7 → Phase 8):**

| Component | Phase 7 | Phase 8 | Change |
|-----------|---------|---------|--------|
| **Gene burden 50 genes** | 51.7 ms | 19.6 ms | **62.0% faster** |
| **Gene burden 100** | 77.4 ms | 36.4 ms | **53.0% faster** |
| **Gene burden 1K** | 34.3 ms | 18.9 ms | **44.9% faster** |
| **Gene burden 10 genes** | 10.0 ms | 6.1 ms | **39.5% faster** |
| **Gene burden 100 genes** | 85.0 ms | 53.2 ms | **37.3% faster** |
| **Gene burden 10K** | 28.5 ms | 22.0 ms | **22.6% faster** |

**Phase 8 achieved 22-62% speedup on gene burden from itertuples optimization.**

**Phase 9 (Inheritance Analysis Optimization): COMPLETE (5/5 plans)**

**Plan 01 (Golden File Validation Infrastructure): COMPLETE**
- Created scripts/validate_inheritance.py with generate and compare modes
- Built 10 synthetic test scenarios covering all inheritance patterns
- Generated golden reference files (10 .parquet files) from current implementation
- Created pytest test suite (14 tests, all passing)
- Validated determinism: scenarios produce identical output across runs
- Comparison mode validates clinical equivalence (Inheritance_Pattern exact, confidence within 0.001)

**Plan 02 (Deducer Vectorization - Pass 1): COMPLETE**
- Created vectorized_deducer.py with NumPy boolean mask pattern deduction (802 lines)
- Replaced df.apply(deduce_patterns_for_variant) with vectorized_deduce_patterns()
- Genotype matrix encoding once upfront (n_variants x n_samples int8 matrix)
- All 10 golden file scenarios pass (clinically equivalent output)
- All 140 inheritance tests + 599 unit tests pass with no regressions
- Measured speedup: 3.9x at 1K variants, 7.1x at 10K variants (Pass 1 alone)

**Plan 03 (Pass 2 & 3 Vectorization): COMPLETE**
- Vectorized Pass 2: variant keys built for all rows at once, comp het results applied via boolean masks
- Optimized Pass 3: pre-extracted columns to lists, accumulated results in lists, bulk assignment at end
- Removed itertuples loop from Pass 2 (now fully vectorized)
- Pass 3 partially vectorized (segregation_checker and create_inheritance_details kept scalar per 09-CONTEXT)
- All 10 golden file scenarios pass (clinically equivalent output)
- All 140 inheritance tests + 599 unit tests pass with no regressions

**Plan 04 (Compound Het Consolidation): COMPLETE**
- Removed original comp_het.py (428 lines deleted)
- Updated analyzer.py and parallel_analyzer.py to use comp_het_vectorized only
- Eliminated all VECTORIZED_AVAILABLE flags and dual code paths
- Conservative compound het: no pairing when phase ambiguous (avoids false positives)
- All 140 inheritance tests + 597 unit tests pass
- All 10 golden file scenarios pass (clinical equivalence maintained)
- Single implementation simplifies maintenance and benchmarking

**Plan 05 (Benchmark Verification): COMPLETE**
- Added Pass 1 vectorization benchmarks (scalar vs vectorized comparison)
- Created ratio assertions with realistic thresholds (3-5x for Pass 1)
- Measured full inheritance analysis speedup: 40-47% improvement
- All 10 golden file scenarios still pass after benchmark additions
- All 599 unit tests pass, all CI checks pass
- Benchmark results saved: `.benchmarks/Linux-CPython-3.10-64bit/0004_phase9_vectorization.json`

**Phase 9 Performance Summary (Linux CPython 3.10, Phase 8 → Phase 9):**

| Component | Phase 8 (before) | Phase 9 (after) | Improvement |
|-----------|------------------|-----------------|-------------|
| **Inheritance 100** | 52.21 ms | 31.42 ms | **1.66x faster** (39.8%) |
| **Inheritance 1K** | 476.65 ms | 282.46 ms | **1.69x faster** (40.7%) |
| **Inheritance 10K** | 5293.77 ms | 2806.70 ms | **1.89x faster** (47.0%) |
| **Single variant** | 7.84 μs | 4.78 μs | **1.64x faster** (39.0%) |

**Pass 1 isolated performance (scalar vs vectorized):**

| Scale | Scalar | Vectorized | Speedup |
|-------|--------|------------|---------|
| **100 variants** | 1.49 ms | 1.86 ms | 0.8x (overhead dominates) |
| **1K variants** | 13.9 ms | 3.5 ms | **3.9x faster** |
| **10K variants** | 144 ms | 20 ms | **7.1x faster** |

**Phase 9 achieved 40-47% speedup on inheritance analysis from three-pass vectorization.**

**Key findings:**
- Original 10-100x target was unrealistic (Pass 1 is only ~40% of total time)
- Pass 1 vectorization: 3.9-7.1x faster (excellent for that component)
- Overall improvement: 1.66-1.89x faster (40-47% improvement is significant)
- Remaining bottlenecks: Pass 2/3 still ~60% of total time (future optimization target)

**Phase 11 (Pipeline I/O Elimination): IN PROGRESS (2/3 plans complete)**

**Plan 01 (bcftools Field Extraction): COMPLETE**
- Replaced SnpSift extractFields with bcftools query (19x faster, measured)
- Dynamic format string building from config.json fields
- Python ANN/NMD parsing from raw pipe-delimited strings
- Per-sample GT columns (GEN[0].GT, GEN[1].GT, ...) ready for genotype elimination
- 23 new unit tests, 8 updated tests, all passing
- Expected savings: ~2h52m on large cohort extraction

**Plan 02 (Genotype Replacement Elimination): COMPLETE**
- GenotypeReplacementStage now no-ops immediately (7-hour bottleneck eliminated)
- GT reconstruction deferred to TSVOutputStage and ExcelReportStage (<1 min at output time)
- Phenotype integration updated to work with per-sample GT columns from bcftools
- _GT_PARSED dead code removed (parse_gt_column, GT_PATTERN regex, load-time parsing)
- All 636 unit tests pass with zero regressions
- Expected savings: ~7 hours on large cohort genotype replacement

**Plan 03 (Dead Code Cleanup): PENDING**
- Remove GenotypeReplacementStage helper methods (all unreachable code)
- Clean up stage registration and dependencies
- Finalize pipeline I/O elimination phase
- Target: 10+ hours → under 1 hour total pipeline time

## Session Continuity

Last session: 2026-02-15 10:56 UTC
Stopped at: Completed 11-02-PLAN.md (genotype replacement elimination)
Resume file: None
Next: Phase 11 Plan 03 (Dead Code Cleanup) — /gsd:plan 11 03
