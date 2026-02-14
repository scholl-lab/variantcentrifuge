---
phase: 06-benchmark-framework
plan: 01
subsystem: testing
tags: [pytest-benchmark, synthetic-data, memory-tracking, fixtures]
requires: []
provides:
  - pytest-benchmark test framework
  - Synthetic variant/pedigree data generators
  - Memory tracking utilities with budget enforcement
  - Performance test fixtures
affects:
  - 06-02 (Inheritance analysis benchmarks)
  - 06-03 (Compound het benchmarks)
  - 06-04 (Gene burden benchmarks)
  - 06-05 (Genotype replacement benchmarks)
  - 06-06 (Scoring benchmarks)
tech-stack:
  added: [pytest-benchmark==5.2.3]
  patterns: [factory-fixtures, context-managers, synthetic-data-generation]
key-files:
  created:
    - benchmarks/.gitignore
    - tests/performance/__init__.py
    - tests/performance/conftest.py
    - tests/performance/helpers/__init__.py
    - tests/performance/helpers/synthetic_data.py
    - tests/performance/helpers/memory_budgets.py
  modified:
    - pyproject.toml
    - .gitignore
decisions: []
metrics:
  duration: 327 seconds (~5.5 minutes)
  completed: 2026-02-14
---

# Phase 06 Plan 01: Benchmark Framework Foundation Summary

**One-liner:** Created pytest-benchmark framework with reproducible synthetic data generators (variants, pedigrees, scoring configs) and tracemalloc-based memory tracking with warning-only budget enforcement

## Objective

Create the benchmark framework foundation: add pytest-benchmark dependency, build synthetic data generators for genomic DataFrames and pedigrees, create tracemalloc memory budget helpers, and set up the project structure for benchmark result storage.

## What Was Built

### 1. Benchmark Infrastructure

**pytest-benchmark Integration:**
- Added `pytest-benchmark>=5.1.0` to dev dependencies (installed v5.2.3)
- Created `benchmarks/` directory for ephemeral JSON results with `.gitignore` containing `*` and `!.gitignore`
- Added `benchmarks/*.json` to root `.gitignore` as belt-and-suspenders protection
- Made `tests/performance/` a proper Python package with `__init__.py`
- Created `tests/performance/helpers/` subpackage for utilities

**Why these choices:**
- Core pytest-benchmark only (no [histogram] extra) avoids heavy matplotlib dependency
- Gitignored benchmarks directory keeps results local while preserving directory structure
- Package structure enables clean imports and pytest fixture discovery

### 2. Synthetic Data Generators (`tests/performance/helpers/synthetic_data.py`)

**Four generator functions, all with reproducible seeding:**

**`generate_synthetic_variants(n_variants, n_samples, seed=42) -> pd.DataFrame`**
- Produces DataFrames matching variantcentrifuge column structure: `CHROM, POS, REF, ALT, GT, GENE, FILTER, EFFECT, IMPACT`
- Realistic genomic distributions:
  - Chromosomes: Weighted toward chr1-10, natural sort order (1, 2, ..., 22, X, Y)
  - Positions: Random integers 1-250M
  - REF/ALT: Single nucleotide bases, ensures REF ≠ ALT
  - Genotypes: Realistic allele frequency distribution (90% rare variants with 95/4/1% for 0-0/0-1/1-1, 8% uncommon with 80/15/5, 2% common with 50/40/10)
  - Genes: ~n_variants/10 unique genes (avg 10 variants per gene)
  - Filter: 95% PASS, 5% LowQual
  - Effect: Weighted distribution (40% missense, 30% synonymous, 10% frameshift, 5% stop_gained, 15% splice_region)
  - Impact: Deterministically mapped from Effect (HIGH for frameshift/stop_gained, MODERATE for missense, LOW for synonymous, MODIFIER for splice_region)
- GT column format: Comma-separated genotype strings per variant `"0/0,0/1,1/1,0/0,..."` (one per sample)
- All sample names: `SAMPLE_NNNN`, gene names: `GENE_NNNN`

**`generate_synthetic_pedigree(n_samples=3, seed=42) -> dict`**
- Minimum trio: CHILD_001 (affected), FATHER_001 (unaffected), MOTHER_001 (unaffected)
- Additional samples as unrelated individuals (father_id=0, mother_id=0)
- Returns dict matching variantcentrifuge pedigree format: `{sample_id: {sample_id, father_id, mother_id, sex, affected_status}}`

**`generate_synthetic_scoring_config() -> dict`**
- Minimal scoring config with variable mappings (IMPACT → impact_var, EFFECT → effect_var)
- Simple formulas using pd.eval syntax
- Enables scoring benchmarks without reading real config files from disk

**`generate_synthetic_gene_burden_data(n_variants, n_samples, n_genes=50, seed=42) -> tuple[pd.DataFrame, set, set]`**
- Generates DataFrame with GT column in gene_burden format: `"SAMPLE_0001(0/1);SAMPLE_0002(0/0);..."`
- 50/50 split into case_samples and control_samples sets
- Cases have higher variant frequency (70/25/5 for 0-0/0-1/1-1) vs controls (90/8/2)
- Returns (df, case_samples, control_samples)

**Reproducibility:**
- All generators use `numpy.random.default_rng(seed)` for deterministic output
- Same seed → identical DataFrames (verified with `df1.equals(df2)`)

**Privacy compliance:**
- Zero private identifiers - all synthetic (SAMPLE_NNNN, GENE_NNNN, CHILD/FATHER/MOTHER)
- No cohort names, patient IDs, or traceable references

### 3. Memory Tracking Utilities (`tests/performance/helpers/memory_budgets.py`)

**`MemoryTracker` context manager:**
- Wraps `tracemalloc` for peak memory measurement
- Properties: `peak_mb`, `current_mb` (converted from bytes to MB)
- Usage: `with MemoryTracker() as tracker: ...` then check `tracker.peak_mb`

**`warn_if_over_budget(peak_mb, budget_mb, context)` function:**
- Issues `warnings.warn()` if peak exceeds budget
- **Never raises exceptions** - memory violations are warning-only per CONTEXT.md
- Allows benchmarks to complete while flagging concerning usage

**Default budget constants (MB):**
- `INHERITANCE_BUDGET_MB = 512`
- `COMP_HET_BUDGET_MB = 256`
- `GENOTYPE_REPLACEMENT_BUDGET_MB = 1024`
- `GENE_BURDEN_BUDGET_MB = 512`
- `SCORING_BUDGET_MB = 256`

These are initial estimates at 2x expected peak for 10K variants × 100 samples. They will be refined once actual profiling data exists in subsequent benchmark plans.

### 4. Pytest Fixtures (`tests/performance/conftest.py`)

**Five fixtures for all benchmark tests:**

1. **`synthetic_variants`** - Factory fixture returning callable `_generate(n_variants, n_samples, seed=42)`
2. **`synthetic_pedigree`** - Factory fixture returning callable `_generate(n_samples=3, seed=42)`
3. **`synthetic_scoring_config`** - Direct fixture (not factory, always same config)
4. **`synthetic_gene_burden_data`** - Factory fixture returning callable `_generate(n_variants, n_samples, n_genes=50, seed=42)`
5. **`memory_tracker`** - Fresh `MemoryTracker` instance

**Why factory fixtures:**
- Allows benchmark tests to parameterize data size (100/1K/10K/50K variants)
- Enables multiple data generations per test with different seeds
- Clean separation: fixture provides factory, test controls parameters

**Fixture scope:**
- Local to `tests/performance/` (does not pollute root `tests/conftest.py`)
- Root conftest already auto-marks performance tests as `@pytest.mark.slow` and `@pytest.mark.performance`

## Deviations from Plan

None - plan executed exactly as written.

## Verification Results

All verification criteria passed:

✅ `python -c "import pytest_benchmark"` succeeds (v5.2.3 installed)
✅ `ls benchmarks/.gitignore` exists
✅ `ls tests/performance/__init__.py` exists
✅ `ls tests/performance/helpers/__init__.py` exists
✅ `grep pytest-benchmark pyproject.toml` returns match
✅ Synthetic variants shape (100, 9) with expected columns
✅ Reproducibility: `df1.equals(df2)` for same seed
✅ Pedigree generates 10 samples with CHILD_001, FATHER_001, MOTHER_001 + SAMPLE_NNNN
✅ MemoryTracker imports successfully
✅ `pytest tests/performance/ --collect-only` loads conftest without errors (7 tests collected from legacy files)
✅ No private identifiers found in committed code (grep returned only docstring explanations)

## Success Criteria Met

- [x] pytest-benchmark is installed and importable
- [x] Synthetic data fixtures generate DataFrames matching variantcentrifuge column structure at 100/1K/10K/50K variants
- [x] Pedigree fixtures produce trio + cohort structures
- [x] Memory tracker provides tracemalloc-based peak tracking with warning-only budget enforcement
- [x] All identifiers are fully synthetic (SAMPLE_NNNN, GENE_NNNN, CHILD_001, etc.)
- [x] benchmarks/ directory is gitignored for ephemeral JSON results

## Next Phase Readiness

**Plans 06-02 through 06-06 are now fully unblocked:**

All subsequent benchmark plans depend on these artifacts:
- Fixtures available via `synthetic_variants`, `synthetic_pedigree`, `synthetic_scoring_config`, `synthetic_gene_burden_data`
- Memory tracking available via `memory_tracker` and `warn_if_over_budget`
- pytest-benchmark plugin active for `benchmark` fixture
- Reproducible data generation with configurable sizes and seeding

**No blockers or concerns.**

The framework is ready for benchmark implementation in inheritance analysis (06-02), compound het (06-03), gene burden (06-04), genotype replacement (06-05), and scoring (06-06).

## Implementation Notes

**Coexistence with legacy performance code:**
- `tests/performance/` already contained `benchmark_pipeline.py`, `test_streaming_parallel_performance.py`, and `README.md` (pre-GSD legacy scripts)
- These were NOT modified or removed - new pytest-benchmark framework coexists alongside them
- Legacy scripts remain as standalone benchmarking tools, new framework provides pytest-integrated profiling

**Chromosome weight normalization:**
- Initial chromosome weights didn't sum to exactly 1.0 due to floating point precision
- Fixed by normalizing with `chrom_weights / chrom_weights.sum()` to ensure exact probability sum
- This prevents `numpy.random.choice` ValueError while maintaining realistic distribution

**GT column format differences:**
- `generate_synthetic_variants`: Comma-separated `"0/0,0/1,1/1"` (inheritance/scoring format)
- `generate_synthetic_gene_burden_data`: Semicolon-separated with sample IDs `"SAMPLE_0001(0/1);SAMPLE_0002(0/0)"` (gene burden format)
- Different formats match what each analysis module expects

## Files Changed

**Created:**
- `benchmarks/.gitignore` (3 lines)
- `tests/performance/__init__.py` (1 line)
- `tests/performance/conftest.py` (158 lines)
- `tests/performance/helpers/__init__.py` (1 line)
- `tests/performance/helpers/synthetic_data.py` (308 lines)
- `tests/performance/helpers/memory_budgets.py` (102 lines)

**Modified:**
- `pyproject.toml` (+1 line in dev dependencies)
- `.gitignore` (+2 lines for benchmarks)

**Total new code:** 573 lines (helpers + conftest)

## Commits

| Task | Commit | Message |
|------|--------|---------|
| 1 | b3fa220 | chore(06-01): add pytest-benchmark and create performance test structure |
| 2 | 2221843 | feat(06-01): create synthetic data generators and benchmark fixtures |
