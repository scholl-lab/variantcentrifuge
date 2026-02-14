---
phase: 06-benchmark-framework
plan: 03
subsystem: testing
tags: [pytest, benchmarking, ratio-assertions, tracemalloc, memory-profiling]
requires:
  - phase: 06-01
    provides: synthetic data generators, MemoryTracker helper, pytest-benchmark framework
provides:
  - Ratio assertion tests for vectorized vs sequential implementations (zero-flakiness comparisons)
  - Memory budget enforcement tests via tracemalloc (warning-only violations)
affects:
  - 06-02, 06-04, 06-05, 06-06 (can use ratio assertion and memory profiling patterns)
tech-stack:
  added: []
  patterns: [ratio-assertions, within-test-timing, tracemalloc-profiling, warning-only-budgets]
key-files:
  created:
    - tests/performance/benchmark_ratio_assertions.py
    - tests/performance/benchmark_memory_budgets.py
  modified: []
key-decisions:
  - "Ratio assertions use time.perf_counter() directly, not pytest-benchmark fixture (enables within-test comparison)"
  - "Memory budget violations are warnings only, never hard test failures"
  - "Memory tests separate from timing tests to avoid tracemalloc overhead contamination"
patterns-established:
  - "Ratio pattern: Time both implementations in same test, compute speedup = sequential/vectorized, assert > 1.0"
  - "Memory pattern: MemoryTracker context manager → warn_if_over_budget → print peak_mb → correctness assertions only"
duration: 268 seconds (~4.5 minutes)
completed: 2026-02-14
---

# Phase 06 Plan 03: Ratio Assertions & Memory Budgets Summary

**Ratio assertions compare vectorized vs sequential comp_het within same run (4.6x speedup observed), memory profiling via tracemalloc for 5 components with warning-only budget enforcement**

## Performance

- **Duration:** 4.5 minutes (268 seconds)
- **Started:** 2026-02-14T07:08:59Z
- **Completed:** 2026-02-14T07:13:27Z
- **Tasks:** 2
- **Files modified:** 2 created (plus 3 unplanned files committed accidentally)

## Accomplishments
- Ratio assertion tests enable zero-flakiness comparison of vectorized vs sequential implementations by timing both in same test run
- Memory budget tests measure peak memory via tracemalloc for inheritance, comp_het, gene_burden, scoring, and DataFrame I/O
- Speedup measurements show vectorized comp_het is 4.6x faster than sequential on 18-variant gene
- Memory profiling infrastructure enables regression detection without causing hard test failures

## Task Commits

Each task was committed atomically:

1. **Task 1: Create ratio assertion tests** - `64d99ca` (test)
2. **Task 2: Create memory budget tests** - `27dd3ac` (test)

**Note:** Task 2 commit accidentally included 3 additional files (benchmark_comp_het.py, benchmark_genotype_replacement.py, benchmark_inheritance.py) that belong to future plans 06-02, 06-04, 06-05. See Deviations section.

## Files Created/Modified

**Created:**
- `tests/performance/benchmark_ratio_assertions.py` (174 lines) - Ratio assertions comparing vectorized vs sequential compound het detection
- `tests/performance/benchmark_memory_budgets.py` (214 lines) - Memory budget enforcement tests with tracemalloc

**Accidentally committed (not part of this plan):**
- `tests/performance/benchmark_comp_het.py` (164 lines) - Belongs to plan 06-02
- `tests/performance/benchmark_genotype_replacement.py` (174 lines) - Belongs to plan 06-05
- `tests/performance/benchmark_inheritance.py` (167 lines) - Belongs to plan 06-02

## Decisions Made

1. **Ratio assertions use `time.perf_counter()` directly instead of pytest-benchmark fixture**
   - Rationale: Ratio tests need to time BOTH implementations in the same test to ensure identical conditions. The `benchmark` fixture can only time one function per test.
   - Benefit: Zero flakiness - both implementations measured on same machine, same process, same test run. No cross-run variance.

2. **Memory budget violations are warnings only, never hard failures**
   - Rationale: Per CONTEXT.md decision, memory budgets should flag concerning usage but not block development.
   - Implementation: `warn_if_over_budget()` issues `warnings.warn()`, tests have zero `assert` statements on memory values.

3. **Memory tests separated from timing tests**
   - Rationale: tracemalloc adds 5-20% overhead that would contaminate timing measurements.
   - Implementation: Separate files (benchmark_memory_budgets.py vs benchmark_ratio_assertions.py), no `benchmark` fixture in memory tests.

## Deviations from Plan

### Unplanned Files Committed

**[Accidental Inclusion] Three future-plan benchmark files**
- **Found during:** Task 2 commit (git add picked up untracked files)
- **Issue:** `git add tests/performance/benchmark_memory_budgets.py` also staged three other benchmark files that shouldn't exist yet:
  - `tests/performance/benchmark_comp_het.py` (plan 06-02)
  - `tests/performance/benchmark_genotype_replacement.py` (plan 06-05)
  - `tests/performance/benchmark_inheritance.py` (plan 06-02)
- **Root cause:** Files were created at some point (possibly by IDE autocomplete or previous session) and were sitting untracked in working directory.
- **Impact:** Files are already committed in 27dd3ac. They contain substantial code (164-174 lines each) but are NOT part of this plan's scope.
- **Resolution:** Document in SUMMARY. These files will be reviewed/validated when their respective plans execute. If they don't match plan requirements, they'll be rewritten.
- **Committed in:** 27dd3ac (Task 2 commit)

---

**Total deviations:** 1 accidental file inclusion (3 files committed that belong to future plans)
**Impact on plan:** Plan 06-03 objectives fully met (ratio assertions and memory budgets complete). The extra files don't affect this plan's deliverables but create technical debt for plans 06-02 and 06-05 to validate/rewrite as needed.

## Issues Encountered

### Implementation Fixes

**Issue 1: Incorrect column name in inheritance test**
- **Problem:** Test checked for "INHERITANCE_PATTERN" column but analyzer returns "Inheritance_Pattern"
- **Resolution:** Fixed assertion to match actual column name
- **Verification:** Test passes, assertion finds expected column

**Issue 2: Gene burden test used wrong data structure**
- **Problem:** Initial test passed gene burden data with GT column format, but `perform_gene_burden_analysis` expects pre-aggregated data with proband_count, control_count, etc.
- **Resolution:** Created pre-aggregated gene burden DataFrame matching expected schema (50 genes with count columns)
- **Verification:** Test passes, gene burden analysis completes successfully

**Issue 3: Scoring config structure mismatch**
- **Problem:** Initial config had formulas as dict instead of list of dicts
- **Resolution:** Changed to `formulas: [{"score_name": "formula"}, ...]` format matching actual API
- **Verification:** Test passes, scoring applies formulas correctly

## Test Results

**Ratio Assertions (benchmark_ratio_assertions.py):**
- `test_comp_het_vectorized_vs_original_ratio`: Speedup 4.62x (1.58ms sequential → 0.34ms vectorized, 18 variants)
- `test_comp_het_ratio_at_scale[500]`: Speedup 3.03x (1.32ms → 0.44ms, 21 variants)
- `test_comp_het_ratio_at_scale[2000]`: Speedup 3.22x (1.62ms → 0.50ms, 18 variants)

**Memory Budgets (benchmark_memory_budgets.py):**
- Inheritance analysis: 7.5 MB peak (10K variants) — well under 512 MB budget
- Compound het: 1.0 MB peak (5K variants single gene) — well under 256 MB budget
- Gene burden: 0.2 MB peak (50 genes, 100 samples) — well under 512 MB budget
- Scoring: 1.4 MB peak (10K variants, 10 samples) — well under 256 MB budget
- DataFrame I/O: 13.5 MB peak (50K variants, 4.6 MB file) — baseline measurement, no budget

All tests pass with zero warnings (all components well within budget).

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Plans 06-02 through 06-06 can now use:**
- Ratio assertion pattern for comparing implementations (see `_time_function_iterations` helper)
- Memory profiling pattern for budget enforcement (see MemoryTracker + warn_if_over_budget)
- Established that current memory usage is very low (<15 MB peak across all components)

**Note on accidentally committed files:**
- Plans 06-02 and 06-05 will need to validate whether the accidentally committed benchmark files (benchmark_comp_het.py, benchmark_genotype_replacement.py, benchmark_inheritance.py) meet their requirements.
- If files don't match plan specs, those plans should rewrite them.
- If files are satisfactory, those plans can mark them as complete and commit only documentation/SUMMARY updates.

**No blockers.**

---
*Phase: 06-benchmark-framework*
*Completed: 2026-02-14*
