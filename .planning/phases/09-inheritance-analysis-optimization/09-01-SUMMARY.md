---
phase: 09
plan: 01
subsystem: inheritance-analysis-validation
tags: [golden-files, validation, testing, vectorization-prep]
requires:
  - "Phase 8: DataFrame optimization infrastructure"
provides:
  - "Golden reference files for inheritance analysis validation"
  - "Automated comparison script for clinical equivalence testing"
  - "Pytest-based golden file validation suite"
affects:
  - "09-02: Deducer vectorization (uses golden files for validation)"
  - "09-03: Compound het vectorization (uses golden files for validation)"
  - "09-04: Full vectorization validation (uses golden files for validation)"
tech-stack:
  added: []
  patterns:
    - "Golden file testing for clinical correctness"
    - "Deterministic synthetic test data generation"
    - "Parquet-based reference storage"
key-files:
  created:
    - scripts/validate_inheritance.py
    - tests/test_inheritance/test_golden_files.py
    - tests/fixtures/golden/README.md
    - tests/fixtures/golden/*.parquet (10 files)
  modified: []
decisions:
  - decision: "Use Parquet format for golden files"
    rationale: "Preserves exact dtypes and column ordering, more reliable than CSV/TSV"
    alternatives: ["CSV", "JSON", "pickle"]
    impact: "Binary files in git, but small (<10KB each) and critical for correctness"
  - decision: "Generate both .parquet and .tsv files"
    rationale: "Parquet for automated comparison, TSV for human inspection"
    alternatives: ["Parquet only", "TSV only"]
    impact: "TSV files gitignored, parquet files committed"
  - decision: "Build 10 synthetic test scenarios"
    rationale: "Cover all major inheritance patterns: trio, single-sample, extended-family, X-linked, mitochondrial, compound-het, edge cases"
    alternatives: ["Use real VCF data", "Fewer scenarios"]
    impact: "Deterministic, reproducible, fast tests; no dependency on external data"
  - decision: "Compare clinical fields with tolerance"
    rationale: "Inheritance_Pattern must match exactly, confidence within epsilon=0.001, all_patterns as sets"
    alternatives: ["Byte-for-byte comparison", "Looser comparison"]
    impact: "Robust to JSON key ordering, floating-point precision, but strict on clinical output"
metrics:
  duration: "6 minutes"
  completed: "2026-02-14"
---

# Phase 09 Plan 01: Golden File Validation Infrastructure

**One-liner:** Created golden reference files and comparison framework for validating inheritance analysis vectorization preserves clinical correctness

## Objective

Before vectorizing inheritance analysis code (Phase 9 Plans 02-04), establish a validation framework that proves the optimized code produces **clinically equivalent output** to the current implementation.

## What Was Built

### 1. Golden File Generation Script (`scripts/validate_inheritance.py`)

**Two-mode CLI script:**

**Generate mode:** `python scripts/validate_inheritance.py generate`
- Builds 10 deterministic synthetic test scenarios
- Runs `analyze_inheritance()` with current (pre-vectorization) implementation
- Saves results as `.parquet` (exact dtypes) and `.tsv` (human-readable)
- Creates golden reference files in `tests/fixtures/golden/`

**Compare mode:** `python scripts/validate_inheritance.py compare`
- Loads golden `.parquet` files
- Re-runs `analyze_inheritance()` on same input data
- Compares outputs:
  - `Inheritance_Pattern`: exact match
  - `Inheritance_Details` JSON:
    - `primary_pattern`: exact match
    - `all_patterns`: set equality (order-independent)
    - `confidence`: within epsilon=0.001
    - `samples_with_pattern`: sample IDs must match
- Reports pass/fail for each scenario
- Exit code 0 if all pass, 1 if any fail

### 2. Test Scenarios (10 scenarios covering all major patterns)

| Scenario | Description | Key Pattern | Variants |
|----------|-------------|-------------|----------|
| `trio_denovo` | De novo variant | Child het, parents ref | 1 |
| `trio_dominant` | Autosomal dominant | Child het, affected father het | 1 |
| `trio_recessive` | Autosomal recessive | Child hom_alt, both parents het | 1 |
| `trio_denovo_candidate` | De novo with missing GT | One parent missing genotype | 1 |
| `single_sample` | No pedigree data | Het → "unknown", hom_alt → "homozygous" | 2 |
| `extended_family` | Multi-generation pedigree | Dominant segregation 3 generations | 1 |
| `x_linked` | X-linked recessive | Male hemizygous, mother carrier | 1 |
| `mitochondrial` | Mitochondrial inheritance | Maternal transmission (CHROM=MT) | 1 |
| `compound_het` | Compound heterozygous | Trans configuration (one from each parent) | 2 |
| `edge_cases` | Edge case handling | Missing genotypes, single variants | 2 |

**Total:** 13 variants across 10 scenarios

### 3. Pytest Integration (`tests/test_inheritance/test_golden_files.py`)

**5 test functions with comprehensive coverage:**

1. **`test_golden_file_match[scenario_name]`** (parametrized, 10 instances)
   - Loads golden file for scenario
   - Runs `analyze_inheritance()` with current code
   - Compares `Inheritance_Pattern` and `Inheritance_Details`
   - Asserts clinical equivalence

2. **`test_scenario_determinism`**
   - Runs each scenario twice
   - Asserts identical output (proves no randomness)

3. **`test_scenario_data_types`**
   - Validates output structure (columns exist, correct dtypes)
   - Validates `Inheritance_Details` is valid JSON

4. **`test_all_scenarios_covered`**
   - Ensures all 10 expected scenarios are present

5. **`test_golden_files_readable`**
   - Verifies all golden `.parquet` files can be loaded
   - Checks for required columns

**Test markers:** `@pytest.mark.unit` and `@pytest.mark.inheritance`

### 4. Documentation (`tests/fixtures/golden/README.md`)

- Explains purpose and file format
- Documents all 10 scenarios
- Provides usage instructions for generate/compare modes
- Describes comparison criteria (what must match, what has tolerance)
- Maintenance guidelines (when to regenerate, when to investigate failures)

## Verification Results

All verification criteria met:

✅ `python scripts/validate_inheritance.py generate` — Created 10 golden files without errors
✅ `python scripts/validate_inheritance.py compare` — Exits with code 0 (all scenarios pass)
✅ `pytest tests/test_inheritance/test_golden_files.py -v` — All 14 tests pass
✅ `pytest -m unit` — No regressions (599 passed)
✅ Golden files exist for all 10 scenarios (trio, single-sample, extended-family, x-linked, mitochondrial, compound-het, edge-cases)

## Technical Implementation

### Scenario Builder Functions

Each scenario is a pure function returning `(df, pedigree_data, sample_list, scenario_name)`:
- `build_trio_denovo_scenario()`
- `build_trio_dominant_scenario()`
- `build_trio_recessive_scenario()`
- `build_trio_denovo_candidate_scenario()`
- `build_single_sample_scenario()`
- `build_extended_family_scenario()`
- `build_x_linked_scenario()`
- `build_mitochondrial_scenario()`
- `build_compound_het_scenario()`
- `build_edge_cases_scenario()`

Builders are defined in `scripts/validate_inheritance.py` and imported by `test_golden_files.py` for data reuse.

### Comparison Logic

**Stable sort key:** `(CHROM, POS, REF, ALT)` — ensures row-order independence

**Field-by-field comparison:**
```python
# Exact match
assert golden_pattern == result_pattern

# Set equality (order-independent)
assert set(golden_all_patterns) == set(result_all_patterns)

# Floating-point tolerance
assert abs(golden_confidence - result_confidence) <= 0.001

# Sample ID set equality
assert set(golden_sample_ids) == set(result_sample_ids)
```

### Determinism Verification

`test_scenario_determinism` runs each scenario twice and compares outputs:
- No random seed dependencies
- No timestamp dependencies
- No filesystem ordering dependencies
- Fully reproducible across runs

## Deviations from Plan

None — Plan executed exactly as written.

## Challenges Encountered

1. **TSV files gitignored:** `.gitignore` excludes `*.tsv` files.
   - **Solution:** Committed only `.parquet` files (source of truth), TSV files are regenerated locally.

2. **Path handling for imports:** Test file needs to import from `scripts/`.
   - **Solution:** Added `sys.path.insert(0, str(project_root / "scripts"))` to enable imports.

## Performance

- **Generate mode:** ~2 seconds (creates 10 golden files)
- **Compare mode:** ~2 seconds (validates 10 scenarios)
- **Pytest run:** ~1.8 seconds (14 tests)
- **Total plan execution:** 6 minutes (including git operations, verification)

## Next Phase Readiness

**Phase 9 Plan 02 (Deducer Vectorization):**
- ✅ Golden files ready for validation
- ✅ Comparison script available via `python scripts/validate_inheritance.py compare`
- ✅ Pytest integration for CI/CD
- ✅ Documentation explains when to regenerate vs investigate failures

**Future plans (09-03, 09-04):**
- Can reuse same golden files for compound het and full vectorization validation
- If vectorized code changes clinical behavior (intentionally), regenerate golden files after verifying correctness

## Example Usage

```bash
# Generate golden files (one-time, or after fixing clinical bugs)
python scripts/validate_inheritance.py generate

# Validate after making changes to inheritance code
python scripts/validate_inheritance.py compare

# Run as part of test suite
pytest tests/test_inheritance/test_golden_files.py -v

# Check determinism
pytest tests/test_inheritance/test_golden_files.py::test_scenario_determinism -v
```

## Files Changed

**Created:**
- `scripts/validate_inheritance.py` (925 lines)
- `tests/test_inheritance/test_golden_files.py` (260 lines)
- `tests/fixtures/golden/README.md` (documentation)
- `tests/fixtures/golden/*.parquet` (10 golden reference files)

**Modified:** None

## Test Coverage

- **New tests added:** 14 (5 functions, 10 parametrized instances)
- **Test markers:** `unit`, `inheritance`
- **All tests pass:** ✅ 14/14

## Conclusion

Successfully created a robust validation framework for inheritance analysis vectorization. The golden file infrastructure provides:

1. **Confidence:** Can prove vectorized code preserves clinical correctness
2. **Speed:** Fast validation (~2 seconds) enables rapid iteration
3. **Reproducibility:** Deterministic tests prevent false failures
4. **Documentation:** Clear guidelines for maintenance and troubleshooting

**Ready to proceed with Phase 9 Plan 02 (Deducer Vectorization).**
