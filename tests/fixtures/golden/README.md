# Golden Reference Files for Inheritance Analysis

This directory contains golden reference files used to validate inheritance analysis vectorization.
These files capture the output of the **current (pre-vectorization) implementation** and serve as
the source of truth for ensuring clinical equivalence after optimization.

## Purpose

Before vectorizing the inheritance analysis code, we need a way to prove that the optimized code
produces **clinically equivalent output**. These golden files enable automated comparison between:
- The original implementation (used to generate these files)
- The vectorized implementation (compared against these files)

## File Format

Each scenario has two files:

1. **`{scenario_name}.parquet`** - Binary reference file that preserves exact dtypes
   - Used for automated comparison in tests
   - Contains full DataFrame with all columns including `Inheritance_Pattern` and `Inheritance_Details`

2. **`{scenario_name}.tsv`** - Human-readable reference file
   - Tab-separated values for manual inspection
   - Same content as parquet, but easier to review

## Test Scenarios

The golden files cover the following inheritance patterns:

| Scenario | Description | Key Pattern |
|----------|-------------|-------------|
| `trio_denovo` | Trio with de novo variant | Child het, parents ref |
| `trio_dominant` | Autosomal dominant | Child het, affected father het |
| `trio_recessive` | Autosomal recessive | Child hom_alt, both parents het |
| `trio_denovo_candidate` | De novo candidate | Child het, one parent missing GT |
| `single_sample` | Single sample (no pedigree) | Het → "unknown", hom_alt → "homozygous" |
| `extended_family` | Multi-generation pedigree | Dominant segregation across 3 generations |
| `x_linked` | X-linked recessive | Male proband hemizygous, mother carrier |
| `mitochondrial` | Mitochondrial inheritance | Maternal transmission (CHROM=MT) |
| `compound_het` | Compound heterozygous | Two het variants in trans configuration |
| `edge_cases` | Edge case handling | Missing genotypes, single variants |

## Generating Golden Files

Golden files are generated using the current (pre-vectorization) code:

```bash
python scripts/validate_inheritance.py generate
```

This will:
1. Build synthetic test DataFrames for each scenario
2. Run `analyze_inheritance()` with the current implementation
3. Save results to `.parquet` and `.tsv` files in this directory

**IMPORTANT:** Only regenerate golden files if the current implementation is known to be correct
and you intentionally want to update the reference output.

## Comparing Against Golden Files

After making changes to inheritance analysis code, validate against golden files:

```bash
python scripts/validate_inheritance.py compare
```

This will:
1. Load each golden `.parquet` file
2. Re-run `analyze_inheritance()` on the same input data
3. Compare outputs for:
   - `Inheritance_Pattern` (must match exactly)
   - `Inheritance_Details` JSON fields:
     - `primary_pattern` (exact match)
     - `all_patterns` (set equality)
     - `confidence` (within epsilon=0.001)
     - `samples_with_pattern` (sample IDs must match)
4. Report any mismatches

Exit code:
- `0` = All scenarios pass
- `1` = At least one scenario fails

## Comparison Criteria

The comparison is **clinically focused**, not byte-for-byte identical:

### Must Match Exactly
- Inheritance pattern classification
- Primary pattern selection
- Set of all possible patterns
- Sample IDs with each pattern

### Allowed Variance
- Confidence score (±0.001)
- Formatting of Inheritance_Details JSON (order of keys, whitespace)
- Non-clinical metadata fields

### Stable Sort Key
All comparisons use a stable sort key: `(CHROM, POS, REF, ALT)`

This ensures row-order independence between implementations.

## Integration with Tests

Golden file validation is also integrated into pytest via `tests/test_inheritance/test_golden_files.py`:

```bash
# Run golden file tests
pytest tests/test_inheritance/test_golden_files.py -v

# Run only determinism check
pytest tests/test_inheritance/test_golden_files.py -v -k determinism
```

## Scenario Determinism

All scenarios use **deterministic synthetic data** (no random values, no timestamps).
This ensures:
- Golden files can be regenerated reliably
- Comparisons are reproducible
- CI/CD pipelines get consistent results

To verify determinism:
```bash
pytest tests/test_inheritance/test_golden_files.py::test_scenario_determinism -v
```

## Version Control

Golden files **are committed to git** because:
1. They document expected behavior
2. They enable historical comparison (detecting regressions)
3. They allow CI/CD validation without regeneration

Parquet files are binary but small (<10KB each).

## Maintenance

- **When to regenerate:** Only when the current implementation is intentionally updated
  (e.g., fixing a clinical bug, adding a new pattern type)
- **When to investigate:** If comparison fails, determine whether:
  - The new code has a bug (fix it)
  - The golden files are outdated (regenerate after verifying correctness)
- **Adding scenarios:** Edit `scripts/validate_inheritance.py` to add new `build_*_scenario()` functions

## References

- Implementation: `scripts/validate_inheritance.py`
- Test integration: `tests/test_inheritance/test_golden_files.py`
- Source analyzer: `variantcentrifuge/inheritance/analyzer.py`
- Planning docs: `.planning/phases/09-inheritance-analysis-optimization/09-01-PLAN.md`
