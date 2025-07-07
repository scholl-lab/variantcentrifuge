# Regression Testing Framework

This directory contains the regression testing framework for validating that the new refactored pipeline produces identical results to the original pipeline.

## Overview

The regression testing ensures that the refactoring maintains 100% compatibility with the original pipeline by:
- Running both pipelines with identical inputs
- Comparing outputs byte-by-byte
- Validating statistics and metadata
- Testing various configurations and edge cases

## Components

### 1. `test_regression_suite.py`
The main pytest test suite that:
- Defines test configurations covering various use cases
- Runs both old and new pipelines
- Validates output consistency
- Reports any differences

### 2. `generate_baseline_outputs.py`
Script to generate baseline outputs using the old pipeline:
```bash
./generate_baseline_outputs.py test_data.vcf.gz --output-dir baseline_outputs
```

### 3. `compare_outputs.py`
Utility for detailed comparison of pipeline outputs:
```bash
./compare_outputs.py old_outputs/ new_outputs/ --verbose
```

## Running Regression Tests

### Prerequisites
1. Install VariantCentrifuge with development dependencies
2. Ensure bcftools, snpEff, SnpSift, and bedtools are in PATH
3. Have a test VCF file available

### Quick Start

1. Generate baseline outputs:
```bash
python tests/regression/generate_baseline_outputs.py test_cohort.vcf.gz
```

2. Run regression tests:
```bash
pytest tests/regression/test_regression_suite.py -v -m regression
```

3. Compare specific outputs:
```bash
python tests/regression/compare_outputs.py baseline_outputs/outputs/ test_outputs/ --verbose
```

## Test Configurations

The regression suite tests the following scenarios:

1. **Basic Gene Extraction**
   - Single gene (BRCA1)
   - Multiple genes from file

2. **Filtering Presets**
   - Rare coding variants
   - Pathogenic variants
   - Complex filter combinations

3. **Advanced Features**
   - Variant scoring
   - Inheritance analysis with PED files
   - Custom BED annotations
   - Final filtering on computed columns

4. **Edge Cases**
   - Empty results
   - Large gene lists
   - Complex filter expressions

## Interpreting Results

### Success
All tests pass when outputs are identical:
```
✓ Regression test passed: single_gene_basic
✓ Regression test passed: gene_with_scoring
...
```

### Failures
Differences are reported with details:
```
TSV outputs differ for single_gene_rare_coding:
Column 'Score' has 5 differences
  Row 10: '0.85' → '0.86'
  Row 25: '0.92' → '0.93'
```

### Common Issues

1. **Column Order**: The new pipeline may output columns in a different order. This is acceptable if all data is present.

2. **Variant ID**: The `Variant_ID` column may differ due to processing order changes. This column is ignored in comparisons.

3. **Floating Point**: Small differences in scores due to floating point precision are acceptable within tolerance (1e-6).

## Adding New Tests

To add a new regression test:

1. Add configuration to `REGRESSION_TESTS` in `test_regression_suite.py`:
```python
RegressionTestConfig(
    name="my_new_test",
    gene_name="GENE1",
    preset="rare",
    extra_args=["--extra-flag", "value"],
    description="Test description",
)
```

2. Regenerate baselines:
```bash
./generate_baseline_outputs.py test_data.vcf.gz
```

3. Run the new test:
```bash
pytest tests/regression/test_regression_suite.py::test_pipeline_output_match[my_new_test] -v
```

## Troubleshooting

### Missing Tools
If tests fail with "command not found":
```bash
# Check tool availability
which bcftools snpEff SnpSift bedtools

# Install missing tools with conda/mamba
mamba install -c bioconda bcftools snpeff snpsift bedtools
```

### File Not Found
Ensure test data exists:
- Place test VCF in `tests/data/test_cohort.vcf.gz`
- Or specify custom location with pytest fixtures

### Debugging Failures
1. Check logs in test output directory
2. Use `--verbose` flag for detailed output
3. Compare specific files with `compare_outputs.py`
4. Examine intermediate files in `tmp_path` directories

## Performance Considerations

- Tests run with `--threads 1` for deterministic results
- Full regression suite may take several minutes
- Use pytest markers to run specific tests:
  ```bash
  pytest -m "regression and not slow"
  ```

## CI Integration

The regression tests can be integrated into CI pipelines:

```yaml
# Example GitHub Actions workflow
- name: Run Regression Tests
  run: |
    pytest tests/regression/test_regression_suite.py \
      -v -m regression \
      --tb=short \
      --junit-xml=regression-results.xml
```