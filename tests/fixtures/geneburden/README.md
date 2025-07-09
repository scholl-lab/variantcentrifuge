# Gene Burden Test Fixtures

This directory contains comprehensive test data and infrastructure for testing all aspects of VariantCentrifuge's gene burden analysis functionality.

## Overview

The gene burden test framework provides complete coverage for all possible ways to specify case/control groups in VariantCentrifuge:

### ✅ **Direct Sample Specification**
- `--case-samples` + `--control-samples`
- `--case-samples-file` + `--control-samples-file`

### ✅ **Phenotype-based Classification**
- `--case-phenotypes` + `--control-phenotypes`
- `--case-phenotypes-file` + `--control-phenotypes-file`
- Integration with `--phenotype-file` + `--phenotype-sample-column` + `--phenotype-value-column`

### ✅ **Advanced Scenarios**
- Multiple phenotype file formats (CSV, TSV)
- Alternative column names
- Multiple HPO terms per sample
- Error handling for missing data
- Edge case testing

## Files Generated

### Core Scripts
- **`generate_comprehensive_test_data.py`** - Main test data generator
- **`verify_test_data.py`** - Verification script for generated data
- **`comprehensive_test_data/`** - Generated test dataset directory

### Test Data Structure
```
comprehensive_test_data/
├── test_data.vcf.gz              # Main VCF with controlled genotypes
├── case_samples.txt              # Case sample IDs
├── control_samples.txt           # Control sample IDs
├── test_genes.txt               # All test genes
├── phenotypes_basic.csv         # GCKD-style phenotype file
├── phenotypes_extended.csv      # Multiple HPO terms per sample
├── phenotypes_alt_columns.csv   # Alternative column names
├── case_hpo_terms.txt          # HPO terms for cases
├── control_hpo_terms.txt       # HPO terms for controls
├── run_comprehensive_tests.sh   # Test runner script
└── README.md                    # Detailed usage instructions
```

## Quick Start

### 1. Generate Test Data
```bash
cd tests/fixtures/geneburden
python generate_comprehensive_test_data.py --output-dir comprehensive_test_data
```

### 2. Verify Test Data
```bash
python verify_test_data.py
```

### 3. Run Tests
```bash
# Pytest integration
pytest tests/test_gene_burden_comprehensive.py -v

# Or manual testing
cd comprehensive_test_data
./run_comprehensive_tests.sh
```

## Test Scenarios

### Scenario 1: Direct Sample Lists
```bash
variantcentrifuge \
  --vcf-file test_data.vcf.gz \
  --gene-file test_genes.txt \
  --case-samples CASE_001,CASE_002,CASE_003 \
  --control-samples CTRL_001,CTRL_002,CTRL_003 \
  --perform-gene-burden \
  --use-new-pipeline
```

### Scenario 2: Sample Files
```bash
variantcentrifuge \
  --vcf-file test_data.vcf.gz \
  --gene-file test_genes.txt \
  --case-samples-file case_samples.txt \
  --control-samples-file control_samples.txt \
  --perform-gene-burden \
  --use-new-pipeline
```

### Scenario 3: HPO-based Classification
```bash
variantcentrifuge \
  --vcf-file test_data.vcf.gz \
  --gene-file test_genes.txt \
  --phenotype-file phenotypes_basic.csv \
  --phenotype-sample-column SampleID \
  --phenotype-value-column identifier \
  --case-phenotypes HP:0000113,HP:0000003,HP:0000107 \
  --control-phenotypes HP:0000001,HP:0032101 \
  --perform-gene-burden \
  --use-new-pipeline
```

### Scenario 4: HPO Term Files
```bash
variantcentrifuge \
  --vcf-file test_data.vcf.gz \
  --gene-file test_genes.txt \
  --phenotype-file phenotypes_basic.csv \
  --phenotype-sample-column SampleID \
  --phenotype-value-column identifier \
  --case-phenotypes-file case_hpo_terms.txt \
  --control-phenotypes-file control_hpo_terms.txt \
  --perform-gene-burden \
  --use-new-pipeline
```

## Expected Results

The test dataset is designed with controlled probabilities to ensure predictable results:

### ✅ **Disease Genes** (should show enrichment)
- **PKD1, PKD2** - Polycystic kidney disease genes
- **BRCA1, BRCA2** - Breast cancer genes
- Expected: OR > 1.0, p-values < 0.05

### ✅ **Control Genes** (should show no enrichment)
- **TTN, OBSCN, MUC16** - Large genes with many benign variants
- Expected: OR ≈ 1.0, p-values > 0.05

## Data Composition

- **100 samples**: 40 cases + 60 controls
- **19 variants** across 8 genes
- **Controlled genotype distributions**:
  - 75% cases have pathogenic disease gene variants
  - 10% controls have pathogenic disease gene variants
  - 40% all samples have control gene variants

## HPO Terms Used

### Case Phenotypes (Disease-related)
- `HP:0000113` - Polycystic kidney dysplasia
- `HP:0000003` - Multicystic kidney dysplasia
- `HP:0000107` - Renal cyst
- `HP:0000822` - Hypertension
- `HP:0003774` - Stage 5 chronic kidney disease

### Control Phenotypes (Normal/Healthy)
- `HP:0000001` - All (root term)
- `HP:0032101` - Normal phenotype

## Integration with Test Suite

### Pytest Markers
```bash
# Run all gene burden tests
pytest -m gene_burden

# Run slow integration tests
pytest -m "gene_burden and slow"

# Run specific test class
pytest tests/test_gene_burden_comprehensive.py::TestGeneBurdenDirectSamples -v
```

### Test Classes
- **`TestGeneBurdenDirectSamples`** - Direct sample specification tests
- **`TestGeneBurdenPhenotypeBased`** - Phenotype-based classification tests
- **`TestGeneBurdenValidation`** - Result validation tests
- **`TestGeneBurdenErrorHandling`** - Error handling and edge cases

## Maintenance

### Regenerating Test Data
If you need to update the test data:

```bash
cd tests/fixtures/geneburden
python generate_comprehensive_test_data.py --output-dir comprehensive_test_data --seed 42
```

### Adding New Test Scenarios
1. Modify `generate_comprehensive_test_data.py` to add new data patterns
2. Add corresponding test methods in `test_gene_burden_comprehensive.py`
3. Update this README with new scenarios

### Custom Base VCF
To use a real VCF file as the base for test data generation:

```bash
python generate_comprehensive_test_data.py \
  --output-dir comprehensive_test_data \
  --base-vcf /path/to/annotated.vcf.gz \
  --seed 42
```

## Troubleshooting

### Test Data Not Found
```bash
# Generate the test data first
cd tests/fixtures/geneburden
python generate_comprehensive_test_data.py --output-dir comprehensive_test_data
```

### VCF Index Issues
```bash
# Recreate the tabix index
cd comprehensive_test_data
tabix -p vcf test_data.vcf.gz
```

### Verification Failures
```bash
# Run verification to identify issues
python verify_test_data.py comprehensive_test_data
```

## Related Files

- **`../../test_gene_burden_comprehensive.py`** - Main pytest test file
- **`../../../testing/tools/test_data_generators/`** - Original test data generators
- **`../../../testing/phenotypes/`** - Real phenotype file examples

---

This comprehensive test framework ensures that all gene burden analysis functionality is thoroughly tested across all supported input methods and edge cases.