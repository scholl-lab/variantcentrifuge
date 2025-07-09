# Comprehensive Gene Burden Test Dataset

Generated on: 2025-07-09T14:18:13.467555

## Overview

This test dataset provides comprehensive coverage for testing all gene burden analysis 
functionality in VariantCentrifuge, including all possible ways to specify case/control groups.

## Dataset Composition

- **Total Samples**: 100 (40 cases, 60 controls)
- **Disease Genes**: BRCA1, BRCA2, PKD1, PKD2
- **Control Genes**: MUC16, OBSCN, TTN
- **Variants**: 22+ variants using **real hg19 genomic coordinates**

### Gene Coordinates (hg19)
- **PKD1**: chr16:2,138,710-2,185,899 (polycystic kidney disease)
- **PKD2**: chr4:88,928,819-88,998,929 (polycystic kidney disease) 
- **BRCA1**: chr17:41,196,311-41,277,500 (breast cancer)
- **BRCA2**: chr13:32,889,610-32,973,805 (breast cancer)
- **TTN**: chr2:179,390,715-179,695,529 (titin, large control gene)
- **OBSCN**: chr1:228,395,830-228,566,577 (obscurin, large control gene)
- **MUC16**: chr19:8,959,519-9,092,018 (mucin, control gene)

## Files Generated

### Core Data Files
- `test_data.vcf.gz` - Main VCF file with controlled genotype distributions
- `test_data.vcf.gz.tbi` - Tabix index

### Sample Specification Files
- `case_samples.txt` - Complete list of case sample IDs
- `control_samples.txt` - Complete list of control sample IDs
- `case_samples_subset.txt` - Subset of case samples (for testing)
- `control_samples_subset.txt` - Subset of control samples (for testing)

### Phenotype Files
- `phenotypes_basic.csv` - Basic phenotype file (GCKD-like format)
- `phenotypes_extended.csv` - Extended with multiple HPO terms per sample
- `phenotypes.tsv` - Simple TSV format
- `phenotypes_alt_columns.csv` - Alternative column names
- `phenotypes_incomplete.csv` - Missing some samples (for error testing)

### HPO Term Files
- `case_hpo_terms.txt` - HPO terms defining case phenotypes
- `control_hpo_terms.txt` - HPO terms defining control phenotypes  
- `mixed_hpo_terms.txt` - Mixed terms for edge case testing

### Gene Lists
- `test_genes.txt` - All test genes
- `disease_genes.txt` - Disease-associated genes only
- `control_genes.txt` - Control genes only

### Test Infrastructure
- `run_comprehensive_tests.sh` - Main test runner script
- `test_1_direct_samples.sh` through `test_5_alt_columns.sh` - Individual test scripts
- `test_configurations.json` - Test scenario configurations
- `dataset_statistics.json` - Detailed dataset statistics

## Test Scenarios

### 1. Direct Sample Specification
```bash
--case-samples CASE_001,CASE_002,CASE_003
--control-samples CTRL_001,CTRL_002,CTRL_003
```

### 2. Sample File Specification  
```bash
--case-samples-file case_samples.txt
--control-samples-file control_samples.txt
```

### 3. HPO-based Classification
```bash
--phenotype-file phenotypes_basic.csv
--phenotype-sample-column SampleID
--phenotype-value-column identifier
--case-phenotypes HP:0000113,HP:0000003,HP:0000107
--control-phenotypes HP:0000001,HP:0032101
```

### 4. HPO Term Files
```bash
--phenotype-file phenotypes_basic.csv
--phenotype-sample-column SampleID
--phenotype-value-column identifier
--case-phenotypes-file case_hpo_terms.txt
--control-phenotypes-file control_hpo_terms.txt
```

### 5. Alternative Column Names
```bash
--phenotype-file phenotypes_alt_columns.csv
--phenotype-sample-column sample_name
--phenotype-value-column hpo_id
--case-phenotypes HP:0000113,HP:0000003,HP:0000107
--control-phenotypes HP:0000001,HP:0032101
```

## Expected Results

When running gene burden analysis on this dataset:

- **Disease genes** (PKD1, PKD2, BRCA1, BRCA2) should show **significant enrichment** in cases
- **Control genes** (TTN, OBSCN, MUC16) should show **no significant enrichment**
- P-values for disease genes should be < 0.05 (before multiple testing correction)
- Odds ratios for disease genes should be > 1.0

## Running Tests

1. **Run all tests**:
   ```bash
   ./run_comprehensive_tests.sh
   ```

2. **Run individual tests**:
   ```bash
   ./test_1_direct_samples.sh
   ./test_2_sample_files.sh
   # etc.
   ```

3. **Custom test**:
   ```bash
   variantcentrifuge \
     --vcf-file test_data.vcf.gz \
     --gene-file test_genes.txt \
     --case-samples-file case_samples.txt \
     --control-samples-file control_samples.txt \
     --perform-gene-burden \
     --preset rare,coding \
     --use-new-pipeline \
     --output-file results.tsv
   ```

## HPO Terms Used

### Case Phenotypes (Disease-related)
- HP:0000113: Polycystic kidney dysplasia
- HP:0000003: Multicystic kidney dysplasia
- HP:0000107: Renal cyst
- HP:0000822: Hypertension
- HP:0003774: Stage 5 chronic kidney disease

### Control Phenotypes (Normal/Healthy)
- HP:0000001: All
- HP:0032101: Normal phenotype

### Mixed Phenotypes (For edge case testing)
- HP:0001626: Abnormality of the cardiovascular system
- HP:0000478: Abnormality of the eye
- HP:0001574: Abnormality of the integument

## Validation

Check `dataset_statistics.json` for detailed statistics on:
- Genotype distribution per gene and sample type
- Variant counts and classifications
- Sample composition
- Test scenario configurations

The dataset is designed with controlled probabilities to ensure reliable test results
while maintaining realistic genetic variant distributions.
