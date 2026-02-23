# Enhanced Gene Burden Analysis Test Data

This directory contains enhanced test data and utilities for comprehensive gene burden analysis testing in VariantCentrifuge, featuring **realistic annotations sampled from real genomic data**.

## Overview

The enhanced gene burden test framework provides complete coverage for all possible ways to specify case/control groups in VariantCentrifuge, using **authentic SnpEff annotations** for maximum realism:

### ✅ **Key Enhancements**
- **Real Annotation Sampling**: Authentic SnpEff annotations from local anonymized templates
- **Realistic Pathogenicity**: Actual dbNSFP scores (CADD, SIFT, PolyPhen2, REVEL)
- **Diverse Variant Effects**: 10 unique effect types from real genomic data
- **Controlled Test Outcomes**: Maintains probabilistic distributions for reliable testing

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

### Directory Structure
```
geneburden/
├── generate_enhanced_test_data.py # Enhanced test data generator
├── README.md                      # This documentation
├── templates/                     # Core template files
│   └── annotation_templates.txt   # Anonymized annotation patterns
├── output/                        # Generated test datasets (gitignored)
│   ├── enhanced_test_data.vcf.gz  # VCF with real SnpEff annotations
│   ├── case_samples.txt           # Case sample IDs
│   ├── control_samples.txt        # Control sample IDs
│   ├── test_genes.txt            # All test genes
│   ├── phenotypes_*.csv          # Various phenotype file formats
│   ├── *_hpo_terms.txt          # HPO term files
│   ├── enhanced_dataset_statistics.json # Detailed statistics
│   └── README.md                  # Dataset-specific documentation
└── temp/                          # Temporary files (gitignored)
```

## Quick Start

### 1. Generate Enhanced Test Data
```bash
cd tests/fixtures/geneburden
python generate_enhanced_test_data.py
# Or specify custom output directory:
# python generate_enhanced_test_data.py --output-dir custom_output
```

### 2. Run Gene Burden Analysis
```bash
cd output
python -m variantcentrifuge.cli \
  --vcf-file enhanced_test_data.vcf.gz \
  --gene-file test_genes.txt \
  --case-samples-file case_samples.txt \
  --control-samples-file control_samples.txt \
  --perform-gene-burden \
  --preset high_or_moderate \
  --output-file results.tsv
```

## Test Scenarios

### Scenario 1: Direct Sample Lists
```bash
python -m variantcentrifuge.cli \
  --vcf-file enhanced_test_data.vcf.gz \
  --gene-file test_genes.txt \
  --case-samples CASE_001,CASE_002,CASE_003 \
  --control-samples CTRL_001,CTRL_002,CTRL_003 \
  --perform-gene-burden \
  --preset high_or_moderate
```

### Scenario 2: Sample Files
```bash
python -m variantcentrifuge.cli \
  --vcf-file enhanced_test_data.vcf.gz \
  --gene-file test_genes.txt \
  --case-samples-file case_samples.txt \
  --control-samples-file control_samples.txt \
  --perform-gene-burden \
  --preset high_or_moderate
```

### Scenario 3: HPO-based Classification
```bash
python -m variantcentrifuge.cli \
  --vcf-file enhanced_test_data.vcf.gz \
  --gene-file test_genes.txt \
  --phenotype-file phenotypes_basic.csv \
  --phenotype-sample-column SampleID \
  --phenotype-value-column identifier \
  --case-phenotypes HP:0000113,HP:0000003,HP:0000107 \
  --control-phenotypes HP:0000001,HP:0032101 \
  --perform-gene-burden \
  --preset high_or_moderate
```

### Scenario 4: HPO Term Files
```bash
python -m variantcentrifuge.cli \
  --vcf-file enhanced_test_data.vcf.gz \
  --gene-file test_genes.txt \
  --phenotype-file phenotypes_basic.csv \
  --phenotype-sample-column SampleID \
  --phenotype-value-column identifier \
  --case-phenotypes-file case_hpo_terms.txt \
  --control-phenotypes-file control_hpo_terms.txt \
  --perform-gene-burden \
  --preset high_or_moderate
```

## Expected Results

The enhanced test dataset uses realistic annotations with controlled probabilities:

### ✅ **Disease Genes** (should show enrichment)
- **PKD1, PKD2** - Polycystic kidney disease genes
- **BRCA1, BRCA2** - Breast cancer genes
- Expected: OR > 1.0, p-values < 0.05 (pathogenic variants enriched in cases)

### ✅ **Control Genes** (should show no enrichment)
- **TTN, OBSCN, MUC16** - Large genes with many benign variants
- Expected: OR ≈ 1.0, p-values > 0.05

## Enhanced Data Composition

- **100 samples**: 40 cases + 60 controls
- **25 variants** across 7 genes with **real annotations**
- **10 effect types**: `missense_variant`, `frameshift_variant`, `disruptive_inframe_deletion`, etc.
- **Authentic dbNSFP scores**: CADD phred scores, SIFT predictions, PolyPhen2 scores
- **Pathogenicity-aware genotype distributions**:
  - High CADD scores (>20) and deleterious predictions drive case enrichment
  - Moderate/low impact variants distributed more evenly

## Real Annotation Features

### Authentic Variant Effects
- **High Impact**: `frameshift_variant`, `stop_gained`
- **Moderate Impact**: `missense_variant`, `disruptive_inframe_deletion`
- **Low Impact**: `synonymous_variant`, `splice_region_variant`

### Real Pathogenicity Scores
- **CADD phred scores**: 15.45, 25.6, 33.0 (from actual variants)
- **SIFT predictions**: D (deleterious), T (tolerated)
- **PolyPhen2 predictions**: D (damaging), B (benign)
- **REVEL scores**: 0.311, 0.659 (machine learning pathogenicity)

### Population Data
- **gnomAD frequencies**: Real allele frequencies from population data
- **GERP scores**: Evolutionary conservation scores
- **Clinical significance**: Actual ClinVar annotations where available

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

## Maintenance

### Regenerating Enhanced Test Data
```bash
cd tests/fixtures/geneburden
python generate_enhanced_test_data.py \
  --seed 42
```

### Using External Annotation Sources (Optional)
```bash
python generate_enhanced_test_data.py \
  --annotation-source /path/to/your/annotated_variants.txt \
  --seed 42
```

## Validation

Check the detailed statistics to verify annotation quality:
```bash
cat output/enhanced_dataset_statistics.json | jq '.annotation_statistics'
```

Expected output shows:
- 10 unique effect types from real data
- High and moderate impact effects available
- Rich dbNSFP annotation fields preserved

## Troubleshooting

### Script Dependencies
```bash
# Verify local annotation templates exist
ls -la templates/annotation_templates.txt
```

### VCF Index Issues
```bash
cd output
tabix -p vcf enhanced_test_data.vcf.gz
```

### Annotation Quality Check
```bash
# Verify authentic annotations are present
cd output
zcat enhanced_test_data.vcf.gz | grep -v "^#" | head -1 | cut -f8
```

---

This enhanced test framework provides **realistic genomic annotations** while maintaining controlled statistical properties for reliable gene burden analysis testing.
