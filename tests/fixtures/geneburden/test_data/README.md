# Enhanced Gene Burden Test Dataset with Real Annotations

Generated on: 2025-07-09T14:54:52.026042

## Overview

This enhanced test dataset provides comprehensive coverage for testing gene burden analysis 
functionality in VariantCentrifuge, featuring **realistic annotations sampled from real genomic data**.

### Key Enhancements

- **Real Annotation Sampling**: Authentic SnpEff annotations extracted from /mnt/c/development/scholl-lab/variantcentrifuge/tests/fixtures/geneburden/../../../testdata_snpef.txt
- **Realistic Pathogenicity Prediction**: Uses actual dbNSFP scores (CADD, SIFT, PolyPhen2, REVEL)
- **Diverse Variant Effects**: 10 unique effect types from real data
- **Controlled Test Outcomes**: Maintains probabilistic genotype distributions for reliable testing

## Dataset Composition

- **Total Samples**: 100 (40 cases, 60 controls)
- **Disease Genes**: BRCA1, BRCA2, PKD1, PKD2
- **Control Genes**: MUC16, OBSCN, TTN
- **Annotation Source**: Real variants from /mnt/c/development/scholl-lab/variantcentrifuge/tests/fixtures/geneburden/../../../testdata_snpef.txt

## Enhanced Features

### Realistic Annotations
- Authentic ANN fields from SnpEff annotations
- Real dbNSFP pathogenicity predictions
- Actual genomic conservation scores (GERP, PhastCons)
- Population frequency data (gnomAD)

### Effect Type Diversity
Available effect types from real data:
- disruptive_inframe_deletion
- frameshift_variant
- intron_variant
- missense_variant
- missense_variant&splice_region_variant
- splice_donor_variant&intron_variant
- splice_region_variant&intron_variant
- stop_gained
- synonymous_variant
- upstream_gene_variant

### Pathogenicity-Aware Genotype Generation
Genotype probabilities are adjusted based on:
- Variant impact level (HIGH, MODERATE, LOW, MODIFIER)
- CADD phred scores (>20 = likely pathogenic)
- SIFT predictions (D = deleterious)
- PolyPhen2 predictions (D = damaging)

## Files Generated

### Core Data Files
- `enhanced_test_data.vcf.gz` - VCF with realistic annotations
- `enhanced_test_data.vcf.gz.tbi` - Tabix index

### Test Infrastructure
- `run_enhanced_tests.sh` - Test runner with real annotation validation
- `enhanced_dataset_statistics.json` - Detailed statistics including annotation metrics

[Rest of files same as comprehensive dataset...]

## Expected Results

When running gene burden analysis on this enhanced dataset:

- **Disease genes** should show **significant enrichment** in cases
- **High-impact variants** in disease genes should have stronger association
- **Annotation-based filtering** should work with realistic dbNSFP scores
- **Effect prioritization** should properly rank variant consequences

## Validation

The enhanced dataset provides:
- Realistic annotation diversity for testing annotation-based filters
- Authentic pathogenicity scores for testing scoring algorithms
- Controlled genotype distributions ensuring reliable statistical results
- Comprehensive test coverage across all gene burden analysis methods

Check `enhanced_dataset_statistics.json` for detailed metrics on annotation sampling and effect distributions.
