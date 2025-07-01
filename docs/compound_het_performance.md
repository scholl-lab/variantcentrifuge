# Compound Heterozygous Analysis - Performance Optimization

## Overview

VariantCentrifuge includes an optimized vectorized implementation of compound heterozygous analysis that provides significant performance improvements for large-scale genomic analyses.

## Features

### Performance Improvements
- **10-50x faster** for genes with many variants
- **60x memory reduction** using numeric genotype encoding
- **Efficient handling** of genes with 3+ heterozygous variants

### When to Use

The vectorized implementation is enabled by default and is recommended for:
- Large multi-sample VCF files
- Genes with many variants (>10 variants)
- Production pipelines requiring optimal performance

## Usage

### Command Line

The vectorized implementation is enabled by default:
```bash
variantcentrifuge --vcf input.vcf --ped family.ped --calculate-inheritance
```

To use the original implementation (for debugging/comparison):
```bash
variantcentrifuge --vcf input.vcf --ped family.ped --calculate-inheritance --no-vectorized-comp-het
```

### Performance Comparison

| Gene Variants | Original Time | Vectorized Time | Speedup |
|---------------|---------------|-----------------|---------|
| 10            | 0.01s         | 0.001s          | 10x     |
| 50            | 0.19s         | 0.004s          | 48x     |
| 100           | 1.23s         | 0.017s          | 72x     |
| 500           | 31.0s         | 0.55s           | 56x     |

## Multiple Variant Handling

When a sample has 3+ heterozygous variants in the same gene:

### Example: 5 Heterozygous Variants
- Total possible pairs: 10
- If father carries 3 variants and mother carries 2 variants
- Trans-configured pairs: 3 × 2 = 6 true compound hets

### Prioritization

When multiple compound het pairs exist, they are prioritized by:

1. **Variant Impact**: HIGH > MODERATE > LOW
2. **Allele Frequency**: Rarer variants prioritized
3. **Functional Scores**: CADD, REVEL, SpliceAI
4. **Inheritance Confidence**: Clear trans > ambiguous

## Technical Details

### Genotype Encoding
```
./. → -1 (missing)
0/0 → 0  (reference)
0/1 → 1  (heterozygous)
1/1 → 2  (homozygous alt)
```

### Memory Usage
- Original: ~20-30 bytes per genotype (string)
- Vectorized: 1 byte per genotype (int8)

### Algorithm Complexity
- Original: O(n²) with high constant factor
- Vectorized: O(n²) with low constant factor + vectorization

## Implementation Notes

The vectorized implementation:
- Uses NumPy arrays for efficient computation
- Processes all samples in batch
- Caches encoded genotypes
- Eliminates repeated DataFrame operations
- Maintains full compatibility with original output