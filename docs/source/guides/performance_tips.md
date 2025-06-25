# Performance Optimization

Tips for optimizing VariantCentrifuge performance with large datasets.

## Hardware Recommendations

- **Memory:** 16GB+ RAM for large VCF files
- **Storage:** SSD for intermediate files
- **CPU:** Multi-core for parallel processing

## Large Dataset Strategies

### Chromosome-based Processing

```bash
# Split by chromosome
for chr in {1..22} X Y; do
    bcftools view -r chr${chr} large_cohort.vcf.gz | \
    variantcentrifuge \
        --gene-file genes.txt \
        --vcf-file /dev/stdin \
        --output-file chr${chr}_results.tsv &
done
wait

# Merge results
cat chr*_results.tsv > combined_results.tsv
```

### Memory Management

```bash
# Increase Java heap size
export JAVA_OPTS="-Xmx16g"

# Use streaming for large files
bcftools view large.vcf.gz | variantcentrifuge --vcf-file /dev/stdin ...
```