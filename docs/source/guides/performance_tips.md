# Performance Optimization

Tips for optimizing VariantCentrifuge performance with large datasets.

## Hardware Recommendations

- **Memory:** 16GB+ RAM for large VCF files
- **Storage:** SSD for intermediate files
- **CPU:** Multi-core for parallel processing

## Performance Features

### bcftools Pre-filtering

The `--bcftools-prefilter` option allows you to apply fast bcftools filters during the initial variant extraction step, significantly reducing the amount of data processed by subsequent steps:

```bash
# Filter out low-quality and common variants early
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file large_cohort.vcf.gz \
  --bcftools-prefilter 'FILTER="PASS" && QUAL>20 && INFO/AC<10' \
  --preset rare,coding \
  --output-file results.tsv
```

**Benefits:**
- Reduces memory usage by filtering early in the pipeline
- Speeds up SnpSift filtering and field extraction
- Particularly effective for large cohort VCFs

**bcftools Filter Syntax:**
- `FILTER="PASS"` - Only PASS variants
- `QUAL>20` - Quality score threshold
- `INFO/AC<10` - Allele count less than 10
- `INFO/AF<0.01` - Allele frequency less than 1%
- Combine with `&&` (AND) or `||` (OR)

### Optimized Filtering Workflow

For maximum performance with large datasets:

```bash
# Combine bcftools pre-filter with late filtering
variantcentrifuge \
  --gene-file large_gene_list.txt \
  --vcf-file cohort_1000_samples.vcf.gz \
  --bcftools-prefilter 'FILTER="PASS" && INFO/AC<20' \
  --late-filtering \
  --preset rare,coding \
  --scoring-config-path scoring/my_scores \
  --filters "score > 0.5" \
  --output-file high_score_rare.tsv
```

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

## Storage Optimization

### Archive Results

For easy storage and transfer of analysis results, use the `--archive-results` option to create a compressed archive:

```bash
# Create timestamped archive after analysis
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --html-report \
  --xlsx \
  --archive-results \
  --output-file results.tsv

# Result: variantcentrifuge_results_input_20241217_143052.tar.gz
```

**Benefits:**
- Automatic compression reduces storage requirements by 50-80%
- Timestamped archives for version tracking
- Single file for easy transfer and backup
- Archive placed outside output directory to avoid recursion

### Intermediate File Management

```bash
# Delete intermediate files automatically (default behavior)
variantcentrifuge ... # intermediates deleted after success

# Keep intermediates for debugging
variantcentrifuge --keep-intermediates ...

# Combine with archiving for complete workflow
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --archive-results \  # Archive everything
  --output-file results.tsv  # Then clean intermediates
```