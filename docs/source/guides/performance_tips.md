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

## Checkpoint and Resume

VariantCentrifuge includes a checkpoint system that tracks pipeline progress and allows resumption after interruptions. This is particularly useful for long-running analyses on large datasets.

### Basic Checkpoint Usage

Enable checkpointing with the `--enable-checkpoint` flag:

```bash
# Run with checkpoint tracking
variantcentrifuge \
  --gene-file large_gene_list.txt \
  --vcf-file large_cohort.vcf.gz \
  --enable-checkpoint \
  --threads 8 \
  --output-file results.tsv
```

If the pipeline is interrupted (e.g., system crash, job timeout), resume from the last checkpoint:

```bash
# Resume from checkpoint
variantcentrifuge \
  --gene-file large_gene_list.txt \
  --vcf-file large_cohort.vcf.gz \
  --enable-checkpoint \
  --resume \
  --threads 8 \
  --output-file results.tsv
```

### How Checkpointing Works

The checkpoint system:
- Tracks each major pipeline step (gene BED creation, variant extraction, filtering, etc.)
- Saves state to `.variantcentrifuge_state.json` in the output directory
- Records input/output files, parameters, and completion status for each step
- Validates configuration hasn't changed between runs
- Skips already-completed steps when resuming

### Checkpoint Features

#### Parallel Processing Support
Checkpoints work seamlessly with parallel processing (`--threads`):
- Tracks individual chunk processing in parallel runs
- Properly handles the merge step after parallel processing
- Maintains thread-safe state updates

#### File Integrity Checking
Optional file checksum validation ensures data integrity:

```bash
# Enable checksum validation (slower but more reliable)
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --enable-checkpoint \
  --checkpoint-checksum \
  --output-file results.tsv
```

#### Status Inspection
Check the status of a previous run without resuming:

```bash
# Show checkpoint status
variantcentrifuge \
  --show-checkpoint-status \
  --output-dir previous_analysis/
```

### Tracked Pipeline Steps

The checkpoint system tracks these major steps:
1. **Gene BED creation** - Converting gene names to genomic regions
2. **Parallel processing** (if using `--threads`):
   - BED file splitting
   - Individual chunk processing
   - Chunk merging
3. **TSV sorting** - Sorting extracted variants by gene
4. **Genotype replacement** - Converting genotypes to sample IDs
5. **Phenotype integration** - Adding phenotype data
6. **Variant analysis** - Computing statistics and scores
7. **Final output** - Writing results and reports

### Best Practices

1. **Use for long-running analyses**: Particularly beneficial for:
   - Large cohort VCFs (>100 samples)
   - Extensive gene lists (>100 genes)
   - Complex scoring and annotation pipelines

2. **Combine with job schedulers**: Ideal for HPC environments:
   ```bash
   #!/bin/bash
   #SBATCH --time=24:00:00
   #SBATCH --mem=32G
   
   variantcentrifuge \
     --gene-file all_genes.txt \
     --vcf-file cohort.vcf.gz \
     --enable-checkpoint \
     --resume \  # Safe to use even on first run
     --threads 16 \
     --output-file results.tsv
   ```

3. **Monitor progress**: The log shows which steps are skipped:
   ```
   INFO: Skipping gene BED creation (already completed)
   INFO: Skipping chunk merging (already completed)
   INFO: Resuming from genotype replacement...
   ```

### Limitations

- Configuration must remain identical between runs (same filters, fields, etc.)
- Pipeline version must match (no resume across VariantCentrifuge updates)
- Intermediate files must not be manually modified
- Output directory structure must remain intact

### Troubleshooting Checkpoint Issues

1. **"Configuration has changed, cannot resume"**
   - Ensure all command-line arguments match the original run
   - Check that config files haven't been modified

2. **"Output file validation failed"**
   - An intermediate file may have been corrupted or deleted
   - Remove `.variantcentrifuge_state.json` to start fresh
   - Use `--checkpoint-checksum` for better validation

3. **"Pipeline version mismatch"**
   - VariantCentrifuge was updated between runs
   - Complete the analysis with the original version or start fresh