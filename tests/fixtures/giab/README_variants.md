# GATK Variant Calling

Script to call variants using GATK HaplotypeCaller on GIAB BAM files with the human_g1k_v37_decoy reference genome.

## Overview

This script performs multi-sample variant calling on the extracted GIAB HG19 exome BAM files using GATK HaplotypeCaller:

1. **Multi-sample variant calling** - Calls variants on all samples together in a single command
2. **Targeted regions** - Only calls variants in the gene regions from the BED file
3. **Simple approach** - Direct VCF output without GVCF intermediate files

## Prerequisites

Before running variant calling, ensure you have:

```bash
# 1. GIAB BAM files
python download_giab_bams.py

# 2. Reference genome
python download_reference_genome.py

# 3. GATK installed and in PATH
gatk --version
```

## Usage

```bash
# Basic variant calling (all samples)
python call_variants_gatk.py

# Custom directories
python call_variants_gatk.py --data-dir data --reference-dir reference --output-dir variants

# More memory for large samples
python call_variants_gatk.py --java-memory 8g
```

## What it does

### Step 1: Preparation
- Validates GATK and samtools installation
- Finds BAM files in data directory
- Checks reference genome files
- Creates FASTA index (`.fai`) if missing
- Creates sequence dictionary (`.dict`) if missing
- Sorts BED file by chromosome and position
- Converts BED file to GATK interval list format

### Step 2: Multi-sample variant calling
- Runs GATK HaplotypeCaller with all samples in a single command
- Uses multiple `-I` inputs for each BAM file
- Directly outputs a multi-sample VCF file
- Processes all samples together for optimal family-based calling

## Output Files

```
variants/
├── gene_regions_sorted.bed           # Sorted BED file
├── gene_regions.interval_list        # GATK interval list
└── HG002-trio_haplotypecaller.vcf.gz # Multi-sample VCF
```

## GATK Parameters

The script uses a simple GATK HaplotypeCaller command:

```bash
gatk --java-options -Xmx4g HaplotypeCaller \
  -R reference/human_g1k_v37_decoy.fasta \
  -I data/HG002.HG19.exome.bam \
  -I data/HG003.HG19.exome.bam \
  -I data/HG004.HG19.exome.bam \
  -L variants/gene_regions.interval_list \
  -O variants/HG002-trio_haplotypecaller.vcf.gz
```

## Performance

**Typical runtime for trio:**
- Small exome regions (~8 genes): 10-30 minutes total
- Memory usage: 4-8GB RAM
- Output size: 2-10 MB VCF

**Optimization tips:**
- Use `--java-memory 8g` for better performance
- Use SSD storage for faster I/O
- Multi-sample calling is more efficient than individual samples

## Gene Regions Covered

The script calls variants in these 8 genes:
- **BRCA1** - Breast cancer susceptibility
- **TP53** - Tumor suppressor gene
- **GRIN2A** - Glutamate receptor
- **NAA10** - N-terminal acetyltransferase
- **PKD2** - Polycystic kidney disease
- **PKD1** - Polycystic kidney disease (primary)
- **TTN** - Titin (largest human gene)
- **COL4A5** - Collagen (X-linked Alport syndrome)

## Integration with VariantCentrifuge

The output VCF file is ready for analysis with VariantCentrifuge:

```bash
# Use multi-sample VCF
variantcentrifuge --input variants/HG002-trio_haplotypecaller.vcf.gz --output analysis/
```

## Troubleshooting

### GATK Installation Issues
```bash
# Check GATK is in PATH
gatk --version

# Common install methods
conda install -c bioconda gatk4
# or
wget https://github.com/broadinstitute/gatk/releases/download/4.x.x/gatk-4.x.x.zip
```

### Memory Issues
```bash
# Increase Java memory
python call_variants_gatk.py --java-memory 8g

# Or reduce parallel threads
python call_variants_gatk.py --threads 1
```

### Missing Reference Files
```bash
# Download reference genome
python download_reference_genome.py

# Check files exist
ls reference/human_g1k_v37_decoy.fasta*
```

### BAM File Issues
```bash
# Check BAM files exist and are indexed
ls data/*.bam data/*.bai

# Re-run BAM extraction if needed
python download_giab_bams.py
```

## Expected Results

**Variant counts (approximate):**
- Total variants: 150-600 across all samples in gene regions
- SNPs: 120-450 across all samples
- Indels: 30-150 across all samples
- High-quality variants: 90-300 across all samples

**File sizes:**
- Multi-sample VCF: 2-10 MB
- Processing time: 10-30 minutes total

## Quality Metrics

The script uses GATK HaplotypeCaller default quality settings for reliable variant calling.

For additional quality control, consider running GATK's Variant Quality Score Recalibration (VQSR) on larger datasets.

## Next Steps

After variant calling, you can:
1. **Analyze with VariantCentrifuge** - Run the full analysis pipeline
2. **Annotate variants** - Add functional annotations
3. **Filter variants** - Apply quality and frequency filters
4. **Inheritance analysis** - Study family-based variant patterns
5. **Benchmark results** - Compare against GIAB truth sets