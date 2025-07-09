# GIAB Test Data Infrastructure

Complete pipeline for downloading, processing, and analyzing GIAB (Genome in a Bottle) test data for genomic variant analysis.

## Overview

This directory contains three scripts that create a complete genomic testing workflow:

1. **`download_giab_bams.py`** - Extract GIAB HG19 exome data from remote sources
2. **`download_reference_genome.py`** - Download human_g1k_v37_decoy reference genome
3. **`call_variants_gatk.py`** - Call variants using GATK HaplotypeCaller

## Quick Start

```bash
# 1. Download GIAB BAM files (HG002, HG003, HG004 trio)
python download_giab_bams.py

# 2. Download reference genome
python download_reference_genome.py

# 3. Call variants
python call_variants_gatk.py
```

## Prerequisites

### Required Tools
- **Python 3.6+** - Script execution
- **samtools** - BAM manipulation and indexing
- **wget** - File downloads from FTP servers
- **GATK 4.x** - Variant calling
- **Internet connection** - Access to remote data sources

### Installation
```bash
# Conda/Mamba installation
mamba install -c bioconda samtools gatk4 wget

# Or individual installations
# samtools: https://github.com/samtools/samtools
# GATK: https://github.com/broadinstitute/gatk
# wget: usually pre-installed on Linux/macOS
```

## Workflow Details

### Step 1: GIAB BAM Download (`download_giab_bams.py`)

Downloads and extracts exome regions from GIAB trio samples using remote streaming.

**Data Source:** Oslo University Hospital HG19 exome BAMs
- **HG002** - Ashkenazi son (proband)
- **HG003** - Ashkenazi father
- **HG004** - Ashkenazi mother

**Gene Regions (8 genes):**
- **BRCA1** (17:41,194,312-41,279,381) - Breast cancer susceptibility
- **TP53** (17:7,571,720-7,590,868) - Tumor suppressor
- **GRIN2A** (16:9,849,257-10,265,950) - Glutamate receptor
- **NAA10** (X:153,198,978-153,205,254) - N-terminal acetyltransferase
- **PKD2** (4:88,999,644-89,070,785) - Polycystic kidney disease
- **PKD1** (16:2,136,709-2,187,899) - Polycystic kidney disease (primary)
- **TTN** (2:179,388,716-179,674,150) - Titin (largest human gene)
- **COL4A5** (X:107,681,068-107,942,775) - Collagen (X-linked Alport syndrome)

**Usage:**
```bash
# Basic download
python download_giab_bams.py

# Skip MD5 verification (faster)
python download_giab_bams.py --skip-md5

# Custom output directory
python download_giab_bams.py --output-dir /path/to/data
```

**Features:**
- **Remote streaming** - Extracts regions directly without downloading full BAMs
- **Automatic indexing** - Creates .bai files for extracted BAMs
- **Smart caching** - Skips existing files
- **Progress tracking** - Real-time download and extraction status
- **Clean extraction** - Removes temporary files after processing

**Output:**
```
data/
├── HG002.HG19.exome.bam     # ~50-100 MB per sample
├── HG002.HG19.exome.bam.bai
├── HG003.HG19.exome.bam
├── HG003.HG19.exome.bam.bai
├── HG004.HG19.exome.bam
└── HG004.HG19.exome.bam.bai
```

### Step 2: Reference Genome Download (`download_reference_genome.py`)

Downloads the standard GRCh37/HG19 reference genome with decoy sequences.

**Data Source:** Broad Institute FTP server
- **human_g1k_v37_decoy.fasta** (~3GB) - Reference sequence
- **human_g1k_v37_decoy.fasta.fai** (~3KB) - FASTA index
- **human_g1k_v37_decoy.dict** (~4KB) - Sequence dictionary

**Usage:**
```bash
# Download and decompress all files
python download_reference_genome.py

# Keep compressed files after decompression
python download_reference_genome.py --keep-compressed

# Download compressed files only
python download_reference_genome.py --compressed-only

# Custom output directory
python download_reference_genome.py --output-dir /path/to/reference
```

**Features:**
- **Resume support** - Continues interrupted downloads
- **Integrity checks** - Validates file sizes and content
- **Automatic decompression** - Extracts .gz files
- **Progress tracking** - Real-time download progress

**Output:**
```
reference/
├── human_g1k_v37_decoy.fasta      # ~3GB decompressed reference
├── human_g1k_v37_decoy.fasta.fai  # ~3KB index
├── human_g1k_v37_decoy.dict       # ~4KB dictionary
├── human_g1k_v37_decoy.fasta.gz   # ~800MB compressed (optional)
├── human_g1k_v37_decoy.fasta.fai.gz # ~1KB compressed (optional)
└── human_g1k_v37_decoy.dict.gz    # ~1KB compressed (optional)
```

### Step 3: Variant Calling (`call_variants_gatk.py`)

Performs multi-sample variant calling using GATK HaplotypeCaller.

**Approach:** Single GATK command with multiple input BAMs for optimal family-based calling.

**Usage:**
```bash
# Basic variant calling
python call_variants_gatk.py

# More memory for better performance
python call_variants_gatk.py --java-memory 8g

# Custom directories
python call_variants_gatk.py --data-dir data --reference-dir reference --output-dir variants
```

**GATK Command:**
```bash
gatk --java-options -Xmx4g HaplotypeCaller \
  -R reference/human_g1k_v37_decoy.fasta \
  -I data/HG002.HG19.exome.bam \
  -I data/HG003.HG19.exome.bam \
  -I data/HG004.HG19.exome.bam \
  -L variants/gene_regions.interval_list \
  -O variants/HG002-trio_haplotypecaller.vcf.gz
```

**Process:**
1. **Validation** - Checks GATK and samtools installation
2. **File discovery** - Finds BAM files and reference genome
3. **BED sorting** - Sorts gene regions by chromosome and position
4. **Interval conversion** - Creates GATK interval list from BED file
5. **Multi-sample calling** - Runs GATK with all samples together

**Output:**
```
variants/
├── gene_regions_sorted.bed           # Sorted BED file
├── gene_regions.interval_list        # GATK interval list
└── HG002-trio_haplotypecaller.vcf.gz # Multi-sample VCF (~2-10 MB)
```

## Complete File Structure

After running all three scripts:

```
tests/fixtures/giab/
├── data/                             # GIAB BAM files
│   ├── HG002.HG19.exome.bam
│   ├── HG002.HG19.exome.bam.bai
│   ├── HG003.HG19.exome.bam
│   ├── HG003.HG19.exome.bam.bai
│   ├── HG004.HG19.exome.bam
│   └── HG004.HG19.exome.bam.bai
├── reference/                        # Reference genome files
│   ├── human_g1k_v37_decoy.fasta
│   ├── human_g1k_v37_decoy.fasta.fai
│   └── human_g1k_v37_decoy.dict
├── variants/                         # Variant calling outputs
│   ├── gene_regions_sorted.bed
│   ├── gene_regions.interval_list
│   └── HG002-trio_haplotypecaller.vcf.gz
├── logs/                            # Execution logs
│   ├── giab_download_*.log
│   ├── reference_download_*.log
│   └── variant_calling_*.log
├── download_giab_bams.py            # Script 1: BAM download
├── download_reference_genome.py     # Script 2: Reference download
├── call_variants_gatk.py            # Script 3: Variant calling
├── gene_regions_GRCh37.bed         # Gene coordinates (HG19)
└── README.md                        # This file
```

## Performance Expectations

### Download Times
- **GIAB BAMs:** 5-15 minutes (depending on network)
- **Reference genome:** 10-30 minutes (depending on network)

### Processing Times
- **Variant calling:** 10-30 minutes (trio on 8 genes)

### Disk Usage
- **Total:** ~3.5GB
- **GIAB BAMs:** ~300 MB (extracted regions only)
- **Reference:** ~3.2 GB (decompressed)
- **Variants:** ~10 MB

### Memory Usage
- **BAM download:** <1GB RAM
- **Reference download:** <1GB RAM
- **Variant calling:** 4-8GB RAM (configurable)

## Integration with VariantCentrifuge

The output VCF is ready for analysis with VariantCentrifuge:

```bash
# Analyze variants with inheritance patterns
variantcentrifuge \
  --input tests/fixtures/giab/variants/HG002-trio_haplotypecaller.vcf.gz \
  --output giab_analysis/ \
  --ped tests/fixtures/giab/trio.ped
```

## Troubleshooting

### Common Issues

**1. Tool Not Found**
```bash
# Check tool availability
gatk --version
samtools --version
wget --version

# Install missing tools
mamba install -c bioconda gatk4 samtools wget
```

**2. Download Failures**
```bash
# Check network connectivity
ping ftp.broadinstitute.org
ping depot.galaxyproject.org

# Resume interrupted downloads (automatic)
python download_giab_bams.py
python download_reference_genome.py
```

**3. Memory Issues**
```bash
# Increase Java memory for GATK
python call_variants_gatk.py --java-memory 8g

# Check available RAM
free -h
```

**4. Disk Space**
```bash
# Check available space (~4GB needed)
df -h .

# Clean up if needed
rm -rf data/temp/  # Remove temporary files
```

**5. GATK Errors**
```bash
# Check GATK installation
gatk --list

# Verify reference files exist
ls -la reference/human_g1k_v37_decoy.*

# Re-download reference if corrupted
rm -rf reference/
python download_reference_genome.py
```

## Advanced Usage

### Custom Gene Lists
Edit `gene_regions_GRCh37.bed` to target different genes:
```bash
# Format: chr    start    end    gene_name
1    43042295    43127364    BRCA1
17   7571720     7590868     TP53
```

### Parallel Processing
For faster BAM extraction with multiple regions:
```bash
# The script automatically optimizes for parallel extraction
python download_giab_bams.py  # Already optimized
```

### Quality Control
```bash
# Check BAM quality
samtools flagstat data/HG002.HG19.exome.bam

# Validate VCF
gatk ValidateVariants -V variants/HG002-trio_haplotypecaller.vcf.gz -R reference/human_g1k_v37_decoy.fasta

# Basic variant stats
bcftools stats variants/HG002-trio_haplotypecaller.vcf.gz
```

## Data Sources and Citations

### GIAB Consortium
- **Website:** https://www.nist.gov/programs-projects/genome-bottle
- **Paper:** Zook, J.M., et al. (2014). "Integrating human sequence data sets provides a resource of benchmark SNP and indel genotypes." Nature Biotechnology.

### Reference Genome
- **Source:** Broad Institute GATK Resource Bundle
- **Build:** GRCh37/hg19 with 1000 Genomes Project decoy sequences
- **Citation:** Genome Reference Consortium. "Genome Reference Consortium Human Build 37 (GRCh37)."

### Oslo University Hospital
- **Data hosting:** HG19 exome BAM files
- **Access:** Public FTP server for research use

## License and Usage

These scripts are designed for research and educational purposes. The GIAB data is publicly available for research use. Please cite appropriate sources when using this data in publications.

## Support

For issues with these scripts:
1. Check the troubleshooting section above
2. Review log files in the `logs/` directory
3. Verify tool installations and versions
4. Ensure adequate disk space and memory

The scripts include comprehensive logging and error handling to help diagnose issues.