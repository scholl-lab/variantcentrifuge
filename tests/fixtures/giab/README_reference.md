# Reference Genome Download

Script to download the human_g1k_v37_decoy reference genome from the Broad Institute FTP server.

## Overview

This downloads the standard reference genome used in variant calling pipelines:
- **human_g1k_v37_decoy.fasta** - Reference genome sequence (~3GB)
- **human_g1k_v37_decoy.fasta.fai** - FASTA index file (~3KB)

The reference includes decoy sequences to improve alignment accuracy and reduce false positive variant calls.

## Usage

```bash
# Download and decompress reference files
python download_reference_genome.py

# Download to specific directory
python download_reference_genome.py --output-dir /path/to/reference

# Keep compressed files after decompression
python download_reference_genome.py --keep-compressed

# Download compressed files only (don't decompress)
python download_reference_genome.py --compressed-only
```

## What it does

- Downloads from Broad Institute FTP server: `ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/`
- Downloads compressed files: `human_g1k_v37_decoy.fasta.gz` and `human_g1k_v37_decoy.fasta.fai.gz`
- Automatically decompresses files (unless `--compressed-only` specified)
- Performs integrity checks on downloaded files
- Shows download progress with progress bars

## Output Files

```
reference/
├── human_g1k_v37_decoy.fasta      # ~3GB decompressed reference
├── human_g1k_v37_decoy.fasta.fai  # ~3KB index file
├── human_g1k_v37_decoy.fasta.gz   # ~800MB compressed (optional)
└── human_g1k_v37_decoy.fasta.fai.gz # ~1KB compressed (optional)
```

## Features

- **Resume support**: Uses `wget -c` to resume interrupted downloads
- **Progress tracking**: Real-time download progress bars
- **Smart caching**: Skips downloads if files already exist
- **Integrity checks**: Verifies file sizes and content
- **Flexible compression**: Option to keep or remove compressed files
- **Detailed logging**: Saves logs with timestamps

## Reference Genome Details

- **Build**: GRCh37/hg19 with 1000 Genomes Project decoy sequences
- **Chromosomes**: 25 main chromosomes (1-22, X, Y, MT) plus decoy contigs
- **Size**: ~3.2 billion base pairs
- **Format**: FASTA with consistent chromosome naming (1, 2, 3... not chr1, chr2, chr3)
- **Use cases**: Variant calling, alignment, reference-based analysis

## Requirements

- **wget** - For downloading files from FTP server
- **Python 3.6+** - For script execution
- **~4GB disk space** - For decompressed files
- **Internet connection** - Access to Broad Institute FTP server

## Troubleshooting

### FTP Access Issues
If you experience FTP connection problems:
- Check firewall settings
- Try from a different network
- Consider using the Google Cloud bucket alternative (if available)

### Download Failures
- Use `--keep-compressed` to preserve partial downloads
- Check available disk space (~4GB needed)
- Verify network stability for large file downloads

### File Corruption
The script performs basic integrity checks:
- File size validation (~3GB for FASTA)
- Index file sequence count verification
- Non-empty file checks

## Integration with Variant Calling

Use the downloaded reference with common tools:

```bash
# GATK
gatk HaplotypeCaller -R reference/human_g1k_v37_decoy.fasta -I sample.bam

# BWA alignment
bwa index reference/human_g1k_v37_decoy.fasta
bwa mem reference/human_g1k_v37_decoy.fasta reads.fastq

# Samtools
samtools faidx reference/human_g1k_v37_decoy.fasta chr1:1000-2000
```

## Important Notes

⚠️ **FTP Deprecation**: The Broad Institute FTP server was scheduled for deprecation. This script provides access to the historical b37 reference, but for production use, consider:
- Google Cloud bucket versions
- NCBI reference genomes  
- Terra workspace preconfigured references

The human_g1k_v37_decoy reference remains a standard in genomics pipelines and is widely used for reproducibility with existing analyses.