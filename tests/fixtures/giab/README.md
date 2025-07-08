# GIAB Test Data

Simple script to download GIAB BAM files and extract specific gene regions.

## Usage

```bash
# Download all samples for GRCh38
python download_giab_bams.py --assembly GRCh38

# Download all samples for GRCh37
python download_giab_bams.py --assembly GRCh37

# Download specific sample
python download_giab_bams.py --assembly GRCh38 --samples HG002

# Keep full BAM files after extraction
python download_giab_bams.py --assembly GRCh38 --keep-full

# Use more threads for faster processing
python download_giab_bams.py --assembly GRCh38 --threads 8
```

## What it does

- Downloads full GIAB BAM files (~100GB each) in parallel
- Extracts only gene regions from `gene_regions_GRCh38.bed` or `gene_regions_GRCh37.bed`
- Creates small BAM files (~2-5 MB each) in the `data/` directory
- Removes full BAM files after extraction (unless --keep-full is used)

## Gene regions

The BED files contain these genes:
- BRCA1
- TP53  
- GRIN2A
- NAA10
- PKD2

## Output

```
data/
├── HG002.GRCh38.bam
├── HG002.GRCh38.bam.bai
├── HG003.GRCh38.bam
├── HG003.GRCh38.bam.bai
├── HG004.GRCh38.bam
└── HG004.GRCh38.bam.bai
```

## Features

- **Resume support**: Uses `wget -c` to resume interrupted downloads
- **MD5 verification**: Verifies downloaded files against known checksums
- **Smart caching**: Skips files that already exist with correct MD5
- **Detailed logging**: Saves logs to `logs/` directory with timestamps
- **Progress tracking**: Saves download status to `download_status.json`
- **Parallel downloads**: Downloads multiple samples simultaneously
- **Integrity checks**: Verifies BAM files with `samtools quickcheck`

## Files Created

- `data/` - Extracted BAM files (~2-5 MB each)
- `logs/` - Download logs with timestamps
- `download_status.json` - Track download progress and MD5 checksums
- `*.extraction_info.json` - Metadata for each extracted BAM

## Requirements

- samtools
- wget