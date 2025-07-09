# GIAB Test Data

Script to extract specific gene regions directly from remote GIAB HG19 exome BAM files using samtools streaming.

## Usage

```bash
# Extract all available samples (HG002, HG003, HG004)
python download_giab_bams.py

# Extract specific sample
python download_giab_bams.py --samples HG002

# Use more threads for faster processing
python download_giab_bams.py --threads 8
```

## What it does

- Streams data directly from remote GIAB HG19 exome BAM files (no local download needed)
- Extracts only gene regions from `gene_regions_GRCh37.bed`
- Creates small BAM files (~1-2 MB each) in the `data/` directory
- Uses samtools with individual region specifications for remote file compatibility
- Automatically sorts and indexes extracted BAM files

## Gene regions

The BED files contain these genes:
- **BRCA1** - Breast cancer gene
- **TP53** - Tumor suppressor gene  
- **GRIN2A** - Glutamate receptor gene
- **NAA10** - N-terminal acetyltransferase
- **PKD2** - Polycystic kidney disease gene
- **PKD1** - Polycystic kidney disease gene (primary)
- **TTN** - Titin gene (largest human gene)
- **COL4A5** - Collagen gene (X-linked Alport syndrome)

## Output

```
data/
├── HG002.HG19.exome.bam
├── HG002.HG19.exome.bam.bai
├── HG003.HG19.exome.bam
├── HG003.HG19.exome.bam.bai
├── HG004.HG19.exome.bam
└── HG004.HG19.exome.bam.bai
```

## Features

- **No downloads**: Streams data directly from remote FTP servers
- **Smart caching**: Skips extraction if valid BAM already exists
- **Detailed logging**: Saves logs to `logs/` directory with timestamps
- **Parallel processing**: Extracts multiple samples simultaneously
- **Integrity checks**: Verifies BAM files with `samtools quickcheck`
- **Performance tracking**: Records extraction time and file sizes

## Files Created

- `data/` - Extracted BAM files (~1-2 MB each from exome data)
- `logs/` - Extraction logs with timestamps
- `*.extraction_info.json` - Metadata for each extracted BAM including timing

## Available Samples

- **HG002**: Son from Ashkenazi trio
- **HG003**: Father from Ashkenazi trio  
- **HG004**: Mother from Ashkenazi trio

## Requirements

- samtools (only dependency needed)