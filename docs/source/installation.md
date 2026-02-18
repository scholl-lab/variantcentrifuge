# Installation

## Dependencies

### External Tools

VariantCentrifuge requires several bioinformatics tools to be installed and available in your `PATH`:

- **bcftools** - For variant extraction and manipulation
- **snpEff** - For generating gene BED files and functional annotations
- **SnpSift** - For filtering and field extraction
- **bedtools** (specifically `sortBed`) - For sorting BED files

#### Installing via mamba/conda

```bash
mamba create -y -n annotation bcftools snpsift snpeff bedtools
mamba activate annotation
```

Ensure these tools are in your `PATH` before running VariantCentrifuge.

### Python Requirements

VariantCentrifuge requires **Python 3.10+** and the following Python packages:

- `pandas` - For XLSX conversion and data handling
- `jinja2` - For HTML template rendering
- `openpyxl` - For XLSX creation
- `scipy` - For Fisher exact test in variant analysis
- `statsmodels` - For multiple testing correction in gene burden analysis

## Installation Methods

### Method 1: Install from PyPI (Recommended)

```bash
pip install variantcentrifuge
```

### Method 2: Install from Source

1. **Clone the repository:**
   ```bash
   git clone https://github.com/scholl-lab/variantcentrifuge.git
   cd variantcentrifuge
   ```

2. **Set up a virtual environment (recommended):**
   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install the package:**
   ```bash
   pip install .
   ```

   For development (editable install):
   ```bash
   pip install -e .
   ```

### Method 3: Using conda environment

```bash
# Clone the repository
git clone https://github.com/scholl-lab/variantcentrifuge.git
cd variantcentrifuge

# Create and activate conda environment
mamba env create -f conda/environment.yml
mamba activate annotation

# Install the package
pip install -e .
```

### Method 4: Docker (Recommended for Quick Setup)

The Docker image bundles all external tools (bcftools, snpEff, SnpSift, bedtools) and Python dependencies into a single container. No local installation of bioinformatics tools is required.

```bash
# Pull the latest image from GitHub Container Registry
docker pull ghcr.io/scholl-lab/variantcentrifuge:latest

# Verify the installation
docker run --rm ghcr.io/scholl-lab/variantcentrifuge:latest --version
```

Place your VCF files in a local directory and mount it into the container:

```bash
docker run --rm -v ./data:/data \
  ghcr.io/scholl-lab/variantcentrifuge:latest \
  --gene-name BRCA1 \
  --vcf-file /data/input.vcf.gz \
  --output-file /data/output.tsv
```

A `docker-compose.yml` is included in the repository for convenience:

```yaml
services:
  variantcentrifuge:
    image: ghcr.io/scholl-lab/variantcentrifuge:latest
    volumes:
      - ./data:/data
      # Mount snpEff databases (download once, reuse)
      # - /path/to/snpeff_data:/snpeff_data:ro
      # Override built-in scoring configs
      # - ./my_scoring:/app/scoring:ro
```

```bash
docker compose run --rm variantcentrifuge \
  --gene-name BRCA1 --vcf-file /data/input.vcf.gz --output-file /data/output.tsv
```

The image runs as a non-root user, includes built-in scoring models at `/app/scoring/`, and is signed with cosign for supply chain security.

## Verification

Verify your installation by running:

```bash
variantcentrifuge --version
```

And check that external tools are available:

```bash
bcftools --version
snpEff -version
java -jar $SNPSIFT_JAR
bedtools --version
```

## Troubleshooting

### Common Issues

1. **External tools not found:**
   - Ensure all external tools are installed and in your PATH
   - For conda installations, activate the correct environment

2. **Permission errors:**
   - Use `--user` flag with pip: `pip install --user variantcentrifuge`
   - Or use a virtual environment

3. **Version conflicts:**
   - Use a fresh virtual environment or conda environment
   - Update pip: `pip install --upgrade pip`

### Environment Variables

You may need to set environment variables for some tools:

```bash
export SNPEFF_JAR=/path/to/snpEff.jar
export SNPSIFT_JAR=/path/to/SnpSift.jar
```

## Next Steps

Once installed, proceed to the [Usage Guide](usage.md) to learn how to configure and run VariantCentrifuge.
