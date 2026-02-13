# Configuration

This directory holds the user-editable configuration files for the
variantcentrifuge Snakemake workflow.

## Files

| File | Purpose |
|------|---------|
| `config_vc.yaml` | Main workflow configuration (paths, parameters, IGV, container) |
| `samples.tsv` | Tab-separated sample sheet listing samples and their VCF basenames |

## Quick start

1. Copy the template files and rename them (or edit in place):

   ```bash
   cp config/config_vc.yaml config/my_config.yaml
   cp config/samples.tsv    config/my_samples.tsv
   ```

2. Replace every `EDIT_ME` placeholder with your actual paths.

3. Populate `samples.tsv` with one row per sample. The `vcf_basename` column
   should contain the VCF filename **without** the `.vcf.gz` extension. The
   workflow resolves each VCF as `{paths.vcf_folder}/{vcf_basename}.vcf.gz`.

## config_vc.yaml reference

| Section | Key | Description |
|---------|-----|-------------|
| `paths` | `samples` | Path to the samples TSV file |
| `paths` | `vcf_folder` | Directory containing `.vcf.gz` input files |
| `paths` | `output_folder` | Root output directory (subdirs created per sample) |
| `paths` | `log_subdir` | Subdirectory name for rule logs |
| `variantcentrifuge` | `config_file` | Path to the variantcentrifuge JSON config |
| `variantcentrifuge` | `genes` | Gene filter (`"all"`, gene name, or gene-list file path) |
| `variantcentrifuge` | `log_level` | Log verbosity: `DEBUG`, `INFO`, `WARN`, `ERROR` |
| `variantcentrifuge` | `threads` | Threads per variantcentrifuge job |
| `variantcentrifuge` | `field_profile` | Field profile name or `null` for default |
| `variantcentrifuge` | `presets` | List of filter preset names |
| `igv` | `enabled` | Enable IGV report generation |
| `igv` | `reference` | IGV genome reference |
| `container` | `image` | Singularity/Apptainer image URI (empty string to disable) |

## Installing variantcentrifuge into the Snakemake conda env

The Snakemake conda environment (`workflow/envs/variantcentrifuge.yml`) provides
all external tool dependencies but does **not** install variantcentrifuge itself
(it is not published on PyPI). After Snakemake creates the conda env you must
install variantcentrifuge manually:

```bash
# Find the conda env prefix created by Snakemake:
snakemake --list-conda-envs
# Example output: .snakemake/conda/a1b2c3d4

# Option A: install from a local checkout
conda run -p .snakemake/conda/<hash> pip install /path/to/variantcentrifuge

# Option B: install directly from GitHub
conda run -p .snakemake/conda/<hash> pip install git+https://github.com/scholl-lab/variantcentrifuge.git@v0.12.0
```

Replace `.snakemake/conda/<hash>` with the env prefix path shown by
`snakemake --list-conda-envs`.

## snpEff reference name

The `"reference"` field in `config.json` must exactly match the name of the
local snpEff database. Common names include `GRCh37.p13`, `GRCh38.p14`,
`hg19`, and `hg38`. To check which databases are available:

```bash
snpEff databases | grep -i grch38
```

## Somatic / tumor-normal presets

The **legacy** `mutect2_TvsN` and `mutect2_TvsN_pass` presets assume the
standard Mutect2 sample order: **GEN[0] = tumor, GEN[1] = normal**. Verify
with:

```bash
bcftools query -l your.vcf.gz
```

The parameterized `somatic` / `somatic_pass` presets accept explicit sample
indices via `--tumor-sample-index` and `--normal-sample-index`. The CLI
defaults are `normal_idx=0` and `tumor_idx=1` (the opposite of Mutect2's
convention). For standard Mutect2 VCFs, pass:

```bash
variantcentrifuge ... --presets somatic_pass --tumor-sample-index 0 --normal-sample-index 1
```

## samples.tsv columns

| Column | Type | Description |
|--------|------|-------------|
| `sample` | string | Unique sample identifier |
| `vcf_basename` | string | VCF filename stem (without `.vcf.gz`) |
