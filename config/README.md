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

## samples.tsv columns

| Column | Type | Description |
|--------|------|-------------|
| `sample` | string | Unique sample identifier |
| `vcf_basename` | string | VCF filename stem (without `.vcf.gz`) |
