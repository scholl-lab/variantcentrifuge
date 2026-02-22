# Frequently Asked Questions

## General Questions

### What is VariantCentrifuge?

VariantCentrifuge is a comprehensive Python-based command-line tool designed for filtering, analyzing, and interpreting genetic variants from VCF (Variant Call Format) files. It's specifically built for clinical geneticists, bioinformaticians, and researchers who need to identify disease-causing variants efficiently.

### What file formats does VariantCentrifuge support?

VariantCentrifuge supports:
- VCF files (Variant Call Format)
- Compressed VCF files (VCF.gz)
- PED files for family/pedigree information
- TSV/CSV files for phenotype data
- BED files for custom region annotation
- JSON files for gene metadata

### Can I use VariantCentrifuge for clinical variant interpretation?

Yes, VariantCentrifuge includes features specifically designed for clinical variant analysis:
- ACMG classification integration
- ClinVar annotation support
- Pathogenicity prediction scores
- Clinical significance filtering
- Professional reporting suitable for clinical review

### What are the system requirements?

- **Operating System:** Linux, macOS, or Windows (via WSL)
- **Python:** Version 3.7 or higher
- **External Tools:** bcftools, snpEff, SnpSift, bedtools (must be in PATH)
- **Memory:** Minimum 8GB RAM recommended (16GB+ for large cohorts)
- **Storage:** Varies by dataset size; intermediate files can be compressed

## Installation and Setup

### How do I install VariantCentrifuge?

The quickest method is Docker, which includes all external tools:

```bash
docker pull ghcr.io/scholl-lab/variantcentrifuge:latest
docker run --rm -v ./data:/data ghcr.io/scholl-lab/variantcentrifuge:latest \
  --gene-name BRCA1 --vcf-file /data/input.vcf.gz --output-file /data/output.tsv
```

For a native installation, use conda/mamba:

```bash
# Create environment with dependencies
mamba create -y -n variantcentrifuge bcftools snpsift snpeff bedtools
mamba activate variantcentrifuge

# Install VariantCentrifuge
pip install variantcentrifuge
```

For development installation:
```bash
git clone https://github.com/scholl-lab/variantcentrifuge.git
cd variantcentrifuge
pip install -e .
```

### Why do I get "command not found" errors?

This typically means the required external tools (bcftools, snpEff, SnpSift, bedtools) are not installed or not in your PATH. Ensure all dependencies are installed and accessible:

```bash
# Check if tools are available
which bcftools snpEff SnpSift bedtools
```

## Usage and Features

### How do I filter variants for a specific gene?

Basic gene filtering example:

```bash
variantcentrifuge \
  --gene-name BRCA1 \
  --vcf-file input.vcf.gz \
  --output-file brca1_variants.tsv
```

For multiple genes:
```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --output-file output.tsv
```

### What filter presets are available?

VariantCentrifuge includes 20+ preset filter combinations:
- `rare`: Variants with AF < 1%
- `super_rare`: Variants with AF < 0.1%
- `coding`: Variants in coding regions
- `high_impact`: High impact variants only
- `pathogenic`: ClinVar pathogenic/likely pathogenic
- `rare,coding,high_impact`: Combined filters

Use `--preset` to apply:
```bash
variantcentrifuge --preset rare,coding,pathogenic ...
```

### How does compound heterozygous detection work?

VariantCentrifuge automatically detects compound heterozygous variants when:
1. A PED file is provided with family information
2. Two or more heterozygous variants exist in the same gene
3. The variants follow expected inheritance patterns

The optimized algorithm is 10-50x faster than traditional approaches and handles large gene panels efficiently.

### Can I analyze multiple samples or families?

Yes! VariantCentrifuge supports:
- Single sample analysis
- Family/trio analysis with PED files
- Case-control cohort studies
- Multi-sample VCFs

Example with family analysis:
```bash
variantcentrifuge \
  --gene-file panel.txt \
  --vcf-file family.vcf.gz \
  --ped family.ped \
  --inheritance-mode full \
  --output-file results.tsv
```

## Performance and Optimization

### How can I improve processing speed for large VCF files?

Several optimization strategies:

1. **Pre-filter with bcftools:**
   ```bash
   variantcentrifuge --bcftools-filter "QUAL>=30" ...
   ```

2. **Use chunked processing:**
   ```bash
   variantcentrifuge --chunks 1000 ...
   ```

3. **Compress intermediate files:**
   ```bash
   variantcentrifuge --gzip-intermediates ...
   ```

4. **Focus on specific regions:**
   ```bash
   variantcentrifuge --interval-expansion 0 ...
   ```

### What's the maximum VCF file size VariantCentrifuge can handle?

There's no hard limit, but performance considerations apply:
- Files <1GB: Process normally
- Files 1-10GB: Use chunking (`--chunks 5000`)
- Files >10GB: Pre-filter with bcftools and use chunking
- Memory usage scales with number of samples and variants per chunk

## Troubleshooting

### Why are some variants missing from my output?

Common reasons:
1. Variants filtered by quality or coverage thresholds
2. Gene name mismatches (check gene symbol aliases)
3. Variants outside gene boundaries (adjust `--interval-expansion`)
4. Filter criteria too stringent

Debug with:
```bash
variantcentrifuge --no-filtering --keep-intermediates ...
```

### How do I debug inheritance pattern detection?

Use verbose inheritance output:
```bash
variantcentrifuge --inheritance-mode full --ped family.ped ...
```

This provides detailed JSON output showing:
- Genotypes for all family members
- Segregation analysis results
- Confidence scores
- Reasons for pattern assignment

### Can I resume an interrupted analysis?

Yes, with checkpoint support:

```bash
# Enable checkpointing
variantcentrifuge --enable-checkpoint ...

# Resume if interrupted
variantcentrifuge --enable-checkpoint --resume ...
```

## Advanced Features

### How do I create custom scoring formulas?

Create a scoring configuration directory:

```bash
scoring/
├── my_score/
│   ├── variable_assignment_config.json
│   └── formula_config.json
```

Then apply:
```bash
variantcentrifuge --scoring-config-path scoring/my_score ...
```

### Can I integrate custom annotations?

Yes, multiple ways:

1. **BED file regions:**
   ```bash
   variantcentrifuge --annotate-bed custom_regions.bed ...
   ```

2. **Gene lists:**
   ```bash
   variantcentrifuge --annotate-gene-list important_genes.txt ...
   ```

3. **JSON metadata:**
   ```bash
   variantcentrifuge --annotate-json-genes metadata.json \
     --json-gene-mapping '{"identifier":"symbol","dataFields":["panel","moi"]}' ...
   ```

### How do I perform gene burden analysis?

For case-control studies:

```bash
variantcentrifuge \
  --perform-gene-burden \
  --case-samples-file cases.txt \
  --control-samples-file controls.txt \
  --output-file burden_analysis.tsv
```

This performs Fisher's exact test with multiple testing correction.

### How do I switch between dbNSFP versions?

Use the `--field-profile` option. The default profile (`dbnsfp4`) uses separate gnomAD exomes/genomes fields. For dbNSFP v5.x with joint gnomAD 4.1 fields, use `dbnsfp5`:

```bash
variantcentrifuge --field-profile dbnsfp5 --preset rare,coding ...
```

List available profiles:
```bash
variantcentrifuge --list-field-profiles
```

See the [Field Profiles](configuration.md#field-profiles) documentation for details on adding custom profiles.

### Can I run VariantCentrifuge without installing bioinformatics tools locally?

Yes, the Docker image includes bcftools, snpEff, SnpSift, bedtools, and all Python dependencies. Mount your data directory and run:

```bash
docker run --rm -v ./data:/data ghcr.io/scholl-lab/variantcentrifuge:latest \
  --gene-name BRCA1 --vcf-file /data/input.vcf.gz --output-file /data/output.tsv
```

See the [Docker installation guide](installation.md#method-4-docker-recommended-for-quick-setup) for complete setup instructions.

## Association Testing

### How do I run a gene burden test with covariate adjustment?

Use `--perform-association` with `--covariate-file` instead of (or in addition to) `--perform-gene-burden`:

```bash
variantcentrifuge \
  --gene-file target_genes.txt \
  --vcf-file cohort.vcf.gz \
  --perform-association \
  --association-tests logistic_burden \
  --covariate-file covariates.tsv \
  --covariates age,sex,batch \
  --case-samples-file cases.txt \
  --control-samples-file controls.txt \
  --output-file burden_results.tsv
```

The covariate file should be a TSV or CSV with sample IDs in the first column and covariate values in subsequent columns (header required). Logistic burden regression fits a weighted burden score as the predictor with covariates, reporting beta coefficients and standard errors. For quantitative traits, use `--trait-type quantitative` with `--association-tests linear_burden`. See the [Association Testing Guide](guides/association_testing.md) for full details.

### What is ACAT-O and why should I use `acat_o_corrected_p_value`?

ACAT-O (Aggregated Cauchy Association Test — Omnibus) combines p-values from multiple association tests per gene into a single omnibus p-value using the Cauchy combination method (Liu and Xie, 2020). It is robust to the unknown correlation between tests and does not require simulation.

When you run multiple tests (e.g., `--association-tests fisher,logistic_burden,skat_o`), each test may capture different genetic architectures: Fisher detects dominant carrier effects, burden tests detect directional cumulative effects, and SKAT-O detects heterogeneous effects without sign constraints. ACAT-O combines evidence from all active tests into a single per-gene p-value. A single Benjamini-Hochberg FDR correction is applied to `acat_o_p_value` across all genes to produce `acat_o_corrected_p_value`. This is the primary significance measure to report; individual test p-values serve for diagnostic signal decomposition only. See the [Association Testing Guide](guides/association_testing.md) for the test selection reference table.

### Do I need R installed for SKAT or COAST?

No. Since v0.15.0, the default backend for both SKAT (`--skat-backend python`) and COAST (`--coast-backend python`) is a pure Python implementation using numpy and scipy. No R or rpy2 installation is required.

The R backends (`--skat-backend r`, `--coast-backend r`) are deprecated as of v0.15.0 and may be removed in a future release. They require rpy2, the SKAT R package, and (for COAST) the AllelicSeries R package. If you have existing workflows relying on the R backends, migrate to the Python backends — the p-values agree within 10% on real data. See the [Association Testing Guide](guides/association_testing.md) for backend migration details.
