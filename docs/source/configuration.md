# Configuration

VariantCentrifuge uses a JSON configuration file (`config.json`) to set default parameters and define reusable filter presets. This allows for flexible, reproducible analysis workflows.

## Configuration File Location

VariantCentrifuge looks for configuration files in the following order:

1. File specified with `--config` command-line option
2. `config.json` in the current working directory
3. Default configuration included with the package

## Configuration Structure

### Required Keys

These parameters must be provided either in the config file or via command-line arguments:

- **`reference`** (`str`) - Reference genome database for snpEff (e.g., "GRCh37.75", "GRCh38.99")
- **`filters`** (`str`) - SnpSift filter expression to select variants
- **`fields_to_extract`** (`str`) - Space-separated list of fields to extract via SnpSift

### Optional Keys

Configuration options with default values:

- **`interval_expand`** (`int`) - Number of bases to expand around genes. *Default: 0*
- **`add_chr`** (`bool`) - Add "chr" prefix to chromosome names. *Default: true*
- **`debug_level`** (`str`) - Logging level: "DEBUG", "INFO", "WARN", "ERROR". *Default: "INFO"*
- **`no_stats`** (`bool`) - Skip statistics computation. *Default: false*
- **`perform_gene_burden`** (`bool`) - Perform gene burden analysis. *Default: false*
- **`gene_burden_mode`** (`str`) - "samples" or "alleles". *Default: "alleles"*
- **`correction_method`** (`str`) - "fdr" or "bonferroni" for multiple testing correction. *Default: "fdr"*

### IGV Integration

- **`igv_enabled`** (`bool`) - Enable IGV.js integration. *Default: false*
- **`bam_mapping_file`** (`str`) - Path to TSV/CSV file mapping sample IDs to BAM files (required if igv_enabled=true)
- **`igv_reference`** (`str`) - Genome reference identifier for IGV (e.g., 'hg19', 'hg38')
- **`igv_fasta`** (`str`) - Path to local FASTA file for IGV reports (takes precedence over igv_reference)
- **`igv_ideogram`** (`str`) - Path to local ideogram file for IGV visualization
- **`igv_flanking`** (`int`) - Flanking bases around variants for IGV. *Default: 50*

## Filter Presets

The configuration system supports predefined filter presets that can be combined and reused. Presets are defined in the `presets` section:

```json
{
  "presets": {
    "rare": "(((dbNSFP_gnomAD_exomes_AF[0] < 0.0001) | (na dbNSFP_gnomAD_exomes_AC[0])) & ((dbNSFP_gnomAD_genomes_AF[0] < 0.0001) | (na dbNSFP_gnomAD_genomes_AC[0])))",
    "coding": "((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE'))",
    "pathogenic": "(((dbNSFP_clinvar_clnsig =~ '[Pp]athogenic') & !(dbNSFP_clinvar_clnsig =~ '[Cc]onflicting')) | ((ClinVar_CLNSIG =~ '[Pp]athogenic') & !(ClinVar_CLNSIG =~ '[Cc]onflicting')))"
  }
}
```

### Using Presets

Presets can be used with the `--preset` command-line option:

```bash
# Single preset
variantcentrifuge --preset rare --gene-name BRCA1 --vcf-file input.vcf

# Multiple presets (combined with AND)
variantcentrifuge --preset rare,coding --gene-name BRCA1 --vcf-file input.vcf

# Combine presets with custom filters
variantcentrifuge --preset rare --filters "QUAL >= 100" --gene-name BRCA1 --vcf-file input.vcf
```

### Built-in Presets

The default configuration includes many useful presets:

#### Rarity Filters
- **`super_rare`** - AC ≤ 2 in gnomAD exomes and genomes
- **`rare`** - AF < 0.0001 in gnomAD exomes and genomes  
- **`1percent`** - AF < 0.001 in gnomAD exomes and genomes
- **`5percent`** - AF < 0.05 in gnomAD exomes and genomes

#### Impact Filters
- **`high`** - High impact variants only
- **`moderate`** - Moderate impact variants only
- **`coding`** - High OR moderate impact variants
- **`high_or_moderate_or_low`** - Any protein-coding impact

#### Clinical Significance
- **`pathogenic`** - ClinVar pathogenic variants
- **`not_benign`** - Exclude ClinVar benign variants
- **`pathogenic_or_rare`** - Pathogenic OR rare variants

#### Quality Filters
- **`not_artefact`** - Quality-based artifact filtering
- **`mutect2_TvsN_pass`** - Tumor vs Normal Mutect2 filters with PASS only

## Example Configurations

### Basic Research Configuration

```json
{
  "reference": "GRCh37.75",
  "filters": "",
  "fields_to_extract": "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P dbNSFP_gnomAD_exomes_AC GEN[*].GT",
  "interval_expand": 1000,
  "add_chr": false,
  "perform_gene_burden": true
}
```

### Clinical Analysis Configuration

```json
{
  "reference": "GRCh38.99",
  "filters": "((dbNSFP_gnomAD_exomes_AF[0] < 0.01) & ((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE')))",
  "fields_to_extract": "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P dbNSFP_CADD_phred dbNSFP_REVEL_score ClinVar_CLNSIG GEN[*].GT",
  "add_chr": true,
  "igv_enabled": true,
  "igv_reference": "hg38"
}
```

### Somatic Variant Configuration

```json
{
  "reference": "GRCh38.99", 
  "filters": "(GEN[0].AF < 0.03) & (GEN[1].AF >= 0.05) & (FILTER = 'PASS')",
  "fields_to_extract": "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P GEN[*].GT GEN[*].AF GEN[*].DP",
  "presets": {
    "somatic_quality": "(GEN[0].AF < 0.03) & (GEN[1].AF >= 0.05) & (GEN[*].DP >= 50)",
    "cosmic_or_rare": "(((dbNSFP_gnomAD_exomes_AC[0] <= 2) | (na dbNSFP_gnomAD_exomes_AC[0])) | (exists ID & ID =~ 'COS'))"
  }
}
```

## External Database Links

Configure links to external databases that will be added to HTML reports:

```json
{
  "links": {
    "SpliceAI": "https://spliceailookup.broadinstitute.org/#variant={CHROM}-{POS}-{REF}-{ALT}&hg=19",
    "Franklin": "https://franklin.genoox.com/clinical-db/variant/snp/{CHROM}-{POS}-{REF}-{ALT}-hg19",
    "Varsome": "https://varsome.com/variant/hg19/{CHROM}-{POS}-{REF}-{ALT}",
    "gnomAD": "https://gnomad.broadinstitute.org/variant/{CHROM}-{POS}-{REF}-{ALT}",
    "ClinVar": "https://www.ncbi.nlm.nih.gov/clinvar/?term={CHROM}-{POS}-{REF}-{ALT}"
  }
}
```

## HTML Report Customization

Control which columns are hidden by default in HTML reports:

```json
{
  "html_report_default_hidden_columns": [
    "QUAL", "AC", "FEATUREID", "AA_POS", "AA_LEN"
  ],
  "html_report_truncate_settings": {
    "default_max_width_px": 120,
    "columns_for_hover_expand": ["HGVS_P", "HGVS_C", "EFFECT"],
    "column_specific_max_widths_px": {
      "GT": 250
    }
  }
}
```

## Configuration Validation

VariantCentrifuge validates your configuration and provides helpful error messages:

- Missing required keys will trigger clear error messages
- Invalid preset references will be caught at startup
- Filter expression syntax is validated when possible

## Best Practices

1. **Use version control** - Store your configuration files in version control alongside your analysis scripts
2. **Environment-specific configs** - Use different configurations for development, testing, and production
3. **Document custom presets** - Add comments explaining complex filter logic
4. **Test filters** - Validate filter expressions on small datasets before large analyses
5. **Modular approach** - Use presets to build complex filters from simpler components