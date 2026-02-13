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

### Variant Scoring

- **`scoring_config_path`** (`str`) - Path to directory containing scoring configuration files (variable_assignment_config.json and formula_config.json)

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

## Field Profiles

Field profiles handle differences between annotation database versions — for example, dbNSFP v4.x and v5.x use different gnomAD field names. Instead of maintaining separate config files, you define profiles once and switch between them at runtime.

### How It Works

Presets can contain parameterized templates using `{{fragment_name:param}}` syntax. At resolution time, the fragment is looked up in the active profile and expanded with the parameter value.

For example, the `rare` preset:

```json
"rare": "({{gnomad_af_below:0.0001}})"
```

expands to different SnpSift expressions depending on the active profile:

- **dbnsfp4**: `(((dbNSFP_gnomAD_exomes_AF[0] < 0.0001) | ...) & ((dbNSFP_gnomAD_genomes_AF[0] < 0.0001) | ...))`
- **dbnsfp5**: `(((dbNSFP_gnomAD4.1_joint_AF[0] < 0.0001) | (na dbNSFP_gnomAD4.1_joint_AC[0])))`

Presets without `{{...}}` patterns pass through unchanged (fully backward compatible).

### Selecting a Profile

Use the `--field-profile` CLI option:

```bash
# Use dbNSFP v5.x field names
variantcentrifuge --field-profile dbnsfp5 --preset rare,coding ...

# List available profiles
variantcentrifuge --list-field-profiles
```

The default profile is set by `default_field_profile` in `config.json` (ships as `dbnsfp4`).

### Profile Configuration Structure

Profiles are defined in the `field_profiles` section of `config.json`:

```json
{
  "default_field_profile": "dbnsfp4",
  "field_profiles": {
    "dbnsfp4": {
      "description": "dbNSFP 4.x (separate gnomAD exomes/genomes fields)",
      "fragments": {
        "gnomad_af_below": "((dbNSFP_gnomAD_exomes_AF[0] < {0}) | (na dbNSFP_gnomAD_exomes_AC[0])) & ((dbNSFP_gnomAD_genomes_AF[0] < {0}) | (na dbNSFP_gnomAD_genomes_AC[0]))",
        "gnomad_af_below_or": "((dbNSFP_gnomAD_exomes_AF[0] < {0}) | (na dbNSFP_gnomAD_exomes_AC[0])) | ((dbNSFP_gnomAD_genomes_AF[0] < {0}) | (na dbNSFP_gnomAD_genomes_AC[0]))",
        "gnomad_ac_atmost": "((dbNSFP_gnomAD_exomes_AC[0] <= {0}) | (na dbNSFP_gnomAD_exomes_AC[0])) & ((dbNSFP_gnomAD_genomes_AC[0] <= {0}) | (na dbNSFP_gnomAD_genomes_AC[0]))"
      },
      "fields_to_extract": "dbNSFP_gnomAD_exomes_AF dbNSFP_gnomAD_genomes_AF dbNSFP_gnomAD_exomes_AC dbNSFP_gnomAD_genomes_AC",
      "hidden_columns": ["dbNSFP_gnomAD_exomes_AC", "dbNSFP_gnomAD_genomes_AC"]
    },
    "dbnsfp5": {
      "description": "dbNSFP 5.x (joint gnomAD 4.1 fields)",
      "fragments": {
        "gnomad_af_below": "((dbNSFP_gnomAD4.1_joint_AF[0] < {0}) | (na dbNSFP_gnomAD4.1_joint_AC[0]))",
        "gnomad_af_below_or": "((dbNSFP_gnomAD4.1_joint_AF[0] < {0}) | (na dbNSFP_gnomAD4.1_joint_AC[0]))",
        "gnomad_ac_atmost": "((dbNSFP_gnomAD4.1_joint_AC[0] <= {0}) | (na dbNSFP_gnomAD4.1_joint_AC[0]))"
      },
      "fields_to_extract": "dbNSFP_gnomAD4.1_joint_AF dbNSFP_gnomAD4.1_joint_AC",
      "hidden_columns": ["dbNSFP_gnomAD4.1_joint_AC"]
    }
  }
}
```

Each profile can define:
- **`fragments`** — Named templates with `{0}` placeholders for parameter substitution
- **`fields_to_extract`** — Profile-specific fields appended to the base field list
- **`hidden_columns`** — Columns hidden by default in HTML reports

### Adding Custom Profiles

To support a new annotation version, add a new profile to `field_profiles` in your config file with appropriate fragments. No code changes are needed — just define the fragment templates for your field naming convention.

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

## Variant Scoring Configuration

VariantCentrifuge supports custom variant scoring through configuration files. The scoring system requires two JSON files in a directory:

### Variable Assignment Configuration

`variable_assignment_config.json` maps VCF annotation fields to formula variables:

```json
{
  "variables": {
    "dbNSFP_gnomAD_exomes_AF": "gnomade_variant|default:0.0",
    "dbNSFP_gnomAD_genomes_AF": "gnomadg_variant|default:0.0",
    "dbNSFP_CADD_phred": "cadd_phred_variant|default:0.0",
    "ANN[0].EFFECT": "consequence_terms_variant|default:''",
    "ANN[0].IMPACT": "impact_variant|default:''"
  }
}
```

Each mapping specifies:
- **Key**: The column name in your VCF data
- **Value**: The variable name for the formula, with optional default value

### Formula Configuration

`formula_config.json` contains the scoring formulas:

```json
{
  "formulas": [
    {
      "score_name": "formula_expression"
    }
  ]
}
```

Formulas use pandas eval syntax and can include:
- Mathematical operations
- Conditional logic (`(condition) * value`)
- String operations (`.str.contains()`)
- Any valid pandas expression

### Example: Nephro Variant Score

A complete scoring configuration for kidney disease variants:

```json
// variable_assignment_config.json
{
  "variables": {
    "dbNSFP_gnomAD_exomes_AF": "gnomade_variant|default:0.0",
    "dbNSFP_gnomAD_genomes_AF": "gnomadg_variant|default:0.0",
    "dbNSFP_CADD_phred": "cadd_phred_variant|default:0.0",
    "ANN[0].EFFECT": "consequence_terms_variant|default:''",
    "ANN[0].IMPACT": "impact_variant|default:''"
  }
}

// formula_config.json
{
  "formulas": [
    {
      "nephro_variant_score": "1 / (1 + 2.718281828459045 ** (-((-36.30796) + ((gnomade_variant - 0.00658) / 0.05959) * (-309.33539) + ((gnomadg_variant - 0.02425) / 0.11003) * (-2.54581) + (((consequence_terms_variant == 'missense_variant') * 1.0 - 0.24333) / 0.42909) * (-1.14313) + ((cadd_phred_variant - 12.47608) / 11.78359) * 2.68520 + ((((impact_variant == 'HIGH') * 4 + (impact_variant == 'MODERATE') * 3 + (impact_variant == 'LOW') * 2 + (impact_variant == 'MODIFIER') * 1) - 2.49999) / 1.11804) * 3.14822)))"
    }
  ]
}
```

### Using Scoring

To apply scoring to your analysis:

```bash
variantcentrifuge \
  --gene-name GENE \
  --vcf-file input.vcf.gz \
  --scoring-config-path /path/to/scoring/config/dir \
  --output-file scored_variants.tsv
```

The scoring module will:
1. Load the configuration files from the specified directory
2. Map VCF columns to formula variables
3. Apply the formula to calculate scores
4. Add score columns to the output

## Best Practices

1. **Use version control** - Store your configuration files in version control alongside your analysis scripts
2. **Environment-specific configs** - Use different configurations for development, testing, and production
3. **Document custom presets** - Add comments explaining complex filter logic
4. **Test filters** - Validate filter expressions on small datasets before large analyses
5. **Modular approach** - Use presets to build complex filters from simpler components
6. **Test scoring formulas** - Validate scoring on known variants before production use
