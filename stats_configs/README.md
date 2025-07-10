# Statistics Configuration Examples

This directory contains example statistics configuration files for VariantCentrifuge. These are **optional templates** that demonstrate how to create custom statistics for your variant analysis workflows.

## Overview

VariantCentrifuge automatically generates basic statistics using built-in defaults. These configuration files allow you to:

- Define custom dataset-level statistics (e.g., pathogenic variant counts)
- Create gene-level aggregations (e.g., burden analysis per gene)
- Generate grouped statistics with pivot tables
- Focus on specific analysis domains (clinical, population genetics, etc.)

## Available Configurations

### `clinical_focus_stats.json`
**Focus**: Clinical interpretation and pathogenicity assessment

**Key Statistics**:
- Pathogenic/Likely Pathogenic (P/LP) variant counts
- Variants of Uncertain Significance (VUS) counts  
- ACMG classification coverage
- Mean CADD scores for pathogenic variants
- ClinVar annotation density per gene
- Clinical significance distribution by gene
- Inheritance pattern correlations

**Use Case**: Clinical genetics labs, diagnostic workflows, pathogenicity assessment

### `example_custom_stats.json`
**Focus**: Comprehensive example showcasing all capabilities

**Key Statistics**:
- Ultra-rare variant analysis (AF < 0.0001)
- Impact type distributions and cross-tabulations
- CADD, REVEL, and conservation score statistics
- Gene burden metrics and deleteriousness predictions
- Allele frequency stratification
- Complex groupings and pivot tables

**Use Case**: Research workflows, method development, learning the configuration format

## Usage

### Basic Usage (Built-in Defaults)
```bash
# Uses automatic default statistics
variantcentrifuge --vcf input.vcf --gene-file genes.txt --xlsx
```

### Using Example Configurations
```bash
# Clinical focus
variantcentrifuge --vcf input.vcf --gene-file genes.txt \
  --stats-config stats_configs/clinical_focus_stats.json --xlsx

# Comprehensive example
variantcentrifuge --vcf input.vcf --gene-file genes.txt \
  --stats-config stats_configs/example_custom_stats.json --xlsx
```

### Creating Custom Configurations
```bash
# Copy and modify an example
cp stats_configs/clinical_focus_stats.json my_custom_stats.json
# Edit my_custom_stats.json to fit your needs
variantcentrifuge --vcf input.vcf --gene-file genes.txt \
  --stats-config my_custom_stats.json --xlsx
```

## Configuration Format

Statistics configurations use JSON format with three main sections:

```json
{
  "stats_version": "1.0",
  "description": "Configuration description",
  
  "dataset_stats": [
    {
      "name": "total_variants",
      "expression": "len(df)",
      "description": "Total number of variants"
    }
  ],
  
  "gene_stats": [
    {
      "name": "variant_count", 
      "expression": "len(group_df)",
      "groupby": "GENE",
      "description": "Number of variants per gene"
    }
  ],
  
  "grouped_stats": [
    {
      "name": "impact_by_gene",
      "expression": "size()",
      "groupby": ["GENE", "IMPACT"],
      "output_format": "pivot",
      "description": "Gene × Impact cross-tabulation"
    }
  ]
}
```

### Expression Context

Expressions are evaluated with these available variables:
- `df`: Full variant DataFrame
- `group_df`: Current group DataFrame (for gene_stats/grouped_stats)
- `pd`: Pandas library
- `np`: NumPy library  
- `len()`, `size()`: Length/size functions

### Required Columns

Use `required_columns` to ensure statistics only run when expected columns exist:

```json
{
  "name": "pathogenic_count",
  "expression": "df['ClinicalSignificance'].str.contains('athogenic').sum()",
  "required_columns": ["ClinicalSignificance"],
  "description": "Count of pathogenic variants"
}
```

## Output Integration

Statistics are automatically included in:
- **Excel Reports**: "Statistics" sheet with all computed metrics
- **TSV Files**: `--stats-output-file` for programmatic analysis
- **HTML Reports**: Summary statistics sections

## Best Practices

1. **Start with examples**: Copy and modify existing configurations
2. **Test expressions**: Validate with small datasets first
3. **Document clearly**: Use descriptive names and descriptions
4. **Handle missing data**: Use `required_columns` and null checks
5. **Keep it focused**: Create domain-specific configurations rather than overly complex ones

## Common Use Cases

### Population Genetics
Focus on allele frequencies, Hardy-Weinberg equilibrium, population stratification

### Functional Impact
Emphasize conservation scores, deleteriousness predictions, protein effects

### Gene Burden Analysis
Case/control comparisons, statistical tests, burden metrics

### Quality Control
Coverage metrics, batch effects, technical artifacts

## Notes

- These files are **examples only** - the tool works without them using built-in defaults
- Configurations are loaded fresh for each analysis
- Invalid expressions are logged but don't crash the pipeline
- Missing columns cause statistics to be skipped with warnings