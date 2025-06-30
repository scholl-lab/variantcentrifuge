# Variant Scoring Guide

This guide covers how to create and use custom variant scoring models in VariantCentrifuge.

## Overview

VariantCentrifuge's scoring system allows you to apply custom mathematical models to score variants based on their annotations. This is useful for:

- Prioritizing variants based on multiple criteria
- Implementing published scoring algorithms
- Creating disease-specific variant rankings
- Combining multiple evidence types into a single score

## How Scoring Works

The scoring system:
1. Loads configuration files that define variable mappings and formulas
2. Maps VCF annotation columns to formula variables
3. Applies the formula to calculate scores
4. Adds score columns to the output

## Configuration Structure

A scoring configuration consists of two JSON files in a directory:

### 1. Variable Assignment Configuration

`variable_assignment_config.json` maps VCF columns to formula variables:

```json
{
  "variables": {
    "VCF_COLUMN_NAME": "formula_variable|default:value"
  }
}
```

### 2. Formula Configuration

`formula_config.json` defines the scoring formulas:

```json
{
  "formulas": [
    {
      "score_name": "formula_expression"
    }
  ]
}
```

## Creating a Custom Scoring Model

### Step 1: Identify Required Annotations

First, determine which VCF annotations you need for your scoring model:

- Allele frequencies (e.g., `dbNSFP_gnomAD_exomes_AF`)
- Pathogenicity scores (e.g., `dbNSFP_CADD_phred`, `dbNSFP_REVEL_score`)
- Variant effects (e.g., `ANN[0].EFFECT`, `ANN[0].IMPACT`)
- Conservation scores (e.g., `dbNSFP_phyloP100way_vertebrate`)

### Step 2: Create Variable Mappings

Create a directory for your scoring configuration:

```bash
mkdir -p scoring/my_custom_score
```

Create `variable_assignment_config.json`:

```json
{
  "variables": {
    "dbNSFP_gnomAD_exomes_AF": "af_exomes|default:0.0",
    "dbNSFP_CADD_phred": "cadd_score|default:0.0",
    "dbNSFP_REVEL_score": "revel_score|default:0.0",
    "ANN[0].IMPACT": "impact|default:'MODIFIER'"
  }
}
```

### Step 3: Define Your Formula

Create `formula_config.json`:

```json
{
  "formulas": [
    {
      "my_variant_score": "(cadd_score / 40) * 0.4 + (revel_score * 0.3) + ((1 - af_exomes) * 0.3)"
    }
  ]
}
```

### Step 4: Test Your Scoring

```bash
variantcentrifuge \
  --gene-name BRCA1 \
  --vcf-file test.vcf.gz \
  --scoring-config-path scoring/my_custom_score \
  --output-file scored_test.tsv
```

## Formula Syntax

Formulas use pandas eval syntax with these capabilities:

### Basic Operations
- Arithmetic: `+`, `-`, `*`, `/`, `**`
- Comparisons: `==`, `!=`, `<`, `>`, `<=`, `>=`
- Logical: `&` (and), `|` (or), `~` (not)

### Conditional Logic
```python
# If-then-else using multiplication
"(impact == 'HIGH') * 1.0 + (impact == 'MODERATE') * 0.5"

# Multiple conditions
"((af < 0.01) & (cadd > 20)) * 1.0"
```

### String Operations
```python
# String matching
"consequence.str.contains('missense') * 0.5"

# Exact match
"(consequence == 'stop_gained') * 1.0"
```

### Mathematical Functions
```python
# Exponential (for logistic regression)
"1 / (1 + 2.718281828459045 ** (-linear_combination))"

# Min/Max bounds
"cadd_score.clip(0, 40) / 40"
```

## Example: Logistic Regression Model

Here's a complete example implementing a logistic regression model:

### Variable Mappings
```json
{
  "variables": {
    "dbNSFP_gnomAD_exomes_AF": "gnomad_af|default:0.0",
    "dbNSFP_CADD_phred": "cadd|default:0.0",
    "ANN[0].EFFECT": "effect|default:''",
    "ANN[0].IMPACT": "impact|default:''"
  }
}
```

### Logistic Formula
```json
{
  "formulas": [
    {
      "pathogenicity_score": "1 / (1 + 2.718281828459045 ** (-(2.5 + (cadd - 15) * 0.1 + (gnomad_af * -50) + ((effect == 'missense_variant') * 0.5) + ((impact == 'HIGH') * 2))))"
    }
  ]
}
```

## Advanced Techniques

### Multiple Scores

You can calculate multiple scores in one configuration:

```json
{
  "formulas": [
    {
      "conservation_score": "(phylop + phastcons) / 2"
    },
    {
      "pathogenicity_score": "(cadd / 40) * 0.5 + (revel * 0.5)"
    },
    {
      "combined_score": "conservation_score * 0.3 + pathogenicity_score * 0.7"
    }
  ]
}
```

### Handling Missing Data

Use default values in variable mappings:

```json
{
  "variables": {
    "dbNSFP_CADD_phred": "cadd|default:10.0",
    "ClinVar_CLNSIG": "clinvar|default:'not_provided'"
  }
}
```

### Categorical Variables

Convert categories to numeric values:

```json
{
  "formulas": [
    {
      "impact_numeric": "(impact == 'HIGH') * 4 + (impact == 'MODERATE') * 3 + (impact == 'LOW') * 2 + (impact == 'MODIFIER') * 1"
    }
  ]
}
```

## Real-World Example: Nephro Variant Score

The included `scoring/nephro_variant_score` configuration implements a sophisticated model for kidney disease variants:

```bash
# Use the nephro variant score
variantcentrifuge \
  --gene-file kidney_genes.txt \
  --vcf-file patient.vcf.gz \
  --scoring-config-path scoring/nephro_variant_score \
  --preset rare,coding \
  --html-report \
  --output-file kidney_analysis.tsv
```

This model:
- Uses logistic regression with multiple predictors
- Handles gnomAD frequencies from exomes and genomes
- Incorporates CADD scores
- Considers variant consequences (missense, frameshift, etc.)
- Accounts for predicted impact levels

## Best Practices

1. **Validate Your Formulas**
   - Test on known pathogenic and benign variants
   - Check score distributions make sense
   - Verify handling of missing data

2. **Document Your Model**
   - Include comments explaining the rationale
   - Document expected score ranges
   - Provide interpretation guidelines

3. **Version Control**
   - Track scoring configurations in git
   - Tag releases of scoring models
   - Document changes between versions

4. **Performance Considerations**
   - Keep formulas reasonably simple
   - Avoid deeply nested conditions
   - Test on large datasets

## Troubleshooting

### Common Issues

1. **"Column not found" errors**
   - Check VCF has the required annotations
   - Verify column names match exactly
   - Use appropriate default values

2. **Formula syntax errors**
   - Test formulas incrementally
   - Check parentheses are balanced
   - Verify string operations use proper syntax

3. **Unexpected scores**
   - Check default values are appropriate
   - Verify numeric conversions
   - Test edge cases (0, 1, missing values)

### Debugging Tips

1. Keep intermediate TSV files to inspect data:
   ```bash
   variantcentrifuge --keep-intermediates ...
   ```

2. Start with simple formulas and build complexity:
   ```json
   {"test_score": "cadd"}  // Start simple
   {"test_score": "cadd / 40"}  // Add normalization
   {"test_score": "(cadd / 40) * 0.5 + (revel * 0.5)"}  // Combine
   ```

3. Use the Python REPL to test pandas expressions:
   ```python
   import pandas as pd
   df = pd.DataFrame({'cadd': [10, 20, 30]})
   df.eval("cadd / 40")
   ```

## Integration with Other Features

Scoring works seamlessly with other VariantCentrifuge features:

- **Filtering**: Apply filters before or after scoring
- **Gene Burden**: Use scores in burden testing
- **Reports**: Scores appear in HTML reports
- **IGV**: Visualize high-scoring variants

## Contributing Scoring Models

If you develop a useful scoring model, consider contributing it:

1. Create a well-documented configuration
2. Include example use cases
3. Provide validation results
4. Submit via GitHub pull request

## Further Reading

- [Configuration Guide](../configuration.md) - General configuration options
- [API Reference](../api/scoring.md) - Technical documentation
- [pandas.eval documentation](https://pandas.pydata.org/docs/reference/api/pandas.eval.html) - Formula syntax reference