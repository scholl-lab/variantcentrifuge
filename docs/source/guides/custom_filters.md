# Custom Filter Development

Learn how to create custom filter expressions for specialized analysis needs. VariantCentrifuge supports multiple filtering approaches at different stages of the pipeline.

## Filtering Options Overview

### 1. SnpSift Filters (`--filters`)
Applied early in the pipeline on VCF data. Use SnpSift's Java-like expression syntax.

### 2. bcftools Pre-filters (`--bcftools-prefilter`)
Applied during variant extraction for performance. Uses bcftools expression syntax.

### 3. Final Filters (`--final-filter`)
Applied at the end of the pipeline on the final TSV. Uses pandas query syntax.

## SnpSift Filter Syntax

SnpSift uses a Java-like expression syntax:

```bash
# Basic comparisons
"QUAL >= 30"
"AC[0] <= 5"
"AF[0] < 0.01"

# Field existence
"exists ClinVar_CLNSIG"
"na gnomAD_exomes_AF"

# String matching
"ClinVar_CLNSIG =~ '[Pp]athogenic'"
"ANN[ANY].IMPACT has 'HIGH'"

# Complex combinations
"((QUAL >= 30) & (AC[0] <= 5)) | (ClinVar_CLNSIG =~ '[Pp]athogenic')"
```

## bcftools Pre-filter Syntax

bcftools uses a different syntax optimized for speed:

```bash
# Basic filters
--bcftools-prefilter 'FILTER="PASS"'
--bcftools-prefilter 'QUAL>30'
--bcftools-prefilter 'INFO/AC<10'
--bcftools-prefilter 'INFO/AF<0.01'

# Combined filters
--bcftools-prefilter 'FILTER="PASS" && QUAL>30 && INFO/AC<10'
--bcftools-prefilter 'INFO/AF<0.001 || INFO/AC<5'
```

## Final Filter Syntax (pandas query)

The `--final-filter` option uses pandas query syntax, allowing filtering on any column including computed values:

```bash
# Numeric comparisons
--final-filter 'inheritance_score > 0.8'
--final-filter 'CADD_phred >= 20'
--final-filter 'Inheritance_Confidence > 0.9'

# String equality
--final-filter 'IMPACT == "HIGH"'
--final-filter 'Inheritance_Pattern == "de_novo"'

# String contains
--final-filter 'Custom_Annotation.str.contains("cancer_panel")'
--final-filter 'GENE.str.startswith("BRCA")'

# IN operator
--final-filter 'Inheritance_Pattern in ["de_novo", "compound_heterozygous"]'
--final-filter 'GENE in ["BRCA1", "BRCA2", "TP53"]'

# Complex expressions
--final-filter '(IMPACT == "HIGH" and inheritance_score > 0.7) or CADD_phred > 30'
--final-filter 'Inheritance_Pattern == "de_novo" and Inheritance_Confidence > 0.8'
```

## Filtering Strategy Examples

### Performance-Optimized Filtering
```bash
# Use bcftools pre-filter to reduce data early
variantcentrifuge \
  --bcftools-prefilter 'FILTER="PASS" && INFO/AC<20' \
  --preset rare,coding \
  --final-filter 'inheritance_score > 0.5' \
  ...
```

### Score-Based Filtering
```bash
# Use late filtering or final filter for computed columns
variantcentrifuge \
  --scoring-config-path scoring/my_model \
  --late-filtering \
  --filters "my_score > 0.7" \
  ...

# Or use final filter
variantcentrifuge \
  --scoring-config-path scoring/my_model \
  --final-filter 'my_score > 0.7 and IMPACT != "LOW"' \
  ...
```

### Inheritance-Based Filtering
```bash
# Filter for high-confidence de novo variants
variantcentrifuge \
  --ped family.ped \
  --inheritance-mode columns \
  --final-filter 'Inheritance_Pattern == "de_novo" and Inheritance_Confidence > 0.9' \
  ...
```

## Example Custom Filter Presets

```json
{
  "presets": {
    "custom_rare_high_impact": "(((gnomAD_exomes_AF < 0.0001) | (na gnomAD_exomes_AF)) & ((ANN[ANY].IMPACT has 'HIGH') | ((ANN[ANY].IMPACT has 'MODERATE') & (dbNSFP_CADD_phred >= 25))))"
  }
}
```
