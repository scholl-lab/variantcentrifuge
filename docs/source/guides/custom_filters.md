# Custom Filter Development

Learn how to create custom SnpSift filter expressions for specialized analysis needs.

## Filter Syntax

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

## Example Custom Filters

```json
{
  "presets": {
    "custom_rare_high_impact": "(((gnomAD_exomes_AF < 0.0001) | (na gnomAD_exomes_AF)) & ((ANN[ANY].IMPACT has 'HIGH') | ((ANN[ANY].IMPACT has 'MODERATE') & (dbNSFP_CADD_phred >= 25))))"
  }
}
```