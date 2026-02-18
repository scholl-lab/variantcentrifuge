# Privacy and Pseudonymization

VariantCentrifuge includes comprehensive pseudonymization features to help you share genomic data while protecting participant privacy.

## Overview

The pseudonymization feature replaces original sample identifiers with consistent, non-identifiable pseudonyms throughout all output files. This enables data sharing for publication or collaboration while maintaining participant privacy.

## Quick Start

Basic pseudonymization with sequential IDs:
```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file cohort.vcf \
  --output-file results.tsv \
  --pseudonymize \
  --pseudonymize-table mapping.tsv
```

## Naming Schemas

### Sequential Schema (Default)

Simple numbered identifiers:
```bash
--pseudonymize-schema sequential \
--pseudonymize-prefix STUDY_A
# Generates: STUDY_A_001, STUDY_A_002, etc.
```

### Categorical Schema

Groups samples by phenotype or other metadata:
```bash
--pseudonymize-schema categorical \
--pseudonymize-category-field phenotype
# Generates: CASE_001, CASE_002, CONTROL_001, etc.
```

### Anonymous Schema

Letter-number combinations or hash-based IDs:
```bash
--pseudonymize-schema anonymous
# Generates: A001, A002, B001, etc.
```

### Custom Schema

Define your own pattern:
```bash
--pseudonymize-schema custom \
--pseudonymize-pattern "{study}_{phenotype}_{index:04d}"
# Generates: STUDY1_CASE_0001, etc.
```

## Security Best Practices

**Critical Security Information**

### 1. Mapping Table Security

The pseudonymization table contains the key to re-identify your data:
- Store it separately from pseudonymized outputs
- Never share it with the pseudonymized data unless authorized
- Consider encrypting the mapping file
- Keep it in a secure location with restricted access

### 2. Archive Behavior

The `--archive-results` flag does NOT include the mapping table:
- Mapping tables are saved in the parent directory by design
- This prevents accidental sharing of identifying information
- Always verify archive contents before sharing

### 3. PED File Handling

Use `--pseudonymize-ped` to create shareable pedigree files:
- Original PED files should never be shared
- Family relationships are preserved with pseudonymized IDs
- The pseudonymized PED is saved alongside the mapping table

## Integration with Other Features

### With Inheritance Analysis

Pseudonymization automatically handles inheritance pattern outputs:
- Sample IDs in `Inheritance_Samples` column are replaced
- Compound heterozygous partner information is preserved
- Family relationships remain interpretable

### With HTML Reports

The HTML report will use pseudonymized IDs throughout:
- Sample columns in variant tables
- IGV integration will show pseudonymized labels
- Download links will use pseudonymized filenames

### With Excel Output

Excel files automatically use pseudonymized data:
- All sheets will contain pseudonymized IDs
- Links and references remain functional

## Example Workflows

### Publication-Ready Analysis

For a case-control study ready for publication:
```bash
variantcentrifuge \
  --gene-file cancer_genes.txt \
  --vcf-file cohort.vcf \
  --output-file publication_results.tsv \
  --ped cohort.ped \
  --pseudonymize \
  --pseudonymize-schema categorical \
  --pseudonymize-table ../SECURE/id_mapping.tsv \
  --pseudonymize-ped \
  --html-report \
  --archive-results

# Results:
# 1. ./publication_results.tsv - pseudonymized variant data
# 2. ./report/ - HTML report with pseudonymized IDs
# 3. ../SECURE/id_mapping.tsv - mapping table (NOT in archive)
# 4. ../SECURE/id_mapping_pedigree.ped - pseudonymized PED
# 5. ../variantcentrifuge_results_*.tar.gz - shareable archive
```

### Family Study with Custom IDs

For trio or family studies:
```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file families.vcf \
  --ped families.ped \
  --output-file results.tsv \
  --pseudonymize \
  --pseudonymize-schema custom \
  --pseudonymize-pattern "FAM{family}_{role}_{index:02d}" \
  --pseudonymize-table family_mapping.tsv
```

### Multi-Site Collaboration

For multi-site studies with complete de-identification:
```bash
variantcentrifuge \
  --gene-file shared_genes.txt \
  --vcf-file site_data.vcf \
  --output-file site_results.tsv \
  --pseudonymize \
  --pseudonymize-schema anonymous \
  --pseudonymize-table ../restricted/site_mapping.tsv \
  --xlsx \
  --html-report
```

## Reversing Pseudonymization

For internal use, you can reverse the process using the mapping table:

```python
import pandas as pd

# Load the mapping table
mapping = pd.read_csv('mapping.tsv', sep='\t')
reverse_map = dict(zip(mapping['pseudonym_id'], mapping['original_id']))

# Load pseudonymized results
results = pd.read_csv('pseudonymized_results.tsv', sep='\t')

# Function to reverse pseudonymization in GT column
def reverse_gt(gt_value):
    if pd.isna(gt_value):
        return gt_value
    for pseudo, orig in reverse_map.items():
        gt_value = gt_value.replace(pseudo + '(', orig + '(')
    return gt_value

# Apply to GT column
results['GT'] = results['GT'].apply(reverse_gt)

# Save de-pseudonymized results
results.to_csv('original_results.tsv', sep='\t', index=False)
```

## Technical Details

### What Gets Pseudonymized

- Sample IDs in the GT (genotype) column
- Sample IDs in inheritance analysis columns
- Sample IDs in PED files (if --pseudonymize-ped is used)
- Any future columns that contain sample identifiers

### Consistency Guarantees

- The same original ID always maps to the same pseudonym within a run
- Pseudonyms are deterministic (sorted before assignment)
- Duplicate handling ensures all pseudonyms are unique

### Performance Considerations

- Pseudonymization adds minimal overhead
- Processing happens in-memory before final output
- No impact on analysis accuracy or completeness

## Troubleshooting

### Common Issues

1. **Missing mapping table path**
   ```
   Error: --pseudonymize-table is required when using --pseudonymize
   ```
   Solution: Always specify where to save the mapping table

2. **Custom schema pattern errors**
   ```
   Warning: Missing placeholder in pattern: 'study'
   ```
   Solution: Ensure all placeholders in your pattern have corresponding metadata

3. **Mapping table in wrong location**
   - The mapping table is saved in the parent directory of your output
   - This is intentional for security
   - Check one directory level up from your results

### Validation Checklist

Before sharing pseudonymized data:
- [ ] Verify mapping table is NOT in the results directory
- [ ] Check that all sample IDs are replaced in output files
- [ ] Confirm PED file is pseudonymized (if using pedigrees)
- [ ] Test that results archive doesn't contain mapping table
- [ ] Document which schema was used for reproducibility

## See Also

- [Usage Guide](../usage.md) — CLI reference for all pseudonymization flags
- [Configuration Guide](../configuration.md) — Config file options
- [Rare Disease Workflow](../guides/rare_disease_workflow.md) — Family analysis with pseudonymization
