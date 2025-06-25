# Cancer Genomics Analysis

This guide covers somatic variant analysis workflows for cancer genomics using VariantCentrifuge.

## Tumor-Normal Analysis

### Configuration for Somatic Variants

```json
{
  "reference": "GRCh38.99",
  "presets": {
    "somatic_quality": "(GEN[0].AF < 0.03) & (GEN[1].AF >= 0.05) & (GEN[*].DP >= 50)",
    "cosmic_or_rare": "(((gnomAD_exomes_AC <= 2) | (na gnomAD_exomes_AC)) | (exists ID & ID =~ 'COS'))"
  }
}
```

### Tumor-Only Analysis

For cases without matched normal samples:

```bash
variantcentrifuge \
    --gene-file oncogenes_tsg.txt \
    --vcf-file tumor_only.vcf.gz \
    --preset super_rare,coding \
    --filters "(GEN[*].AF >= 0.05) & (GEN[*].DP >= 30)" \
    --output-file tumor_variants.tsv
```