# Cancer Genomics Analysis

This guide covers somatic variant analysis workflows for cancer genomics using VariantCentrifuge, including tumor-normal paired analysis, tumor-only analysis, and loss of heterozygosity (LOH) detection.

## Overview

VariantCentrifuge supports several cancer analysis modes through built-in somatic presets and configurable tumor-normal parameters:

- **Tumor-normal paired** — Compare tumor and matched normal samples to identify somatic variants
- **Tumor-only** — Analyze tumor samples without a matched normal
- **LOH detection** — Identify loss of heterozygosity events
- **Germline shared** — Find variants present in both tumor and normal

## Tumor-Normal Paired Analysis

### Prerequisites

Your VCF must contain both tumor and normal samples, typically produced by a somatic variant caller such as Mutect2, Strelka2, or VarDict.

### Basic Workflow

```bash
variantcentrifuge \
  --gene-file oncogenes_tsg.txt \
  --vcf-file tumor_normal.vcf.gz \
  --preset somatic,coding \
  --tumor-sample-index 1 \
  --normal-sample-index 0 \
  --html-report \
  --output-file somatic_variants.tsv
```

### Configuring Sample Indices

VCF files list samples in a specific order. Use `--tumor-sample-index` and `--normal-sample-index` to identify which sample is which:

| Flag | Default | Description |
|------|---------|-------------|
| `--tumor-sample-index` | `1` | 0-based index of the tumor sample |
| `--normal-sample-index` | `0` | 0-based index of the normal sample |

To check sample order in your VCF:
```bash
bcftools query -l tumor_normal.vcf.gz
```

### Quality Thresholds

Fine-tune somatic variant calling quality with these parameters:

| Flag | Default | Description |
|------|---------|-------------|
| `--tumor-dp-min` | `20` | Minimum read depth in tumor |
| `--normal-dp-min` | `20` | Minimum read depth in normal |
| `--tumor-af-min` | `0.05` | Minimum allele frequency in tumor (5%) |
| `--normal-af-max` | `0.03` | Maximum allele frequency in normal (3%) |

These values are substituted into preset filter expressions via template variables (`{tumor_idx}`, `{normal_idx}`, `{tumor_dp_min}`, etc.).

### Strict Somatic Filtering

For high-confidence somatic calls with stringent thresholds:

```bash
variantcentrifuge \
  --gene-file cancer_panel.txt \
  --vcf-file tumor_normal.vcf.gz \
  --preset somatic_strict,coding \
  --tumor-sample-index 1 \
  --normal-sample-index 0 \
  --tumor-dp-min 30 \
  --tumor-af-min 0.10 \
  --normal-af-max 0.01 \
  --html-report \
  --xlsx \
  --output-file strict_somatic.tsv
```

## Tumor-Only Analysis

When a matched normal sample is not available, use the `tumor_only` preset which filters on population frequency databases to remove common germline variants:

```bash
variantcentrifuge \
  --gene-file oncogenes_tsg.txt \
  --vcf-file tumor_only.vcf.gz \
  --preset tumor_only,coding \
  --tumor-sample-index 0 \
  --tumor-dp-min 30 \
  --tumor-af-min 0.05 \
  --html-report \
  --output-file tumor_only_variants.tsv
```

:::{note}
Tumor-only analysis has a higher false-positive rate for somatic calls because rare germline variants may pass population frequency filters. Consider using stricter frequency thresholds.
:::

## Loss of Heterozygosity (LOH)

Detect LOH events where the normal sample is heterozygous but the tumor shows allelic imbalance:

```bash
variantcentrifuge \
  --gene-file tsg_panel.txt \
  --vcf-file tumor_normal.vcf.gz \
  --preset loh \
  --tumor-sample-index 1 \
  --normal-sample-index 0 \
  --output-file loh_events.tsv
```

## Available Somatic Presets

VariantCentrifuge ships with several somatic-focused presets:

| Preset | Description |
|--------|-------------|
| `somatic` | Standard somatic filter: tumor AF >= threshold, normal AF <= threshold, minimum depth |
| `somatic_pass` | Somatic filter restricted to PASS variants |
| `somatic_strict` | Stringent somatic filter with higher depth and AF requirements |
| `loh` | Loss of heterozygosity: normal is het (0/1), tumor shows imbalance |
| `germline_shared` | Variants present in both tumor and normal (shared germline) |
| `tumor_only` | Tumor-only mode: filters by population frequency without matched normal |
| `mutect2_TvsN_pass` | Mutect2-specific: PASS filter for tumor-normal pairs |
| `mutect2_TvsN` | Mutect2-specific: includes non-PASS variants |
| `mutect2_To_pass` | Mutect2-specific: PASS filter for tumor-only |

Combine with impact presets for targeted analysis: `--preset somatic,coding` or `--preset somatic,high`.

## Custom Somatic Configuration

For non-standard somatic calling pipelines, define custom presets in your config file:

```json
{
  "reference": "GRCh38.99",
  "presets": {
    "my_somatic": "(GEN[{normal_idx}].AF < {normal_af_max}) & (GEN[{tumor_idx}].AF >= {tumor_af_min}) & (GEN[{tumor_idx}].DP >= {tumor_dp_min}) & (GEN[{normal_idx}].DP >= {normal_dp_min})",
    "cosmic_hotspot": "((exists ID) & (ID =~ 'COS'))",
    "somatic_rare": "((dbNSFP_gnomAD_exomes_AC[0] <= 2) | (na dbNSFP_gnomAD_exomes_AC[0]))"
  }
}
```

Template variables (`{tumor_idx}`, `{normal_idx}`, `{tumor_dp_min}`, etc.) are expanded at runtime using the CLI flags.

## Worked Example: Comprehensive Cancer Panel

```bash
# Step 1: Somatic variant analysis with scoring
variantcentrifuge \
  --gene-file comprehensive_cancer_panel.txt \
  --vcf-file tumor_normal.vcf.gz \
  --preset somatic,coding \
  --tumor-sample-index 1 \
  --normal-sample-index 0 \
  --tumor-dp-min 30 \
  --tumor-af-min 0.05 \
  --scoring-config-path scoring/nephro_candidate_score \
  --final-filter 'IMPACT in ["HIGH", "MODERATE"]' \
  --html-report \
  --xlsx \
  --output-file cancer_panel_results.tsv
```

## Interpretation Guidelines

### Variant Prioritization for Cancer

1. **Tier 1 — Known oncogenic**: Variants in COSMIC hotspots or known driver mutations
2. **Tier 2 — Likely oncogenic**: HIGH impact variants in known cancer genes with low population frequency
3. **Tier 3 — Uncertain significance**: MODERATE impact variants requiring further evidence
4. **Tier 4 — Likely passenger**: Common variants or LOW impact changes

### Key Columns for Cancer Analysis

- **IMPACT** — Functional impact (HIGH, MODERATE, LOW, MODIFIER)
- **gnomAD AF** — Population allele frequency (lower = more likely somatic)
- **CADD** — Combined Annotation Dependent Depletion score
- **ClinVar** — Clinical significance from ClinVar database
- **Tumor AF** — Variant allele frequency in tumor (from genotype field)
- **Normal AF** — Variant allele frequency in normal (should be ~0 for true somatic)

## Best Practices

1. **Verify sample order** in your VCF before running — swapped indices silently produce wrong results
2. **Use PASS variants** (`somatic_pass` preset) for clinical reporting; include non-PASS for research
3. **Adjust depth thresholds** based on your sequencing coverage (WGS ~30x, WES ~100x, panel ~500x)
4. **Combine with population filters** — true somatic variants are absent from population databases
5. **Review LOH alongside somatic** — LOH in tumor suppressors is a common second hit
6. **Use `--bcftools-prefilter`** on large WGS VCFs to speed up analysis

## See Also

- [Configuration Guide](../configuration.md) — Preset definitions and custom config
- [Custom Filters](custom_filters.md) — Writing custom filter expressions
- [Variant Scoring](variant_scoring.md) — Configuring scoring models
- [Usage Guide](../usage.md) — Complete CLI reference
