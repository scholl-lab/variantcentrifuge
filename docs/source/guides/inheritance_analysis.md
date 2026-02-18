# Inheritance Analysis

VariantCentrifuge includes a comprehensive Mendelian inheritance analysis system for identifying disease-causing variants in family studies. It supports de novo, autosomal dominant, autosomal recessive, X-linked, compound heterozygous, and mitochondrial inheritance patterns.

## Overview

The inheritance analysis pipeline uses a three-pass approach:

1. **Deduction** — Determine the most likely inheritance pattern for each variant based on genotypes and family structure
2. **Compound heterozygous detection** — Identify pairs of heterozygous variants in the same gene inherited from different parents (trans configuration)
3. **Prioritization** — Rank patterns by clinical significance when multiple patterns are possible

## Quick Start

### Family Trio Analysis

```bash
variantcentrifuge \
  --gene-file disease_genes.txt \
  --vcf-file trio.vcf.gz \
  --ped family.ped \
  --inheritance-mode columns \
  --preset rare,coding \
  --html-report \
  --output-file trio_results.tsv
```

### Singleton Analysis (No PED File)

Without a PED file, all samples are treated as affected singletons:

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file proband.vcf.gz \
  --inheritance-mode simple \
  --preset rare,coding \
  --output-file singleton_results.tsv
```

## PED File Format

VariantCentrifuge uses the standard PLINK PED format (tab-separated, 6 columns):

```
#Family_ID  Individual_ID  Father_ID  Mother_ID  Sex  Affected
FAM001      proband        father     mother     1    2
FAM001      father         0          0          1    1
FAM001      mother         0          0          2    1
```

| Column | Values |
|--------|--------|
| Family_ID | Family identifier (groups related individuals) |
| Individual_ID | Must match sample IDs in VCF |
| Father_ID | Father's Individual_ID, or `0` if unknown/unavailable |
| Mother_ID | Mother's Individual_ID, or `0` if unknown/unavailable |
| Sex | `1` = male, `2` = female, `0` = unknown |
| Affected | `2` = affected, `1` = unaffected, `0` = unknown |

:::{tip}
Sample IDs in the PED file must exactly match the sample names in the VCF file. Use `bcftools query -l your.vcf.gz` to check sample names.
:::

## Inheritance Patterns

### Supported Patterns

| Pattern | Description | Typical Genotypes |
|---------|-------------|-------------------|
| `de_novo` | Not present in either parent | Proband: 0/1, Parents: 0/0 |
| `autosomal_dominant` (AD) | Present in one affected parent | Proband: 0/1, Affected parent: 0/1, Unaffected parent: 0/0 |
| `autosomal_recessive` (AR) | Both parents are carriers | Proband: 1/1, Parents: 0/1 |
| `x_linked_recessive` (XLR) | X-linked, male affected | Male proband: 1 (hemizygous), Mother: 0/1 (carrier) |
| `x_linked_dominant` (XLD) | X-linked, both sexes affected | Proband: 0/1 or 1/1 on chrX |
| `compound_heterozygous` | Two het variants in same gene from different parents | Proband: 0/1 at two loci, each parent carries one |
| `mitochondrial` | Mitochondrial (chrM) inheritance | Proband: variant on chrM |
| `unknown` | Pattern cannot be determined | Missing data or ambiguous genotypes |

### Pattern Prioritization

When multiple patterns are possible, they are ranked by clinical significance:

1. de_novo (highest priority)
2. compound_heterozygous
3. autosomal_recessive / x_linked_recessive
4. x_linked_dominant
5. autosomal_dominant
6. unknown (lowest priority)

## Output Modes

The `--inheritance-mode` flag controls how inheritance results appear in the output:

### Simple Mode (`--inheritance-mode simple`)

Adds a single `Inheritance_Pattern` column with the pattern name:

```
GENE    POS     Inheritance_Pattern
BRCA1   41276   de_novo
PKD1    2138    autosomal_recessive
```

### Columns Mode (`--inheritance-mode columns`)

Expands inheritance results into separate columns for easier filtering:

```
GENE    POS     Inheritance_Pattern    Inheritance_Confidence    Inheritance_Samples    Inheritance_Details
BRCA1   41276   de_novo                0.95                      proband(0/1)           Father: 0/0, Mother: 0/0
PKD1    2138    autosomal_recessive    0.90                      proband(1/1)           Father: 0/1, Mother: 0/1
```

### Full Mode (`--inheritance-mode full`)

Outputs the complete inheritance analysis as a JSON object in a single column. Useful for programmatic downstream analysis.

## Compound Heterozygous Detection

Compound heterozygous variants are two heterozygous variants in the same gene inherited from different parents (one from each parent — trans configuration).

### How It Works

1. For each gene, identify samples with 2+ heterozygous variants
2. Check parental genotypes to confirm trans configuration (one variant from father, one from mother)
3. When parents are unavailable, variants are marked as `compound_het_possible`
4. Multiple compound het pairs in a gene are ranked by variant impact and allele frequency

### Performance

The default vectorized implementation is 10-50x faster than the original for genes with many variants:

| Variants per Gene | Original | Vectorized | Speedup |
|-------------------|----------|------------|---------|
| 10 | 0.01s | 0.001s | 10x |
| 50 | 0.19s | 0.004s | 48x |
| 100 | 1.23s | 0.017s | 72x |
| 500 | 31.0s | 0.55s | 56x |

To use the original implementation (for debugging or comparison):
```bash
--no-vectorized-comp-het
```

## Integration with Scoring

Inheritance patterns integrate with the variant scoring system. The built-in `inheritance_score` model assigns scores based on clinical significance:

| Pattern | Score |
|---------|-------|
| de_novo | 0.95 |
| autosomal_recessive / compound_het / homozygous | 0.80 |
| x_linked_recessive | 0.70 |
| x_linked_dominant | 0.50 |
| autosomal_dominant | 0.40 |
| compound_het_possible | 0.40 |
| unknown | 0.10 |

Use inheritance in scoring formulas:

```bash
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file trio.vcf.gz \
  --ped family.ped \
  --inheritance-mode columns \
  --scoring-config-path scoring/nephro_candidate_score \
  --output-file scored_trio.tsv
```

## Filtering on Inheritance

Use `--final-filter` to filter results by inheritance pattern:

```bash
# Only de novo and compound het variants
--final-filter 'Inheritance_Pattern in ["de_novo", "compound_heterozygous"]'

# High-confidence inheritance calls
--final-filter 'Inheritance_Confidence > 0.8'

# Recessive patterns only
--final-filter 'Inheritance_Pattern in ["autosomal_recessive", "compound_heterozygous", "x_linked_recessive"]'
```

## Memory Management

For large multi-sample VCFs, inheritance analysis memory usage is automatically managed:

- System memory is detected (supports SLURM, PBS, cgroups, and bare metal)
- Chunk sizes are calculated based on available memory
- Use `--max-memory-gb` to set explicit limits
- Use `--force-inheritance-processing` to override memory safety checks
- `--memory-safety-factor` (default 0.92) controls how much of detected memory to use

## Worked Example: Rare Disease Trio

```bash
# Full trio analysis with inheritance, scoring, and reports
variantcentrifuge \
  --gene-file intellectual_disability_genes.txt \
  --vcf-file trio_annotated.vcf.gz \
  --ped trio.ped \
  --inheritance-mode columns \
  --preset rare,coding \
  --scoring-config-path scoring/nephro_candidate_score \
  --final-filter 'Inheritance_Pattern != "unknown" and IMPACT in ["HIGH", "MODERATE"]' \
  --html-report \
  --xlsx \
  --output-file trio_id_analysis.tsv
```

## Troubleshooting

### No inheritance patterns detected
- Verify PED file sample IDs match VCF sample names exactly
- Check that the PED file has correct family relationships
- Ensure the VCF contains genotype data (GT field) for all family members

### All patterns show as "unknown"
- Missing genotype data (./.) for key family members
- Parent samples not present in VCF
- Incorrect sex assignments (affects X-linked analysis)

### Compound het not detected
- Need at least 2 heterozygous variants per gene per affected sample
- Parent genotypes required to confirm trans configuration (without parents, `compound_het_possible` is reported)
- Check that gene names are consistent between variants

## See Also

- [Rare Disease Workflow](rare_disease_workflow.md) — End-to-end rare disease analysis
- [Variant Scoring](variant_scoring.md) — Integrating inheritance with scoring models
- [Performance Tips](performance_tips.md) — Memory management for large datasets
- [Usage Guide](../usage.md) — Complete CLI reference
