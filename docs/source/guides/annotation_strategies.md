# Annotation Strategies

This guide provides recommended strategies for annotating your VCF files before using VariantCentrifuge. Proper annotation is crucial for effective variant filtering and analysis.

## Overview

VariantCentrifuge relies on annotated VCF fields to filter and analyze variants. The quality and comprehensiveness of your annotations directly impact the effectiveness of your analysis. This guide covers:

- Recommended annotation tools and databases
- Step-by-step annotation workflows
- Database selection guidelines
- Quality control and validation
- Common pitfalls and solutions

## Recommended Tools

### Core Annotation Tools

#### SnpEff
**Purpose:** Functional effect prediction and gene annotation

**Key Features:**
- Predicts variant effects (HIGH, MODERATE, LOW, MODIFIER impact)
- Provides gene, transcript, and protein annotations
- Supports multiple genome builds and species
- Generates comprehensive ANN fields

**Installation:**
```bash
# Via conda/mamba
mamba install -c bioconda snpeff

# Download databases
snpEff download GRCh38.99  # or your preferred version
```

**Basic Usage:**
```bash
java -jar snpEff.jar GRCh38.99 input.vcf > annotated.vcf
```

#### SnpSift
**Purpose:** Database annotation and field extraction

**Key Features:**
- Adds population frequency data (gnomAD, ExAC, 1000G)
- Clinical significance annotations (ClinVar)
- Pathogenicity predictions (CADD, REVEL, SIFT, PolyPhen)
- Flexible field extraction and filtering

**Installation:**
```bash
# Usually bundled with SnpEff
mamba install -c bioconda snpsift
```

#### bcftools
**Purpose:** VCF manipulation and annotation

**Key Features:**
- High-performance VCF processing
- Custom annotation from BED/VCF files
- Header manipulation and normalization
- Variant normalization and filtering

**Installation:**
```bash
mamba install -c bioconda bcftools
```

### Database Selection

#### Population Frequency Databases

**gnomAD (Genome Aggregation Database)**
- **Current version:** v4.1 (recommended)
- **Coverage:** 125,748 exomes, 15,708 genomes
- **Best for:** General population frequency filtering
- **Fields:** `gnomAD_exomes_AF`, `gnomAD_genomes_AF`, `gnomAD_*_AC`

**ExAC (Exome Aggregation Consortium)**
- **Status:** Legacy (replaced by gnomAD)
- **Use only if:** gnomAD unavailable for your analysis

**1000 Genomes Project**
- **Best for:** Population-specific frequency data
- **Limitation:** Smaller sample size than gnomAD

#### Clinical Databases

**ClinVar**
- **Purpose:** Clinical significance annotations
- **Updated:** Monthly
- **Fields:** `ClinVar_CLNSIG`, `ClinVar_CLNDN`
- **Critical for:** Clinical variant interpretation

**HGMD (Human Gene Mutation Database)**
- **Purpose:** Known disease mutations
- **License:** Commercial (free version available)
- **Best for:** Known pathogenic variant identification

#### Pathogenicity Prediction

**CADD (Combined Annotation Dependent Depletion)**
- **Range:** 0-99 (higher = more deleterious)
- **Threshold:** ≥20 (top 1% most deleterious)
- **Best for:** Genome-wide pathogenicity scoring

**REVEL (Rare Exome Variant Ensemble Learner)**
- **Range:** 0-1 (higher = more pathogenic)
- **Threshold:** ≥0.5 for likely pathogenic
- **Best for:** Missense variant prediction

**dbNSFP**
- **Contains:** Multiple prediction scores (SIFT, PolyPhen, CADD, REVEL, etc.)
- **Advantage:** Comprehensive collection
- **Best for:** Comparative scoring

## End-to-End Annotation Workflows

### Workflow 1: Comprehensive Research Annotation

This workflow provides comprehensive annotation suitable for most research applications.

#### Step 1: Variant Normalization
```bash
# Normalize variants (split multi-allelic, left-align)
bcftools norm -m-both -f reference.fasta input.vcf.gz > normalized.vcf
```

#### Step 2: Functional Annotation with SnpEff
```bash
# Annotate with SnpEff
java -Xmx8g -jar snpEff.jar GRCh38.99 normalized.vcf > snpeff_annotated.vcf

# Alternative with custom options
java -Xmx8g -jar snpEff.jar -v -stats snpeff_summary.html GRCh38.99 normalized.vcf > snpeff_annotated.vcf
```

#### Step 3: Database Annotation with SnpSift
```bash
# Add gnomAD frequencies
java -jar SnpSift.jar annotate gnomad.exomes.r4.0.sites.vcf.gz snpeff_annotated.vcf > gnomad_annotated.vcf

# Add ClinVar annotations
java -jar SnpSift.jar annotate clinvar.vcf.gz gnomad_annotated.vcf > clinvar_annotated.vcf

# Add dbNSFP predictions
java -jar SnpSift.jar dbnsfp -db dbNSFP4.1a.txt.gz clinvar_annotated.vcf > final_annotated.vcf
```

#### Step 4: Quality Control
```bash
# Check annotation completeness
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t%INFO/gnomAD_exomes_AF\n' final_annotated.vcf | head -10

# Generate annotation statistics
java -jar SnpSift.jar extractFields final_annotated.vcf CHROM POS REF ALT "ANN[0].IMPACT" "gnomAD_exomes_AF" > annotation_check.tsv
```

### Workflow 2: Clinical Diagnostic Annotation

Optimized for clinical variant interpretation with emphasis on known pathogenic variants.

#### Step 1: Core Annotation
```bash
# SnpEff with clinical focus
java -Xmx8g -jar snpEff.jar -canon GRCh38.99 input.vcf > clinical_snpeff.vcf
```

#### Step 2: Clinical Database Priority
```bash
# ClinVar (highest priority for clinical interpretation)
java -jar SnpSift.jar annotate -info "CLNSIG,CLNDN,CLNREVSTAT" clinvar.vcf.gz clinical_snpeff.vcf > clinical_clinvar.vcf

# HGMD (if available)
java -jar SnpSift.jar annotate hgmd.vcf.gz clinical_clinvar.vcf > clinical_hgmd.vcf

# gnomAD for population frequency
java -jar SnpSift.jar annotate -info "AF,AF_popmax" gnomad.exomes.vcf.gz clinical_hgmd.vcf > clinical_annotated.vcf
```

#### Step 3: Pathogenicity Scores
```bash
# Add CADD and REVEL via dbNSFP
java -jar SnpSift.jar dbnsfp -f CADD_phred,REVEL_score,SIFT_score,Polyphen2_HDIV_score -db dbNSFP4.1a.txt.gz clinical_annotated.vcf > clinical_final.vcf
```

### Workflow 3: Cancer Genomics Annotation

Focused on somatic variant annotation with cancer-relevant databases.

#### Step 1: Somatic-Specific Annotation
```bash
# SnpEff with all transcripts (cancer analysis often needs multiple isoforms)
java -Xmx8g -jar snpEff.jar -canon -v GRCh38.99 somatic.vcf > cancer_snpeff.vcf
```

#### Step 2: Cancer Databases
```bash
# COSMIC (Catalogue of Somatic Mutations in Cancer)
java -jar SnpSift.jar annotate cosmic.vcf.gz cancer_snpeff.vcf > cancer_cosmic.vcf

# OncoKB annotations (if available)
# Custom annotation scripts may be needed

# Population frequencies (for germline contamination assessment)
java -jar SnpSift.jar annotate gnomad.genomes.vcf.gz cancer_cosmic.vcf > cancer_annotated.vcf
```

## VariantCentrifuge Configuration Examples

### Example 1: Rare Disease Analysis

```json
{
  "reference": "GRCh38.99",
  "filters": "",
  "fields_to_extract": "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P dbNSFP_CADD_phred dbNSFP_REVEL_score gnomAD_exomes_AF gnomAD_genomes_AF ClinVar_CLNSIG GEN[*].GT",
  "presets": {
    "rare_pathogenic": "(((gnomAD_exomes_AF < 0.001) | (na gnomAD_exomes_AF)) & ((gnomAD_genomes_AF < 0.001) | (na gnomAD_genomes_AF))) & ((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE') | (ClinVar_CLNSIG =~ '[Pp]athogenic'))"
  }
}
```

### Example 2: Population Genetics Study

```json
{
  "reference": "GRCh38.99",
  "fields_to_extract": "CHROM POS REF ALT ID ANN[0].GENE ANN[0].IMPACT gnomAD_exomes_AF gnomAD_exomes_AC gnomAD_genomes_AF gnomAD_genomes_AC GEN[*].GT GEN[*].DP",
  "presets": {
    "coding_variants": "((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE'))",
    "high_quality": "(GEN[*].DP >= 20) & (QUAL >= 30)"
  }
}
```

### Example 3: Clinical Exome Analysis

```json
{
  "reference": "GRCh38.99",
  "fields_to_extract": "CHROM POS REF ALT ANN[0].GENE ANN[0].FEATUREID ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P dbNSFP_CADD_phred dbNSFP_REVEL_score gnomAD_exomes_AF ClinVar_CLNSIG ClinVar_CLNDN GEN[*].GT",
  "presets": {
    "clinical_significance": "(ClinVar_CLNSIG =~ '[Pp]athogenic') | (ClinVar_CLNSIG =~ 'VUS')",
    "rare_coding": "(((gnomAD_exomes_AF < 0.01) | (na gnomAD_exomes_AF)) & ((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE')))"
  }
}
```

## Quality Control and Validation

### Annotation Completeness Check

```bash
# Check annotation rates
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\t%INFO/gnomAD_exomes_AF\n' annotated.vcf | \
awk 'BEGIN{total=0; ann=0; gnomad=0}
     {total++; if($3!="") ann++; if($4!="") gnomad++}
     END{print "Total variants:", total; print "SnpEff annotated:", ann, "(" ann/total*100 "%)"; print "gnomAD annotated:", gnomad, "(" gnomad/total*100 "%)"}'
```

### Validate Database Versions

```bash
# Check database versions in VCF headers
bcftools view -h annotated.vcf | grep "##INFO.*database"
```

### Common Issues and Solutions

#### Issue: Low annotation rates
**Symptoms:** Many variants lack annotations
**Causes:** Reference genome mismatch, outdated databases
**Solutions:**
- Verify genome build consistency (GRCh37 vs GRCh38)
- Update annotation databases
- Check chromosome naming (chr1 vs 1)

#### Issue: Inconsistent population frequencies
**Symptoms:** Unexpected frequency distributions
**Causes:** Population stratification, database version differences
**Solutions:**
- Use population-specific frequency data
- Validate against known common variants
- Check for ancestry-specific databases

#### Issue: Missing pathogenicity predictions
**Symptoms:** Empty CADD/REVEL scores
**Causes:** Variant type limitations, database coverage
**Solutions:**
- Use multiple prediction tools
- Consider variant type-specific predictors
- Manual curation for critical variants

## Performance Optimization

### Memory Management
```bash
# Increase Java heap size for large VCFs
java -Xmx16g -jar snpEff.jar ...

# Use streaming processing for very large files
bcftools annotate -a database.vcf.gz ...
```

### Parallel Processing
```bash
# Split VCF by chromosome for parallel annotation
for chr in {1..22} X Y; do
    bcftools view -r chr${chr} input.vcf.gz | \
    java -jar snpEff.jar GRCh38.99 /dev/stdin > chr${chr}_annotated.vcf &
done
wait

# Merge results
bcftools concat chr*_annotated.vcf | bcftools sort > final_annotated.vcf.gz
```

### Database Preparation
```bash
# Create tabix indexes for faster annotation
bgzip database.vcf
tabix -p vcf database.vcf.gz

# Prepare custom annotation files
bcftools annotate --check-ref e -a custom_annotations.vcf.gz ...
```

## Custom Gene Annotations with VariantCentrifuge

In addition to standard VCF annotations, VariantCentrifuge provides built-in functionality to add custom annotations during analysis. These annotations are applied after variant extraction and are included in the final output.

### JSON Gene Annotations

The `--annotate-json-genes` feature allows you to integrate structured gene metadata from JSON files directly into your variant analysis.

#### JSON File Format

Create a JSON file containing an array of gene objects:

```json
[
  {
    "gene_symbol": "BRCA1",
    "panel": "HereditaryCancer",
    "inheritance": "AD",
    "function": "DNA repair",
    "disease_association": "Breast/Ovarian cancer",
    "actionability": "Tier1"
  },
  {
    "gene_symbol": "TP53",
    "panel": "HereditaryCancer",
    "inheritance": "AD",
    "function": "Tumor suppressor",
    "disease_association": "Li-Fraumeni syndrome",
    "actionability": "Tier1"
  },
  {
    "gene_symbol": "MLH1",
    "panel": "Lynch",
    "inheritance": "AD",
    "function": "DNA mismatch repair",
    "disease_association": "Lynch syndrome",
    "actionability": "Tier1"
  }
]
```

#### Field Mapping Configuration

Use the `--json-gene-mapping` parameter to specify how JSON fields map to annotations:

```json
{
  "identifier": "gene_symbol",
  "dataFields": ["panel", "inheritance", "actionability"]
}
```

- `identifier`: The JSON field containing the gene symbol
- `dataFields`: Array of fields to include in the Custom_Annotation column

#### Usage Example

```bash
variantcentrifuge \
  --gene-file cancer_genes.txt \
  --vcf-file patient.vcf.gz \
  --annotate-json-genes gene_metadata.json \
  --json-gene-mapping '{"identifier":"gene_symbol","dataFields":["panel","inheritance","actionability"]}' \
  --output-file annotated_variants.tsv
```

This will add annotations like `panel=HereditaryCancer;inheritance=AD;actionability=Tier1` to the Custom_Annotation column for variants in matching genes.

#### Adding JSON Annotations as Separate Columns

For a more structured output suitable for direct analysis, you can add each JSON data field as its own column instead of bundling them into the `Custom_Annotation` field. This is enabled with the `--json-genes-as-columns` flag.

**Usage Example:**

```bash
variantcentrifuge \
  --gene-name GENE \
  --vcf-file input.vcf.gz \
  --annotate-json-genes gene_data.json \
  --json-gene-mapping '{"identifier":"symbol","dataFields":["ngs","actionability"]}' \
  --json-genes-as-columns \
  --output-file variants_with_columns.tsv
```

This command will produce a TSV with two new columns, `ngs` and `actionability`, populated with the corresponding values from `gene_data.json`.

### BED File Annotations

Annotate variants with genomic regions using BED files:

```bash
variantcentrifuge \
  --gene-name GENE \
  --vcf-file input.vcf.gz \
  --annotate-bed hotspots.bed \
  --annotate-bed regulatory_regions.bed \
  --output-file output.tsv
```

### Gene List Annotations

Check if variants affect genes in custom lists:

```bash
variantcentrifuge \
  --gene-file all_genes.txt \
  --vcf-file input.vcf.gz \
  --annotate-gene-list actionable_genes.txt \
  --annotate-gene-list drug_targets.txt \
  --output-file output.tsv
```

### Combined Annotation Strategy

For comprehensive analysis, combine multiple annotation sources:

```bash
variantcentrifuge \
  --gene-file disease_genes.txt \
  --vcf-file patient.vcf.gz \
  --annotate-bed known_hotspots.bed \
  --annotate-gene-list clinically_actionable.txt \
  --annotate-json-genes gene_database.json \
  --json-gene-mapping '{"identifier":"symbol","dataFields":["panel","evidence","notes"]}' \
  --preset rare,coding \
  --html-report \
  --output-file comprehensive_analysis.tsv
```

### Integration with Scoring

Custom annotations can be used in variant scoring formulas. For example, if you annotate with `actionability` levels, you can create scoring formulas that prioritize Tier1 actionable variants.

## Best Practices Summary

1. **Plan your annotation strategy** based on your analysis goals
2. **Use the latest database versions** when possible
3. **Validate annotation completeness** before analysis
4. **Document your annotation workflow** for reproducibility
5. **Test on small datasets** before processing large cohorts
6. **Keep database versions consistent** across related analyses
7. **Consider population-specific databases** for diverse cohorts
8. **Backup original VCF files** before annotation
9. **Monitor resource usage** for large-scale annotations
10. **Validate critical variants manually** when needed
11. **Use custom annotations** to integrate project-specific gene metadata
12. **Combine annotation sources** for comprehensive variant characterization

By following these annotation strategies, you'll ensure that VariantCentrifuge has access to high-quality, comprehensive variant annotations for effective filtering and analysis.
