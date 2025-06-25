# Rare Disease Analysis Workflow

This guide provides a comprehensive workflow for analyzing rare disease variants using VariantCentrifuge, from initial data preparation through final interpretation.

## Overview

Rare disease variant analysis requires careful filtering to identify potentially causative variants while minimizing false positives. This workflow focuses on:

- Identifying rare, high-impact variants
- Prioritizing variants in known disease genes
- Leveraging clinical databases for variant interpretation
- Generating comprehensive reports for clinical review

## Prerequisites

- Annotated VCF files (see [Annotation Strategies](annotation_strategies.md))
- Disease gene lists (e.g., OMIM, ClinGen, custom panels)
- Patient phenotype information
- Family information (if available)

## Step-by-Step Workflow

### Step 1: Prepare Gene Lists

```bash
# Create disease-specific gene lists
echo "BRCA1\nBRCA2\nTP53\nATM\nCHEK2" > breast_cancer_genes.txt

# Or use comprehensive panels
wget https://ftp.clinicalgenome.org/ClinGen_gene_curation_list.tsv
awk -F'\t' '$3 == "Breast cancer" {print $1}' ClinGen_gene_curation_list.tsv > clinical_breast_genes.txt
```

### Step 2: Configure Analysis

Create a rare disease-specific configuration:

```json
{
  "reference": "GRCh38.99",
  "filters": "",
  "fields_to_extract": "CHROM POS REF ALT ANN[0].GENE ANN[0].FEATUREID ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P dbNSFP_CADD_phred dbNSFP_REVEL_score gnomAD_exomes_AF gnomAD_genomes_AF ClinVar_CLNSIG ClinVar_CLNDN HGMD_CLASS GEN[*].GT GEN[*].DP",
  "interval_expand": 20,
  "perform_gene_burden": true,
  "presets": {
    "rare_pathogenic": "(((gnomAD_exomes_AF < 0.001) | (na gnomAD_exomes_AF)) & ((gnomAD_genomes_AF < 0.001) | (na gnomAD_genomes_AF))) & ((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE') | ((ClinVar_CLNSIG =~ '[Pp]athogenic') & !(ClinVar_CLNSIG =~ '[Cc]onflicting')))",
    "ultra_rare": "(((gnomAD_exomes_AC <= 2) | (na gnomAD_exomes_AC)) & ((gnomAD_genomes_AC <= 2) | (na gnomAD_genomes_AC)))",
    "high_confidence": "(GEN[*].DP >= 20) & (GEN[*].GQ >= 20) & (QUAL >= 30)",
    "protein_affecting": "((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE'))",
    "clinical_significance": "((ClinVar_CLNSIG =~ '[Pp]athogenic') | (ClinVar_CLNSIG =~ 'VUS') | (ClinVar_CLNSIG =~ '[Ll]ikely'))"
  }
}
```

### Step 3: Run Initial Analysis

```bash
# Single patient analysis
variantcentrifuge \
    --config rare_disease_config.json \
    --gene-file disease_genes.txt \
    --vcf-file patient_001.vcf.gz \
    --preset rare_pathogenic,high_confidence \
    --output-file patient_001_rare_variants.tsv \
    --html-report \
    --xlsx
```

### Step 4: Inheritance Pattern Analysis

For family-based analysis:

```bash
# Trio analysis (proband + parents)
variantcentrifuge \
    --config rare_disease_config.json \
    --gene-file disease_genes.txt \
    --vcf-file family_trio.vcf.gz \
    --preset ultra_rare,protein_affecting \
    --samples-file trio_samples.txt \
    --phenotype-file family_phenotypes.tsv \
    --phenotype-sample-column "sample_id" \
    --phenotype-value-column "affected_status" \
    --output-file trio_analysis.tsv \
    --perform-gene-burden \
    --html-report
```

Trio samples file (`trio_samples.txt`):
```
proband_001	Proband
father_001	Father  
mother_001	Mother
```

Family phenotypes file (`family_phenotypes.tsv`):
```
sample_id	affected_status	age	sex
phenotype_details
proband_001	affected	25	F	intellectual disability, seizures
father_001	unaffected	55	M	normal
mother_001	unaffected	52	F	normal
```

### Step 5: Multi-Gene Panel Analysis

```bash
# Comprehensive gene panel analysis
variantcentrifuge \
    --config rare_disease_config.json \
    --gene-file comprehensive_panel.txt \
    --vcf-file patient_001.vcf.gz \
    --preset clinical_significance \
    --filters "(GEN[*].DP >= 15) & (QUAL >= 20)" \
    --output-file patient_001_panel.tsv \
    --igv \
    --bam-mapping-file patient_bams.tsv \
    --igv-reference hg38 \
    --html-report
```

## Interpretation Guidelines

### Variant Prioritization

1. **Pathogenic/Likely Pathogenic in ClinVar**
   - Highest priority for known disease variants
   - Review evidence and conflicting interpretations

2. **Loss-of-Function in Disease Genes**
   - Nonsense, frameshift, splice site variants
   - High confidence if gene is known for haploinsufficiency

3. **Missense with Strong Predictions**
   - CADD > 20, REVEL > 0.7
   - Affecting conserved domains
   - Novel or ultra-rare (AC ≤ 2)

4. **Splice Region Variants**
   - Within 2bp of exon boundaries
   - Novel splice predictions
   - Experimental validation recommended

### Filtering Strategy

```bash
# Tier 1: Known pathogenic variants
variantcentrifuge \
    --preset clinical_significance \
    --filters "(ClinVar_CLNSIG =~ '[Pp]athogenic') & !(ClinVar_CLNSIG =~ '[Cc]onflicting')" \
    --gene-file disease_genes.txt \
    --vcf-file patient.vcf.gz \
    --output-file tier1_pathogenic.tsv

# Tier 2: High-impact rare variants
variantcentrifuge \
    --preset ultra_rare,protein_affecting \
    --filters "(ANN[ANY].IMPACT has 'HIGH') & (GEN[*].DP >= 20)" \
    --gene-file disease_genes.txt \
    --vcf-file patient.vcf.gz \
    --output-file tier2_high_impact.tsv

# Tier 3: Moderate impact with strong predictions
variantcentrifuge \
    --preset rare_pathogenic \
    --filters "(ANN[ANY].IMPACT has 'MODERATE') & ((dbNSFP_CADD_phred >= 25) | (dbNSFP_REVEL_score >= 0.7))" \
    --gene-file disease_genes.txt \
    --vcf-file patient.vcf.gz \
    --output-file tier3_moderate_predicted.tsv
```

## Case Studies

### Case 1: Intellectual Disability

**Patient:** 8-year-old with developmental delay and seizures

```bash
# Analysis focusing on neurodevelopmental genes
variantcentrifuge \
    --config rare_disease_config.json \
    --gene-file intellectual_disability_genes.txt \
    --vcf-file patient_ID001.vcf.gz \
    --preset ultra_rare,protein_affecting \
    --output-file ID001_neurodevelopmental.tsv \
    --html-report
```

**Key findings:**
- De novo nonsense variant in *SCN1A* (Dravet syndrome)
- Ultra-rare missense variant in *STXBP1* with CADD=28
- Multiple VUS requiring functional studies

### Case 2: Cardiomyopathy

**Patient:** 35-year-old with hypertrophic cardiomyopathy

```bash
# Cardiomyopathy gene panel analysis
variantcentrifuge \
    --config rare_disease_config.json \
    --gene-file cardiomyopathy_genes.txt \
    --vcf-file patient_CM002.vcf.gz \
    --preset clinical_significance \
    --filters "(ClinVar_CLNSIG =~ '[Pp]athogenic|VUS') | ((gnomAD_exomes_AF < 0.0001) & (ANN[ANY].IMPACT has 'HIGH|MODERATE'))" \
    --output-file CM002_cardiomyopathy.tsv \
    --igv \
    --html-report
```

**Key findings:**
- Pathogenic variant in *MYBPC3* (known HCM gene)
- Family segregation analysis recommended
- Cascade screening for relatives

## Quality Control

### Sample Quality Metrics

```bash
# Check coverage and quality metrics
variantcentrifuge \
    --preset high_confidence \
    --filters "(GEN[*].DP >= 20) & (GEN[*].GQ >= 20)" \
    --fields "CHROM POS REF ALT GEN[*].DP GEN[*].GQ GEN[*].AD" \
    --gene-file disease_genes.txt \
    --vcf-file patient.vcf.gz \
    --output-file quality_check.tsv
```

### Validation Requirements

1. **Sanger sequencing** for pathogenic/likely pathogenic variants
2. **Family segregation** when family samples available
3. **Functional studies** for novel VUS in critical genes
4. **CNV analysis** for genes with known deletion/duplication syndromes

## Reporting Templates

### Clinical Report Structure

1. **Patient Information**
   - Demographics and phenotype
   - Family history
   - Indication for testing

2. **Methods**
   - Sequencing platform and coverage
   - Analysis pipeline and filters
   - Gene panel composition

3. **Results**
   - Pathogenic/likely pathogenic variants
   - Variants of uncertain significance
   - Negative findings in relevant genes

4. **Interpretation**
   - Clinical significance
   - Inheritance pattern
   - Recommendations

### Automated Report Generation

```bash
# Generate clinical report with key findings
variantcentrifuge \
    --config clinical_report_config.json \
    --gene-file clinical_panel.txt \
    --vcf-file patient.vcf.gz \
    --preset clinical_significance \
    --phenotype-file patient_phenotype.tsv \
    --output-file clinical_report.tsv \
    --html-report \
    --xlsx
```

## Best Practices

1. **Use validated gene panels** for specific conditions
2. **Apply appropriate frequency thresholds** based on disease prevalence
3. **Consider inheritance patterns** in filtering strategy
4. **Validate critical findings** with orthogonal methods
5. **Document analysis parameters** for reproducibility
6. **Regular database updates** for clinical annotations
7. **Multidisciplinary review** of complex cases
8. **Genetic counseling** for positive findings

## Common Pitfalls

1. **Over-reliance on prediction tools** without experimental validation
2. **Ignoring population-specific frequencies** in diverse cohorts
3. **Inadequate coverage assessment** in critical gene regions
4. **Missing structural variants** in single nucleotide variant analysis
5. **Insufficient phenotype documentation** for variant interpretation

By following this workflow, you can systematically analyze rare disease variants and generate clinically actionable reports using VariantCentrifuge.