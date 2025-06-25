# Cohort Analysis Guide

This guide explains how to perform cohort-level variant analysis using VariantCentrifuge, including aggregating results from multiple samples and generating comprehensive cohort reports.

## Overview

VariantCentrifuge supports both single-sample and cohort-level analysis workflows. While the main CLI processes individual VCF files, the cohort analysis tools allow you to:

- Aggregate results from multiple single-sample analyses
- Identify recurrently mutated genes across the cohort  
- Generate interactive HTML reports for cohort-wide visualization
- Apply dynamic filtering and explore variant patterns
- Compute cohort-level statistics and gene burden analyses

## Cohort Analysis Workflow

### Step 1: Single-Sample Analysis

First, run VariantCentrifuge on each sample in your cohort:

```bash
# Example: Process multiple samples
for sample in APA12 APA13 APA15 APA18; do
    variantcentrifuge \
        --gene-file cancer_genes.txt \
        --vcf-file "${sample}_TvsN.filtered.annotated.vcf.gz" \
        --output-file "output/${sample}_variants.tsv" \
        --xlsx \
        --html-report
done
```

This creates individual analysis results for each sample that can later be aggregated.

### Step 2: Cohort Report Generation

Use the `create_cohort_report.py` script to aggregate individual sample results into a comprehensive cohort report.

#### Prerequisites

- Python 3 with pandas and Jinja2 installed
- Collection of VariantCentrifuge TSV output files
- Consistent analysis parameters across all samples

#### Usage

```bash
python scripts/create_cohort_report.py \
  --input-pattern 'output/**/*.vc_analysis.tsv' \
  --output-dir 'cohort_analysis' \
  --report-name 'Cancer Gene Cohort Analysis' \
  --sample-regex 'output/([^/]+)_variants'
```

#### Command-Line Arguments

- **`--input-pattern`** (Required): Glob pattern to find input TSV files
- **`--output-dir`** (Required): Directory for output files (created if needed)
- **`--report-name`** (Optional): Custom report title (default: "Variant Cohort Analysis")
- **`--sample-regex`** (Required): Regex to extract sample IDs from file paths

#### Sample ID Extraction

The `--sample-regex` parameter uses a Python regular expression with one capturing group to extract sample identifiers:

```bash
# For file: output/variantcentrifuge_analysis/APA12_TvsN.filtered.annotated/APA12_TvsN.filtered.annotated.vc_analysis.tsv
# Regex: 'variantcentrifuge_analysis/([^/]+)/'
# Captures: APA12_TvsN.filtered.annotated

# For file: results/sample_APA15/variants.tsv  
# Regex: 'sample_([^/]+)/'
# Captures: APA15
```

### Step 3: Interactive Report Features

The generated cohort report provides:

#### Dynamic Filtering
- **Text filters** for gene names, variant effects, and sample IDs
- **Range sliders** for numeric fields (allele frequency, CADD scores, etc.)
- **Multi-select filters** for categorical data (impact levels, clinical significance)
- **Real-time updates** of statistics and visualizations

#### Interactive Visualizations
- **Gene frequency bar chart** showing most commonly mutated genes
- **Clickable gene bars** that filter the main variant table
- **Summary statistics dashboard** that updates with applied filters
- **Exportable data** in CSV/Excel formats

#### High-Performance Table
- **Pagination** for large datasets
- **Column sorting** and visibility controls
- **Global search** across all fields
- **Export functionality** for filtered results

## Advanced Cohort Analysis

### Gene Burden Testing

For statistical gene burden analysis across the cohort:

```bash
# Run gene burden analysis for each sample
variantcentrifuge \
    --gene-file target_genes.txt \
    --vcf-file cohort_merged.vcf.gz \
    --perform-gene-burden \
    --phenotype-file sample_phenotypes.tsv \
    --phenotype-sample-column "sample_id" \
    --phenotype-value-column "case_control_status" \
    --output-file cohort_gene_burden.tsv
```

### Custom Cohort Configurations

Create cohort-specific configuration files:

```json
{
  "reference": "GRCh38.99",
  "filters": "",
  "fields_to_extract": "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P dbNSFP_gnomAD_exomes_AF ClinVar_CLNSIG GEN[*].GT",
  "presets": {
    "cohort_rare": "(((dbNSFP_gnomAD_exomes_AF < 0.001) | (na dbNSFP_gnomAD_exomes_AF)) & ((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE')))",
    "cohort_quality": "(GEN[*].DP >= 20) & (QUAL >= 30)"
  },
  "perform_gene_burden": true,
  "html_report_default_hidden_columns": [
    "QUAL", "FEATUREID", "AA_POS", "AA_LEN"
  ]
}
```

### Batch Processing with Snakemake

For large cohorts, use the included Snakemake workflow:

```bash
# Configure analysis in snakemake/config_vc.yaml
# Then run the workflow
cd snakemake
snakemake --cores 8 --use-conda
```

## Example Workflows

### Workflow 1: Cancer Cohort Analysis

```bash
# Step 1: Process tumor-normal pairs
for sample in TCGA_*; do
    variantcentrifuge \
        --config cancer_config.json \
        --gene-file oncogenes_tsg.txt \
        --vcf-file "${sample}_somatic.vcf.gz" \
        --preset mutect2_TvsN,coding \
        --output-file "results/${sample}_somatic_variants.tsv" \
        --html-report \
        --igv \
        --bam-mapping-file bam_files.tsv
done

# Step 2: Generate cohort report
python scripts/create_cohort_report.py \
    --input-pattern 'results/*_somatic_variants.tsv' \
    --output-dir 'tcga_cohort_report' \
    --report-name 'TCGA Somatic Variants - Oncogenes & TSGs' \
    --sample-regex 'results/([^_]+)_somatic'
```

### Workflow 2: Rare Disease Cohort

```bash
# Step 1: Process affected individuals
for sample in patient_*; do
    variantcentrifuge \
        --config rare_disease_config.json \
        --gene-file disease_genes.txt \
        --vcf-file "${sample}_filtered.vcf.gz" \
        --preset rare,coding,not_benign \
        --phenotype-file patient_phenotypes.tsv \
        --phenotype-sample-column "sample_id" \
        --phenotype-value-column "affected_status" \
        --output-file "results/${sample}_rare_variants.tsv" \
        --perform-gene-burden \
        --xlsx
done

# Step 2: Aggregate and analyze
python scripts/create_cohort_report.py \
    --input-pattern 'results/*_rare_variants.tsv' \
    --output-dir 'rare_disease_cohort' \
    --report-name 'Rare Disease Cohort Analysis' \
    --sample-regex 'results/(patient_[^_]+)_rare'
```

### Workflow 3: Population Genetics Study

```bash
# Step 1: Process population samples
for pop in EUR AFR EAS SAS; do
    for sample in ${pop}_*; do
        variantcentrifuge \
            --config population_config.json \
            --gene-file population_genes.txt \
            --vcf-file "${sample}.vcf.gz" \
            --preset 5percent,coding \
            --output-file "results/${pop}/${sample}_variants.tsv"
    done
done

# Step 2: Generate population-specific reports
for pop in EUR AFR EAS SAS; do
    python scripts/create_cohort_report.py \
        --input-pattern "results/${pop}/*_variants.tsv" \
        --output-dir "cohort_reports/${pop}" \
        --report-name "${pop} Population Variant Analysis" \
        --sample-regex "results/${pop}/([^_]+)_variants"
done
```

## Output Structure

The cohort analysis generates the following output structure:

```
cohort_analysis/
├── cohort_report.html          # Interactive HTML report
└── data/
    ├── variants.json          # Complete variant dataset
    ├── summary_stats.json     # Cohort-wide statistics
    └── gene_frequencies.json  # Gene mutation frequencies
```

### Cohort Report Features

#### Summary Dashboard
- Total variants and samples
- Most frequently mutated genes
- Impact distribution
- Quality metrics

#### Interactive Variant Table
- All variants from all samples
- Sample-specific annotations
- Filterable by any column
- Exportable in multiple formats

#### Gene-Level Analysis
- Mutation frequency per gene
- Sample distribution per gene
- Statistical significance testing
- Visual gene ranking

## Best Practices

### Data Consistency
1. **Use identical configurations** across all samples in a cohort
2. **Apply consistent quality filters** before VariantCentrifuge analysis
3. **Ensure sample naming consistency** for proper aggregation
4. **Validate annotation completeness** across all input files

### Performance Optimization
1. **Process samples in parallel** when possible
2. **Use appropriate hardware resources** for large cohorts
3. **Implement checkpointing** for long-running analyses
4. **Monitor disk space** for intermediate files

### Quality Control
1. **Validate sample identities** before aggregation
2. **Check for batch effects** across processing groups
3. **Compare individual vs cohort statistics** for consistency
4. **Review outlier samples** that may indicate processing issues

### Report Customization
1. **Customize report names** to reflect study context
2. **Hide irrelevant columns** in the configuration
3. **Adjust filtering ranges** based on your data distribution
4. **Document analysis parameters** for reproducibility

## Troubleshooting

### Common Issues

**Empty cohort report**
- Check file path patterns and regex
- Verify TSV file formats and headers
- Ensure sample ID extraction is working

**Performance issues with large cohorts**
- Consider splitting into smaller batches
- Increase system memory allocation
- Use pagination in the HTML report

**Inconsistent results across samples**
- Verify identical analysis configurations
- Check for annotation version differences
- Validate reference genome consistency

### Debugging Commands

```bash
# Test sample ID extraction
python -c "
import re
pattern = 'results/([^/]+)_variants'
test_path = 'results/sample_001_variants.tsv'
match = re.search(pattern, test_path)
print('Sample ID:', match.group(1) if match else 'No match')
"

# Validate TSV file structure
head -1 results/*_variants.tsv | grep -c "CHROM"  # Should equal number of files

# Check variant counts per sample
for file in results/*_variants.tsv; do
    echo "$(basename $file): $(tail -n +2 $file | wc -l) variants"
done
```

By following this guide, you can effectively perform cohort-level variant analysis and generate comprehensive, interactive reports for your research or clinical studies.