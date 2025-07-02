# Usage Guide

## Basic Usage

The most basic command to run VariantCentrifuge:

```bash
variantcentrifuge \
  --gene-name BICC1 \
  --vcf-file path/to/your.vcf \
  --output-file output.tsv
```

## Command Line Options

### Required Arguments

- `--vcf-file` - Input VCF file (can be compressed with gzip)
- `--output-file` - Output TSV file path

### Gene Selection (choose one)

- `--gene-name GENE` - Single gene name
- `--gene-file GENES.TXT` - File containing multiple genes (one per line)

### Configuration

- `--config CONFIG_FILE` - Load custom parameters from JSON config file
- `--reference REFERENCE` - snpEff reference database (overrides config)
- `--filters "FILTER_EXPRESSION"` - Custom SnpSift filters (overrides config)
- `--fields "FIELD_LIST"` - Custom fields to extract (overrides config)

### Input/Output Options

- `--samples-file SAMPLES.TXT` - Sample ID mapping for genotype replacement
- `--phenotype-file PHENO.TSV` - Phenotype data file
- `--phenotype-sample-column` - Column name for sample IDs in phenotype file
- `--phenotype-value-column` - Column name for phenotype values
- `--xlsx` - Convert final output TSV to XLSX format
- `--keep-intermediates` - Retain intermediate files after successful run

### Analysis Options

- `--perform-gene-burden` - Run gene burden analysis
- `--html-report` - Generate interactive HTML report
- `--igv` - Enable IGV.js integration (requires additional options)
- `--bam-mapping-file` - TSV/CSV file mapping sample IDs to BAM files
- `--igv-reference` - Genome reference for IGV (e.g., 'hg19', 'hg38')

### Scoring Options

- `--scoring-config-path` - Path to scoring configuration directory containing variable_assignment_config.json and formula_config.json

### Annotation Options

- `--annotate-bed BED_FILE` - Annotate variants with genomic regions from BED files (can specify multiple)
- `--annotate-gene-list GENE_LIST` - Check if variants affect genes in custom gene lists (can specify multiple)
- `--annotate-json-genes JSON_FILE` - Annotate variants with gene information from JSON file
- `--json-gene-mapping MAPPING` - Specify JSON field mapping for gene annotations (required with --annotate-json-genes)
- `--json-genes-as-columns` - Output JSON gene data as separate columns instead of appending to Custom_Annotation column

### Other Options

- `--version` - Show version and exit
- `--help` - Show help message

## Examples

### Basic Gene Analysis

```bash
variantcentrifuge \
  --gene-name BRCA1 \
  --vcf-file samples.vcf.gz \
  --output-file brca1_variants.tsv
```

### Multiple Genes with Custom Filters

```bash
variantcentrifuge \
  --gene-file cancer_genes.txt \
  --vcf-file samples.vcf.gz \
  --filters "(( dbNSFP_gnomAD_exomes_AC[0] <= 2 ) | ( na dbNSFP_gnomAD_exomes_AC[0] )) & ((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE'))" \
  --output-file cancer_variants.tsv \
  --xlsx
```

### Comprehensive Analysis with Reports

```bash
variantcentrifuge \
  --gene-name BRCA1 \
  --vcf-file samples.vcf.gz \
  --samples-file sample_mapping.txt \
  --phenotype-file patient_data.tsv \
  --phenotype-sample-column "sample_id" \
  --phenotype-value-column "disease_status" \
  --perform-gene-burden \
  --html-report \
  --xlsx \
  --output-file brca1_analysis.tsv
```

### IGV Integration

```bash
variantcentrifuge \
  --gene-name TP53 \
  --vcf-file samples.vcf.gz \
  --igv \
  --bam-mapping-file bam_files.tsv \
  --igv-reference hg38 \
  --html-report \
  --output-file tp53_variants.tsv
```

### Variant Scoring

```bash
# Apply custom scoring model to variants
variantcentrifuge \
  --gene-file kidney_genes.txt \
  --vcf-file patient.vcf.gz \
  --scoring-config-path scoring/nephro_variant_score \
  --preset rare,coding \
  --html-report \
  --output-file scored_variants.tsv
```

### Custom Annotations

```bash
# Annotate with JSON gene information
variantcentrifuge \
  --gene-name BRCA1 \
  --vcf-file samples.vcf.gz \
  --annotate-json-genes gene_metadata.json \
  --json-gene-mapping '{"identifier":"gene_symbol","dataFields":["panel","inheritance","function"]}' \
  --output-file annotated_variants.tsv

# Multiple annotation sources
variantcentrifuge \
  --gene-file cancer_genes.txt \
  --vcf-file samples.vcf.gz \
  --annotate-bed cancer_hotspots.bed \
  --annotate-gene-list actionable_genes.txt \
  --annotate-json-genes gene_panels.json \
  --json-gene-mapping '{"identifier":"symbol","dataFields":["panel_name","evidence_level"]}' \
  --html-report \
  --output-file multi_annotated.tsv
```

## Input File Formats

### VCF Files

- Standard VCF format (v4.0 or later)
- Can be compressed with gzip (.vcf.gz)
- Should be annotated with snpEff for optimal functionality

### Gene Files

Text file with one gene name per line:

```
BRCA1
BRCA2
TP53
ATM
```

### Sample Mapping Files

Tab-separated file for genotype replacement:

```
original_id	new_id
sample_001	Patient_A
sample_002	Patient_B
sample_003	Control_001
```

### Phenotype Files

Tab or comma-separated file with sample information:

```
sample_id	disease_status	age	sex
Patient_A	case	45	F
Patient_B	case	52	M
Control_001	control	48	F
```

### BAM Mapping Files

For IGV integration, provide a mapping from sample IDs to BAM file paths:

```
sample_id	bam_path
Patient_A	/path/to/patient_a.bam
Patient_B	/path/to/patient_b.bam
Control_001	/path/to/control_001.bam
```

### JSON Gene Files

For gene annotation, provide a JSON file containing an array of gene objects:

```json
[
  {
    "gene_symbol": "BRCA1",
    "panel": "HereditaryCancer",
    "inheritance": "AD",
    "function": "DNA repair"
  },
  {
    "gene_symbol": "TP53",
    "panel": "HereditaryCancer",
    "inheritance": "AD",
    "function": "Tumor suppressor"
  }
]
```

The `--json-gene-mapping` parameter specifies:
- `identifier`: The field containing the gene symbol (e.g., "gene_symbol")
- `dataFields`: Array of fields to include as annotations (e.g., ["panel", "inheritance", "function"])

## Output Files

### Main Output

- **TSV file** - Tab-separated variant table with extracted fields
- **XLSX file** - Excel format (if `--xlsx` specified)
- **Metadata file** - Analysis parameters and tool versions

### Optional Outputs

- **HTML report** - Interactive variant browser (if `--html-report` specified)
- **IGV reports** - Individual variant visualization (if `--igv` specified)
- **Gene burden results** - Statistical analysis (if `--perform-gene-burden` specified)

## Configuration

See the [Configuration Guide](configuration.md) for detailed information about setting up configuration files and customizing VariantCentrifuge behavior.

## Troubleshooting

### Common Issues

1. **No variants found:**
   - Check that your VCF file contains variants in the specified gene regions
   - Verify gene names are correct and match your reference annotation
   - Review filter expressions - they may be too restrictive

2. **External tool errors:**
   - Ensure all required tools are installed and in PATH
   - Check that snpEff database matches your VCF reference
   - Verify file permissions and disk space

3. **Memory issues:**
   - Large VCF files may require more memory
   - Consider filtering your VCF file beforehand to reduce size
   - Use `--keep-intermediates` to debug intermediate file sizes

### Getting Help

- Use `variantcentrifuge --help` for command-line options
- Check the [API Reference](api/index.md) for detailed function documentation
- Report issues on [GitHub](https://github.com/scholl-lab/variantcentrifuge/issues)