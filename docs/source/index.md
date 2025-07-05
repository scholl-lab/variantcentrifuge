# VariantCentrifuge: Clinical Variant Analysis Tool for VCF Files

**VariantCentrifuge** is a powerful Python-based command-line tool for filtering, analyzing, and interpreting genetic variants from VCF (Variant Call Format) files. Designed for clinical geneticists, bioinformaticians, and researchers, it streamlines the identification of disease-causing variants through sophisticated filtering, Mendelian inheritance analysis, and comprehensive reporting.

## Why Choose VariantCentrifuge?

VariantCentrifuge addresses the critical challenges in clinical variant interpretation by providing:

### 🧬 Advanced Variant Filtering
- **Gene-Centric Analysis:** Focus on specific genes of interest with automatic region extraction
- **Rare Variant Detection:** Identify clinically relevant rare variants using population frequency databases
- **Impact-Based Filtering:** Prioritize high and moderate impact variants using SnpEff annotations
- **Custom Filter Presets:** 20+ pre-configured filters for common clinical scenarios

### 🔬 Clinical Interpretation Features
- **Mendelian Inheritance Analysis:** Automatically detect de novo, recessive, dominant, and compound heterozygous patterns
- **ACMG Classification Support:** Integrate pathogenicity predictions and clinical significance
- **Phenotype Integration:** Link variants to patient phenotypes for enhanced interpretation
- **Custom Variant Scoring:** Apply configurable scoring algorithms for variant prioritization

### 📊 Comprehensive Analysis Tools
- **Gene Burden Testing:** Statistical analysis for case-control studies using Fisher's exact test
- **Compound Heterozygous Detection:** Optimized algorithm (10-50x faster) for identifying compound het variants
- **Multi-Sample Support:** Analyze families, trios, and cohorts with pedigree files
- **Performance Optimization:** Handle large VCF files with chunked processing and parallel execution

### 📈 Professional Reporting
- **Interactive HTML Reports:** Sortable tables, filtering, and search functionality
- **IGV.js Integration:** Visualize variants directly in the browser with genomic context
- **Excel Export:** Generate formatted Excel workbooks for clinical review
- **Detailed Statistics:** Comprehensive variant and gene-level statistics

## Quick Start

```bash
# Basic usage
variantcentrifuge \
  --gene-name BICC1 \
  --vcf-file input.vcf.gz \
  --output-file output.tsv

# With custom filters and Excel output
variantcentrifuge \
  --gene-name BICC1 \
  --vcf-file input.vcf.gz \
  --filters "((dbNSFP_gnomAD_exomes_AC[0] <= 2) & (ANN[ANY].IMPACT has 'HIGH'))" \
  --xlsx
```

## Documentation Contents

```{toctree}
:maxdepth: 2
:caption: User Guide

installation
usage
configuration
guides/index
faq
```

```{toctree}
:maxdepth: 2
:caption: API Reference

api/index
```

```{toctree}
:maxdepth: 1
:caption: Development

development
contributing
changelog
```

## External Dependencies

VariantCentrifuge requires these bioinformatics tools to be installed and available in your PATH:

- **bcftools** - VCF manipulation
- **snpEff** - Functional annotation and BED file generation  
- **SnpSift** - Variant filtering and field extraction
- **bedtools** - BED file operations

## Getting Help

- **Issues:** Report bugs and request features on [GitHub Issues](https://github.com/scholl-lab/variantcentrifuge/issues)
- **Discussions:** Join the conversation by [creating an issue](https://github.com/scholl-lab/variantcentrifuge/issues)
- **Documentation:** Browse the complete documentation on this site

## License

This project is licensed under the [MIT License](https://github.com/scholl-lab/variantcentrifuge/blob/main/LICENSE).

## Indices and tables

* {ref}`genindex`
* {ref}`modindex`
* {ref}`search`