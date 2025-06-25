# Welcome to VariantCentrifuge!

**VariantCentrifuge** is a Python-based command-line tool designed to filter, extract, and refine genetic variant data (VCF files) based on genes of interest, rarity criteria, and impact annotations. Built with modularity and extensibility in mind, VariantCentrifuge replaces the complexity of traditional Bash/R pipelines with a cleaner, maintainable Python codebase.

## Key Features

- **Gene-Centric Filtering:** Extract variants from regions defined by genes of interest
- **Rare Variant Identification:** Apply custom filters to isolate rare and moderate/high-impact variants
- **Flexible Field Extraction:** Easily specify which fields to extract from the VCF
- **Genotype Replacement:** Replace genotype fields with corresponding sample IDs
- **Phenotype Integration:** Integrate phenotype data for enhanced variant analysis
- **Gene List Annotation:** Annotate variant outputs with gene membership information
- **Comprehensive Analysis:** Perform gene burden analyses and variant-level statistics
- **Rich Reporting:** Generate interactive HTML reports with IGV.js integration

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
- **Discussions:** Join the conversation on [GitHub Discussions](https://github.com/scholl-lab/variantcentrifuge/discussions)
- **Documentation:** Browse the complete documentation on this site

## License

This project is licensed under the [MIT License](https://github.com/scholl-lab/variantcentrifuge/blob/main/LICENSE).

## Indices and tables

* {ref}`genindex`
* {ref}`modindex`
* {ref}`search`