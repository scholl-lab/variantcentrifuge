# API Reference

This section provides auto-generated documentation from the VariantCentrifuge source code.

## Core Modules

```{toctree}
:maxdepth: 1

cli
pipeline
config
```

## Analysis Components

```{toctree}
:maxdepth: 1

analyze_variants
gene_burden
stats
```

## Data Processing

```{toctree}
:maxdepth: 1

filters
extractor
replacer
converter
```

## Utilities and Helpers

```{toctree}
:maxdepth: 1

utils
helpers
validators
```

## Report Generation

```{toctree}
:maxdepth: 1

generate_html_report
generate_igv_report
links
```

## Data Integration

```{toctree}
:maxdepth: 1

phenotype
phenotype_filter
gene_bed
```

## Module Overview

### Core Pipeline

- **{doc}`cli`** - Command-line interface and argument parsing
- **{doc}`pipeline`** - Main workflow orchestration and coordination
- **{doc}`config`** - Configuration file loading and validation

### Analysis & Statistics

- **{doc}`analyze_variants`** - Variant-level analysis and gene burden testing
- **{doc}`gene_burden`** - Statistical methods for gene burden analysis
- **{doc}`stats`** - Summary statistics and data aggregation

### Data Processing

- **{doc}`filters`** - SnpSift-based variant filtering operations
- **{doc}`extractor`** - VCF field extraction utilities
- **{doc}`replacer`** - Genotype replacement with sample IDs
- **{doc}`converter`** - TSV to Excel conversion and formatting

### Utilities

- **{doc}`utils`** - Common utility functions and external tool integration
- **{doc}`helpers`** - Data manipulation and processing helpers
- **{doc}`validators`** - Input validation and sanity checking

### Reporting

- **{doc}`generate_html_report`** - Interactive HTML report generation
- **{doc}`generate_igv_report`** - IGV.js integration for genomic visualization
- **{doc}`links`** - External database link generation

### Data Integration

- **{doc}`phenotype`** - Phenotype data loading and integration
- **{doc}`phenotype_filter`** - Phenotype-based filtering operations
- **{doc}`gene_bed`** - Gene coordinate processing and BED file generation