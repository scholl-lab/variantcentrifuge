# VariantCentrifuge

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://scholl-lab.github.io/variantcentrifuge/)
[![GitHub Issues](https://img.shields.io/github/issues/scholl-lab/variantcentrifuge)](https://github.com/scholl-lab/variantcentrifuge/issues)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Python Version](https://img.shields.io/badge/python-3.7+-blue.svg)](https://python.org)

**VariantCentrifuge** is a Python-based command-line tool designed to filter, extract, and refine genetic variant data (VCF files) based on genes of interest, rarity criteria, and impact annotations. Built with modularity and extensibility in mind, VariantCentrifuge replaces the complexity of traditional Bash/R pipelines with a cleaner, maintainable Python codebase.

## ✨ Key Features

- **🎯 Gene-Centric Filtering:** Extract variants from regions defined by genes of interest
- **🔍 Rare Variant Identification:** Apply custom filters to isolate rare and moderate/high-impact variants
- **⚙️ Flexible Configuration:** JSON-based configuration with reusable filter presets
- **📊 Interactive Reports:** Generate HTML reports with sortable tables and IGV.js integration
- **🧬 Gene Burden Analysis:** Perform statistical analysis with Fisher's exact test
- **🔗 Clinical Integration:** ClinVar, gnomAD, and external database annotations
- **👥 Cohort Analysis:** Aggregate results from multiple samples with interactive visualizations
- **📦 Results Archiving:** Automatically create compressed archives of complete analysis results
- **💾 Space Optimization:** Gzip compression for intermediate files to reduce disk usage

## 🚀 Quick Start

### Installation

```bash
# Install from PyPI (recommended)
pip install variantcentrifuge

# Or install from source
git clone https://github.com/scholl-lab/variantcentrifuge.git
cd variantcentrifuge
pip install .
```

### Basic Usage

```bash
# Analyze variants in a single gene
variantcentrifuge \\\n  --gene-name BRCA1 \\\n  --vcf-file input.vcf.gz \\\n  --output-file brca1_variants.tsv

# Use predefined filters for rare, coding variants
variantcentrifuge \\\n  --gene-file cancer_genes.txt \\\n  --vcf-file input.vcf.gz \\\n  --preset rare,coding \\\n  --html-report \\\n  --xlsx

# Analysis with space optimization and result archiving
variantcentrifuge \\\n  --gene-file cancer_genes.txt \\\n  --vcf-file input.vcf.gz \\\n  --preset rare,coding \\\n  --gzip-intermediates \\\n  --archive-results \\\n  --html-report
```

## 📋 Prerequisites

**External Tools** (must be in PATH):
- `bcftools` - VCF manipulation
- `snpEff` - Functional annotation
- `SnpSift` - Variant filtering and field extraction
- `bedtools` - BED file operations

**Install via conda:**
```bash
mamba create -y -n variantcentrifuge bcftools snpsift snpeff bedtools
mamba activate variantcentrifuge
```

## 📖 Documentation

**📚 [Complete Documentation](https://scholl-lab.github.io/variantcentrifuge/)**

### Quick Links

- **[Installation Guide](https://scholl-lab.github.io/variantcentrifuge/installation.html)** - Detailed setup instructions
- **[Usage Guide](https://scholl-lab.github.io/variantcentrifuge/usage.html)** - Command-line options and examples
- **[Configuration](https://scholl-lab.github.io/variantcentrifuge/configuration.html)** - Filter presets and customization
- **[API Reference](https://scholl-lab.github.io/variantcentrifuge/api/)** - Developer documentation

### Practical Guides

- **[Annotation Strategies](https://scholl-lab.github.io/variantcentrifuge/guides/annotation_strategies.html)** - VCF annotation best practices
- **[Cohort Analysis](https://scholl-lab.github.io/variantcentrifuge/guides/cohort_analysis.html)** - Multi-sample analysis workflows
- **[Rare Disease Workflow](https://scholl-lab.github.io/variantcentrifuge/guides/rare_disease_workflow.html)** - Clinical variant analysis
- **[Cancer Analysis](https://scholl-lab.github.io/variantcentrifuge/guides/cancer_analysis.html)** - Somatic variant workflows

## 🎯 Use Cases

### Rare Disease Analysis
```bash
variantcentrifuge \\\n  --gene-file disease_genes.txt \\\n  --vcf-file patient.vcf.gz \\\n  --preset rare_pathogenic,high_confidence \\\n  --phenotype-file patient_data.tsv \\\n  --html-report \\\n  --output-file rare_disease_analysis.tsv
```

### Cancer Genomics
```bash
variantcentrifuge \\\n  --gene-file oncogenes_tsg.txt \\\n  --vcf-file tumor_normal.vcf.gz \\\n  --preset mutect2_TvsN,coding \\\n  --igv \\\n  --bam-mapping-file bam_files.tsv \\\n  --html-report
```

### Population Genetics
```bash
variantcentrifuge \\\n  --gene-file population_genes.txt \\\n  --vcf-file cohort.vcf.gz \\\n  --preset 5percent,coding \\\n  --perform-gene-burden \\\n  --html-report
```

## 🏗️ Example Configuration

```json
{
  "reference": "GRCh38.99",
  "filters": "",
  "fields_to_extract": "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P gnomAD_exomes_AF ClinVar_CLNSIG GEN[*].GT",
  "presets": {
    "rare": "(((gnomAD_exomes_AF < 0.0001) | (na gnomAD_exomes_AF)) & ((gnomAD_genomes_AF < 0.0001) | (na gnomAD_genomes_AF)))",
    "coding": "((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE'))",
    "pathogenic": "((ClinVar_CLNSIG =~ '[Pp]athogenic') & !(ClinVar_CLNSIG =~ '[Cc]onflicting'))"
  }
}
```

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](https://scholl-lab.github.io/variantcentrifuge/contributing.html) for details.

- **🐛 Bug Reports:** [GitHub Issues](https://github.com/scholl-lab/variantcentrifuge/issues)
- **💡 Feature Requests:** [GitHub Issues](https://github.com/scholl-lab/variantcentrifuge/issues)
- **💬 Discussions:** [GitHub Discussions](https://github.com/scholl-lab/variantcentrifuge/discussions)

## 📄 License

This project is licensed under the [MIT License](LICENSE).

## 🙏 Acknowledgments

- Built upon the rich ecosystem of bioinformatics tools (snpEff, SnpSift, bcftools, bedtools)
- Special thanks to contributors and the open-source community
- Inspired by prior Bash/R pipelines for variant filtering

---

**📖 [View Full Documentation](https://scholl-lab.github.io/variantcentrifuge/) | 🚀 [Get Started](https://scholl-lab.github.io/variantcentrifuge/installation.html) | 💬 [Join Discussion](https://github.com/scholl-lab/variantcentrifuge/discussions)**