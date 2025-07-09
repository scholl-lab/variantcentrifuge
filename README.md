# VariantCentrifuge

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://scholl-lab.github.io/variantcentrifuge/)
[![GitHub Issues](https://img.shields.io/github/issues/scholl-lab/variantcentrifuge)](https://github.com/scholl-lab/variantcentrifuge/issues)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Python Version](https://img.shields.io/badge/python-3.7+-blue.svg)](https://python.org)

**VariantCentrifuge** is a production-ready Python-based command-line tool designed to filter, extract, and refine genetic variant data (VCF files) based on genes of interest, rarity criteria, and impact annotations. Built with modularity and extensibility in mind, VariantCentrifuge replaces the complexity of traditional Bash/R pipelines with a cleaner, maintainable Python codebase.

> **🏆 NEW: Stage-Based Pipeline Architecture** - Now featuring a completely refactored modular pipeline with 36 specialized stages for enhanced performance, maintainability, and parallel execution capabilities. Enable with `--use-new-pipeline` flag.

## ✨ Key Features

- **🎯 Gene-Centric Filtering:** Extract variants from regions defined by genes of interest
- **🔍 Advanced Filtering Pipeline:**
  - **bcftools pre-filtering** for performance optimization on large VCFs
  - **SnpSift filtering** for complex variant selection
  - **Final filtering** with pandas query syntax on computed columns (scores, inheritance patterns)
- **🏗️ Modern Architecture:** Choose between classic monolithic pipeline or new stage-based architecture
  - **Stage-Based Pipeline** (`--use-new-pipeline`): 36 modular stages with parallel execution
  - **Classic Pipeline** (default): Battle-tested monolithic implementation
- **⚙️ Flexible Configuration:** JSON-based configuration with reusable filter presets
- **📊 Interactive Reports:** Generate HTML reports with sortable tables and IGV.js integration
- **🧬 Gene Burden Analysis:** Perform statistical analysis with Fisher's exact test
- **🔗 Clinical Integration:** ClinVar, gnomAD, and external database annotations
- **👥 Cohort Analysis:** Aggregate results from multiple samples with interactive visualizations
- **🎨 Custom Scoring:** Apply configurable variant scoring models without code changes
- **⚡ Performance Optimized:** Parallel processing, memory-efficient streaming, and chunked analysis

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

# Apply custom variant scoring
variantcentrifuge \\\n  --gene-name GENE \\\n  --vcf-file input.vcf.gz \\\n  --scoring-config-path scoring/nephro_variant_score \\\n  --output-file scored_variants.tsv

# Performance-optimized filtering with final filter on scores
variantcentrifuge \\\n  --gene-file genes.txt \\\n  --vcf-file large_cohort.vcf.gz \\\n  --bcftools-prefilter 'FILTER="PASS" && INFO/AC<10' \\\n  --preset rare,coding \\\n  --scoring-config-path scoring/my_model \\\n  --final-filter 'score > 0.8 and IMPACT == "HIGH"' \\\n  --output-file high_priority_variants.tsv
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

### Variant Scoring
```bash
variantcentrifuge \\\n  --gene-file kidney_genes.txt \\\n  --vcf-file patient.vcf.gz \\\n  --preset rare,coding \\\n  --scoring-config-path scoring/nephro_variant_score \\\n  --html-report \\\n  --output-file scored_kidney_variants.tsv
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
- Inspired by prior Bash/R pipelines for variant filtering

---

**📖 [View Full Documentation](https://scholl-lab.github.io/variantcentrifuge/) | 🚀 [Get Started](https://scholl-lab.github.io/variantcentrifuge/installation.html) | 💬 [Join Discussion](https://github.com/scholl-lab/variantcentrifuge/discussions)**