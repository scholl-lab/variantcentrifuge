# VariantCentrifuge

[![CI](https://img.shields.io/github/actions/workflow/status/scholl-lab/variantcentrifuge/test.yml?branch=main&label=CI)](https://github.com/scholl-lab/variantcentrifuge/actions/workflows/test.yml)
[![Docker](https://img.shields.io/github/actions/workflow/status/scholl-lab/variantcentrifuge/docker.yml?branch=main&label=Docker)](https://github.com/scholl-lab/variantcentrifuge/actions/workflows/docker.yml)
[![Docs](https://img.shields.io/github/actions/workflow/status/scholl-lab/variantcentrifuge/docs.yml?branch=main&label=Docs)](https://github.com/scholl-lab/variantcentrifuge/actions/workflows/docs.yml)
[![Release](https://img.shields.io/github/v/release/scholl-lab/variantcentrifuge)](https://github.com/scholl-lab/variantcentrifuge/releases)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)](https://python.org)
[![License](https://img.shields.io/github/license/scholl-lab/variantcentrifuge)](LICENSE)

A command-line tool for filtering, extracting, and prioritizing genetic variants from VCF files.
VariantCentrifuge combines gene-centric region extraction, multi-tier filtering (bcftools, SnpSift, pandas), inheritance analysis, and configurable scoring into a single reproducible pipeline.

## Features

- Gene-centric variant extraction using gene names or BED regions
- Three-tier filtering: bcftools prefilter, SnpSift expressions, pandas final filter
- Inheritance pattern analysis (de novo, AD, AR, X-linked, compound het)
- Configurable variant scoring models
- Gene burden analysis with Fisher's exact test
- **Modular association testing**: SKAT-O, COAST allelic series, logistic/linear burden regression with covariates and PCA, ACAT-O omnibus combination — see the [Association Testing Guide](docs/source/guides/association_testing.md)
- Interactive HTML reports with column-level filtering, semantic color badges, summary dashboard, accessible design (WCAG 2.1 AA), and PDF export
- ClinVar, gnomAD, and SpliceAI annotation links
- Cohort aggregation across multiple samples
- Field profiles for switching annotation database versions (e.g., dbNSFP v4/v5)
- Docker image with all bioinformatics dependencies included
- Stage-based pipeline architecture with parallel execution

## Installation

**Docker** (recommended -- all tools included):

```bash
docker pull ghcr.io/scholl-lab/variantcentrifuge:latest
```

The production image contains an immutable runtime environment. `pip` is used
only while building the image and is removed before distribution; package
installation belongs in a development environment.

**pip:**

```bash
pip install variantcentrifuge
```

**From source:**

```bash
git clone https://github.com/scholl-lab/variantcentrifuge.git
cd variantcentrifuge && pip install .
```

External tools (bcftools, snpEff, SnpSift, bedtools) must be in PATH when not using Docker.
Install via conda: `mamba create -y -n vc bcftools snpsift snpeff bedtools`

## Quick Start

```bash
# Filter rare coding variants in a gene list
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --preset rare,coding \
  --html-report \
  --xlsx

# Score and filter with a custom model
variantcentrifuge \
  --gene-file genes.txt \
  --vcf-file input.vcf.gz \
  --preset rare,coding \
  --scoring-config-path scoring/nephro_variant_score \
  --final-filter 'score > 0.8 and IMPACT == "HIGH"' \
  --output-file results.tsv
```

## Snakemake Workflow

A Snakemake 8+ workflow for batch-processing multiple VCFs on HPC clusters (SLURM, PBS) is included under `workflow/`, with cluster profiles in `profiles/` and sample configuration in `config/`.
See `scripts/run_snakemake.sh` for the auto-detecting launcher.

## Documentation

Full documentation: [scholl-lab.github.io/variantcentrifuge](https://scholl-lab.github.io/variantcentrifuge/)

- [Installation Guide](https://scholl-lab.github.io/variantcentrifuge/installation.html)
- [Usage Guide](https://scholl-lab.github.io/variantcentrifuge/usage.html)
- [Configuration](https://scholl-lab.github.io/variantcentrifuge/configuration.html)
- [API Reference](https://scholl-lab.github.io/variantcentrifuge/api/)

## Contributing

Contributions are welcome. Please see the [Contributing Guide](https://scholl-lab.github.io/variantcentrifuge/contributing.html) for details.

- [Bug Reports and Feature Requests](https://github.com/scholl-lab/variantcentrifuge/issues)
- [Discussions](https://github.com/scholl-lab/variantcentrifuge/discussions)

## How to cite

If you use VariantCentrifuge in your research, cite the software and record the
exact release or commit used in your analysis. [CITATION.cff](CITATION.cff)
provides machine-readable citation metadata for GitHub's **Cite this repository**
menu. Until an archive DOI is available, cite the corresponding GitHub release
or commit URL. For example, for the latest published release:

> Popp, B. (2026). VariantCentrifuge (version 0.17.12) [Computer software].
> https://github.com/scholl-lab/variantcentrifuge/releases/tag/v0.17.12

The current development version is **0.17.13**, as recorded in
`variantcentrifuge/version.py`, `CITATION.cff`, and [.zenodo.json](.zenodo.json).
It has no release date yet; the latest published release is **v0.17.12**
(2026-07-19). For a development checkout, include its full commit hash and URL
instead of citing it as a published release.

**Zenodo archive pending:** metadata is prepared, but no DOI has been minted.
See the [release and archival guide](docs/source/citation.md) for the maintainer
steps to archive a release and add the resulting DOI badge.

## License

This project is licensed under the [MIT License](LICENSE).
