# Changelog

All notable changes to VariantCentrifuge will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

_No unreleased changes._

## [0.15.0] - 2026-02-22

### Added
- **Modular association testing framework** with pluggable test architecture:
  - Fisher's exact test with carrier/allele modes, odds ratio, and 95% CI
  - Logistic burden regression for binary traits with Firth fallback for separation
  - Linear burden regression for quantitative traits (OLS)
  - SKAT and SKAT-O with pure Python backend (numpy/scipy, thread-safe); R backend deprecated
  - COAST allelic series test (BMV/DMV/PTV categories); pure Python backend default
  - ACAT-O per-gene omnibus via Cauchy combination; ACAT-V per-variant score within SKAT
- **Covariate system**: TSV/CSV covariate files with auto-detected delimiter and categorical encoding
- **PCA integration**: PLINK `.eigenvec`, AKT output, and generic TSV formats; optional AKT subprocess
- **Variant weights**: Beta(MAF), uniform, CADD, REVEL, and combined functional weight schemes
- **Diagnostics**: Genomic inflation factor (lambda_GC), QQ data TSV, optional matplotlib QQ plot
- **JSON config**: `"association"` section in config.json with validation and CLI override precedence
- **Single FDR strategy**: Benjamini-Hochberg correction on ACAT-O p-values only (not per-test)
- **Association testing guide**: Comprehensive user documentation with examples and troubleshooting

## [0.14.0] - 2026-02-18

### Added
- **Modern HTML report** — complete UX overhaul of the individual variant report across 5 phases:
  - **JS Stack Modernization**: DataTables v2, Chart.js (65KB vs Plotly 3.5MB, 98% reduction), Tippy.js tooltips, all assets vendored for offline single-file reports
  - **Summary dashboard**: metric cards (total variants, genes, samples, impact breakdown, top genes), impact distribution and inheritance pattern charts
  - **Semantic color badges**: IMPACT (red/orange/amber/gray), ClinVar (red-to-green severity), inheritance patterns (de novo=red, compound het=purple, AD=blue, AR=green, X-linked=teal)
  - **Table redesign**: sticky GENE column (FixedColumns), dark header, expandable row detail panels, content density toggle (Compact/Regular/Relaxed with localStorage), intelligent column widths, zebra striping
  - **Column-level filtering**: noUiSlider range sliders for numeric columns (POS, gnomAD AF, CADD), categorical dropdowns (IMPACT, ClinVar, Inheritance), text search (GENE), removable filter chips, "Include missing values" toggle, reactive chart updates
  - **Unified toolbar**: all controls in one 28px row — entries/page, filters, missing, search, density, show/hide columns, PDF export
  - **Accessibility (WCAG 2.1 AA)**: skip-link, ARIA roles/labels, keyboard-accessible tooltips, SVG icons with screen-reader text, chart data table fallbacks, 4.5:1 contrast ratios on all badges
  - **Print/PDF support**: @media print stylesheet hiding interactive controls, PDF export via browser print dialog
  - Report metadata footer with filter criteria, VCF source, reference genome, version, and date
  - Expanded summary.json with inheritance distribution, top genes, and sample count
  - Loading skeleton with shimmer animation during DataTable initialization
- 107 HTML report-specific tests (structure, behavior, assets, accessibility, print)
- Unified resource auto-detection across pipeline modes

## [0.13.1] - 2026-02-16

### Fixed
- Resource auto-detection across pipeline modes

## [0.13.0] - 2026-02-16

### Added
- **Performance optimization** across 7 phases:
  - Benchmark framework with timing instrumentation
  - Vectorized genotype replacement (Pandas-native operations)
  - DataFrame optimization with sanitized column names for itertuples compatibility
  - Inheritance analysis optimization with vectorized deduction
  - Output stage optimization
  - Pipeline I/O elimination (reduced intermediate file writes)
  - Parallelization and chunking with memory-aware chunk sizing
- Memory management system with SLURM/PBS/cgroup detection
- Golden file infrastructure for inheritance validation
- Performance benchmark test suite

## [0.12.0] - 2025-12-01

### Added
- **Stage-based pipeline architecture** (`pipeline_core/`) with modular stages, dependency graph, topological sort, and parallel execution via ThreadPoolExecutor/ProcessPoolExecutor
- **Restructured Snakemake workflow** matching lab pipeline conventions (Issue #68):
  - Snakemake 8+ with `min_version("8.0")` and executor plugins
  - Schema-validated config (`config/config_vc.yaml`) and sample sheet (`config/samples.tsv`)
  - Profile layering: `profiles/default/` (resources) + `profiles/{bih,charite,local}/` (executor)
  - Auto-detecting launcher script (`scripts/run_snakemake.sh`) for BIH, Charite, and local
  - Singularity/Apptainer container support via `container:` directive with conda fallback
- **Docker image** on GHCR (`ghcr.io/scholl-lab/variantcentrifuge`) with all bioinformatics tools pre-installed:
  - Multi-stage build using micromamba for minimal image size
  - CI/CD pipeline with automated builds, Trivy security scanning, and cosign image signing
  - `docker-compose.yml` with volume mount patterns for data, snpEff databases, and custom scoring configs
  - Non-root container for security best practices
- **Field profiles** for annotation database version compatibility (`--field-profile`):
  - Built-in `dbnsfp4` and `dbnsfp5` profiles for seamless gnomAD field switching
  - Template syntax `{{fragment:param}}` for parameterized filter presets
  - `--list-field-profiles` to show available profiles
  - Custom profiles configurable in `config.json` without code changes
- **Inheritance analysis** with three-pass pipeline (deduction, compound het, prioritization):
  - Supported patterns: de novo, AD, AR, X-linked (XLR/XLD), compound heterozygous, mitochondrial
  - PED file integration for family-based analysis (`--ped`, `--inheritance-mode`)
  - Vectorized compound het detection (10-50x faster than original)
  - Segregation analysis with Fisher's exact test
  - Pattern prioritization by clinical significance
- **bcftools pre-filtering** (`--bcftools-prefilter`) for early variant filtering during extraction
- **Final filtering** (`--final-filter`) using pandas query syntax on any column including computed scores
- **Sample pseudonymization** for privacy-preserving data sharing (Issue #34):
  - Multiple naming schemas: sequential, categorical, anonymous, and custom patterns
  - Consistent pseudonym mapping across all output formats (TSV, Excel, HTML)
  - PED file pseudonymization support (`--pseudonymize-ped`)
  - Secure mapping table storage in parent directory
- **Checkpoint and resume system** for robust pipeline execution:
  - Automatic pipeline state tracking with `.variantcentrifuge_state.json`
  - Resume capability after interruptions (`--enable-checkpoint` and `--resume`)
  - Optional file checksum validation (`--checkpoint-checksum`)
  - Interactive resume point selection (`--interactive-resume`)
  - Thread-safe state updates for parallel chunk processing
- Unified annotation system supporting BED files, gene lists, and JSON gene data
- JSON gene annotation feature with flexible field mapping (`--annotate-json-genes` and `--json-gene-mapping`)
- Comprehensive Sphinx documentation with modern Furo theme
- GitHub Actions workflow for automated documentation deployment
- Tumor-normal filtering presets (`somatic`, `loh`, `tumor_only`) with configurable sample indices and thresholds
- VCF annotation inspection (`--show-vcf-annotations`) for field discovery
- Genotype filtering (`--genotype-filter`) with per-gene override support
- Transcript-level filtering (`--transcript-list`, `--transcript-file`)
- ClinVar PM5 annotation support (`--clinvar-pm5-lookup`)

### Changed
- Documentation migrated from README to structured Sphinx documentation
- Enhanced filtering with three-stage approach (bcftools pre-filter, SnpSift filter, final filter)

### Fixed
- Numeric type conversion in final filtering to handle mixed data types correctly
- **Gene burden analysis edge cases** (Issue #31): Improved handling of infinite and zero odds ratios

## [0.5.0] - 2024-08-01

### Added
- Interactive HTML report generation with sortable tables
- IGV.js integration for genomic visualization
- Cohort analysis and reporting functionality
- Gene burden analysis with Fisher's exact test
- Phenotype data integration and filtering
- Preset filter system for common analysis workflows
- Excel output format support
- External database links (SpliceAI, Franklin, Varsome, gnomAD, ClinVar)

### Changed
- Improved command-line interface with better argument organization
- Enhanced configuration system with JSON-based presets
- Optimized variant filtering pipeline
- Better error handling and user feedback

### Fixed
- VCF header normalization for indexed fields
- Sample identification in cohort reports
- IGV report generation with proper FASTA handling

## [0.4.0] - 2024-06-01

### Added
- Comprehensive test suite with pytest
- Pre-commit hooks for code quality
- Gene list annotation functionality
- Variant statistics and metadata generation

### Changed
- Modular code architecture with clear separation of concerns
- Improved logging and debugging capabilities
- Enhanced VCF processing pipeline

### Fixed
- Gene BED file generation edge cases
- Genotype replacement functionality

## [0.3.0] - 2024-04-01

### Added
- Gene-centric variant filtering
- SnpSift integration for field extraction
- Basic HTML report generation
- Phenotype integration capabilities

### Changed
- Migrated from Bash/R to Python-based pipeline
- Improved error handling and validation

## [0.2.0] - 2024-02-01

### Added
- Initial Python CLI implementation
- VCF processing with bcftools integration
- Basic filtering capabilities
- Configuration file support

### Changed
- Complete rewrite from shell scripts to Python

## [0.1.0] - 2024-01-01

### Added
- Initial release
- Basic variant filtering functionality
- Shell script-based pipeline
- Simple output generation

---

## Legend

- **Added** for new features
- **Changed** for changes in existing functionality
- **Deprecated** for soon-to-be removed features
- **Removed** for now removed features
- **Fixed** for any bug fixes
- **Security** for vulnerability fixes
