# Changelog

All notable changes to VariantCentrifuge will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive Sphinx documentation with modern Furo theme
- GitHub Actions workflow for automated documentation deployment
- API reference documentation with autodoc
- User guides for installation, usage, and configuration
- Development and contributing guides
- Annotation strategy guides for VCF preprocessing
- Unified annotation system supporting BED files, gene lists, and JSON gene data
- JSON gene annotation feature with flexible field mapping (`--annotate-json-genes` and `--json-gene-mapping`)
- Custom annotation integration in the pipeline workflow
- **bcftools pre-filtering** (`--bcftools-prefilter`) for early variant filtering during extraction, significantly improving performance on large VCFs
- **Final filtering** (`--final-filter`) using pandas query syntax, allowing filtering on any column including computed scores and inheritance patterns
- Comprehensive test suite for new filtering features
- **Sample pseudonymization** for privacy-preserving data sharing (Issue #34):
  - Multiple naming schemas: sequential, categorical, anonymous, and custom patterns
  - Consistent pseudonym mapping across all output formats (TSV, Excel, HTML)
  - Automatic handling of genotype and inheritance columns
  - PED file pseudonymization support (`--pseudonymize-ped`)
  - Secure mapping table storage in parent directory
  - Comprehensive test suite and documentation

### Changed
- Documentation migrated from README to structured Sphinx documentation
- Improved configuration documentation with preset examples
- Enhanced filtering capabilities with three-stage filtering approach (bcftools pre-filter, SnpSift filter, final filter)

### Fixed
- Documentation structure and navigation
- Numeric type conversion in final filtering to handle mixed data types correctly
- **Gene burden analysis edge cases** (Issue #31): Improved handling of infinite and zero odds ratios by:
  - Detecting structural zeros (e.g., entire row/column zero) and returning NaN appropriately
  - Applying continuity correction (default 0.5) to zero cells to calculate meaningful confidence intervals
  - Using score method as primary CI calculation (more robust for sparse data)
  - Removing arbitrary fallback bounds in favor of statistically sound methods

## [0.5.0] - 2024-XX-XX

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

## [0.4.0] - 2024-XX-XX

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

## [0.3.0] - 2024-XX-XX

### Added
- Gene-centric variant filtering
- SnpSift integration for field extraction
- Basic HTML report generation
- Phenotype integration capabilities

### Changed
- Migrated from Bash/R to Python-based pipeline
- Improved error handling and validation

## [0.2.0] - 2024-XX-XX

### Added
- Initial Python CLI implementation
- VCF processing with bcftools integration
- Basic filtering capabilities
- Configuration file support

### Changed
- Complete rewrite from shell scripts to Python

## [0.1.0] - 2024-XX-XX

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