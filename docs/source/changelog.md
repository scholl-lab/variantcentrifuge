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

### Changed
- Documentation migrated from README to structured Sphinx documentation
- Improved configuration documentation with preset examples

### Fixed
- Documentation structure and navigation

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