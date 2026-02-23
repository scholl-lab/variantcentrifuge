# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Tests
pytest                              # All tests
pytest -m unit                      # Unit tests only
pytest -m integration               # Integration tests only
pytest -m "not slow"                # Skip slow tests
make test-fast                      # Non-slow, non-integration tests

# Code quality
make lint                           # Ruff linter
make format                         # Ruff format + fix
make typecheck                      # Mypy type checker
make ci-check                       # All CI checks locally

# Development
uv pip install -e ".[dev]"          # Install with dev deps
make clean                          # Remove build artifacts
```

## Important Rules

- **NEVER execute `variantcentrifuge` directly** — it processes large genomic datasets requiring specific input files and can run for hours.
- **NEVER commit without explicit user permission** — always ask first.
- Version is managed in `variantcentrifuge/version.py` (`__version__`). All other references import from there.

## Architecture

### Pipeline Architecture
The tool uses a stage-based pipeline architecture invoked through `variantcentrifuge.cli:main`:
- **`pipeline.py`** — Builds and runs the stage-based pipeline using `PipelineRunner`, `PipelineContext`, and registered stages from `pipeline_core/`
- **`pipeline_core/`** — Framework with 40+ modular stages, dependency graph analysis, parallel execution

### Pipeline Core (pipeline_core/)
- **`runner.py`** — PipelineRunner: dependency graph analysis, topological sort, parallel execution via ThreadPoolExecutor/ProcessPoolExecutor, checkpoint integration
- **`context.py`** — PipelineContext: thread-safe dataclass carrying args, config, workspace, data artifacts, completed stages between all stages
- **`stage.py`** — Abstract base Stage class all stages inherit from; declares dependencies, parallel-safety, estimated runtime
- **`workspace.py`** — Manages all intermediate/output file paths and directory structure
- **`error_handling.py`** — Custom exceptions, error propagation

### Stages (stages/)
Four files, each containing multiple Stage subclasses:

| File | Purpose | Examples |
|------|---------|---------|
| `setup_stages.py` | Config/data loading | ConfigurationLoadingStage, PedigreeLoadingStage |
| `processing_stages.py` | VCF processing | BCFToolsPrefilterStage, SnpSiftFilterStage, GenotypeReplacementStage |
| `analysis_stages.py` | Variant analysis | InheritanceAnalysisStage, VariantScoringStage, GeneBurdenAnalysisStage |
| `output_stages.py` | Report generation | HTMLReportStage, ExcelReportStage, IGVReportStage |

`stage_registry.py` handles central stage registration, discovery, and dependency validation.

### Data Flow
1. **Setup** — Load config, create gene BED from gene names, validate VCF, load pedigree/phenotype
2. **Processing** — bcftools prefilter → SnpSift filter → field extraction → genotype replacement (vectorized or sequential)
3. **Analysis** — Inheritance pattern deduction → compound het detection → scoring → statistics → gene burden
4. **Output** — TSV → Excel/HTML/IGV reports, optional pseudonymization and archiving

### Inheritance Analysis (inheritance/)
Three-pass analysis: deduction → compound het → prioritization
- **`analyzer.py`** — Main orchestrator with vectorized and original implementations
- **`deducer.py`** — Pattern deduction (de novo, AD, AR, X-linked, mitochondrial)
- **`comp_het.py`** / **`comp_het_vectorized.py`** — Compound heterozygous detection (original and vectorized)
- **`segregation_checker.py`** — Segregation analysis with Fisher's exact test
- **`prioritizer.py`** — Ranks patterns by clinical significance
- **`parallel_analyzer.py`** — Multi-gene parallel analysis with chunking

### Memory Management (memory/)
`inheritance_memory_manager.py` — Detects system memory (supports SLURM, PBS, cgroups), estimates dataset memory needs, calculates chunk sizes with safety thresholds. Important for HPC environments.

### Checkpoint/Resume System
`checkpoint.py` provides AtomicFileOperation for safe writes, PipelineState tracking, and file validation with checksums. `interactive_resume.py` handles interactive checkpoint selection.

### External Tool Dependencies (must be in PATH)
- `bcftools` — VCF manipulation and prefiltering
- `snpEff` / `SnpSift` — Variant annotation and field extraction
- `bedtools` — BED file operations
- `igv-reports` — IGV.js report generation (optional)

### Configuration
- **`variantcentrifuge/config.json`** — Main config: filter presets (rare, coding, pathogenic, etc.), field extraction definitions, link templates (SpliceAI, Franklin, Varsome, gnomAD, ClinVar), HTML report settings
- **`scoring/*/`** — Scoring models, each with `formula_config.json` and `variable_assignment_config.json` (e.g., `nephro_candidate_score/`, `inheritance_score/`)
- **`stats_configs/*.json`** — Custom statistics configurations

### Snakemake Workflow
`snakemake/variantcentrifuge.smk` — Batch processing of multiple VCFs with HPC cluster support (SLURM, PBS). Configured via `config_vc.yaml`.

## Testing Structure

- `tests/unit/` — Fast unit tests (pipeline core, stages, scoring logic, determinism)
- `tests/integration/` — Full pipeline tests with mocked external tools
- `tests/test_inheritance/` — Inheritance-specific tests with own `conftest.py` (pedigree fixtures, variant DataFrames)
- `tests/mocks/` — Mock implementations for external tools (bcftools, snpEff, etc.)
- `tests/fixtures/` — Test data (gene burden data, GIAB samples)
- `tests/performance/` — Pipeline benchmarks

Pytest markers: `unit`, `integration`, `slow`, `inheritance`, `comp_het`, `segregation`, `deducer`, `prioritizer`, `replacer`, `analyze`, `phenotype`, `helpers`

## Code Style
- Ruff: 100 char line length, Python 3.10+ target (formatting + linting)
- mypy: gradual type checking (ignore_missing_imports)
- Pre-commit hooks enforce ruff plus trailing whitespace, end-of-file fixer, YAML/JSON/TOML checks
