# Resume System

VariantCentrifuge includes a robust checkpoint and resume system that allows you to recover from interrupted pipeline runs, significantly reducing computational waste and improving the user experience for long-running analyses.

## Overview

The resume system automatically saves pipeline state to disk as stages complete, allowing you to:

- **Resume interrupted pipelines** without losing progress
- **Skip completed stages** and continue from where you left off
- **Recover from system failures** or accidental interruptions
- **Optimize development workflows** by resuming from specific stages

## Quick Start

### Basic Resume

```bash
# Run initial pipeline
variantcentrifuge --gene-name PKD1 --vcf-file input.vcf.gz --enable-checkpoint

# If interrupted, resume with the same command + --resume
variantcentrifuge --gene-name PKD1 --vcf-file input.vcf.gz --enable-checkpoint --resume
```

### Resume from Specific Stage

```bash
# Resume from a specific stage (advanced usage)
variantcentrifuge --gene-name PKD1 --vcf-file input.vcf.gz --enable-checkpoint --resume-from dataframe_loading
```

## How It Works

### Checkpoint State File

The resume system creates a `.variantcentrifuge_state.json` file in your output directory that tracks:

- **Stage completion status** (`pending`, `running`, `finalizing`, `completed`, `failed`)
- **Execution times** for performance monitoring
- **Input/output files** with optional checksums for validation
- **Configuration hash** to detect parameter changes
- **Pipeline version** to ensure compatibility

### Stage Lifecycle

Each pipeline stage follows this checkpoint lifecycle:

1. **Start**: Stage begins execution, marked as `running`
2. **Process**: Stage performs its main work
3. **Finalize** (optional): Stage enters `finalizing` state while moving temp files to final locations
4. **Complete**: Stage finishes, marked as `completed` with output files recorded
5. **Resume action**: On resume, completed stages are restored, skipped, or
   recomputed according to their explicit checkpoint resume policy

### Atomic File Operations

For reliability, stages use atomic file operations to prevent partial file corruption:

```python
from variantcentrifuge.checkpoint import AtomicFileOperation

# Safe file creation
with AtomicFileOperation('/path/to/final/output.tsv') as temp_path:
    # Write to temporary file
    with open(temp_path, 'w') as f:
        f.write("data")
    # File is atomically moved to final location on success
```

## Configuration

### Enabling Checkpoints

```bash
# Enable checkpoint system
variantcentrifuge --enable-checkpoint [other options]

# Enable with checksum validation (RECOMMENDED FOR PRODUCTION)
variantcentrifuge --enable-checkpoint --checkpoint-checksum [other options]

# Enable with custom output directory
variantcentrifuge --enable-checkpoint --output-dir /path/to/output [other options]
```

### Checksum Validation (Production Recommendation)

For maximum reliability, especially in production environments, enable checksum validation:

```bash
# Production-grade checkpoint validation
variantcentrifuge --enable-checkpoint --checkpoint-checksum [other options]
```

**Benefits of checksum validation:**
- **File integrity verification**: Detects corrupted or truncated files
- **Reliable recovery**: Only recovers stages with verified complete outputs
- **Production safety**: Prevents silent failures from partial files

**Without checksum validation:**
- Faster checkpoint operations (size/time-based validation only)
- Interrupted stages are always re-executed for safety
- Suitable for development or when performance is critical

### Resume Options

| Option | Description | Example |
|--------|-------------|---------|
| `--resume` | Resume from last completed stage | `--resume` |
| `--resume-from STAGE` | Restart from specific stage (re-execute stage and all subsequent stages) | `--resume-from dataframe_loading` |
| `--enable-checkpoint` | Enable checkpoint system | `--enable-checkpoint` |
| `--checkpoint-checksum` | Enable checksum validation (recommended for production) | `--checkpoint-checksum` |

## Stage Types and Resume Behavior

Completed checkpoint entries are not all equivalent after a process restart. Some
stages only produce durable files, while others mutate in-memory DataFrames that no
longer exist in the resumed process. VariantCentrifuge uses an explicit resume policy
for each stage:

- `restore`: the stage may be skipped only after `_handle_checkpoint_skip()` rebuilds
  the downstream context state from durable artifacts.
- `skip`: the stage has no downstream runtime state to restore, or is an intentional
  no-op in the current configuration.
- `recompute`: the completed checkpoint entry is not trusted after process restart; the
  stage and downstream stages are re-executed.

The conservative default is `recompute`.

### Durable File Stages

Stages that produce durable output files can be restored only when those files exist and
validate. Missing expected artifacts abort resume instead of logging a warning and
continuing.

Example: `parallel_complete_processing` restores the merged extracted TSV and marks the
constituent extraction/filtering stages complete only after the merged TSV validates:

```python
# Example: ParallelCompleteProcessingStage
def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
    if not merged_tsv.exists() or not self._validate_chunk_tsv(merged_tsv):
        raise RuntimeError("Cannot restore parallel_complete_processing from checkpoint")

    context.extracted_tsv = merged_tsv
    context.data = merged_tsv
    context.config["parallel_vcf_processing_complete"] = True
    context.mark_complete("variant_extraction")
    context.mark_complete("field_extraction")
    context.mark_complete("snpsift_filtering")
    return context
```

`GenotypeReplacementStage` is a Phase-11 no-op. Its resume skip is therefore a clean
no-op and does not require a `genotype_replaced.tsv` file.

### Restorable Memory Stages

Some memory-based stages can rebuild their runtime state from durable TSV artifacts.
`dataframe_loading` reloads the current DataFrame when full loading is appropriate, or
defers to `chunked_analysis` for large files. `chunked_analysis` restores
`context.current_dataframe` from `chunked_analysis_results.tsv(.gz)`, including
header-only empty results.

```python
# Example: DataFrameLoadingStage
def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
    input_file = self._find_input_file(context)
    if not input_file or not Path(input_file).exists():
        raise RuntimeError("Cannot restore dataframe_loading from checkpoint")

    if self._should_use_chunks(context, input_file) and not context.config.get(
        "dataframe_chunked_analysis_complete", False
    ):
        context.config["use_chunked_processing"] = True
        return context

    df, rename_map = load_optimized_dataframe(str(input_file), sep="\t")
    context.current_dataframe = df
    context.column_rename_map = rename_map
    return context
```

The old `chunked_processing_complete` meaning has been split:

- `parallel_vcf_processing_complete`: VCF chunk extraction/filter/merge is complete.
- `dataframe_chunked_analysis_complete`: DataFrame chunked analysis has produced a
  restorable DataFrame TSV.

The parallel VCF flag does not suppress DataFrame chunking.

### Volatile Memory Stages

Stages such as custom annotation, inheritance analysis, variant analysis, scoring,
final filtering, pseudonymization, and ClinVar PM5 mutate DataFrames in memory without
durable post-stage DataFrame artifacts. After process restart, these stages are
recomputed even if the checkpoint says they were completed.

This recompute behavior must clear the persisted checkpoint completion for the first
volatile stage and all downstream stages, because both `PipelineRunner` and
`Stage.__call__()` can skip completed checkpoint stages. Clearing only in-memory
`PipelineContext.completed_stages` is not sufficient.

### Composite Stages

Stages that represent multiple operations (e.g., `parallel_complete_processing`):

```python
# Example: ParallelCompleteProcessingStage
def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
    """Mark constituent stages as complete when composite stage is skipped."""
    # Restore output file
    context.extracted_tsv = merged_tsv

    # Mark constituent stages as complete
    context.mark_complete("variant_extraction")
    context.mark_complete("field_extraction")
    context.mark_complete("snpsift_filtering")

    return context
```

### Required Outputs and Empty Results

Requested analysis and output stages fail when required runtime state is missing. For
example, requested association, gene burden, statistics, TSV output, and Excel output
do not mark their checkpoint stage complete when `context.current_dataframe` or the
required TSV artifact is absent.

An empty DataFrame is different from a missing DataFrame. Valid zero-row results produce
header-only requested outputs, including final TSV, association sidecar, and gene burden
sidecar files.

## Resume Validation

The resume system performs several validation checks:

### Configuration Consistency

```python
# Configuration hash validation
current_hash = self._hash_configuration(configuration)
stored_hash = self.state.get("configuration_hash")
if current_hash != stored_hash:
    logger.warning("Configuration has changed, cannot resume")
    return False
```

### Pipeline Version Compatibility

```python
# Version compatibility check
if self.state.get("pipeline_version") != pipeline_version:
    logger.warning("Pipeline version mismatch, cannot resume")
    return False
```

### File Existence Validation

```python
# Output file validation
for step_name, step_info in self.state["steps"].items():
    if step_info.status == "completed":
        for file_info in step_info.output_files:
            if not file_info.validate():
                logger.warning(f"Output file validation failed: {file_info.path}")
                return False
```

## Stale State Recovery

The system handles interrupted pipeline runs with robust recovery logic that prioritizes data integrity:

### Improved Recovery Strategy

**Safety-First Approach**: For maximum reliability, interrupted stages (`running` or `finalizing`) are always re-executed unless checksum validation confirms complete files.

```python
# Enhanced recovery logic
if step_info.status in ("running", "finalizing"):
    logger.warning(f"Found stale {step_info.status} stage, marking for re-execution")

    # Only attempt recovery with checksum validation enabled
    if self.enable_checksum and step_info.output_files:
        logger.info("Attempting checksum-based recovery")

        # Validate all output files with checksums
        all_valid = True
        for file_info in step_info.output_files:
            if not file_info.validate(calculate_checksum=True):
                all_valid = False
                break

        if all_valid:
            # Files verified - mark as completed
            step_info.status = "completed"
            return True

    # Mark for re-execution (safe default)
    step_info.status = "failed"
    step_info.error = "Pipeline was interrupted - stage will be re-executed for safety"
```

### Recovery Modes

**1. Checksum Validation Mode** (`--checkpoint-checksum`):
- Validates file integrity using SHA256 checksums
- Recovers stages only when files are verified as complete
- **Recommended for production pipelines**

**2. Conservative Mode** (default):
- Always re-executes interrupted stages for safety
- Faster checkpoint operations (no checksum calculation)
- Prevents potential issues from partial/corrupt files

### Finalizing State

The new `finalizing` state indicates stages that are moving temporary files to final locations:

```python
# Stage progression with finalizing state
context.checkpoint_state.start_step("example_stage")          # running
context.checkpoint_state.finalize_step("example_stage")       # finalizing
context.checkpoint_state.complete_step("example_stage")       # completed
```

This intermediate state ensures that:
- Partial files are never considered complete
- Atomic file operations are properly tracked
- Recovery logic can distinguish between main processing and file finalization

## Advanced Usage

### Resume from Specific Stage

```bash
# Resume from dataframe_loading stage
variantcentrifuge --gene-name PKD1 --vcf-file input.vcf.gz \
    --enable-checkpoint --resume-from dataframe_loading

# Resume from variant_analysis stage
variantcentrifuge --gene-name PKD1 --vcf-file input.vcf.gz \
    --enable-checkpoint --resume-from variant_analysis
```

`--resume-from STAGE` uses restart semantics: the target stage and all subsequent stages
are re-executed. The target can be a previously failed or incomplete stage. Before the
target, the runner restores only safe durable/restorable prerequisites. If the requested
target depends on volatile in-memory state that cannot be restored, resume fails with an
earlier safe stage to use instead, such as `--resume-from dataframe_loading`.

For issue #101-style runs, fixing resume state restoration can correctly expose the
original `association_analysis` failure again. That association categorical genotype
error is separate from checkpoint resume correctness.

### Development Workflow

```bash
# Initial run with checkpoint
variantcentrifuge --gene-name PKD1 --vcf-file input.vcf.gz --enable-checkpoint

# Modify scoring configuration and resume from scoring stage
variantcentrifuge --gene-name PKD1 --vcf-file input.vcf.gz --enable-checkpoint \
    --scoring-config-path new_scoring --resume-from variant_scoring

# Resume from report generation
variantcentrifuge --gene-name PKD1 --vcf-file input.vcf.gz --enable-checkpoint \
    --resume-from tsv_output
```

## Pipeline State Summary

The resume system provides detailed state information:

```
Pipeline State Summary:
  State file: /path/to/output/.variantcentrifuge_state.json
  Pipeline version: 0.5.0
  Started: 2025-07-14 12:59:31

Steps:
  ✓ parallel_complete_processing (22.0s)
  ✓ genotype_replacement (2.9s)
  ✓ phenotype_integration (0.1s)
  ✓ dataframe_loading (0.1s)
  ◐ custom_annotation (finalizing)
  → inheritance_analysis (running)
  ✗ variant_scoring
    Error: Pipeline was interrupted - stage will be re-executed for safety
```

## Performance Benefits

### Time Savings

- **Skip completed stages**: Avoid re-running expensive operations
- **Parallel stage optimization**: Resume from parallel processing points
- **Memory efficiency**: Restore only necessary data structures

### Development Efficiency

- **Iterative development**: Test changes without full pipeline runs
- **Parameter tuning**: Resume from analysis stages with new parameters
- **Report generation**: Re-generate reports without reprocessing data

## Troubleshooting

### Common Issues

#### Configuration Mismatch

```
Cannot resume - configuration or version mismatch
```

**Solution**: Ensure all parameters match the original run, or start fresh without `--resume`.

#### Missing Output Files

```
Output file validation failed for step 'stage_name': /path/to/file.tsv
```

**Solution**: Check that intermediate files weren't manually deleted. Remove checkpoint file to start fresh.

#### Missing Required Checkpoint Artifact

```
Cannot restore chunked_analysis from checkpoint: chunked output missing
```

**Solution**: Resume from an earlier durable stage that can regenerate the missing
artifact, or start fresh. The pipeline intentionally aborts rather than reporting
success with missing requested outputs.

#### Corrupted Checkpoint

```
Checkpoint appears corrupted (all stages were stale)
```

**Solution**: Delete the `.variantcentrifuge_state.json` file and start fresh.

### Manual Checkpoint Management

```bash
# Remove checkpoint file to start fresh
rm /path/to/output/.variantcentrifuge_state.json

# Check checkpoint status
variantcentrifuge --enable-checkpoint --resume --dry-run
```

## Best Practices

### 1. Production Environments

```bash
# Production-grade configuration with checksum validation
variantcentrifuge --vcf-file large_cohort.vcf.gz \
    --enable-checkpoint --checkpoint-checksum \
    --output-dir /reliable/path [options]
```

**Production recommendations:**
- Always use `--checkpoint-checksum` for file integrity verification
- Use dedicated output directories with sufficient disk space
- Monitor checkpoint file sizes (larger with checksums enabled)
- Consider backup strategies for checkpoint files in critical workflows

### 2. Development Environments

```bash
# Development configuration optimized for speed
variantcentrifuge --vcf-file test_data.vcf.gz \
    --enable-checkpoint \
    --resume-from analysis_stage [options]
```

**Development recommendations:**
- Checksum validation optional (faster iteration)
- Use `--resume-from` for rapid testing of specific stages
- Clean checkpoint files when changing major parameters

### 3. Consistent Output Directories

```bash
# Use same output directory for resume
variantcentrifuge --output-dir /consistent/path --enable-checkpoint [options]
```

### 4. Parameter Validation

```bash
# Verify parameters before resuming
variantcentrifuge --resume --dry-run [options]
```

### 5. File Safety Best Practices

```python
# Use atomic file operations in custom stages
from variantcentrifuge.checkpoint import AtomicFileOperation

def write_output_safely(data, output_path):
    with AtomicFileOperation(output_path) as temp_path:
        # Write to temp file
        with open(temp_path, 'w') as f:
            f.write(data)
        # Automatically moved to final location on success
```

## Technical Implementation

### State Management

The `PipelineState` class manages checkpoint state:

```python
class PipelineState:
    def __init__(self, output_dir: str, enable_checksum: bool = False):
        self.state_file = os.path.join(output_dir, ".variantcentrifuge_state.json")
        self.state = {
            "version": "1.0",
            "pipeline_version": None,
            "start_time": None,
            "steps": {},
            "configuration_hash": None
        }
```

### Stage Integration

Stages integrate with the checkpoint system through:

```python
class Stage(ABC):
    def __call__(self, context: PipelineContext) -> PipelineContext:
        # Stage.__call__ has its own checkpoint skip gate. The runner clears
        # persisted completion first for recompute stages so this gate cannot
        # skip volatile in-memory work after process restart.
        if context.checkpoint_state.should_skip_step(self.name):
            context.mark_complete(self.name)
            if hasattr(self, "_handle_checkpoint_skip"):
                context = self._handle_checkpoint_skip(context)

            return context

        # Execute stage normally
        result = self._process(context)

        # Mark complete with output files
        context.checkpoint_state.complete_step(
            self.name,
            output_files=self.get_output_files(context)
        )

        return result
```

## Future Enhancements

### Planned Features

- **Incremental updates**: Resume with modified input files
- **Parallel resume**: Resume multiple branches simultaneously
- **Cloud storage**: Store checkpoints in cloud storage
- **Checkpoint compression**: Reduce checkpoint file size
- **Visual resume interface**: GUI for selecting resume points

### API Extensions

```python
# Future API for programmatic resume control
from variantcentrifuge.checkpoint import PipelineCheckpoint

checkpoint = PipelineCheckpoint("/path/to/output")
if checkpoint.can_resume():
    pipeline.resume_from(checkpoint.get_last_stage())
```

## Related Documentation

- [Configuration Guide](configuration.md) - Pipeline configuration options
- [Performance Tips](guides/performance_tips.md) - Optimization strategies
- [Development Guide](development.md) - Development workflow
- [API Reference](api/checkpoint.md) - Checkpoint API documentation
