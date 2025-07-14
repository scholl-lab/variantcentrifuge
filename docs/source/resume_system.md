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

- **Stage completion status** (`completed`, `running`, `failed`)
- **Execution times** for performance monitoring
- **Input/output files** for validation
- **Configuration hash** to detect parameter changes
- **Pipeline version** to ensure compatibility

### Stage Lifecycle

Each pipeline stage follows this checkpoint lifecycle:

1. **Start**: Stage begins execution, marked as `running`
2. **Process**: Stage performs its work
3. **Complete**: Stage finishes, marked as `completed` with output files recorded
4. **Skip**: On resume, completed stages are skipped and their state is restored

## Configuration

### Enabling Checkpoints

```bash
# Enable checkpoint system
variantcentrifuge --enable-checkpoint [other options]

# Enable with custom output directory
variantcentrifuge --enable-checkpoint --output-dir /path/to/output [other options]
```

### Resume Options

| Option | Description | Example |
|--------|-------------|---------|
| `--resume` | Resume from last completed stage | `--resume` |
| `--resume-from STAGE` | Resume from specific stage | `--resume-from dataframe_loading` |
| `--enable-checkpoint` | Enable checkpoint system | `--enable-checkpoint` |

## Stage Types and Resume Behavior

### File-Based Stages

Stages that produce output files (e.g., `parallel_complete_processing`, `genotype_replacement`):

```python
# Example: GenotypeReplacementStage
def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
    """Restore file paths when stage is skipped."""
    output_tsv = context.workspace.get_intermediate_path(
        f"{context.workspace.base_name}.genotype_replaced.tsv.gz"
    )
    if output_tsv.exists():
        context.genotype_replaced_tsv = output_tsv
        context.data = output_tsv
    return context
```

### Memory-Based Stages

Stages that work with DataFrames (e.g., `dataframe_loading`, `custom_annotation`):

```python
# Example: DataFrameLoadingStage
def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
    """Restore DataFrame from TSV file when stage is skipped."""
    input_file = self._find_input_file(context)
    df = pd.read_csv(input_file, sep="\t", dtype=str, compression="gzip")
    context.current_dataframe = df
    return context
```

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

The system handles interrupted pipeline runs by detecting and cleaning up stale states:

### Stale Running Stages

```python
# Cleanup logic for interrupted stages
if step_info.status == "running":
    # Check if output files exist
    if step_info.output_files:
        all_outputs_exist = all(
            os.path.exists(file_info.path) 
            for file_info in step_info.output_files
        )
        if all_outputs_exist:
            # Stage actually completed
            step_info.status = "completed"
            return True
    
    # Determine if stage was interrupted
    if this_was_last_stage_started:
        step_info.status = "failed"
        step_info.error = "Pipeline was interrupted"
```

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
  ✓ custom_annotation (0.3s)
  ✗ inheritance_analysis (19.7s)
    Error: Pipeline was interrupted
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

### 1. Use Checkpoints for Long-Running Analyses

```bash
# Always enable checkpoints for large datasets
variantcentrifuge --vcf-file large_cohort.vcf.gz --enable-checkpoint [options]
```

### 2. Consistent Output Directories

```bash
# Use same output directory for resume
variantcentrifuge --output-dir /consistent/path --enable-checkpoint [options]
```

### 3. Parameter Validation

```bash
# Verify parameters before resuming
variantcentrifuge --resume --dry-run [options]
```

### 4. Development Workflow

```bash
# Enable checkpoints during development
variantcentrifuge --enable-checkpoint --resume-from analysis_stage [options]
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
        # Check if stage should be skipped
        if context.checkpoint_state.should_skip_step(self.name):
            context.mark_complete(self.name)
            
            # Restore state if needed
            if hasattr(self, '_handle_checkpoint_skip'):
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