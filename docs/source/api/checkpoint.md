# checkpoint

```{eval-rst}
.. automodule:: variantcentrifuge.checkpoint
   :members:
   :undoc-members:
   :show-inheritance:
```

## Overview

The checkpoint module provides pipeline state tracking and resume functionality for VariantCentrifuge. It enables robust handling of long-running analyses by saving progress at key pipeline steps and allowing resumption after interruptions.

## Key Classes

### PipelineState

The main class for managing pipeline checkpoints. It tracks:
- Pipeline steps and their completion status
- Input/output files for each step
- Configuration hash to ensure consistency
- Pipeline version for compatibility checks

### CheckpointContext

A context manager that automatically tracks step execution:
- Records start and end times
- Captures input and output files
- Handles errors gracefully
- Saves state on successful completion

### FileInfo

Tracks file metadata for validation:
- File path, size, and modification time
- Optional checksum calculation
- Validation to ensure files haven't changed

### StepInfo

Stores information about each pipeline step:
- Step name and status (pending, running, completed, failed)
- Execution times
- Input/output files
- Step-specific parameters
- Error information if failed

## Usage Example

```python
from variantcentrifuge.checkpoint import PipelineState, CheckpointContext

# Initialize pipeline state
pipeline_state = PipelineState(output_dir, enable_checksum=False)
pipeline_state.initialize(config, pipeline_version)

# Check if we can resume from a previous run
if pipeline_state.load() and pipeline_state.can_resume(config, pipeline_version):
    logger.info("Resuming from checkpoint")

# Use checkpoint context for a pipeline step
if not pipeline_state.should_skip_step("variant_extraction"):
    with CheckpointContext(pipeline_state, "variant_extraction") as ctx:
        # Add input files
        ctx.add_input_file(input_vcf)
        
        # Perform the step
        output_file = extract_variants(input_vcf, output_path)
        
        # Add output files
        ctx.add_output_file(output_file)
        
        # State is automatically saved on successful completion
```

## Parallel Processing Support

The checkpoint system is designed to work with parallel processing:
- Thread-safe state updates
- Tracks individual chunks in parallel runs
- Proper ordering of parallel step completion
- Handles worker failures gracefully

## State File Format

The checkpoint state is stored as JSON in `.variantcentrifuge_state.json`:

```json
{
  "version": "1.0",
  "pipeline_version": "0.5.0",
  "start_time": 1234567890.0,
  "configuration_hash": "abc123...",
  "steps": {
    "gene_bed_creation": {
      "name": "gene_bed_creation",
      "status": "completed",
      "start_time": 1234567890.0,
      "end_time": 1234567891.0,
      "input_files": [],
      "output_files": [
        {
          "path": "output/genes.bed",
          "size": 12345,
          "mtime": 1234567891.0,
          "checksum": null
        }
      ],
      "parameters": {},
      "error": null
    }
  }
}
```