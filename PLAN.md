# Plan: Fix Parallel Stage Resume Issue

## Problem Analysis

The issue is that parallel extraction stages (like `ParallelVariantExtractionStage` and `ParallelCompleteProcessingStage`) are completely repeating all work during resume, even when individual chunks or substeps have already been completed. This happens because:

1. **Missing Checkpoint Skip Logic**: `ParallelVariantExtractionStage` lacks `_handle_checkpoint_skip` method
2. **No Granular Checkpoint Tracking**: Individual chunks within parallel stages aren't tracked in the checkpoint system
3. **Chunk-Level Work Repetition**: The parallel processing code blindly re-executes all chunks without checking if intermediate files exist

## Root Cause

**Current Logic Flow:**
1. Pipeline resumes and marks `parallel_complete_processing` as complete if checkpoint exists
2. But `ParallelVariantExtractionStage` doesn't have checkpoint skip logic
3. Individual chunks aren't checkpointed, so all bcftools/snpsift work repeats
4. Even if chunk outputs exist, the parallel processing code recreates them

## Solution Strategy

### Main Step Resume Control
- User-facing resume functionality will control main pipeline steps only
- Users can resume from major stages like `variant_extraction`, `analysis`, `report_generation`
- Main step resume provides clear, predictable control points

### Automatic Substep Resume
- General resume system will automatically detect and reuse existing substep files
- Chunk-level files (bcftools outputs, filtered VCFs, extracted TSVs) will be validated and reused
- No user intervention needed for substep-level resumption

## Implementation Plan

### 1. **Add Checkpoint Skip Logic to ParallelVariantExtractionStage**
- Implement `_handle_checkpoint_skip()` method
- Restore expected output files and context state
- Mark constituent stages as complete

### 2. **Implement Automatic Substep File Detection**
- Create chunk-level file existence checking before processing
- Skip individual chunks that have already been processed
- Only process missing/incomplete chunks
- Validate existing files (size, timestamp, content integrity)

### 3. **Enhanced Parallel Processing Logic**
- Modify `_process_chunks_parallel` to check for existing chunk outputs
- Skip chunks where final TSV files already exist
- Implement robust chunk validation
- Handle mixed completion scenarios (some chunks done, others missing)

### 4. **Improve ParallelCompleteProcessingStage Resume**
- Add chunk-level output validation in `_handle_checkpoint_skip`
- Better detection of partial completion states
- Smarter reconstruction of context state from existing files

### 5. **Main Step Resume Controls**
- Expose only major pipeline stages for user-controlled resume
- Provide clear stage names and descriptions
- Maintain simple, predictable resume points

## Architecture

```
Main Step Resume (User Control)
├── gene_bed_creation
├── variant_extraction  ← User resumes here
│   ├── chunk_0.vcf.gz     ← Auto-detected/reused
│   ├── chunk_1.vcf.gz     ← Auto-detected/reused
│   └── chunk_N.vcf.gz     ← Auto-detected/reused
├── filtering_and_extraction
│   ├── chunk_0.tsv.gz     ← Auto-detected/reused
│   ├── chunk_1.tsv.gz     ← Auto-detected/reused
│   └── merged.tsv.gz      ← Auto-detected/reused
├── analysis              ← User resumes here
└── report_generation     ← User resumes here
```

## Expected Benefits

- **Intuitive User Experience**: Simple main-step resume controls
- **Massive Performance Improvement**: Automatic reuse of completed substeps
- **Robust Resumption**: Handle partial parallel stage completion gracefully
- **Better Resource Usage**: Avoid redundant bcftools/snpsift operations
- **No User Complexity**: Automatic substep detection requires no user knowledge

## Files to Modify

1. `variantcentrifuge/stages/processing_stages.py` - Add checkpoint skip logic and chunk validation
2. `variantcentrifuge/checkpoint.py` - Enhanced substep file detection (if needed)
3. `variantcentrifuge/pipeline_core/runner.py` - Improved parallel stage handling
4. `variantcentrifuge/cli.py` - Main step resume controls
5. `tests/unit/stages/test_parallel_complete_processing.py` - Add resumption tests

## Implementation Priority

1. **High Priority**: Automatic substep file detection and reuse
2. **Medium Priority**: Enhanced parallel stage checkpoint skip logic
3. **Low Priority**: Main step resume UI improvements

This approach provides the best of both worlds: simple user controls for main steps while automatically optimizing substep resumption behind the scenes.