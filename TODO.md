# TODO: Parallel Stage Resume Implementation

## High Priority Tasks

### 1. Fix ParallelVariantExtractionStage Resume ✅ COMPLETED
- [x] Add `_handle_checkpoint_skip()` method to `ParallelVariantExtractionStage`
- [x] Implement output file detection and context restoration
- [x] Test that stage properly skips when already completed

### 2. Implement Automatic Substep File Detection ✅ COMPLETED  
- [x] Add chunk file existence checking in `_process_chunks_parallel` methods
- [x] Create `_validate_existing_chunk()` method for file integrity checks
- [x] Skip processing for valid existing chunk files
- [x] Handle mixed completion scenarios (some chunks done, others missing)

### 3. Enhance ParallelCompleteProcessingStage ✅ COMPLETED
- [x] Improve `_handle_checkpoint_skip()` to validate chunk-level outputs
- [x] Add logic to detect partial completion states
- [x] Reconstruct context state from existing intermediate files
- [x] Test with various partial completion scenarios

## Medium Priority Tasks

### 4. Improve Chunk Validation Logic ✅ COMPLETED
- [x] Add file size validation for chunk outputs
- [x] Implement timestamp-based staleness detection
- [x] Add basic content integrity checks (header validation for TSV files)
- [x] Handle corrupted/partial chunk files gracefully

### 5. Enhanced Error Handling ✅ COMPLETED
- [x] Add robust error handling for chunk validation failures
- [x] Implement fallback to full processing when validation fails
- [x] Add detailed logging for chunk skip/reuse decisions
- [x] Handle edge cases (empty files, permission issues)

### 6. Main Step Resume Controls (CLI) ✅ COMPLETED
- [x] Define clear main step resume points in CLI
- [x] Add `--resume-from` option with main stage names
- [x] Provide clear stage descriptions and estimated progress
- [x] Add validation for resume point names

## Low Priority Tasks

### 7. Performance Optimizations
- [ ] Implement parallel chunk validation
- [ ] Add caching for file metadata checks
- [ ] Optimize file I/O for large numbers of chunks
- [ ] Add progress reporting for chunk validation

### 8. Enhanced Testing ✅ COMPLETED
- [x] Add unit tests for chunk validation logic
- [x] Test partial parallel stage completion scenarios
- [x] Add integration tests for mixed resume states
- [x] Test performance with large numbers of chunks
- [x] Add regression tests for resume functionality

### 9. Documentation Updates
- [ ] Update resume system documentation with substep behavior
- [ ] Document main step resume controls
- [ ] Add troubleshooting guide for resume issues
- [ ] Update CLI help with resume examples

## Implementation Details

### Chunk Validation Strategy
```python
def _validate_existing_chunk(self, chunk_path: Path) -> bool:
    """Validate that an existing chunk file is complete and valid."""
    if not chunk_path.exists():
        return False
    
    # Size check (must be > 0)
    if chunk_path.stat().st_size == 0:
        return False
    
    # Basic content validation for TSV files
    if chunk_path.suffix == '.tsv' or chunk_path.name.endswith('.tsv.gz'):
        return self._validate_tsv_header(chunk_path)
    
    return True
```

### Chunk Skip Logic
```python
def _process_chunks_parallel(self, context, bed_chunks):
    """Process chunks, skipping existing valid outputs."""
    chunks_to_process = []
    existing_outputs = []
    
    for i, chunk_bed in enumerate(bed_chunks):
        expected_output = self._get_expected_chunk_output(context, i)
        
        if self._validate_existing_chunk(expected_output):
            logger.info(f"Reusing existing chunk {i}: {expected_output}")
            existing_outputs.append(expected_output)
        else:
            logger.info(f"Processing missing/invalid chunk {i}")
            chunks_to_process.append((i, chunk_bed))
    
    # Process only missing chunks
    new_outputs = self._process_missing_chunks(context, chunks_to_process)
    
    # Combine existing and new outputs
    all_outputs = existing_outputs + new_outputs
    return sorted(all_outputs)
```

## Testing Strategy

### Unit Tests
- Test chunk validation with various file states
- Test mixed completion scenarios
- Test error handling for corrupted files

### Integration Tests
- Test full pipeline resume with partial chunk completion
- Test resume after various interruption points
- Test performance with many chunks

### Manual Testing
- Test with real VCF files and interruptions
- Verify chunk reuse saves time
- Test main step resume controls

## Success Criteria

1. **Automatic Substep Reuse**: Parallel stages automatically detect and reuse valid chunk files
2. **Performance Improvement**: Resume time reduced by 50-90% when chunks exist
3. **Robust Validation**: Handle corrupted/partial files gracefully
4. **Simple User Interface**: Main step resume controls are intuitive
5. **Backward Compatibility**: Existing resume functionality continues to work
6. **Comprehensive Testing**: All edge cases covered by tests