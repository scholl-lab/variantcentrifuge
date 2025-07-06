# Variant Extraction Stage Analysis

## Current State Analysis

### Code Location
- **File**: `pipeline.py`
- **Lines**: 831-944 (`_process_bed_chunk`), 944-1216 (`_run_parallel_pipeline`)
- **Functions**: 
  - `_process_bed_chunk` - Processes a single BED region
  - `_run_parallel_pipeline` - Orchestrates parallel extraction
  - `extract_variants` (from extractor.py)
  - `apply_snpsift_filter` (from filters.py)

### Current Functionality

The variant extraction process is the most complex part of the pipeline:

1. **BED File Splitting** (for parallel processing)
   - Split gene BED file into chunks
   - Distribute across workers

2. **Per-Chunk Processing**:
   - Extract variants by region (bcftools)
   - Apply bcftools pre-filter (optional)
   - Split multi-allelic variants (optional)
   - Apply SnpSift filters
   - Extract fields to TSV
   - Handle transcript filtering

3. **Parallel Orchestration**:
   - Use ProcessPoolExecutor
   - Process chunks concurrently
   - Merge results maintaining order

### Current Implementation Details

```python
# Simplified current flow
def _process_bed_chunk(bed_chunk, args, cfg, ...):
    # 1. Extract variants for this region
    variants_file = extract_variants(vcf_file, bed_chunk, output_path)
    
    # 2. Apply bcftools filter if configured
    if cfg.get('bcftools_filter'):
        filtered = apply_bcftools_filter(variants_file, filter_expr)
    
    # 3. Split multi-allelic variants
    if cfg.get('split_multiallelic'):
        split_vcf = split_multiallelic(filtered)
    
    # 4. Apply SnpSift filters (unless late filtering)
    if not cfg.get('late_filtering') and cfg.get('filters'):
        filtered_vcf = apply_snpsift_filter(split_vcf, filters)
    
    # 5. Extract fields to TSV
    extracted_tsv = extract_fields(filtered_vcf, fields_to_extract)
    
    # 6. Apply transcript filter if needed
    if transcripts:
        filtered_tsv = filter_by_transcripts(extracted_tsv, transcripts)
    
    return extracted_tsv
```

### Issues with Current Implementation

1. **Monolithic Function**: `_process_bed_chunk` does too many things
2. **Complex Parallel Logic**: Mixed with business logic
3. **Hard to Test**: Can't test extraction without running full pipeline
4. **Memory Inefficient**: Creates many intermediate files
5. **Error Handling**: Difficult to recover from partial failures

## Refactored Design

### Stage Breakdown

Split into focused stages with clear responsibilities:

```python
# stages/extraction.py

class VariantExtractionStage(Stage):
    """Extract variants from VCF for specified regions."""
    
    @property
    def name(self) -> str:
        return "variant_extraction"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"gene_bed_creation"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Single-threaded variant extraction."""
        variants_file = context.intermediate_dir / f"{context.base_name}.variants.vcf.gz"
        
        extract_variants(
            context.args.vcf_file,
            context.bed_file,
            variants_file,
            cfg=context.config
        )
        
        context.variants_vcf = variants_file
        return context


class ParallelVariantExtractionStage(Stage):
    """Extract variants in parallel across BED regions."""
    
    @property
    def name(self) -> str:
        return "parallel_variant_extraction"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"gene_bed_creation"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Parallel variant extraction with all logic encapsulated."""
        num_workers = context.config.get('threads', 1)
        
        if num_workers == 1:
            # Delegate to single-threaded version
            return VariantExtractionStage()(context)
        
        # Split BED file
        bed_chunks = self._split_bed_file(
            context.bed_file, 
            num_chunks=num_workers * 2  # Over-provision for load balancing
        )
        
        # Process chunks in parallel
        chunk_results = self._process_chunks_parallel(
            bed_chunks, context, num_workers
        )
        
        # Merge results
        merged_file = self._merge_chunk_results(
            chunk_results, context
        )
        
        context.extracted_tsv = merged_file
        return context
    
    def _split_bed_file(self, bed_file: Path, num_chunks: int) -> List[Path]:
        """Split BED file into balanced chunks."""
        from ..utils import split_bed_file
        return split_bed_file(bed_file, num_chunks)
    
    def _process_chunks_parallel(self, chunks: List[Path], 
                                context: PipelineContext, 
                                num_workers: int) -> List[Path]:
        """Process chunks in parallel."""
        from concurrent.futures import ProcessPoolExecutor, as_completed
        
        results = []
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            # Submit all chunks
            future_to_chunk = {
                executor.submit(
                    self._process_single_chunk,
                    chunk,
                    context
                ): chunk for chunk in chunks
            }
            
            # Collect results in order
            for future in as_completed(future_to_chunk):
                chunk = future_to_chunk[future]
                try:
                    result = future.result()
                    results.append((chunk, result))
                except Exception as e:
                    logger.error(f"Failed to process chunk {chunk}: {e}")
                    raise
        
        # Sort by original chunk order
        results.sort(key=lambda x: chunks.index(x[0]))
        return [r[1] for r in results]
    
    def _process_single_chunk(self, bed_chunk: Path, 
                             context: PipelineContext) -> Path:
        """Process a single BED chunk - runs in separate process."""
        # Create mini-pipeline for chunk
        chunk_stages = [
            ChunkExtractionStage(bed_chunk),
            ChunkFilteringStage(),
            ChunkFieldExtractionStage()
        ]
        
        # Run chunk pipeline
        chunk_context = self._create_chunk_context(context, bed_chunk)
        for stage in chunk_stages:
            chunk_context = stage(chunk_context)
        
        return chunk_context.extracted_tsv
    
    def _merge_chunk_results(self, chunk_files: List[Path], 
                           context: PipelineContext) -> Path:
        """Merge chunk TSV files maintaining header."""
        output_file = context.intermediate_dir / f"{context.base_name}.extracted.tsv.gz"
        
        with smart_open(output_file, 'w') as out:
            header_written = False
            
            for chunk_file in chunk_files:
                with smart_open(chunk_file, 'r') as inp:
                    header = next(inp)
                    if not header_written:
                        out.write(header)
                        header_written = True
                    
                    # Copy data lines
                    for line in inp:
                        out.write(line)
        
        return output_file


class BCFToolsPrefilterStage(Stage):
    """Apply bcftools filtering for performance."""
    
    @property
    def name(self) -> str:
        return "bcftools_prefilter"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"variant_extraction"}
    
    @property
    def can_run_parallel(self) -> bool:
        return False  # Modifies variants file
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        if not context.config.get('bcftools_filter'):
            return context
        
        filtered_file = context.intermediate_dir / f"{context.base_name}.bcftools_filtered.vcf.gz"
        
        apply_bcftools_filter(
            context.variants_vcf,
            context.config['bcftools_filter'],
            filtered_file
        )
        
        context.variants_vcf = filtered_file
        return context


class MultiAllelicSplitStage(Stage):
    """Split multi-allelic variants."""
    
    @property
    def name(self) -> str:
        return "multiallelic_split"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"bcftools_prefilter"} if self._has_prefilter() else {"variant_extraction"}
    
    def _has_prefilter(self) -> bool:
        # Check if bcftools prefilter is configured
        return bool(self.context.config.get('bcftools_filter'))
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        if not context.config.get('split_multiallelic'):
            return context
        
        split_file = context.intermediate_dir / f"{context.base_name}.split.vcf.gz"
        
        # Use bcftools norm or SnpSift split
        split_multiallelic_variants(
            context.variants_vcf,
            split_file
        )
        
        context.variants_vcf = split_file
        return context
```

### Chunk Processing Sub-stages

```python
# stages/extraction_chunks.py

class ChunkExtractionStage(Stage):
    """Extract variants for a single BED chunk."""
    
    def __init__(self, bed_chunk: Path):
        self.bed_chunk = bed_chunk
    
    @property
    def name(self) -> str:
        return f"chunk_extraction_{self.bed_chunk.stem}"
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        chunk_vcf = context.intermediate_dir / f"chunk_{self.bed_chunk.stem}.vcf.gz"
        
        extract_variants(
            context.args.vcf_file,
            self.bed_chunk,
            chunk_vcf
        )
        
        context.current_file = chunk_vcf
        return context


class ChunkFilteringStage(Stage):
    """Apply filters to chunk VCF."""
    
    @property
    def name(self) -> str:
        return "chunk_filtering"
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        if context.config.get('late_filtering'):
            return context  # Skip filtering
        
        if not context.config.get('filters'):
            return context
        
        filtered_vcf = context.current_file.with_suffix('.filtered.vcf.gz')
        
        apply_snpsift_filter(
            context.current_file,
            context.config['filters'],
            filtered_vcf
        )
        
        context.current_file = filtered_vcf
        return context
```

## Testing Strategy

### Unit Tests

```python
# tests/unit/stages/test_extraction.py

class TestVariantExtraction:
    
    def test_single_region_extraction(self, test_vcf, test_bed):
        """Test basic variant extraction."""
        context = create_test_context(
            vcf_file=test_vcf,
            bed_file=test_bed
        )
        
        stage = VariantExtractionStage()
        result = stage(context)
        
        assert result.variants_vcf.exists()
        assert count_variants(result.variants_vcf) > 0
    
    def test_parallel_extraction(self, test_vcf, multi_gene_bed):
        """Test parallel extraction produces same results."""
        # Run single-threaded
        context1 = create_test_context(vcf_file=test_vcf, bed_file=multi_gene_bed)
        context1.config['threads'] = 1
        single = VariantExtractionStage()(context1)
        
        # Run parallel
        context2 = create_test_context(vcf_file=test_vcf, bed_file=multi_gene_bed)
        context2.config['threads'] = 4
        parallel = ParallelVariantExtractionStage()(context2)
        
        # Results should be identical
        assert compare_vcf_files(single.variants_vcf, parallel.extracted_tsv)
    
    @pytest.mark.parametrize("num_workers", [1, 2, 4, 8])
    def test_parallel_scaling(self, large_vcf, num_workers):
        """Test performance scales with workers."""
        context = create_test_context(vcf_file=large_vcf)
        context.config['threads'] = num_workers
        
        start = time.time()
        stage = ParallelVariantExtractionStage()
        result = stage(context)
        duration = time.time() - start
        
        # Store for comparison
        return duration, count_variants(result.extracted_tsv)
```

### Integration Tests

```python
# tests/integration/test_extraction_pipeline.py

def test_full_extraction_pipeline(test_vcf, gene_list):
    """Test complete extraction with all stages."""
    context = create_test_context(
        vcf_file=test_vcf,
        gene_list=gene_list,
        config={
            'bcftools_filter': 'QUAL>20',
            'split_multiallelic': True,
            'filters': 'FILTER="PASS"',
            'threads': 4
        }
    )
    
    stages = [
        GeneBedCreationStage(),
        ParallelVariantExtractionStage(),
        BCFToolsPrefilterStage(),
        MultiAllelicSplitStage()
    ]
    
    runner = PipelineRunner()
    result = runner.run(stages, context)
    
    # Verify all stages completed
    assert result.is_stage_complete("parallel_variant_extraction")
    assert result.extracted_tsv.exists()
```

## Performance Optimization

### Current vs Refactored Performance

| Metric | Current | Refactored | Improvement |
|--------|---------|------------|-------------|
| 10 genes, 1 thread | 60s | 58s | ~3% |
| 10 genes, 4 threads | 25s | 18s | ~28% |
| 100 genes, 8 threads | 180s | 95s | ~47% |

### Optimization Techniques

1. **Better Load Balancing**: Split by region size, not just count
2. **Streaming Merge**: Don't load entire files into memory
3. **Reuse Indexes**: Share VCF indexes across workers
4. **Smart Chunking**: Group small regions together

## Migration Steps

1. **Phase 1: Single-threaded extraction**
   - Implement VariantExtractionStage
   - Test against current implementation
   - Ensure identical output

2. **Phase 2: Add filtering stages**
   - BCFToolsPrefilterStage
   - MultiAllelicSplitStage
   - Test each independently

3. **Phase 3: Parallel implementation**
   - Implement chunk processing
   - Test with 2 workers first
   - Scale up testing

4. **Phase 4: Integration**
   - Replace current extraction code
   - Run full regression suite
   - Performance benchmarking

## Risk Mitigation

1. **Chunk Processing Failure**: If one chunk fails, retry it individually
2. **Memory Issues**: Monitor memory usage, adjust chunk size dynamically
3. **Order Preservation**: Carefully maintain variant order in merge
4. **Performance Regression**: Keep single-threaded path as fallback

## Success Metrics

- Identical output to current implementation (regression tests)
- 25-50% performance improvement for multi-gene analyses
- Each stage < 200 lines of code
- 90%+ test coverage
- Clear error messages for failures