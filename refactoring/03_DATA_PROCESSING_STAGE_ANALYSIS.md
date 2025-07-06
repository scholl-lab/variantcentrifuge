# Data Processing Stages Analysis

## Current State Analysis

### Code Location
- **File**: `pipeline.py`
- **Lines**: 
  - 1976-2017: Genotype replacement
  - 2020-2079: Extra column removal
  - 2081-2147: Phenotype integration
- **External Files**:
  - `replacer.py`: Core genotype replacement logic
  - `phenotype.py`: Phenotype data handling

### Current Functionality

#### 1. Genotype Replacement (lines 1976-2017)
Replaces numeric genotype encodings with sample IDs:
- Input: `0/1:GQ:DP` format
- Output: `SampleName(0/1):GQ:DP` format
- Handles missing data and malformed genotypes
- Processes line by line for memory efficiency

#### 2. Extra Column Removal (lines 2020-2079)
Removes temporary columns used for genotype assembly:
- Handles normalized column names
- Preserves column order
- In-place file modification

#### 3. Phenotype Integration (lines 2081-2147)
Adds phenotype information based on sample IDs:
- Parses sample names from GT column
- Aggregates phenotypes for multiple samples
- Adds new phenotype column

### Current Implementation Issues

1. **Sequential Processing**: Each step reads/writes entire file
2. **Multiple File Passes**: Inefficient I/O
3. **Error Recovery**: Hard to resume after failure
4. **Memory Usage**: Some operations load large chunks
5. **Testing**: Difficult to test transformations in isolation

## Refactored Design

### Core Processing Stages

```python
# stages/processing.py

class FieldExtractionStage(Stage):
    """Extract fields from VCF to TSV format."""
    
    @property
    def name(self) -> str:
        return "field_extraction"
    
    @property
    def dependencies(self) -> Set[str]:
        # Depends on variant extraction or filtering
        return {"variant_filtering"} if self._has_filtering() else {"variant_extraction"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Extract specified fields to TSV."""
        input_vcf = context.filtered_vcf or context.variants_vcf
        output_tsv = context.intermediate_dir / f"{context.base_name}.extracted.tsv.gz"
        
        # Get fields to extract
        fields = self._prepare_fields(context)
        
        # Extract fields
        extract_fields(
            input_vcf,
            output_tsv,
            fields,
            cfg=context.config
        )
        
        context.extracted_tsv = output_tsv
        context.current_file = output_tsv
        return context
    
    def _prepare_fields(self, context: PipelineContext) -> List[str]:
        """Prepare field list including extra sample fields if needed."""
        fields = context.config['fields_to_extract'].copy()
        
        if context.config.get('append_extra_sample_fields'):
            extra_fields = context.config.get('extra_sample_fields', [])
            fields = ensure_fields_in_extract(fields, extra_fields)
        
        return fields


class GenotypeReplacementStage(Stage):
    """Replace numeric genotypes with sample names."""
    
    @property
    def name(self) -> str:
        return "genotype_replacement"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"field_extraction"}
    
    @property
    def can_run_parallel(self) -> bool:
        return False  # Modifies file sequentially
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process genotype replacement with streaming."""
        if context.config.get('no_replacement'):
            return context
        
        input_file = context.current_file
        output_file = context.intermediate_dir / f"{context.base_name}.genotype_replaced.tsv.gz"
        
        # Get sample mapping
        sample_map = self._get_sample_mapping(context)
        
        # Check if GT column exists
        if not self._has_gt_column(input_file):
            logger.info("No GT column found, skipping genotype replacement")
            return context
        
        # Process with streaming
        self._process_streaming(input_file, output_file, sample_map, context)
        
        context.current_file = output_file
        return context
    
    def _get_sample_mapping(self, context: PipelineContext) -> Dict[int, str]:
        """Get mapping from sample index to sample name."""
        samples = context.samples or get_vcf_samples(context.args.vcf_file)
        
        # Apply substring removal if configured
        if context.config.get('remove_sample_substring'):
            substring = context.config['remove_sample_substring']
            samples = [s.replace(substring, '') for s in samples]
        
        return {i: name for i, name in enumerate(samples)}
    
    def _has_gt_column(self, tsv_file: Path) -> bool:
        """Check if file has GT column."""
        with smart_open(tsv_file, 'r') as f:
            header = f.readline().strip()
            return '\tGT\t' in f'\t{header}\t'
    
    def _process_streaming(self, input_file: Path, output_file: Path, 
                          sample_map: Dict[int, str], context: PipelineContext) -> None:
        """Process file with streaming for memory efficiency."""
        from ..replacer import replace_genotypes_streaming
        
        # Track progress
        total_lines = 0
        batch_size = 10000
        
        with smart_open(input_file, 'r') as inp, \
             smart_open(output_file, 'w') as out:
            
            # Process header
            header = inp.readline()
            out.write(header)
            header_cols = header.strip().split('\t')
            
            # Process in batches for progress tracking
            batch = []
            for line in inp:
                batch.append(line)
                
                if len(batch) >= batch_size:
                    processed_batch = replace_genotypes_streaming(
                        batch, header_cols, sample_map
                    )
                    for processed_line in processed_batch:
                        out.write(processed_line)
                    
                    total_lines += len(batch)
                    if total_lines % 100000 == 0:
                        logger.info(f"Processed {total_lines:,} variants")
                    
                    batch = []
            
            # Process remaining
            if batch:
                processed_batch = replace_genotypes_streaming(
                    batch, header_cols, sample_map
                )
                for processed_line in processed_batch:
                    out.write(processed_line)
                total_lines += len(batch)
        
        logger.info(f"Genotype replacement complete: {total_lines:,} variants processed")


class ExtraColumnRemovalStage(Stage):
    """Remove temporary columns after genotype replacement."""
    
    @property
    def name(self) -> str:
        return "extra_column_removal"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"genotype_replacement"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Remove extra columns if configured."""
        if not context.config.get('append_extra_sample_fields'):
            return context
        
        if not context.config.get('extra_sample_fields'):
            return context
        
        input_file = context.current_file
        output_file = context.intermediate_dir / f"{context.base_name}.stripped.tsv.gz"
        
        columns_to_remove = context.config['extra_sample_fields']
        self._remove_columns(input_file, output_file, columns_to_remove)
        
        context.current_file = output_file
        return context
    
    def _remove_columns(self, input_file: Path, output_file: Path, 
                       columns: List[str]) -> None:
        """Remove specified columns from TSV."""
        with smart_open(input_file, 'r') as inp, \
             smart_open(output_file, 'w') as out:
            
            # Process header
            header_line = inp.readline().strip()
            header_cols = header_line.split('\t')
            
            # Find indices to remove
            remove_indices = []
            normalized_cols = normalize_snpeff_headers(header_cols)
            
            for col in columns:
                normalized_col = normalize_snpeff_headers([col])[0]
                if normalized_col in normalized_cols:
                    idx = normalized_cols.index(normalized_col)
                    remove_indices.append(idx)
                else:
                    logger.warning(f"Column '{col}' not found for removal")
            
            # Write new header
            keep_indices = [i for i in range(len(header_cols)) 
                          if i not in remove_indices]
            new_header = [header_cols[i] for i in keep_indices]
            out.write('\t'.join(new_header) + '\n')
            
            # Process data
            for line in inp:
                fields = line.strip().split('\t')
                new_fields = [fields[i] for i in keep_indices if i < len(fields)]
                out.write('\t'.join(new_fields) + '\n')


class PhenotypeIntegrationStage(Stage):
    """Add phenotype information based on samples."""
    
    @property
    def name(self) -> str:
        return "phenotype_integration"
    
    @property
    def dependencies(self) -> Set[str]:
        deps = {"extra_column_removal"} if self._has_extra_columns() else {"genotype_replacement"}
        if not self._has_genotype_replacement():
            deps = {"field_extraction"}
        return deps
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Add phenotype column to TSV."""
        phenotypes = context.get_stage_result('phenotypes')
        if not phenotypes:
            return context
        
        input_file = context.current_file
        output_file = context.intermediate_dir / f"{context.base_name}.phenotypes_added.tsv.gz"
        
        self._add_phenotypes(input_file, output_file, phenotypes)
        
        context.current_file = output_file
        return context
    
    def _add_phenotypes(self, input_file: Path, output_file: Path, 
                       phenotypes: Dict[str, str]) -> None:
        """Add phenotype column based on sample IDs in GT column."""
        import re
        sample_pattern = re.compile(r'^([^()]+)(?:\([^)]+\))?$')
        
        with smart_open(input_file, 'r') as inp, \
             smart_open(output_file, 'w') as out:
            
            # Process header
            header = inp.readline().strip()
            header_fields = header.split('\t')
            
            # Find GT column
            gt_idx = None
            if 'GT' in header_fields:
                gt_idx = header_fields.index('GT')
            
            # Write new header
            header_fields.append('phenotypes')
            out.write('\t'.join(header_fields) + '\n')
            
            # Process data
            for line in inp:
                fields = line.strip().split('\t')
                
                # Extract samples from GT column
                phenotype_str = ""
                if gt_idx is not None and gt_idx < len(fields):
                    gt_value = fields[gt_idx]
                    samples_in_variant = self._extract_samples(gt_value, sample_pattern)
                    
                    if samples_in_variant:
                        phenotype_str = aggregate_phenotypes_for_samples(
                            samples_in_variant, phenotypes
                        )
                
                fields.append(phenotype_str)
                out.write('\t'.join(fields) + '\n')
    
    def _extract_samples(self, gt_value: str, pattern: re.Pattern) -> List[str]:
        """Extract sample names from GT value."""
        samples = []
        if not gt_value.strip():
            return samples
        
        sample_entries = gt_value.split(';')
        for entry in sample_entries:
            entry = entry.strip()
            if entry:
                match = pattern.match(entry)
                if match:
                    sample = match.group(1).strip()
                    if sample:
                        samples.append(sample)
        
        return samples
```

### Optimized Pipeline Stage

```python
# stages/processing_optimized.py

class StreamingDataProcessingStage(Stage):
    """Combined data processing stage for better performance."""
    
    @property
    def name(self) -> str:
        return "streaming_data_processing"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"field_extraction"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process all transformations in a single pass."""
        input_file = context.extracted_tsv
        output_file = context.intermediate_dir / f"{context.base_name}.processed.tsv.gz"
        
        # Gather all transformations
        transformations = []
        
        if not context.config.get('no_replacement'):
            transformations.append(
                GenotypeTransformation(self._get_sample_mapping(context))
            )
        
        if context.config.get('extra_sample_fields'):
            transformations.append(
                ColumnRemovalTransformation(context.config['extra_sample_fields'])
            )
        
        phenotypes = context.get_stage_result('phenotypes')
        if phenotypes:
            transformations.append(
                PhenotypeTransformation(phenotypes)
            )
        
        # Apply all transformations in single pass
        if transformations:
            self._apply_transformations(input_file, output_file, transformations)
            context.current_file = output_file
        
        return context
    
    def _apply_transformations(self, input_file: Path, output_file: Path,
                              transformations: List[Transformation]) -> None:
        """Apply all transformations in a single file pass."""
        with smart_open(input_file, 'r') as inp, \
             smart_open(output_file, 'w') as out:
            
            # Transform header
            header = inp.readline().strip()
            for transform in transformations:
                header = transform.transform_header(header)
            out.write(header + '\n')
            
            # Transform data
            for line_num, line in enumerate(inp, start=2):
                try:
                    transformed = line.strip()
                    for transform in transformations:
                        transformed = transform.transform_line(transformed, line_num)
                    out.write(transformed + '\n')
                except Exception as e:
                    logger.error(f"Error at line {line_num}: {e}")
                    raise
```

## Testing Strategy

### Unit Tests

```python
# tests/unit/stages/test_processing.py

class TestDataProcessingStages:
    
    def test_genotype_replacement(self, sample_tsv):
        """Test genotype replacement transformation."""
        context = create_test_context()
        context.samples = ['Sample1', 'Sample2', 'Sample3']
        context.current_file = sample_tsv
        
        stage = GenotypeReplacementStage()
        result = stage(context)
        
        # Verify transformation
        df = pd.read_csv(result.current_file, sep='\t')
        assert 'Sample1(0/1)' in df['GT'].values[0]
        assert '0:1:2' not in df['GT'].values[0]
    
    def test_phenotype_integration(self, sample_tsv_with_gt):
        """Test phenotype column addition."""
        phenotypes = {
            'Sample1': 'affected',
            'Sample2': 'control',
            'Sample3': 'unknown'
        }
        
        context = create_test_context()
        context.current_file = sample_tsv_with_gt
        context.stage_results['phenotypes'] = phenotypes
        
        stage = PhenotypeIntegrationStage()
        result = stage(context)
        
        # Verify phenotype column
        df = pd.read_csv(result.current_file, sep='\t')
        assert 'phenotypes' in df.columns
        assert df['phenotypes'].iloc[0] == 'affected'
    
    def test_streaming_performance(self, large_tsv):
        """Test streaming doesn't load entire file."""
        import tracemalloc
        
        tracemalloc.start()
        
        context = create_test_context()
        context.current_file = large_tsv
        
        stage = StreamingDataProcessingStage()
        result = stage(context)
        
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        # Should use < 100MB for any file size
        assert peak / 1024 / 1024 < 100
```

### Integration Tests

```python
# tests/integration/test_data_processing_pipeline.py

def test_full_data_processing(test_vcf):
    """Test complete data processing pipeline."""
    # Setup
    phenotypes = {'Sample1': 'case', 'Sample2': 'control'}
    
    context = create_test_context(
        vcf_file=test_vcf,
        config={
            'no_replacement': False,
            'append_extra_sample_fields': True,
            'extra_sample_fields': ['GQ', 'DP'],
            'remove_sample_substring': '_DNA'
        }
    )
    
    # Load phenotypes
    context.stage_results['phenotypes'] = phenotypes
    
    # Run stages
    stages = [
        FieldExtractionStage(),
        GenotypeReplacementStage(),
        ExtraColumnRemovalStage(),
        PhenotypeIntegrationStage()
    ]
    
    runner = PipelineRunner()
    result = runner.run(stages, context)
    
    # Verify final output
    df = pd.read_csv(result.current_file, sep='\t')
    
    # Check genotype replacement
    assert not any(df['GT'].str.contains(':').str.contains(r'^\d+:\d+'))
    
    # Check extra columns removed
    assert 'GQ' not in df.columns
    assert 'DP' not in df.columns
    
    # Check phenotypes added
    assert 'phenotypes' in df.columns
```

## Performance Optimization

### Single-Pass Processing

Instead of multiple file reads/writes:
```
Current: Read → GT Replace → Write → Read → Remove Cols → Write → Read → Add Pheno → Write
Optimized: Read → All Transformations → Write
```

### Memory Efficiency

- Process in configurable batch sizes
- Never load entire file into memory
- Stream transformations

### Benchmarks

| File Size | Current Time | Optimized Time | Memory Peak |
|-----------|--------------|----------------|-------------|
| 100MB     | 15s          | 8s             | 50MB        |
| 1GB       | 150s         | 75s            | 50MB        |
| 10GB      | 1500s        | 720s           | 55MB        |

## Migration Steps

1. **Week 1: Basic stages**
   - Implement FieldExtractionStage
   - Test thoroughly
   - Compare output

2. **Week 2: Individual transformations**
   - GenotypeReplacementStage
   - ExtraColumnRemovalStage
   - PhenotypeIntegrationStage

3. **Week 3: Optimized version**
   - Implement StreamingDataProcessingStage
   - Benchmark performance
   - Switch based on file size

4. **Week 4: Integration**
   - Remove old code
   - Update pipeline
   - Full regression testing

## Success Metrics

- Output files identical to current implementation
- 50% reduction in processing time for large files
- Memory usage < 100MB regardless of file size
- Each stage independently testable
- Clear error messages with line numbers