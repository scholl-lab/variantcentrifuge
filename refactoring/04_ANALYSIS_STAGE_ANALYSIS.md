# Analysis Stages Analysis

## Current State Analysis

### Code Location
- **File**: `pipeline.py`
- **Lines**: 
  - 2216-2296: Main variant analysis
  - 2410-2433: Inheritance output processing
  - Chunked processing: 271-546
- **External Files**:
  - `analyze_variants.py`: Core analysis logic
  - `inheritance/`: Inheritance pattern analysis
  - `scoring.py`: Variant scoring
  - `stats.py`: Statistics generation

### Current Functionality

The analysis phase is the most complex part of the pipeline:

#### 1. Variant Analysis
- Loads TSV into DataFrame
- Applies custom annotations
- Calculates inheritance patterns
- Applies scoring formulas
- Generates statistics

#### 2. Chunked Processing
For large files (>100MB):
- Processes data in gene-based chunks
- Maintains memory efficiency
- Handles scoring and inheritance per chunk

#### 3. Gene Burden Analysis
Optional statistical analysis:
- Case vs control comparisons
- Fisher's exact test
- Multiple testing correction

### Current Issues

1. **Memory Usage**: Loads entire file for small datasets
2. **Complex Flow**: Chunked vs regular processing paths
3. **Tight Coupling**: Analysis steps are interdependent
4. **Hard to Parallelize**: Sequential processing of analysis
5. **Testing Difficulty**: Can't test components independently

## Refactored Design

### Core Analysis Stages

```python
# stages/analysis.py

class DataFrameLoadingStage(Stage):
    """Load TSV data into DataFrame for analysis."""
    
    @property
    def name(self) -> str:
        return "dataframe_loading"
    
    @property
    def dependencies(self) -> Set[str]:
        # Depends on data processing completion
        return {"phenotype_integration"} or {"streaming_data_processing"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load data based on file size."""
        file_size_mb = context.current_file.stat().st_size / (1024 * 1024)
        
        if file_size_mb > 100 and not context.config.get('no_chunked_processing'):
            # Don't load - will use chunked processing
            context.use_chunked_processing = True
            logger.info(f"File size {file_size_mb:.1f}MB - will use chunked processing")
        else:
            # Load into DataFrame
            df = pd.read_csv(
                context.current_file,
                sep='\t',
                dtype=str,
                keep_default_na=False,
                compression='gzip' if str(context.current_file).endswith('.gz') else None
            )
            context.current_dataframe = df
            context.use_chunked_processing = False
            logger.info(f"Loaded {len(df)} variants into memory")
        
        return context


class CustomAnnotationStage(Stage):
    """Apply unified custom annotations."""
    
    @property
    def name(self) -> str:
        return "custom_annotation"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"dataframe_loading"}
    
    @property
    def can_run_parallel(self) -> bool:
        return True  # If using chunked processing
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply annotations to data."""
        custom_features = context.get_stage_result('custom_features')
        
        if not custom_features:
            # Add empty annotation column for consistency
            if context.current_dataframe is not None:
                context.current_dataframe['Custom_Annotation'] = ''
            return context
        
        if context.use_chunked_processing:
            # Annotations will be applied per chunk
            context.custom_features = custom_features
        else:
            # Apply to loaded DataFrame
            df = annotate_dataframe_with_features(
                context.current_dataframe,
                custom_features
            )
            context.current_dataframe = df
        
        return context


class InheritanceAnalysisStage(Stage):
    """Calculate Mendelian inheritance patterns."""
    
    @property
    def name(self) -> str:
        return "inheritance_analysis"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"custom_annotation"}
    
    @property
    def can_run_parallel(self) -> bool:
        return True  # Can process genes independently
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Analyze inheritance patterns."""
        if not context.config.get('calculate_inheritance'):
            return context
        
        pedigree_data = context.get_stage_result('pedigree_data')
        if not pedigree_data:
            logger.info("No pedigree data - skipping inheritance analysis")
            return context
        
        if context.use_chunked_processing:
            # Will be handled in chunked processor
            context.pedigree_data = pedigree_data
        else:
            # Analyze full DataFrame
            df = self._analyze_inheritance(
                context.current_dataframe,
                pedigree_data,
                context.config
            )
            context.current_dataframe = df
        
        return context
    
    def _analyze_inheritance(self, df: pd.DataFrame, 
                           pedigree_data: Dict[str, Any],
                           config: Dict[str, Any]) -> pd.DataFrame:
        """Perform inheritance analysis on DataFrame."""
        from ..inheritance.analyzer import analyze_inheritance_patterns
        
        # Group by gene for compound het detection
        gene_groups = df.groupby('Gene_Name', sort=False)
        
        results = []
        for gene, gene_df in gene_groups:
            analyzed_df = analyze_inheritance_patterns(
                gene_df,
                pedigree_data,
                config
            )
            results.append(analyzed_df)
        
        return pd.concat(results, ignore_index=True)


class VariantScoringStage(Stage):
    """Apply scoring formulas to variants."""
    
    @property
    def name(self) -> str:
        return "variant_scoring"
    
    @property
    def dependencies(self) -> Set[str]:
        # Scoring may use inheritance results
        deps = {"inheritance_analysis"}
        if not self._has_inheritance():
            deps = {"custom_annotation"}
        return deps
    
    @property
    def can_run_parallel(self) -> bool:
        return True  # Scoring is independent per variant
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply scoring configuration."""
        scoring_config = context.get_stage_result('scoring_config')
        if not scoring_config:
            return context
        
        if context.use_chunked_processing:
            # Store for chunk processor
            context.scoring_config = scoring_config
        else:
            # Apply to DataFrame
            df = self._apply_scoring(
                context.current_dataframe,
                scoring_config
            )
            context.current_dataframe = df
        
        return context
    
    def _apply_scoring(self, df: pd.DataFrame, 
                      scoring_config: Dict[str, Any]) -> pd.DataFrame:
        """Apply scoring formulas."""
        from ..scoring import apply_scoring_config
        
        return apply_scoring_config(df, scoring_config)


class StatisticsGenerationStage(Stage):
    """Generate summary statistics."""
    
    @property
    def name(self) -> str:
        return "statistics_generation"
    
    @property
    def dependencies(self) -> Set[str]:
        # Run after all analysis
        return {"variant_scoring"}
    
    @property
    def can_run_parallel(self) -> bool:
        return True  # Statistics are independent
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate statistics."""
        if context.config.get('no_stats'):
            return context
        
        stats_output = context.config.get('stats_output_file')
        if not stats_output:
            stats_output = context.intermediate_dir / f"{context.base_name}.statistics.tsv"
        
        if context.use_chunked_processing:
            # Aggregate statistics from chunks
            stats = self._aggregate_chunk_statistics(context)
        else:
            # Calculate from DataFrame
            stats = self._calculate_statistics(
                context.current_dataframe,
                context.config
            )
        
        # Write statistics
        self._write_statistics(stats, stats_output)
        context.statistics = stats
        
        return context
    
    def _calculate_statistics(self, df: pd.DataFrame, 
                            config: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate comprehensive statistics."""
        from ..stats import generate_variant_statistics
        
        return generate_variant_statistics(df, config)


class ChunkedAnalysisStage(Stage):
    """Process large files in memory-efficient chunks."""
    
    @property
    def name(self) -> str:
        return "chunked_analysis"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"dataframe_loading"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process file in gene-based chunks."""
        if not context.use_chunked_processing:
            return context
        
        input_file = context.current_file
        output_file = context.intermediate_dir / f"{context.base_name}.analyzed.tsv.gz"
        
        # Process chunks
        self._process_chunks(
            input_file,
            output_file,
            context,
            chunksize=context.config.get('chunk_size', 10000)
        )
        
        context.current_file = output_file
        return context
    
    def _process_chunks(self, input_file: Path, output_file: Path,
                       context: PipelineContext, chunksize: int) -> None:
        """Process file in chunks grouped by gene."""
        from ..analyze_variants import read_tsv_in_gene_chunks
        
        chunk_stats = []
        
        with smart_open(output_file, 'w') as out:
            header_written = False
            
            for chunk_df in read_tsv_in_gene_chunks(input_file, chunksize):
                # Apply all analysis steps to chunk
                if context.custom_features:
                    chunk_df = annotate_dataframe_with_features(
                        chunk_df, context.custom_features
                    )
                
                if context.pedigree_data:
                    chunk_df = self._analyze_inheritance(
                        chunk_df, context.pedigree_data, context.config
                    )
                
                if context.scoring_config:
                    chunk_df = self._apply_scoring(
                        chunk_df, context.scoring_config
                    )
                
                # Collect statistics
                if not context.config.get('no_stats'):
                    chunk_stats.append(
                        self._calculate_chunk_statistics(chunk_df)
                    )
                
                # Write chunk
                if not header_written:
                    chunk_df.to_csv(out, sep='\t', index=False, mode='w')
                    header_written = True
                else:
                    chunk_df.to_csv(out, sep='\t', index=False, mode='a', header=False)
        
        # Aggregate statistics
        if chunk_stats:
            context.statistics = self._merge_statistics(chunk_stats)
```

### Parallel Analysis for Independent Stages

```python
# stages/analysis_parallel.py

class ParallelAnalysisOrchestrator(Stage):
    """Orchestrate parallel analysis stages."""
    
    @property
    def name(self) -> str:
        return "parallel_analysis_orchestrator"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"dataframe_loading"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        """Run independent analyses in parallel."""
        if context.use_chunked_processing:
            # Chunked processing handles its own parallelization
            return ChunkedAnalysisStage()(context)
        
        # Split DataFrame by gene for parallel processing
        gene_groups = self._split_by_gene(context.current_dataframe)
        
        # Process genes in parallel
        with ThreadPoolExecutor(max_workers=context.config.get('threads', 4)) as executor:
            futures = []
            
            for gene, gene_df in gene_groups:
                future = executor.submit(
                    self._analyze_gene,
                    gene_df,
                    context
                )
                futures.append((gene, future))
            
            # Collect results
            results = []
            for gene, future in futures:
                try:
                    result_df = future.result()
                    results.append(result_df)
                except Exception as e:
                    logger.error(f"Failed to analyze gene {gene}: {e}")
                    raise
        
        # Merge results
        context.current_dataframe = pd.concat(results, ignore_index=True)
        return context
```

## Testing Strategy

### Unit Tests

```python
# tests/unit/stages/test_analysis.py

class TestAnalysisStages:
    
    def test_inheritance_analysis(self, variants_df, pedigree_data):
        """Test inheritance pattern detection."""
        context = create_test_context()
        context.current_dataframe = variants_df
        context.stage_results['pedigree_data'] = pedigree_data
        context.config['calculate_inheritance'] = True
        
        stage = InheritanceAnalysisStage()
        result = stage(context)
        
        # Check inheritance columns added
        df = result.current_dataframe
        assert 'Inheritance_Pattern' in df.columns
        assert 'Inheritance_Details' in df.columns
        
        # Verify patterns detected
        assert 'de_novo' in df['Inheritance_Pattern'].values
    
    def test_scoring_application(self, variants_df, scoring_config):
        """Test variant scoring."""
        context = create_test_context()
        context.current_dataframe = variants_df
        context.stage_results['scoring_config'] = scoring_config
        
        stage = VariantScoringStage()
        result = stage(context)
        
        # Check score columns added
        df = result.current_dataframe
        assert 'VariantScore' in df.columns
        assert df['VariantScore'].dtype == 'float64'
    
    def test_chunked_processing(self, large_tsv_file):
        """Test chunked processing for large files."""
        context = create_test_context()
        context.current_file = large_tsv_file
        context.use_chunked_processing = True
        
        stage = ChunkedAnalysisStage()
        result = stage(context)
        
        # Verify output exists and is complete
        assert result.current_file.exists()
        
        # Count lines
        with smart_open(large_tsv_file) as f:
            input_lines = sum(1 for _ in f) - 1  # Subtract header
        
        with smart_open(result.current_file) as f:
            output_lines = sum(1 for _ in f) - 1
        
        assert input_lines == output_lines
```

### Performance Tests

```python
# tests/performance/test_analysis_performance.py

def test_parallel_vs_sequential_analysis(large_variants_df):
    """Compare parallel vs sequential processing."""
    # Sequential
    context1 = create_test_context()
    context1.current_dataframe = large_variants_df.copy()
    context1.config['threads'] = 1
    
    start = time.time()
    sequential_result = run_analysis_stages(context1)
    sequential_time = time.time() - start
    
    # Parallel
    context2 = create_test_context()
    context2.current_dataframe = large_variants_df.copy()
    context2.config['threads'] = 8
    
    start = time.time()
    parallel_result = ParallelAnalysisOrchestrator()(context2)
    parallel_time = time.time() - start
    
    # Verify same results
    pd.testing.assert_frame_equal(
        sequential_result.current_dataframe.sort_values(['CHROM', 'POS']),
        parallel_result.current_dataframe.sort_values(['CHROM', 'POS'])
    )
    
    # Check speedup
    speedup = sequential_time / parallel_time
    assert speedup > 2.0  # At least 2x faster with 8 threads
```

## Optimization Strategies

### 1. Smart Chunking
- Group by gene to keep related variants together
- Dynamic chunk size based on available memory
- Preserve variant order

### 2. Parallel Processing
- Independent gene analysis
- Parallel scoring calculation
- Concurrent statistics generation

### 3. Memory Management
- Stream processing for large files
- Lazy loading of annotations
- Efficient data structures

### Performance Benchmarks

| Dataset | Variants | Current Time | Optimized Time | Memory Usage |
|---------|----------|--------------|----------------|--------------|
| Small   | 10K      | 5s           | 4s             | 100MB        |
| Medium  | 100K     | 45s          | 20s            | 500MB        |
| Large   | 1M       | 450s         | 150s           | 1GB          |
| XLarge  | 10M      | 4500s        | 900s           | 2GB (chunked)|

## Migration Plan

### Week 1: Basic Analysis Stages
1. Implement DataFrameLoadingStage
2. Implement CustomAnnotationStage
3. Test with small datasets

### Week 2: Complex Analysis
1. Implement InheritanceAnalysisStage
2. Implement VariantScoringStage
3. Test inheritance patterns

### Week 3: Performance Optimization
1. Implement ChunkedAnalysisStage
2. Implement ParallelAnalysisOrchestrator
3. Benchmark performance

### Week 4: Integration
1. Remove old analysis code
2. Update pipeline integration
3. Full regression testing

## Risk Mitigation

1. **Result Consistency**: Sort DataFrames before comparison in tests
2. **Memory Limits**: Monitor memory usage, fall back to chunked processing
3. **Parallel Bugs**: Extensive testing of thread safety
4. **Performance Regression**: Keep sequential path as baseline

## Success Metrics

- Identical analysis results (inheritance, scores, stats)
- 3-5x performance improvement for multi-core systems
- Memory usage scales with chunk size, not file size
- Each analysis component independently testable
- Clear separation of concerns