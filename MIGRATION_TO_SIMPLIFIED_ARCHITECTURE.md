# Migration Guide to Simplified Architecture

## Overview

This guide shows how to migrate from the current monolithic `pipeline.py` to the simplified Stage + PipelineContext architecture.

## Migration Principles

1. **One Stage per Logical Operation**: Each major function becomes a Stage
2. **Context Carries All State**: No global variables or shared state
3. **Dependencies as Strings**: Simple stage name references
4. **Complexity Hidden**: Parallel/chunked logic stays inside stages

## Step-by-Step Migration Examples

### Example 1: Migrating Configuration Loading

**Current Code (pipeline.py lines 1387-1416):**
```python
def run_pipeline(args):
    # Load configuration
    config = load_config(args.config) if args.config else {}
    
    # Merge with CLI args
    if args.reference:
        config['reference'] = args.reference
    if args.threads:
        config['threads'] = args.threads
    # ... more merging ...
```

**Migrated Stage:**
```python
class ConfigurationLoadingStage(Stage):
    @property
    def name(self) -> str:
        return "configuration_loading"
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        # Load base config
        if context.args.config:
            config = load_config(context.args.config)
        else:
            config = {}
        
        # Merge CLI args (exact same logic)
        for key, value in vars(context.args).items():
            if value is not None and key != 'config':
                config[key] = value
        
        # Store in context
        context.config = config
        return context
```

### Example 2: Migrating Parallel Variant Extraction

**Current Code (pipeline.py lines 1840-1980):**
```python
def _run_parallel_pipeline(args, normalized_genes, bed_file, config):
    # Complex parallel extraction logic
    bed_chunks = split_bed_file(bed_file, args.threads)
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for i, chunk in enumerate(bed_chunks):
            future = executor.submit(_process_bed_chunk, args, chunk, i)
            futures.append(future)
        
        results = [f.result() for f in futures]
    
    # Merge results
    merged_vcf = merge_vcfs(results)
```

**Migrated Stage:**
```python
class ParallelVariantExtractionStage(Stage):
    @property
    def name(self) -> str:
        return "variant_extraction"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"gene_bed_creation"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        # Same logic, but encapsulated
        bed_chunks = self._split_bed_file(
            context.gene_bed_file, 
            context.config.get('threads', 1)
        )
        
        with ProcessPoolExecutor(max_workers=context.config.get('threads', 1)) as executor:
            futures = []
            for i, chunk in enumerate(bed_chunks):
                future = executor.submit(
                    self._process_bed_chunk,
                    context.args.vcf_file,
                    chunk,
                    context.workspace.get_temp_path(f"chunk_{i}.vcf")
                )
                futures.append(future)
            
            results = [f.result() for f in futures]
        
        # Merge and store result
        merged_vcf = self._merge_vcfs(
            results,
            context.workspace.get_intermediate_path("variants.vcf.gz")
        )
        
        context.extracted_vcf = merged_vcf
        context.data = merged_vcf  # Update primary data artifact
        return context
    
    def _split_bed_file(self, bed_file: Path, n_chunks: int) -> List[Path]:
        # Move the split logic here
        pass
    
    def _process_bed_chunk(self, vcf: Path, bed: Path, output: Path) -> Path:
        # Move the chunk processing here
        pass
    
    def _merge_vcfs(self, vcf_files: List[Path], output: Path) -> Path:
        # Move the merge logic here
        pass
```

### Example 3: Migrating Analysis with DataFrame

**Current Code (pipeline.py lines 2342-2445):**
```python
# Inheritance analysis
if args.ped_file:
    logger.info("Calculating inheritance patterns...")
    
    # Load data
    df = pd.read_csv(current_file, sep='\t', dtype=str)
    
    # Apply inheritance
    df = apply_inheritance_analysis(df, ped_data, config)
    
    # Save
    output_file = os.path.join(intermediate_dir, f"{base_name}.inheritance.tsv")
    df.to_csv(output_file, sep='\t', index=False)
```

**Migrated Stage:**
```python
class InheritanceAnalysisStage(Stage):
    @property
    def name(self) -> str:
        return "inheritance_analysis"
    
    @property
    def dependencies(self) -> Set[str]:
        return {"dataframe_loading", "pedigree_loading"}
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        if not context.pedigree_data:
            # Skip if no pedigree
            return context
        
        logger.info("Calculating inheritance patterns...")
        
        # Work with DataFrame from context
        df = context.current_dataframe
        
        # Apply same analysis logic
        df = apply_inheritance_analysis(
            df, 
            context.pedigree_data, 
            context.config
        )
        
        # Update context
        context.current_dataframe = df
        
        # Also save to file if needed
        if context.config.get('save_intermediates'):
            output_file = context.workspace.get_intermediate_path(
                f"{context.workspace.base_name}.inheritance.tsv"
            )
            df.to_csv(output_file, sep='\t', index=False)
        
        return context
```

## Common Migration Patterns

### 1. Function → Stage

```python
# Before
def extract_fields_from_vcf(vcf_file, fields, output_file):
    # ... implementation ...

# After
class FieldExtractionStage(Stage):
    name = "field_extraction"
    dependencies = {"variant_filtering"}
    
    def _process(self, context):
        # Same implementation, but using context
        self._extract_fields(
            context.filtered_vcf,
            context.config['fields_to_extract'],
            context.workspace.get_intermediate_path('extracted.tsv')
        )
```

### 2. Global State → Context Attributes

```python
# Before
samples = []  # Global
pedigree_data = None  # Global

# After
# Everything in context
context.vcf_samples = []
context.pedigree_data = None
```

### 3. File Paths → Workspace

```python
# Before
output_file = os.path.join(output_dir, f"{base_name}.filtered.vcf")

# After
output_file = context.workspace.get_output_path(".filtered", ".vcf")
```

### 4. Conditional Logic → Optional Stages

```python
# Before
if args.scoring_config:
    apply_scoring(df, scoring_config)

# After (in main pipeline)
stages = [
    # ... other stages ...
    VariantScoringStage() if args.scoring_config else None,
]
stages = [s for s in stages if s is not None]
```

## Testing During Migration

For each migrated stage:

1. **Create Unit Test First**
```python
def test_new_stage():
    # Test the stage in isolation
    context = create_test_context()
    stage = NewStage()
    result = stage(context)
    assert_expected_behavior(result)
```

2. **Run Side-by-Side Comparison**
```python
def test_stage_matches_original():
    # Run original function
    original_output = original_function(test_input)
    
    # Run new stage
    context = create_context_from_input(test_input)
    stage = NewStage()
    result = stage(context)
    
    # Compare outputs
    assert files_are_identical(original_output, result.data)
```

3. **Integration Test**
```python
def test_stage_in_pipeline():
    # Run with previous and next stages
    stages = [PreviousStage(), NewStage(), NextStage()]
    runner = PipelineRunner()
    result = runner.run(stages, initial_context)
    assert_pipeline_works(result)
```

## Gradual Migration Strategy

### Phase 1: Core Infrastructure
1. Create `pipeline/` package structure
2. Implement PipelineContext, Stage, Workspace, Runner
3. Test infrastructure with dummy stages

### Phase 2: Independent Stages
Start with stages that have no dependencies:
- ConfigurationLoadingStage
- PhenotypeLoadingStage  
- ScoringConfigLoadingStage

### Phase 3: Processing Pipeline
Migrate in dependency order:
- GeneBedCreationStage
- VariantExtractionStage
- FilteringStages
- FieldExtractionStage

### Phase 4: Analysis Pipeline
- DataFrameLoadingStage
- All analysis stages
- Keep chunked logic internal

### Phase 5: Output Pipeline
- All output stages
- Parallel report generation

### Phase 6: Cleanup
- Remove old pipeline.py code
- Update imports
- Final testing

## Validation Checklist

For each migrated stage:
- [ ] All logic preserved exactly
- [ ] Unit tests pass
- [ ] Integration tests pass
- [ ] Output identical to original
- [ ] Performance equal or better
- [ ] Code under 200 lines
- [ ] Dependencies declared correctly

## Benefits After Migration

1. **pipeline.py**: 2,831 lines → ~100 lines
2. **Each stage**: Self-contained, testable unit
3. **Parallel execution**: Automatic based on dependencies
4. **Easy debugging**: Clear stage boundaries
5. **Simple testing**: Consistent patterns

The migration preserves all functionality while dramatically improving maintainability.