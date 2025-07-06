# Testing Strategy with Simplified Architecture

## Overview

The simplified Stage + PipelineContext architecture makes testing significantly easier and more consistent. This document shows how the testing strategy adapts to the new design.

## Core Testing Pattern

With the unified architecture, every stage follows the same testing pattern:

```python
def test_any_stage():
    # 1. Create context with required state
    context = PipelineContext(
        args=Mock(relevant_args=value),
        config={'relevant_config': 'value'},
        workspace=Mock()
    )
    
    # 2. Set up prerequisite state
    context.data = "path/to/input.file"
    context.mark_complete('dependency_stage')
    
    # 3. Execute stage
    stage = StageUnderTest()
    result = stage(context)
    
    # 4. Assert on returned context
    assert result.data == "path/to/output.file"
    assert result.is_complete('stage_under_test')
```

## Simplified Test Examples

### 1. Testing Configuration Stage

```python
def test_configuration_loading():
    # Given
    args = Mock(
        config='config.json',
        threads=8,
        reference=None
    )
    context = PipelineContext(args=args, config={}, workspace=Mock())
    
    # When
    stage = ConfigurationLoadingStage()
    result = stage(context)
    
    # Then
    assert result.config['threads'] == 8  # CLI override
    assert result.config['reference'] == 'GRCh37'  # From config file
```

### 2. Testing Processing Stage with Mocked External Tool

```python
def test_variant_extraction():
    # Given
    context = PipelineContext(
        args=Mock(vcf_file='input.vcf'),
        config={'threads': 4},
        workspace=Mock(get_temp_path=lambda x: f"/tmp/{x}")
    )
    context.gene_bed_file = Path("genes.bed")
    context.mark_complete('gene_bed_creation')
    
    # When
    with patch('subprocess.run') as mock_run:
        mock_run.return_value.returncode = 0
        
        stage = ParallelVariantExtractionStage()
        result = stage(context)
    
    # Then
    assert result.data == result.extracted_vcf
    assert mock_run.call_count >= 1  # Called for extraction
```

### 3. Testing Analysis Stage with DataFrame

```python
def test_inheritance_analysis():
    # Given
    test_df = pd.DataFrame({
        'CHROM': ['chr1'],
        'POS': [100],
        'GT': ['child(0/1);father(0/0);mother(0/0)']  # De novo
    })
    
    context = PipelineContext(
        args=Mock(),
        config={'inheritance_mode': 'simple'},
        workspace=Mock()
    )
    context.current_dataframe = test_df
    context.mark_complete('dataframe_loading')
    
    # When
    stage = InheritanceAnalysisStage()
    result = stage(context)
    
    # Then
    assert 'Inheritance_Pattern' in result.current_dataframe.columns
    assert result.current_dataframe.iloc[0]['Inheritance_Pattern'] == 'de_novo'
```

## Integration Testing Simplified

Testing stage sequences becomes trivial:

```python
def test_processing_pipeline():
    # Create initial context
    context = PipelineContext(
        args=test_args,
        config=test_config,
        workspace=Workspace(tmp_path, "test")
    )
    
    # Run stages in sequence
    stages = [
        GeneBedCreationStage(),
        VariantExtractionStage(),
        FieldExtractionStage(),
        GenotypeReplacementStage()
    ]
    
    for stage in stages:
        context = stage(context)
    
    # Assert final state
    assert context.data.endswith('.genotype_replaced.tsv')
    assert Path(context.data).exists()
```

## Mocking Strategies

### 1. Mock External Tools

```python
@pytest.fixture
def mock_bcftools():
    with patch('subprocess.run') as mock:
        mock.return_value.returncode = 0
        mock.return_value.stdout = "##fileformat=VCFv4.2\n"
        yield mock
```

### 2. Mock File System

```python
@pytest.fixture
def mock_workspace(tmp_path):
    workspace = Mock()
    workspace.output_dir = tmp_path / "output"
    workspace.intermediate_dir = tmp_path / "intermediate"
    workspace.get_temp_path = lambda x: tmp_path / "temp" / x
    return workspace
```

### 3. Mock Complex Dependencies

```python
@pytest.fixture
def context_with_complete_analysis(mock_workspace):
    context = PipelineContext(
        args=Mock(),
        config=full_config,
        workspace=mock_workspace
    )
    
    # Mark all analysis stages complete
    for stage in ['dataframe_loading', 'custom_annotation', 
                  'inheritance_analysis', 'variant_scoring']:
        context.mark_complete(stage)
    
    # Add analysis results
    context.current_dataframe = analyzed_dataframe
    context.statistics = computed_stats
    
    return context
```

## Testing Benefits with Simplified Architecture

### 1. Consistent Pattern
Every stage test follows the same pattern:
- Create context
- Set prerequisites
- Execute stage
- Assert results

### 2. Easy Dependency Management
```python
# Dependencies are just strings to mark complete
context.mark_complete('dependency1')
context.mark_complete('dependency2')
```

### 3. Isolated Testing
Each stage can be tested completely in isolation:
```python
# No need to run previous stages
context.data = "mocked/previous/output.tsv"
context.mark_complete('all_previous_stages')
```

### 4. Simple Integration Tests
```python
# Just run stages in sequence
runner = PipelineRunner()
result = runner.run([stage1, stage2, stage3], initial_context)
```

## Parallel Execution Testing

The simplified architecture makes parallel testing straightforward:

```python
def test_parallel_report_generation():
    # Given
    context = create_context_with_output_ready()
    
    # When
    with patch('concurrent.futures.ThreadPoolExecutor') as mock_executor:
        stage = ParallelReportGenerationStage()
        result = stage(context)
    
    # Then
    assert mock_executor.called
    assert all(report in result.report_paths 
              for report in ['excel', 'html', 'igv'])
```

## Performance Testing Pattern

```python
def test_stage_performance(benchmark):
    # Create context with large dataset
    context = create_large_dataset_context()
    
    # Benchmark stage execution
    stage = ChunkedAnalysisStage()
    result = benchmark(stage, context)
    
    # Assert performance
    assert benchmark.stats['mean'] < 5.0  # seconds
```

## Summary

The simplified architecture reduces testing complexity by:

1. **Unified Interface**: All stages test the same way
2. **Explicit Dependencies**: Just mark stages complete
3. **Clear Data Flow**: Context in → Context out
4. **Easy Mocking**: Standard patterns for all mocks
5. **Simple Integration**: Chain stages naturally

This makes the comprehensive test suite of 376 tests much easier to implement and maintain.