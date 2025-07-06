# Configuration and Setup Stage Analysis

## Current State Analysis

### Code Location
- **File**: `pipeline.py`
- **Lines**: 1387-1516 (configuration loading)
- **Functions**: Multiple inline configuration blocks in `run_pipeline`

### Current Functionality

The configuration loading currently handles:

1. **Phenotype Configuration** (lines 1392-1411)
   - Load phenotype file with sample-value mappings
   - Store in config dict
   
2. **Scoring Configuration** (lines 1412-1423)
   - Load scoring JSON configurations
   - Validate scoring formulas
   
3. **Pedigree Data** (lines 1424-1452)
   - Load PED file if provided
   - Create singleton pedigrees if not
   
4. **Custom Annotations** (lines 1443-1452)
   - Validate annotation configuration
   - Load BED files, gene lists, JSON gene data
   
5. **Phenotype Terms** (lines 1453-1466)
   - Parse case/control phenotype terms
   - Load from files or command line
   
6. **Sample Lists** (lines 1467-1480)
   - Load case/control sample lists
   - Parse from files or command line

### Dependencies
- `load_phenotypes()` from phenotype.py
- `read_scoring_config()` from scoring.py
- `read_pedigree()` from ped_reader.py
- `validate_annotation_config()` and `load_custom_features()` from annotator.py
- `load_terms_from_file()` - local helper function

### Issues with Current Implementation

1. **Mixed Concerns**: Configuration loading is intertwined with validation and setup
2. **No Error Recovery**: One failed config stops entire pipeline
3. **Scattered Logic**: Configuration happens throughout the function
4. **Hard to Test**: Can't test configuration loading in isolation
5. **No Parallelization**: Everything loads sequentially

## Refactored Design

### Stage Breakdown

We'll create multiple focused stages that can run in parallel where possible:

```python
# stages/setup.py

class ConfigurationLoadingStage(Stage):
    """Load and validate base configuration."""
    
    @property
    def name(self) -> str:
        return "config_loading"
    
    @property
    def can_run_parallel(self) -> bool:
        return True  # No dependencies
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        # Just configuration merging and validation
        # No file I/O here
        return context


class PhenotypeLoadingStage(Stage):
    """Load phenotype data from files."""
    
    @property
    def name(self) -> str:
        return "phenotype_loading"
    
    @property
    def can_run_parallel(self) -> bool:
        return True  # Can run parallel with other loaders
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        if not context.args.phenotype_file:
            return context
            
        phenotypes = load_phenotypes(
            context.args.phenotype_file,
            context.args.phenotype_sample_column,
            context.args.phenotype_value_column
        )
        
        context.stage_results['phenotypes'] = phenotypes
        return context


class ScoringConfigLoadingStage(Stage):
    """Load variant scoring configuration."""
    
    @property  
    def name(self) -> str:
        return "scoring_config_loading"
    
    @property
    def can_run_parallel(self) -> bool:
        return True
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        if not context.config.get('scoring_config_path'):
            return context
            
        scoring_config = read_scoring_config(
            context.config['scoring_config_path']
        )
        
        context.stage_results['scoring_config'] = scoring_config
        return context


class PedigreeLoadingStage(Stage):
    """Load pedigree data for inheritance analysis."""
    
    @property
    def name(self) -> str:
        return "pedigree_loading"
    
    @property
    def can_run_parallel(self) -> bool:
        return True
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        if not context.config.get('calculate_inheritance'):
            return context
            
        if context.config.get('ped_file'):
            pedigree_data = read_pedigree(context.config['ped_file'])
        else:
            # Will populate later when we have samples
            pedigree_data = {}
            
        context.stage_results['pedigree_data'] = pedigree_data
        return context


class AnnotationConfigLoadingStage(Stage):
    """Load and validate custom annotation sources."""
    
    @property
    def name(self) -> str:
        return "annotation_config_loading"
    
    @property
    def can_run_parallel(self) -> bool:
        return True
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        # Validate configuration
        errors = validate_annotation_config(context.config)
        if errors:
            raise ValueError(f"Annotation config errors: {errors}")
            
        # Load features
        custom_features = load_custom_features(context.config)
        context.stage_results['custom_features'] = custom_features
        return context


class SampleConfigLoadingStage(Stage):
    """Load case/control sample configurations."""
    
    @property
    def name(self) -> str:
        return "sample_config_loading"
    
    @property
    def can_run_parallel(self) -> bool:
        return True
    
    def _process(self, context: PipelineContext) -> PipelineContext:
        # Load case samples
        case_samples = []
        if context.args.case_samples:
            case_samples = [s.strip() for s in context.args.case_samples.split(",")]
        case_samples += load_terms_from_file(context.args.case_samples_file, logger)
        
        # Load control samples
        control_samples = []
        if context.args.control_samples:
            control_samples = [s.strip() for s in context.args.control_samples.split(",")]
        control_samples += load_terms_from_file(context.args.control_samples_file, logger)
        
        # Load phenotype terms
        case_terms = []
        if context.args.case_phenotypes:
            case_terms = [t.strip() for t in context.args.case_phenotypes.split(",")]
        case_terms += load_terms_from_file(context.args.case_phenotypes_file, logger)
        
        control_terms = []
        if context.args.control_phenotypes:
            control_terms = [t.strip() for t in context.args.control_phenotypes.split(",")]
        control_terms += load_terms_from_file(context.args.control_phenotypes_file, logger)
        
        context.stage_results['case_samples'] = case_samples
        context.stage_results['control_samples'] = control_samples
        context.stage_results['case_phenotypes'] = case_terms
        context.stage_results['control_phenotypes'] = control_terms
        
        return context
```

### Parallel Execution Plan

These stages can all run in parallel since they have no dependencies:

```
┌─────────────────────┐  ┌──────────────────────┐  ┌─────────────────────┐
│ ConfigurationLoading│  │  PhenotypeLoading    │  │  ScoringConfigLoading│
└─────────────────────┘  └──────────────────────┘  └─────────────────────┘
         │                         │                          │
         └─────────────────────────┴──────────────────────────┘
                                   │
┌─────────────────────┐  ┌──────────────────────┐  ┌─────────────────────┐
│  PedigreeLoading    │  │AnnotationConfigLoading│ │ SampleConfigLoading │
└─────────────────────┘  └──────────────────────┘  └─────────────────────┘
```

## Testing Strategy

### Unit Tests

```python
# tests/unit/stages/test_setup.py

class TestConfigurationStages:
    
    def test_phenotype_loading_stage(self, tmp_path):
        # Create test phenotype file
        phenotype_file = tmp_path / "phenotypes.txt"
        phenotype_file.write_text("sample1\taffected\nsample2\tcontrol\n")
        
        # Create context
        context = create_test_context(
            phenotype_file=str(phenotype_file),
            phenotype_sample_column="0",
            phenotype_value_column="1"
        )
        
        # Run stage
        stage = PhenotypeLoadingStage()
        result = stage(context)
        
        # Verify
        phenotypes = result.get_stage_result('phenotypes')
        assert phenotypes['sample1'] == 'affected'
        assert phenotypes['sample2'] == 'control'
    
    def test_parallel_loading(self):
        # Test that multiple stages can run concurrently
        context = create_test_context()
        
        stages = [
            PhenotypeLoadingStage(),
            ScoringConfigLoadingStage(),
            PedigreeLoadingStage(),
            AnnotationConfigLoadingStage(),
            SampleConfigLoadingStage()
        ]
        
        runner = ParallelPipelineRunner(max_parallel=5)
        result = runner.run(stages, context)
        
        # All stages should complete
        assert all(result.is_stage_complete(s.name) for s in stages)
```

### Integration Tests

```python
# tests/integration/test_configuration_loading.py

def test_full_configuration_loading(test_data_dir):
    """Test loading all configuration types together."""
    args = create_test_args(
        phenotype_file=test_data_dir / "phenotypes.txt",
        scoring_config_path=test_data_dir / "scoring_config",
        ped_file=test_data_dir / "family.ped",
        annotate_bed=[test_data_dir / "regions.bed"],
        case_samples_file=test_data_dir / "cases.txt"
    )
    
    context = PipelineContext(
        config=load_config(),
        args=args,
        start_time=datetime.now(),
        output_dir=Path("test_output"),
        intermediate_dir=Path("test_output/intermediate"),
        base_name="test"
    )
    
    # Run all config stages
    stages = [
        ConfigurationLoadingStage(),
        PhenotypeLoadingStage(),
        ScoringConfigLoadingStage(),
        PedigreeLoadingStage(),
        AnnotationConfigLoadingStage(),
        SampleConfigLoadingStage()
    ]
    
    runner = ParallelPipelineRunner()
    result = runner.run(stages, context)
    
    # Verify all loaded correctly
    assert result.get_stage_result('phenotypes') is not None
    assert result.get_stage_result('scoring_config') is not None
    assert result.get_stage_result('pedigree_data') is not None
```

## Migration Steps

1. **Create stage files**
   ```bash
   mkdir -p variantcentrifuge/stages
   touch variantcentrifuge/stages/__init__.py
   touch variantcentrifuge/stages/setup.py
   ```

2. **Implement stages one by one**
   - Start with ConfigurationLoadingStage (no file I/O)
   - Test in isolation
   - Add remaining stages

3. **Update pipeline.py**
   - Remove configuration loading code
   - Add stages to pipeline

4. **Test thoroughly**
   ```bash
   # After each stage implementation
   python testing/test_scripts/regression_test.py
   ```

## Rollback Plan

If any stage fails:
1. The entire configuration loading can be reverted by removing stage imports
2. Original inline code remains until all stages are proven
3. Git reset to last working commit if needed

## Success Metrics

- Configuration loading time reduced by 30-50% through parallelization
- Each stage independently testable
- Clear error messages for configuration issues
- No change to output files (regression tests pass)