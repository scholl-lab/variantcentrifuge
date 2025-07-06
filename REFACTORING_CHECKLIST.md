# VariantCentrifuge Pipeline Refactoring Checklist

## Pre-Refactoring Setup

### Environment Preparation
- [ ] Create feature branch: `git checkout -b refactor/complete-pipeline-redesign`
- [ ] Set up test environment: `mamba activate annotation`
- [ ] Verify all tools available: `which variantcentrifuge bcftools snpEff SnpSift bedtools`
- [ ] Create baseline test outputs: `./testing/test_scripts/run_comprehensive_test.sh baseline testing/test_output_baseline`
- [ ] Generate checksums: `find testing/test_output_baseline -type f -name "*.tsv" | xargs md5sum > baseline_checksums.txt`
- [ ] Create regression test script: `testing/test_scripts/regression_test.py`
- [ ] Verify regression test works: `python testing/test_scripts/regression_test.py`

### Directory Structure
- [ ] Create pipeline package: `mkdir -p variantcentrifuge/pipeline`
- [ ] Create stages package: `mkdir -p variantcentrifuge/stages`
- [ ] Create test directories: `mkdir -p tests/unit/stages tests/integration/stages`
- [ ] Add __init__.py files to all new packages

## Phase 1: Core Infrastructure (Week 1)

### Day 1-2: Core Abstractions
- [ ] Implement `pipeline/context.py` with PipelineContext dataclass
- [ ] Run regression test ✓
- [ ] Implement `pipeline/stage.py` with Stage abstract class
- [ ] Run regression test ✓
- [ ] Implement `pipeline/runner.py` with PipelineRunner
- [ ] Run regression test ✓
- [ ] Create unit tests for core abstractions
- [ ] All tests pass ✓

### Day 3-4: Testing Infrastructure
- [ ] Set up pytest configuration for new structure
- [ ] Create mock factories in `tests/mocks/`
- [ ] Create integration test helpers
- [ ] Create performance benchmarking framework
- [ ] Document testing procedures

### Day 5: First Configuration Stage
- [ ] Implement ConfigurationLoadingStage
- [ ] Create unit tests for ConfigurationLoadingStage
- [ ] Run regression test ✓
- [ ] Integrate into pipeline with feature flag
- [ ] Test with feature flag OFF ✓
- [ ] Test with feature flag ON ✓

## Phase 2: Configuration Stages (Week 1-2)

### Configuration Loading Stages
- [ ] Implement PhenotypeLoadingStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement ScoringConfigLoadingStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement PedigreeLoadingStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement AnnotationConfigLoadingStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement SampleConfigLoadingStage
  - [ ] Unit tests
  - [ ] Regression test ✓

### Parallel Configuration Loading
- [ ] Test all config stages running in parallel
- [ ] Verify no race conditions
- [ ] Performance benchmark
- [ ] Full regression test ✓

## Phase 3: Processing Stages (Week 2)

### Gene Processing
- [ ] Implement GeneBedCreationStage
  - [ ] Extract logic from pipeline.py
  - [ ] Unit tests
  - [ ] Regression test ✓

### Variant Extraction
- [ ] Implement VariantExtractionStage (single-threaded)
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement ParallelVariantExtractionStage
  - [ ] Chunk processing logic
  - [ ] Parallel execution
  - [ ] Result merging
  - [ ] Unit tests
  - [ ] Regression test ✓
  - [ ] Performance test (verify speedup)

### Filtering Stages
- [ ] Implement BCFToolsPrefilterStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement MultiAllelicSplitStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement SnpSiftFilterStage
  - [ ] Unit tests
  - [ ] Regression test ✓

### Data Processing
- [ ] Implement FieldExtractionStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement GenotypeReplacementStage
  - [ ] Streaming implementation
  - [ ] Unit tests
  - [ ] Memory usage test
  - [ ] Regression test ✓
- [ ] Implement ExtraColumnRemovalStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement PhenotypeIntegrationStage
  - [ ] Unit tests
  - [ ] Regression test ✓

### Optimized Processing
- [ ] Implement StreamingDataProcessingStage
  - [ ] Single-pass implementation
  - [ ] Performance benchmarks
  - [ ] Memory usage tests
  - [ ] Regression test ✓

## Phase 4: Analysis Stages (Week 3)

### Data Loading
- [ ] Implement DataFrameLoadingStage
  - [ ] Size-based decision logic
  - [ ] Unit tests
  - [ ] Regression test ✓

### Analysis Components
- [ ] Implement CustomAnnotationStage
  - [ ] DataFrame mode
  - [ ] Chunked mode support
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement InheritanceAnalysisStage
  - [ ] Pedigree integration
  - [ ] Compound het detection
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement VariantScoringStage
  - [ ] Formula evaluation
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement StatisticsGenerationStage
  - [ ] Statistics calculation
  - [ ] Chunk aggregation
  - [ ] Unit tests
  - [ ] Regression test ✓

### Chunked Processing
- [ ] Implement ChunkedAnalysisStage
  - [ ] Gene-based chunking
  - [ ] Memory efficiency tests
  - [ ] Large file tests
  - [ ] Regression test ✓

### Parallel Analysis
- [ ] Implement ParallelAnalysisOrchestrator
  - [ ] Gene-level parallelization
  - [ ] Performance benchmarks
  - [ ] Regression test ✓

## Phase 5: Output Stages (Week 4)

### Post-Processing
- [ ] Implement VariantIdentifierStage
  - [ ] DataFrame mode
  - [ ] Streaming mode
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement FinalFilteringStage
  - [ ] Pandas query support
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement PseudonymizationStage
  - [ ] Sample extraction
  - [ ] Mapping generation
  - [ ] Unit tests
  - [ ] Regression test ✓

### Output Generation
- [ ] Implement TSVOutputStage
  - [ ] File output
  - [ ] Stdout support
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement MetadataGenerationStage
  - [ ] Tool version capture
  - [ ] Unit tests
  - [ ] Regression test ✓

### Report Generation
- [ ] Implement ExcelReportStage
  - [ ] Multi-sheet support
  - [ ] Formatting
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement HTMLReportStage
  - [ ] JSON generation
  - [ ] HTML generation
  - [ ] Unit tests
  - [ ] Regression test ✓
- [ ] Implement IGVReportStage
  - [ ] BAM integration
  - [ ] Local FASTA support
  - [ ] Unit tests
  - [ ] Regression test ✓

### Parallel Output
- [ ] Implement ParallelReportGenerationStage
  - [ ] Concurrent report generation
  - [ ] Performance tests
  - [ ] Regression test ✓

## Phase 6: Integration (Week 5)

### Pipeline Integration
- [ ] Update main pipeline.py to use stages
- [ ] Remove feature flags
- [ ] Remove old implementation code
- [ ] Update imports throughout codebase
- [ ] Full regression test suite ✓

### Performance Testing
- [ ] Benchmark single-gene analysis
- [ ] Benchmark 10-gene analysis (1 vs 8 cores)
- [ ] Benchmark 100-gene analysis (1 vs 8 cores)
- [ ] Benchmark large file processing (1GB+)
- [ ] Document performance improvements

### Edge Case Testing
- [ ] Test with missing optional inputs
- [ ] Test with malformed data
- [ ] Test checkpoint/resume scenarios
- [ ] Test error recovery
- [ ] Test memory limits

### Documentation
- [ ] Update docstrings for all new modules
- [ ] Update README with new architecture
- [ ] Create developer guide
- [ ] Update user documentation
- [ ] Create migration guide

## Phase 7: Final Validation

### Comprehensive Testing
- [ ] Run full test suite: `pytest tests/`
- [ ] Run integration tests: `pytest tests/integration/`
- [ ] Run performance tests: `pytest tests/performance/`
- [ ] Run regression tests with all test scenarios
- [ ] Verify all outputs identical to baseline

### Code Quality
- [ ] Run linting: `flake8 variantcentrifuge/`
- [ ] Run formatting: `black variantcentrifuge/`
- [ ] Check import sorting: `isort variantcentrifuge/`
- [ ] Verify test coverage > 85%
- [ ] Code review all changes

### Final Steps
- [ ] Create pull request
- [ ] Team code review
- [ ] Address review comments
- [ ] Final regression test ✓
- [ ] Merge to main branch
- [ ] Tag release
- [ ] Update changelog

## Post-Refactoring

### Monitoring
- [ ] Monitor for any issues in first week
- [ ] Gather performance metrics from users
- [ ] Document any lessons learned
- [ ] Plan future improvements

### Cleanup
- [ ] Remove any deprecated code
- [ ] Archive old documentation
- [ ] Update CI/CD pipelines
- [ ] Celebrate! 🎉

## Important Reminders

⚠️ **Run regression test after EVERY change**
⚠️ **Commit only working code**
⚠️ **Keep feature branch up to date with main**
⚠️ **Document any deviations from plan**
⚠️ **Ask for help if stuck**

## Progress Tracking

Start Date: ____________
Target End Date: ____________

### Weekly Progress
- Week 1: ___% complete
- Week 2: ___% complete
- Week 3: ___% complete
- Week 4: ___% complete
- Week 5: ___% complete

### Blockers
1. ________________________________
2. ________________________________
3. ________________________________

### Notes
_____________________________________
_____________________________________
_____________________________________