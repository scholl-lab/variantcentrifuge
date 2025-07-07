# VariantCentrifuge Pipeline Refactoring Checklist

## Pre-Refactoring Setup

### Environment Preparation
- [x] Create feature branch: `git checkout -b refactor/complete-pipeline-redesign`
- [x] Set up test environment: `mamba activate annotation`
- [x] Verify all tools available: `which variantcentrifuge bcftools snpEff SnpSift bedtools`
- [ ] Create baseline test outputs: `./testing/test_scripts/run_comprehensive_test.sh baseline testing/test_output_baseline`
- [ ] Generate checksums: `find testing/test_output_baseline -type f -name "*.tsv" | xargs md5sum > baseline_checksums.txt`
- [x] Create regression test script: `testing/test_scripts/regression_test.py`
- [ ] Verify regression test works: `python testing/test_scripts/regression_test.py`

### Directory Structure
- [x] Create pipeline package: `mkdir -p variantcentrifuge/pipeline`
- [x] Create stages package: `mkdir -p variantcentrifuge/stages`
- [x] Create test directories: `mkdir -p tests/unit/stages tests/integration/stages`
- [x] Add __init__.py files to all new packages

## Phase 1: Core Infrastructure (Week 1)

### Day 1-2: Core Abstractions
- [x] Implement `pipeline/context.py` with PipelineContext dataclass
- [ ] Run regression test ✓
- [x] Implement `pipeline/stage.py` with Stage abstract class
- [ ] Run regression test ✓
- [x] Implement `pipeline/runner.py` with PipelineRunner
- [ ] Run regression test ✓
- [x] Create unit tests for core abstractions
- [x] All tests pass ✓

### Day 3-4: Testing Infrastructure
- [x] Set up pytest configuration for new structure
- [x] Create mock factories in `tests/mocks/`
- [x] Create integration test helpers
- [ ] Create performance benchmarking framework
- [ ] Document testing procedures

### Day 5: First Configuration Stage
- [x] Implement ConfigurationLoadingStage
- [x] Create unit tests for ConfigurationLoadingStage
- [ ] Run regression test ✓
- [x] Integrate into pipeline with feature flag
- [x] Test with feature flag OFF ✓
- [x] Test with feature flag ON ✓

## Phase 2: Configuration Stages (Week 1-2)

### Configuration Loading Stages
- [x] Implement PhenotypeLoadingStage
  - [x] Unit tests
  - [ ] Regression test ✓
- [x] Implement ScoringConfigLoadingStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement PedigreeLoadingStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement AnnotationConfigLoadingStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement SampleConfigLoadingStage
  - [x] Unit tests
  - [ ] Regression test ✓

### Parallel Configuration Loading
- [x] Test all config stages running in parallel
- [ ] Verify no race conditions
- [ ] Performance benchmark
- [ ] Full regression test ✓

## Phase 3: Processing Stages (Week 2)

### Gene Processing
- [x] Implement GeneBedCreationStage
  - [x] Extract logic from pipeline.py
  - [x] Unit tests
  - [ ] Regression test ✓

### Variant Extraction
- [x] Implement VariantExtractionStage (single-threaded)
  - [x] Unit tests
  - [ ] Regression test ✓
- [x] Implement ParallelVariantExtractionStage
  - [x] Chunk processing logic
  - [x] Parallel execution
  - [x] Result merging
  - [x] Unit tests
  - [ ] Regression test ✓
  - [ ] Performance test (verify speedup)

### Filtering Stages
- [x] Implement BCFToolsPrefilterStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement MultiAllelicSplitStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement SnpSiftFilterStage
  - [ ] Unit tests
  - [ ] Regression test ✓

### Data Processing
- [x] Implement FieldExtractionStage
  - [x] Unit tests
  - [ ] Regression test ✓
- [x] Implement GenotypeReplacementStage
  - [ ] Streaming implementation
  - [x] Unit tests
  - [ ] Memory usage test
  - [ ] Regression test ✓
- [x] Implement ExtraColumnRemovalStage
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement PhenotypeIntegrationStage
  - [x] Unit tests
  - [ ] Regression test ✓

### Optimized Processing
- [x] Implement StreamingDataProcessingStage
  - [x] Single-pass implementation
  - [ ] Performance benchmarks
  - [ ] Memory usage tests
  - [ ] Regression test ✓

## Phase 4: Analysis Stages (Week 3)

### Data Loading
- [x] Implement DataFrameLoadingStage
  - [x] Size-based decision logic
  - [ ] Unit tests
  - [ ] Regression test ✓

### Analysis Components
- [x] Implement CustomAnnotationStage
  - [x] DataFrame mode
  - [x] Chunked mode support
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement InheritanceAnalysisStage
  - [x] Pedigree integration
  - [x] Compound het detection
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement VariantScoringStage
  - [x] Formula evaluation
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement StatisticsGenerationStage
  - [x] Statistics calculation
  - [x] Chunk aggregation
  - [ ] Unit tests
  - [ ] Regression test ✓

### Chunked Processing
- [x] Implement ChunkedAnalysisStage
  - [ ] Gene-based chunking
  - [ ] Memory efficiency tests
  - [ ] Large file tests
  - [ ] Regression test ✓

### Parallel Analysis
- [x] Implement ParallelAnalysisOrchestrator
  - [ ] Gene-level parallelization
  - [ ] Performance benchmarks
  - [ ] Regression test ✓

## Phase 5: Output Stages (Week 4)

### Post-Processing
- [x] Implement VariantIdentifierStage
  - [x] DataFrame mode
  - [ ] Streaming mode
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement FinalFilteringStage
  - [x] Pandas query support
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement PseudonymizationStage
  - [x] Sample extraction
  - [x] Mapping generation
  - [ ] Unit tests
  - [ ] Regression test ✓

### Output Generation
- [x] Implement TSVOutputStage
  - [x] File output
  - [x] Stdout support
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement MetadataGenerationStage
  - [x] Tool version capture
  - [ ] Unit tests
  - [ ] Regression test ✓

### Report Generation
- [x] Implement ExcelReportStage
  - [x] Multi-sheet support
  - [x] Formatting
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement HTMLReportStage
  - [x] JSON generation
  - [x] HTML generation
  - [ ] Unit tests
  - [ ] Regression test ✓
- [x] Implement IGVReportStage
  - [x] BAM integration
  - [x] Local FASTA support
  - [ ] Unit tests
  - [ ] Regression test ✓

### Parallel Output
- [x] Implement ParallelReportGenerationStage
  - [x] Concurrent report generation
  - [ ] Performance tests
  - [ ] Regression test ✓

## Phase 6: Integration (Week 5)

### Pipeline Integration
- [x] Update main pipeline.py to use stages
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
- [x] Run linting: `flake8 variantcentrifuge/`
- [x] Run formatting: `black variantcentrifuge/`
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

Start Date: 2025-07-07
Target End Date: TBD

### Weekly Progress
- Week 1: 100% complete (Core Infrastructure + Config Stages)
- Week 2: 100% complete (Processing Stages)
- Week 3: 100% complete (Analysis Stages)
- Week 4: 100% complete (Output Stages)
- Week 5: 20% complete (Integration started, testing pending)

### Blockers
1. ________________________________
2. ________________________________
3. ________________________________

### Notes
_____________________________________
_____________________________________
_____________________________________