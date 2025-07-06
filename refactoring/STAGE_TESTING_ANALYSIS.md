# Comprehensive Testing Analysis for Pipeline Refactoring

## Overview

This document provides a detailed testing analysis for each stage in the pipeline refactoring, identifying existing tests that can be reused, tests that need adaptation, and new tests that must be created.

## Testing Categories

1. **Reusable Tests**: Existing tests that can be used with minimal changes
2. **Adaptable Tests**: Tests that need modification to work with the new architecture
3. **New Tests Required**: Tests that must be created from scratch
4. **Integration Points**: Tests that verify stage interactions

---

## Stage 1: Configuration Loading Stages

### 1.1 ConfigurationLoadingStage

**Current Code Location**: `pipeline.py` lines 1387-1516

**Existing Tests**:
- `test_cli.py`: Tests configuration loading from CLI args
- No specific unit tests for configuration merging

**Tests to Adapt**:
```python
# From test_cli.py - needs adaptation for new stage
def test_config_loading():
    # OLD: Tests args passed to run_pipeline
    # NEW: Test ConfigurationLoadingStage directly
    
    # Adapt to:
    context = PipelineContext(config={}, args=args, ...)
    stage = ConfigurationLoadingStage()
    result = stage(context)
    assert result.config['reference'] == 'GRCh37'
```

**New Tests Required**:
```python
# tests/unit/stages/test_configuration_loading.py

def test_config_merging_priority():
    """Test that CLI args override config file settings."""
    base_config = {'threads': 1, 'reference': 'hg19'}
    args = Namespace(threads=8, reference=None)
    
    context = PipelineContext(config=base_config, args=args)
    stage = ConfigurationLoadingStage()
    result = stage(context)
    
    assert result.config['threads'] == 8  # CLI override
    assert result.config['reference'] == 'hg19'  # Config retained

def test_default_values_applied():
    """Test default values are set when not provided."""
    context = PipelineContext(config={}, args=Namespace())
    stage = ConfigurationLoadingStage()
    result = stage(context)
    
    assert 'fields_to_extract' in result.config
    assert 'no_stats' in result.config
    assert result.config['no_stats'] is False

def test_path_resolution():
    """Test relative paths are resolved to absolute."""
    args = Namespace(output_dir='./output')
    context = PipelineContext(config={}, args=args)
    stage = ConfigurationLoadingStage()
    result = stage(context)
    
    assert Path(result.output_dir).is_absolute()
```

### 1.2 PhenotypeLoadingStage

**Current Code Location**: `pipeline.py` lines 1392-1411

**Existing Tests**:
- No direct tests for phenotype loading
- Integration tests in `test_pipeline_integration.py` include phenotypes

**Tests to Adapt**:
None directly applicable

**New Tests Required**:
```python
# tests/unit/stages/test_phenotype_loading.py

def test_phenotype_file_loading(tmp_path):
    """Test loading phenotypes from TSV file."""
    phenotype_file = tmp_path / "phenotypes.tsv"
    phenotype_file.write_text("sample1\taffected\nsample2\tcontrol\n")
    
    args = Namespace(
        phenotype_file=str(phenotype_file),
        phenotype_sample_column='0',
        phenotype_value_column='1'
    )
    context = PipelineContext(config={}, args=args)
    stage = PhenotypeLoadingStage()
    result = stage(context)
    
    phenotypes = result.get_stage_result('phenotypes')
    assert phenotypes['sample1'] == 'affected'
    assert phenotypes['sample2'] == 'control'

def test_phenotype_missing_file():
    """Test graceful handling of missing phenotype file."""
    args = Namespace(phenotype_file='nonexistent.tsv')
    context = PipelineContext(config={}, args=args)
    stage = PhenotypeLoadingStage()
    
    with pytest.raises(FileNotFoundError):
        stage(context)

def test_phenotype_malformed_file(tmp_path):
    """Test error handling for malformed phenotype file."""
    phenotype_file = tmp_path / "bad.tsv"
    phenotype_file.write_text("only_one_column\n")
    
    args = Namespace(
        phenotype_file=str(phenotype_file),
        phenotype_sample_column='0',
        phenotype_value_column='1'
    )
    context = PipelineContext(config={}, args=args)
    stage = PhenotypeLoadingStage()
    
    with pytest.raises(ValueError, match="Invalid phenotype file format"):
        stage(context)

def test_no_phenotype_file():
    """Test stage completes when no phenotype file provided."""
    args = Namespace(phenotype_file=None)
    context = PipelineContext(config={}, args=args)
    stage = PhenotypeLoadingStage()
    result = stage(context)
    
    assert result.get_stage_result('phenotypes') is None
```

### 1.3 ScoringConfigLoadingStage

**Current Code Location**: `pipeline.py` lines 1412-1423

**Existing Tests**:
- `test_scoring.py`: Tests scoring configuration loading
- `test_scoring_integration.py`: Integration tests

**Tests to Adapt**:
```python
# From test_scoring.py
def test_read_scoring_config():
    # OLD: Direct function call
    config = read_scoring_config(config_path)
    
    # NEW: Through stage
    context = PipelineContext(
        config={'scoring_config_path': config_path},
        args=Namespace()
    )
    stage = ScoringConfigLoadingStage()
    result = stage(context)
    config = result.get_stage_result('scoring_config')
```

**New Tests Required**:
```python
# tests/unit/stages/test_scoring_config_loading.py

def test_scoring_config_validation():
    """Test validation of scoring configuration."""
    # Test missing required fields
    bad_config_path = create_bad_scoring_config()
    context = PipelineContext(
        config={'scoring_config_path': bad_config_path},
        args=Namespace()
    )
    stage = ScoringConfigLoadingStage()
    
    with pytest.raises(ValueError, match="Invalid scoring configuration"):
        stage(context)

def test_scoring_config_caching():
    """Test that scoring config is cached in stage results."""
    context = PipelineContext(
        config={'scoring_config_path': 'scoring/test_config'},
        args=Namespace()
    )
    stage = ScoringConfigLoadingStage()
    result = stage(context)
    
    config = result.get_stage_result('scoring_config')
    assert 'variable_mapping' in config
    assert 'formulas' in config
```

### 1.4 PedigreeLoadingStage

**Current Code Location**: `pipeline.py` lines 1424-1452

**Existing Tests**:
- `test_ped_reader.py`: Tests PED file parsing
- `test_inheritance/`: Extensive inheritance tests

**Tests to Adapt**:
```python
# From test_ped_reader.py
def test_read_pedigree():
    # OLD: Direct function test
    pedigree = read_pedigree(ped_file)
    
    # NEW: Through stage
    context = PipelineContext(
        config={'ped_file': ped_file, 'calculate_inheritance': True},
        args=Namespace()
    )
    stage = PedigreeLoadingStage()
    result = stage(context)
    pedigree = result.get_stage_result('pedigree_data')
```

**New Tests Required**:
```python
# tests/unit/stages/test_pedigree_loading.py

def test_pedigree_without_inheritance_calculation():
    """Test pedigree not loaded when inheritance not requested."""
    context = PipelineContext(
        config={'ped_file': 'family.ped', 'calculate_inheritance': False},
        args=Namespace()
    )
    stage = PedigreeLoadingStage()
    result = stage(context)
    
    assert result.get_stage_result('pedigree_data') is None

def test_singleton_pedigree_creation():
    """Test creation of singleton pedigrees when no PED file."""
    context = PipelineContext(
        config={'calculate_inheritance': True},
        args=Namespace()
    )
    stage = PedigreeLoadingStage()
    result = stage(context)
    
    # Should create empty dict, populated later
    assert result.get_stage_result('pedigree_data') == {}

def test_malformed_ped_file(tmp_path):
    """Test error handling for malformed PED file."""
    bad_ped = tmp_path / "bad.ped"
    bad_ped.write_text("only three columns\n")
    
    context = PipelineContext(
        config={'ped_file': str(bad_ped), 'calculate_inheritance': True},
        args=Namespace()
    )
    stage = PedigreeLoadingStage()
    
    with pytest.raises(ValueError, match="Invalid PED file"):
        stage(context)
```

### 1.5 AnnotationConfigLoadingStage

**Current Code Location**: `pipeline.py` lines 1443-1452

**Existing Tests**:
- `test_annotator.py`: Tests annotation configuration
- `test_gene_list_annotation.py`: Gene list annotation tests

**Tests to Adapt**:
```python
# From test_annotator.py
def test_validate_annotation_config():
    # OLD: Direct validation
    errors = validate_annotation_config(cfg)
    
    # NEW: Through stage (validation happens in stage)
    context = PipelineContext(config=cfg, args=Namespace())
    stage = AnnotationConfigLoadingStage()
    # Should raise if invalid
    with pytest.raises(ValueError):
        stage(context)
```

**New Tests Required**:
```python
# tests/unit/stages/test_annotation_config_loading.py

def test_bed_file_loading(tmp_path):
    """Test loading BED file annotations."""
    bed_file = tmp_path / "regions.bed"
    bed_file.write_text("chr1\t1000\t2000\tRegion1\n")
    
    context = PipelineContext(
        config={'annotate_bed': [str(bed_file)]},
        args=Namespace()
    )
    stage = AnnotationConfigLoadingStage()
    result = stage(context)
    
    features = result.get_stage_result('custom_features')
    assert 'regions_by_chrom' in features
    assert 'chr1' in features['regions_by_chrom']

def test_gene_list_loading(tmp_path):
    """Test loading gene list annotations."""
    gene_list = tmp_path / "genes.txt"
    gene_list.write_text("BRCA1\nTP53\nEGFR\n")
    
    context = PipelineContext(
        config={'annotate_gene_list': [str(gene_list)]},
        args=Namespace()
    )
    stage = AnnotationConfigLoadingStage()
    result = stage(context)
    
    features = result.get_stage_result('custom_features')
    assert 'gene_lists' in features
    assert 'BRCA1' in features['gene_lists']['genes']

def test_json_gene_data_loading(tmp_path):
    """Test loading JSON gene annotations."""
    json_file = tmp_path / "genes.json"
    json_file.write_text('{"BRCA1": {"panel": "cancer", "inheritance": "AD"}}')
    
    context = PipelineContext(
        config={
            'annotate_json_genes': str(json_file),
            'json_gene_mapping': {'identifier': 'symbol', 'dataFields': ['panel']}
        },
        args=Namespace()
    )
    stage = AnnotationConfigLoadingStage()
    result = stage(context)
    
    features = result.get_stage_result('custom_features')
    assert 'json_gene_data' in features
    assert features['json_gene_data']['BRCA1']['panel'] == 'cancer'

def test_multiple_annotation_sources(tmp_path):
    """Test loading multiple annotation types together."""
    # Create all annotation types
    bed_file = tmp_path / "regions.bed"
    bed_file.write_text("chr1\t1000\t2000\tRegion1\n")
    
    gene_list = tmp_path / "genes.txt"
    gene_list.write_text("BRCA1\n")
    
    context = PipelineContext(
        config={
            'annotate_bed': [str(bed_file)],
            'annotate_gene_list': [str(gene_list)]
        },
        args=Namespace()
    )
    stage = AnnotationConfigLoadingStage()
    result = stage(context)
    
    features = result.get_stage_result('custom_features')
    assert 'regions_by_chrom' in features
    assert 'gene_lists' in features
```

### 1.6 SampleConfigLoadingStage

**Current Code Location**: `pipeline.py` lines 1467-1480

**Existing Tests**:
- Limited direct tests
- Some coverage in integration tests

**Tests to Adapt**:
None directly applicable

**New Tests Required**:
```python
# tests/unit/stages/test_sample_config_loading.py

def test_case_control_from_cli():
    """Test loading case/control samples from CLI args."""
    args = Namespace(
        case_samples='Sample1,Sample2',
        control_samples='Sample3,Sample4',
        case_samples_file=None,
        control_samples_file=None,
        case_phenotypes='HP:0001234,HP:0005678',
        control_phenotypes=None,
        case_phenotypes_file=None,
        control_phenotypes_file=None
    )
    context = PipelineContext(config={}, args=args)
    stage = SampleConfigLoadingStage()
    result = stage(context)
    
    assert result.get_stage_result('case_samples') == ['Sample1', 'Sample2']
    assert result.get_stage_result('control_samples') == ['Sample3', 'Sample4']
    assert result.get_stage_result('case_phenotypes') == ['HP:0001234', 'HP:0005678']

def test_sample_lists_from_files(tmp_path):
    """Test loading sample lists from files."""
    case_file = tmp_path / "cases.txt"
    case_file.write_text("Sample1\nSample2\n")
    
    control_file = tmp_path / "controls.txt"
    control_file.write_text("Sample3\nSample4\n")
    
    args = Namespace(
        case_samples=None,
        control_samples=None,
        case_samples_file=str(case_file),
        control_samples_file=str(control_file),
        case_phenotypes=None,
        control_phenotypes=None,
        case_phenotypes_file=None,
        control_phenotypes_file=None
    )
    context = PipelineContext(config={}, args=args)
    stage = SampleConfigLoadingStage()
    result = stage(context)
    
    assert result.get_stage_result('case_samples') == ['Sample1', 'Sample2']
    assert result.get_stage_result('control_samples') == ['Sample3', 'Sample4']

def test_combined_cli_and_file_loading(tmp_path):
    """Test combining CLI and file-based sample lists."""
    case_file = tmp_path / "cases.txt"
    case_file.write_text("Sample3\nSample4\n")
    
    args = Namespace(
        case_samples='Sample1,Sample2',
        control_samples=None,
        case_samples_file=str(case_file),
        control_samples_file=None,
        case_phenotypes=None,
        control_phenotypes=None,
        case_phenotypes_file=None,
        control_phenotypes_file=None
    )
    context = PipelineContext(config={}, args=args)
    stage = SampleConfigLoadingStage()
    result = stage(context)
    
    # Should combine both sources
    assert result.get_stage_result('case_samples') == ['Sample1', 'Sample2', 'Sample3', 'Sample4']

def test_empty_sample_lists():
    """Test handling of empty sample lists."""
    args = Namespace(
        case_samples=None,
        control_samples=None,
        case_samples_file=None,
        control_samples_file=None,
        case_phenotypes=None,
        control_phenotypes=None,
        case_phenotypes_file=None,
        control_phenotypes_file=None
    )
    context = PipelineContext(config={}, args=args)
    stage = SampleConfigLoadingStage()
    result = stage(context)
    
    assert result.get_stage_result('case_samples') == []
    assert result.get_stage_result('control_samples') == []
```

---

## Stage 2: Processing Stages

### 2.1 GeneBedCreationStage

**Current Code Location**: `pipeline.py` lines 1563-1580, `gene_bed.py`

**Existing Tests**:
- No direct tests for BED creation
- Integration tests verify BED files created

**Tests to Adapt**:
None directly applicable

**New Tests Required**:
```python
# tests/unit/stages/test_gene_bed_creation.py

def test_single_gene_bed_creation(tmp_path):
    """Test BED creation for single gene."""
    args = Namespace(gene_name='BRCA1', gene_file=None)
    context = PipelineContext(
        config={'reference': 'GRCh37'},
        args=args,
        output_dir=tmp_path,
        base_name='test'
    )
    stage = GeneBedCreationStage()
    result = stage(context)
    
    assert result.bed_file.exists()
    assert result.normalized_genes == 'BRCA1'

def test_gene_list_bed_creation(tmp_path):
    """Test BED creation from gene list file."""
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("BRCA1\nTP53\nEGFR\n")
    
    args = Namespace(gene_name=None, gene_file=str(gene_file))
    context = PipelineContext(
        config={'reference': 'GRCh37'},
        args=args,
        output_dir=tmp_path,
        base_name='test'
    )
    stage = GeneBedCreationStage()
    result = stage(context)
    
    assert result.bed_file.exists()
    assert result.normalized_genes == 'BRCA1,EGFR,TP53'  # Alphabetized

def test_interval_expansion(tmp_path):
    """Test BED creation with interval expansion."""
    args = Namespace(gene_name='BRCA1', gene_file=None)
    context = PipelineContext(
        config={'reference': 'GRCh37', 'interval_expansion': 1000},
        args=args,
        output_dir=tmp_path,
        base_name='test'
    )
    stage = GeneBedCreationStage()
    result = stage(context)
    
    # Verify intervals are expanded
    with open(result.bed_file) as f:
        # Would need to check actual coordinates
        lines = f.readlines()
        assert len(lines) > 0

def test_missing_gene_error():
    """Test error handling for non-existent gene."""
    args = Namespace(gene_name='INVALID_GENE_XYZ', gene_file=None)
    context = PipelineContext(
        config={'reference': 'GRCh37'},
        args=args,
        output_dir=Path('/tmp'),
        base_name='test'
    )
    stage = GeneBedCreationStage()
    
    with pytest.raises(RuntimeError, match="No intervals found"):
        stage(context)

def test_empty_gene_list():
    """Test error handling for empty gene list."""
    args = Namespace(gene_name='', gene_file=None)
    context = PipelineContext(
        config={'reference': 'GRCh37'},
        args=args,
        output_dir=Path('/tmp'),
        base_name='test'
    )
    stage = GeneBedCreationStage()
    
    with pytest.raises(ValueError, match="No genes specified"):
        stage(context)
```

### 2.2 VariantExtractionStage / ParallelVariantExtractionStage

**Current Code Location**: `pipeline.py` lines 831-944, 944-1216

**Existing Tests**:
- `test_pipeline_integration.py`: Tests variant extraction
- `test_chunked_processing.py`: Tests chunked processing

**Tests to Adapt**:
```python
# From test_pipeline_integration.py
def test_variant_extraction():
    # OLD: Full pipeline test
    # NEW: Test extraction stage directly
    
    context = PipelineContext(
        config={'threads': 1},
        args=Namespace(vcf_file='test.vcf'),
        bed_file=Path('genes.bed')
    )
    stage = VariantExtractionStage()
    result = stage(context)
    
    assert result.variants_vcf.exists()
```

**New Tests Required**:
```python
# tests/unit/stages/test_variant_extraction.py

def test_single_region_extraction(test_vcf, test_bed):
    """Test basic variant extraction for single region."""
    context = PipelineContext(
        config={},
        args=Namespace(vcf_file=test_vcf),
        bed_file=test_bed,
        intermediate_dir=Path('/tmp/intermediate'),
        base_name='test'
    )
    stage = VariantExtractionStage()
    result = stage(context)
    
    assert result.variants_vcf.exists()
    assert count_variants(result.variants_vcf) > 0

def test_parallel_extraction_consistency(test_vcf, multi_region_bed):
    """Test parallel extraction produces same results as single-threaded."""
    # Single-threaded
    context1 = PipelineContext(
        config={'threads': 1},
        args=Namespace(vcf_file=test_vcf),
        bed_file=multi_region_bed,
        intermediate_dir=Path('/tmp/intermediate1'),
        base_name='test'
    )
    single_result = VariantExtractionStage()(context1)
    
    # Multi-threaded
    context2 = PipelineContext(
        config={'threads': 4},
        args=Namespace(vcf_file=test_vcf),
        bed_file=multi_region_bed,
        intermediate_dir=Path('/tmp/intermediate2'),
        base_name='test'
    )
    parallel_result = ParallelVariantExtractionStage()(context2)
    
    # Compare variant counts
    assert count_variants(single_result.variants_vcf) == count_variants(parallel_result.extracted_tsv)

def test_chunk_processing_order_preservation(test_vcf, multi_region_bed):
    """Test that parallel processing preserves variant order."""
    context = PipelineContext(
        config={'threads': 4},
        args=Namespace(vcf_file=test_vcf),
        bed_file=multi_region_bed,
        intermediate_dir=Path('/tmp/intermediate'),
        base_name='test'
    )
    
    stage = ParallelVariantExtractionStage()
    # Mock the chunk processing to verify order
    original_process = stage._process_single_chunk
    chunk_order = []
    
    def mock_process(chunk, ctx):
        chunk_order.append(chunk.name)
        return original_process(chunk, ctx)
    
    stage._process_single_chunk = mock_process
    result = stage(context)
    
    # Verify chunks processed (order may vary)
    assert len(chunk_order) > 1
    
    # Verify output maintains correct order
    variants = read_variants(result.extracted_tsv)
    assert variants_are_sorted(variants)

def test_extraction_with_missing_vcf():
    """Test error handling for missing VCF file."""
    context = PipelineContext(
        config={},
        args=Namespace(vcf_file='nonexistent.vcf'),
        bed_file=Path('genes.bed')
    )
    stage = VariantExtractionStage()
    
    with pytest.raises(FileNotFoundError):
        stage(context)

def test_extraction_with_empty_bed(test_vcf, tmp_path):
    """Test extraction with empty BED file."""
    empty_bed = tmp_path / "empty.bed"
    empty_bed.touch()
    
    context = PipelineContext(
        config={},
        args=Namespace(vcf_file=test_vcf),
        bed_file=empty_bed
    )
    stage = VariantExtractionStage()
    result = stage(context)
    
    assert result.variants_vcf.exists()
    assert count_variants(result.variants_vcf) == 0
```

### 2.3 BCFToolsPrefilterStage

**Current Code Location**: Uses `apply_bcftools_filter` from filters.py

**Existing Tests**:
- `test_filters.py`: Has some bcftools tests
- Integration tests use bcftools filtering

**Tests to Adapt**:
```python
# From test_filters.py
def test_bcftools_filter():
    # OLD: Direct function call
    result = apply_bcftools_filter(vcf, filter_expr, output)
    
    # NEW: Through stage
    context = PipelineContext(
        config={'bcftools_filter': 'QUAL>20'},
        variants_vcf=vcf_path
    )
    stage = BCFToolsPrefilterStage()
    result = stage(context)
```

**New Tests Required**:
```python
# tests/unit/stages/test_bcftools_prefilter.py

def test_bcftools_quality_filter(test_variants_vcf):
    """Test filtering by quality score."""
    context = PipelineContext(
        config={'bcftools_filter': 'QUAL>30'},
        variants_vcf=test_variants_vcf,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = BCFToolsPrefilterStage()
    result = stage(context)
    
    assert result.variants_vcf != test_variants_vcf  # New file created
    assert count_variants(result.variants_vcf) < count_variants(test_variants_vcf)

def test_bcftools_complex_filter(test_variants_vcf):
    """Test complex bcftools filter expression."""
    context = PipelineContext(
        config={'bcftools_filter': '(QUAL>30) && (DP>10) && (AF<0.01)'},
        variants_vcf=test_variants_vcf,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = BCFToolsPrefilterStage()
    result = stage(context)
    
    # Verify filter applied
    variants = read_variants(result.variants_vcf)
    for v in variants:
        assert v.qual > 30
        assert v.info['DP'] > 10
        assert v.info['AF'] < 0.01

def test_no_bcftools_filter():
    """Test stage skips when no filter configured."""
    context = PipelineContext(
        config={},  # No bcftools_filter
        variants_vcf=Path('test.vcf')
    )
    stage = BCFToolsPrefilterStage()
    result = stage(context)
    
    assert result.variants_vcf == Path('test.vcf')  # Unchanged

def test_invalid_bcftools_expression(test_variants_vcf):
    """Test error handling for invalid filter expression."""
    context = PipelineContext(
        config={'bcftools_filter': 'INVALID SYNTAX [['},
        variants_vcf=test_variants_vcf
    )
    stage = BCFToolsPrefilterStage()
    
    with pytest.raises(RuntimeError, match="bcftools filter failed"):
        stage(context)
```

### 2.4 FieldExtractionStage

**Current Code Location**: Uses `extract_fields` from extractor.py

**Existing Tests**:
- Limited direct tests
- Covered in integration tests

**Tests to Adapt**:
None directly applicable

**New Tests Required**:
```python
# tests/unit/stages/test_field_extraction.py

def test_basic_field_extraction(test_filtered_vcf):
    """Test extraction of basic fields."""
    context = PipelineContext(
        config={'fields_to_extract': ['CHROM', 'POS', 'REF', 'ALT', 'QUAL']},
        filtered_vcf=test_filtered_vcf,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = FieldExtractionStage()
    result = stage(context)
    
    assert result.extracted_tsv.exists()
    
    # Verify TSV structure
    df = pd.read_csv(result.extracted_tsv, sep='\t')
    assert list(df.columns) == ['CHROM', 'POS', 'REF', 'ALT', 'QUAL']

def test_info_field_extraction(test_vcf_with_annotations):
    """Test extraction of INFO fields."""
    context = PipelineContext(
        config={'fields_to_extract': ['CHROM', 'POS', 'ANN[*].GENE', 'ANN[*].IMPACT']},
        filtered_vcf=test_vcf_with_annotations,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = FieldExtractionStage()
    result = stage(context)
    
    df = pd.read_csv(result.extracted_tsv, sep='\t')
    assert 'ANN[*].GENE' in df.columns
    assert 'ANN[*].IMPACT' in df.columns

def test_sample_field_extraction(test_multisample_vcf):
    """Test extraction of sample-specific fields."""
    context = PipelineContext(
        config={
            'fields_to_extract': ['CHROM', 'POS', 'GEN[*].GT', 'GEN[*].DP'],
            'append_extra_sample_fields': True,
            'extra_sample_fields': ['GQ']
        },
        filtered_vcf=test_multisample_vcf,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = FieldExtractionStage()
    result = stage(context)
    
    df = pd.read_csv(result.extracted_tsv, sep='\t')
    assert 'GEN[*].GT' in df.columns
    assert 'GEN[*].DP' in df.columns
    assert 'GEN[*].GQ' in df.columns  # Extra field added

def test_missing_field_handling(test_vcf):
    """Test handling of non-existent fields."""
    context = PipelineContext(
        config={'fields_to_extract': ['CHROM', 'POS', 'NONEXISTENT_FIELD']},
        filtered_vcf=test_vcf,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = FieldExtractionStage()
    result = stage(context)
    
    df = pd.read_csv(result.extracted_tsv, sep='\t')
    assert 'NONEXISTENT_FIELD' in df.columns
    assert df['NONEXISTENT_FIELD'].isna().all()  # All NA values
```

### 2.5 GenotypeReplacementStage

**Current Code Location**: `pipeline.py` lines 1976-2017, uses replacer.py

**Existing Tests**:
- `test_replacer.py`: Comprehensive genotype replacement tests
- `test_replacer_missing_allele.py`: Edge case tests

**Tests to Adapt**:
```python
# From test_replacer.py
def test_replace_genotypes():
    # OLD: Direct function call
    result = replace_genotypes(input_file, output_file, sample_names)
    
    # NEW: Through stage
    context = PipelineContext(
        samples=['Sample1', 'Sample2'],
        current_file=input_tsv,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = GenotypeReplacementStage()
    result = stage(context)
```

**New Tests Required**:
```python
# tests/unit/stages/test_genotype_replacement.py

def test_streaming_genotype_replacement(large_tsv_with_genotypes):
    """Test memory-efficient streaming replacement."""
    context = PipelineContext(
        config={'no_replacement': False},
        samples=['Sample1', 'Sample2', 'Sample3'],
        current_file=large_tsv_with_genotypes,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    
    # Monitor memory usage
    import tracemalloc
    tracemalloc.start()
    
    stage = GenotypeReplacementStage()
    result = stage(context)
    
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    # Should use minimal memory regardless of file size
    assert peak / 1024 / 1024 < 100  # Less than 100MB
    
    # Verify replacements
    df = pd.read_csv(result.current_file, sep='\t')
    assert 'Sample1(0/1)' in df['GT'].iloc[0]

def test_sample_substring_removal():
    """Test removal of substring from sample names."""
    context = PipelineContext(
        config={
            'no_replacement': False,
            'remove_sample_substring': '_DNA'
        },
        samples=['Sample1_DNA', 'Sample2_DNA'],
        current_file=create_test_tsv_with_gt(),
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = GenotypeReplacementStage()
    result = stage(context)
    
    df = pd.read_csv(result.current_file, sep='\t')
    assert 'Sample1(0/1)' in df['GT'].iloc[0]  # _DNA removed
    assert 'Sample1_DNA' not in df['GT'].iloc[0]

def test_no_gt_column_handling(tsv_without_gt):
    """Test handling of TSV without GT column."""
    context = PipelineContext(
        config={'no_replacement': False},
        samples=['Sample1', 'Sample2'],
        current_file=tsv_without_gt,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = GenotypeReplacementStage()
    result = stage(context)
    
    # Should return unchanged
    assert result.current_file == tsv_without_gt

def test_malformed_genotype_handling():
    """Test handling of malformed genotype values."""
    malformed_tsv = create_tsv_with_malformed_gt()
    context = PipelineContext(
        config={'no_replacement': False},
        samples=['Sample1', 'Sample2'],
        current_file=malformed_tsv,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = GenotypeReplacementStage()
    result = stage(context)
    
    # Should handle gracefully
    df = pd.read_csv(result.current_file, sep='\t')
    # Malformed entries should be preserved or marked
```

### 2.6 PhenotypeIntegrationStage

**Current Code Location**: `pipeline.py` lines 2081-2147

**Existing Tests**:
- Limited direct tests
- Some integration test coverage

**Tests to Adapt**:
None directly applicable

**New Tests Required**:
```python
# tests/unit/stages/test_phenotype_integration.py

def test_phenotype_column_addition(tsv_with_replaced_genotypes):
    """Test adding phenotype column based on samples."""
    phenotypes = {
        'Sample1': 'affected',
        'Sample2': 'control',
        'Sample3': 'unknown'
    }
    
    context = PipelineContext(
        config={},
        current_file=tsv_with_replaced_genotypes,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    context.stage_results['phenotypes'] = phenotypes
    
    stage = PhenotypeIntegrationStage()
    result = stage(context)
    
    df = pd.read_csv(result.current_file, sep='\t')
    assert 'phenotypes' in df.columns
    assert df.columns[-1] == 'phenotypes'  # Added at end
    
    # Verify phenotype mapping
    assert 'affected' in df['phenotypes'].iloc[0]

def test_multi_sample_phenotype_aggregation():
    """Test phenotype aggregation for multi-sample variants."""
    phenotypes = {
        'Sample1': 'affected',
        'Sample2': 'affected',
        'Sample3': 'control'
    }
    
    # Create TSV with multi-sample variant
    tsv_content = "CHROM\tPOS\tGT\nchr1\t100\tSample1(0/1);Sample2(0/1);Sample3(0/0)\n"
    test_file = create_temp_file(tsv_content)
    
    context = PipelineContext(
        config={},
        current_file=test_file,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    context.stage_results['phenotypes'] = phenotypes
    
    stage = PhenotypeIntegrationStage()
    result = stage(context)
    
    df = pd.read_csv(result.current_file, sep='\t')
    # Should aggregate phenotypes
    assert 'affected,affected,control' in df['phenotypes'].iloc[0]

def test_no_phenotype_data():
    """Test stage completes when no phenotype data available."""
    context = PipelineContext(
        config={},
        current_file=Path('test.tsv'),
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    # No phenotypes in stage_results
    
    stage = PhenotypeIntegrationStage()
    result = stage(context)
    
    # Should return unchanged
    assert result.current_file == Path('test.tsv')

def test_complex_genotype_format_parsing():
    """Test parsing of complex genotype formats."""
    phenotypes = {'Sample1': 'affected'}
    
    # Test various GT formats
    tsv_content = """CHROM\tPOS\tGT
chr1\t100\tSample1(0/1:20:99)
chr1\t200\tSample1(1|1)
chr1\t300\tSample1(.)
chr1\t400\t
"""
    test_file = create_temp_file(tsv_content)
    
    context = PipelineContext(
        config={},
        current_file=test_file,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    context.stage_results['phenotypes'] = phenotypes
    
    stage = PhenotypeIntegrationStage()
    result = stage(context)
    
    df = pd.read_csv(result.current_file, sep='\t')
    assert df['phenotypes'].iloc[0] == 'affected'
    assert df['phenotypes'].iloc[1] == 'affected'
    assert df['phenotypes'].iloc[2] == 'affected'
    assert df['phenotypes'].iloc[3] == ''  # Empty GT
```

---

## Stage 3: Analysis Stages

### 3.1 DataFrameLoadingStage

**Current Code Location**: `pipeline.py` lines 2216-2232

**Existing Tests**:
- `test_chunked_processing.py`: Tests size-based decisions

**Tests to Adapt**:
```python
# From test_chunked_processing.py
def test_chunked_vs_memory_decision():
    # OLD: Embedded in pipeline
    # NEW: Test stage decision logic
    
    # Small file - load to memory
    small_context = PipelineContext(
        config={'no_chunked_processing': False},
        current_file=small_file  # <100MB
    )
    stage = DataFrameLoadingStage()
    result = stage(small_context)
    assert result.current_dataframe is not None
    assert not result.use_chunked_processing
```

**New Tests Required**:
```python
# tests/unit/stages/test_dataframe_loading.py

def test_size_based_loading_decision():
    """Test decision between memory and chunked processing."""
    # Create files of different sizes
    small_file = create_tsv_file(size_mb=50)
    large_file = create_tsv_file(size_mb=150)
    
    # Small file - should load to memory
    context1 = PipelineContext(
        config={},
        current_file=small_file,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    stage = DataFrameLoadingStage()
    result1 = stage(context1)
    
    assert result1.current_dataframe is not None
    assert not result1.use_chunked_processing
    
    # Large file - should use chunked
    context2 = PipelineContext(
        config={},
        current_file=large_file,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    result2 = stage(context2)
    
    assert result2.current_dataframe is None
    assert result2.use_chunked_processing

def test_force_no_chunked_processing():
    """Test forcing memory loading regardless of size."""
    large_file = create_tsv_file(size_mb=150)
    
    context = PipelineContext(
        config={'no_chunked_processing': True},
        current_file=large_file
    )
    stage = DataFrameLoadingStage()
    result = stage(context)
    
    assert result.current_dataframe is not None
    assert not result.use_chunked_processing

def test_compressed_file_loading():
    """Test loading gzipped TSV files."""
    compressed_file = create_compressed_tsv()
    
    context = PipelineContext(
        config={},
        current_file=compressed_file
    )
    stage = DataFrameLoadingStage()
    result = stage(context)
    
    assert result.current_dataframe is not None
    assert len(result.current_dataframe) > 0

def test_malformed_tsv_handling():
    """Test error handling for malformed TSV."""
    bad_tsv = create_malformed_tsv()
    
    context = PipelineContext(
        config={},
        current_file=bad_tsv
    )
    stage = DataFrameLoadingStage()
    
    with pytest.raises(pd.errors.ParserError):
        stage(context)
```

### 3.2 CustomAnnotationStage

**Current Code Location**: Uses annotator.py

**Existing Tests**:
- `test_annotator.py`: Comprehensive annotation tests
- `test_gene_list_annotation.py`: Gene list specific tests

**Tests to Adapt**:
```python
# From test_annotator.py
def test_unified_annotation():
    # OLD: Direct function call
    df = annotate_dataframe_with_features(df, features)
    
    # NEW: Through stage
    context = PipelineContext(
        current_dataframe=df
    )
    context.stage_results['custom_features'] = features
    stage = CustomAnnotationStage()
    result = stage(context)
```

**New Tests Required**:
```python
# tests/unit/stages/test_custom_annotation.py

def test_annotation_with_all_feature_types():
    """Test applying all annotation types together."""
    df = create_test_variants_df()
    
    features = {
        'regions_by_chrom': {
            'chr1': [{'start': 1000, 'end': 2000, 'name': 'Region1'}]
        },
        'gene_lists': {
            'cancer_genes': {'BRCA1', 'TP53'}
        },
        'json_gene_data': {
            'BRCA1': {'panel': 'hereditary', 'inheritance': 'AD'}
        }
    }
    
    context = PipelineContext(
        current_dataframe=df,
        use_chunked_processing=False
    )
    context.stage_results['custom_features'] = features
    
    stage = CustomAnnotationStage()
    result = stage(context)
    
    assert 'Custom_Annotation' in result.current_dataframe.columns
    
    # Check various annotations
    annotations = result.current_dataframe['Custom_Annotation']
    assert any('Region=Region1' in str(a) for a in annotations)
    assert any('InGeneList=cancer_genes' in str(a) for a in annotations)
    assert any('panel=hereditary' in str(a) for a in annotations)

def test_chunked_annotation_setup():
    """Test annotation setup for chunked processing."""
    features = {'regions_by_chrom': {'chr1': []}}
    
    context = PipelineContext(
        current_dataframe=None,
        use_chunked_processing=True
    )
    context.stage_results['custom_features'] = features
    
    stage = CustomAnnotationStage()
    result = stage(context)
    
    # Should store features for chunk processor
    assert result.custom_features == features

def test_empty_annotation_column():
    """Test adding empty annotation column when no features."""
    df = create_test_variants_df()
    
    context = PipelineContext(
        current_dataframe=df,
        use_chunked_processing=False
    )
    # No custom_features in stage_results
    
    stage = CustomAnnotationStage()
    result = stage(context)
    
    assert 'Custom_Annotation' in result.current_dataframe.columns
    assert result.current_dataframe['Custom_Annotation'].eq('').all()

def test_annotation_performance(large_variants_df):
    """Test annotation performance with large dataset."""
    features = create_complex_features()
    
    context = PipelineContext(
        current_dataframe=large_variants_df,
        use_chunked_processing=False
    )
    context.stage_results['custom_features'] = features
    
    import time
    start = time.time()
    
    stage = CustomAnnotationStage()
    result = stage(context)
    
    duration = time.time() - start
    
    # Should complete in reasonable time
    variants_per_second = len(large_variants_df) / duration
    assert variants_per_second > 1000  # At least 1000 variants/second
```

### 3.3 InheritanceAnalysisStage

**Current Code Location**: Uses inheritance/ module

**Existing Tests**:
- `test_inheritance/`: Extensive inheritance testing
- Multiple integration tests

**Tests to Adapt**:
```python
# From test_inheritance/test_analyzer.py
def test_inheritance_analysis():
    # OLD: Direct module call
    df = analyze_inheritance_patterns(df, pedigree, config)
    
    # NEW: Through stage
    context = PipelineContext(
        config={'calculate_inheritance': True},
        current_dataframe=df
    )
    context.stage_results['pedigree_data'] = pedigree
    stage = InheritanceAnalysisStage()
    result = stage(context)
```

**New Tests Required**:
```python
# tests/unit/stages/test_inheritance_analysis.py

def test_inheritance_with_trio_data():
    """Test inheritance analysis with trio (parents + child)."""
    df = create_trio_variants_df()
    pedigree = create_trio_pedigree()
    
    context = PipelineContext(
        config={
            'calculate_inheritance': True,
            'inheritance_mode': 'simple'
        },
        current_dataframe=df,
        use_chunked_processing=False
    )
    context.stage_results['pedigree_data'] = pedigree
    
    stage = InheritanceAnalysisStage()
    result = stage(context)
    
    assert 'Inheritance_Pattern' in result.current_dataframe.columns
    assert 'Inheritance_Details' in result.current_dataframe.columns
    
    # Check for de novo detection
    patterns = result.current_dataframe['Inheritance_Pattern']
    assert 'de_novo' in patterns.values

def test_compound_heterozygous_detection():
    """Test detection of compound heterozygous variants."""
    df = create_comp_het_variants_df()  # Two variants in same gene
    pedigree = create_trio_pedigree()
    
    context = PipelineContext(
        config={
            'calculate_inheritance': True,
            'no_vectorized_comp_het': False  # Use optimized version
        },
        current_dataframe=df,
        use_chunked_processing=False
    )
    context.stage_results['pedigree_data'] = pedigree
    
    stage = InheritanceAnalysisStage()
    result = stage(context)
    
    patterns = result.current_dataframe['Inheritance_Pattern']
    assert 'compound_heterozygous' in patterns.values

def test_inheritance_without_pedigree():
    """Test handling when no pedigree data available."""
    df = create_test_variants_df()
    
    context = PipelineContext(
        config={'calculate_inheritance': True},
        current_dataframe=df,
        use_chunked_processing=False
    )
    # No pedigree_data in stage_results
    
    stage = InheritanceAnalysisStage()
    result = stage(context)
    
    # Should skip analysis
    assert result.current_dataframe is df  # Unchanged

def test_inheritance_mode_columns():
    """Test different inheritance output modes."""
    df = create_trio_variants_df()
    pedigree = create_trio_pedigree()
    
    # Test columns mode
    context = PipelineContext(
        config={
            'calculate_inheritance': True,
            'inheritance_mode': 'columns'
        },
        current_dataframe=df,
        use_chunked_processing=False
    )
    context.stage_results['pedigree_data'] = pedigree
    
    stage = InheritanceAnalysisStage()
    result = stage(context)
    
    # Should have expanded columns
    assert 'Inheritance_Pattern' in result.current_dataframe.columns
    assert 'Inheritance_Confidence' in result.current_dataframe.columns
    assert 'Inheritance_Description' in result.current_dataframe.columns
    assert 'Inheritance_Samples' in result.current_dataframe.columns

def test_gene_based_inheritance_grouping():
    """Test that inheritance is analyzed per gene."""
    # Create variants in different genes
    df = pd.DataFrame({
        'Gene_Name': ['GENE1', 'GENE1', 'GENE2', 'GENE2'],
        'CHROM': ['chr1'] * 4,
        'POS': [100, 200, 300, 400],
        'GT': ['child(0/1);mother(0/1);father(0/0)'] * 4
    })
    
    pedigree = create_trio_pedigree()
    
    context = PipelineContext(
        config={'calculate_inheritance': True},
        current_dataframe=df,
        use_chunked_processing=False
    )
    context.stage_results['pedigree_data'] = pedigree
    
    # Mock the analysis to track gene grouping
    analysis_calls = []
    original_analyze = InheritanceAnalysisStage._analyze_inheritance
    
    def mock_analyze(self, df, pedigree, config):
        analysis_calls.append(df['Gene_Name'].unique().tolist())
        return original_analyze(self, df, pedigree, config)
    
    InheritanceAnalysisStage._analyze_inheritance = mock_analyze
    
    stage = InheritanceAnalysisStage()
    result = stage(context)
    
    # Should be called once per gene
    assert len(analysis_calls) == 2
    assert ['GENE1'] in analysis_calls
    assert ['GENE2'] in analysis_calls
```

### 3.4 VariantScoringStage

**Current Code Location**: Uses scoring.py

**Existing Tests**:
- `test_scoring.py`: Basic scoring tests
- `test_scoring_integration.py`: Full pipeline scoring
- `test_scoring_inheritance.py`: Scoring with inheritance

**Tests to Adapt**:
```python
# From test_scoring.py
def test_apply_scoring():
    # OLD: Direct function
    df = apply_scoring_config(df, scoring_config)
    
    # NEW: Through stage
    context = PipelineContext(
        current_dataframe=df
    )
    context.stage_results['scoring_config'] = scoring_config
    stage = VariantScoringStage()
    result = stage(context)
```

**New Tests Required**:
```python
# tests/unit/stages/test_variant_scoring.py

def test_basic_scoring_formula():
    """Test application of basic scoring formula."""
    df = pd.DataFrame({
        'CADD_PHRED': [20.0, 30.0, 10.0],
        'gnomAD_AF': [0.001, 0.0001, 0.01]
    })
    
    scoring_config = {
        'variable_mapping': {
            'cadd': {'field': 'CADD_PHRED', 'default': 0},
            'af': {'field': 'gnomAD_AF', 'default': 1}
        },
        'formulas': [
            {
                'name': 'SimpleScore',
                'formula': 'cadd * (1 - af)'
            }
        ]
    }
    
    context = PipelineContext(
        current_dataframe=df,
        use_chunked_processing=False
    )
    context.stage_results['scoring_config'] = scoring_config
    
    stage = VariantScoringStage()
    result = stage(context)
    
    assert 'SimpleScore' in result.current_dataframe.columns
    
    # Verify calculations
    scores = result.current_dataframe['SimpleScore']
    assert abs(scores.iloc[0] - 19.98) < 0.01  # 20 * (1 - 0.001)
    assert abs(scores.iloc[1] - 29.997) < 0.01  # 30 * (1 - 0.0001)

def test_multiple_scoring_formulas():
    """Test application of multiple scoring formulas."""
    df = create_test_variants_df()
    
    scoring_config = {
        'variable_mapping': create_variable_mapping(),
        'formulas': [
            {'name': 'Score1', 'formula': 'impact_score * rarity'},
            {'name': 'Score2', 'formula': 'cadd + conservation'},
            {'name': 'CombinedScore', 'formula': 'Score1 * 0.5 + Score2 * 0.5'}
        ]
    }
    
    context = PipelineContext(
        current_dataframe=df,
        use_chunked_processing=False
    )
    context.stage_results['scoring_config'] = scoring_config
    
    stage = VariantScoringStage()
    result = stage(context)
    
    # All scores should be added
    assert 'Score1' in result.current_dataframe.columns
    assert 'Score2' in result.current_dataframe.columns
    assert 'CombinedScore' in result.current_dataframe.columns

def test_scoring_with_missing_fields():
    """Test scoring when required fields are missing."""
    df = pd.DataFrame({
        'CHROM': ['chr1'],
        'POS': [100]
        # Missing scoring fields
    })
    
    scoring_config = {
        'variable_mapping': {
            'cadd': {'field': 'CADD_PHRED', 'default': 10}
        },
        'formulas': [
            {'name': 'Score', 'formula': 'cadd * 2'}
        ]
    }
    
    context = PipelineContext(
        current_dataframe=df,
        use_chunked_processing=False
    )
    context.stage_results['scoring_config'] = scoring_config
    
    stage = VariantScoringStage()
    result = stage(context)
    
    # Should use default value
    assert result.current_dataframe['Score'].iloc[0] == 20  # 10 * 2

def test_scoring_formula_error_handling():
    """Test handling of invalid scoring formulas."""
    df = create_test_variants_df()
    
    scoring_config = {
        'variable_mapping': {},
        'formulas': [
            {'name': 'BadScore', 'formula': 'undefined_variable * 2'}
        ]
    }
    
    context = PipelineContext(
        current_dataframe=df,
        use_chunked_processing=False
    )
    context.stage_results['scoring_config'] = scoring_config
    
    stage = VariantScoringStage()
    
    with pytest.raises(NameError):
        stage(context)

def test_scoring_with_inheritance_data():
    """Test scoring that uses inheritance results."""
    df = pd.DataFrame({
        'CADD_PHRED': [20.0],
        'Inheritance_Pattern': ['de_novo'],
        'Inheritance_Confidence': [0.95]
    })
    
    scoring_config = {
        'variable_mapping': {
            'cadd': {'field': 'CADD_PHRED', 'default': 0},
            'is_denovo': {
                'field': 'Inheritance_Pattern',
                'transform': 'lambda x: 1 if x == "de_novo" else 0'
            },
            'inheritance_conf': {'field': 'Inheritance_Confidence', 'default': 0}
        },
        'formulas': [
            {
                'name': 'InheritanceScore',
                'formula': 'cadd * is_denovo * inheritance_conf'
            }
        ]
    }
    
    context = PipelineContext(
        current_dataframe=df,
        use_chunked_processing=False
    )
    context.stage_results['scoring_config'] = scoring_config
    
    stage = VariantScoringStage()
    result = stage(context)
    
    assert abs(result.current_dataframe['InheritanceScore'].iloc[0] - 19.0) < 0.01
```

### 3.5 StatisticsGenerationStage

**Current Code Location**: Uses stats.py and stats_engine.py

**Existing Tests**:
- `test_stats_engine.py`: Statistics engine tests
- `test_stats_integration.py`: Integration tests

**Tests to Adapt**:
```python
# From test_stats_integration.py
def test_statistics_generation():
    # OLD: Called within pipeline
    # NEW: Test stage directly
    
    context = PipelineContext(
        config={'no_stats': False, 'stats_output_file': 'stats.tsv'},
        current_dataframe=df
    )
    stage = StatisticsGenerationStage()
    result = stage(context)
```

**New Tests Required**:
```python
# tests/unit/stages/test_statistics_generation.py

def test_basic_statistics_calculation():
    """Test calculation of basic variant statistics."""
    df = create_test_variants_df(num_variants=100)
    
    context = PipelineContext(
        config={'no_stats': False},
        current_dataframe=df,
        use_chunked_processing=False,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    
    stage = StatisticsGenerationStage()
    result = stage(context)
    
    assert result.statistics is not None
    assert 'total_variants' in result.statistics
    assert 'variants_by_impact' in result.statistics
    assert 'variants_by_gene' in result.statistics
    
    # Verify file written
    stats_file = Path('/tmp/test.statistics.tsv')
    assert stats_file.exists()

def test_statistics_with_custom_config():
    """Test statistics with custom configuration."""
    df = create_test_variants_df()
    
    stats_config = {
        'dataset_stats': [
            {
                'name': 'high_impact_count',
                'expression': 'len(df[df["Impact"] == "HIGH"])'
            }
        ],
        'gene_stats': [
            {
                'name': 'mean_cadd_per_gene',
                'expression': 'df.groupby("Gene_Name")["CADD_PHRED"].mean()'
            }
        ]
    }
    
    context = PipelineContext(
        config={
            'no_stats': False,
            'stats_config': stats_config
        },
        current_dataframe=df,
        use_chunked_processing=False,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    
    stage = StatisticsGenerationStage()
    result = stage(context)
    
    assert 'high_impact_count' in result.statistics
    assert 'mean_cadd_per_gene' in result.statistics

def test_chunked_statistics_aggregation():
    """Test statistics aggregation from chunked processing."""
    chunk_stats = [
        {'total_variants': 100, 'high_impact': 10},
        {'total_variants': 150, 'high_impact': 20},
        {'total_variants': 50, 'high_impact': 5}
    ]
    
    context = PipelineContext(
        config={'no_stats': False},
        use_chunked_processing=True,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    
    # Mock chunk statistics
    context.chunk_statistics = chunk_stats
    
    stage = StatisticsGenerationStage()
    result = stage(context)
    
    # Should aggregate correctly
    assert result.statistics['total_variants'] == 300
    assert result.statistics['high_impact'] == 35

def test_statistics_disabled():
    """Test that statistics generation can be disabled."""
    df = create_test_variants_df()
    
    context = PipelineContext(
        config={'no_stats': True},
        current_dataframe=df
    )
    
    stage = StatisticsGenerationStage()
    result = stage(context)
    
    assert result.statistics == {}
    # No file should be written

def test_statistics_error_handling():
    """Test handling of errors in statistics calculation."""
    df = create_test_variants_df()
    
    stats_config = {
        'dataset_stats': [
            {
                'name': 'bad_stat',
                'expression': '1/0'  # Division by zero
            }
        ]
    }
    
    context = PipelineContext(
        config={
            'no_stats': False,
            'stats_config': stats_config
        },
        current_dataframe=df,
        use_chunked_processing=False,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    
    stage = StatisticsGenerationStage()
    
    # Should handle error gracefully
    result = stage(context)
    assert 'bad_stat' not in result.statistics
    # Other stats should still be calculated
```

---

## Stage 4: Output Stages

### 4.1 VariantIdentifierStage

**Current Code Location**: `pipeline.py` lines 2342-2399

**Existing Tests**:
- No direct tests for VAR_ID generation
- Verified in integration tests

**Tests to Adapt**:
None directly applicable

**New Tests Required**:
```python
# tests/unit/stages/test_variant_identifier.py

def test_variant_id_generation():
    """Test generation of unique variant identifiers."""
    df = pd.DataFrame({
        'CHROM': ['chr1', 'chr1', 'chr2'],
        'POS': [100, 200, 300],
        'REF': ['A', 'G', 'C'],
        'ALT': ['T', 'C', 'G']
    })
    
    context = PipelineContext(
        current_dataframe=df,
        use_chunked_processing=False
    )
    
    stage = VariantIdentifierStage()
    result = stage(context)
    
    assert 'VAR_ID' in result.current_dataframe.columns
    assert result.current_dataframe.columns[0] == 'VAR_ID'
    
    # Check format
    var_ids = result.current_dataframe['VAR_ID']
    assert all(var_ids.str.match(r'^var_\d{4}_[a-f0-9]{4}$'))
    
    # Check uniqueness
    assert len(var_ids.unique()) == len(df)

def test_variant_id_streaming(large_tsv_file):
    """Test streaming VAR_ID addition for large files."""
    context = PipelineContext(
        current_file=large_tsv_file,
        use_chunked_processing=True,
        intermediate_dir=Path('/tmp'),
        base_name='test'
    )
    
    # Monitor memory
    import tracemalloc
    tracemalloc.start()
    
    stage = VariantIdentifierStage()
    result = stage(context)
    
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    # Should use minimal memory
    assert peak / 1024 / 1024 < 100  # Less than 100MB
    
    # Verify IDs added
    with open(result.current_file) as f:
        header = f.readline()
        assert header.startswith('VAR_ID\t')

def test_variant_id_consistency():
    """Test that same variant gets same hash portion."""
    df1 = pd.DataFrame({
        'CHROM': ['chr1'],
        'POS': [100],
        'REF': ['A'],
        'ALT': ['T']
    })
    
    df2 = df1.copy()
    
    context1 = PipelineContext(current_dataframe=df1, use_chunked_processing=False)
    context2 = PipelineContext(current_dataframe=df2, use_chunked_processing=False)
    
    stage = VariantIdentifierStage()
    result1 = stage(context1)
    result2 = stage(context2)
    
    # Hash portion should be identical
    id1 = result1.current_dataframe['VAR_ID'].iloc[0]
    id2 = result2.current_dataframe['VAR_ID'].iloc[0]
    
    assert id1.split('_')[2] == id2.split('_')[2]  # Same hash

def test_missing_required_columns():
    """Test handling when required columns are missing."""
    df = pd.DataFrame({
        'CHROM': ['chr1'],
        'POS': [100]
        # Missing REF and ALT
    })
    
    context = PipelineContext(
        current_dataframe=df,
        use_chunked_processing=False
    )
    
    stage = VariantIdentifierStage()
    result = stage(context)
    
    # Should still generate IDs
    assert 'VAR_ID' in result.current_dataframe.columns
```

### 4.2 FinalFilteringStage

**Current Code Location**: `pipeline.py` lines 2484-2510

**Existing Tests**:
- `test_final_filter_simple.py`: Basic final filtering
- `test_final_filter_integration.py`: Integration tests

**Tests to Adapt**:
```python
# From test_final_filter_simple.py
def test_final_filter():
    # OLD: Through CLI
    # NEW: Test stage directly
    
    df = create_test_df()
    context = PipelineContext(
        config={'final_filter': 'Score > 0.5'},
        current_dataframe=df
    )
    stage = FinalFilteringStage()
    result = stage(context)
```

**New Tests Required**:
```python
# tests/unit/stages/test_final_filtering.py

def test_pandas_query_filtering():
    """Test filtering with pandas query syntax."""
    df = pd.DataFrame({
        'Score': [0.3, 0.6, 0.8, 0.4],
        'Impact': ['HIGH', 'MODERATE', 'HIGH', 'LOW'],
        'AF': [0.001, 0.01, 0.0001, 0.1]
    })
    
    context = PipelineContext(
        config={'final_filter': '(Score > 0.5) & (AF < 0.01)'},
        current_dataframe=df,
        use_chunked_processing=False
    )
    
    stage = FinalFilteringStage()
    result = stage(context)
    
    assert len(result.current_dataframe) == 2
    assert all(result.current_dataframe['Score'] > 0.5)
    assert all(result.current_dataframe['AF'] < 0.01)

def test_complex_filtering_expression():
    """Test complex filtering with multiple conditions."""
    df = create_test_variants_df()
    
    filter_expr = '''
    ((Impact == "HIGH") | (Score > 0.8)) & 
    (Inheritance_Pattern.isin(["de_novo", "autosomal_recessive"])) &
    (AF < 0.01)
    '''
    
    context = PipelineContext(
        config={'final_filter': filter_expr},
        current_dataframe=df,
        use_chunked_processing=False
    )
    
    stage = FinalFilteringStage()
    result = stage(context)
    
    # Verify all conditions met
    filtered = result.current_dataframe
    for _, row in filtered.iterrows():
        assert (row['Impact'] == 'HIGH' or row['Score'] > 0.8)
        assert row['Inheritance_Pattern'] in ['de_novo', 'autosomal_recessive']
        assert row['AF'] < 0.01

def test_no_final_filter():
    """Test stage passes through when no filter configured."""
    df = create_test_variants_df()
    
    context = PipelineContext(
        config={},  # No final_filter
        current_dataframe=df,
        use_chunked_processing=False
    )
    
    stage = FinalFilteringStage()
    result = stage(context)
    
    # Should be unchanged
    pd.testing.assert_frame_equal(result.current_dataframe, df)

def test_filtering_from_chunked_file():
    """Test loading file to memory for filtering."""
    context = PipelineContext(
        config={'final_filter': 'Score > 0.5'},
        current_file=Path('chunked_results.tsv'),
        use_chunked_processing=True
    )
    
    stage = FinalFilteringStage()
    result = stage(context)
    
    # Should load to DataFrame and filter
    assert result.current_dataframe is not None
    assert not result.use_chunked_processing  # Now in memory

def test_empty_filter_result_warning(caplog):
    """Test warning when filter removes all variants."""
    df = pd.DataFrame({
        'Score': [0.1, 0.2, 0.3]
    })
    
    context = PipelineContext(
        config={'final_filter': 'Score > 0.9'},  # Nothing passes
        current_dataframe=df,
        use_chunked_processing=False
    )
    
    stage = FinalFilteringStage()
    result = stage(context)
    
    assert len(result.current_dataframe) == 0
    assert "No variants passed final filter" in caplog.text
```

### 4.3 PseudonymizationStage

**Current Code Location**: `pipeline.py` lines 2511-2562, uses pseudonymizer.py

**Existing Tests**:
- `test_pseudonymizer.py`: Comprehensive pseudonymization tests

**Tests to Adapt**:
```python
# From test_pseudonymizer.py
def test_pseudonymization():
    # OLD: Direct function call
    result_lines, mapper = apply_pseudonymization(lines, samples, config)
    
    # NEW: Through stage
    context = PipelineContext(
        config={'pseudonymize': True, 'pseudonymize_schema': 'sequential'},
        current_dataframe=df
    )
    stage = PseudonymizationStage()
    result = stage(context)
```

**New Tests Required**:
```python
# tests/unit/stages/test_pseudonymization.py

def test_sequential_pseudonymization():
    """Test sequential pseudonymization schema."""
    df = pd.DataFrame({
        'GT': ['Sample1(0/1);Sample2(0/0)', 'Sample3(1/1)'],
        'Other': ['data1', 'data2']
    })
    
    context = PipelineContext(
        config={
            'pseudonymize': True,
            'pseudonymize_schema': 'sequential'
        },
        current_dataframe=df,
        output_dir=Path('/tmp/output')
    )
    
    stage = PseudonymizationStage()
    result = stage(context)
    
    # Check pseudonymization applied
    assert 'SAMPLE_001' in result.current_dataframe['GT'].iloc[0]
    assert 'SAMPLE_002' in result.current_dataframe['GT'].iloc[0]
    assert 'SAMPLE_003' in result.current_dataframe['GT'].iloc[1]
    
    # Check mapping saved
    assert context.pseudonymizer is not None

def test_categorical_pseudonymization():
    """Test categorical pseudonymization with metadata."""
    df = pd.DataFrame({
        'GT': ['Case1(0/1);Control1(0/0)', 'Case2(1/1)']
    })
    
    metadata = {
        'Case1': {'status': 'case'},
        'Case2': {'status': 'case'},
        'Control1': {'status': 'control'}
    }
    
    context = PipelineContext(
        config={
            'pseudonymize': True,
            'pseudonymize_schema': 'categorical',
            'pseudonymize_metadata': metadata
        },
        current_dataframe=df,
        output_dir=Path('/tmp/output')
    )
    
    stage = PseudonymizationStage()
    result = stage(context)
    
    # Check categorical naming
    assert 'CASE_001' in result.current_dataframe['GT'].iloc[0]
    assert 'CONTROL_001' in result.current_dataframe['GT'].iloc[0]
    assert 'CASE_002' in result.current_dataframe['GT'].iloc[1]

def test_sample_extraction_from_gt():
    """Test extraction of sample names from GT column."""
    df = pd.DataFrame({
        'GT': [
            'Sample1(0/1:20:99);Sample2(0/0)',
            'Sample3(1|1)',
            'Sample4(./.)',
            ''  # Empty
        ]
    })
    
    context = PipelineContext(
        config={'pseudonymize': True},
        current_dataframe=df,
        output_dir=Path('/tmp/output')
    )
    
    stage = PseudonymizationStage()
    samples = stage._extract_samples(context)
    
    assert samples == {'Sample1', 'Sample2', 'Sample3', 'Sample4'}

def test_pseudonymization_file_outputs(tmp_path):
    """Test generation of mapping and PED files."""
    df = create_test_df_with_samples()
    pedigree_data = create_test_pedigree()
    
    context = PipelineContext(
        config={
            'pseudonymize': True,
            'pseudonymize_table': 'mapping.tsv',
            'pseudonymize_ped': True,
            'ped_file': 'original.ped'
        },
        current_dataframe=df,
        output_dir=tmp_path / 'output'
    )
    context.stage_results['pedigree_data'] = pedigree_data
    
    stage = PseudonymizationStage()
    result = stage(context)
    
    # Check files created in parent directory
    parent_dir = tmp_path
    assert (parent_dir / 'mapping.tsv').exists()
    assert (parent_dir / 'mapping.ped').exists()

def test_no_pseudonymization():
    """Test stage skips when not configured."""
    df = create_test_df_with_samples()
    
    context = PipelineContext(
        config={},  # No pseudonymize flag
        current_dataframe=df
    )
    
    stage = PseudonymizationStage()
    result = stage(context)
    
    # Should be unchanged
    pd.testing.assert_frame_equal(result.current_dataframe, df)
```

### 4.4 TSVOutputStage

**Current Code Location**: `pipeline.py` lines 2563-2592

**Existing Tests**:
- Integration tests verify output
- No specific unit tests

**Tests to Adapt**:
None directly applicable

**New Tests Required**:
```python
# tests/unit/stages/test_tsv_output.py

def test_file_output(tmp_path):
    """Test writing TSV to file."""
    df = create_test_variants_df()
    
    context = PipelineContext(
        args=Namespace(output_file='results.tsv'),
        current_dataframe=df,
        output_dir=tmp_path
    )
    
    stage = TSVOutputStage()
    result = stage(context)
    
    assert result.final_output_path == tmp_path / 'results.tsv'
    assert result.final_output_path.exists()
    
    # Verify content
    written_df = pd.read_csv(result.final_output_path, sep='\t')
    pd.testing.assert_frame_equal(written_df, df)

def test_stdout_output(capsys):
    """Test writing TSV to stdout."""
    df = pd.DataFrame({
        'CHROM': ['chr1'],
        'POS': [100]
    })
    
    context = PipelineContext(
        args=Namespace(output_file='stdout'),
        current_dataframe=df
    )
    
    stage = TSVOutputStage()
    result = stage(context)
    
    assert result.final_output_path is None
    
    captured = capsys.readouterr()
    assert 'CHROM\tPOS' in captured.out
    assert 'chr1\t100' in captured.out

def test_streaming_file_copy(large_tsv_file, tmp_path):
    """Test copying from current_file when no DataFrame."""
    context = PipelineContext(
        args=Namespace(output_file='results.tsv'),
        current_file=large_tsv_file,
        current_dataframe=None,
        output_dir=tmp_path
    )
    
    stage = TSVOutputStage()
    result = stage(context)
    
    assert result.final_output_path.exists()
    
    # Should be identical copy
    assert filecmp.cmp(large_tsv_file, result.final_output_path)

def test_compressed_output(tmp_path):
    """Test handling of compressed output."""
    df = create_test_variants_df()
    
    context = PipelineContext(
        args=Namespace(output_file='results.tsv.gz'),
        current_dataframe=df,
        output_dir=tmp_path
    )
    
    stage = TSVOutputStage()
    result = stage(context)
    
    assert result.final_output_path.suffix == '.gz'
    
    # Verify can read compressed
    written_df = pd.read_csv(result.final_output_path, sep='\t', compression='gzip')
    assert len(written_df) == len(df)
```

### 4.5 Report Generation Stages

**Current Code Location**: Various report generators

**Existing Tests**:
- `test_igv_report_generation.py`: IGV tests
- Integration tests for reports

**Tests to Adapt**:
Report generation tests need significant adaptation

**New Tests Required**:
```python
# tests/unit/stages/test_report_stages.py

def test_excel_report_generation(tmp_path):
    """Test Excel report with multiple sheets."""
    # Create test files
    tsv_file = tmp_path / 'results.tsv'
    create_test_tsv(tsv_file)
    
    metadata_file = tmp_path / 'metadata.tsv'
    create_metadata_file(metadata_file)
    
    context = PipelineContext(
        args=Namespace(xlsx=True),
        final_output_path=tsv_file,
        metadata_path=metadata_file,
        output_dir=tmp_path,
        base_name='test'
    )
    
    stage = ExcelReportStage()
    result = stage(context)
    
    assert result.excel_path.exists()
    assert result.excel_path.suffix == '.xlsx'
    
    # Verify sheets
    with pd.ExcelFile(result.excel_path) as xls:
        assert 'Variants' in xls.sheet_names
        assert 'Metadata' in xls.sheet_names

def test_html_report_generation(tmp_path):
    """Test HTML report generation."""
    tsv_file = tmp_path / 'results.tsv'
    create_test_tsv(tsv_file)
    
    context = PipelineContext(
        args=Namespace(html_report=True),
        final_output_path=tsv_file,
        output_dir=tmp_path,
        config={}
    )
    
    stage = HTMLReportStage()
    result = stage(context)
    
    assert result.html_report_path.exists()
    assert result.html_report_path.name == 'index.html'
    
    # Check supporting files
    report_dir = result.html_report_path.parent
    assert (report_dir / 'variants.json').exists()
    assert (report_dir / 'summary.json').exists()

def test_igv_report_validation():
    """Test IGV configuration validation."""
    context = PipelineContext(
        config={'igv_enabled': True},
        # Missing required config
        final_output_path=Path('results.tsv')
    )
    
    stage = IGVReportStage()
    
    with pytest.raises(ValueError, match="requires --bam-mapping-file"):
        stage(context)

def test_metadata_generation(tmp_path):
    """Test metadata file generation."""
    context = PipelineContext(
        config={'version': '1.0.0', 'reference': 'hg38'},
        args=Namespace(),
        start_time=datetime.now(),
        output_dir=tmp_path,
        base_name='test',
        final_output_path=tmp_path / 'results.tsv'
    )
    
    stage = MetadataGenerationStage()
    result = stage(context)
    
    assert result.metadata_path.exists()
    
    # Verify content
    metadata = pd.read_csv(result.metadata_path, sep='\t')
    assert 'Tool' in metadata['Parameter'].values
    assert 'Version' in metadata['Parameter'].values
    assert 'Run_start_time' in metadata['Parameter'].values

def test_parallel_report_generation():
    """Test parallel generation of multiple reports."""
    # This would test ParallelReportGenerationStage
    # but requires more complex setup
    pass
```

---

## Integration Test Requirements

### Stage Interaction Tests

```python
# tests/integration/stages/test_stage_interactions.py

def test_configuration_to_processing_flow():
    """Test data flow from configuration stages to processing."""
    # Load all configurations
    config_stages = [
        ConfigurationLoadingStage(),
        PhenotypeLoadingStage(),
        ScoringConfigLoadingStage(),
        PedigreeLoadingStage()
    ]
    
    # Then processing stages that depend on config
    processing_stages = [
        GeneBedCreationStage(),
        VariantExtractionStage()
    ]
    
    # Run full sequence and verify data flows correctly

def test_processing_to_analysis_flow():
    """Test data flow from processing to analysis stages."""
    # Test that extracted TSV properly flows to analysis
    
def test_analysis_to_output_flow():
    """Test data flow from analysis to output stages."""
    # Test scored and annotated data reaches output correctly

def test_parallel_stage_execution():
    """Test that parallel stages execute correctly."""
    # Test configuration stages run in parallel
    # Test report stages run in parallel
```

### Full Pipeline Tests

```python
# tests/integration/test_refactored_pipeline.py

def test_complete_refactored_pipeline():
    """Test entire pipeline with all stages."""
    # This is the ultimate test - run all stages
    # and compare with original pipeline output
    
def test_checkpoint_resume_with_stages():
    """Test checkpoint/resume works with stage architecture."""
    # Run partial pipeline
    # Kill and resume
    # Verify stages properly skip/resume
```

## Summary

This comprehensive testing analysis identifies:

1. **Existing tests to adapt**: ~40 test functions
2. **New tests required**: ~150 test functions
3. **Key testing areas**:
   - Stage isolation and dependencies
   - Parallel execution correctness
   - Memory efficiency for large files
   - Error handling and edge cases
   - Integration between stages

Each stage needs 5-10 dedicated unit tests plus integration tests to verify interactions. The refactoring will result in much better test coverage and maintainability.