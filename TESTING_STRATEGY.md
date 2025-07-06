# Testing Strategy for Refactored VariantCentrifuge

## Executive Summary

This document outlines a comprehensive testing strategy for the refactored VariantCentrifuge pipeline. The strategy emphasizes testability, maintainability, and confidence in the refactored codebase through multiple layers of testing.

## Testing Philosophy

### Core Principles

1. **Test Pyramid Approach**
   - Many unit tests (fast, isolated)
   - Moderate integration tests (component interactions)
   - Few end-to-end tests (full pipeline validation)

2. **Test-Driven Refactoring**
   - Write tests for existing functionality first
   - Refactor with confidence
   - Ensure no regression

3. **Continuous Testing**
   - Tests run on every commit
   - Fast feedback loop
   - Automated test execution

## Testing Layers

### 1. Unit Testing

#### Purpose
Test individual components in isolation with mocked dependencies.

#### Coverage Goals
- **Target**: 85% code coverage per module
- **Critical paths**: 100% coverage
- **Error handling**: 100% coverage

#### Testing Framework
```python
# Using pytest with plugins
pytest==7.4.0
pytest-cov==4.1.0
pytest-mock==3.11.1
pytest-asyncio==0.21.0
pytest-xdist==3.3.1  # Parallel test execution
pytest-benchmark==4.0.0  # Performance testing
```

#### Unit Test Structure

```python
# tests/unit/processors/test_gene_processor.py

import pytest
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path

from variantcentrifuge.processors.gene_processor import GeneProcessor
from variantcentrifuge.exceptions import ProcessingError

class TestGeneProcessor:
    """Unit tests for GeneProcessor."""
    
    @pytest.fixture
    def processor(self):
        """Create GeneProcessor instance for testing."""
        return GeneProcessor()
    
    @pytest.fixture
    def mock_config(self):
        """Create mock configuration."""
        return {
            'reference': 'hg38',
            'snpeff_path': '/usr/local/bin/snpEff',
            'interval_expansion': 0
        }
    
    def test_normalize_genes_single_gene(self, processor):
        """Test normalizing a single gene name."""
        result = processor.normalize_genes("brca1", None)
        assert result == "BRCA1"
    
    def test_normalize_genes_from_file(self, processor, tmp_path):
        """Test normalizing genes from a file."""
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("BRCA1\nTP53\nEGFR\n")
        
        result = processor.normalize_genes(None, gene_file)
        assert result == "BRCA1,TP53,EGFR"
    
    @patch('subprocess.run')
    def test_create_bed_file_success(self, mock_run, processor, tmp_path):
        """Test successful BED file creation."""
        mock_run.return_value = Mock(returncode=0, stdout="Success")
        
        result = processor.create_bed_file(
            "BRCA1",
            "hg38",
            tmp_path,
            interval_expansion=100
        )
        
        assert result.exists()
        mock_run.assert_called_once()
    
    @patch('subprocess.run')
    def test_create_bed_file_failure(self, mock_run, processor, tmp_path):
        """Test BED file creation failure handling."""
        mock_run.return_value = Mock(
            returncode=1,
            stderr="Gene not found in database"
        )
        
        with pytest.raises(ProcessingError) as exc_info:
            processor.create_bed_file("INVALID_GENE", "hg38", tmp_path)
        
        assert "Gene not found" in str(exc_info.value)
    
    @pytest.mark.parametrize("genes,expected", [
        ("BRCA1", "BRCA1"),
        ("brca1,tp53", "BRCA1,TP53"),
        ("  EGFR  ,  KRAS  ", "EGFR,KRAS"),
        ("", ""),
    ])
    def test_normalize_genes_various_inputs(self, processor, genes, expected):
        """Test gene normalization with various inputs."""
        result = processor.normalize_genes(genes, None)
        assert result == expected
```

#### Mocking Strategy

```python
# tests/unit/mocks.py

from unittest.mock import Mock, MagicMock
from typing import Dict, Any, List
import pandas as pd

class MockFactory:
    """Factory for creating consistent mocks."""
    
    @staticmethod
    def create_mock_vcf_reader():
        """Create mock VCF reader."""
        mock = MagicMock()
        mock.samples = ['Sample1', 'Sample2', 'Sample3']
        mock.__iter__.return_value = [
            Mock(CHROM='chr1', POS=12345, REF='A', ALT=['G']),
            Mock(CHROM='chr2', POS=67890, REF='C', ALT=['T'])
        ]
        return mock
    
    @staticmethod
    def create_mock_dataframe(num_rows: int = 10) -> pd.DataFrame:
        """Create mock variant dataframe."""
        data = {
            'CHROM': [f'chr{i%23 or 23}' for i in range(num_rows)],
            'POS': list(range(1000, 1000 + num_rows)),
            'REF': ['A'] * num_rows,
            'ALT': ['G'] * num_rows,
            'Gene': [f'GENE{i%5}' for i in range(num_rows)],
            'Impact': ['HIGH', 'MODERATE', 'LOW'][i%3] for i in range(num_rows)],
            'GT': ['Sample1(0/1)'] * num_rows
        }
        return pd.DataFrame(data)
    
    @staticmethod
    def create_mock_config(overrides: Dict[str, Any] = None) -> Dict[str, Any]:
        """Create mock configuration."""
        config = {
            'reference': 'hg38',
            'fields_to_extract': ['CHROM', 'POS', 'REF', 'ALT'],
            'filters': 'FILTER="PASS"',
            'no_stats': False,
            'perform_gene_burden': False
        }
        if overrides:
            config.update(overrides)
        return config
```

### 2. Integration Testing

#### Purpose
Test interactions between components and with external tools.

#### Integration Test Categories

1. **Component Integration**
   ```python
   # tests/integration/test_processor_integration.py
   
   class TestProcessorIntegration:
       """Test integration between processors."""
       
       def test_gene_to_variant_flow(self, tmp_path):
           """Test flow from gene processor to variant processor."""
           # Setup
           gene_processor = GeneProcessor()
           variant_processor = VariantProcessor()
           
           # Create test VCF
           test_vcf = create_test_vcf(tmp_path / "test.vcf")
           
           # Process genes
           bed_file = gene_processor.create_bed_file(
               "BRCA1,TP53", "hg38", tmp_path
           )
           
           # Extract variants
           variants = variant_processor.extract_variants_by_region(
               test_vcf, bed_file, tmp_path / "output.vcf"
           )
           
           # Verify
           assert variants.exists()
           assert count_variants(variants) > 0
   ```

2. **External Tool Integration**
   ```python
   # tests/integration/test_external_tools.py
   
   @pytest.mark.external_tools
   class TestExternalToolIntegration:
       """Test integration with external bioinformatics tools."""
       
       def test_bcftools_integration(self, test_vcf):
           """Test bcftools command execution."""
           from variantcentrifuge.utils import run_command
           
           result = run_command([
               "bcftools", "view", "-H", str(test_vcf)
           ])
           
           assert result.returncode == 0
           assert len(result.stdout.splitlines()) > 0
       
       @pytest.mark.skipif(not shutil.which("snpEff"), 
                          reason="snpEff not available")
       def test_snpeff_integration(self, test_vcf):
           """Test snpEff annotation."""
           # Test snpEff functionality
           pass
   ```

3. **Database Integration**
   ```python
   # tests/integration/test_annotation_integration.py
   
   class TestAnnotationIntegration:
       """Test annotation system integration."""
       
       def test_unified_annotation_flow(self, sample_variants_df):
           """Test complete annotation workflow."""
           # Setup annotation sources
           bed_file = create_test_bed_file()
           gene_list = ['BRCA1', 'TP53', 'EGFR']
           json_data = {'BRCA1': {'panel': 'cancer', 'inheritance': 'AD'}}
           
           # Create annotation config
           custom_features = {
               'regions_by_chrom': load_bed_to_intervals(bed_file),
               'gene_lists': {'cancer_genes': gene_list},
               'json_gene_data': json_data
           }
           
           # Apply annotations
           annotated_df = annotate_dataframe_with_features(
               sample_variants_df, custom_features
           )
           
           # Verify annotations
           assert 'Custom_Annotation' in annotated_df.columns
           assert annotated_df['Custom_Annotation'].notna().any()
   ```

### 3. End-to-End Testing

#### Purpose
Validate complete pipeline functionality with real-world scenarios.

#### E2E Test Scenarios

```python
# tests/e2e/test_full_pipeline.py

@pytest.mark.e2e
class TestFullPipeline:
    """End-to-end pipeline tests."""
    
    def test_single_gene_analysis(self, test_data_dir, tmp_path):
        """Test complete pipeline for single gene."""
        args = create_test_args(
            vcf_file=test_data_dir / "sample.vcf",
            gene_name="BRCA1",
            output_dir=tmp_path,
            config={
                'preset': 'rare,coding',
                'calculate_inheritance': True
            }
        )
        
        # Run pipeline
        run_pipeline(args, load_config(), datetime.now())
        
        # Verify outputs
        assert (tmp_path / "sample_BRCA1.final.tsv").exists()
        assert (tmp_path / "sample_BRCA1.metadata.tsv").exists()
        
        # Verify content
        df = pd.read_csv(tmp_path / "sample_BRCA1.final.tsv", sep='\t')
        assert len(df) > 0
        assert 'VAR_ID' in df.columns
        assert 'Inheritance_Pattern' in df.columns
    
    def test_multi_sample_cohort_analysis(self, cohort_vcf, tmp_path):
        """Test pipeline with multi-sample cohort."""
        args = create_test_args(
            vcf_file=cohort_vcf,
            gene_file="cancer_genes.txt",
            output_dir=tmp_path,
            config={
                'perform_gene_burden': True,
                'case_samples_file': 'cases.txt',
                'control_samples_file': 'controls.txt'
            }
        )
        
        # Run pipeline
        run_pipeline(args, load_config(), datetime.now())
        
        # Verify gene burden output
        burden_file = tmp_path / "cohort_cancer_genes.gene_burden.tsv"
        assert burden_file.exists()
        
        burden_df = pd.read_csv(burden_file, sep='\t')
        assert 'p_value' in burden_df.columns
        assert 'odds_ratio' in burden_df.columns
```

### 4. Performance Testing

#### Purpose
Ensure refactored code meets performance requirements.

#### Performance Benchmarks

```python
# tests/performance/test_pipeline_performance.py

@pytest.mark.benchmark
class TestPipelinePerformance:
    """Performance benchmarks for pipeline components."""
    
    def test_variant_processing_speed(self, benchmark, large_vcf):
        """Benchmark variant processing speed."""
        processor = VariantProcessor()
        
        def process_variants():
            return processor.extract_fields(
                large_vcf, 
                ['CHROM', 'POS', 'REF', 'ALT'],
                "output.tsv"
            )
        
        result = benchmark(process_variants)
        
        # Assert performance requirements
        assert benchmark.stats['mean'] < 10.0  # Less than 10 seconds
    
    def test_memory_usage(self, large_tsv_file):
        """Test memory usage stays within limits."""
        import tracemalloc
        
        tracemalloc.start()
        
        # Process large file
        processor = ChunkProcessor()
        for chunk in processor.read_tsv_in_chunks(large_tsv_file, chunksize=10000):
            # Process chunk
            pass
        
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        # Assert memory usage under 2GB
        assert peak / 1024 / 1024 / 1024 < 2.0
```

### 5. Regression Testing

#### Purpose
Ensure refactored code produces identical results to original.

#### Regression Test Suite

```python
# tests/regression/test_output_regression.py

class TestOutputRegression:
    """Regression tests comparing old and new implementations."""
    
    def test_variant_extraction_regression(self, test_vcf, reference_output):
        """Compare variant extraction with reference output."""
        # Run new implementation
        new_output = run_new_pipeline(test_vcf)
        
        # Compare with reference
        assert compare_tsv_files(
            new_output, 
            reference_output,
            ignore_columns=['VAR_ID', 'timestamp']
        )
    
    def test_scoring_regression(self, scored_reference):
        """Ensure scoring produces same results."""
        # Test scoring consistency
        pass
```

## Test Data Management

### Test Data Categories

1. **Minimal Test Data**
   ```
   tests/data/minimal/
   ├── single_variant.vcf
   ├── single_gene.bed
   └── single_sample.ped
   ```

2. **Functional Test Data**
   ```
   tests/data/functional/
   ├── multi_gene.vcf
   ├── cohort_10_samples.vcf
   ├── complex_genotypes.vcf
   └── inheritance_test_cases/
   ```

3. **Performance Test Data**
   ```
   tests/data/performance/
   ├── large_cohort_100_samples.vcf.gz
   ├── whole_exome_variants.vcf.gz
   └── stress_test_1M_variants.vcf.gz
   ```

### Test Data Generation

```python
# tests/fixtures/data_generators.py

def generate_test_vcf(num_variants: int, num_samples: int, 
                     output_path: Path) -> Path:
    """Generate test VCF with specified characteristics."""
    # VCF generation logic
    pass

def generate_test_pedigree(family_structure: str, 
                          output_path: Path) -> Path:
    """Generate test pedigree file."""
    # Pedigree generation logic
    pass
```

## Testing Infrastructure

### Continuous Integration

```yaml
# .github/workflows/tests.yml

name: Tests

on: [push, pull_request]

jobs:
  unit-tests:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, 3.10, 3.11]
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}
    
    - name: Install dependencies
      run: |
        pip install -e ".[test]"
    
    - name: Run unit tests
      run: |
        pytest tests/unit -v --cov=variantcentrifuge --cov-report=xml
    
    - name: Upload coverage
      uses: codecov/codecov-action@v3
  
  integration-tests:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Install bioinformatics tools
      run: |
        conda install -c bioconda bcftools snpeff snpsift bedtools
    
    - name: Run integration tests
      run: |
        pytest tests/integration -v -m "not slow"
  
  e2e-tests:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Full environment setup
      run: |
        # Complete setup script
    
    - name: Run E2E tests
      run: |
        pytest tests/e2e -v
```

### Test Execution Strategy

```bash
# Run all tests
pytest

# Run only unit tests
pytest tests/unit

# Run tests in parallel
pytest -n auto

# Run with coverage
pytest --cov=variantcentrifuge --cov-report=html

# Run specific test categories
pytest -m "not slow and not external_tools"

# Run performance tests
pytest tests/performance --benchmark-only

# Run regression tests
pytest tests/regression --regression-baseline=v1.0
```

## Test Quality Metrics

### Coverage Requirements

1. **Overall Coverage**: > 80%
2. **Core Modules**: > 90%
3. **Error Handling**: 100%
4. **New Code**: > 95%

### Test Execution Time

1. **Unit Tests**: < 30 seconds
2. **Integration Tests**: < 5 minutes
3. **E2E Tests**: < 15 minutes
4. **Full Test Suite**: < 30 minutes

### Test Maintainability

1. **DRY Principle**: Shared fixtures and utilities
2. **Clear Names**: Descriptive test names
3. **Documentation**: Docstrings for complex tests
4. **Independence**: Tests don't depend on execution order

## Testing Best Practices

### 1. Test Organization

```python
# Good test organization
tests/
├── conftest.py          # Shared fixtures
├── unit/
│   ├── conftest.py      # Unit-specific fixtures
│   ├── processors/
│   ├── analyzers/
│   └── reporters/
├── integration/
│   ├── conftest.py      # Integration fixtures
│   └── workflows/
├── e2e/
│   └── scenarios/
└── fixtures/
    ├── data/
    └── mocks/
```

### 2. Fixture Design

```python
# conftest.py - Shared fixtures

@pytest.fixture(scope="session")
def test_data_dir():
    """Path to test data directory."""
    return Path(__file__).parent / "data"

@pytest.fixture
def temp_workspace(tmp_path):
    """Create temporary workspace with structure."""
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    (workspace / "input").mkdir()
    (workspace / "output").mkdir()
    return workspace

@pytest.fixture
def mock_external_tools(monkeypatch):
    """Mock external tool calls."""
    def mock_run_command(cmd, *args, **kwargs):
        if cmd[0] == "bcftools":
            return Mock(returncode=0, stdout="Mocked output")
        return Mock(returncode=1, stderr="Tool not found")
    
    monkeypatch.setattr(
        "variantcentrifuge.utils.run_command",
        mock_run_command
    )
```

### 3. Assertion Patterns

```python
# Clear, specific assertions

# Bad
assert result

# Good
assert result is not None, "Expected result but got None"
assert len(result) == 5, f"Expected 5 items, got {len(result)}"

# Using pytest features
with pytest.raises(ValueError, match="Invalid gene name"):
    processor.normalize_genes("123invalid")

# Approximate comparisons
assert result == pytest.approx(0.95, rel=0.01)

# Custom assertions
def assert_valid_vcf(vcf_path):
    """Assert that file is a valid VCF."""
    assert vcf_path.exists(), f"VCF file not found: {vcf_path}"
    with open(vcf_path) as f:
        first_line = f.readline()
        assert first_line.startswith("##fileformat=VCF")
```

## Test Documentation

### Test Case Documentation

```python
class TestInheritanceAnalyzer:
    """Test inheritance pattern analysis.
    
    These tests validate the inheritance analyzer's ability to:
    1. Detect de novo variants
    2. Identify compound heterozygous pairs
    3. Calculate segregation scores
    4. Handle edge cases in pedigree data
    """
    
    def test_de_novo_detection(self, trio_data):
        """Test de novo variant detection in parent-child trio.
        
        Given:
            - Child has heterozygous variant
            - Both parents are homozygous reference
            - Variant has high quality scores
        
        Expect:
            - Inheritance pattern marked as 'de_novo'
            - High confidence score (>0.9)
        """
        # Test implementation
```

### Test Result Reporting

```python
# pytest.ini configuration

[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_classes = ["Test*"]
python_functions = ["test_*"]
markers = [
    "slow: marks tests as slow",
    "external_tools: requires external bioinformatics tools",
    "benchmark: performance benchmark tests",
    "e2e: end-to-end tests",
    "regression: regression tests"
]
addopts = [
    "--strict-markers",
    "--verbose",
    "--tb=short",
    "--junit-xml=test-results.xml"
]
```

## Conclusion

This comprehensive testing strategy ensures:

1. **Confidence**: Refactored code works correctly
2. **Quality**: High code quality standards
3. **Performance**: Meets performance requirements
4. **Maintainability**: Tests are easy to understand and maintain
5. **Documentation**: Tests serve as living documentation

By following this strategy, the refactoring process will be safer, more reliable, and result in a more robust codebase.