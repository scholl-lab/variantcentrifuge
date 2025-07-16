# Inheritance Analysis Test Suite

This directory contains comprehensive tests for the enhanced inheritance analysis system in VariantCentrifuge.

## Test Organization

### Unit Tests

- **`test_enhanced_deducer.py`** - Tests for pattern deduction logic (Pass 1)
  - Single sample patterns
  - De novo detection
  - Dominant/recessive patterns
  - X-linked patterns
  - Mitochondrial patterns

- **`test_enhanced_comp_het.py`** - Tests for compound heterozygous analysis (Pass 2)
  - Trans/cis configuration detection
  - Single sample compound het detection
  - Missing data handling

- **`test_segregation_checker.py`** - Tests for segregation analysis
  - Pattern-specific segregation checks
  - Confidence score calculation
  - Incomplete penetrance handling

- **`test_enhanced_prioritizer.py`** - Tests for pattern prioritization (Pass 3)
  - Priority scoring with segregation
  - Variant/sample info adjustments
  - Confidence calculation
  - Pattern categorization

### Genotype Format Tests

- **`test_colon_separated_genotypes.py`** - Tests for colon-separated GT format
  - Handling GEN[*].GT extraction format
  - Performance with many samples
  - Different separators
  - Missing data handling

- **`test_gen_star_gt_format.py`** - Documentation tests for GEN[*].GT format
  - Regression tests for compound het detection
  - Format expectations
  - Pipeline behavior documentation

### Integration Tests

- **`test_inheritance_integration.py`** - Full workflow tests
  - Complete 3-pass analysis
  - Single sample vs family analysis
  - Summary generation
  - Pattern filtering

- **`test_variant_linker_compatibility.py`** - Compatibility tests
  - Ensures behavior matches variant-linker tool
  - Complex pattern scenarios
  - Edge cases

## Running Tests

Run all inheritance tests:
```bash
pytest tests/test_inheritance/
```

Run specific test file:
```bash
pytest tests/test_inheritance/test_enhanced_deducer.py
```

Run with verbose output:
```bash
pytest -v tests/test_inheritance/
```

Run with coverage:
```bash
pytest --cov=variantcentrifuge.inheritance tests/test_inheritance/
```

## Test Fixtures

The `conftest.py` file provides common test fixtures:
- `sample_pedigree_data` - Basic trio pedigree
- `single_sample_pedigree` - Single affected sample
- `multi_generation_pedigree` - Extended family
- `sample_variant_df` - Basic variant DataFrame
- `compound_het_df` - Variants for compound het testing
- `x_linked_variants` - X-linked test data
- `mitochondrial_variants` - Mitochondrial test data

## Test Coverage

The test suite covers:
1. **Pattern Detection** - All Mendelian inheritance patterns
2. **Compound Het Analysis** - Including phase determination
3. **Segregation Analysis** - Pattern validation in families
4. **Prioritization** - Choosing best pattern with confidence
5. **Edge Cases** - Missing data, single samples, errors
6. **Integration** - Full pipeline functionality
7. **Compatibility** - Matching variant-linker behavior
