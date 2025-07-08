# Regression Test Issues Summary

## Overview

The regression tests are failing due to several configuration and argument parsing issues when comparing the old and new pipeline implementations. Out of 20 tests, all 20 failed with various errors.

## Error Categories

### 1. Unrecognized Arguments (6 failures)
**Issue**: The test construction is incorrectly passing file paths as extra positional arguments instead of as option values.

**Affected Tests**:
- `gene_with_inheritance`: `variantcentrifuge: error: unrecognized arguments: test_family.ped`
- `gene_with_annotations`: `variantcentrifuge: error: unrecognized arguments: test_regions.bed`

**Root Cause**: In `test_regression_suite.py`, the `extra_args` list contains:
```python
# Incorrect
"--ped", "test_family.ped",  # This becomes two separate args
"--annotate-bed", "test_regions.bed"  # This becomes two separate args

# Should be single arguments:
"--ped=test_family.ped"
# OR passed as part of cmd.extend() properly
```

### 2. TypeError: sequence item 3: expected str instance, NoneType found (7 failures)
**Issue**: The test command construction is including `None` values in the command list.

**Affected Tests**:
- `test_basic_gene_extraction`
- `test_with_filtering`
- `test_parallel_extraction`
- `test_with_presets[rare]`
- `test_with_presets[coding]`
- `test_with_presets[pathogenic]`

**Root Cause**: The `logger.info()` call is trying to join a list containing `None` values.

### 3. Incorrect Field Extraction Flag
**Issue**: Tests are using `--extract` but the CLI expects `--fields` or `-e`.

**Evidence**: In `cli.py` line 157:
```python
format_group.add_argument("-e", "--fields", help="Fields to extract with SnpSift extractFields")
```

But regression tests use `--extract` which doesn't exist.

### 4. Filter Expression Format Issues
**Issue**: New pipeline is failing on rare preset filter with gnomAD fields.

**Error**: 
```
Command failed: SnpSift filter ((((dbNSFP_gnomAD_exomes_AF[0] < 0.0001) | (na dbNSFP_gnomAD_exomes_AC[0])) & ((dbNSFP_gnomAD_genomes_AF[0] < 0.0001) | (na dbNSFP_gnomAD_genomes_AC[0]))))
```

**Possible Cause**: Missing fields in test VCF or incorrect filter syntax for test data.

### 5. Missing Test Files
**Issue**: Several tests reference files that don't exist in the expected locations.

**Missing Files**:
- `test_genes.txt`
- `test_family.ped`
- `test_regions.bed`
- `scoring/inheritance_score`

### 6. --replace-genotypes Flag Removed
**Issue**: Test uses `--replace-genotypes` which appears to have been removed or renamed.

**Error**: `pytest: error: unrecognized arguments: --replace-genotypes`

## Detailed Fix Plan

### Phase 1: Fix Test Infrastructure (Priority: High)

#### 1.1 Fix Argument Construction in test_regression_suite.py
```python
# In RegressionTestConfig definitions, change:
extra_args=["--ped", "test_family.ped"]
# To:
extra_args=["--ped=tests/fixtures/test_family.ped"]

# OR fix in the run() method to handle multi-arg options properly
```

#### 1.2 Create Missing Test Files
Create the following test files in `tests/fixtures/`:
- `test_genes.txt` - List of test genes (BRCA1, TP53, etc.)
- `test_family.ped` - Simple trio family structure
- `test_regions.bed` - BED file with test regions

#### 1.3 Fix Field Extraction Arguments
Replace all instances of `--extract` with `--fields` in:
- `tests/regression/test_pipeline_regression.py`
- `tests/regression/test_regression_suite.py`
- `tests/regression/test_baseline_comparison.py`

### Phase 2: Fix Pipeline Compatibility Issues (Priority: High)

#### 2.1 Update Filter Handling
Check if the test VCF has the required annotation fields for filters:
- `dbNSFP_gnomAD_exomes_AF`
- `dbNSFP_gnomAD_genomes_AF`
- `dbNSFP_gnomAD_exomes_AC`
- `dbNSFP_gnomAD_genomes_AC`

If not, either:
- Update test VCF to include these fields
- Modify test filters to use available fields
- Add field existence checks in the pipeline

#### 2.2 Handle Removed/Renamed Arguments
Search for and update any uses of `--replace-genotypes`:
- Check if it was renamed to `--no-replacement` (inverse logic)
- Update tests accordingly

### Phase 3: Fix Command Construction Logic (Priority: Medium)

#### 3.1 Fix None Values in Commands
In `test_pipeline_regression.py`, ensure all command parts are strings:
```python
# Add validation before logging
cmd = [str(x) for x in cmd if x is not None]
logger.info("Running: %s", " ".join(cmd))
```

#### 3.2 Update Path Resolution
Ensure all test file paths are resolved correctly:
```python
# In test files, use proper path resolution
test_file = Path(__file__).parent.parent / "fixtures" / "test_family.ped"
```

### Phase 4: Validate Test Data (Priority: Medium)

#### 4.1 Create Comprehensive Test VCF
Ensure test VCF (`test_regression_variants.GRCh37.annotated.vcf.gz`) contains:
- All required annotation fields
- Multiple genes for testing
- Various variant types
- Sample genotypes for inheritance testing

#### 4.2 Create Scoring Configuration
Create a minimal scoring configuration at `scoring/inheritance_score/` or update path to existing config.

### Phase 5: Update Regression Test Suite (Priority: Low)

#### 5.1 Add Debug Output
Add more detailed logging to understand failures:
```python
# Log actual command being run
logger.debug(f"Command parts: {cmd}")
# Log working directory
logger.debug(f"Working directory: {os.getcwd()}")
```

#### 5.2 Add Test Fixtures Setup
Create a `conftest.py` or setup method to ensure all test files exist before running tests.

## Implementation Order

1. **Immediate Fixes** (Fix blocking issues):
   - Fix `--extract` → `--fields` globally
   - Fix argument construction for `--ped` and `--annotate-bed`
   - Remove None values from command construction

2. **Test Data Setup**:
   - Create missing test fixture files
   - Verify test VCF has required fields

3. **Pipeline Updates**:
   - Handle missing annotation fields gracefully
   - Update any renamed arguments

4. **Validation**:
   - Run tests individually to verify fixes
   - Add comprehensive logging for debugging

## Typical Working Command Reference

Based on the provided example, a working command should look like:
```bash
variantcentrifuge \
  --vcf-file testing/gckd_all.GRCh37.annotated.vcf.gz \
  --reference GRCh37.75 \
  --gene-name "PKD1,PKD2,CEP290,HNF1B,UMOD,REN,NPHP1,COL4A5,COL4A4,COL4A3" \
  --preset not_artefact \
  --preset 1percent \
  --preset high_or_moderate_or_high_splice_prediction_or_high_prediction \
  --output-file gckd_all.scored.tsv \
  --output-dir testing/test_output \
  --scoring-config-path scoring/nephro_candidate_score \
  --annotate-json-genes testing/ncs/gene_info_summary.json \
  --json-gene-mapping '{"identifier":"symbol","dataFields":["ngs"]}' \
  --inheritance-mode simple \
  --log-level DEBUG \
  --threads 8 \
  --json-genes-as-columns \
  --bcftools-prefilter 'FILTER="PASS" & INFO/AC[0] < 10' \
  --keep-intermediates \
  --enable-checkpoint \
  --log-file testing/test_output/run.log \
  --use-new-pipeline
```

Key differences from test commands:
- Uses `--bcftools-prefilter` instead of `--bcftools-filter`
- Uses proper quoting for filter expressions
- All file paths are properly specified
- No `--extract` or `--fields` in this example (may use defaults)