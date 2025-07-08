# Regression Test Fix Implementation Plan

## Step-by-Step Fix Implementation

### Step 1: Fix Immediate Blocking Issues

#### 1.1 Fix --extract to --fields
**Files to update**:
- `tests/regression/test_pipeline_regression.py`
- `tests/regression/test_regression_suite.py` 
- `tests/regression/test_baseline_comparison.py`

**Changes**:
```python
# Replace all instances of:
"--extract", fields_str
# With:
"--fields", fields_str
```

#### 1.2 Fix Argument Construction
**File**: `tests/regression/test_regression_suite.py`

**Current** (lines 74-90):
```python
RegressionTestConfig(
    name="gene_with_inheritance",
    gene_name="PKD1",
    extra_args=[
        "--filters",
        "FILTER = 'PASS'",
        "--ped",
        "test_family.ped",  # WRONG: becomes separate argument
        "--inheritance-mode",
        "columns",
    ],
),
```

**Fixed**:
```python
RegressionTestConfig(
    name="gene_with_inheritance",
    gene_name="PKD1",
    extra_args=[
        "--filters",
        "FILTER = 'PASS'",
        "--ped", "tests/fixtures/test_family.ped",  # Keep as two items
        "--inheritance-mode",
        "columns",
    ],
),
```

**Alternative Fix in run() method**:
Update how paths are resolved for relative files.

#### 1.3 Fix None Values in Commands
**File**: `tests/regression/test_pipeline_regression.py`

Add validation before command execution:
```python
# Before subprocess.run() or logging
cmd = [str(arg) for arg in cmd if arg is not None]
```

### Step 2: Create Missing Test Fixtures

#### 2.1 Create test_genes.txt
**File**: `tests/fixtures/test_genes.txt`
```
BRCA1
BRCA2
TP53
PKD1
CFTR
LDLR
APOB
```

#### 2.2 Create test_family.ped
**File**: `tests/fixtures/test_family.ped`
```
FAM001	father	0	0	1	1
FAM001	mother	0	0	2	1
FAM001	child	father	mother	1	2
```

#### 2.3 Create test_regions.bed
**File**: `tests/fixtures/test_regions.bed`
```
chr17	41196312	41277500	BRCA1	0	+
chr13	32889611	32973805	BRCA2	0	+
chr17	7571720	7590868	TP53	0	+
```

#### 2.4 Create or Link Scoring Config
Either create minimal config or update test to use existing:
```bash
# Option 1: Link to existing
ln -s scoring/inheritance_score tests/fixtures/scoring/inheritance_score

# Option 2: Create minimal config
mkdir -p tests/fixtures/scoring/inheritance_score
# Create variable_assignment_config.json and formula_config.json
```

### Step 3: Fix Filter Compatibility

#### 3.1 Check Test VCF Fields
```bash
# Check what fields are available in test VCF
bcftools view -h tests/fixtures/test_regression_variants.GRCh37.annotated.vcf.gz | grep "##INFO"
```

#### 3.2 Update Filter Expressions
If gnomAD fields are missing, update config.json or test filters:
```python
# Option 1: Use simpler filters for tests
"rare": "AF < 0.01"  # Instead of complex gnomAD filters

# Option 2: Add fallback in pipeline
"(dbNSFP_gnomAD_exomes_AF[0] < 0.0001) | (! exists dbNSFP_gnomAD_exomes_AF)"
```

### Step 4: Fix Removed Arguments

#### 4.1 Update --replace-genotypes
Check if functionality moved to:
- `--no-replacement` (inverse logic)
- Always enabled by default

Update tests accordingly:
```python
# If removed entirely, remove from tests
# If renamed, update to new name
```

### Step 5: Update Path Resolution

#### 5.1 Fix Relative Paths in Tests
**File**: `tests/regression/test_regression_suite.py`

Add path resolution in `run()` method:
```python
def run(self, vcf_file: str, output_dir: Path, config: RegressionTestConfig):
    # ... existing code ...
    
    # Fix relative paths in extra_args
    fixed_args = []
    i = 0
    while i < len(config.extra_args):
        arg = config.extra_args[i]
        if arg in ["--ped", "--annotate-bed", "--scoring-config-path", "--gene-file"]:
            # Next arg is a file path
            if i + 1 < len(config.extra_args):
                file_path = config.extra_args[i + 1]
                if not Path(file_path).is_absolute():
                    # Make absolute relative to project root
                    abs_path = project_root / file_path
                    if abs_path.exists():
                        fixed_args.extend([arg, str(abs_path)])
                    else:
                        # Try relative to tests/fixtures
                        abs_path = project_root / "tests" / "fixtures" / file_path
                        fixed_args.extend([arg, str(abs_path)])
                    i += 2
                    continue
        fixed_args.append(arg)
        i += 1
    
    cmd.extend(fixed_args)
```

### Step 6: Add Better Error Handling

#### 6.1 Add Debug Logging
```python
# In test files, add comprehensive logging
logger.debug(f"Working directory: {os.getcwd()}")
logger.debug(f"Command parts: {[repr(x) for x in cmd]}")
logger.debug(f"Environment: {os.environ.get('PATH', '')}")
```

#### 6.2 Add Pre-flight Checks
```python
def check_test_environment():
    """Verify test environment is properly set up."""
    required_files = [
        "tests/fixtures/test_regression_variants.GRCh37.annotated.vcf.gz",
        "tests/fixtures/test_genes.txt",
        "tests/fixtures/test_family.ped",
        "tests/fixtures/test_regions.bed"
    ]
    
    for file in required_files:
        if not Path(file).exists():
            raise FileNotFoundError(f"Required test file missing: {file}")
    
    # Check for required tools
    for tool in ["bcftools", "snpEff", "SnpSift"]:
        if not shutil.which(tool):
            pytest.skip(f"Required tool '{tool}' not found in PATH")
```

### Step 7: Create Test Helper Script

Create a script to verify fixes work:
```python
#!/usr/bin/env python3
# test_single_regression.py
import sys
import subprocess
from pathlib import Path

def test_single_command():
    """Test a single variantcentrifuge command."""
    cmd = [
        "variantcentrifuge",
        "--vcf-file", "tests/fixtures/test_regression_variants.GRCh37.annotated.vcf.gz",
        "--gene-name", "BRCA1",
        "--filters", "FILTER = 'PASS'",
        "--fields", "CHROM POS REF ALT",
        "--output-file", "/tmp/test_output.tsv",
        "--use-new-pipeline"
    ]
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    print(f"Return code: {result.returncode}")
    print(f"Stdout: {result.stdout}")
    print(f"Stderr: {result.stderr}")
    
    return result.returncode == 0

if __name__ == "__main__":
    success = test_single_command()
    sys.exit(0 if success else 1)
```

## Execution Order

1. **Quick Fixes** (5 minutes):
   - Fix `--extract` → `--fields` globally
   - Add None value filtering
   - Create missing fixture files

2. **Test Individual Components** (10 minutes):
   - Run helper script to test basic command
   - Verify each component works independently

3. **Fix Path Issues** (15 minutes):
   - Update path resolution logic
   - Ensure all files are found correctly

4. **Fix Filter Issues** (20 minutes):
   - Check VCF annotation fields
   - Update filters or add fallbacks

5. **Run Full Test Suite** (5 minutes):
   - Execute regression tests
   - Debug any remaining issues

## Verification Commands

```bash
# Test basic functionality
python test_single_regression.py

# Test specific regression test
pytest tests/regression/test_regression_suite.py::TestPipelineRegression::test_pipeline_output_match[single_gene_basic] -xvs

# Run all regression tests
pytest tests/regression/ -xvs

# Check available VCF fields
bcftools view -h tests/fixtures/test_regression_variants.GRCh37.annotated.vcf.gz | grep "##INFO"
```

## Expected Outcomes

After implementing these fixes:
1. All 20 regression tests should pass
2. Old and new pipeline outputs should match
3. No unrecognized argument errors
4. No None value TypeErrors
5. All test files properly resolved

## Rollback Plan

If fixes cause issues:
1. Revert changes to test files
2. Mark regression tests as skipped temporarily
3. Focus on unit tests for individual components
4. Create simplified regression tests that work