# VariantCentrifuge Refactoring Test Protocol

## Critical Testing Requirements

**GOLDEN RULE**: No functionality can be broken. Every refactoring step MUST produce identical output to the original code.

## Test Environment Setup

### 1. Activate Correct Environment
```bash
mamba activate annotation
```

### 2. Verify Environment
```bash
# Check that all tools are available
which variantcentrifuge
which bcftools
which snpEff
which SnpSift
which bedtools
```

## Baseline Test Creation

### Step 1: Create Baseline Output Directory
```bash
# Create directories
mkdir -p testing/test_output_baseline
mkdir -p testing/test_output_current
mkdir -p testing/test_scripts
```

### Step 2: Run Comprehensive Baseline Test
```bash
#!/bin/bash
# Save as: testing/test_scripts/run_comprehensive_test.sh

# Set test name
TEST_NAME="${1:-comprehensive}"
OUTPUT_DIR="${2:-testing/test_output_current}"

# Clear output directory
rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Run the comprehensive test
variantcentrifuge \
  --vcf-file testing/gckd_all.GRCh37.annotated.vcf.gz \
  --reference GRCh37.75 \
  --gene-name "PKD1,PKD2,CEP290,HNF1B,UMOD,REN,NPHP1,COL4A5,COL4A4,COL4A3" \
  --preset not_artefact \
  --preset 1percent \
  --preset high_or_moderate_or_high_splice_prediction_or_high_prediction \
  --output-file gckd_all.scored.tsv \
  --output-dir "${OUTPUT_DIR}" \
  --scoring-config-path scoring/nephro_candidate_score \
  --annotate-json-genes testing/ncs/gene_info_summary.json \
  --json-gene-mapping '{"identifier":"symbol","dataFields":["ngs"]}' \
  --inheritance-mode simple \
  --log-level DEBUG \
  --threads 8 \
  --json-genes-as-columns \
  --bcftools-prefilter 'FILTER="PASS" & INFO/AC[0] < 10' \
  --enable-checkpoint \
  --log-file "${OUTPUT_DIR}/run.log" \
  --keep-intermediates

# Check exit code
if [ $? -eq 0 ]; then
    echo "SUCCESS: Pipeline completed successfully"
else
    echo "FAILURE: Pipeline failed"
    exit 1
fi
```

### Step 3: Create Initial Baseline
```bash
# Make script executable
chmod +x testing/test_scripts/run_comprehensive_test.sh

# Run baseline test
./testing/test_scripts/run_comprehensive_test.sh baseline testing/test_output_baseline

# Create checksums for all output files
find testing/test_output_baseline -type f -name "*.tsv" -o -name "*.xlsx" -o -name "*.json" | \
  sort | xargs md5sum > testing/baseline_checksums.txt

# Create file listing with sizes
find testing/test_output_baseline -type f -exec ls -la {} \; > testing/baseline_files.txt
```

## Automated Regression Testing

### Create Regression Test Script
```python
#!/usr/bin/env python3
# Save as: testing/test_scripts/regression_test.py

import subprocess
import filecmp
import difflib
import hashlib
import sys
from pathlib import Path
import json

class RegressionTester:
    def __init__(self, baseline_dir="testing/test_output_baseline", 
                 test_dir="testing/test_output_current"):
        self.baseline_dir = Path(baseline_dir)
        self.test_dir = Path(test_dir)
        self.failures = []
        
    def run_pipeline(self):
        """Run the pipeline with current code."""
        print("Running pipeline with current code...")
        cmd = ["./testing/test_scripts/run_comprehensive_test.sh", 
               "current", str(self.test_dir)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILURE: Pipeline execution failed")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            return False
        
        print("Pipeline completed successfully")
        return True
    
    def compare_file_checksums(self, file1, file2):
        """Compare files using MD5 checksums."""
        def get_checksum(filepath):
            hash_md5 = hashlib.md5()
            with open(filepath, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_md5.update(chunk)
            return hash_md5.hexdigest()
        
        return get_checksum(file1) == get_checksum(file2)
    
    def show_diff(self, file1, file2, max_lines=50):
        """Show differences between two text files."""
        try:
            with open(file1, 'r') as f1, open(file2, 'r') as f2:
                lines1 = f1.readlines()
                lines2 = f2.readlines()
                
                diff = list(difflib.unified_diff(
                    lines1[:max_lines], 
                    lines2[:max_lines],
                    fromfile=str(file1),
                    tofile=str(file2),
                    lineterm=''
                ))
                
                if diff:
                    print("\n".join(diff[:100]))  # Show first 100 lines of diff
                    if len(diff) > 100:
                        print(f"... and {len(diff) - 100} more lines")
                        
        except Exception as e:
            print(f"Could not show diff: {e}")
    
    def compare_json_files(self, file1, file2):
        """Compare JSON files semantically."""
        try:
            with open(file1, 'r') as f1, open(file2, 'r') as f2:
                data1 = json.load(f1)
                data2 = json.load(f2)
                return data1 == data2
        except:
            return False
    
    def compare_outputs(self):
        """Compare all output files between baseline and current."""
        print("\nComparing outputs...")
        
        # Get all files from baseline
        baseline_files = set()
        for f in self.baseline_dir.rglob("*"):
            if f.is_file() and not f.name.endswith('.log'):
                baseline_files.add(f.relative_to(self.baseline_dir))
        
        # Check each baseline file
        for rel_path in sorted(baseline_files):
            baseline_file = self.baseline_dir / rel_path
            test_file = self.test_dir / rel_path
            
            if not test_file.exists():
                self.failures.append(f"MISSING: {rel_path}")
                continue
            
            # Special handling for different file types
            if rel_path.suffix == '.json':
                if not self.compare_json_files(baseline_file, test_file):
                    self.failures.append(f"JSON DIFF: {rel_path}")
            elif rel_path.suffix in ['.tsv', '.txt']:
                if not filecmp.cmp(baseline_file, test_file, shallow=False):
                    self.failures.append(f"CONTENT DIFF: {rel_path}")
                    self.show_diff(baseline_file, test_file)
            else:
                # Binary comparison for other files
                if not self.compare_file_checksums(baseline_file, test_file):
                    self.failures.append(f"CHECKSUM DIFF: {rel_path}")
        
        # Check for extra files in test output
        test_files = set()
        for f in self.test_dir.rglob("*"):
            if f.is_file() and not f.name.endswith('.log'):
                test_files.add(f.relative_to(self.test_dir))
        
        extra_files = test_files - baseline_files
        for f in extra_files:
            self.failures.append(f"EXTRA FILE: {f}")
        
        return len(self.failures) == 0
    
    def print_summary(self):
        """Print test summary."""
        print("\n" + "="*70)
        if self.failures:
            print(f"REGRESSION TEST FAILED - {len(self.failures)} issues found:")
            for failure in self.failures:
                print(f"  - {failure}")
            print("="*70)
        else:
            print("REGRESSION TEST PASSED - All outputs match baseline!")
            print("="*70)
    
    def run(self):
        """Run complete regression test."""
        # Run pipeline
        if not self.run_pipeline():
            return False
        
        # Compare outputs
        success = self.compare_outputs()
        
        # Print summary
        self.print_summary()
        
        return success

if __name__ == "__main__":
    tester = RegressionTester()
    success = tester.run()
    sys.exit(0 if success else 1)
```

### Make Scripts Executable
```bash
chmod +x testing/test_scripts/regression_test.py
```

## Testing Protocol for Each Refactoring Step

### Before ANY Code Change

1. **Verify Clean State**
   ```bash
   # Ensure no uncommitted changes
   git status
   
   # Run regression test to ensure current state is good
   python testing/test_scripts/regression_test.py
   ```

2. **Create Feature Branch**
   ```bash
   git checkout -b refactor/step-X-description
   ```

### After EACH Code Change

1. **Run Quick Smoke Test** (if you have one)
   ```bash
   # Quick test with single gene
   variantcentrifuge \
     --vcf-file testing/small_test.vcf \
     --gene-name BRCA1 \
     --output-dir testing/smoke_test
   ```

2. **Run Full Regression Test**
   ```bash
   python testing/test_scripts/regression_test.py
   ```

3. **If Test Fails**
   ```bash
   # Check specific differences
   diff -u testing/test_output_baseline/gckd_all.scored.tsv \
           testing/test_output_current/gckd_all.scored.tsv | head -100
   
   # Check intermediate files if needed
   diff -u testing/test_output_baseline/intermediate/gckd_all.extracted.tsv \
           testing/test_output_current/intermediate/gckd_all.extracted.tsv
   ```

4. **Commit ONLY If Tests Pass**
   ```bash
   git add -p  # Review changes carefully
   git commit -m "refactor: [description of change]"
   ```

## Additional Test Scenarios

### Test 1: Single Gene, Minimal Options
```bash
variantcentrifuge \
  --vcf-file testing/test_data.vcf \
  --gene-name BRCA1 \
  --output-dir testing/test_single_gene
```

### Test 2: Multi-threading
```bash
variantcentrifuge \
  --vcf-file testing/test_data.vcf \
  --gene-file testing/gene_list.txt \
  --threads 4 \
  --output-dir testing/test_parallel
```

### Test 3: All Features Enabled
```bash
variantcentrifuge \
  --vcf-file testing/test_data.vcf \
  --gene-file testing/gene_list.txt \
  --preset rare,coding,pathogenic \
  --xlsx \
  --html-report \
  --perform-gene-burden \
  --case-samples-file testing/cases.txt \
  --control-samples-file testing/controls.txt \
  --calculate-inheritance \
  --ped testing/family.ped \
  --output-dir testing/test_full_features
```

### Test 4: Checkpoint/Resume
```bash
# First run with checkpoint
variantcentrifuge \
  --vcf-file testing/test_data.vcf \
  --gene-file testing/gene_list.txt \
  --enable-checkpoint \
  --output-dir testing/test_checkpoint

# Kill process halfway through (Ctrl+C)

# Resume
variantcentrifuge \
  --vcf-file testing/test_data.vcf \
  --gene-file testing/gene_list.txt \
  --enable-checkpoint \
  --resume \
  --output-dir testing/test_checkpoint
```

## Performance Testing

### Create Performance Benchmark
```bash
#!/bin/bash
# Save as: testing/test_scripts/performance_test.sh

echo "Running performance benchmark..."

# Time the baseline
echo "Baseline performance:"
time ./testing/test_scripts/run_comprehensive_test.sh baseline testing/perf_baseline

# Time the current code
echo "Current performance:"
time ./testing/test_scripts/run_comprehensive_test.sh current testing/perf_current

# Compare file sizes (should be identical)
echo "Output file sizes:"
ls -la testing/perf_baseline/*.tsv
ls -la testing/perf_current/*.tsv
```

## Emergency Rollback Procedure

If something goes wrong:

1. **Immediate Rollback**
   ```bash
   git checkout main
   git branch -D refactor/step-X-description
   ```

2. **Verify Rollback**
   ```bash
   python testing/test_scripts/regression_test.py
   ```

3. **Document Issue**
   - Note what change caused the failure
   - Save diff output for analysis
   - Update refactoring plan with lessons learned

## Test Maintenance

### Weekly Tasks
1. Re-run baseline to ensure it still works
2. Update test data if needed
3. Add new test cases for any bugs found

### Before Major Refactoring Steps
1. Create new baseline with current main branch
2. Run all test scenarios
3. Document current performance metrics

## Success Criteria Checklist

Before considering ANY refactoring step complete:

- [ ] Regression test passes (identical output)
- [ ] No performance degradation (±5%)
- [ ] All test scenarios pass
- [ ] Code review completed
- [ ] Documentation updated
- [ ] Changes committed with clear message

Remember: **Patient, methodical testing prevents disasters!**