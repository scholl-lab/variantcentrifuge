# VariantCentrifuge Refactoring - Quick Action Checklist

## 🚀 Immediate Actions (Day 1)

### 1. Install Tools
```bash
mamba create -n vc-test python=3.8
mamba activate vc-test
mamba install -c bioconda bcftools=1.17 snpeff=5.1 snpsift=5.1 bedtools=2.31
pip install -e .
```

### 2. Generate Test Data
```bash
# Create test VCF if needed
python tests/fixtures/create_test_vcf.py

# Or use existing test data
ls tests/data/*.vcf.gz
```

## 🧪 Validation (Days 2-3)

### 3. Run Regression Tests
```bash
# Generate baseline with old pipeline
python tests/regression/generate_baseline.py

# Run regression suite
pytest tests/regression/test_regression_suite.py -v

# Check for differences
python tests/regression/analyze_differences.py
```

### 4. Performance Testing
```bash
# Run benchmarks
python tests/performance/benchmark_pipeline.py

# View results
python tests/performance/visualize_results.py
```

## 📝 Documentation (Day 4)

### 5. Update Docs
- [ ] Update `docs/source/architecture.md` with Stage-based design
- [ ] Add `docs/source/migration.md` for users
- [ ] Update `README.md` with new architecture note
- [ ] Update `CHANGELOG.md` with refactoring details

## ✅ Final Steps (Day 5)

### 6. Production Testing
```bash
# Test with real production data
variantcentrifuge --gene-file production_genes.txt \
                 --vcf-file production.vcf.gz \
                 --use-new-pipeline \
                 --output-dir test_new/

# Compare with old pipeline
variantcentrifuge --gene-file production_genes.txt \
                 --vcf-file production.vcf.gz \
                 --output-dir test_old/

# Diff outputs
diff -r test_old/ test_new/
```

### 7. Sign-off Checklist
- [ ] All regression tests pass
- [ ] Performance improvement confirmed (≥50%)
- [ ] Documentation reviewed and approved
- [ ] Production data processed successfully
- [ ] Stakeholder approval received

## 🎯 Definition of Done

The refactoring is complete when:
1. ✅ New pipeline produces identical outputs to old pipeline
2. ✅ Performance improvement ≥50% on multi-core systems
3. ✅ All tests pass (unit, integration, regression)
4. ✅ Documentation is updated
5. ✅ Successfully processes production workloads

## 📞 Contact

Questions? Check:
- Technical details: `refactoring/FINAL_COMPLETION_PLAN.md`
- Architecture: `refactoring/README.md`
- Test coverage: `refactoring/TESTING_COVERAGE_REVIEW.md`