#!/bin/bash
# Comprehensive Gene Burden Analysis Test Suite
# Generated on 2025-07-09T14:18:13.405087

set -e

echo "=== VariantCentrifuge Gene Burden Test Suite ==="
echo "Test data directory: /mnt/c/development/scholl-lab/variantcentrifuge/tests/fixtures/geneburden/comprehensive_test_data"
echo

# Test 1: Direct sample specification
echo "Test 1: Direct sample specification"
variantcentrifuge \
    --vcf-file test_data.vcf.gz \
    --gene-file test_genes.txt \
    --case-samples CASE_001,CASE_002,CASE_003,CASE_004,CASE_005 \
    --control-samples CTRL_001,CTRL_002,CTRL_003,CTRL_004,CTRL_005,CTRL_006,CTRL_007 \
    --perform-gene-burden \
    --preset rare,coding \
    --output-dir test1_direct_samples \
    --output-file test1_results.tsv \
    --use-new-pipeline

echo "✓ Test 1 completed"
echo

# Test 2: Sample files
echo "Test 2: Sample file specification"
variantcentrifuge \
    --vcf-file test_data.vcf.gz \
    --gene-file test_genes.txt \
    --case-samples-file case_samples.txt \
    --control-samples-file control_samples.txt \
    --perform-gene-burden \
    --preset rare,coding \
    --output-dir test2_sample_files \
    --output-file test2_results.tsv \
    --use-new-pipeline

echo "✓ Test 2 completed"
echo

# Test 3: HPO-based classification
echo "Test 3: HPO-based classification"
variantcentrifuge \
    --vcf-file test_data.vcf.gz \
    --gene-file test_genes.txt \
    --phenotype-file phenotypes_basic.csv \
    --phenotype-sample-column SampleID \
    --phenotype-value-column identifier \
    --case-phenotypes HP:0000113,HP:0000003,HP:0000107,HP:0000822,HP:0003774 \
    --control-phenotypes HP:0000001,HP:0032101 \
    --perform-gene-burden \
    --preset rare,coding \
    --output-dir test3_hpo_classification \
    --output-file test3_results.tsv \
    --use-new-pipeline

echo "✓ Test 3 completed"
echo

# Test 4: HPO term files
echo "Test 4: HPO term files"
variantcentrifuge \
    --vcf-file test_data.vcf.gz \
    --gene-file test_genes.txt \
    --phenotype-file phenotypes_basic.csv \
    --phenotype-sample-column SampleID \
    --phenotype-value-column identifier \
    --case-phenotypes-file case_hpo_terms.txt \
    --control-phenotypes-file control_hpo_terms.txt \
    --perform-gene-burden \
    --preset rare,coding \
    --output-dir test4_hpo_files \
    --output-file test4_results.tsv \
    --use-new-pipeline

echo "✓ Test 4 completed"
echo

# Test 5: Alternative column names
echo "Test 5: Alternative column names"
variantcentrifuge \
    --vcf-file test_data.vcf.gz \
    --gene-file test_genes.txt \
    --phenotype-file phenotypes_alt_columns.csv \
    --phenotype-sample-column sample_name \
    --phenotype-value-column hpo_id \
    --case-phenotypes HP:0000113,HP:0000003,HP:0000107,HP:0000822,HP:0003774 \
    --control-phenotypes HP:0000001,HP:0032101 \
    --perform-gene-burden \
    --preset rare,coding \
    --output-dir test5_alt_columns \
    --output-file test5_results.tsv \
    --use-new-pipeline

echo "✓ Test 5 completed"
echo

echo "=== All tests completed! ==="
echo "Check individual test directories for results:"
echo "- test1_direct_samples/"
echo "- test2_sample_files/"  
echo "- test3_hpo_classification/"
echo "- test4_hpo_files/"
echo "- test5_alt_columns/"
echo

echo "Expected results:"
echo "- Disease genes (PKD1, PKD2, BRCA1, BRCA2) should show enrichment in cases"
echo "- Control genes (TTN, OBSCN, MUC16) should show no significant enrichment"
