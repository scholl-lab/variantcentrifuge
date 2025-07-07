#!/bin/bash
# Test script to compare old and new pipeline outputs

set -e

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate annotation

echo "=== Testing VariantCentrifuge Pipeline Comparison ==="
echo "Testing with BRCA1 gene..."

# Create output directories
mkdir -p test_output_old test_output_new

# Test parameters
VCF="testing/all_merged.ann.dbnsfp.BRCA1.vcf.gz"
GENE="BRCA1"
REFERENCE="GRCh37.75"
FIELDS="CHROM POS REF ALT Gene_Name FILTER AF_gnomAD"

echo ""
echo "1. Running OLD pipeline..."
time variantcentrifuge \
  --vcf-file "$VCF" \
  --gene-name "$GENE" \
  --output-file test_old.tsv \
  --output-dir test_output_old \
  --reference "$REFERENCE" \
  --preset pass \
  --preset rare \
  -e "$FIELDS" \
  --log-level INFO

echo ""
echo "2. Running NEW pipeline..."
time variantcentrifuge \
  --vcf-file "$VCF" \
  --gene-name "$GENE" \
  --output-file test_new.tsv \
  --output-dir test_output_new \
  --reference "$REFERENCE" \
  --preset pass \
  --preset rare \
  -e "$FIELDS" \
  --log-level INFO \
  --use-new-pipeline

echo ""
echo "3. Comparing outputs..."

# Check if files exist
if [ ! -f "test_output_old/test_old.tsv" ]; then
    echo "ERROR: Old pipeline output not found!"
    exit 1
fi

if [ ! -f "test_output_new/test_new.tsv" ]; then
    echo "ERROR: New pipeline output not found!"
    exit 1
fi

# Compare line counts
OLD_LINES=$(wc -l < test_output_old/test_old.tsv)
NEW_LINES=$(wc -l < test_output_new/test_new.tsv)

echo "Old pipeline: $OLD_LINES lines"
echo "New pipeline: $NEW_LINES lines"

# Show first few lines of each
echo ""
echo "Old pipeline output (first 5 lines):"
head -5 test_output_old/test_old.tsv

echo ""
echo "New pipeline output (first 5 lines):"
head -5 test_output_new/test_new.tsv

# Detailed comparison
echo ""
echo "Checking differences..."
if diff -q test_output_old/test_old.tsv test_output_new/test_new.tsv > /dev/null; then
    echo "SUCCESS: Outputs are identical!"
else
    echo "WARNING: Outputs differ. Showing first 10 differences:"
    diff test_output_old/test_old.tsv test_output_new/test_new.tsv | head -20 || true
fi

echo ""
echo "=== Test Complete ==="