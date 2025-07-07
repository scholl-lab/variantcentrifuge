#!/bin/bash
# Test the new pipeline architecture

set -e

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate annotation

echo "=== Testing New Pipeline Architecture ==="

# Create output directory
mkdir -p test_output_new

# Test parameters
VCF="testing/all_merged.ann.dbnsfp.BRCA1.vcf.gz"
GENE="BRCA1"
REFERENCE="GRCh37.75"
FIELDS="CHROM POS REF ALT Gene_Name FILTER AF_gnomAD"

echo "Running new pipeline with debug logging..."
python -m variantcentrifuge.cli \
  --vcf-file "$VCF" \
  --gene-name "$GENE" \
  --output-file test_new.tsv \
  --output-dir test_output_new \
  --reference "$REFERENCE" \
  --preset pass \
  --preset rare \
  -e "$FIELDS" \
  --log-level DEBUG \
  --use-new-pipeline \
  2>&1 | tee new_pipeline_debug.log

echo ""
echo "Checking output..."
if [ -f "test_output_new/test_new.tsv" ]; then
    echo "SUCCESS: Output file created"
    echo "Line count: $(wc -l < test_output_new/test_new.tsv)"
    echo ""
    echo "First 5 lines:"
    head -5 test_output_new/test_new.tsv
else
    echo "ERROR: Output file not found!"
    exit 1
fi