# Fix for Missing VAR_ID and Custom_Annotation in New Pipeline

## Issues Identified

1. **VAR_ID column missing** - The unique variant identifiers were not appearing in the new pipeline output
2. **Custom_Annotation column missing** - Custom annotations from BED files, gene lists, and JSON data were lost
3. **Count columns are present** - These are working correctly (proband_count, control_count, etc.)
4. **ID field shows "NA"** - This is expected behavior for empty fields (SnpSift configuration)

## Root Cause

The `VariantAnalysisStage` creates a completely new DataFrame from the `analyze_variants` function output, discarding any columns that were added by previous stages:

1. `CustomAnnotationStage` runs and adds `Custom_Annotation` column
2. `VariantAnalysisStage` runs and creates a new DataFrame, losing `Custom_Annotation`
3. `VariantIdentifierStage` was supposed to run before output, but the DataFrame it would modify had already been replaced

## Fixes Applied

### 1. Fixed VariantAnalysisStage Column Preservation

Updated `analysis_stages.py` to properly merge columns from the original DataFrame:
- Uses CHROM, POS, REF, ALT as key columns for merging
- Preserves any columns that exist in the original DataFrame but not in the analysis output
- Handles cases where analyze_variants might filter or reorder rows
- Properly positions VAR_ID as first column and Custom_Annotation after GT column

### 2. Fixed Stage Ordering

Updated `pipeline_refactored.py` to run stages in the correct order:
- Moved `VariantIdentifierStage` to run AFTER `VariantAnalysisStage`
- This ensures VAR_ID is added to the final DataFrame that will be output

## Code Changes

### analysis_stages.py (VariantAnalysisStage._process)
```python
# Old: Simple assignment that loses columns
analysis_df.insert(0, "VAR_ID", df["VAR_ID"].values)

# New: Proper merge that handles filtered/reordered data
key_cols = ["CHROM", "POS", "REF", "ALT"]
if all(col in df.columns for col in key_cols) and all(col in analysis_df.columns for col in key_cols):
    original_only_cols = [col for col in df.columns if col not in analysis_df.columns]
    if original_only_cols:
        merge_df = df[key_cols + original_only_cols].copy()
        analysis_df = pd.merge(analysis_df, merge_df, on=key_cols, how="left")
```

### pipeline_refactored.py (build_pipeline_stages)
```python
# Old order:
stages.append(VariantAnalysisStage())
stages.append(StatisticsGenerationStage())
stages.append(GeneBurdenAnalysisStage())
stages.append(VariantIdentifierStage())  # Too late!

# New order:
stages.append(VariantAnalysisStage())
stages.append(VariantIdentifierStage())  # Right after analysis
stages.append(StatisticsGenerationStage())
stages.append(GeneBurdenAnalysisStage())
```

## Expected Output

After these fixes, the new pipeline output should include:
1. **VAR_ID** as the first column with format `var_XXXX_YYYY`
2. **Custom_Annotation** column after GT (when annotations are used)
3. All count columns (already working)
4. Same column order as the old pipeline

## Testing

To verify the fix works:
```bash
# Run with new pipeline
variantcentrifuge --use-new-pipeline \
  --vcf-file test.vcf \
  --gene-file genes.txt \
  --annotate-bed regions.bed \
  --output-file output.tsv

# Check output has VAR_ID and Custom_Annotation
head -1 output.tsv | grep -E "VAR_ID|Custom_Annotation"
```