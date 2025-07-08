# Analysis Stages Review: DataFrame Modifications and Parallel Safety

## Summary of Findings

### Stages That Modify context.current_dataframe

1. **DataFrameLoadingStage** (line 31)
   - **Modifies DataFrame**: YES - Sets `context.current_dataframe = df` (line 123)
   - **parallel_safe**: Not explicitly set (defaults to False)
   - **Type of modification**: Initial loading, not in-place modification

2. **CustomAnnotationStage** (line 143)
   - **Modifies DataFrame**: YES - Sets `context.current_dataframe = df` (line 198)
   - **parallel_safe**: Not explicitly set (defaults to False)
   - **Type of modification**: Replaces entire DataFrame with annotated version

3. **InheritanceAnalysisStage** (line 202)
   - **Modifies DataFrame**: YES - Sets `context.current_dataframe = df` (line 273)
   - **parallel_safe**: Not explicitly set (defaults to False)
   - **Type of modification**: Replaces entire DataFrame with inheritance analysis results

4. **VariantScoringStage** (line 277)
   - **Modifies DataFrame**: YES - Sets `context.current_dataframe = df` (line 325)
   - **parallel_safe**: Not explicitly set (defaults to False)
   - **Type of modification**: Replaces entire DataFrame with scored version

5. **StatisticsGenerationStage** (line 329)
   - **Modifies DataFrame**: NO - Only reads DataFrame, doesn't modify it
   - **parallel_safe**: TRUE (line 350-353) - Explicitly marked as safe
   - **Type of modification**: None - read-only computation

6. **VariantAnalysisStage** (line 411)
   - **Modifies DataFrame**: YES - Sets `context.current_dataframe = analysis_df` (line 511)
   - **parallel_safe**: Not explicitly set (defaults to False)
   - **Type of modification**: Creates new DataFrame from analyze_variants results

7. **GeneBurdenAnalysisStage** (line 521)
   - **Modifies DataFrame**: NO - Only reads DataFrame for analysis
   - **parallel_safe**: Not explicitly set (defaults to False)
   - **Type of modification**: None - stores results separately in context.gene_burden_results

8. **ChunkedAnalysisStage** (line 596)
   - **Modifies DataFrame**: YES - Sets `context.current_dataframe = pd.concat(...)` (line 699)
   - **parallel_safe**: Not explicitly set (defaults to False)
   - **Type of modification**: Builds new DataFrame from processed chunks

9. **ParallelAnalysisOrchestrator** (line 851)
   - **Modifies DataFrame**: YES - Sets `context.current_dataframe = pd.concat(...)` (line 978)
   - **parallel_safe**: FALSE (line 870-874) - Explicitly marked as not safe
   - **Type of modification**: Manages its own parallelism, rebuilds DataFrame

## Critical Issue Identified

**StatisticsGenerationStage is marked as parallel_safe = True, but this could be problematic if:**
- It runs in a separate process where it cannot access the modified DataFrame from other stages
- The DataFrame modifications from other stages haven't been completed yet

## Recommendations

1. **Remove parallel_safe = True from StatisticsGenerationStage** - Even though it's read-only, it depends on having access to the current DataFrame which may not be available in a parallel process.

2. **All stages that modify context.current_dataframe should have parallel_safe = False** (which is the default).

3. **Consider creating a separate read-only context for parallel stages** that contains a snapshot of the DataFrame at the time of parallel execution.

## Code Pattern Analysis

All DataFrame-modifying stages follow this pattern:
```python
df = context.current_dataframe  # Get reference
df = some_operation(df)         # Create new DataFrame
context.current_dataframe = df  # Replace in context
```

This is NOT in-place modification (which would be like `df.insert()` or `df['new_col'] = ...`), but rather replacing the entire DataFrame reference. However, this still cannot work correctly in parallel processes because:
- Each process gets a copy of the context
- Changes to context.current_dataframe in a subprocess won't propagate back to the main process
- Other parallel stages won't see these changes

## Conclusion

Only **StatisticsGenerationStage** has `parallel_safe = True`, and while it's read-only, it should probably be changed to `False` to ensure it has access to the complete, modified DataFrame from previous stages.