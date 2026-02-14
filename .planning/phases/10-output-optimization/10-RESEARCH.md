# Phase 10: Output Optimization - Research

**Researched:** 2026-02-14
**Domain:** Excel generation and genotype parsing optimization
**Confidence:** HIGH

## Summary

Phase 10 focuses on accelerating Excel generation (2-5x) and eliminating redundant GT column parsing across the pipeline. Current implementation uses openpyxl exclusively for Excel generation, which is 2-5x slower than xlsxwriter for initial writes. GT parsing happens repeatedly (in converter.py for IGV links, in JSON generation, and potentially in gene burden/statistics stages) using regex pattern `r"([^()]+)\(([^)]+)\)"` to extract Sample(GT) tuples.

The standard approach is a **two-pass Excel strategy**: xlsxwriter for fast initial write, then openpyxl for finalization (hyperlinks, freeze panes, auto-filters). GT pre-parsing should happen once at DataFrame load time in `dataframe_optimizer.py:load_optimized_dataframe()`, storing parsed data as additional columns for downstream consumption.

**Primary recommendation:** Use xlsxwriter via `pd.ExcelWriter(engine='xlsxwriter')` for initial write, then `openpyxl.load_workbook()` to add hyperlinks/formatting. Pre-parse GT column during DataFrame load into structured cache columns (e.g., `_GT_PARSED` dict column).

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| xlsxwriter | 3.x | Fast Excel write | Industry standard for write-only Excel generation; 2-5x faster than openpyxl for large DataFrames |
| openpyxl | 3.1.x | Excel read/modify/format | De facto standard for reading and modifying existing .xlsx files; required for hyperlinks/formatting |
| pandas | 2.x+ | DataFrame Excel I/O | Integrates seamlessly with both xlsxwriter and openpyxl via ExcelWriter engine parameter |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pyarrow | 15.x+ | Fast CSV read (already in use) | Already used in Phase 8 for 5-15x faster DataFrame loading |
| NumPy | 1.x | Structured arrays for GT cache | If GT cache uses NumPy structured arrays instead of dict columns |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| xlsxwriter+openpyxl | openpyxl only | Simpler but 2-5x slower for initial write on large datasets |
| xlsxwriter+openpyxl | xlsxwriter only | Faster but cannot add hyperlinks to existing files (xlsxwriter is write-only) |

**Installation:**
```bash
# Already installed - verify versions
pip show xlsxwriter openpyxl pandas
```

## Architecture Patterns

### Recommended Project Structure
```
variantcentrifuge/
├── converter.py                # Excel generation functions
│   ├── convert_to_excel()      # MODIFY: Use xlsxwriter engine
│   ├── finalize_excel_file()   # KEEP: openpyxl for hyperlinks/formatting
│   └── append_tsv_as_sheet()   # MODIFY: Use xlsxwriter for speed
├── dataframe_optimizer.py      # DataFrame loading optimizations
│   ├── load_optimized_dataframe()  # ADD: GT pre-parsing logic
│   └── parse_gt_column()           # NEW: GT parsing function
└── stages/
    └── output_stages.py        # ExcelReportStage
        └── ExcelReportStage._process()  # MODIFY: Pass df to avoid disk read
```

### Pattern 1: Two-Pass Excel Generation (xlsxwriter → openpyxl)
**What:** Write initial Excel file with xlsxwriter for speed, then load with openpyxl to add features xlsxwriter can't handle (hyperlinks to existing cells, openpyxl-specific formatting)

**When to use:** Any Excel generation where speed matters (>1000 rows) AND hyperlinks/advanced formatting are required

**Example:**
```python
# Source: https://xlsxwriter.readthedocs.io/working_with_pandas.html
# Phase 1: Fast write with xlsxwriter
with pd.ExcelWriter('output.xlsx', engine='xlsxwriter') as writer:
    df.to_excel(writer, sheet_name='Results', index=False)

    # Can add multiple sheets here (Metadata, Statistics, Gene Burden)
    metadata_df.to_excel(writer, sheet_name='Metadata', index=False)
    stats_df.to_excel(writer, sheet_name='Statistics', index=False)

# Phase 2: Finalize with openpyxl (hyperlinks, freeze panes, auto-filters)
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

wb = load_workbook('output.xlsx')
ws = wb['Results']

# Freeze top row
ws.freeze_panes = "A2"

# Auto-filter
max_col_letter = get_column_letter(ws.max_column)
ws.auto_filter.ref = f"A1:{max_col_letter}1"

# Add hyperlinks
for row_idx in range(2, ws.max_row + 1):
    cell = ws.cell(row=row_idx, column=5)  # Example column
    cell.hyperlink = "https://example.com"
    cell.style = "Hyperlink"

wb.save('output.xlsx')
```

### Pattern 2: GT Pre-Parsing at DataFrame Load Time
**What:** Parse GT column once during `load_optimized_dataframe()`, store parsed data as additional DataFrame columns for all downstream stages to use

**When to use:** When GT parsing is needed by multiple stages (converter.py, statistics, gene burden) - eliminates redundant regex operations

**Example:**
```python
# In dataframe_optimizer.py
def parse_gt_column(df: pd.DataFrame) -> pd.DataFrame:
    """Parse GT column into structured data for downstream use.

    GT format: "Sample1(0/1);Sample2(1/1);Sample3(0/0)"
    Parsed format: [{"sample": "Sample1", "gt": "0/1", "allele1": 0, "allele2": 1}, ...]
    """
    import re
    pattern = re.compile(r"([^()]+)\(([^)]+)\)")

    def parse_row(gt_value):
        if pd.isna(gt_value) or not gt_value:
            return []

        parsed = []
        for entry in str(gt_value).split(";"):
            entry = entry.strip()
            m = pattern.match(entry)
            if m:
                sample_id = m.group(1).strip()
                genotype = m.group(2).strip()

                # Parse alleles
                separator = "|" if "|" in genotype else "/"
                parts = genotype.split(separator)
                allele1 = None if parts[0] == "." else int(parts[0])
                allele2 = None if len(parts) > 1 and parts[1] != "." else None if len(parts) > 1 else allele1
                if len(parts) > 1 and parts[1] != ".":
                    allele2 = int(parts[1])

                parsed.append({
                    "sample": sample_id,
                    "gt": genotype,
                    "allele1": allele1,
                    "allele2": allele2
                })

        return parsed

    # Add parsed GT as new column
    df['_GT_PARSED'] = df['GT'].apply(parse_row)

    return df

# In load_optimized_dataframe(), after loading:
if 'GT' in df.columns and not df.empty:
    logger.info("Pre-parsing GT column for downstream stages")
    df = parse_gt_column(df)
```

### Pattern 3: In-Memory DataFrame Pass-Through for Excel
**What:** Reuse in-memory DataFrame (`context.variants_df`) for Excel generation instead of re-reading from disk

**When to use:** When DataFrame fits in memory (Phase 8's `should_use_memory_passthrough()` already determines this)

**Example:**
```python
# In ExcelReportStage._process()
# CURRENT: Always reads from disk
excel_df = None
if context.variants_df is not None:
    excel_df = context.variants_df.copy()  # Use in-memory copy
    logger.info("Using in-memory DataFrame for Excel (skipping disk read)")

# Then pass to convert_to_excel
xlsx_file = convert_to_excel(str(input_file), context.config, df=excel_df)
```

### Anti-Patterns to Avoid
- **Opening workbook multiple times:** Open once with openpyxl, make all changes, save once. Each open/save cycle is expensive.
- **Row-by-row writes with openpyxl:** Use xlsxwriter or pandas bulk write instead. Openpyxl row-by-row is 10-100x slower.
- **Mixing engines mid-process:** Write entire file with one engine (xlsxwriter), then finalize with another (openpyxl). Don't interleave.
- **Re-parsing GT for each stage:** Parse once at load time, store in DataFrame, consume from cache.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Fast Excel write | Custom XML generation | xlsxwriter via pandas | xlsxwriter handles all OOXML complexity, compression, relationships; 2-5x faster than openpyxl |
| Hyperlinks in existing Excel | String manipulation of XML | openpyxl's cell.hyperlink API | Hyperlinks require relationship IDs, rels files, proper escaping - openpyxl handles all of this |
| GT string parsing | Hand-rolled string splits | Centralized parse function with caching | GT format has edge cases (missing values, phasing); parse once, cache result |
| DataFrame to Excel | Writing TSV then converting | pandas.to_excel() with xlsxwriter engine | pandas optimizes DataFrame→Excel conversion; handles types, NA values, formatting |

**Key insight:** Excel OOXML format is complex (ZIP with XML files, relationships, shared strings). Libraries handle edge cases (special chars, formulas, merged cells, 1-based indexing). Custom solutions miss these and break on real-world data.

## Common Pitfalls

### Pitfall 1: Using openpyxl for Initial Write
**What goes wrong:** Excel generation takes 2-5x longer than necessary for large datasets (10K+ rows)
**Why it happens:** openpyxl is optimized for read/modify workflows, not write-only bulk operations. It maintains in-memory object model which adds overhead.
**How to avoid:** Use xlsxwriter for initial write via `pd.ExcelWriter(engine='xlsxwriter')`, then load with openpyxl only for finalization features
**Warning signs:** Excel generation taking >1 second per 1000 rows; memory usage spiking during write

### Pitfall 2: Not Closing ExcelWriter Before Finalization
**What goes wrong:** openpyxl fails to load the file or loads corrupted data
**Why it happens:** xlsxwriter uses context manager (`with` statement) or requires explicit `.close()` to flush and finalize the ZIP archive. If not closed, file is incomplete.
**How to avoid:** Always use `with pd.ExcelWriter(...) as writer:` pattern, or call `writer.close()` explicitly before opening with openpyxl
**Warning signs:** openpyxl raises `zipfile.BadZipFile` or `InvalidFileException`

### Pitfall 3: Openpyxl Column Indexing (1-Based)
**What goes wrong:** Hyperlinks or formatting applied to wrong columns (off-by-one errors)
**Why it happens:** openpyxl uses 1-based indexing for rows and columns (Excel convention), while pandas/Python use 0-based
**How to avoid:** When iterating over DataFrame column indices, add 1 when passing to openpyxl: `ws.cell(row=row_idx, column=col_idx+1)`
**Warning signs:** First column missing formatting; hyperlinks appearing in wrong cells

### Pitfall 4: Regex Pattern Recompilation
**What goes wrong:** GT parsing is slower than necessary due to repeated regex compilation
**Why it happens:** `re.compile()` called inside loops or per-row functions instead of once at module/function level
**How to avoid:** Compile pattern once (module-level or function-level constant), reuse for all rows
**Warning signs:** Profiling shows significant time in `re.compile()` or `sre_compile.compile()`

### Pitfall 5: Hyperlink Cell Value Replacement
**What goes wrong:** Cell values disappear when hyperlinks are added, or hyperlinks don't display as links
**Why it happens:** Setting `cell.hyperlink` doesn't automatically set display value or style. Need to set both `cell.value` and `cell.style = "Hyperlink"`
**How to avoid:** Always set value, hyperlink, and style together:
```python
cell.value = "Display Text"
cell.hyperlink = "https://example.com"
cell.style = "Hyperlink"
```
**Warning signs:** Cells showing formula instead of link text; links not blue/underlined

### Pitfall 6: GT Cache Column Name Collision
**What goes wrong:** Pre-parsed GT cache column (`_GT_PARSED`) collides with existing column or gets written to output
**Why it happens:** Cache column not marked as internal/temporary; not removed before output
**How to avoid:** Use underscore prefix (`_GT_PARSED`) to indicate internal column; remove in TSVOutputStage before writing
**Warning signs:** TSV/Excel output contains unexpected `_GT_PARSED` column; DataFrame has duplicate columns

## Code Examples

Verified patterns from official sources and current codebase:

### Current Excel Generation (Baseline)
```python
# Source: variantcentrifuge/converter.py (lines 26-77)
def convert_to_excel(tsv_file: str, cfg: dict, df: pd.DataFrame | None = None) -> str:
    """Convert TSV to XLSX with openpyxl engine (current implementation)."""
    if df is not None:
        logger.debug("Using provided DataFrame (skipping disk read)")
    else:
        df = pd.read_csv(tsv_file, sep="\t", na_values="NA", low_memory=False)

    xlsx_file = os.path.splitext(tsv_file)[0] + ".xlsx"
    df.to_excel(xlsx_file, index=False, sheet_name="Results")  # Uses openpyxl by default
    return xlsx_file
```

### Optimized Excel Generation (Two-Pass Strategy)
```python
# Source: https://xlsxwriter.readthedocs.io/working_with_pandas.html
def convert_to_excel_optimized(tsv_file: str, cfg: dict, df: pd.DataFrame | None = None) -> str:
    """Convert TSV to XLSX with xlsxwriter for speed."""
    if df is not None:
        logger.debug("Using provided DataFrame (skipping disk read)")
    else:
        df = pd.read_csv(tsv_file, sep="\t", na_values="NA", low_memory=False)

    xlsx_file = os.path.splitext(tsv_file)[0] + ".xlsx"

    # Phase 1: Fast write with xlsxwriter
    with pd.ExcelWriter(xlsx_file, engine='xlsxwriter') as writer:
        df.to_excel(writer, sheet_name='Results', index=False)
        # Note: finalize_excel_file() will add hyperlinks/formatting via openpyxl

    return xlsx_file
```

### Current GT Parsing (Baseline)
```python
# Source: variantcentrifuge/converter.py (lines 286, 290-299)
# Pattern compiled INSIDE loop - inefficient
pattern = re.compile(r"([^()]+)\(([^)]+)\)")  # Compiled once per variant row

sample_entries = gt_value.split(";")
for entry in sample_entries:
    entry = entry.strip()
    m = pattern.match(entry)
    if m:
        sample_id = m.group(1).strip()
        genotype = m.group(2).strip()
        # Process sample_id and genotype
```

### Optimized GT Parsing (Pre-Parse + Cache)
```python
# NEW: In dataframe_optimizer.py
# Compile pattern ONCE at module level
GT_PATTERN = re.compile(r"([^()]+)\(([^)]+)\)")

def parse_gt_entry(gt_value: str) -> list[dict]:
    """Parse GT column entry into structured data.

    Args:
        gt_value: GT string like "Sample1(0/1);Sample2(1/1)"

    Returns:
        List of dicts: [{"sample": "Sample1", "gt": "0/1", "is_variant": True}, ...]
    """
    if pd.isna(gt_value) or not gt_value:
        return []

    parsed = []
    for entry in str(gt_value).split(";"):
        entry = entry.strip()
        if not entry:
            continue

        m = GT_PATTERN.match(entry)
        if m:
            sample_id = m.group(1).strip()
            genotype = m.group(2).strip()

            # Determine if variant (not 0/0 or ./.)
            is_variant = genotype not in ["0/0", "./."]

            parsed.append({
                "sample": sample_id,
                "gt": genotype,
                "is_variant": is_variant
            })

    return parsed

def add_gt_cache(df: pd.DataFrame) -> pd.DataFrame:
    """Add pre-parsed GT cache column to DataFrame."""
    if 'GT' not in df.columns:
        return df

    logger.info("Pre-parsing GT column for downstream stages")
    df['_GT_PARSED'] = df['GT'].apply(parse_gt_entry)

    return df
```

### Using GT Cache in Downstream Stages
```python
# In converter.py finalize_excel_file()
# CURRENT: Re-parses GT for every row
pattern = re.compile(r"([^()]+)\(([^)]+)\)")
for entry in gt_value.split(";"):
    m = pattern.match(entry)
    # ... process

# OPTIMIZED: Use pre-parsed cache
if '_GT_PARSED' in df.columns:
    # Get pre-parsed data from DataFrame row
    parsed_gt = row._GT_PARSED  # Access via namedtuple attribute
    for entry in parsed_gt:
        sample_id = entry["sample"]
        genotype = entry["gt"]
        is_variant = entry["is_variant"]
        # ... process
```

### Batch Hyperlink Addition (openpyxl)
```python
# Source: variantcentrifuge/converter.py (lines 238-252) + openpyxl best practices
# CURRENT: Iterates all URL columns, all rows
for row_idx in range(2, ws.max_row + 1):
    for col_idx in url_columns:
        cell = ws.cell(row=row_idx, column=col_idx)
        if cell.value and str(cell.value).startswith("http"):
            cell.hyperlink = str(cell.value)
            cell.value = "Link"  # Display text
            cell.style = "Hyperlink"

# OPTIMIZED: Same approach but minimize cell access
# (openpyxl already caches cells; main optimization is using pre-parsed GT to reduce processing)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| openpyxl-only Excel generation | xlsxwriter + openpyxl two-pass | 2020+ | 2-5x faster for write-heavy workflows; industry standard split |
| Per-stage GT parsing | Pre-parse once at load time | This phase (Phase 10) | Eliminates redundant regex; faster downstream processing |
| Disk read for Excel | In-memory DataFrame pass-through | Phase 8 (implemented), Phase 10 (full utilization) | Eliminates redundant I/O for small-medium datasets |
| Manual dtype specification | Auto-detect categorical columns | Phase 8 | 50-75% memory reduction; faster load via PyArrow |

**Deprecated/outdated:**
- **xlwt/xlrd:** Old libraries for .xls (Excel 97-2003) format. Deprecated in favor of openpyxl/xlsxwriter for .xlsx (modern Excel 2007+)
- **Row-by-row openpyxl writes:** Slow for large datasets. Replaced by pandas bulk write with xlsxwriter, then openpyxl for finalization
- **Global regex pattern in loops:** Pre-2015 pattern. Now compile once, reuse.

## Open Questions

Things that couldn't be fully resolved:

1. **GT Cache Storage Format**
   - What we know: Dict column (`list[dict]`) works but may not be most efficient; NumPy structured arrays possible but more complex
   - What's unclear: Actual performance difference between dict column vs structured array for GT cache; memory footprint comparison
   - Recommendation: Start with dict column for simplicity (easier to work with, JSON-serializable), profile, optimize to structured array only if needed

2. **Hyperlink Batch Performance**
   - What we know: Current code iterates all rows/columns sequentially; openpyxl doesn't have bulk hyperlink API
   - What's unclear: Whether batching hyperlink writes (e.g., collecting all hyperlinks, applying in single pass) provides meaningful speedup
   - Recommendation: Keep current row-by-row approach for hyperlinks (most time is in initial write, not finalization); profile to confirm

3. **GT Cache Column Cleanup**
   - What we know: `_GT_PARSED` column should not appear in final output; TSVOutputStage restores original column names
   - What's unclear: Whether to drop cache columns explicitly or rely on column name restoration filtering them out
   - Recommendation: Explicitly drop cache columns before output for clarity: `df = df.drop(columns=[c for c in df.columns if c.startswith('_')])`

## Sources

### Primary (HIGH confidence)
- [xlsxwriter documentation - Working with Pandas](https://xlsxwriter.readthedocs.io/working_with_pandas.html) - Official integration guide
- [pandas.ExcelWriter documentation](https://pandas.pydata.org/docs/reference/api/pandas.ExcelWriter.html) - Engine parameter usage
- [openpyxl documentation - Performance](https://openpyxl.readthedocs.io/en/3.1/performance.html) - Performance comparison with xlsxwriter
- variantcentrifuge/converter.py (lines 26-498) - Current Excel generation implementation
- variantcentrifuge/dataframe_optimizer.py (lines 232-343) - Current DataFrame loading with optimizations
- variantcentrifuge/stages/output_stages.py (lines 553-692) - ExcelReportStage implementation

### Secondary (MEDIUM confidence)
- [Openpyxl vs XlsxWriter comparison](https://hive.blog/python/@geekgirl/openpyxl-vs-xlsxwriter-the-ultimate-showdown-for-excel-automation) - Performance benchmarks
- [Real Python - openpyxl guide](https://realpython.com/openpyxl-excel-spreadsheets-python/) - Freeze panes, auto-filters, hyperlinks examples
- [Medium - Excel optimization case study](https://mass-software-solutions.medium.com/optimizing-excel-report-generation-from-openpyxl-to-xlsmerge-processing-52-columns-200k-rows-5b5a03ecbcd4) - 200K row optimization (9min → 3min with xlsxwriter)

### Tertiary (LOW confidence)
- WebSearch results on xlsxwriter best practices (2026) - General ecosystem patterns

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - xlsxwriter and openpyxl are industry standard, well-documented, actively maintained
- Architecture: HIGH - Two-pass pattern is established best practice; GT pre-parsing pattern is straightforward extension of Phase 8
- Pitfalls: MEDIUM - Based on codebase review and documentation; some edge cases may exist

**Research date:** 2026-02-14
**Valid until:** 2026-04-14 (60 days - stable ecosystem for Excel libraries)
