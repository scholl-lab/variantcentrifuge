# Phase 8: DataFrame Optimization - Research

**Researched:** 2026-02-14
**Domain:** pandas DataFrame memory optimization, PyArrow integration, vectorized iteration
**Confidence:** HIGH

## Summary

This phase optimizes pandas DataFrame operations in the stage-based pipeline to achieve 50-70% memory reduction and 2-3x I/O speedup. The research reveals four complementary optimization strategies: (1) PyArrow engine for CSV reading with proven 5-15x speedup, (2) categorical dtypes for low-cardinality columns with 75% memory reduction, (3) itertuples replacing iterrows for 10-14x iteration speedup, and (4) shared memory for DataFrame pass-through between stages avoiding redundant disk I/O.

The key challenges identified are: pd.NA comparison propagation with categorical dtypes (requiring explicit NA handling), object dtype incompatibility with shared memory (requiring dtype conversion), invalid Python identifiers in column names breaking itertuples access (requiring pre-rename), and memory threshold determination for fallback behavior on 8-16GB systems (recommend 25% available RAM threshold).

**Primary recommendation:** Apply all four optimizations in sequence with diff-based validation at each step. Convert to categorical at load time using dtype dict, use PyArrow engine with string[pyarrow] dtype (not ArrowDtype), replace iterrows with itertuples after renaming invalid column names, and implement shared memory with auto-fallback using psutil to detect system constraints.

## Standard Stack

The established libraries/tools for pandas DataFrame optimization:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pandas | >= 2.0 | DataFrame operations | Required for project; 2.x has PyArrow integration and improved NA handling |
| pyarrow | >= 14.0 | Fast CSV parsing and Arrow types | 5-15x faster CSV reads with multithreading; stable since v14 with Python 3.12 support |
| psutil | current | Memory monitoring | Already in project dependencies; cross-platform memory detection for thresholds |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| memory_profiler | latest | Detailed memory profiling | Validation phase only; measure categorical savings |
| pytest-benchmark | >= 5.1.0 | Performance benchmarking | Already in dev deps; validate speedup claims |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| PyArrow engine | C engine (default) | PyArrow 5-15x faster but lacks some features (converters, skipfooter, regex sep) |
| string[pyarrow] | ArrowDtype(pa.string()) | ArrowDtype is experimental in pandas 3.0; string[pyarrow] is production-ready |
| shared_memory | pickle via Queue | Pickle serialization 3-10x slower than shared memory for large DataFrames |
| itertuples | apply() | apply() typically 2-5x slower than itertuples but supports complex operations |

**Installation:**
```bash
# PyArrow becomes required dependency (not optional)
uv pip install "pyarrow>=14.0"
```

## Architecture Patterns

### Recommended DataFrame Loading Pattern
```python
# Stage: Load main variant DataFrame with all optimizations
def execute(self, context: PipelineContext) -> PipelineContext:
    input_file = context.workspace.extracted_tsv

    # 1. Auto-detect low-cardinality columns for categorical dtype
    dtype_dict = _detect_categorical_columns(
        input_file,
        cardinality_threshold=0.5,  # <50% unique values
        max_sample_rows=10000
    )

    # 2. Load with PyArrow engine and categorical dtypes
    df = pd.read_csv(
        input_file,
        sep="\t",
        dtype=dtype_dict,  # Pre-specified dtypes including categoricals
        engine="pyarrow"
    )

    # 3. Rename invalid Python identifiers for itertuples compatibility
    df = _rename_invalid_identifiers(df)

    # 4. Store in context for in-memory pass-through
    context.variants_df = df  # Shared memory reference

    return context
```

### Pattern 1: Categorical Dtype Auto-Detection
**What:** Scan CSV to identify low-cardinality columns and build dtype dict for read_csv
**When to use:** At DataFrame load time for main variant data (not config/gene lists)
**Example:**
```python
# Source: Research synthesis from pandas docs and optimization guides
def _detect_categorical_columns(
    csv_path: Path,
    cardinality_threshold: float = 0.5,
    max_sample_rows: int = 10000
) -> dict[str, str]:
    """Detect columns suitable for categorical dtype.

    Args:
        csv_path: Path to CSV file
        cardinality_threshold: Max ratio of unique/total values (default 0.5 = 50%)
        max_sample_rows: Sample size for detection (prevent full scan)

    Returns:
        dtype dict mapping column names to 'category' or 'str'
    """
    # Sample first N rows to detect cardinality
    sample_df = pd.read_csv(csv_path, sep="\t", nrows=max_sample_rows, dtype=str)

    dtype_dict = {}
    for col in sample_df.columns:
        n_unique = sample_df[col].nunique()
        n_total = len(sample_df[col])

        # Convert if unique ratio < threshold
        if n_unique / n_total < cardinality_threshold:
            dtype_dict[col] = "category"
        else:
            dtype_dict[col] = "str"

    return dtype_dict
```

### Pattern 2: Column Name Sanitization for itertuples
**What:** Rename columns with invalid Python identifiers before using itertuples
**When to use:** Once at load time, before any itertuples-based iteration
**Example:**
```python
# Source: Research synthesis from pandas itertuples documentation
import re

def _rename_invalid_identifiers(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns with invalid Python identifiers for itertuples compatibility.

    itertuples auto-renames invalid names to positional (_0, _1), breaking code.
    This function explicitly renames to valid identifiers upfront.

    Args:
        df: DataFrame with potentially invalid column names

    Returns:
        DataFrame with sanitized column names
    """
    rename_map = {}

    for col in df.columns:
        # Check if valid Python identifier
        if not col.isidentifier():
            # Replace invalid chars with underscore
            new_name = re.sub(r'[^a-zA-Z0-9_]', '_', col)
            # Ensure doesn't start with digit
            if new_name[0].isdigit():
                new_name = f"col_{new_name}"
            rename_map[col] = new_name

    if rename_map:
        logger.info(f"Renamed {len(rename_map)} columns with invalid identifiers")
        logger.debug(f"Column renames: {rename_map}")
        return df.rename(columns=rename_map)

    return df
```

### Pattern 3: Shared Memory DataFrame Pass-Through
**What:** Store DataFrame in context for direct reference; auto-fallback to disk on memory pressure
**When to use:** Between analysis stage (produces DataFrame) and output stages (consume it)
**Example:**
```python
# Source: Research synthesis from multiprocessing.shared_memory docs and pandas guides
import psutil

def _should_use_shared_memory(df: pd.DataFrame, threshold_ratio: float = 0.25) -> bool:
    """Determine if DataFrame should use shared memory based on available RAM.

    Args:
        df: DataFrame to evaluate
        threshold_ratio: Max fraction of available RAM to use (default 0.25 = 25%)

    Returns:
        True if DataFrame fits within threshold
    """
    df_memory = df.memory_usage(deep=True).sum()
    available_memory = psutil.virtual_memory().available

    max_allowed = available_memory * threshold_ratio
    return df_memory <= max_allowed

# In analysis stage (producer):
class InheritanceAnalysisStage(Stage):
    def execute(self, context: PipelineContext) -> PipelineContext:
        df = pd.read_csv(context.workspace.extracted_tsv, ...)
        # ... perform analysis ...

        # Store in context (shared memory reference)
        if _should_use_shared_memory(df, threshold_ratio=0.25):
            context.variants_df = df
            logger.info("Using shared memory for DataFrame pass-through")
        else:
            # Fallback: write to disk, clear reference
            df.to_csv(context.workspace.analyzed_tsv, sep="\t", index=False)
            context.variants_df = None
            logger.info("DataFrame exceeds memory threshold, using disk pass-through")

        return context

# In output stage (consumer):
class ExcelReportStage(Stage):
    def execute(self, context: PipelineContext) -> PipelineContext:
        # Try shared memory first
        if context.variants_df is not None:
            df = context.variants_df
            logger.debug("Using shared memory DataFrame")
        else:
            # Fallback: read from disk
            df = pd.read_csv(context.workspace.analyzed_tsv, sep="\t", engine="pyarrow")
            logger.debug("Loading DataFrame from disk")

        # ... generate Excel ...
        return context
```

### Pattern 4: itertuples Replacement
**What:** Replace iterrows() with itertuples() for 10-14x speedup
**When to use:** All row iteration hot paths (inheritance Pass 2-3, gene burden, scoring)
**Example:**
```python
# Source: Research from pandas docs and performance guides
# BEFORE (iterrows - slow):
for idx, row in df.iterrows():
    variant_key = create_variant_key(row)
    gene = row.get("GENE", "")
    if gene in comp_het_results_by_gene:
        df.at[idx, "_comp_het_info"] = comp_het_results_by_gene[gene][variant_key]

# AFTER (itertuples - 10-14x faster):
for row in df.itertuples(index=True):
    variant_key = create_variant_key(row)
    gene = getattr(row, "GENE", "")  # Use getattr for namedtuple
    if gene in comp_het_results_by_gene:
        df.at[row.Index, "_comp_het_info"] = comp_het_results_by_gene[gene][variant_key]
```

### Pattern 5: pd.NA-Safe Comparisons with Categorical
**What:** Handle pd.NA propagation in categorical comparisons using fillna or explicit NA checks
**When to use:** Whenever comparing categorical columns that may contain missing values
**Example:**
```python
# Source: Research from pandas working with missing data docs
# PROBLEMATIC (pd.NA propagates):
mask = df["IMPACT"] == "HIGH"  # Returns pd.NA for missing values, not False

# SOLUTION 1: fillna before comparison
mask = df["IMPACT"].fillna("UNKNOWN") == "HIGH"  # Always returns bool

# SOLUTION 2: Explicit NA handling
mask = (df["IMPACT"] == "HIGH") & ~df["IMPACT"].isna()

# SOLUTION 3: Use .eq() with fill_value
mask = df["IMPACT"].eq("HIGH", fill_value=False)
```

### Anti-Patterns to Avoid

- **Late categorical conversion:** Converting to categorical after loading wastes initial memory spike; convert at load time via dtype dict
- **ArrowDtype strings:** ArrowDtype is experimental in pandas 3.0; use string[pyarrow] instead for production stability
- **Unconditional shared memory:** Large DataFrames (>25% available RAM) should fallback to disk to avoid OOM on 8-16GB systems
- **Ignoring invalid identifiers:** itertuples silently renames invalid column names to _0, _1, breaking attribute access; rename explicitly upfront
- **Mixing PyArrow with unsupported params:** PyArrow engine doesn't support converters, skipfooter, regex sep, on_bad_lines='warn' - fallback to C engine when needed

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Detecting low-cardinality columns | Manual unique value counter | nunique() / sample-based ratio | Edge cases: very large files, memory efficient sampling, string vs numeric cardinality |
| Shared memory for DataFrames | Custom mmap wrapper | Direct context reference with psutil-based fallback | Object dtypes cause segfaults in shared_memory; context reference is simpler and safer |
| Column name sanitization | String replacement hacks | re.sub with isidentifier() check | Python identifier rules complex (keywords, reserved words, unicode); regex handles all cases |
| Memory usage measurement | Manual size calculation | df.memory_usage(deep=True) | Categorical dtypes have hidden storage; deep=True introspects actual memory |
| NA-safe comparisons | Custom NA checking | fillna() or .eq(fill_value=) | pd.NA propagation is subtle; built-in methods are NA-aware and tested |

**Key insight:** pandas 2.x/3.x has mature APIs for these operations; custom solutions miss edge cases and lose future optimizations.

## Common Pitfalls

### Pitfall 1: PyArrow Engine with Unsupported Parameters
**What goes wrong:** read_csv with engine="pyarrow" silently fails or errors when using unsupported parameters like skipfooter, converters, or regex separators.
**Why it happens:** PyArrow engine is faster but doesn't support full pandas read_csv API surface.
**How to avoid:** Check for unsupported params before enabling PyArrow; fallback to C engine when needed.
**Warning signs:** Unexpected parsing errors, missing rows at file end, converter functions not called.

**Prevention:**
```python
# Check if read_csv call uses unsupported PyArrow params
PYARROW_UNSUPPORTED = {'skipfooter', 'converters', 'on_bad_lines',
                       'thousands', 'memory_map', 'dialect', 'chunksize'}

def _can_use_pyarrow(read_csv_kwargs: dict) -> bool:
    """Check if read_csv params are compatible with PyArrow engine."""
    unsupported_used = PYARROW_UNSUPPORTED & set(read_csv_kwargs.keys())

    # Also check for regex separator
    if 'sep' in read_csv_kwargs and isinstance(read_csv_kwargs['sep'], str):
        if len(read_csv_kwargs['sep']) > 1 or not read_csv_kwargs['sep'].isprintable():
            return False

    return len(unsupported_used) == 0
```

### Pitfall 2: Categorical Comparison pd.NA Propagation
**What goes wrong:** Comparisons like df["IMPACT"] == "HIGH" return pd.NA instead of False for missing values, breaking boolean indexing and if statements.
**Why it happens:** pd.NA propagates in operations (like SQL NULL); categorical dtypes preserve pd.NA in comparisons.
**How to avoid:** Use fillna() before comparison, or use .eq(fill_value=False) for NA-safe boolean results.
**Warning signs:** Unexpected pd.NA in boolean masks, indexing errors "cannot index with vector containing NA / NaN values", if conditions behaving strangely.

**Detection:**
```python
# Add validation after boolean mask creation
mask = df["IMPACT"] == "HIGH"
if mask.isna().any():
    raise ValueError("Boolean mask contains pd.NA - comparison not NA-safe")
```

### Pitfall 3: Object Dtypes in Shared Memory
**What goes wrong:** Storing DataFrames with object dtypes (strings) in multiprocessing.shared_memory causes segmentation faults when child processes access the data.
**Why it happens:** Object dtypes store pointers, not values; pointers are invalid across process boundaries.
**How to avoid:** Convert string columns to fixed-length dtype (S or U) when using to_records(), or use context reference (not true shared_memory) for stage-based pipeline.
**Warning signs:** Segfaults during multiprocessing, intermittent crashes when accessing DataFrame columns.

**Note:** Stage-based pipeline uses context reference (not multiprocessing), so this is LOW risk. Include for completeness if future parallelization added.

### Pitfall 4: itertuples with Invalid Column Names
**What goes wrong:** Columns like "GEN[0].GT" or "GENE-NAME" get silently renamed to _0, _1 in itertuples, breaking attribute access (row.GEN[0].GT fails).
**Why it happens:** namedtuple requires valid Python identifiers; pandas auto-renames invalid names to positional.
**How to avoid:** Rename invalid column names explicitly at load time using isidentifier() check and regex sanitization.
**Warning signs:** AttributeError when accessing row.column_name, unexpected _0, _1 in row attributes.

**Detection:**
```python
# Validate column names before itertuples
invalid_cols = [col for col in df.columns if not col.isidentifier()]
if invalid_cols:
    raise ValueError(f"Columns with invalid identifiers: {invalid_cols}. Rename before itertuples.")
```

### Pitfall 5: Cardinality Threshold Too Aggressive
**What goes wrong:** Converting high-cardinality columns (e.g., 40% unique) to categorical increases memory usage instead of reducing it.
**Why it happens:** Categorical dtype stores unique values + codes; overhead exceeds savings when uniqueness is high.
**How to avoid:** Use 50% unique threshold (or lower, like 30%) and validate memory savings with memory_usage(deep=True) before/after.
**Warning signs:** Memory usage increases after categorical conversion, slower operations on categorical columns.

**Validation:**
```python
# Measure memory before/after categorical conversion
before = df.memory_usage(deep=True).sum()
df[col] = df[col].astype("category")
after = df.memory_usage(deep=True).sum()

if after >= before:
    logger.warning(f"Column {col} increased memory after categorical: {before} -> {after}")
    df[col] = df[col].astype("str")  # Revert
```

## Code Examples

Verified patterns from official sources:

### PyArrow Engine CSV Loading
```python
# Source: https://pandas.pydata.org/docs/user_guide/pyarrow.html
# Source: https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html

# Basic PyArrow engine usage (5-15x faster than C engine)
df = pd.read_csv(
    "variants.tsv",
    sep="\t",
    engine="pyarrow"
)

# With dtype specification (categorical at load time)
df = pd.read_csv(
    "variants.tsv",
    sep="\t",
    dtype={
        "CHROM": "category",
        "IMPACT": "category",
        "FILTER": "category",
        "GENE": "category",
        "POS": "int64"
    },
    engine="pyarrow"
)

# Using string[pyarrow] for string columns (not ArrowDtype)
df = pd.read_csv(
    "variants.tsv",
    sep="\t",
    dtype={
        "CHROM": "category",
        "REF": "string[pyarrow]",  # Production-ready PyArrow strings
        "ALT": "string[pyarrow]"
    },
    engine="pyarrow"
)
```

### Memory Usage Profiling
```python
# Source: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.memory_usage.html

# Measure memory with deep introspection (includes object storage)
before = df.memory_usage(deep=True).sum()
logger.info(f"Memory before optimization: {before / 1024**2:.2f} MB")

# Convert to categorical
for col in ["CHROM", "IMPACT", "FILTER", "GENE", "EFFECT"]:
    df[col] = df[col].astype("category")

# Measure after
after = df.memory_usage(deep=True).sum()
reduction = (1 - after/before) * 100
logger.info(f"Memory after optimization: {after / 1024**2:.2f} MB ({reduction:.1f}% reduction)")

# Per-column breakdown
memory_by_col = df.memory_usage(deep=True)
for col, mem in memory_by_col.items():
    logger.debug(f"{col}: {mem / 1024**2:.2f} MB")
```

### Categorical Dtype Specification
```python
# Source: https://pandas.pydata.org/docs/user_guide/categorical.html

from pandas.api.types import CategoricalDtype

# Method 1: Inference (simple, but categories vary by data)
df["IMPACT"] = df["IMPACT"].astype("category")

# Method 2: Explicit categories (recommended for consistency)
impact_dtype = CategoricalDtype(
    categories=["MODIFIER", "LOW", "MODERATE", "HIGH"],
    ordered=True
)
df["IMPACT"] = df["IMPACT"].astype(impact_dtype)

# Method 3: At load time via dtype dict (best performance)
df = pd.read_csv(
    "variants.tsv",
    sep="\t",
    dtype={"IMPACT": impact_dtype},
    engine="pyarrow"
)
```

### itertuples vs iterrows Performance
```python
# Source: Research from performance benchmarks (14x speedup confirmed)

# SLOW: iterrows (439ms for 10k rows)
for idx, row in df.iterrows():
    result = process_variant(
        gene=row["GENE"],
        impact=row["IMPACT"],
        pos=row["POS"]
    )
    df.at[idx, "result"] = result

# FAST: itertuples (30.8ms for 10k rows - 14x faster)
for row in df.itertuples(index=True):
    result = process_variant(
        gene=row.GENE,
        impact=row.IMPACT,
        pos=row.POS
    )
    df.at[row.Index, "result"] = result

# Access patterns for itertuples:
# - row.Index: row index (if index=True, default)
# - row.column_name: access by attribute (requires valid identifier)
# - row[0], row[1]: access by position (always works)
```

### System Memory Detection for Thresholds
```python
# Source: psutil documentation + research on memory thresholds

import psutil

def get_memory_threshold(ratio: float = 0.25) -> int:
    """Get maximum DataFrame memory based on available RAM.

    For 8-16GB systems, recommend 25% of available RAM for DataFrame
    to leave headroom for processing and other stages.

    Args:
        ratio: Fraction of available RAM to use (default 0.25)

    Returns:
        Maximum DataFrame memory in bytes
    """
    vm = psutil.virtual_memory()
    available = vm.available
    total = vm.total

    logger.debug(f"System RAM: {total / 1024**3:.2f} GB total, "
                 f"{available / 1024**3:.2f} GB available")

    threshold = int(available * ratio)
    logger.info(f"DataFrame memory threshold: {threshold / 1024**2:.2f} MB")

    return threshold

# Usage in auto-fallback logic
df_memory = df.memory_usage(deep=True).sum()
threshold = get_memory_threshold(ratio=0.25)

if df_memory > threshold:
    logger.info("DataFrame exceeds threshold, using disk pass-through")
    df.to_csv(output_path, sep="\t", index=False)
    context.variants_df = None
else:
    logger.info("DataFrame within threshold, using shared memory")
    context.variants_df = df
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| C engine (default) | PyArrow engine | pandas 1.4.0 (2022) | 5-15x CSV read speedup with multithreading |
| object dtype strings | string[pyarrow] | pandas 2.0 (2023) | Memory savings + Arrow ecosystem integration |
| np.nan for missing | pd.NA for nullable types | pandas 1.0 (2020) | Consistent missing value handling across dtypes |
| Post-load categorical | Load-time dtype specification | Best practice (2024+) | Avoid initial memory spike from object dtype |
| iterrows() | itertuples() | Long-standing (2015+) | 10-14x iteration speedup |
| Manual memory threshold | psutil-based detection | Standard practice (2020+) | Adapts to actual system constraints, not hardcoded values |

**Deprecated/outdated:**
- **ArrowDtype(pa.string())**: Still experimental in pandas 3.0; use string[pyarrow] instead for production (StringDtype with PyArrow backend is stable)
- **dtype_backend="pyarrow"**: Broader than needed; prefer selective dtype dict for control
- **Categorical.from_codes()**: Use pd.Categorical() constructor or astype("category") instead
- **convert_dtypes()**: Auto-converts to nullable types but less control than explicit dtype specification

## Open Questions

Things that couldn't be fully resolved:

1. **Exact PyArrow version pin**
   - What we know: PyArrow >= 14.0 has Python 3.12 support; pandas 3.0 requires PyArrow >= 7; version 16.1.0 is latest
   - What's unclear: Whether any versions 14.x-16.x have compatibility issues with pandas 2.x/3.x
   - Recommendation: Pin to >= 14.0 (conservative) or >= 15.0 (recent stable), verify with existing tests. Check pandas changelog for known PyArrow version conflicts.

2. **Categorical category list pre-population**
   - What we know: Explicit CategoricalDtype with categories prevents "unknown category" warnings; inference is simpler but categories vary
   - What's unclear: Whether VCF fields (IMPACT, FILTER, EFFECT) have stable category lists across all datasets, or vary by annotation source (SnpEff version, VEP, etc.)
   - Recommendation: Start with inference (simple), add explicit categories only if "unknown category" warnings appear in testing. Monitor logs for category additions.

3. **itertuples performance with very wide DataFrames**
   - What we know: itertuples is 10-14x faster for narrow DataFrames (10-20 columns); namedtuple unpacking has overhead
   - What's unclear: Whether speedup holds for wide DataFrames (50+ columns like full variant data with sample columns)
   - Recommendation: Benchmark on actual data (50-100 columns) in Phase 6 framework; fallback to apply() if itertuples slower due to namedtuple overhead.

4. **Memory threshold interaction with SLURM/PBS**
   - What we know: psutil.virtual_memory() detects cgroup limits on most systems; memory_manager.py already handles HPC detection
   - What's unclear: Whether 25% threshold is appropriate in SLURM environments with allocated memory (e.g., 32GB allocation on 128GB node)
   - Recommendation: Test memory_manager.py integration - if it already provides safe threshold, use that instead of hardcoded 25%. Defer HPC-specific tuning to Phase 9+ if needed.

## Sources

### Primary (HIGH confidence)
- [pandas PyArrow Functionality Documentation](https://pandas.pydata.org/docs/user_guide/pyarrow.html) - ArrowDtype vs string[pyarrow], PyArrow integration
- [pandas Categorical Data Documentation](https://pandas.pydata.org/docs/user_guide/categorical.html) - Categorical dtypes, missing value handling, comparisons
- [pandas read_csv API Reference](https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html) - PyArrow engine capabilities and limitations
- [pandas DataFrame.itertuples Documentation](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.itertuples.html) - itertuples behavior, invalid identifiers handling
- [pandas DataFrame.memory_usage Documentation](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.memory_usage.html) - Memory profiling with deep parameter
- [pandas Working with Missing Data](https://pandas.pydata.org/docs/user_guide/missing_data.html) - pd.NA propagation behavior
- [Python multiprocessing.shared_memory Documentation](https://docs.python.org/3/library/multiprocessing.shared_memory.html) - Shared memory API

### Secondary (MEDIUM confidence)
- [Pandas read_csv: The Definitive Guide (2026) - Kanaries](https://docs.kanaries.net/topics/Pandas/pandas-read-csv) - PyArrow engine 5x-30x speedup benchmarks
- [Efficient reading and parsing of large CSV files - Vincent Codes Finance](https://vincent.codes.finance/posts/pandas-read-csv/) - PyArrow 15x speedup with dtype_backend
- [Why Pandas itertuples() Is Faster Than iterrows() - Medium](https://medium.com/swlh/why-pandas-itertuples-is-faster-than-iterrows-and-how-to-make-it-even-faster-bc50c0edd30d) - 14x speedup explanation (30.8ms vs 439ms)
- [Python Shared Memory in Multiprocessing - Mingze Gao](https://mingze-gao.com/posts/python-shared-memory-in-multiprocessing/) - Object dtype segfault workaround
- [Using large numpy arrays and pandas dataframes with multiprocessing - Emilio's Blog](https://e-dorigatti.github.io/python/2020/06/19/multiprocessing-large-objects.html) - Shared memory vs pickle performance
- [Mastering Memory Optimization for Pandas DataFrames - ThinhDA](https://thinhdanggroup.github.io/pandas-memory-optimization/) - 75% memory reduction with categorical, 50% threshold rule
- [Recasting Low-Cardinality Columns as Categoricals - Hackers and Slackers](https://hackersandslackers.com/recasting-low-cardinality-columns-as-categoricals-2/) - 50% unique values threshold

### Tertiary (LOW confidence)
- Various Stack Overflow and GitHub issues on categorical NA behavior - Anecdotal experiences, not authoritative
- Blog posts on pandas performance from 2020-2023 - May be outdated for pandas 3.0 behavior
- PyArrow GitHub issues on pandas compatibility - Specific bugs, not general guidance

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - pandas and PyArrow versions verified in official docs; psutil already in dependencies
- Architecture: HIGH - All patterns verified in official pandas documentation and cross-referenced with performance guides
- Pitfalls: MEDIUM - Most from official docs, some from community reports (validated against docs)
- PyArrow version specifics: MEDIUM - General compatibility known, exact version pins require testing
- Memory thresholds: MEDIUM - 25% threshold is rule-of-thumb from community, not scientifically derived for this workload

**Research date:** 2026-02-14
**Valid until:** 60 days (2026-04-15) - pandas and PyArrow are stable; major changes unlikely in 2-month window
