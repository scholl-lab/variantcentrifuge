# Phase 9: Inheritance Analysis Optimization - Research

**Researched:** 2026-02-14
**Domain:** Vectorization of complex genomic inheritance pattern analysis
**Confidence:** HIGH

## Summary

This research analyzes the current three-pass inheritance analysis architecture (pattern deduction, compound heterozygous detection, prioritization) to identify vectorization opportunities for achieving 10-100x speedup. The primary bottleneck is Pass 1 (df.apply with deduce_patterns_for_variant) which applies complex conditional logic per-variant via Python-level row iteration.

The standard approach for vectorizing complex conditional logic in pandas/NumPy genomic analysis is:
1. **Boolean mask operations** with bitwise operators (&, |, ~) on NumPy arrays
2. **np.select()** for multi-branch conditionals (replaces nested if/elif chains)
3. **Categorical dtype encoding** for genotypes (already completed in Phase 8)
4. **Gene-level groupby().apply()** with vectorized inner functions

**Primary recommendation:** Vectorize pass-by-pass with golden file validation between each pass. Start with Pass 1 using np.select() for pattern branches, followed by Pass 2 (compound het already has vectorized implementation), then Pass 3 (prioritization logic).

## Standard Stack

The established libraries/tools for vectorizing genomic analysis:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| NumPy | 1.24+ | Array operations, boolean masking, np.select() | Compiled C operations 10-100x faster than Python loops |
| pandas | 2.0+ | DataFrame operations, groupby with observed=True | Industry standard for tabular genomic data |
| PyArrow | 14.0+ | Categorical dtype backend | 50-75% memory reduction (already implemented Phase 8) |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pytest-benchmark | 4.0+ | Performance measurement | Verify 10-100x speedup claims |
| scikit-allel | 1.3+ | Mendelian error arrays | Reference for genomic vectorization patterns |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| np.select() | Nested np.where() | np.where() becomes unreadable with >3 branches |
| Boolean masks | Custom ufuncs | ufuncs add complexity, masks are more maintainable |
| pandas groupby | Manual chunking | groupby with observed=True is cleaner, well-tested |

**Installation:**
```bash
# All dependencies already present in project
# PyArrow backend enabled in Phase 8
```

## Architecture Patterns

### Current Three-Pass Structure
```
Pass 1: Per-Variant Pattern Deduction (86 lines analyzer.py)
├── df.apply(lambda row: deduce_patterns_for_variant())  ← BOTTLENECK
├── deduce_patterns_for_variant() calls:
│   ├── deduce_patterns_for_sample() per sample
│   │   ├── check_dominant_pattern()
│   │   ├── check_recessive_pattern()
│   │   ├── check_x_linked_patterns()
│   │   └── check_mitochondrial_pattern()
│   └── Returns list of pattern strings per variant
└── Result: df['_inheritance_patterns'] = list of patterns

Pass 2: Compound Heterozygous Analysis (128-134 lines analyzer.py)
├── Per-gene groupby with observed=True
├── analyze_gene_for_compound_het_vectorized() ← ALREADY VECTORIZED
│   ├── Encode genotypes as int8 NumPy arrays
│   ├── Boolean masks for heterozygous variants
│   ├── Parent-of-origin determination via array ops
│   └── Returns dict[variant_key][sample_id] = comp_het_info
└── Result: df['_comp_het_info'] applied via itertuples ← VECTORIZE THIS

Pass 3: Prioritization and Finalization (138-187 lines analyzer.py)
├── Per-variant itertuples loop ← VECTORIZE THIS
├── Collect all patterns (deduced + compound het)
├── calculate_segregation_score() per pattern
├── prioritize_patterns() using PATTERN_PRIORITY dict
├── create_inheritance_details() builds JSON
└── Result: df['Inheritance_Pattern'] and df['Inheritance_Details']
```

### Pattern 1: Vectorized Pattern Deduction
**What:** Replace df.apply(deduce_patterns_for_variant) with NumPy boolean mask operations
**When to use:** Pass 1 deduction of de novo, AD, AR, X-linked, mitochondrial patterns
**Example:**
```python
# Current (scalar): analyzer.py line 86-88
df["_inheritance_patterns"] = df.apply(
    lambda row: deduce_patterns_for_variant(row.to_dict(), pedigree_data, sample_list),
    axis=1
)

# Vectorized approach using boolean masks:
def vectorized_deduce_patterns(df, pedigree_data, sample_list):
    """
    Vectorized inheritance pattern deduction using NumPy boolean masks.

    Key insight: Most patterns are per-sample checks applied to ALL variants.
    Instead of looping over variants, create boolean masks for each pattern type.
    """
    # Initialize result array (list of patterns per variant)
    n_variants = len(df)
    patterns_per_variant = [[] for _ in range(n_variants)]

    # Create genotype encoding matrix (variants x samples)
    # -1=missing, 0=ref, 1=het, 2=hom_alt
    gt_matrix = encode_all_genotypes(df, sample_list)

    # For each sample, deduce patterns using vectorized operations
    for sample_id in sample_list:
        sample_idx = sample_list.index(sample_id)
        sample_gts = gt_matrix[:, sample_idx]  # All variants for this sample

        # Get parent genotypes if available
        father_id, mother_id = get_parents(sample_id, pedigree_data)
        if father_id and mother_id:
            father_idx = sample_list.index(father_id) if father_id in sample_list else None
            mother_idx = sample_list.index(mother_id) if mother_id in sample_list else None

            if father_idx is not None and mother_idx is not None:
                father_gts = gt_matrix[:, father_idx]
                mother_gts = gt_matrix[:, mother_idx]

                # De novo mask: child has variant, both parents ref
                de_novo_mask = (
                    (sample_gts > 0) &
                    (father_gts == 0) &
                    (mother_gts == 0)
                )

                # Recessive mask: child hom_alt, both parents het
                recessive_mask = (
                    (sample_gts == 2) &
                    (father_gts == 1) &
                    (mother_gts == 1)
                )

                # Dominant mask: child het, at least one parent het
                dominant_mask = (
                    (sample_gts == 1) &
                    ((father_gts > 0) | (mother_gts > 0))
                )

                # Apply masks to patterns_per_variant
                for i in np.where(de_novo_mask)[0]:
                    patterns_per_variant[i].append("de_novo")
                for i in np.where(recessive_mask)[0]:
                    patterns_per_variant[i].append("autosomal_recessive")
                for i in np.where(dominant_mask)[0]:
                    patterns_per_variant[i].append("autosomal_dominant")

    # X-linked and mitochondrial patterns require chromosome filtering
    # Use np.select for chromosome-specific logic
    chrom_array = df['CHROM'].values
    x_mask = np.isin(chrom_array, ['X', 'CHRX', '23'])
    mt_mask = np.isin(chrom_array, ['MT', 'M', 'CHRM', 'CHRMT'])

    # Apply X-linked and MT patterns only to relevant variants
    if x_mask.any():
        vectorized_x_linked_patterns(df[x_mask], gt_matrix[x_mask],
                                     patterns_per_variant, pedigree_data, sample_list)
    if mt_mask.any():
        vectorized_mitochondrial_patterns(df[mt_mask], gt_matrix[mt_mask],
                                          patterns_per_variant, pedigree_data, sample_list)

    return patterns_per_variant
```

### Pattern 2: np.select() for Multi-Branch Conditionals
**What:** Use np.select() instead of nested if/elif chains for pattern classification
**When to use:** When deduce_patterns_for_sample has 5+ conditional branches (de novo, AD, AR, X-linked, MT)
**Example:**
```python
# Instead of nested if/elif in deducer.py:
def classify_pattern_vectorized(sample_gts, father_gts, mother_gts, affected_mask, chrom):
    """
    Use np.select to classify inheritance patterns for all variants at once.

    Replaces 5-level if/elif chain with single np.select call.
    """
    # Define condition list (evaluated in order, first match wins)
    conditions = [
        # De novo: child variant, both parents ref
        (sample_gts > 0) & (father_gts == 0) & (mother_gts == 0),

        # Recessive: child hom_alt, both parents het, affected
        (sample_gts == 2) & (father_gts == 1) & (mother_gts == 1) & affected_mask,

        # Dominant: child het, at least one parent het, affected
        (sample_gts == 1) & ((father_gts > 0) | (mother_gts > 0)) & affected_mask,

        # X-linked (only on X chromosome): handled separately
        (chrom == 'X') & (sample_gts > 0),

        # Mitochondrial (only on MT): maternal transmission
        (chrom == 'MT') & (sample_gts > 0) & (mother_gts > 0),
    ]

    # Define choice list (pattern names)
    choices = [
        'de_novo',
        'autosomal_recessive',
        'autosomal_dominant',
        'x_linked',
        'mitochondrial',
    ]

    # Select pattern based on conditions
    # Default to 'unknown' if no conditions match
    patterns = np.select(conditions, choices, default='unknown')

    return patterns
```

### Pattern 3: Vectorized Compound Het Application
**What:** Replace itertuples loop (lines 128-134 analyzer.py) with vectorized column operations
**When to use:** Pass 2 result application
**Example:**
```python
# Current (scalar): analyzer.py lines 128-134
for row in df.itertuples(index=True):
    variant_key = create_variant_key(row)
    gene = getattr(row, "GENE", "")
    if gene in comp_het_results_by_gene and variant_key in comp_het_results_by_gene[gene]:
        df.at[row.Index, "_comp_het_info"] = comp_het_results_by_gene[gene][variant_key]

# Vectorized approach:
def apply_comp_het_results_vectorized(df, comp_het_results_by_gene):
    """
    Apply compound het results using vectorized operations.

    Key insight: Create variant keys for all variants at once,
    then use dict lookup on entire arrays.
    """
    # Create variant keys vectorized
    variant_keys = (
        df['CHROM'].astype(str) + ':' +
        df['POS'].astype(str) + ':' +
        df['REF'].astype(str) + '>' +
        df['ALT'].astype(str)
    )

    # Create gene array
    genes = df['GENE'].values

    # Initialize _comp_het_info column
    df['_comp_het_info'] = None

    # For each gene with comp het results, apply in batch
    for gene, gene_results in comp_het_results_by_gene.items():
        # Boolean mask for variants in this gene
        gene_mask = genes == gene

        # Get variant keys for this gene
        gene_variant_keys = variant_keys[gene_mask]

        # Apply results where variant_key exists in gene_results
        for variant_key, comp_het_info in gene_results.items():
            # Find matching variants
            match_mask = gene_variant_keys == variant_key
            if match_mask.any():
                # Get original indices
                original_indices = df.index[gene_mask][match_mask]
                df.loc[original_indices, '_comp_het_info'] = [comp_het_info] * len(original_indices)

    return df
```

### Anti-Patterns to Avoid
- **Scalar genotype parsing in loops:** Always encode genotypes to NumPy arrays first, then operate on arrays
- **Row-wise pattern accumulation:** Build pattern lists via boolean masks, not append in loops
- **Nested np.where() chains:** Use np.select() for >2 conditions
- **Missing genotype special cases in vectorized code:** Encode missing as -1, handle via mask operations

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Genotype encoding | String parsing per variant | encode_genotypes() with GENOTYPE_ENCODING dict | comp_het_vectorized.py already has int8 encoding, 60x memory reduction |
| Pattern priority sorting | Custom sort key function | PATTERN_PRIORITY dict from prioritizer.py | Well-tested clinical priority order, matches variant-linker |
| Segregation p-values | Custom Fisher's exact | scipy.stats.fisher_exact | segregation_checker.py already uses this |
| Compound het phase determination | Manual parent-of-origin | find_potential_partners_vectorized() | Handles trans/cis/ambiguous correctly |
| Boolean mask combination | Manual array indexing | Bitwise operators (&, \|, ~) | Vectorized in NumPy, avoids index errors |

**Key insight:** comp_het_vectorized.py (lines 1-449) provides the reference pattern for vectorization: genotype encoding, boolean masks, int8 arrays, parent-of-origin via NumPy ops. Pass 1 should follow the same architecture.

## Common Pitfalls

### Pitfall 1: Logical operators (and/or) instead of bitwise (&/|)
**What goes wrong:** `and` evaluates entire array as single boolean, `&` operates element-wise
**Why it happens:** Python intuition doesn't transfer to NumPy arrays
**How to avoid:** Always use bitwise operators (&, |, ~) for boolean array operations
**Warning signs:** ValueError: "The truth value of an array with more than one element is ambiguous"
**Example:**
```python
# WRONG - uses logical operator
mask = (sample_gts > 0) and (father_gts == 0)  # ValueError!

# CORRECT - uses bitwise operator
mask = (sample_gts > 0) & (father_gts == 0)  # Element-wise AND
```

### Pitfall 2: Mixed categorical and object dtypes in boolean operations
**What goes wrong:** Categorical columns don't broadcast correctly with object columns
**Why it happens:** Phase 8 introduced categorical dtypes for low-cardinality columns
**How to avoid:** Convert to .values or ensure consistent dtype before boolean ops
**Warning signs:** TypeError: "Cannot perform operation with categorical dtype"
**Example:**
```python
# WRONG - mixing categorical CHROM with string literal
x_mask = df['CHROM'] == 'X'  # May fail if CHROM is categorical

# CORRECT - use .isin() for categorical-safe comparison
x_mask = df['CHROM'].isin(['X', 'CHRX', '23'])
```

### Pitfall 3: Applying vectorized function per-row via df.apply()
**What goes wrong:** df.apply() calls Python interpreter per row, negating vectorization benefits
**Why it happens:** Habit of using df.apply() as default pattern
**How to avoid:** Call vectorized function ONCE on entire DataFrame, not per-row
**Warning signs:** Vectorized function is fast but overall time is slow
**Example:**
```python
# WRONG - calls vectorized function per row
df['patterns'] = df.apply(lambda row: vectorized_deduce(row), axis=1)

# CORRECT - call once on entire DataFrame
df['patterns'] = vectorized_deduce(df, pedigree_data, sample_list)
```

### Pitfall 4: Genotype encoding inside hot loop
**What goes wrong:** Parsing "0/1" strings to integers millions of times
**Why it happens:** Not recognizing genotype encoding as preprocessing step
**How to avoid:** Encode ALL genotypes to int8 array ONCE before pattern deduction
**Warning signs:** String split/parse operations show up in profiler
**Example:**
```python
# WRONG - parsing per variant per sample
for variant in variants:
    for sample in samples:
        gt = variant[sample]
        if gt == "0/1": ...  # String comparison in hot loop

# CORRECT - encode once, operate on integers
gt_matrix = encode_all_genotypes(df, sample_list)  # Once
for sample_idx in range(len(sample_list)):
    sample_gts = gt_matrix[:, sample_idx]
    het_mask = sample_gts == 1  # Integer comparison
```

### Pitfall 5: Creating pattern lists per-variant instead of using boolean indexing
**What goes wrong:** Appending to lists in loops defeats vectorization
**Why it happens:** Scalar thinking (build result per item) vs vector thinking (filter arrays)
**How to avoid:** Create boolean mask for each pattern, use np.where() to get indices
**Warning signs:** List comprehensions or append in loops
**Example:**
```python
# WRONG - building lists per variant
patterns_per_variant = []
for i in range(len(df)):
    patterns = []
    if de_novo_condition(i):
        patterns.append('de_novo')
    if recessive_condition(i):
        patterns.append('autosomal_recessive')
    patterns_per_variant.append(patterns)

# CORRECT - boolean masks + indexing
de_novo_mask = vectorized_de_novo_check(df, pedigree_data)
recessive_mask = vectorized_recessive_check(df, pedigree_data)

patterns_per_variant = [[] for _ in range(len(df))]
for i in np.where(de_novo_mask)[0]:
    patterns_per_variant[i].append('de_novo')
for i in np.where(recessive_mask)[0]:
    patterns_per_variant[i].append('autosomal_recessive')
```

### Pitfall 6: Ignoring groupby observed=True with categorical dtypes
**What goes wrong:** groupby on categorical columns iterates over ALL categories (including unused), 3500x slowdown
**Why it happens:** Phase 8 added categorical dtypes, old groupby code doesn't specify observed=True
**How to avoid:** ALWAYS use .groupby(col, observed=True) with categorical columns
**Warning signs:** Groupby is mysteriously slow after Phase 8 PyArrow loading
**Example:**
```python
# WRONG - iterates over all possible GENE values (thousands)
for gene, gene_df in df.groupby('GENE'):  # Slow with categorical

# CORRECT - only iterates over observed genes in data
for gene, gene_df in df.groupby('GENE', observed=True):  # Fast
```

## Code Examples

Verified patterns from comp_het_vectorized.py and current implementation:

### Example 1: Genotype Encoding (comp_het_vectorized.py lines 36-66)
```python
# Source: comp_het_vectorized.py encode_genotypes()
GENOTYPE_ENCODING = {
    "./.": -1, "0/0": 0, "0/1": 1, "1/0": 1,
    "1/1": 2, "0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2,
}

def encode_genotypes(genotype_series: pd.Series) -> np.ndarray:
    """Encode genotypes as int8 for vectorized operations."""
    gt_strings = genotype_series.fillna("./.").astype(str)
    encoded = np.zeros(len(gt_strings), dtype=np.int8)

    for gt_str, code in GENOTYPE_ENCODING.items():
        mask = gt_strings == gt_str
        encoded[mask] = code

    # Handle unrecognized as missing
    unrecognized = ~gt_strings.isin(GENOTYPE_ENCODING.keys())
    if unrecognized.any():
        encoded[unrecognized] = -1

    return encoded
```

### Example 2: Boolean Mask Pattern Detection (comp_het_vectorized.py lines 129-134)
```python
# Source: comp_het_vectorized.py analyze_gene_for_compound_het_vectorized()
sample_genotypes = genotype_matrix.get(sample_id, np.array([]))

# Find heterozygous variants (encoded as 1)
het_mask = sample_genotypes == 1
het_indices = np.where(het_mask)[0]

# Skip if fewer than 2 heterozygous variants
if len(het_indices) < 2:
    continue
```

### Example 3: Parent-of-Origin Determination (comp_het_vectorized.py lines 252-265)
```python
# Source: comp_het_vectorized.py find_potential_partners_vectorized()
# Extract parent genotypes for heterozygous positions
father_het_gts = father_genotypes[het_indices]
mother_het_gts = mother_genotypes[het_indices]

# Check if each variant is present in each parent
father_has_var = father_het_gts > 0
mother_has_var = mother_het_gts > 0

# Determine origin of each variant
from_father_only = father_has_var & ~mother_has_var
from_mother_only = ~father_has_var & mother_has_var
from_both = father_has_var & mother_has_var
from_neither = ~father_has_var & ~mother_has_var
```

### Example 4: Chromosome-Specific Filtering (should be added to Pass 1)
```python
# Pattern for X-linked and mitochondrial handling
# Use categorical-safe comparison
chrom_array = df['CHROM'].values  # Convert to NumPy array
x_mask = np.isin(chrom_array, ['X', 'CHRX', '23'])
mt_mask = np.isin(chrom_array, ['MT', 'M', 'CHRM', 'CHRMT'])

# Apply pattern detection only to relevant chromosomes
if x_mask.any():
    x_variants = df[x_mask]
    x_patterns = deduce_x_linked_vectorized(x_variants, gt_matrix[x_mask], ...)
```

### Example 5: Creating Variant Keys Vectorized (for Pass 2 application)
```python
# Source: Adapted from create_variant_key pattern
# Instead of per-row string concatenation:
variant_keys = (
    df['CHROM'].astype(str) + ':' +
    df['POS'].astype(str) + ':' +
    df['REF'].astype(str) + '>' +
    df['ALT'].astype(str)
)
# Returns pd.Series of variant keys, vectorized string ops
```

### Example 6: np.select() for Multi-Condition Assignment
```python
# Pattern for replacing if/elif chains with np.select()
# Source: NumPy documentation + best practices research

# Define conditions (evaluated in order)
conditions = [
    (sample_gts > 0) & (father_gts == 0) & (mother_gts == 0),  # De novo
    (sample_gts == 2) & (father_gts == 1) & (mother_gts == 1),  # Recessive
    (sample_gts == 1) & ((father_gts > 0) | (mother_gts > 0)),  # Dominant
]

# Define choices (corresponding values)
choices = ['de_novo', 'autosomal_recessive', 'autosomal_dominant']

# Apply selection (first matching condition wins)
patterns = np.select(conditions, choices, default='unknown')
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| df.apply(axis=1) for all operations | itertuples for hot paths | Phase 8 (2026-02-14) | 10-14x iteration speedup |
| String genotypes in DataFrames | int8 encoded genotype matrices | comp_het_vectorized.py (existing) | 60x memory reduction |
| iterrows for compound het | NumPy boolean masks | comp_het_vectorized.py (existing) | 10-50x speedup |
| Regular dtypes | PyArrow categorical | Phase 8 (2026-02-14) | 50-75% memory reduction |
| All groupby calls | groupby(observed=True) | Phase 7 (2026-02-14) | Prevents 3500x slowdown with categorical |
| Manual GC | gc.collect() after stages | Phase 7 (2026-02-14) | 20-58% collateral speedup on inheritance |

**Deprecated/outdated:**
- **df.apply(axis=1) in hot paths:** Phase 8 established itertuples as standard, Phase 9 goes further to full vectorization
- **Nested np.where() chains:** np.select() is more readable and maintainable for >2 conditions (best practice as of 2026)
- **Row-wise genotype parsing:** Genotype encoding as preprocessing is now standard pattern (comp_het_vectorized.py proves viability)

## Open Questions

Things that couldn't be fully resolved:

1. **Edge case handling strategy for vectorized Pass 1**
   - What we know: Missing genotypes encoded as -1, handled via boolean masks
   - What's unclear: Whether unusual ploidy (e.g., "1/1/1") should be vectorized or fall back to scalar
   - Recommendation: Encode unusual genotypes as -1 (missing), log warning, handle via mask operations. This maintains vectorization while preserving correctness.

2. **Optimal chunking strategy for large pedigrees (>100 samples)**
   - What we know: comp_het_vectorized.py processes all samples for a gene at once
   - What's unclear: At what sample count does memory pressure require chunking?
   - Recommendation: Profile with synthetic 100-sample pedigree. If memory acceptable, no chunking needed (simplicity wins).

3. **Segregation score calculation vectorization (Pass 3)**
   - What we know: calculate_segregation_score() calls per-variant, uses scalar conditionals
   - What's unclear: Whether segregation can be vectorized or if it's inherently per-variant
   - Recommendation: Investigate during Pass 3 planning. May be acceptable to leave scalar if not on critical path.

4. **Pattern list accumulation vs sparse matrix representation**
   - What we know: patterns_per_variant is list[list[str]], can have multiple patterns per variant
   - What's unclear: Whether sparse matrix (variant x pattern boolean matrix) would be more efficient
   - Recommendation: Start with list accumulation (matches current output structure). Profile and optimize if bottleneck.

## Sources

### Primary (HIGH confidence)
- **Codebase analysis:** variantcentrifuge/inheritance/*.py (analyzer.py, deducer.py, comp_het.py, comp_het_vectorized.py, segregation_checker.py, prioritizer.py, parallel_analyzer.py)
- **Benchmark infrastructure:** tests/performance/benchmark_inheritance.py, benchmark_comp_het.py (existing baseline measurements)
- **Test coverage:** tests/test_inheritance/*.py (4665 lines, comprehensive fixtures and test cases)
- **Phase 7-8 optimizations:** .planning/STATE.md (20-58% inheritance speedup from gc.collect + observed=True, PyArrow categorical dtypes)

### Secondary (MEDIUM confidence)
- [Python Loop Replacement: Handling Conditional Logic (PyTorch & NumPy)](https://medium.com/@zmadscientist/python-loop-replacement-pytorch-numpy-optimizations-d2e64ed2f355) — Boolean mask patterns
- [Comparisons, Masks, and Boolean Logic | Python Data Science Handbook](https://jakevdp.github.io/PythonDataScienceHandbook/02.06-boolean-arrays-and-masks.html) — Bitwise vs logical operators
- [Vectorizing Conditional Logic in Pandas](https://medium.com/@connect.hashblock/vectorizing-conditional-logic-in-pandas-the-trick-that-saved-me-hours-821e13a9fb9d) — np.select() usage
- [np.select vs. np.where: How to Choose the Right Function](https://satnamsingh99.medium.com/np-select-vs-np-where-how-to-choose-the-right-function-for-your-numpy-arrays-a0093a96b1d9) — Multi-condition selection
- [Mendelian transmission](http://alimanfoo.github.io/2017/02/14/mendelian-transmission.html) — Genomic vectorization with scikit-allel, mendel_errors arrays
- [How to Speed Up Python Pandas by over 300x](https://www.exxactcorp.com/blog/Deep-Learning/how-to-speed-up-python-pandas-by-over-300x) — Vectorization performance benefits
- [Tutorial: 700x speed improvement on Pandas with Vectorization](https://www.linkedin.com/pulse/tutorial-basic-vectorization-pandas-iterrows-apply-duc-lai-trung-minh-75d4c) — Apply vs vectorization

### Tertiary (LOW confidence)
- General pandas vectorization best practices articles (multiple sources) — Conceptual guidance, not specific to genomics

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — NumPy, pandas, PyArrow already in use, comp_het_vectorized.py proves viability
- Architecture: HIGH — Three-pass structure directly analyzed from analyzer.py source code
- Pitfalls: HIGH — Based on comp_het_vectorized.py implementation experience, Phase 8 categorical dtype issues
- Vectorization patterns: MEDIUM — np.select() and boolean masks are best practices, but genomic-specific application requires validation
- Edge case handling: MEDIUM — Missing genotypes well-understood, unusual ploidy needs empirical testing

**Research date:** 2026-02-14
**Valid until:** 30 days (stable domain, but NumPy 2.0 may introduce new optimization patterns)
