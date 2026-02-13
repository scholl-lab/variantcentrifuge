# Domain Pitfalls: Performance Optimization in Clinical Genomics Pipelines

**Domain:** Clinical genomics variant analysis pipeline
**Focus:** Adding performance optimizations to pandas-heavy bioinformatics system
**Researched:** 2026-02-14
**Context:** Optimizing variantcentrifuge (clinical tool, 1035 tests, Windows + Linux)

---

## Executive Summary

Optimizing a **clinical genomics pipeline** is fundamentally different from optimizing a web service or data analytics tool. The critical constraint is **bitwise output identity** — a single dropped variant, changed filter result, or reordered row can invalidate clinical interpretation. Performance optimizations that introduce subtle semantic changes are **catastrophic failures**, not acceptable tradeoffs.

This document catalogs pitfalls specific to the 11 planned optimizations for variantcentrifuge, with emphasis on:
1. **Clinical correctness risks** (output must be identical)
2. **Cross-platform pitfalls** (Windows + Linux compatibility)
3. **Integration with existing defensive architecture** (dtype=str is intentional)
4. **Subtle pandas behavior changes** (especially pandas 3.0 string dtype defaults)

---

## Critical Pitfalls

### Pitfall 1: Pandas 3.0 String Dtype Breaking Changes with PyArrow Engine

**What goes wrong:**

Starting January 21, 2026, pandas 3.0 changed the default string dtype from `object` to dedicated `string[pyarrow]` when PyArrow is installed. Code that checks `df.dtypes == 'object'` or assumes string columns are object dtype will break. More critically: **missing value semantics change** from `np.nan` (object dtype) to `pd.NA` (PyArrow strings), breaking comparisons.

```python
# Pre-pandas 3.0 with dtype=str
df = pd.read_csv("data.tsv", sep="\t", dtype=str)
df["GENE"].dtype  # → dtype('O')  (object)
df["GENE"].isna().sum()  # Uses np.nan

# Pandas 3.0 with engine="pyarrow"
df = pd.read_csv("data.tsv", sep="\t", engine="pyarrow")
df["GENE"].dtype  # → string[pyarrow]
df["GENE"].isna().sum()  # Uses pd.NA, different comparison semantics
```

**Why it happens:**

The codebase defensively uses `dtype=str` everywhere to avoid type coercion surprises. Adding `engine="pyarrow"` while keeping `dtype=str` works, but removing `dtype=str` to gain memory benefits triggers pandas 3.0's new string inference.

**Consequences:**

- String comparisons that worked with `==` may fail with PyArrow strings (`pd.NA == "some_value"` returns `pd.NA`, not `False`)
- `fillna("")` behavior changes (PyArrow strings preserve `pd.NA`)
- `.str` accessor methods return PyArrow-backed Series (faster but different dtype)
- Genotype string matching code like `if genotype == "0/1"` may fail if genotype is `pd.NA`

**Prevention:**

1. **Explicit dtype retention:**
   ```python
   # Keep dtype=str when using pyarrow engine
   df = pd.read_csv(path, sep="\t", dtype=str, engine="pyarrow")
   # NOT: df = pd.read_csv(path, sep="\t", engine="pyarrow")  # DANGEROUS
   ```

2. **Test with pandas 3.0+ in CI:**
   ```yaml
   strategy:
     matrix:
       pandas-version: ["2.2", "3.0"]
   ```

3. **Explicit missing value checks:**
   ```python
   # BEFORE: if value == "":
   # AFTER: if pd.isna(value) or value == "":
   ```

4. **Regression test all string comparisons:**
   ```python
   def test_genotype_matching_with_missing():
       df = pd.DataFrame({"GT": ["0/1", "1/1", pd.NA]})
       result = match_genotypes(df)
       assert result.tolist() == expected  # Explicit comparison
   ```

**Detection:**

- Tests fail with `TypeError: boolean value of NA is ambiguous`
- Different variant counts between pandas 2.x and 3.0
- `AssertionError` in output identity checks

**Severity:** **CRITICAL** — silently changes clinical output

**Affected optimizations:** #4 (PyArrow engine), #5 (categorical/PyArrow dtypes)

**Sources:**
- [pandas 3.0.0 whatsnew](https://pandas.pydata.org/docs/whatsnew/v3.0.0.html)
- [String dtype breaking changes](https://github.com/pandas-dev/pandas/issues/59328)
- [pandas 3.0 InfoQ coverage](https://www.infoq.com/news/2026/02/pandas-library/)

---

### Pitfall 2: Vectorization Changing Row Processing Order with Non-Deterministic Results

**What goes wrong:**

Replacing `df.apply(axis=1)` or `iterrows()` with vectorized operations can change the order of operations when those operations have side effects or depend on row processing order. For inheritance analysis with compound heterozygous detection, order matters when building gene-level state.

```python
# BEFORE: Sequential row processing
for idx, row in df.iterrows():
    variant_key = create_variant_key(row)  # Uses CHROM:POS:REF:ALT
    gene = row["GENE"]
    # Side effect: updates global lookup
    comp_het_results[gene][variant_key] = analyze_variant(row)

# AFTER: Vectorized (WRONG if order matters)
df["variant_key"] = df["CHROM"] + ":" + df["POS"] + ":" + df["REF"] + ":" + df["ALT"]
# This processes all rows simultaneously — no sequential state
```

**Why it happens:**

Pandas vectorized operations execute in parallel conceptually (though CPython GIL means single-threaded). Code written assuming "process row 1, then row 2" breaks when vectorized.

**Consequences:**

- Compound het detection may miss variant pairs (expects sequential gene analysis)
- Variant prioritization changes if tie-breaking uses implicit row order
- Non-deterministic output across runs if dict insertion order matters

**Prevention:**

1. **Audit for side effects before vectorizing:**
   ```python
   # RED FLAG: Mutation of external state
   for idx, row in df.iterrows():
       global_lookup[row["key"]] = row["value"]  # SIDE EFFECT

   # SAFE: Pure function returning values
   results = df.apply(lambda row: pure_function(row), axis=1)
   ```

2. **Preserve order-dependent logic explicitly:**
   ```python
   # If compound het needs pairs from same gene processed together:
   for gene, gene_df in df.groupby("GENE", sort=False):  # sort=False preserves input order
       # Process gene_df rows in order
   ```

3. **Test determinism:**
   ```python
   @pytest.mark.parametrize("run_number", range(5))
   def test_inheritance_analysis_deterministic(run_number):
       result = analyze_inheritance(df_fixture.copy())
       assert result.equals(expected_output)  # Must be identical every run
   ```

4. **Add checksums to output:**
   ```python
   # In final TSV/Excel, include row hash
   df["_row_checksum"] = df.apply(lambda r: hash(tuple(r)), axis=1)
   ```

**Detection:**

- Different output on repeated runs with same input
- Test failures with "Expected X variants, got Y"
- Compound het pairs missing in vectorized version

**Severity:** **CRITICAL** — breaks clinical correctness through non-determinism

**Affected optimizations:** #1 (vectorize inheritance deduction), #6 (replace iterrows), #11 (full vectorization)

**Sources:**
- [Pandas vectorization pitfalls](https://pythonspeed.com/articles/pandas-vectorization/)
- [Row-wise function application](https://www.architecture-performance.fr/ap_blog/applying-a-row-wise-function-to-a-pandas-dataframe/)

---

### Pitfall 3: Categorical Groupby 3500x Slowdown Without observed=True

**What goes wrong:**

Converting low-cardinality columns (CHROM, GENE, IMPACT) to categorical dtype saves 70-97% memory and speeds up operations **except groupby**. Without `observed=True`, pandas groupby on categorical columns creates groups for **every possible category** (even if not present in data), causing catastrophic slowdown.

```python
df["GENE"] = df["GENE"].astype("category")  # 2500 genes defined

# WITHOUT observed=True (DISASTER)
gene_groups = df.groupby("GENE")  # Creates 2500 groups even if df has 10 variants in 3 genes
# Each iteration processes empty DataFrames for 2497 genes → 3500x slower

# WITH observed=True (CORRECT)
gene_groups = df.groupby("GENE", observed=True)  # Only 3 groups
```

**Why it happens:**

Pandas defaults to `observed=False` for backward compatibility (pre-3.0). This is correct for statistical analysis (you want to see zero-count categories), but catastrophic for bioinformatics pipelines where you only care about genes with variants.

**Consequences:**

- Gene burden analysis on 2500-gene panel with 50 variants takes **hours instead of seconds**
- Memory explosion from creating empty DataFrames for every unobserved category
- Statistical aggregations produce nonsensical results (mean of empty group = NaN)

**Prevention:**

1. **ALWAYS use observed=True with categorical groupby:**
   ```python
   # Search and replace across entire codebase
   .groupby("GENE")  # WRONG
   .groupby("GENE", observed=True, sort=False)  # CORRECT
   ```

2. **Add lint rule:**
   ```python
   # Custom ruff rule or grep check in CI
   if "groupby(" in line and "observed=" not in line and column_is_categorical:
       raise Error("Missing observed=True")
   ```

3. **Regression test:**
   ```python
   def test_gene_burden_with_categorical_dtype():
       df["GENE"] = df["GENE"].astype("category")  # Force categorical
       start = time.time()
       result = compute_gene_burden(df)
       elapsed = time.time() - start
       assert elapsed < 5.0, f"Took {elapsed}s, likely missing observed=True"
   ```

4. **Note pandas 3.0 changes default:**
   In pandas 3.0+, `observed=True` becomes default. Add explicit `observed=True` now for forward compatibility and performance.

**Detection:**

- Tests hang or timeout on categorical columns
- Memory usage spikes to gigabytes for small DataFrames
- Warning: "FutureWarning: The default of observed=False is deprecated"

**Severity:** **CRITICAL** — makes optimization counterproductive

**Affected optimizations:** #3 (add observed=True), #5 (categorical dtypes)

**Sources:**
- [Categorical groupby slowdown](https://github.com/pandas-dev/pandas/issues/32976)
- [pandas 3.0 observed default change](https://pandas.pydata.org/docs/whatsnew/v3.0.0.html)
- [Categorical dtype gotchas](https://towardsdatascience.com/staying-sane-while-adopting-pandas-categorical-datatypes-78dbd19dcd8a)

---

### Pitfall 4: Dead Code Removal Breaking Validation Logic

**What goes wrong:**

Lines 220-249 in `gene_burden.py` parse GT strings but results are never used (line 251 says "Use the existing aggregated counts for now"). This looks like dead code to remove. **However:** the parsing loop may serve as validation that GT format is correct, even if results aren't used. Removing it could allow malformed data to pass silently.

```python
# Lines 220-249: Looks like dead code
for _, row in gene_df.iterrows():
    gt_value = str(row.get("GT", ""))
    for sample_entry in gt_value.split(";"):
        # Parse but results discarded
        sample_name = sample_entry.split("(")[0]
        genotype = ...
        # Sets not used!

# Line 251: Uses different columns entirely
p_var_count = int(gene_df["proband_variant_count"].sum())
```

**Why it happens:**

Code evolves. Original implementation manually parsed GT, later implementation used pre-computed columns, but the parsing loop remained as a "safety check" or was forgotten.

**Consequences if removed without verification:**

- Malformed GT strings that would have raised `IndexError` during parsing now pass silently
- No validation that GT format matches expected pattern
- Downstream crashes when other code expects valid GT format

**Prevention:**

1. **Confirm dead code is truly unused:**
   ```python
   # Instrument before removing
   removed_code_reached = False
   for _, row in gene_df.iterrows():
       removed_code_reached = True  # Set flag
       # ... parsing code ...
   if removed_code_reached:
       logger.warning("Code marked for removal was executed!")
   ```

2. **Extract validation if needed:**
   ```python
   # Replace dead loop with explicit validation
   def validate_gt_format(gt_string: str) -> None:
       for entry in gt_string.split(";"):
           if "(" not in entry or ")" not in entry:
               raise ValueError(f"Malformed GT: {entry}")

   # Apply once, fail fast
   df["GT"].apply(validate_gt_format)
   ```

3. **Add regression test for malformed input:**
   ```python
   def test_gene_burden_rejects_malformed_gt():
       df = make_gene_df()
       df.loc[0, "GT"] = "INVALID_FORMAT"  # No parentheses
       with pytest.raises(ValueError, match="Malformed GT"):
           compute_gene_burden(df)
   ```

4. **Check git blame and PR history:**
   Understand **why** code was written before removing it.

**Detection:**

- Tests pass after removal but production crashes on edge cases
- Different error messages (no ValueError on malformed GT)

**Severity:** **HIGH** — removes implicit validation, allows bad data

**Affected optimizations:** #2 (remove dead code)

**Sources:**
- Clinical validation requires explicit data quality checks
- [AMP/CAP NGS validation guidelines](https://www.sciencedirect.com/science/article/pii/S1525157817303732)

---

### Pitfall 5: xlsxwriter Incompatibility with Existing openpyxl Post-Processing

**What goes wrong:**

The pipeline uses openpyxl to **finalize** Excel files after initial write: adding hyperlinks, freeze panes, auto-filters. xlsxwriter is write-only and cannot open existing files for modification. Switching to xlsxwriter for initial write breaks the finalization step.

```python
# Current workflow (openpyxl)
df.to_excel(xlsx_file, engine="openpyxl")  # Write
finalize_excel_file(xlsx_file, config)     # Open with openpyxl, add formatting

# After switching to xlsxwriter (BROKEN)
df.to_excel(xlsx_file, engine="xlsxwriter")  # Write
finalize_excel_file(xlsx_file, config)       # openpyxl.load_workbook() fails!
# xlsxwriter doesn't support reading, openpyxl can't modify xlsxwriter files in-place
```

**Why it happens:**

xlsxwriter optimizes for write performance by streaming to disk. It never builds an in-memory workbook object, so cannot be reopened for editing.

**Consequences:**

- Hyperlinks missing (clinical users rely on SpliceAI/ClinVar links)
- No freeze panes (usability regression)
- No auto-filters (filtering by gene/impact broken)

**Prevention:**

1. **Do all formatting during initial write:**
   ```python
   # Use xlsxwriter for everything, including formatting
   with pd.ExcelWriter(xlsx_file, engine="xlsxwriter") as writer:
       df.to_excel(writer, sheet_name="Results", index=False)
       workbook = writer.book
       worksheet = writer.sheets["Results"]

       # Add formatting via xlsxwriter API
       worksheet.freeze_panes(1, 0)
       worksheet.autofilter(0, 0, len(df), len(df.columns)-1)

       # Add hyperlinks
       for row_idx in range(len(df)):
           if df.loc[row_idx, "SpliceAI_URL"]:
               worksheet.write_url(row_idx+1, col_idx, df.loc[row_idx, "SpliceAI_URL"])
   ```

2. **Hybrid approach (complex but works):**
   ```python
   # Write data with xlsxwriter (fast)
   df.to_excel(xlsx_file, engine="xlsxwriter", index=False)

   # Reopen with openpyxl for formatting (slower but familiar)
   wb = openpyxl.load_workbook(xlsx_file)
   # ... add formatting ...
   wb.save(xlsx_file)
   ```

3. **Benchmark both engines before committing:**
   ```python
   @pytest.mark.benchmark
   def test_excel_generation_speed(benchmark, large_df):
       benchmark(large_df.to_excel, "test.xlsx", engine="xlsxwriter")
       # Compare with openpyxl baseline
   ```

4. **Test Excel output structure:**
   ```python
   def test_excel_has_hyperlinks():
       wb = openpyxl.load_workbook(output_xlsx)
       ws = wb["Results"]
       assert ws["B2"].hyperlink is not None  # SpliceAI column
   ```

**Detection:**

- Missing hyperlinks in Excel output
- Missing freeze panes and auto-filters
- Error: "openpyxl.utils.exceptions.InvalidFileException"

**Severity:** **HIGH** — breaks clinical usability features

**Affected optimizations:** #8 (xlsxwriter for Excel)

**Sources:**
- [xlsxwriter vs openpyxl comparison](https://medium.com/@badr.t/excel-file-writing-showdown-pandas-xlsxwriter-and-openpyxl-29ff5bcb4fcd)
- [xlsxwriter alternatives](https://xlsxwriter.readthedocs.io/alternatives.html)

---

## High Severity Pitfalls

### Pitfall 6: Type Coercion Surprises When Removing dtype=str

**What goes wrong:**

The codebase uses `dtype=str` defensively to prevent pandas from inferring types incorrectly. Removing it enables optimizations (categoricals, PyArrow) but risks **silent type changes** that break string comparisons.

```python
# WITH dtype=str (current)
df = pd.read_csv("data.tsv", sep="\t", dtype=str)
df["POS"].dtype  # → dtype('O'), string "12345"
df["POS"] == "12345"  # → True

# WITHOUT dtype=str (RISKY)
df = pd.read_csv("data.tsv", sep="\t")
df["POS"].dtype  # → dtype('int64'), integer 12345
df["POS"] == "12345"  # → False! Type mismatch
```

**Why it happens:**

Pandas infers numeric columns as int64/float64. Code expecting strings breaks.

**Consequences:**

- Filters comparing to string literals fail (`df[df["POS"] == "12345"]` returns empty)
- Genotype matching breaks (`"0/1" == 0` is False)
- JSON serialization changes (int vs string)

**Prevention:**

1. **Explicit dtype mapping:**
   ```python
   GENOMIC_DTYPES = {
       "CHROM": "category",
       "POS": "Int64",  # Nullable int, not string
       "REF": "string[pyarrow]",
       "ALT": "string[pyarrow]",
       "GT": str,  # Keep as object/string
   }
   df = pd.read_csv(path, sep="\t", dtype=GENOMIC_DTYPES)
   ```

2. **Test type assumptions:**
   ```python
   def test_pos_is_numeric():
       df = load_variants()
       assert pd.api.types.is_integer_dtype(df["POS"])

   def test_gt_is_string():
       df = load_variants()
       assert pd.api.types.is_string_dtype(df["GT"])
   ```

3. **Gradual migration:**
   Convert one hot-path file at a time, test extensively before moving to next.

**Detection:**

- Failing tests with type mismatches
- Empty filter results
- `TypeError: '>' not supported between instances of 'str' and 'int'`

**Severity:** **HIGH** — silent correctness bugs

**Affected optimizations:** #5 (categorical dtypes), #4 (PyArrow engine)

**Sources:**
- [Pandas dtype pitfalls](https://notes.dsc80.com/content/02/data-types.html)

---

### Pitfall 7: Subprocess Pipe Fusion Breaking Error Handling

**What goes wrong:**

Current implementation calls `SnpSift filter`, checks return code, then calls `bgzip`. Fusing into a pipe makes error detection harder — if SnpSift fails, bgzip may still succeed, hiding the error.

```python
# BEFORE: Sequential with error checking
result = subprocess.run(["SnpSift", "filter", expr, input_vcf], check=True, capture_output=True)
if result.returncode != 0:
    raise RuntimeError(f"SnpSift failed: {result.stderr}")
subprocess.run(["bgzip", temp_vcf], check=True)

# AFTER: Piped (RISKY)
snpsift = subprocess.Popen(["SnpSift", "filter", expr, input_vcf], stdout=PIPE)
bgzip = subprocess.Popen(["bgzip", "-c"], stdin=snpsift.stdout, stdout=output_file)
bgzip.wait()  # But what if SnpSift failed?
```

**Why it happens:**

Unix pipe semantics: if early process fails but writes partial output, later process succeeds. Return code is only from final process.

**Consequences:**

- Truncated VCF files with no error
- Silent filter failures (output contains unfiltered variants)
- Non-reproducible results depending on where pipe breaks

**Prevention:**

1. **Check all return codes:**
   ```python
   snpsift.stdout.close()  # Allow SIGPIPE
   bgzip_rc = bgzip.wait()
   snpsift_rc = snpsift.wait()

   if snpsift_rc != 0:
       raise RuntimeError(f"SnpSift failed with code {snpsift_rc}")
   if bgzip_rc != 0:
       raise RuntimeError(f"bgzip failed with code {bgzip_rc}")
   ```

2. **Validate output:**
   ```python
   # After pipe completes, verify VCF is valid
   result = subprocess.run(["bcftools", "view", "-h", output_vcf], capture_output=True)
   if result.returncode != 0:
       raise RuntimeError("Output VCF is malformed")
   ```

3. **Use `set -e -o pipefail` if using shell:**
   ```python
   subprocess.run(
       "SnpSift filter 'expr' input.vcf | bgzip -c > output.vcf.gz",
       shell=True, check=True, executable="/bin/bash",
       env={**os.environ, "BASH_ENV": "set -e -o pipefail"}
   )
   ```

4. **Test pipe failure handling:**
   ```python
   @pytest.mark.integration
   def test_pipe_detects_snpsift_failure(tmp_path):
       malformed_vcf = tmp_path / "bad.vcf"
       malformed_vcf.write_text("NOT A VCF")

       with pytest.raises(RuntimeError, match="SnpSift failed"):
           apply_filter_piped(malformed_vcf, "1 == 1", config)
   ```

**Detection:**

- Truncated output files
- Different variant counts between piped and sequential
- Silent test failures

**Severity:** **HIGH** — hides tool failures

**Affected optimizations:** #7 (pipe fusion)

**Sources:**
- [Subprocess pipe chaining](https://rednafi.com/python/unix_style_pipeline_with_subprocess/)
- [bcftools pipe streaming](https://samtools.github.io/bcftools/bcftools.html)

---

### Pitfall 8: Cython/Rust Extension Breaking Cross-Platform Compatibility

**What goes wrong:**

Adding Cython or Rust extensions for genotype replacement requires C compiler on user machines and increases build complexity. Windows users face MSVC vs MinGW issues, conda environments break, and CI needs per-platform builds.

**Why it happens:**

Binary extensions are platform-specific. PyPI wheels must be built for Windows, macOS, Linux across Python 3.10, 3.11, 3.12.

**Consequences:**

- `pip install` fails on Windows with "error: Microsoft Visual C++ 14.0 is required"
- Conda environment cannot resolve dependencies
- Different behavior between platforms (Rust and Cython may have subtle differences)
- Maintenance burden for cross-compilation

**Prevention:**

1. **Use pure-Python optimization first:**
   ```python
   # Vectorized string operations often sufficient
   # 10-50x speedup without compilation
   ```

2. **If extension needed, provide wheels:**
   ```yaml
   # .github/workflows/wheels.yml
   - uses: pypa/cibuildwheel@v2.16
     env:
       CIBW_BUILD: "cp310-* cp311-* cp312-*"
       CIBW_SKIP: "*-win32 *-manylinux_i686"
   ```

3. **Make extension optional:**
   ```python
   try:
       from variantcentrifuge._genotype_replace_ext import replace_fast
       USE_NATIVE = True
   except ImportError:
       logger.warning("Native extension unavailable, using pure Python")
       USE_NATIVE = False
   ```

4. **Test on all platforms in CI:**
   ```yaml
   strategy:
     matrix:
       os: [ubuntu-latest, windows-latest, macos-latest]
       python-version: ["3.10", "3.12"]
   ```

5. **Document build requirements:**
   ```markdown
   ## Building from source (optional)

   Native extensions require:
   - Windows: Visual Studio 2022 with C++ tools
   - Linux: gcc 9+
   - macOS: Xcode command line tools
   ```

**Detection:**

- Build failures on clean Windows machine
- ImportError on fresh conda environment
- Different output between platforms

**Severity:** **HIGH** — breaks installation, platform compatibility

**Affected optimizations:** #10 (Cython/Rust extension)

**Sources:**
- [Cython cross-platform](https://cython-devel.python.narkive.com/Vy84tUt0/cython-cross-platform-extensions)
- [Rust vs Cython comparison](https://pythonspeed.com/articles/rust-cython-python-extensions/)

---

## Medium Severity Pitfalls

### Pitfall 9: Passing DataFrame Through Pipeline Breaking Checkpoint/Resume

**What goes wrong:**

Current pipeline writes TSV to disk between stages, enabling checkpoint/resume. Passing DataFrame in memory breaks this — if pipeline crashes after 2-hour inheritance analysis, must restart from scratch.

```python
# BEFORE: TSV on disk enables checkpointing
analyze_variants(input_vcf) → writes annotated.tsv
analyze_inheritance(annotated.tsv) → writes inheritance.tsv  # Can resume here!

# AFTER: In-memory DataFrame (LOST)
df = analyze_variants(input_vcf)  # 2 hours
df = analyze_inheritance(df)      # Crashes here → restart from beginning
```

**Why it happens:**

Optimization removes redundant disk I/O but loses crash recovery.

**Consequences:**

- Long-running jobs (>1 hour) cannot recover from errors
- OOM kills lose all intermediate results
- Debugging harder (no intermediate files to inspect)

**Prevention:**

1. **Hybrid: memory + checkpoint writes:**
   ```python
   df = analyze_variants(input_vcf)
   # Keep in memory BUT also checkpoint
   df.to_csv(checkpoint_path, sep="\t", index=False)
   context.data["variants_df"] = df  # In memory for next stage
   ```

2. **Use pipeline checkpoint system:**
   ```python
   # In runner.py
   if stage.supports_checkpoint and config.checkpoints_enabled:
       context.data["stage_df"].to_parquet(checkpoint_path)  # Fast binary format
   ```

3. **Memory-map DataFrames for crash recovery:**
   ```python
   # Use parquet with memory mapping
   df.to_parquet(checkpoint_path)
   df = pd.read_parquet(checkpoint_path, memory_map=True)  # Maps without full load
   ```

4. **Test crash recovery:**
   ```python
   def test_pipeline_resumes_after_crash(tmp_path):
       # Run pipeline to stage 3
       run_partial(stages=3)

       # Simulate crash, restart from checkpoint
       result = run_from_checkpoint(stage=3)
       assert result.equals(run_full())  # Same result
   ```

**Detection:**

- Long jobs restart from beginning after errors
- No intermediate files in output directory
- User complaints about lost progress

**Severity:** **MEDIUM** — usability regression for long jobs

**Affected optimizations:** #9 (pass DataFrame through pipeline)

---

### Pitfall 10: High-Cardinality Categorical Dtype Increasing Memory

**What goes wrong:**

Categorical dtype saves memory when cardinality is low (<10% of row count). For high-cardinality columns like sample IDs or variant keys, categorical uses **more** memory than object dtype.

```python
# LOW cardinality: CHROM has 25 values, 50K rows → 97% memory savings
df["CHROM"] = df["CHROM"].astype("category")  # GOOD

# HIGH cardinality: 50K unique variant keys for 50K rows → memory INCREASE
df["variant_key"] = df["variant_key"].astype("category")  # BAD
```

**Why it happens:**

Categorical stores: `categories array + codes array`. If categories ≈ row count, you store the data twice.

**Consequences:**

- Memory usage increases instead of decreases
- Slower operations (category lookup overhead)
- Misleading optimization results

**Prevention:**

1. **Check cardinality before converting:**
   ```python
   def should_categorize(series: pd.Series, threshold: float = 0.5) -> bool:
       n_unique = series.nunique()
       n_total = len(series)
       cardinality = n_unique / n_total
       return cardinality < threshold  # Only if <50% unique

   if should_categorize(df["GENE"]):
       df["GENE"] = df["GENE"].astype("category")
   ```

2. **Profile memory before/after:**
   ```python
   before = df.memory_usage(deep=True).sum()
   df["GENE"] = df["GENE"].astype("category")
   after = df.memory_usage(deep=True).sum()
   assert after < before, f"Categorical increased memory: {after} > {before}"
   ```

3. **Document which columns to categorize:**
   ```python
   # In GENOMIC_DTYPES mapping
   "CHROM": "category",      # ~25 values
   "IMPACT": "category",     # 4 values: HIGH, MODERATE, LOW, MODIFIER
   "GENE": "string[pyarrow]",  # Can be 2500+ unique, use PyArrow instead
   ```

**Detection:**

- Memory usage increases after optimization
- Slower groupby performance
- Profiler shows category overhead

**Severity:** **MEDIUM** — counterproductive optimization

**Affected optimizations:** #5 (categorical dtypes)

**Sources:**
- [Categorical dtype memory](https://pandas.pydata.org/docs/user_guide/categorical.html)
- [Categorical pitfalls](https://pbpython.com/pandas_dtypes_cat.html)

---

### Pitfall 11: Itertuples Name Collision with Reserved Attributes

**What goes wrong:**

`itertuples()` returns namedtuples with column names as attributes. If a column is named `count`, `index`, or other reserved names, access via `.count` returns the namedtuple method, not the column value.

```python
df = pd.DataFrame({"GENE": ["BRCA1"], "count": [5]})

for row in df.itertuples():
    print(row.count)  # → <bound method count of ...>, NOT 5!
```

**Why it happens:**

Namedtuples use `_make`, `_asdict`, `_replace` internally. Column names shadow these.

**Consequences:**

- Silent bugs where `.count` returns a method
- Type errors if code expects integer but gets callable
- Hard to debug (looks correct in code review)

**Prevention:**

1. **Use positional access or rename:**
   ```python
   for row in df.itertuples(index=False):
       gene = row[0]  # Positional
       count = row[1]

   # OR rename parameter
   for row in df.itertuples(index=False, name="Variant"):
       print(row.GENE, row.count)  # Now safe
   ```

2. **Check for reserved names:**
   ```python
   RESERVED = {"count", "index", "_fields", "_make", "_asdict", "_replace"}
   assert not (set(df.columns) & RESERVED), "Column names conflict with namedtuple"
   ```

3. **Test with real column names:**
   ```python
   def test_itertuples_with_real_columns():
       df = load_real_data()
       for row in df.itertuples(index=False):
           # Access each column by name
           _ = row.GENE, row.POS, row.REF  # Will fail if names conflict
   ```

**Detection:**

- Type errors: "expected int, got method"
- Wrong values silently used

**Severity:** **MEDIUM** — subtle correctness bug

**Affected optimizations:** #6 (replace iterrows with itertuples)

**Sources:**
- [itertuples documentation](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.itertuples.html)

---

### Pitfall 12: Regression Test False Negatives from Floating Point Comparison

**What goes wrong:**

When verifying "output identical before/after optimization," using `.equals()` or `==` for floating-point columns fails due to rounding differences in optimized code paths.

```python
# Before optimization
df["AF"] = pd.read_csv(...)["AF"]  # AF = 0.12345678901234567

# After optimization (PyArrow engine)
df["AF"] = pd.read_csv(..., engine="pyarrow")["AF"]  # AF = 0.12345678901234566

df_before.equals(df_after)  # → False! Difference in 16th decimal place
```

**Why it happens:**

Different parsing engines (C vs PyArrow) have subtly different floating-point rounding. Scientifically identical, but bitwise different.

**Consequences:**

- Optimization rejected due to "different output"
- Manual review needed for every floating-point column
- False alarm on clinically irrelevant differences

**Prevention:**

1. **Fuzzy comparison for floats:**
   ```python
   def dataframes_clinically_equal(df1, df2, rtol=1e-9) -> bool:
       # Exact match for strings, categoricals
       string_cols = df1.select_dtypes(include=["object", "string", "category"]).columns
       pd.testing.assert_frame_equal(df1[string_cols], df2[string_cols], check_exact=True)

       # Fuzzy match for numerics
       numeric_cols = df1.select_dtypes(include=[np.number]).columns
       pd.testing.assert_frame_equal(df1[numeric_cols], df2[numeric_cols], rtol=rtol)
   ```

2. **Round before comparison:**
   ```python
   # Round to clinical significance (e.g., allele frequency to 6 decimals)
   df1["AF"] = df1["AF"].round(6)
   df2["AF"] = df2["AF"].round(6)
   assert df1.equals(df2)
   ```

3. **Use checksums on sorted data:**
   ```python
   # Hash row tuples after rounding
   def row_hash(row):
       # Round floats to 10 decimals, then hash
       rounded = tuple(round(x, 10) if isinstance(x, float) else x for x in row)
       return hash(rounded)

   df["_hash"] = df.apply(row_hash, axis=1)
   assert df1["_hash"].equals(df2["_hash"])
   ```

**Detection:**

- Test failures on floating-point columns
- Diffs show changes in 10th+ decimal place
- Allele frequencies differ by <0.000001

**Severity:** **MEDIUM** — false negatives in testing

**Affected optimizations:** #4 (PyArrow engine), #5 (numeric dtypes)

**Sources:**
- [Pandas testing utilities](https://pandas.pydata.org/docs/reference/api/pandas.testing.assert_frame_equal.html)

---

## Low Severity Pitfalls

### Pitfall 13: Benchmark Flakiness on GitHub Actions

**What goes wrong:**

GitHub Actions runners have ~2.66% performance variance. A 5% performance gate fails 45% of the time even when code is unchanged.

**Prevention:**

Use ratio assertions and performance canaries (see issue assessment). Test relative performance (vectorized vs sequential) rather than absolute thresholds.

**Severity:** **LOW** — testing infrastructure annoyance

**Affected optimizations:** All (validation)

**Sources:**
- [GitHub Actions performance stability](https://aakinshin.net/posts/github-actions-perf-stability/)
- [Practical performance tests](https://solidean.com/blog/2025/practical-performance-tests/)

---

### Pitfall 14: xlsx File Size Increase with xlsxwriter

**What goes wrong:**

xlsxwriter produces larger files than openpyxl for identical data (~10-20% larger) due to different compression settings.

**Prevention:**

Acceptable tradeoff for 2-5x speed increase. Document expected file size change.

**Severity:** **LOW** — cosmetic issue

**Affected optimizations:** #8 (xlsxwriter)

---

## Summary: Risk Matrix

| Optimization | Critical Risk | High Risk | Medium Risk | Recommended Mitigation Priority |
|--------------|--------------|-----------|-------------|-------------------------------|
| #1: Vectorize `apply(axis=1)` | Order-dependent logic (#2) | Type coercion (#6) | Floating-point diffs (#12) | HIGH — test determinism |
| #2: Remove dead code | Dead code as validation (#4) | - | - | HIGH — confirm truly unused |
| #3: Add `observed=True` | 3500x slowdown without (#3) | - | - | CRITICAL — trivial fix, huge impact |
| #4: PyArrow engine | Pandas 3.0 dtype changes (#1) | Type coercion (#6) | Float precision (#12) | HIGH — test with pandas 3.0 |
| #5: Categorical dtypes | Pandas 3.0 + groupby (#3) | Type coercion (#6) | High cardinality (#10) | HIGH — cardinality checks |
| #6: Replace `iterrows` | Order-dependent logic (#2) | - | Name collision (#11) | MEDIUM — test column names |
| #7: Pipe fusion | - | Error handling (#7) | - | HIGH — validate pipe errors |
| #8: xlsxwriter | - | openpyxl incompatibility (#5) | File size (#14) | HIGH — test hyperlinks |
| #9: Pass DataFrame | - | - | Checkpoint loss (#9) | MEDIUM — hybrid approach |
| #10: Cython/Rust | - | Cross-platform (#8) | - | HIGH — wheels + optional |
| #11: Full vectorization | Order-dependent logic (#2) | Type coercion (#6) | - | CRITICAL — extensive testing |

---

## Validation Protocol

Before merging any optimization:

1. **Output identity check:**
   ```python
   df_before = run_pipeline_baseline(test_vcf)
   df_after = run_pipeline_optimized(test_vcf)
   assert dataframes_clinically_equal(df_before, df_after)
   ```

2. **All 1035 tests pass** (no skips, no warnings)

3. **Cross-platform validation:** Windows + Linux CI

4. **Benchmark shows improvement:**
   ```python
   assert time_after < time_before * 0.9  # At least 10% faster
   ```

5. **Memory profile shows improvement or no regression:**
   ```python
   assert memory_after <= memory_before * 1.1  # Allow 10% tolerance
   ```

6. **Manual review of 5 random output files** for hyperlinks, formatting, clinical correctness

---

## Key References

**Pandas 3.0 Breaking Changes:**
- [pandas 3.0.0 whatsnew](https://pandas.pydata.org/docs/whatsnew/v3.0.0.html)
- [String dtype breaking changes issue](https://github.com/pandas-dev/pandas/issues/59328)
- [pandas 3.0 InfoQ coverage](https://www.infoq.com/news/2026/02/pandas-library/)

**Performance Optimization:**
- [Pandas vectorization pitfalls](https://pythonspeed.com/articles/pandas-vectorization/)
- [PyArrow CSV benchmark](https://pythonspeed.com/articles/pandas-read-csv-fast/)
- [Categorical groupby slowdown](https://github.com/pandas-dev/pandas/issues/32976)
- [Categorical dtype gotchas](https://towardsdatascience.com/staying-sane-while-adopting-pandas-categorical-datatypes-78dbd19dcd8a)

**Clinical Validation:**
- [AMP/CAP NGS bioinformatics validation guidelines](https://www.sciencedirect.com/science/article/pii/S1525157817303732)
- [Assembling and validating bioinformatic pipelines](https://meridian.allenpress.com/aplm/article/144/9/1118/427496)
- [AI-powered medical software validation 2026](https://medium.com/@sginsbourg/ai-powered-medical-software-validation-in-2026-from-bottleneck-to-competitive-advantage-13f3535a65dd)

**Tools and Libraries:**
- [xlsxwriter vs openpyxl comparison](https://medium.com/@badr.t/excel-file-writing-showdown-pandas-xlsxwriter-and-openpyxl-29ff5bcb4fcd)
- [Cython cross-platform extensions](https://cython-devel.python.narkive.com/Vy84tUt0/cython-cross-platform-extensions)
- [Subprocess pipe chaining](https://rednafi.com/python/unix_style_pipeline_with_subprocess/)

**Testing:**
- [GitHub Actions performance stability](https://aakinshin.net/posts/github-actions-perf-stability/)
- [Practical CI-friendly performance tests](https://solidean.com/blog/2025/practical-performance-tests/)

---

## Research Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Pandas 3.0 string dtype changes | HIGH | Official docs + verified with pandas team GitHub issues |
| Categorical groupby slowdown | HIGH | Documented GitHub issue with reproducible example |
| Clinical validation requirements | MEDIUM | AMP/CAP guidelines verified, but 2026 AI practices from blog |
| Cross-platform Cython/Rust | MEDIUM | General knowledge + Cython docs, not variantcentrifuge-specific |
| xlsxwriter formatting limitations | HIGH | Official docs + cross-referenced multiple sources |
| Pipeline-specific risks | HIGH | Based on code review of actual implementation |

**Overall confidence:** HIGH — pitfalls are specific to planned optimizations with actionable prevention strategies.
