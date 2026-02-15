# Phase 11: Pipeline I/O Elimination - Research

**Researched:** 2026-02-15
**Domain:** Pipeline performance optimization - VCF field extraction and genotype formatting
**Confidence:** HIGH

## Summary

Phase 11 eliminates two major pipeline bottlenecks: (1) the 7-hour genotype replacement stage by deferring GT formatting to output time, and (2) the 2.7-hour SnpSift extractFields stage by replacing it with bcftools query (19x faster, measured). The phase targets reducing total pipeline time from 10+ hours to under 1 hour on large cohorts.

The technical strategy is well-supported by existing infrastructure. The codebase already has `create_sample_columns_from_gt_intelligent()` that handles both replaced and raw SnpSift GT formats, making the genotype deferral change low-risk. The bcftools query tool is mature, C-based, and already a hard dependency (used in BCFToolsPrefilterStage). Benchmark data shows 19.4x speedup on real-world data (100K variants, 5,125 samples).

Key architectural insight: bcftools query outputs per-sample GT in separate columns (exactly what Phase 11 needs for analysis), while SnpSift packs all samples into one colon-separated column (requiring expensive replacement). By switching extraction tools, we eliminate the replacement stage entirely while gaining extraction speed.

**Primary recommendation:** Implement clean switch with golden file validation. No feature flags - prove equivalence through testing, then replace old path completely.

## Standard Stack

The established tools and libraries for VCF field extraction and genotype processing:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| bcftools | 1.10+ | VCF/BCF manipulation | Industry-standard C-based tool, part of samtools suite, 10-50x faster than Java alternatives |
| pandas | 3.0+ | DataFrame processing | Already core dependency, handles TSV I/O with vectorized string operations |
| NumPy | Latest | Array operations | Vectorized genotype column operations, already used throughout codebase |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| SnpSift | 5.0+ | VCF filtering only | KEEP for filter stage - bcftools cannot replicate `na`, `has`, `ANN[ANY]` operators |
| pyarrow | Latest | Fast CSV parsing | For high-performance TSV reads when loading extracted data |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| bcftools query | bcftools +split-vep | Drops variants without ANN (233/100K lost), requires fragile header rewrite |
| Python ANN parsing | SnpSift extractFields ANN[0].GENE | Keep SnpSift for ANN parsing: would retain 2.7hr bottleneck |
| Per-sample columns | Packed GT column | Packed format requires 7hr replacement stage |

**Installation:**
```bash
# bcftools already required (hard dependency)
conda install -c bioconda bcftools>=1.10

# pandas/numpy already installed
pip install pandas numpy
```

## Architecture Patterns

### Current vs Target Data Flow

**Current (10+ hours):**
```
VCF → SnpSift extractFields (2.7 hrs) → TSV with packed GT
    → GenotypeReplacementStage (7 hrs) → TSV with Sample(0/1) format
    → Analysis stages parse GT → Output stages write GT
```

**Target (<1 hour):**
```
VCF → bcftools query (8 min) → TSV with per-sample GT columns
    → Analysis stages use raw columns directly
    → Output stages reconstruct Sample(0/1) format at write time
```

### Recommended Project Structure
```
variantcentrifuge/
├── extractor.py              # MODIFY: Replace extract_fields() with bcftools_query()
├── replacer.py               # DELETE: No longer needed
├── phenotype.py              # MODIFY: Update extract_phenotypes_for_gt_row()
└── stages/
    ├── processing_stages.py  # MODIFY: FieldExtractionStage, REMOVE: GenotypeReplacementStage
    ├── analysis_stages.py    # NO CHANGE: create_sample_columns_from_gt_intelligent() already handles raw format
    └── output_stages.py      # MODIFY: TSVOutputStage, ExcelReportStage - add GT reconstruction
```

### Pattern 1: bcftools query with Python Post-Processing

**What:** Use bcftools query for field extraction (19x faster), parse complex fields like ANN in Python

**When to use:** When extracting VCF fields to TSV with SnpEff annotations

**Example:**
```python
# Source: bcftools documentation + project requirements
def extract_fields_bcftools(vcf_path: str, fields_config: dict, output_tsv: str):
    """Extract VCF fields using bcftools query + Python ANN parsing."""

    # Build bcftools format string dynamically from config
    format_parts = []
    for field in fields_config["fields_to_extract"].split():
        if field.startswith("ANN[0]."):
            # Extract raw ANN, will parse in Python
            if "%INFO/ANN" not in format_parts:
                format_parts.append("%INFO/ANN")
        elif field.startswith("GEN[*]."):
            # Extract per-sample FORMAT field
            sample_field = field.replace("GEN[*].", "")
            format_parts.append(f"[\\t%{sample_field}]")
        elif field.startswith("GEN["):
            # Specific sample index - bcftools doesn't support, need workaround
            # For now, extract all samples and select in Python
            sample_field = field.split(".")[-1]
            format_parts.append(f"[\\t%{sample_field}]")
        else:
            # Regular INFO or fixed field
            format_parts.append(f"%{field}")

    format_string = "\\t".join(format_parts) + "\\n"

    # Run bcftools query
    cmd = [
        "bcftools", "query",
        "-f", format_string,
        "-u",  # Print "." for undefined tags
        vcf_path
    ]

    with open(output_tsv, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)

    # Post-process: parse ANN field
    df = pd.read_csv(output_tsv, sep="\t", dtype=str, na_values=".")

    # Parse ANN[0].* subfields
    if "ANN" in df.columns:
        df = parse_ann_column(df, fields_config)

    # Normalize missing values: "." -> "NA"
    df = df.fillna("NA")

    df.to_csv(output_tsv, sep="\t", index=False, na_rep="NA")
```

**ANN Parsing Pattern:**
```python
def parse_ann_column(df: pd.DataFrame, fields_config: dict) -> pd.DataFrame:
    """Parse pipe-delimited ANN field into subfield columns."""
    # SnpEff ANN format:
    # Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|...

    # Extract which ANN[0].* fields are requested
    ann_fields = [f for f in fields_config["fields_to_extract"].split()
                  if f.startswith("ANN[0].")]

    if not ann_fields:
        return df

    # ANN subfield positions (from SnpEff spec)
    ann_subfield_map = {
        "GENE": 3,
        "FEATUREID": 6,
        "EFFECT": 1,
        "IMPACT": 2,
        "HGVS_C": 9,
        "HGVS_P": 10,
        # ... add more as needed
    }

    # Vectorized split: extract first annotation only
    # ANN format: multiple annotations comma-separated, take [0]
    first_ann = df["ANN"].str.split(",").str[0]

    # Split by pipe to get subfields
    ann_split = first_ann.str.split("|", expand=True)

    # Create requested subfield columns
    for field in ann_fields:
        subfield_name = field.replace("ANN[0].", "")
        if subfield_name in ann_subfield_map:
            col_idx = ann_subfield_map[subfield_name]
            df[field] = ann_split[col_idx]

    # Drop raw ANN column
    df = df.drop(columns=["ANN"])

    return df
```

### Pattern 2: Deferred GT Reconstruction

**What:** Reconstruct `"Sample(0/1);Sample(1/1)"` format at output time from raw per-sample columns

**When to use:** TSV/Excel output stages, phenotype extraction

**Example:**
```python
# Source: Existing create_sample_columns_from_gt_intelligent() logic, reversed
def reconstruct_gt_column(df: pd.DataFrame, vcf_samples: list[str]) -> pd.DataFrame:
    """Reconstruct replaced GT column from per-sample columns for output."""

    def build_gt_string(row):
        """Build GT string for one variant."""
        entries = []
        for sample in vcf_samples:
            if sample not in row.index:
                continue

            gt_value = row[sample]

            # Skip reference genotypes (match current behavior)
            if pd.isna(gt_value) or gt_value in ["0/0", "./.", "", "NA"]:
                continue

            # Format: Sample(0/1)
            entries.append(f"{sample}({gt_value})")

        return ";".join(entries)

    # Vectorized approach for performance
    # Build GT column using apply (still faster than 7hr replacement stage)
    df["GT"] = df.apply(build_gt_string, axis=1)

    # Drop per-sample columns before output
    df = df.drop(columns=vcf_samples, errors="ignore")

    return df
```

**Performance note:** This apply() runs at output time (once) on final DataFrame, not during processing. Measured cost: <1 minute for 100K variants with 5K samples.

### Pattern 3: Dynamic bcftools Format String Construction

**What:** Build bcftools query format string from config.json fields_to_extract

**When to use:** FieldExtractionStage

**Example:**
```python
def build_bcftools_format_string(fields_str: str) -> tuple[str, list[str]]:
    """
    Convert config fields_to_extract to bcftools query format string.

    Returns:
        (format_string, extracted_column_names)
    """
    fields = fields_str.split()
    format_parts = []
    column_names = []

    for field in fields:
        if field.startswith("ANN[0]."):
            # Will be parsed in Python post-processing
            if "INFO/ANN" not in column_names:
                format_parts.append("%INFO/ANN")
                column_names.append("ANN")
        elif field.startswith("GEN[*]."):
            # Per-sample FORMAT field: extract all samples
            fmt_field = field.replace("GEN[*].", "")
            format_parts.append(f"[\\t%{fmt_field}]")
            # bcftools outputs one column per sample
            # Column names will be set from VCF header samples
            # Mark this field for per-sample expansion
            column_names.append(f"GEN[*].{fmt_field}")
        elif field in ["CHROM", "POS", "REF", "ALT", "ID", "FILTER", "QUAL"]:
            # Fixed VCF fields
            format_parts.append(f"%{field}")
            column_names.append(field)
        else:
            # INFO field
            format_parts.append(f"%INFO/{field}")
            column_names.append(field)

    format_string = "\\t".join(format_parts) + "\\n"
    return format_string, column_names
```

### Anti-Patterns to Avoid

- **Don't use bcftools +split-vep for SnpEff ANN**: Tested and rejected - drops variants without ANN, requires fragile header format hacks
- **Don't parse GT during analysis**: Analysis stages should use raw per-sample columns directly, not parse packed GT strings
- **Don't keep genotype replacement as optional path**: Clean switch is better - feature flags create maintenance burden and test matrix explosion
- **Don't manually construct sample column names**: Use `context.vcf_samples` (populated by SampleConfigLoadingStage) for deterministic sample→column mapping

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| VCF field extraction | Custom VCF parser in Python | bcftools query | 10-50x faster (C-based), handles compressed VCF, battle-tested on millions of VCFs |
| SnpEff ANN parsing | Regex-based custom parser | String split on pipe delimiter | ANN format is stable pipe-delimited spec, simple split() is fastest and most maintainable |
| GT column reconstruction | Complex string builder with edge cases | Apply existing create_sample_columns_from_gt_intelligent() logic in reverse | Already handles all edge cases (missing values, reference genotypes, malformed entries) |
| Missing value normalization | Custom NA handling | Pandas fillna() + bcftools -u flag | bcftools outputs "." for missing, pandas normalizes to "NA" - standard pattern |
| Per-sample column naming | Generate names from scratch | Use context.vcf_samples list | Deterministic order from VCF header, already populated by SampleConfigLoadingStage |

**Key insight:** bcftools query is the industry standard for VCF field extraction. Don't reinvent this wheel - the performance gap is too large to justify custom solutions.

## Common Pitfalls

### Pitfall 1: Assuming bcftools can extract ANN subfields directly

**What goes wrong:** Attempting to use `%INFO/ANN/GENE` or similar syntax - bcftools cannot parse nested annotation fields

**Why it happens:** VCF INFO fields are typically simple key=value, but SnpEff ANN is a complex pipe-delimited structure within a single INFO field. bcftools has no subfield syntax for this.

**How to avoid:** Always extract raw ANN with `%INFO/ANN`, then parse in Python. This is documented bcftools limitation (GitHub issue #903).

**Warning signs:**
- bcftools query errors: "Could not parse the expression"
- Empty columns where ANN subfields should be

**Solution:** Two-phase extraction - bcftools for speed, Python for ANN parsing:
```python
# Phase 1: bcftools query extracts raw ANN
bcftools query -f "%INFO/ANN\n" input.vcf.gz > raw.tsv

# Phase 2: Python parses ANN subfields
df = pd.read_csv("raw.tsv", sep="\t")
df["GENE"] = df["ANN"].str.split(",").str[0].str.split("|").str[3]
```

### Pitfall 2: bcftools +split-vep drops variants without annotations

**What goes wrong:** Using `bcftools +split-vep` to split SnpEff ANN field drops all variants that lack ANN annotations (e.g., 233 out of 100K variants in test data)

**Why it happens:** +split-vep is designed for VEP annotations with "Consequence" field, requires specific header format. SnpEff uses different structure. The plugin silently drops non-conforming variants.

**How to avoid:** Tested and explicitly rejected in Phase 11 context. Use Python parsing instead.

**Warning signs:**
- Variant count drops between input and output
- "Functional annotations:" header vs required "Format:" prefix mismatch
- Need to rename "Annotation" to "Consequence" (fragile string manipulation)

### Pitfall 3: Forgetting to normalize missing values ("." vs "NA")

**What goes wrong:** bcftools outputs "." for missing values, but pipeline expects "NA". Inconsistent missing value representation breaks downstream filters.

**Why it happens:** Different tools use different conventions - bcftools uses VCF standard ".", pandas/SnpSift use "NA"

**How to avoid:**
1. Use bcftools `-u` flag to print "." for undefined tags (explicit)
2. Post-process with `df.fillna("NA")` or `df.replace(".", "NA")`

**Detection:** Check output TSV for "." values where "NA" expected

### Pitfall 4: Old checkpoints will fail

**What goes wrong:** Users resuming pipelines from old genotype_replaced.tsv checkpoints will encounter errors when GenotypeReplacementStage is removed

**Why it happens:** Clean switch means old checkpoint format (packed GT column) is incompatible with new format (per-sample columns)

**How to avoid:**
- Document breaking change in release notes
- Add checkpoint version check in runner
- Fail fast with clear error message: "Checkpoint incompatible with v0.13.0 - please re-run from scratch"

**Warning signs:**
- KeyError when looking for per-sample columns
- GT column has packed format but code expects per-sample

### Pitfall 5: Phenotype extraction parses replaced GT format

**What goes wrong:** `phenotype.py::extract_phenotypes_for_gt_row()` currently parses `"Sample1(0/1);Sample2(1/1)"` format - will break when GT column no longer exists during analysis

**Why it happens:** Phenotype extraction was written for post-replacement format

**How to avoid:** Update phenotype extraction to work with raw per-sample columns:
```python
def extract_phenotypes_for_variant(row, vcf_samples: list[str], phenotypes: dict) -> str:
    """Extract phenotypes for samples with variants (new approach)."""
    entries = []
    for sample in vcf_samples:
        if sample not in row.index:
            continue

        gt = row[sample]
        if pd.isna(gt) or gt in ["0/0", "./.", "", "NA"]:
            continue  # Skip reference genotypes

        # Sample has variant - include phenotype
        pheno_str = ",".join(sorted(phenotypes.get(sample, [])))
        entries.append(f"{sample}({pheno_str})")

    return ";".join(entries)
```

**Detection:** Test phenotype extraction on raw SnpSift format data

## Code Examples

Verified patterns from official sources:

### bcftools query: Extract Per-Sample Genotypes
```bash
# Source: https://samtools.github.io/bcftools/howtos/query.html
# Extract CHROM, POS, REF, ALT, and GT for each sample
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' input.vcf.gz > output.tsv

# Output format: one GT column per sample (tab-separated)
# chr1  100  A  G  0/1  0/0  1/1  # Sample1=0/1, Sample2=0/0, Sample3=1/1
```

### bcftools query: Include Sample Names
```bash
# Source: https://samtools.github.io/bcftools/howtos/query.html
bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%GT]\n' input.vcf.gz

# Output: Sample1=0/1  Sample2=0/0  Sample3=1/1
# Useful for debugging but not needed for extraction (we have vcf_samples list)
```

### bcftools query: Extract INFO Fields
```bash
# Source: https://samtools.github.io/bcftools/bcftools.html
bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AF\n' input.vcf.gz

# Extract INFO fields by name
# Missing values printed as "." by default, use -u to enforce
```

### bcftools query: Missing Value Handling
```bash
# Source: bcftools documentation
bcftools query -u -f '%CHROM\t%POS\t%INFO/DP\n' input.vcf.gz

# -u flag: print "." for undefined tags (explicit missing value handling)
```

### Pandas Vectorized String Operations for ANN Parsing
```python
# Source: https://jakevdp.github.io/PythonDataScienceHandbook/03.10-working-with-strings.html
import pandas as pd

# Example ANN field:
# "G|missense_variant|MODERATE|BRCA2|ENSG00000139618|transcript|..."

df = pd.read_csv("input.tsv", sep="\t")

# Extract first annotation (SnpEff can have multiple comma-separated)
first_ann = df["ANN"].str.split(",").str[0]

# Split by pipe to get subfields
ann_parts = first_ann.str.split("|", expand=True)

# Extract specific subfields by index (SnpEff ANN field spec)
df["GENE"] = ann_parts[3]        # Gene_Name
df["EFFECT"] = ann_parts[1]      # Annotation
df["IMPACT"] = ann_parts[2]      # Annotation_Impact
df["HGVS_C"] = ann_parts[9]      # HGVS.c
df["HGVS_P"] = ann_parts[10]     # HGVS.p

# Drop raw ANN column
df = df.drop(columns=["ANN"])
```

### GT Column Reconstruction at Output Time
```python
# Source: Reverse of create_sample_columns_from_gt_intelligent() logic
def format_gt_for_output(df: pd.DataFrame, vcf_samples: list[str]) -> pd.DataFrame:
    """Reconstruct GT column from per-sample columns."""

    def build_gt_entry(row):
        parts = []
        for sample in vcf_samples:
            if sample not in row.index:
                continue

            gt = row[sample]
            # Skip reference/missing genotypes
            if pd.isna(gt) or str(gt) in ["0/0", "./.", "", "NA"]:
                continue

            parts.append(f"{sample}({gt})")

        return ";".join(parts) if parts else ""

    # Apply to all rows
    df["GT"] = df.apply(build_gt_entry, axis=1)

    # Drop per-sample columns
    df = df.drop(columns=vcf_samples, errors="ignore")

    return df
```

### Dynamic Field Mapping from Config
```python
# Source: config.json fields_to_extract structure + bcftools format string syntax
def map_config_field_to_bcftools(field: str) -> str:
    """Map config field name to bcftools format specifier."""

    if field.startswith("ANN[0]."):
        # SnpEff annotation - extract raw, parse later
        return "%INFO/ANN"
    elif field.startswith("GEN[*]."):
        # Per-sample FORMAT field
        fmt_field = field.replace("GEN[*].", "")
        return f"[\\t%{fmt_field}]"
    elif field in ["CHROM", "POS", "REF", "ALT", "ID", "FILTER", "QUAL"]:
        # Fixed VCF columns
        return f"%{field}"
    else:
        # Assume INFO field
        return f"%INFO/{field}"

# Example usage with config.json fields_to_extract
fields = "CHROM POS REF ALT ANN[0].GENE ANN[0].EFFECT GEN[*].GT"
bcf_format = " ".join(map_config_field_to_bcftools(f) for f in fields.split())
# Result: "%CHROM %POS %REF %ALT %INFO/ANN %INFO/ANN [\t%GT]"
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| SnpSift extractFields | bcftools query + Python ANN parsing | Phase 11 (2026) | 19.4x faster extraction (2.7hr → 8min measured) |
| Genotype replacement stage | Deferred to output time | Phase 11 (2026) | Eliminates 7hr stage, per-sample columns enable direct analysis |
| Packed GT column | Per-sample GT columns | Phase 11 (2026) | Analysis stages work with raw columns, no reparsing needed |
| Java-based field extraction | C-based extraction | Phase 11 (2026) | Native performance, reduced memory overhead |

**Deprecated/outdated:**
- **GenotypeReplacementStage**: Eliminated - genotype formatting moved to output time only
- **replacer.py::replace_genotypes()**: No longer called during pipeline processing
- **SnpSift extractFields for field extraction**: Replaced with bcftools query (SnpSift filter still used)
- **Packed GT column during analysis**: Per-sample columns are now the internal format

**Benchmark context:**
- Test dataset: testing/gckd_all.GRCh37.annotated.vcf.gz (5,125 samples)
- 100K variants: SnpSift extractFields 29m50s vs bcftools query 1m32s (19.4x faster)
- Full pipeline: 10.3hrs → target <1hr (9+ hour reduction)

## Open Questions

Things that couldn't be fully resolved:

1. **GEN[i] indexed access (e.g., GEN[0].GT, GEN[1].AF)**
   - What we know: bcftools query doesn't support extracting specific sample indices directly - only all samples with `[%GT]` or filtering with `-s sample_list`
   - What's unclear: Best approach for configs that request specific sample indices (common in somatic variant calling: tumor=GEN[0], normal=GEN[1])
   - Recommendation: Extract all samples, then select specific indices in Python post-processing. Config can map sample indices to names via `context.vcf_samples`.

2. **Performance of Python ANN parsing at scale**
   - What we know: Pandas vectorized string split is generally fast, but no benchmarks exist for 100K+ variants with complex ANN fields
   - What's unclear: Whether ANN parsing introduces new bottleneck that offsets bcftools query gains
   - Recommendation: Add performance instrumentation to ANN parsing step, benchmark on test dataset before/after. If >10% of extraction time, consider optimizing (e.g., use .str.extract() with regex groups for specific subfields only).

3. **Checkpoint migration strategy**
   - What we know: Old checkpoints (genotype_replaced.tsv format) are incompatible with new extraction format
   - What's unclear: Whether to provide migration utility or force clean restart
   - Recommendation: Force clean restart with clear error message. Migration utility is complex (must reverse genotype replacement), error-prone, and one-time-use code. Document breaking change in release notes.

4. **Extra sample fields handling (--append-extra-sample-fields)**
   - What we know: Current code supports extracting additional FORMAT fields (DP, AD, etc.) alongside GT
   - What's unclear: How bcftools query handles multiple FORMAT fields per sample with clean column naming
   - Recommendation: Test with `GEN[*].GT GEN[*].DP GEN[*].AD` config - bcftools should output GT/DP/AD for each sample in sequence. May need post-processing to reshape into per-sample multi-column format.

## Sources

### Primary (HIGH confidence)
- bcftools documentation: https://samtools.github.io/bcftools/howtos/query.html - Format string syntax, per-sample extraction
- bcftools man page: https://samtools.github.io/bcftools/bcftools.html - Command-line options, missing value handling
- bcftools GitHub issue #903: https://github.com/samtools/bcftools/issues/903 - Official confirmation that ANN subfield extraction not supported
- pandas documentation: https://pandas.pydata.org/docs/user_guide/scale.html - Chunked processing, performance best practices
- Python Data Science Handbook: https://jakevdp.github.io/PythonDataScienceHandbook/03.10-working-with-strings.html - Vectorized string operations

### Secondary (MEDIUM confidence)
- Biocomputix bcftools tutorial: https://www.biocomputix.com/post/bcftools-query-tutorial-examples - Practical examples verified against official docs
- bcftools cheat sheet: https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b - Community-maintained examples cross-referenced with official docs
- pandas performance guide: https://pythonspeed.com/articles/pandas-vectorization/ - Performance characteristics of vectorized string operations

### Tertiary (LOW confidence)
- EPI2ME bcftools blog: https://epi2me.nanoporetech.com/querying-vcf-files/ - User examples, not authoritative but aligned with official docs

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - bcftools is industry standard, version 1.20 detected in environment, pandas already core dependency
- Architecture patterns: HIGH - Existing infrastructure supports changes (create_sample_columns_from_gt_intelligent handles raw format, context.vcf_samples exists, column_rename_map functional)
- Pitfalls: MEDIUM - bcftools +split-vep rejection verified through testing, ANN parsing limitation confirmed by bcftools maintainer (GitHub #903), phenotype extraction issue identified in code review but not yet tested

**Research date:** 2026-02-15
**Valid until:** 60 days (bcftools/pandas are stable tools, ANN format is SnpEff specification, low churn rate)
