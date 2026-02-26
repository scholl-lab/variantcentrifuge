# Association Performance Investigation (2026-02-25)

## Problem

Running COAST association on GCKD cohort (5,125 samples, 5,084 genes, `coding` + `1percent`)
takes hours. After 37 min of association stage execution, only ~44% of genes had been
processed — and the **association engine hadn't even started yet**. The bottleneck is in
genotype matrix construction, not COAST itself.

**Observed:** 207,351 variants x 5,159 columns (5,125 per-sample GT + 34 annotation), 16.4 GB RAM

---

## 1. The Antipattern: Triple Genotype Parsing

The pipeline parses genotype strings **3 separate times**, each module independently and
each discarding its work:

| Parse | Stage | Method | Encoding | Fate |
|-------|-------|--------|----------|------|
| #1 | Inheritance (Pass 1) | `encode_genotypes()` vectorized NumPy | int8 (0/1/2/-1) | **DISCARDED** after compound het |
| #2 | Gene Burden | `.map(_gt_to_dosage)` vectorized | carrier counts only | Counts kept, GTs dropped |
| #3 | Association | `build_genotype_matrix()` **iterrows()** | float64 with imputation | Used for COAST/SKAT |

### Why This Happened (Architecture Archaeology)

1. **Independent module development without cross-module awareness:**
   - Inheritance module (Phase 9, 2026-02-15): Built `encode_genotypes()` and
     `_build_genotype_matrix()` as internal utilities. Self-contained, encapsulated.
   - Association module (Phase 19, 2026-02-20): Built `genotype_matrix.py` from scratch,
     200 lines, unaware of inheritance's existing encoder.
   - Within inheritance, Pass 1→2 matrix reuse was added (commit 55358e0), but cross-stage
     reuse was never considered.

2. **PipelineContext carries metadata, not data artifacts:**
   - `context.current_dataframe` and `context.variants_df` store DataFrames
   - No mechanism for passing computed arrays (matrices) between stages
   - No `stage_results` convention for large intermediate data

3. **Encoding incompatibility (real but solvable):**
   - Inheritance: int8, all variants, -1 for missing (needed for segregation checks)
   - Association: float64, per-gene filtered, MAF-imputed (needed for regression)
   - Different filtering: association applies site missingness filter (>10% missing → drop)
   - These differences are **real** but the shared work (GT string → dosage integer) is
     the expensive part that could be factored out.

4. **GT column lifecycle creates unnecessary work:**
   ```
   bcftools → per-sample GT columns (GEN_0__GT, ..., 5125 cols)
     → InheritanceStage: parses to int8 matrix, DISCARDS
     → VariantAnalysisStage: PACKS per-sample cols → single GT string, DROPS originals
     → GeneBurdenStage: works with packed GT
     → AssociationStage: RECOVERS per-sample cols from context.variants_df cache
       → RE-PARSES each GT string via iterrows() (the bottleneck)
   ```

### The Two Specific Bottlenecks

**Bottleneck A: Per-gene DataFrame filter** — `analysis_stages.py:2578`
```python
gene_df = gt_source_df[gt_source_df["GENE"] == gene_name]  # 207K-row scan per gene
```
O(n_genes x n_rows) = 5,084 x 207,351 ~ **1 billion comparisons**

**Bottleneck B: iterrows() GT parsing** — `genotype_matrix.py:173-180`
```python
for v_idx, (_, row) in enumerate(gene_df.iterrows()):       # O(n_variants)
    for s_idx, col in enumerate(gt_columns_list):            # O(5125 samples)
        dosage, is_multi = parse_gt_to_dosage(str(gt_val))   # Python call per cell
```
~1 billion `parse_gt_to_dosage()` calls total across all genes.

### Why the 500-Gene Run Was Fast

| Run | Genes | Total Rows | Per-gene filter cost | Total parse calls | Time |
|-----|-------|------------|---------------------|-------------------|------|
| 500-gene | ~500 | ~20K | 500 x 20K = 10M | ~100M | ~30s |
| 5K-gene | ~5,084 | ~207K | 5K x 207K = 1B | ~1B+ | ~90 min |

Scales **quadratically** — 10x genes and 10x rows = 100x total work.

---

## 2. Resource Allocation Gap

### What Exists: ResourceManager (Well-Designed)

`variantcentrifuge/memory/resource_manager.py` (378 lines) provides:
- **Memory detection** (5-level priority): CLI `--max-memory-gb` > SLURM > PBS > cgroups > psutil
- **CPU detection** (5-level priority): SLURM > PBS > psutil > os.cpu_count() > fallback(4)
- `auto_workers(task_count, memory_per_task_gb)` — memory+CPU-constrained worker count
- `auto_chunk_size(total_items, num_samples)` — memory-safe chunk sizes
- `should_parallelize(total_items)` — threshold for parallelization (default >=100 items)

### How `--threads` Currently Flows

```
CLI: --threads auto (or N)
  ↓
cli.py:1438-1449: ResourceManager resolves to integer
  ↓
context.config["threads"] = resolved_integer
  ↓
Stages consume independently:
  ├── ParallelCompleteProcessing:  config["threads"] directly        ❌ No ResourceManager
  ├── InheritanceAnalysis:         config["threads"] directly        ⚠️  Sometimes uses RM
  ├── ChunkedAnalysis:             rm.auto_workers() capped          ✅ Uses ResourceManager
  ├── GeneBurdenAnalysis:          rm.auto_workers() capped          ✅ Uses ResourceManager
  ├── AssociationAnalysis:         config["association_workers"]     ❌ Separate parameter!
  └── IGVReport:                   min(config["threads"], 8)         ❌ Hardcoded cap
```

### The Gap: `--threads` Does NOT Flow to Association

- `--association-workers` defaults to **1** (hardcoded in `analysis_stages.py:2205`)
- `--threads` is consumed by bcftools, SnpSift, inheritance, and chunked processing
- Association workers are completely separate — no auto-allocation from `--threads`
- ResourceManager is never consulted for association worker count

### Best Practice: Snakemake/Nextflow Thread Budget Model

Industry-standard workflow managers (Snakemake 9.x, Nextflow) use a **global thread budget**
with per-stage allocation:

- **Global budget**: `--cores N` sets total available cores
- **Per-rule declaration**: Each rule declares `threads: M` (max it can use)
- **Scheduler constraint**: Sum of running threads never exceeds global budget
- **Dynamic allocation**: Rules can use `threads: workflow.cores * 0.75` (fraction of total)
- **Resource-aware**: Memory constraints can further reduce thread count

Our pipeline runs stages **sequentially** (one stage at a time within each dependency level),
so the Snakemake concern about concurrent stages sharing a budget is less relevant. But
the pattern of **declaring resource needs per stage and having a central allocator** is
directly applicable.

**Reference pattern (Snakemake):**
```python
rule bwa_map:
    threads: workflow.cores    # Uses all available cores
    resources:
        mem_mb=8000            # Memory constraint

rule samtools_sort:
    threads: 8                 # Fixed max threads
```

This maps to our pipeline as: each stage should declare its thread needs, and
`ResourceManager.auto_workers()` should be the single gatekeeper, consulted by ALL stages.

Sources:
- [Snakemake Advanced Tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html)
- [Snakemake Rules — threads and resources](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)
- [Nextflow HPC Optimization](https://seqera.io/blog/optimizing-nextflow-for-hpc-and-cloud-at-scale/)
- [Bioinformatics Pipeline Frameworks 2025](https://www.tracer.cloud/resources/bioinformatics-pipeline-frameworks-2025)
- [Reproducible pipelines with workflow managers (Nature Methods)](https://www.nature.com/articles/s41592-021-01254-9)
- [PyExPool: Resource-constrained execution pool](https://github.com/eXascaleInfolab/PyExPool)

---

## 3. Proposed Fix Plan

### Phase 35 Scope: Association Performance & Resource Normalization

#### Fix 1: Vectorize `build_genotype_matrix()` [CRITICAL]

Replace `iterrows()` nested loop with vectorized NumPy/Pandas operations.
Reuse the pattern from `encode_genotypes()` (inheritance) which uses NumPy boolean masks:

```python
# encode_genotypes() pattern (already proven fast):
result = np.full(n, -1, dtype=np.int8)
result[series == '0/0'] = 0
result[series.isin(['0/1', '1/0'])] = 1
result[series == '1/1'] = 2
```

Adapt for association needs (float64 + multi-allelic + NaN for missing):

```python
def _vectorized_parse_gt_column(gt_series: pd.Series) -> tuple[np.ndarray, int]:
    """Parse an entire GT column to float64 dosage array. Returns (dosages, multi_count)."""
    s = gt_series.astype(str).str.replace('|', '/', regex=False)
    dosage = np.full(len(s), np.nan, dtype=np.float64)
    dosage[s == '0/0'] = 0.0
    dosage[s.isin(['0/1', '1/0'])] = 1.0
    dosage[s == '1/1'] = 2.0
    # Multi-allelic: any allele > 1 → het-equivalent dosage 1
    multi_mask = s.str.match(r'.*[2-9].*', na=False)
    dosage[multi_mask] = 1.0
    return dosage, int(multi_mask.sum())
```

Then build matrix column-by-column (outer loop over samples, vectorized inner):
```python
raw = np.full((n_variants, n_samples), np.nan, dtype=np.float64)
multi_count = 0
for s_idx, col in enumerate(gt_columns_list):
    dosages, mc = _vectorized_parse_gt_column(gene_df[col])
    raw[:, s_idx] = dosages
    multi_count += mc
```

**Expected speedup:** 50-100x (eliminates ~1B Python function calls)

**Files:** `association/genotype_matrix.py`

#### Fix 2: GroupBy instead of per-gene filter [CRITICAL]

```python
# Current: O(n_genes x n_rows) — 1 billion comparisons
for gene_data in gene_burden_data:
    gene_df = gt_source_df[gt_source_df["GENE"] == gene_name]

# Proposed: O(n_rows) single pass + O(1) per lookup
gene_groups = gt_source_df.groupby("GENE")
for gene_data in gene_burden_data:
    try:
        gene_df = gene_groups.get_group(gene_name)
    except KeyError:
        gene_df = gt_source_df.iloc[0:0]  # Empty DF
```

**Expected speedup:** ~5000x for the filtering step

**Files:** `stages/analysis_stages.py:2570-2600`

#### Fix 3: Normalize `--threads` → `association_workers` via ResourceManager [HIGH]

Make `association_workers` default to auto-allocation via ResourceManager instead of
hardcoded 1. The `--association-workers` CLI flag becomes an **override**, not the primary
source. This follows the Snakemake pattern where per-rule threads derive from the global
budget unless explicitly overridden.

```python
# In _build_assoc_config_from_context():
raw_workers = _get("association_workers", default=None, nullable=True)
if raw_workers is None or raw_workers == 0:
    # Auto-allocate from --threads via ResourceManager
    rm = ResourceManager(config=cfg)
    # Estimate memory per gene: ~50MB for 5K samples (genotype matrix + SKAT/COAST)
    mem_per_gene = 0.05  # GB
    auto_workers = rm.auto_workers(
        task_count=1000,  # Typical gene panel size
        memory_per_task_gb=mem_per_gene,
    )
    association_workers = max(1, min(auto_workers, int(cfg.get("threads", 1))))
else:
    association_workers = raw_workers  # Explicit override
```

Require consistent ResourceManager usage across ALL stages:

| Stage | Current | Target |
|-------|---------|--------|
| ParallelCompleteProcessing | `config["threads"]` directly | `rm.auto_workers()` |
| InheritanceAnalysis | Sometimes RM | Always `rm.auto_workers()` |
| ChunkedAnalysis | `rm.auto_workers()` | No change (already correct) |
| GeneBurdenAnalysis | `rm.auto_workers()` | No change (already correct) |
| AssociationAnalysis | Hardcoded 1 | `rm.auto_workers()` with threads cap |
| IGVReport | `min(threads, 8)` | `rm.auto_workers()` with hardcoded cap |

**Files:** `stages/analysis_stages.py`, `cli.py`, `association/base.py`

#### Fix 4: Store ResourceManager in PipelineContext [MEDIUM]

Currently each stage creates its own ResourceManager instance. Add a shared instance
to PipelineContext so all stages use identical detection results and safety factors:

```python
# In pipeline_core/context.py:
@dataclass
class PipelineContext:
    ...
    resource_manager: ResourceManager | None = None

# In pipeline.py (setup):
context.resource_manager = ResourceManager(config=context.config)
```

Stages then use `context.resource_manager.auto_workers()` instead of creating local instances.

**Files:** `pipeline_core/context.py`, `pipeline.py`, all stages

#### Fix 5: Eliminate unnecessary GT column drop/recover cycle [MEDIUM]

The current flow drops per-sample GT columns in VariantAnalysisStage, then recovers them
from `context.variants_df` in AssociationStage. This is wasteful — keep per-sample columns
alive through the pipeline and only drop them at output time:

- Move `reconstruct_gt_column()` to output stages (TSVOutputStage, ExcelReportStage) only
- Keep per-sample GT columns in `current_dataframe` through analysis stages
- Eliminates the recovery fallback and the 16 GB duplicate DataFrame

**Files:** `stages/analysis_stages.py` (VariantAnalysis, GeneBurden, Association),
`stages/output_stages.py`

#### Fix 6: Stream genotype matrix construction [LOW]

Currently all 5K+ matrices are built upfront before association engine starts (stored
in `gene_burden_data` list), causing 16+ GB peak memory.

Build matrix per-gene just before test execution, discard after:
- Move `build_genotype_matrix()` call into `engine.run_all()` gene loop
- Or: lazy-build via callback in gene_data dict
- Reduces peak memory from 16+ GB to ~1 GB

**Files:** `stages/analysis_stages.py`, `association/engine.py`

---

## 4. Impact Assessment

For 5K-gene GCKD cohort (5,125 samples):

| Fix | Current Time | After Fix | Speedup |
|-----|-------------|-----------|---------|
| Fix 1 (vectorize GT parse) | ~60 min | ~1 min | ~60x |
| Fix 2 (groupby) | ~30 min | ~0.01 min | ~3000x |
| Fix 3 (auto workers) | sequential COAST | parallel COAST | ~Nx (N=cores) |
| Fix 5 (no drop/recover) | 16 GB duplicate DF | 0 GB | memory |
| Fix 6 (stream matrices) | 16 GB all matrices | ~1 GB | memory |

**Combined estimate:** Association stage drops from 3+ hours to ~10-15 min
(with 8 workers: ~2-5 min for COAST execution + ~1 min for matrix construction).

---

## 5. Files to Modify

| File | Fixes | Changes |
|------|-------|---------|
| `association/genotype_matrix.py` | 1 | Vectorize `build_genotype_matrix()` |
| `stages/analysis_stages.py` | 2,3,5,6 | GroupBy, RM normalization, GT lifecycle |
| `pipeline_core/context.py` | 4 | Add `resource_manager` field |
| `pipeline.py` | 4 | Initialize shared ResourceManager |
| `association/base.py` | 3 | Change `association_workers` default to None |
| `cli.py` | 3 | Update help text for `--association-workers` |

## 6. Related

- Phase 27: Association Performance (GL quadrature + parallelization — done)
- Phase 8: DataFrame Optimization (itertuples migration — missed `genotype_matrix.py`)
- Issue #76: Pipeline performance analysis
- ResourceManager: `variantcentrifuge/memory/resource_manager.py` (well-designed, underused)
