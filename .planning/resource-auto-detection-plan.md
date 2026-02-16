# Plan: Unified Resource Auto-Detection Across Pipeline (IMPLEMENTED)

## Problem

The `ResourceManager` auto-detects CPUs and memory, but only the inheritance
analysis stages used it. The processing stages (bcftools, SnpSift, field
extraction) — which are the biggest bottlenecks — gated on `cfg["threads"]`
which defaulted to 1. Users got single-threaded processing unless they
explicitly passed `--threads N`.

Additionally, SnpSift was invoked with only 1GB JVM memory (conda wrapper
default), and the parallel processing stage used SnpSift extractFields which
crashed on chunks with missing VCF header fields.

## Goals

1. Auto-detect optimal thread count when user doesn't explicitly set `--threads`
2. Pass appropriate JVM memory to SnpSift (`-Xmx`)
3. Switch parallel field extraction from SnpSift to bcftools (faster, robust)
4. Add HPC environment variable detection (SLURM, PBS) to ResourceManager
5. Keep backward compatibility: explicit `--threads N` always honored

## Research Corrections (vs. original plan)

During implementation, web research and `--help` checks revealed:

- **`bcftools query` does NOT support `--threads`** — only `bcftools view`,
  `bcftools index`, etc. do. Removed from plan.
- **`SnpSift filter` does NOT have `-noStats`** — that flag doesn't exist.
  Removed from plan.
- **SnpSift conda wrapper** (`/home/bernt/miniforge3/bin/SnpSift`) is a Python
  script that intercepts `-Xm*` args and passes them to `java`. Default is
  `-Xmx1g`. The `-Xmx` flag must come before the subcommand (e.g.,
  `SnpSift -Xmx5g filter ...`).

## Implementation

### 1. CLI: `--threads` default → `"auto"` (`cli.py`)

Changed `--threads` from `type=int, default=1` to `default="auto"` (string).
Resolution happens once at CLI level so `cfg["threads"]` is always an int
by the time any stage sees it:

```python
if str(args.threads).lower() == "auto":
    from .memory.resource_manager import ResourceManager
    rm = ResourceManager(config=cfg)
    resolved_threads = rm.cpu_cores
    logger.info(f"Auto-detected {resolved_threads} CPU cores for --threads")
    cfg["threads"] = resolved_threads
    args.threads = resolved_threads
else:
    cfg["threads"] = int(args.threads)
    args.threads = int(args.threads)
```

Both `cfg["threads"]` and `args.threads` are set to the resolved int, so
`build_pipeline_stages(args)` and all downstream consumers work without changes.

Help text: `"Use 'auto' to detect available CPU cores (default: auto). Set to 1
to disable parallelization."`

### 2. ResourceManager: HPC CPU detection (`resource_manager.py`)

Enhanced `_detect_cpus()` with priority chain:

1. `SLURM_CPUS_PER_TASK` environment variable
2. `PBS_NUM_PPN` environment variable
3. `psutil.cpu_count(logical=False)` — physical cores
4. `os.cpu_count()` — logical cores
5. Fallback: 4 cores

Added `_get_cpu_source()` for log provenance (e.g., "SLURM", "psutil (physical)").
Updated `__init__` log message and `get_summary()` to include CPU source.

### 3. SnpSift JVM memory tuning (`filters.py`, `extractor.py`)

Added `_snpsift_memory_flag(cfg)` helper in `filters.py`:

```python
def _snpsift_memory_flag(cfg: dict[str, Any]) -> str:
    try:
        from .memory.resource_manager import ResourceManager
        rm = ResourceManager(config=cfg)
        memory_gb = max(2, min(16, int(rm.memory_gb * 0.25)))
    except Exception:
        memory_gb = 4
    return f"-Xmx{memory_gb}g"
```

Applied to:
- `SnpSift filter` command in `filters.py`
- `SnpSift extractFields` command in `extractor.py`

On a 20GB system: `-Xmx5g` (20 * 0.25 = 5), up from the default 1GB.

### 4. Parallel stage: SnpSift → bcftools extraction (`processing_stages.py`)

**Root cause of chunk failures:** `ParallelCompleteProcessingStage._process_single_chunk`
used `extract_fields` (aliased to `extract_fields_snpsift`). When a BED chunk's
filtered VCF lacked certain annotation INFO fields in its header (e.g.,
`dbNSFP_hg19_chr`), SnpSift extractFields crashed with:
`RuntimeException: INFO field 'dbNSFP_hg19_chr' not found in VCF header`

The single-threaded `FieldExtractionStage` already used `extract_fields_bcftools`
which handles missing fields gracefully (outputs `.`/NA). The parallel stage
was inconsistent.

**Fix:**
- Switched `_process_single_chunk` from `extract_fields` to `extract_fields_bcftools`
- Added `vcf_samples` to `worker_config` for proper per-sample column naming
- Removed unused `extract_fields` import

**Empty chunk guard** (`extractor.py`): When `bcftools query` produces empty
output (0 variants after filtering), `extract_fields_bcftools` now writes a
header-only TSV instead of crashing on `pd.read_csv` with an empty file.

### 5. Simplified inheritance auto-detection (`analysis_stages.py`)

**Before:** A heuristic branch used `threads_config == 1` as a proxy for
"user didn't set `--threads`" to decide whether to auto-detect workers.
This was ambiguous (`--threads 1` explicit vs. default) and skipped memory
constraints in the explicit-threads path.

**After:** Single path, no branching:
```python
threads_config = context.config.get("threads", 1)
memory_per_chunk_gb = rm.estimate_memory(chunk_size, num_samples)
memory_capped_workers = rm.auto_workers(
    task_count=max(1, threads_config),
    memory_per_task_gb=memory_per_chunk_gb,
)
max_workers = min(threads_config, memory_capped_workers)
```

`min(threads_config, memory_capped_workers)` ensures:
- Never exceeds user-specified or auto-detected thread count
- Never exceeds memory-safe worker count (OOM guard always active)
- `--threads 1` → `min(1, anything) = 1` → sequential path

**Note:** The gene-boundary-safe chunking in `_read_tsv_in_gene_chunks` is
unchanged — all variants for a gene stay in the same chunk, so compound
heterozygous detection remains safe.

## Testing

### New tests added

- `tests/unit/test_auto_threads.py` (6 tests):
  - SnpSift memory flag: known memory, low clamp, high clamp, auto-detection
  - `--threads auto` resolves to int
  - SLURM env var respected

- `tests/unit/memory/test_resource_manager.py` (7 new tests):
  - CPU detection: SLURM, PBS, invalid values, priority, fallback
  - `get_summary()` includes `cpu_source`

### Existing tests updated

- `tests/test_cli.py`: 5 tests updated for string-type `args.threads`
- `tests/unit/stages/test_parallel_complete_processing.py`: Mock target
  `extract_fields` → `extract_fields_bcftools`
- `tests/unit/test_compression_fix.py`: Same mock target update
- `tests/unit/test_parallel_stage_resume.py`: Added `vcf_samples` to mock context

### Verification

- **1152 tests pass**, CI clean (lint, format, typecheck, test-fast)
- **Subsample test** with 6-gene list on ~5000-sample VCF:
  - Auto-detected 8 CPUs, 20.2GB memory
  - 4 parallel chunks, `bcftools view --threads 2` per worker
  - `SnpSift -Xmx5g` for filter and extractFields
  - Empty chunk (0 variants) handled gracefully
  - All chunks succeeded, 100 variants output in 37.7s

## File Summary

| File | Change |
|------|--------|
| `cli.py` | `--threads` default → `"auto"`, resolve via ResourceManager |
| `resource_manager.py` | SLURM/PBS CPU detection, `_get_cpu_source()` |
| `filters.py` | `_snpsift_memory_flag()` helper, `-Xmx` on SnpSift filter |
| `extractor.py` | `-Xmx` on SnpSift extractFields, empty output guard in bcftools path |
| `processing_stages.py` | Parallel extraction: SnpSift → bcftools, pass `vcf_samples` |
| `analysis_stages.py` | Remove `threads==1` heuristic, always apply memory cap |
| `tests/` | 13 new tests, 3 test files updated for mock targets |
