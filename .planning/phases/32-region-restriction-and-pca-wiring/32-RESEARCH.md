# Phase 32: Region Restriction and PCA Wiring - Research

**Researched:** 2026-02-23
**Domain:** Pipeline stage extension, BED file intersection, PCA subprocess wiring
**Confidence:** HIGH

## Summary

Phase 32 adds two independent features: (1) a BED-based region restriction prefilter that
intersects the gene BED with a user-supplied capture kit BED before variant extraction, and (2)
wiring the existing `PCAComputationStage` stub (removed in Phase 30) back into the pipeline so
`--pca-tool akt` actually runs AKT via subprocess and feeds eigenvectors to `AssociationAnalysisStage`.

**Region restriction** modifies `GeneBedCreationStage` to intersect the generated gene BED with a
`--regions-bed` restriction BED via `bedtools intersect`. The intersected BED is written back to
`context.gene_bed_file` so all downstream stages (single-threaded extraction, parallel chunk
splitting, bcftools `-R`) consume it without modification. Chromosome naming mismatches are a hard
fail with a clear error message.

**PCA wiring** adds a new `PCAComputationStage` to `processing_stages.py` that runs AKT via
subprocess when `--pca akt` is given. The resulting eigenvector file path is stored in
`context.stage_results["pca_computation"]` so `AssociationAnalysisStage` can load it via the
existing `load_pca_file` / `merge_pca_covariates` machinery. The `--pca` flag replaces the separate
`--pca-file` / `--pca-tool` flags, auto-detecting whether the argument is a tool name or a file
path.

**Primary recommendation:** Both features are surgical additions to existing stages. Region
restriction is a pre-step inside `GeneBedCreationStage._process()`. PCA computation is a new stage
inserted into `pipeline.py` before `AssociationAnalysisStage`.

## Standard Stack

### Core (already in use — no new dependencies)

| Library / Tool | Version | Purpose | Why Standard |
|----------------|---------|---------|--------------|
| `bedtools` | any in PATH | BED intersection | Already used in `gene_bed.py` (`bedtools merge`) |
| `subprocess` | stdlib | Run AKT | Used throughout pipeline for bcftools, snpEff |
| `bcftools` | any in PATH | VCF extraction with `-R` | Core extraction tool, already wired |
| `association.pca.load_pca_file` | internal | Load eigenvec files | Already handles PLINK + AKT formats |
| `association.pca.merge_pca_covariates` | internal | Merge PCs into covariates | Already tested and used by AssociationAnalysisStage |

### New External Tool Dependency

| Tool | How Invoked | Output Format |
|------|-------------|---------------|
| `akt pca` | `subprocess.run(["akt", "pca", ...])` | AKT stdout format (headerless TSV: sample\_id PC1 PC2 ...) — already handled by `load_pca_file` |

**Installation:** No new Python packages needed.

## Architecture Patterns

### Recommended Project Structure (unchanged)

The feature slots into the existing stage architecture without restructuring:

```
variantcentrifuge/
├── stages/
│   ├── processing_stages.py   # Add PCAComputationStage here (new class)
│   │                          # Modify GeneBedCreationStage._process() for region restriction
│   └── analysis_stages.py     # Modify AssociationAnalysisStage to read pca_computation result
├── pipeline.py                # Wire PCAComputationStage into build_pipeline_stages()
└── cli.py                     # Add --regions-bed, unify --pca (replace --pca-file/--pca-tool)
```

### Pattern 1: Region Restriction via bedtools intersect in GeneBedCreationStage

**What:** After `get_gene_bed()` produces `context.gene_bed_file`, run `bedtools intersect -a
gene.bed -b restriction.bed > intersected.bed`. Replace `context.gene_bed_file` with the result.
Perform chromosome naming mismatch detection before intersection.

**When to use:** `context.config.get("regions_bed")` is set.

**Why intersect here, not in VariantExtractionStage:** The requirement states "intersection happens
once upstream, before chunk splitting." `GeneBedCreationStage` runs before
`ParallelCompleteProcessingStage._split_bed_file()`, so inserting here guarantees all chunks use
the restricted BED.

**Example:**
```python
# Source: codebase inspection — pattern from gene_bed.py
def _intersect_with_restriction_bed(
    self,
    gene_bed: Path,
    restriction_bed: str,
    context: PipelineContext,
) -> Path:
    # 1. Validate restriction BED exists and is not empty
    if not os.path.exists(restriction_bed):
        raise FileNotFoundError(f"Restriction BED not found: {restriction_bed}")

    # 2. Detect chromosome naming mismatch BEFORE intersection
    gene_chroms = self._read_chromosomes(gene_bed)
    restrict_chroms = self._read_chromosomes(Path(restriction_bed))
    # Mismatch: one set has 'chr' prefix, other does not
    gene_has_chr = any(c.startswith("chr") for c in gene_chroms)
    restrict_has_chr = any(c.startswith("chr") for c in restrict_chroms)
    if gene_has_chr != restrict_has_chr:
        raise ValueError(
            f"Chromosome naming mismatch between gene BED (chr prefix: {gene_has_chr}) "
            f"and restriction BED (chr prefix: {restrict_has_chr}). "
            "Ensure both BED files use the same chromosome naming convention."
        )

    # 3. Run bedtools intersect
    output_path = context.workspace.get_intermediate_path(
        f"{context.workspace.base_name}.restricted.bed"
    )
    cmd = ["bedtools", "intersect", "-a", str(gene_bed), "-b", restriction_bed]
    with open(output_path, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, check=True)

    # 4. Validate output is not empty
    if output_path.stat().st_size == 0:
        raise ValueError(
            "Intersection of gene BED and restriction BED produced no regions. "
            "Check that the BED files overlap and have matching chromosome names."
        )

    # 5. Log summary
    orig_count = self._count_regions(gene_bed)
    new_count = self._count_regions(output_path)
    logger.info(
        f"Region restriction: {new_count} of {orig_count} gene regions retained "
        f"after intersection with {restriction_bed}"
    )

    return output_path
```

**BED validation — strict mode:** The CONTEXT.md requires failing on malformed BED (overlapping,
unsorted, missing columns). The simplest implementation: run `bedtools intersect` and let bedtools
raise on malformed input. Do not add a separate pre-validation loop. Bedtools error output will be
surfaced in the exception. This is sufficient for strict validation.

**Zero-variant genes:** This is a downstream concern. The restriction BED produces fewer regions,
but whether a region has zero variants after bcftools extraction is handled in existing stages.
The requirement says "log at DEBUG level" — this already happens implicitly (bcftools output is
empty, FieldExtractionStage produces empty TSV). No special handling needed.

### Pattern 2: PCA Computation Stage (new Stage class)

**What:** `PCAComputationStage` runs before `AssociationAnalysisStage`. It either:
- Loads a pre-computed eigenvector file (if `--pca` arg is a file path), or
- Runs `akt pca <vcf_file>` via subprocess and captures stdout to a temp file (if `--pca` arg is
  `"akt"`).

Result is stored in `context.stage_results["pca_computation"]` as a dict with key `"pca_file"`.
`AssociationAnalysisStage` checks this result to set `assoc_config.pca_file`.

**When to use:** `context.config.get("pca")` is set (new unified flag).

**parallel_safe:** `False` — AKT subprocess is not thread-safe.

**Example:**
```python
# Source: codebase inspection — pattern from analysis_stages.py subprocess usage
class PCAComputationStage(Stage):
    @property
    def name(self) -> str:
        return "pca_computation"

    @property
    def dependencies(self) -> set[str]:
        return {"sample_config_loading", "configuration_loading"}

    @property
    def parallel_safe(self) -> bool:
        return False  # subprocess, not thread-safe

    def _process(self, context: PipelineContext) -> PipelineContext:
        pca_arg = context.config.get("pca")
        if not pca_arg:
            return context

        n_components = context.config.get("pca_components", 10)

        # Auto-detect: file path vs tool name
        if os.path.isfile(pca_arg):
            # Pre-computed eigenvectors — just store path
            pca_file = pca_arg
            logger.info(f"PCA: using pre-computed eigenvectors from {pca_file}")
        elif pca_arg == "akt":
            pca_file = self._run_akt(context, n_components)
        else:
            raise ValueError(
                f"--pca argument '{pca_arg}' is neither a valid file path nor a "
                "recognized tool name ('akt')."
            )

        # Store result for AssociationAnalysisStage
        context.mark_complete(self.name, result={"pca_file": pca_file})
        return context

    def _run_akt(self, context: PipelineContext, n_components: int) -> str:
        vcf_file = context.config.get("vcf_file")
        if not vcf_file:
            raise ValueError("--pca akt requires a VCF file to be specified")

        output_path = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.pca.eigenvec"
        )

        # Check cache: if output_path exists and VCF hasn't changed, reuse
        if output_path.exists() and output_path.stat().st_size > 0:
            logger.info(f"PCA: reusing cached eigenvectors from {output_path}")
            return str(output_path)

        cmd = ["akt", "pca", vcf_file, "-o", str(output_path), "-n", str(n_components)]
        logger.info(f"Running AKT PCA: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.debug(f"AKT PCA stdout: {result.stdout[:200]}")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"AKT PCA failed (exit code {e.returncode}): {e.stderr}"
            ) from e
        except FileNotFoundError as e:
            raise RuntimeError(
                "AKT not found in PATH. Install akt or use --pca <eigenvec_file> "
                "to supply pre-computed eigenvectors."
            ) from e

        if not output_path.exists() or output_path.stat().st_size == 0:
            raise RuntimeError(
                "AKT PCA completed but produced no output. "
                "Check that the VCF has sufficient samples and variants."
            )

        return str(output_path)
```

**Note on AKT output format:** AKT `pca` outputs to a file (not stdout) when `-o` is given, or
stdout otherwise. The existing `load_pca_file` in `association/pca.py` already handles the AKT
format (`akt_or_generic`). The file format is: no header, first column = sample ID, remaining
columns = PC values. This is confirmed by `test_pca.py::TestDetectPcaFormat::test_akt_format`.

### Pattern 3: Wire PCA Stage into AssociationAnalysisStage

**What:** `AssociationAnalysisStage._process()` already loads `pca_file` from `assoc_config`. The
only change is: if `context.stage_results.get("pca_computation")` contains a `pca_file`, set it
on the config before building `assoc_config`.

**Where:** At the top of `AssociationAnalysisStage._process()`, before `_build_assoc_config_from_context()`.

```python
# In AssociationAnalysisStage._process():
# PCA-02: pick up pca_file from PCAComputationStage result
pca_result = context.get_result("pca_computation")
if pca_result and pca_result.get("pca_file"):
    # Override config so _build_assoc_config_from_context picks it up
    context.config["pca_file"] = pca_result["pca_file"]
    logger.debug(f"Association: using PCA file from pipeline stage: {pca_result['pca_file']}")
```

### Pattern 4: Unified --pca CLI flag

**What:** Replace `--pca-file` and `--pca-tool` with a single `--pca` flag in `cli.py`. The
heuristic: `os.path.exists(value)` → treat as file; `value in {"akt"}` → treat as tool name.

**Backward compatibility:** Both `--pca-file` and `--pca-tool` are currently in the CLI. Phase 32
replaces them with `--pca`. The CONTEXT.md says "Unify with any existing PCA-related flags (check
`--pca-tool` from current CLI)". The old flags can be deprecated (argparse hidden) or removed;
since the CONTEXT.md mentions unification, remove them and replace.

**Config key:** Store unified value as `cfg["pca"] = getattr(args, "pca", None)`. Leave `cfg["pca_file"]`
and `cfg["pca_tool"]` as None unless `PCAComputationStage` resolves them, to avoid breaking
existing JSON config paths that reference `pca_file`.

### Pattern 5: Pipeline Wiring in pipeline.py

`PCAComputationStage` must run before `AssociationAnalysisStage`. Insert it when `pca` is set:

```python
# In build_pipeline_stages():
if hasattr(args, "pca") and args.pca:
    from .stages.processing_stages import PCAComputationStage
    stages.append(PCAComputationStage())
```

The dependency graph handles ordering: `AssociationAnalysisStage` can add a soft dependency on
`pca_computation`, or `PCAComputationStage` can declare no dependency on association-related stages,
and the runner will naturally sequence it before association (which depends on `dataframe_loading`).

**Safest approach:** Add `"pca_computation"` as a soft dependency of `AssociationAnalysisStage`.
This is non-breaking (soft dependencies only apply when the stage is present).

### Anti-Patterns to Avoid

- **Storing eigenvectors in `context.current_dataframe` or any DataFrame field** — genotype matrix
  must never enter PipelineContext (memory constraint from architecture invariants). Only the FILE
  PATH is stored via `stage_results`.
- **Running AKT on a per-chunk basis** — PCA is a global operation on all samples. Run once before
  chunking.
- **Using `--pca` to take multiple values** — keep single-value interface matching CONTEXT.md.
- **Auto-correcting chr-prefix mismatches** — CONTEXT.md explicitly requires hard fail with clear
  error message.
- **Writing intersected BED to a new `context` field** — overwrite `context.gene_bed_file` directly
  so all downstream consumers (VariantExtractionStage, ParallelCompleteProcessingStage) get the
  restricted version without modification.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| BED intersection | Custom Python interval logic | `bedtools intersect` | Already in PATH (bedtools already used in gene_bed.py), handles edge cases |
| Chromosome detection | Custom chr-prefix parser | Read first column of each BED file | Simple string `startswith("chr")` check is sufficient |
| PCA file format detection | Custom parser | `association.pca.load_pca_file` | Already handles PLINK (with/without header) + AKT formats, tested |
| PCA eigenvector loading | Custom numpy parsing | `association.pca.load_pca_file` | Sample alignment, missing sample errors, n_components slicing all handled |
| Merging PCs with covariates | Custom hstack | `association.pca.merge_pca_covariates` | Handles None covariate matrix case |
| Cache key for PCA | Custom hash | File path + mtime check | Simple: if output file exists and non-empty, reuse it |

**Key insight:** The `association/pca.py` module is complete and well-tested. Phase 32's PCA work
is purely wiring: subprocess invocation of AKT + storing the output path so the existing machinery
can load it.

## Common Pitfalls

### Pitfall 1: bedtools intersect -a/-b ordering matters

**What goes wrong:** `bedtools intersect -a restriction.bed -b gene.bed` vs `-a gene.bed -b restriction.bed`
produce the same output regions, but the output takes the coordinates from the `-a` file.
**Why it happens:** bedtools intersect reports intervals from the `-a` file that overlap with `-b`.
**How to avoid:** Always use `-a gene.bed -b restriction.bed`. This ensures the output retains gene
BED coordinates (which may have been expanded by `interval_expand`).
**Warning signs:** Output regions are misaligned with expected gene intervals.

### Pitfall 2: bedtools intersect produces empty output without error

**What goes wrong:** If no regions overlap (e.g., chr-prefix mismatch or truly non-overlapping
BEDs), bedtools returns exit code 0 with empty output.
**Why it happens:** Empty intersection is a valid result for bedtools.
**How to avoid:** Check output file size after intersection. If empty, raise a clear ValueError
with mismatch details. The chromosome mismatch detection (reading first column of each BED)
should catch this before the intersect call.
**Warning signs:** Pipeline produces 0 variants in output without any error.

### Pitfall 3: AKT writes to file only when -o is given, otherwise stdout

**What goes wrong:** Running `akt pca vcf.vcf.gz` without `-o` writes to stdout. If not captured
correctly, the eigenvectors are lost.
**Why it happens:** AKT's default behavior.
**How to avoid:** Always pass `-o <output_path>` to akt. Use `subprocess.run(cmd, check=True)` to
detect failures; do not use `capture_output=True` when AKT writes to file via `-o`.
**Warning signs:** Output file is empty or missing after AKT command succeeds.

### Pitfall 4: Unifying --pca breaks existing --pca-file users

**What goes wrong:** Removing `--pca-file` and `--pca-tool` breaks users who have scripts using
those flags.
**Why it happens:** Breaking change without deprecation path.
**How to avoid:** Either (a) keep old flags as hidden aliases that set `cfg["pca"]` to their value,
or (b) clearly document the flag rename in CHANGELOG. Given this is v0.16.0, option (b) is
acceptable. The CONTEXT.md says "unify", implying replacement not alias.
**Warning signs:** Tests using `--pca-file` or `--pca-tool` break — update those tests.

### Pitfall 5: PCAComputationStage soft dependency ordering

**What goes wrong:** If `PCAComputationStage` is appended to the stage list AFTER
`AssociationAnalysisStage`, the runner may execute association before PCA if it groups them
incorrectly.
**Why it happens:** The runner uses topological sort on declared dependencies. If soft dependency
is not declared, ordering is not guaranteed.
**How to avoid:** Declare `pca_computation` as a soft dependency of `AssociationAnalysisStage`
OR ensure `PCAComputationStage` is added to the list BEFORE `AssociationAnalysisStage` in
`build_pipeline_stages()`. Both work; soft dependency is more explicit.
**Warning signs:** AssociationAnalysisStage runs before PCA file is created.

### Pitfall 6: Parallel chunk processing uses gene_bed_file after restriction

**What goes wrong:** If region restriction modifies `context.gene_bed_file`, parallel chunks
must use the restricted BED. But `ParallelCompleteProcessingStage._split_bed_file()` reads from
`context.gene_bed_file`, so this works automatically IF restriction happens inside
`GeneBedCreationStage` (which runs before parallel processing).
**Why it happens:** Ordering is correct because of declared dependencies.
**How to avoid:** Verify that `ParallelCompleteProcessingStage` lists `gene_bed_creation` in its
dependencies (it does: `return {"gene_bed_creation", "configuration_loading"}`). No additional
work needed.

## Code Examples

### Reading chromosomes from BED file for mismatch detection

```python
# Source: codebase inspection, derived from split_bed_file in utils.py
def _read_chromosomes(self, bed_file: Path) -> set[str]:
    """Return the set of chromosome names (first column) from a BED file."""
    chroms = set()
    with open(bed_file) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            parts = line.split("\t", 1)
            if parts:
                chroms.add(parts[0].strip())
    return chroms
```

### BED region count for INFO logging

```python
def _count_regions(self, bed_file: Path) -> int:
    """Count non-comment lines in a BED file."""
    count = 0
    with open(bed_file) as f:
        for line in f:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                count += 1
    return count
```

### Reading pca_computation result in AssociationAnalysisStage

```python
# Source: context.py - get_result() is thread-safe
# In AssociationAnalysisStage._process() before _build_assoc_config_from_context():
pca_stage_result = context.get_result("pca_computation")
if pca_stage_result and isinstance(pca_stage_result, dict):
    pca_file_from_stage = pca_stage_result.get("pca_file")
    if pca_file_from_stage:
        context.config.setdefault("pca_file", pca_file_from_stage)
        logger.info(f"Association: PCA file from pipeline stage: {pca_file_from_stage}")
```

### AKT invocation pattern

AKT PCA command signature (from AKT docs / bioinformatics community):
```bash
akt pca input.vcf.gz -o output.eigenvec -n 10
```

Where:
- `-o` output file (AKT writes eigenvectors here)
- `-n` number of PCs (matches `--pca-components` / `pca_components`)

The output format matches `akt_or_generic` in `load_pca_file`: no header, tab/space separated,
first column = sample ID, columns 2-N = PC values.

### CLI --pca flag auto-detection

```python
# In build_assoc_config or wherever pca arg is processed:
pca_arg = getattr(args, "pca", None)
if pca_arg:
    if os.path.isfile(pca_arg):
        cfg["pca_file"] = pca_arg       # pre-computed file
        cfg["pca_tool"] = None
    elif pca_arg == "akt":
        cfg["pca_tool"] = "akt"          # tool to invoke
        cfg["pca_file"] = None
    else:
        parser.error(
            f"--pca: '{pca_arg}' is not a valid file path or recognized tool. "
            "Use a file path to pre-computed eigenvectors or 'akt' to compute."
        )
```

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|------------------|--------|
| `--pca-file` + `--pca-tool` (two flags) | `--pca` (single smart flag) | Simpler UX; auto-detects file vs tool |
| PCAComputationStage removed (Phase 30) | Re-implement as proper Stage | AKT wiring now complete |
| No region restriction | `--regions-bed` intersection in GeneBedCreationStage | Global prefilter before chunking |

**Deprecated/outdated:**
- `--pca-file` and `--pca-tool` CLI flags: replaced by unified `--pca`
- Any reference to `PCAComputationStage` from Phase 30 removal notes: Phase 32 brings it back

## Open Questions

1. **AKT -o flag behavior**
   - What we know: AKT is a bioinformatics tool for PCA from WGS data; `-o` is standard for
     output file
   - What's unclear: Exact flags and whether AKT accepts gzipped VCF or requires indexed VCF
   - Recommendation: Test with `akt pca --help` in the environment. Assume `-o` writes eigenvec
     file, `-n` sets PC count. If AKT requires indexed VCF, add indexing check.

2. **Cache invalidation for PCA eigenvectors**
   - What we know: Gene BED uses file-content hash (md5 of gene list + params) for caching
   - What's unclear: Should PCA cache invalidate when VCF changes? Simple mtime check may be
     fragile for large VCFs
   - Recommendation: For Phase 32, use the simple "output file exists and non-empty" check.
     Full cache invalidation is future work. Log a DEBUG message when cache is reused.

3. **VALID_ASSOCIATION_KEYS update**
   - What we know: `VALID_ASSOCIATION_KEYS` in `analysis_stages.py` lists valid JSON config keys
   - What's unclear: Should `"pca"` (unified flag) be added to replace `"pca_file"` / `"pca_tool"`?
   - Recommendation: Keep `"pca_file"` and `"pca_tool"` in `VALID_ASSOCIATION_KEYS` for JSON
     config backward compatibility. Add `"pca"` as a new valid key. The stage resolves either to
     the same internal logic.

## Sources

### Primary (HIGH confidence)

- Codebase: `variantcentrifuge/association/pca.py` — complete PCA file loading logic confirmed
- Codebase: `variantcentrifuge/association/covariates.py` — covariate loading pattern for PCA merge
- Codebase: `variantcentrifuge/stages/processing_stages.py` — GeneBedCreationStage and
  ParallelCompleteProcessingStage structure confirmed
- Codebase: `variantcentrifuge/stages/analysis_stages.py` — AssociationAnalysisStage PCA wiring
  point confirmed (lines 2443-2463)
- Codebase: `variantcentrifuge/pipeline_core/context.py` — `stage_results` dict + `get_result()` +
  `mark_complete()` confirmed as the correct inter-stage communication mechanism
- Codebase: `variantcentrifuge/pipeline.py` — `build_pipeline_stages()` pattern confirmed for
  conditional stage insertion
- Codebase: `variantcentrifuge/gene_bed.py` — `bedtools merge` already used, confirming bedtools
  is in PATH and subprocess pattern is established
- Codebase: `variantcentrifuge/cli.py` — existing `--pca-file` / `--pca-tool` flags confirmed;
  `--regions-bed` does not yet exist
- Codebase: `variantcentrifuge/filters.py` — `extract_variants()` uses `bcftools view -R bed_file`
  confirming BED is passed directly to bcftools

### Secondary (MEDIUM confidence)

- `tests/unit/test_pca.py` — AKT output format (`akt_or_generic`) confirmed via test fixtures
- AKT documentation (general bioinformatics knowledge) — `akt pca -o output -n N input.vcf.gz`
  is the standard command pattern

### Tertiary (LOW confidence — verify in environment)

- AKT exact CLI flags: verify `akt pca --help` before implementing
- Whether AKT requires tabix-indexed VCF: verify in test environment

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all tools already in use in codebase
- Architecture: HIGH — existing patterns confirmed by direct code inspection
- Pitfalls: HIGH — derived from direct code analysis + BED intersection semantics
- AKT subprocess details: MEDIUM — command pattern from bioinformatics community knowledge

**Research date:** 2026-02-23
**Valid until:** 2026-03-23 (stable domain — bedtools and AKT APIs are stable)
