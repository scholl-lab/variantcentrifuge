---
phase: 32-region-restriction-and-pca-wiring
verified: 2026-02-23T21:05:03Z
status: passed
score: 4/4 must-haves verified
re_verification: false
---

# Phase 32: Region Restriction and PCA Wiring Verification Report

**Phase Goal:** Users can restrict variant extraction to callable regions via a BED file, and PCA covariates computed by the pipeline are automatically available to association tests without manual file passing.
**Verified:** 2026-02-23T21:05:03Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Passing `--regions-bed capture_kit.bed` restricts bcftools variant extraction to those regions | VERIFIED | `--regions-bed` CLI flag exists in `cli.py` (line 128), mapped to `cfg["regions_bed"]` (line 1247); `GeneBedCreationStage._intersect_with_restriction_bed()` calls `bedtools intersect -a gene_bed -b restriction_bed` (line 247) and replaces `context.gene_bed_file` with the intersected BED (line 138–141); all downstream extraction stages use `context.gene_bed_file` |
| 2 | Chromosome naming mismatch triggers a clear error message, not silent 0-variant output | VERIFIED | `_intersect_with_restriction_bed()` reads chromosomes from both BEDs, compares `chr`-prefix presence, and raises `ValueError("Chromosome naming mismatch between gene BED and restriction BED: ...")` (lines 225–236) before bedtools runs |
| 3 | When `--pca akt` is set, PCA eigenvectors are computed automatically and appear as covariates in association output without extra steps | VERIFIED | `PCAComputationStage` in `processing_stages.py` (line 1600) runs `akt pca` via subprocess (line 1666), stores `pca_file` in `context.config` (line 1643) and marks stage complete with result dict; `AssociationAnalysisStage` reads result via `context.get_result("pca_computation")` and calls `context.config.setdefault("pca_file", ...)` (lines 2289–2295), then `_build_assoc_config_from_context` picks it up |
| 4 | The restriction BED intersection happens once in `GeneBedCreationStage`, upstream of chunk splitting | VERIFIED | Pipeline appends `GeneBedCreationStage()` at line 181, `ParallelCompleteProcessingStage` / `VariantExtractionStage` at lines 188–191, `ChunkedAnalysisStage` at line 224. The intersection logic is inside `GeneBedCreationStage._process()` after `context.gene_bed_file` is assigned (lines 135–141); all downstream stages receive the already-intersected BED via `context.gene_bed_file` |

**Score:** 4/4 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `variantcentrifuge/cli.py` | `--regions-bed` CLI argument | VERIFIED | Flag at line 128, file-existence validation at line 1482, mapped to `cfg["regions_bed"]` at line 1247 |
| `variantcentrifuge/cli.py` | Unified `--pca` flag | VERIFIED | `--pca` flag with backward-compat `--pca-file`/`--pca-tool` aliases using `dest="pca"`; `cfg["pca"]` at line 1249 |
| `variantcentrifuge/stages/processing_stages.py` | BED intersection in `GeneBedCreationStage` | VERIFIED | `_read_chromosomes()` (line 145), `_count_regions()` (line 169), `_intersect_with_restriction_bed()` (line 191), called in `_process()` (lines 135–141); 260+ lines, fully substantive |
| `variantcentrifuge/stages/processing_stages.py` | `PCAComputationStage` class | VERIFIED | `class PCAComputationStage(Stage)` at line 1600 with `_run_akt` method (line 1647); 82 lines, fully substantive |
| `variantcentrifuge/pipeline.py` | `PCAComputationStage` wired before `AssociationAnalysisStage` | VERIFIED | Imported at line 50, appended at line 289 before `AssociationAnalysisStage` at line 292; `create_stages_from_config` wires it at line 370 |
| `variantcentrifuge/stages/analysis_stages.py` | `pca_computation` soft dependency and result pickup | VERIFIED | `soft_dependencies` returns `{"custom_annotation", "pca_computation"}` (line 2227); `pca_result = context.get_result("pca_computation")` at line 2290; `"pca"` in `VALID_ASSOCIATION_KEYS` at line 1959 |
| `tests/unit/stages/test_processing_stages.py` | Unit tests for region restriction | VERIFIED | `TestGeneBedCreationRegionRestriction` class with 7 tests: `test_read_chromosomes`, `test_read_chromosomes_empty_file`, `test_count_regions`, `test_chr_mismatch_detection`, `test_empty_intersection_raises`, `test_successful_intersection`, `test_missing_restriction_bed_raises` |
| `tests/unit/stages/test_processing_stages.py` | Unit tests for PCA stage | VERIFIED | `TestPCAComputationStage` class with 9 tests covering: no-arg skip, file path, akt tool, cache reuse, invalid arg, AKT not found, AKT failure, no VCF, stage properties |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `cli.py` | `processing_stages.py (GeneBedCreationStage)` | `cfg["regions_bed"]` config key | VERIFIED | `cfg["regions_bed"] = getattr(args, "regions_bed", None)` (cli.py line 1247); `context.config.get("regions_bed")` (processing_stages.py line 136) |
| `processing_stages.py (GeneBedCreationStage)` | `context.gene_bed_file` | Overwrites with intersected BED | VERIFIED | `context.gene_bed_file = self._intersect_with_restriction_bed(...)` (line 138) |
| `cli.py` | `processing_stages.py (PCAComputationStage)` | `cfg["pca"]` config key | VERIFIED | `cfg["pca"] = getattr(args, "pca", None)` (cli.py line 1249); `context.config.get("pca")` (processing_stages.py line 1625) |
| `processing_stages.py (PCAComputationStage)` | `analysis_stages.py (AssociationAnalysisStage)` | `context.stage_results["pca_computation"]["pca_file"]` | VERIFIED | `context.mark_complete(self.name, result={"pca_file": pca_file})` (line 1644); `context.get_result("pca_computation")` (analysis_stages.py line 2290); `context.config.setdefault("pca_file", pca_file_from_stage)` (line 2294) |
| `pipeline.py` | `processing_stages.py (PCAComputationStage)` | Import and append in `build_pipeline_stages` | VERIFIED | Imported at line 50; appended at line 289 when `pca_arg` is set; `AssociationAnalysisStage` appended at line 292 (after PCA stage) |

### Requirements Coverage

No REQUIREMENTS.md entries mapped explicitly to phase 32 were found. The 4 observable truths derived from the stated phase goal are all verified.

### Anti-Patterns Found

No blockers or warnings found. All key methods have real implementations. No TODO/FIXME/placeholder patterns present in the phase-modified files.

### Human Verification Required

None. All behaviors verified structurally:
- `--regions-bed` restriction is enforced at the `context.gene_bed_file` level which all extraction stages consume
- Chromosome mismatch detection logic is clear and tested
- PCAComputationStage subprocess wiring is straightforward and unit-tested with mocks
- Stage ordering in `pipeline.py` is unambiguous (line numbers confirm PCA before Association)

### Gaps Summary

No gaps. All 4 must-have truths are verified with complete artifact, substantive, and wiring checks. All 16 new unit tests (7 region restriction + 9 PCA) pass. Full unit test suite of 1570 tests passes with no regressions.

---
_Verified: 2026-02-23T21:05:03Z_
_Verifier: Claude (gsd-verifier)_
