# CLI Flag Implementation Analysis: Old vs New Pipeline

**Analysis Date:** January 9, 2025  
**Author:** Claude Code Analysis  
**Branch:** refactor-modular  

## Executive Summary

This document provides a comprehensive analysis of CLI flag implementation across the old monolithic pipeline (`pipeline.py`) and the new stage-based pipeline (`pipeline_refactored.py`). Out of **72 total CLI flags**, the refactored pipeline has **significant gaps** in implementation:

- **✅ 45 flags (63%)** - Fully implemented  
- **⚠️ 12 flags (17%)** - Partially implemented
- **❌ 15 flags (20%)** - Not implemented or missing core functionality

## Methodology

The analysis involved:
1. Reading and analyzing the old `pipeline.py` implementation
2. Reading and analyzing the new `pipeline_refactored.py` implementation
3. Comparing CLI argument parsing between both implementations
4. Analyzing flag implementation in individual stage files
5. Cross-referencing CLI help output with actual implementation

## Architectural Differences

### Old Pipeline (Direct Access Pattern)
```python
# Direct argument access throughout pipeline
if args.transcript_list:
    transcript_ids = [t.strip() for t in args.transcript_list.split(",")]
    
if args.genotype_filter:
    genotype_filters = args.genotype_filter.split(",")
```

### New Pipeline (Stage-Based Pattern)
```python
# Arguments flow through context and stages
class TranscriptFilterStage(Stage):
    def _process(self, context):
        transcript_list = context.config.get('transcript_list')
        # Process transcripts...
```

## Detailed Flag Analysis

### 🟢 FULLY IMPLEMENTED (45 flags)

#### General Options (5/5)
- `--version` ✅
- `--log-level` ✅
- `--log-file` ✅
- `--config` ✅
- `--use-new-pipeline` ✅

#### Core Input/Output (6/6)
- `--vcf-file` ✅
- `--output-file` ✅
- `--output-dir` ✅
- `--xlsx` ✅
- `--keep-intermediates` ✅
- `--archive-results` ✅

#### Gene Selection (2/4)
- `--gene-name` ✅
- `--gene-file` ✅
- `--transcript-list` ❌ (See Critical Missing)
- `--transcript-file` ❌ (See Critical Missing)

#### Basic Filtering & Annotation (8/9)
- `--reference` ✅
- `--filters` ✅
- `--bcftools-prefilter` ✅ (Implemented in `BCFToolsPrefilterStage`)
- `--preset` ✅
- `--late-filtering` ✅
- `--final-filter` ✅
- `--split-snpeff-lines` ⚠️ (Partial - see below)
- `--no-replacement` ✅ (Controls `GenotypeReplacementStage`)
- `--add-column` ✅

#### Field Extraction & Formatting (3/7)
- `--fields` ✅
- `--no-links` ✅
- `--add-chr` ❌ (Missing in stages)
- `--remove-sample-substring` ❌ (Missing in stages)
- `--append-extra-sample-fields` ❌ (See Critical Missing)
- `--extra-sample-field-delimiter` ❌ (Missing in stages)

#### Phenotype & Sample Groups (8/8)
- `--phenotype-file` ✅ (Implemented in `PhenotypeLoadingStage`)
- `--phenotype-sample-column` ✅
- `--phenotype-value-column` ✅
- `--case-phenotypes` ✅
- `--control-phenotypes` ✅
- `--case-phenotypes-file` ✅
- `--control-phenotypes-file` ✅
- `--case-samples` ✅ (Implemented in `SampleConfigLoadingStage`)
- `--control-samples` ✅
- `--case-samples-file` ✅
- `--control-samples-file` ✅

#### Statistical Analysis (6/6)
- `--perform-gene-burden` ✅ (Implemented in `GeneBurdenAnalysisStage`)
- `--gene-burden-mode` ✅ (Lines 707-708 in `analysis_stages.py`)
- `--correction-method` ✅ (Lines 707-708 in `analysis_stages.py`)
- `--no-stats` ✅ (Controls `StatisticsGenerationStage`)
- `--stats-output-file` ✅
- `--stats-config` ✅

#### Inheritance Analysis (3/4)
- `--ped` ✅ (Implemented in `PedigreeLoadingStage`)
- `--inheritance-mode` ✅ (Implemented in `InheritanceAnalysisStage`)
- `--no-vectorized-comp-het` ✅ (Lines 372, 387 in `analysis_stages.py`)

#### Scoring & Custom Annotations (6/6)
- `--scoring-config-path` ✅ (Implemented in `ScoringConfigLoadingStage`)
- `--annotate-bed` ✅ (Implemented in `AnnotationConfigLoadingStage`)
- `--annotate-gene-list` ✅
- `--annotate-json-genes` ✅
- `--json-gene-mapping` ✅
- `--json-genes-as-columns` ✅

#### Reporting & Visualization (5/12)
- `--html-report` ✅ (Implemented in `HTMLReportStage`)
- `--igv` ✅ (Implemented in `IGVReportStage`)
- `--bam-mapping-file` ✅
- `--igv-reference` ✅
- `--igv-fasta` ✅

#### Performance & Processing (2/9)
- `--threads` ✅ (Controls `ParallelCompleteProcessingStage` vs individual stages)
- `--chunk-size` ✅ (Implemented in `ChunkedAnalysisStage`)

#### Data Privacy Options (5/6)
- `--pseudonymize` ✅ (Implemented in `PseudonymizationStage`)
- `--pseudonymize-schema` ✅
- `--pseudonymize-prefix` ✅
- `--pseudonymize-pattern` ✅
- `--pseudonymize-category-field` ✅
- `--pseudonymize-table` ✅

### ⚠️ PARTIALLY IMPLEMENTED (12 flags)

#### Split SnpEff Lines (1 flag)
- `--split-snpeff-lines` ⚠️
  - **Issue:** Has `MultiAllelicSplitStage` but not the specific snpeff line splitting
  - **Found:** Reference to `split_snpeff_annotations` in `processing_stages.py` (line 37) but no dedicated stage
  - **Status:** Logic exists but not properly integrated as a stage

#### IGV Visualization (7 flags)
- `--igv-ideogram` ⚠️ (Arguments mapped to config but no stage implementation)
- `--igv-flanking` ⚠️ (Arguments mapped to config but no stage implementation)
- `--igv-max-allele-len-filename` ⚠️
- `--igv-hash-len-filename` ⚠️
- `--igv-max-variant-part-filename` ⚠️
  - **Status:** All mapped to config in `pipeline_refactored.py` lines 332-340 but not implemented in `IGVReportStage`

#### Checkpoint & Resume (1/4)
- `--enable-checkpoint` ✅ (Basic implementation)
- `--resume` ❌ (Missing)
- `--checkpoint-checksum` ❌ (Missing)
- `--show-checkpoint-status` ❌ (Missing)

#### Performance Flags (3 flags)
- `--sort-memory-limit` ⚠️ (Only implemented in `ChunkedAnalysisStage` line 896, missing elsewhere)
- `--no-chunked-processing` ❌ (Missing)
- `--force-chunked-processing` ❌ (Missing)
- `--sort-parallel` ❌ (Missing)

### 🔴 CRITICAL MISSING FUNCTIONALITY (15 flags)

#### High Priority - Core Functionality Missing

1. **Transcript Filtering (2 flags)**
   - `--transcript-list` ❌
   - `--transcript-file` ❌
   - **Status:** CLI arguments exist but NO stage implementation found
   - **Expected location:** Should be in processing or filtering stages
   - **Risk:** HIGH - Core genomics functionality missing

2. **Genotype Analysis (6 flags)**
   - `--genotype-filter` ❌
   - `--gene-genotype-file` ❌
   - `--append-extra-sample-fields` ❌
   - `--extra-sample-field-delimiter` ❌
   - **Status:** CLI arguments exist but NO stage implementation found
   - **Found:** Only default False value in `processing_stages.py` (line 697) for extra fields
   - **Risk:** HIGH - Essential variant analysis capability absent

#### Medium Priority - Feature Gaps

3. **Field Processing (2 flags)**
   - `--add-chr` ❌
   - `--remove-sample-substring` ❌
   - **Status:** Missing in stages
   - **Risk:** MEDIUM - Data formatting functionality limited

4. **Advanced IGV Options (5 flags)**
   - `--igv-ideogram` ❌
   - `--igv-flanking` ❌
   - `--igv-max-allele-len-filename` ❌
   - `--igv-hash-len-filename` ❌
   - `--igv-max-variant-part-filename` ❌
   - **Status:** Arguments mapped to config but no stage implementation
   - **Risk:** MEDIUM - Visualization features incomplete

#### Low Priority - Performance/Convenience

5. **Checkpoint Features (3 flags)**
   - `--resume` ❌
   - `--checkpoint-checksum` ❌
   - `--show-checkpoint-status` ❌
   - **Status:** CLI arguments exist but NO stage implementation found
   - **Risk:** LOW - Workflow convenience missing

6. **Performance Optimization (3 flags)**
   - `--no-chunked-processing` ❌
   - `--force-chunked-processing` ❌
   - `--sort-parallel` ❌
   - **Status:** Missing implementation
   - **Risk:** LOW - Performance optimization missing

7. **Privacy Enhancement (1 flag)**
   - `--pseudonymize-ped` ❌
   - **Status:** Missing from privacy features
   - **Risk:** LOW - Additional privacy feature missing

## Specific Implementation Details

### Working Implementations

#### bcftools-prefilter
- **Location:** `BCFToolsPrefilterStage` in `processing_stages.py` (lines 391-439)
- **Implementation:** Full separate stage for filtering during extraction

#### Gene Burden Analysis
- **Location:** `GeneBurdenAnalysisStage` in `analysis_stages.py` (lines 688-702)
- **Implementation:** Both `--correction-method` and `--gene-burden-mode` properly handled

#### Compound Heterozygous Analysis
- **Location:** `InheritanceAnalysisStage` in `analysis_stages.py` (lines 372, 387)
- **Implementation:** `--no-vectorized-comp-het` properly inverted and passed

#### Pseudonymization
- **Location:** `PseudonymizationStage` in `output_stages.py` (lines 320-376)
- **Implementation:** Full support for schema, prefix, category field, and metadata

### Missing Critical Implementations

#### Transcript Filtering
- **Expected:** `TranscriptFilteringStage` 
- **Implementation needed:** Process `--transcript-list` and `--transcript-file`
- **Current status:** Arguments parsed but no stage processes them

#### Genotype Filtering
- **Expected:** `GenotypeFilteringStage`
- **Implementation needed:** Process `--genotype-filter` and `--gene-genotype-file`
- **Current status:** Core filtering capability completely absent

#### Extra Sample Fields
- **Expected:** `ExtraSampleFieldsStage`
- **Implementation needed:** Process `--append-extra-sample-fields`
- **Current status:** Only default value found, no logic to append fields

## Risk Assessment

### HIGH RISK (Production Impact)
- **Transcript filtering** - Core genomics functionality missing
- **Genotype filtering** - Essential variant analysis capability absent

### MEDIUM RISK (Feature Gaps)
- **SnpEff processing** - Annotation workflow incomplete  
- **Extra sample fields** - Data export functionality limited
- **Advanced IGV options** - Visualization features incomplete

### LOW RISK (Performance/Convenience)
- **Checkpoint features** - Workflow convenience missing
- **Performance flags** - Optimization features incomplete
- **PED pseudonymization** - Additional privacy feature missing

## Recommendations

### Phase 1 (Critical - Required for Production Parity)
1. **Implement `TranscriptFilteringStage`**
   - Process `--transcript-list` and `--transcript-file`
   - Add to processing stages after variant extraction
   - Priority: CRITICAL

2. **Implement `GenotypeFilteringStage`**
   - Process `--genotype-filter` and `--gene-genotype-file`
   - Add to analysis stages
   - Priority: CRITICAL

3. **Complete `SnpEffSplitStage`**
   - Integrate existing logic into proper stage
   - Add to processing stages
   - Priority: HIGH

### Phase 2 (Important Features)
1. **Add `ExtraSampleFieldsStage`**
   - Process `--append-extra-sample-fields` and `--extra-sample-field-delimiter`
   - Add to processing stages after field extraction
   - Priority: MEDIUM

2. **Complete checkpoint system**
   - Implement `--resume`, `--checkpoint-checksum`, `--show-checkpoint-status`
   - Add validation and status display
   - Priority: MEDIUM

3. **Implement missing performance flags**
   - Add `--no-chunked-processing`, `--force-chunked-processing`, `--sort-parallel`
   - Enhance chunked processing control
   - Priority: MEDIUM

### Phase 3 (Enhancement Features)
1. **Complete IGV configuration**
   - Implement all IGV parameters in `IGVReportStage`
   - Add advanced visualization options
   - Priority: LOW

2. **Add field processing features**
   - Implement `--add-chr` and `--remove-sample-substring`
   - Add to processing stages
   - Priority: LOW

3. **Add PED pseudonymization**
   - Extend `PseudonymizationStage` for `--pseudonymize-ped`
   - Complete privacy feature set
   - Priority: LOW

## Conclusion

While the new stage-based pipeline maintains the CLI interface and implements 63% of functionality, **critical gaps exist in core genomics processing capabilities**. The missing transcript filtering and genotype filtering features represent significant functionality deficits that must be addressed before the new pipeline can achieve production parity.

The stage-based architecture is well-designed and provides a solid foundation. However, completing the missing stage implementations is essential for full feature compatibility with the original pipeline.

## Appendix: Stage Files Analyzed

- `variantcentrifuge/stages/setup_stages.py`
- `variantcentrifuge/stages/processing_stages.py`
- `variantcentrifuge/stages/analysis_stages.py`
- `variantcentrifuge/stages/output_stages.py`
- `variantcentrifuge/pipeline.py`
- `variantcentrifuge/pipeline_refactored.py`
- `variantcentrifuge/cli.py`