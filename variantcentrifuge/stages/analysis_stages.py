"""
Analysis stages for variant processing and computation.

This module contains stages that perform analysis on the extracted data:
- DataFrame loading and chunked processing
- Custom annotations (BED, gene lists, JSON)
- Inheritance pattern analysis
- Variant scoring
- Statistics generation
- Gene burden analysis
"""

import logging
from pathlib import Path
from typing import Any, Dict, Optional, Set

import pandas as pd

from ..analyze_variants import analyze_variants
from ..annotator import annotate_dataframe_with_features, load_custom_features
from ..gene_burden import perform_gene_burden_analysis
from ..inheritance import analyze_inheritance
from ..pipeline_core import PipelineContext, Stage
from ..scoring import apply_scoring
from ..stats_engine import StatsEngine
from ..links import add_links_to_table

logger = logging.getLogger(__name__)


class DataFrameLoadingStage(Stage):
    """Load TSV data into DataFrame or prepare for chunked processing."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "dataframe_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load data into DataFrame for analysis"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Always depends on field extraction (the TSV file must exist)
        # We don't add other dependencies here to avoid circular dependencies
        # Instead, we handle the data flow in _process() method
        return {"field_extraction"}

    @property
    def soft_dependencies(self) -> Set[str]:
        """Return the set of stage names that should run before if present."""
        # Prefer to run after all data transformation stages
        return {"genotype_replacement", "phenotype_integration", "extra_column_removal"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load data into DataFrame or prepare chunked processing."""
        # Use the most recent TSV file in the pipeline
        # Priority: extra_columns_removed > phenotypes_added > genotype_replaced > extracted
        if hasattr(context, "extra_columns_removed_tsv") and context.extra_columns_removed_tsv:
            input_file = context.extra_columns_removed_tsv
        elif hasattr(context, "phenotypes_added_tsv") and context.phenotypes_added_tsv:
            input_file = context.phenotypes_added_tsv
        elif hasattr(context, "genotype_replaced_tsv") and context.genotype_replaced_tsv:
            input_file = context.genotype_replaced_tsv
        elif hasattr(context, "extracted_tsv") and context.extracted_tsv:
            input_file = context.extracted_tsv
        else:
            # Fallback to context.data, but verify it's a TSV file
            input_file = context.data
            if input_file and (
                str(input_file).endswith(".vcf.gz") or str(input_file).endswith(".vcf")
            ):
                logger.warning(
                    f"DataFrameLoadingStage received VCF file instead of TSV: {input_file}"
                )
                # Try to find any TSV file in order of preference
                for attr in [
                    "extracted_tsv",
                    "genotype_replaced_tsv",
                    "phenotypes_added_tsv",
                    "extra_columns_removed_tsv",
                ]:
                    if hasattr(context, attr):
                        tsv_file = getattr(context, attr)
                        if tsv_file and tsv_file.exists():
                            logger.info(f"Using {attr} instead: {tsv_file}")
                            input_file = tsv_file
                            break
                else:
                    raise ValueError(
                        f"DataFrameLoadingStage requires a TSV file, but got VCF: {input_file}. "
                        f"No TSV files found in context."
                    )

        # Check if we should use chunked processing
        if self._should_use_chunks(context, input_file):
            logger.info("File too large, will use chunked processing")
            context.config["use_chunked_processing"] = True
            # Don't load full DataFrame - chunked stages will handle it
            return context

        # Load full DataFrame
        logger.info(f"Loading data from {input_file}")
        compression = "gzip" if str(input_file).endswith(".gz") else None

        # Use quoting=3 (QUOTE_NONE) to handle data with quotes
        # Use low_memory=False to avoid dtype inference issues
        df = pd.read_csv(
            input_file,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            na_values=[""],
            compression=compression,
            quoting=3,  # QUOTE_NONE - don't treat quotes specially
            low_memory=False,
            on_bad_lines="warn",  # Warn about problematic lines instead of failing
        )

        context.current_dataframe = df
        logger.info(f"Loaded {len(df)} variants into DataFrame")

        return context

    def _should_use_chunks(self, context: PipelineContext, file_path: Path) -> bool:
        """Determine if chunked processing should be used."""
        # Check explicit chunking request
        if context.config.get("chunks"):
            return True

        # Check file size to determine if chunking is needed
        try:
            size_mb = file_path.stat().st_size / (1024 * 1024)
            threshold_mb = context.config.get("chunk_threshold_mb", 500)
            return size_mb > threshold_mb
        except Exception:
            return False


class CustomAnnotationStage(Stage):
    """Apply custom annotations from BED files, gene lists, and JSON."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "custom_annotation"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Apply custom annotations to variants"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Always depends on dataframe being loaded
        # Annotation config is handled in setup stages if needed
        return {"dataframe_loading"}

    def _has_annotation_config(self) -> bool:
        """Check if annotation config stage exists."""
        return True  # Simplified

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply custom annotations."""
        if not context.annotation_configs:
            logger.debug("No custom annotations configured")
            return context

        # Check if using chunked processing
        if context.config.get("use_chunked_processing"):
            logger.info("Annotations will be applied during chunked processing")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for annotation")
            return context

        # Load annotation features
        # Create config dict for load_custom_features
        annotation_cfg = {
            "annotate_bed_files": context.annotation_configs.get("bed_files", []),
            "annotate_gene_lists": context.annotation_configs.get("gene_lists", []),
            "annotate_json_genes": context.annotation_configs.get("json_genes", []),
            "json_gene_mapping": context.annotation_configs.get("json_mapping", ""),
            "json_genes_as_columns": context.config.get("json_genes_as_columns", False),
        }
        features = load_custom_features(annotation_cfg)

        # Apply annotations
        logger.info(f"Applying {len(features)} custom annotation features")
        df = annotate_dataframe_with_features(df, features)

        context.current_dataframe = df
        return context


class InheritanceAnalysisStage(Stage):
    """Calculate Mendelian inheritance patterns."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "inheritance_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Calculate inheritance patterns"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depend on dataframe_loading as a hard requirement
        # custom_annotation is optional and will run before if present
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> Set[str]:
        """Return the set of stage names that should run before if present."""
        # Prefer to run after custom_annotation if it exists
        return {"custom_annotation"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Calculate inheritance patterns."""
        # Check if we have necessary data
        df = context.current_dataframe
        if df is None and not context.config.get("use_chunked_processing"):
            logger.warning("No DataFrame loaded for inheritance analysis")
            return context

        # Determine inheritance mode
        inheritance_mode = context.config.get("inheritance_mode", "simple")

        # Check if using chunked processing
        if context.config.get("use_chunked_processing"):
            logger.info(
                f"Inheritance analysis ({inheritance_mode} mode) "
                "will be applied during chunked processing"
            )
            return context

        # Get samples info
        vcf_samples = context.vcf_samples or []
        pedigree_data = context.pedigree_data

        if not vcf_samples:
            logger.warning("No VCF samples available for inheritance analysis")
            return context

        logger.info(
            f"Calculating inheritance patterns ({inheritance_mode} mode) "
            f"for {len(vcf_samples)} samples"
        )

        # Apply inheritance analysis
        df = analyze_inheritance(
            df=df,
            sample_list=vcf_samples,
            pedigree_data=pedigree_data,
            use_vectorized_comp_het=not context.config.get("no_vectorized_comp_het", False),
        )

        # Process output based on inheritance mode
        from ..inheritance.analyzer import process_inheritance_output

        df = process_inheritance_output(df, inheritance_mode)

        context.current_dataframe = df
        return context


class VariantScoringStage(Stage):
    """Apply variant scoring formulas."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "variant_scoring"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Calculate variant scores"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depend on dataframe_loading as a hard requirement
        # Other stages (custom_annotation, inheritance_analysis) are optional
        # The pipeline runner will ensure correct ordering based on what stages are present
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> Set[str]:
        """Return the set of stage names that should run before if present."""
        # Prefer to run after annotations and inheritance if they exist
        return {"custom_annotation", "inheritance_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply variant scoring."""
        if not context.scoring_config:
            logger.debug("No scoring configuration loaded")
            return context

        # Check if using chunked processing
        if context.config.get("use_chunked_processing"):
            logger.info("Scoring will be applied during chunked processing")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for scoring")
            return context

        logger.info("Applying variant scoring formulas")

        # Apply scoring
        df = apply_scoring(df, context.scoring_config)

        context.current_dataframe = df
        return context


class StatisticsGenerationStage(Stage):
    """Generate summary statistics."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "statistics_generation"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Generate summary statistics"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Only depends on dataframe being loaded
        # Statistics can be run on whatever columns are available
        return {"dataframe_loading"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - read-only computation on DataFrame

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate statistics."""
        if context.config.get("no_stats"):
            logger.debug("Statistics generation disabled")
            return context

        # Check if using chunked processing
        if context.config.get("use_chunked_processing"):
            logger.info("Statistics will be aggregated during chunked processing")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for statistics")
            return context

        # Get stats config
        stats_config_path = context.config.get("stats_config")
        config = stats_config_path if stats_config_path else {}

        logger.info("Generating summary statistics")

        # Create stats engine and compute statistics
        engine = StatsEngine(config)
        stats = engine.compute(df)

        context.statistics = stats

        # Write stats if output file specified
        stats_output = context.config.get("stats_output_file")
        if stats_output:
            self._write_statistics(stats, stats_output)

        return context

    def _write_statistics(self, stats: Dict[str, Any], output_file: str) -> None:
        """Write statistics to file."""
        with open(output_file, "w") as f:
            # Write dataset statistics
            if "dataset" in stats and not stats["dataset"].empty:
                f.write("## Dataset Statistics\n")
                stats["dataset"].to_csv(f, sep="\t", index=False, header=False)

            # Write gene statistics if present
            if "genes" in stats and not stats["genes"].empty:
                f.write("\n## Gene Statistics\n")
                stats["genes"].to_csv(f, sep="\t", index=False)

            # Write any grouped statistics
            for key, df in stats.items():
                if key not in ["dataset", "genes"] and not df.empty:
                    f.write(f"\n## {key.title()} Statistics\n")
                    df.to_csv(f, sep="\t", index=False)

        logger.info(f"Wrote statistics to {output_file}")


class VariantAnalysisStage(Stage):
    """Run variant-level analysis to add analysis columns."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "variant_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Perform variant-level analysis"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Depends on having a DataFrame with variants
        deps = {"dataframe_loading"}
        return deps

    @property
    def soft_dependencies(self) -> Set[str]:
        """Return the set of stage names that should run before if present."""
        # Run after variant_identifier to preserve VAR_ID column
        return {"variant_identifier"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Run variant analysis."""
        if context.current_dataframe is None:
            logger.warning("No DataFrame available for variant analysis")
            return context

        df = context.current_dataframe

        # Prepare config for analyze_variants
        analysis_config = {
            "vcf_sample_names": context.vcf_samples or [],
            "sample_list": ",".join(context.vcf_samples or []),
            "no_db_links": context.config.get("no_db_links", False),
            "base_name": context.workspace.base_name,
            "reference": context.config.get("reference", "GRCh37"),
        }

        # Create temporary TSV file for analyze_variants
        temp_tsv = context.workspace.get_intermediate_path("temp_for_analysis.tsv")
        df.to_csv(temp_tsv, sep="\t", index=False)

        logger.info("Running variant-level analysis")

        # Run analyze_variants and collect results
        analysis_results = []
        with open(temp_tsv, "r") as inp:
            for line in analyze_variants(inp, analysis_config):
                analysis_results.append(line)

        # Parse results back into DataFrame
        if analysis_results:
            # First line is header
            header = analysis_results[0].split("\t")
            data = [line.split("\t") for line in analysis_results[1:]]

            # Create new DataFrame with analysis results
            analysis_df = pd.DataFrame(data, columns=header)

            # The analyze_variants function may filter or reorder rows, so we need to match them properly
            # We'll use the key columns (CHROM, POS, REF, ALT) to merge back any missing columns
            key_cols = ["CHROM", "POS", "REF", "ALT"]
            
            # Check if we have the key columns in both dataframes
            if all(col in df.columns for col in key_cols) and all(col in analysis_df.columns for col in key_cols):
                # Get columns from original df that are not in analysis_df
                original_only_cols = [col for col in df.columns if col not in analysis_df.columns]
                
                if original_only_cols:
                    # Create a subset of original df with key columns and missing columns
                    merge_df = df[key_cols + original_only_cols].copy()
                    
                    # Merge to get the missing columns
                    analysis_df = pd.merge(analysis_df, merge_df, on=key_cols, how="left")
                    
                    # Reorder columns to put VAR_ID first if it exists
                    if "VAR_ID" in analysis_df.columns:
                        cols = analysis_df.columns.tolist()
                        cols.remove("VAR_ID")
                        cols = ["VAR_ID"] + cols
                        analysis_df = analysis_df[cols]
                    
                    # Put Custom_Annotation after GT if both exist
                    if "Custom_Annotation" in analysis_df.columns and "GT" in analysis_df.columns:
                        cols = analysis_df.columns.tolist()
                        cols.remove("Custom_Annotation")
                        gt_idx = cols.index("GT")
                        cols.insert(gt_idx + 1, "Custom_Annotation")
                        analysis_df = analysis_df[cols]
                    
                    logger.debug(f"Preserved columns from original dataframe: {original_only_cols}")
            else:
                logger.warning("Could not preserve columns from original dataframe - missing key columns for merging")

            context.current_dataframe = analysis_df
            logger.info(f"Variant analysis complete: {len(analysis_df)} variants analyzed")

        # Clean up temp file
        if temp_tsv.exists():
            temp_tsv.unlink()

        return context


class GeneBurdenAnalysisStage(Stage):
    """Perform gene burden analysis for case-control studies."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "gene_burden_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Perform gene burden analysis"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Depends on having a DataFrame with variants
        return {"dataframe_loading", "custom_annotation"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Perform gene burden analysis."""
        if not context.config.get("perform_gene_burden"):
            logger.debug("Gene burden analysis not requested")
            return context

        # Check required inputs
        case_samples = context.config.get("case_samples", [])
        control_samples = context.config.get("control_samples", [])

        if not case_samples or not control_samples:
            logger.warning("Case/control samples not defined for gene burden analysis")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for gene burden analysis")
            return context

        logger.info(
            f"Performing gene burden analysis: {len(case_samples)} cases, "
            f"{len(control_samples)} controls"
        )

        # Run analysis with proper config
        burden_config = {
            "gene_burden_mode": context.config.get("gene_burden_mode", "alleles"),
            "correction_method": context.config.get("correction_method", "fdr"),
            "confidence_interval_method": context.config.get(
                "confidence_interval_method", "normal_approx"
            ),
            "confidence_interval_alpha": context.config.get("confidence_interval_alpha", 0.05),
        }

        # First analyze variants to get counts
        analyzed_df = analyze_variants(
            df=df,
            case_samples=case_samples,
            control_samples=control_samples,
            gene_column=context.config.get("gene_column", "GENE"),
        )

        # Then perform gene burden analysis
        burden_results = perform_gene_burden_analysis(df=analyzed_df, cfg=burden_config)

        context.gene_burden_results = burden_results

        # Write results if output specified
        burden_output = context.config.get("gene_burden_output")
        if burden_output:
            burden_results.to_csv(burden_output, sep="\t", index=False)
            logger.info(f"Wrote gene burden results to {burden_output}")

        return context


class ChunkedAnalysisStage(Stage):
    """Process large files in chunks with all analysis steps."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "chunked_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Process large files in memory-efficient chunks"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Depends on configuration but not full DataFrame loading
        deps = {"field_extraction", "phenotype_integration"}
        if self._has_configs():
            deps.update({"annotation_config_loading", "scoring_config_loading", "pedigree_loading"})
        return deps

    def _has_configs(self) -> bool:
        """Check if config stages exist."""
        return True  # Simplified

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process file in chunks."""
        if not context.config.get("use_chunked_processing"):
            logger.debug("Chunked processing not enabled")
            return context

        chunk_size = context.config.get("chunks", 10000)

        # Determine input file
        if context.extra_columns_removed_tsv:
            input_file = context.extra_columns_removed_tsv
        elif context.phenotypes_added_tsv:
            input_file = context.phenotypes_added_tsv
        elif context.genotype_replaced_tsv:
            input_file = context.genotype_replaced_tsv
        elif context.extracted_tsv:
            input_file = context.extracted_tsv
        else:
            logger.error("No TSV file available for chunked processing")
            return context

        logger.info(f"Processing {input_file} in chunks of {chunk_size} variants")

        # Import required modules for chunk processing
        from variantcentrifuge.helpers import find_column
        from variantcentrifuge.analyze_variants import analyze_variants
        from variantcentrifuge.scoring import score_dataframe

        # Find gene column
        gene_col_name = None
        with open(input_file, "r") as f:
            header = f.readline().strip().split("\t")
            for gene_col in ["GENE", "Gene", "gene", "SYMBOL", "ANN[*].GENE", "EFF[*].GENE"]:
                if gene_col in header:
                    gene_col_name = gene_col
                    break

        if not gene_col_name:
            logger.warning("No GENE column found, falling back to non-chunked processing")
            return context

        # Sort file by gene column if not already sorted
        sorted_file = self._ensure_sorted_by_gene(input_file, gene_col_name, context)

        # Process chunks
        output_chunks = []
        chunk_num = 0

        for chunk_df in self._read_tsv_in_gene_chunks(sorted_file, gene_col_name, chunk_size):
            logger.info(f"Processing chunk {chunk_num + 1} with {len(chunk_df)} variants")

            # Apply scoring if configured
            if context.scoring_config:
                chunk_df = score_dataframe(chunk_df, context.scoring_config)

            # Apply analysis
            if not context.config.get("skip_analysis", False):
                analysis_cfg = {
                    "add_variant_id": context.config.get("add_variant_id", True),
                    "perform_gene_burden": context.config.get("perform_gene_burden", False),
                    "custom_gene_list": context.config.get("custom_gene_list"),
                    "case_samples_file": context.config.get("case_samples_file"),
                    "control_samples_file": context.config.get("control_samples_file"),
                    "vcf_samples": context.vcf_samples,
                    "aggregate_operation": context.config.get("aggregate_operation", "max"),
                    "aggregate_column": context.config.get("aggregate_column"),
                    "sample_values": context.config.get("sample_values", {}),
                }

                chunk_df = analyze_variants(chunk_df, context.pedigree_data or {}, analysis_cfg)

            output_chunks.append(chunk_df)
            chunk_num += 1

        # Combine all chunks
        if output_chunks:
            logger.info(f"Combining {len(output_chunks)} processed chunks")
            context.current_dataframe = pd.concat(output_chunks, ignore_index=True)
            context.config["chunked_processing_complete"] = True
            logger.info(
                f"Chunked processing completed: {len(context.current_dataframe)} total variants"
            )
        else:
            logger.warning("No chunks were processed")

        return context

    def _ensure_sorted_by_gene(
        self, input_file: Path, gene_col: str, context: PipelineContext
    ) -> Path:
        """Ensure the TSV file is sorted by gene column."""
        import subprocess
        import tempfile

        # Check if already sorted
        if context.config.get("assume_sorted", False):
            return input_file

        sorted_file = context.workspace.intermediate / f"{input_file.stem}.sorted.tsv"

        if sorted_file.exists():
            logger.info(f"Using existing sorted file: {sorted_file}")
            return sorted_file

        logger.info(f"Sorting {input_file} by {gene_col} column")

        # Find column index (1-based for sort command)
        with open(input_file, "r") as f:
            header = f.readline().strip().split("\t")
            try:
                col_idx = header.index(gene_col) + 1  # 1-based index
            except ValueError:
                logger.error(f"Column {gene_col} not found in header")
                return input_file

        # Use system sort for memory efficiency
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
            # Write header to temp file
            with open(input_file, "r") as f:
                tmp.write(f.readline())
            tmp_path = tmp.name

        # Sort data (skip header)
        sort_cmd = [
            "sort",
            "-t",
            "\t",  # Tab delimiter
            "-k",
            str(col_idx),  # Sort by gene column
            "--stable",  # Stable sort
            "-T",
            str(context.workspace.temp),  # Use temp directory
        ]

        memory_limit = context.config.get("sort_memory_limit", "1G")
        if memory_limit:
            sort_cmd.extend(["-S", memory_limit])

        # Sort command: tail -n +2 input | sort ... >> tmp_file
        with open(tmp_path, "a") as out_f:
            tail_proc = subprocess.Popen(
                ["tail", "-n", "+2", str(input_file)], stdout=subprocess.PIPE
            )
            sort_proc = subprocess.Popen(sort_cmd, stdin=tail_proc.stdout, stdout=out_f)
            tail_proc.stdout.close()
            sort_proc.wait()

        # Move temp file to final location
        import shutil

        shutil.move(tmp_path, sorted_file)

        logger.info(f"Sorting complete: {sorted_file}")
        return sorted_file

    def _read_tsv_in_gene_chunks(self, filepath: Path, gene_column: str, chunksize: int):
        """Read TSV file in gene-aware chunks."""
        import pandas as pd

        logger.info(
            f"Reading {filepath} in gene-aware chunks (column: {gene_column}, size: {chunksize})"
        )

        # Create reader for chunked reading
        reader = pd.read_csv(
            filepath,
            sep="\t",
            chunksize=chunksize,
            low_memory=False,
            dtype=str,  # Read all as strings to preserve data
        )

        gene_buffer = pd.DataFrame()
        chunks_yielded = 0

        for chunk in reader:
            # Combine with buffer
            if gene_buffer.empty:
                gene_buffer = chunk
            else:
                gene_buffer = pd.concat([gene_buffer, chunk], ignore_index=True)

            # Check buffer size
            if len(gene_buffer) > chunksize * 10:
                logger.warning(
                    f"Gene buffer has grown to {len(gene_buffer)} rows. "
                    f"This may indicate a gene with too many variants or unsorted data."
                )

            # Find gene boundaries
            if len(gene_buffer) >= chunksize:
                # Find where genes change
                gene_values = gene_buffer[gene_column].values
                gene_change_indices = [0]
                current_gene = gene_values[0]

                for i in range(1, len(gene_values)):
                    if gene_values[i] != current_gene:
                        gene_change_indices.append(i)
                        current_gene = gene_values[i]

                # Find split point after chunksize
                split_index = None
                for idx in gene_change_indices:
                    if idx >= chunksize:
                        split_index = idx
                        break

                if split_index is None:
                    # All rows belong to same gene or last gene extends beyond chunksize
                    # Keep accumulating unless we're at the last chunk
                    continue
                else:
                    # Yield complete genes
                    yield_df = gene_buffer.iloc[:split_index].copy()
                    yield yield_df
                    chunks_yielded += 1

                    # Keep remaining for next iteration
                    gene_buffer = gene_buffer.iloc[split_index:].reset_index(drop=True)

        # Yield any remaining data
        if not gene_buffer.empty:
            yield gene_buffer
            chunks_yielded += 1

        logger.info(f"Completed reading {chunks_yielded} gene-aware chunks")


class ParallelAnalysisOrchestrator(Stage):
    """Orchestrate parallel analysis of independent genes."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "parallel_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Run analysis in parallel by gene"

    @property
    def dependencies(self) -> Set[str]:
        """Return the set of stage names this stage depends on."""
        # Can run after basic data loading
        return {"dataframe_loading", "custom_annotation"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        # This stage manages its own parallelism
        return False

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Run analysis in parallel by gene."""
        if not context.config.get("parallel_analysis"):
            logger.debug("Parallel analysis not enabled")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for parallel analysis")
            return context

        threads = context.config.get("threads", 1)
        if threads <= 1:
            logger.debug("Thread count <= 1, skipping parallel analysis")
            return context

        gene_column = context.config.get("gene_column", "GENE")

        # Find the gene column if not specified
        if gene_column not in df.columns:
            for col in ["GENE", "Gene", "gene", "SYMBOL", "ANN[*].GENE", "EFF[*].GENE"]:
                if col in df.columns:
                    gene_column = col
                    break
            else:
                logger.warning("No gene column found, cannot perform parallel analysis by gene")
                return context

        unique_genes = df[gene_column].unique()
        num_genes = len(unique_genes)

        logger.info(f"Running parallel analysis for {num_genes} genes using {threads} threads")

        # Import required modules
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from variantcentrifuge.analyze_variants import analyze_variants
        from variantcentrifuge.scoring import score_dataframe
        import pickle
        import tempfile

        # Split DataFrame by gene
        gene_groups = df.groupby(gene_column)

        # Prepare analysis configuration
        analysis_cfg = {
            "add_variant_id": context.config.get("add_variant_id", True),
            "perform_gene_burden": context.config.get("perform_gene_burden", False),
            "custom_gene_list": context.config.get("custom_gene_list"),
            "case_samples_file": context.config.get("case_samples_file"),
            "control_samples_file": context.config.get("control_samples_file"),
            "vcf_samples": context.vcf_samples,
            "aggregate_operation": context.config.get("aggregate_operation", "max"),
            "aggregate_column": context.config.get("aggregate_column"),
            "sample_values": context.config.get("sample_values", {}),
        }

        # Track successful and failed genes
        successful_results = []
        failed_genes = []

        # Process genes in parallel
        with ProcessPoolExecutor(max_workers=threads) as executor:
            # Submit tasks
            futures = {}
            for gene_name, gene_df in gene_groups:
                future = executor.submit(
                    self._process_gene_group,
                    gene_name,
                    gene_df,
                    context.pedigree_data,
                    context.scoring_config,
                    analysis_cfg,
                    context.config.get("skip_analysis", False),
                )
                futures[future] = gene_name

            # Process results as they complete
            try:
                for future in as_completed(futures):
                    gene_name = futures[future]
                    try:
                        result_df = future.result()
                        if result_df is not None:
                            successful_results.append(result_df)
                            logger.debug(f"Successfully processed gene: {gene_name}")
                    except Exception as e:
                        logger.error(f"Failed to process gene {gene_name}: {str(e)}")
                        failed_genes.append(gene_name)
                        # Don't cancel other tasks, continue processing

            except KeyboardInterrupt:
                logger.warning("Parallel processing interrupted by user")
                # Cancel remaining futures
                for future in futures:
                    future.cancel()
                raise

        # Merge results
        if successful_results:
            logger.info(
                f"Merging results from {len(successful_results)} successfully processed genes"
            )
            context.current_dataframe = pd.concat(successful_results, ignore_index=True)

            # Restore original column order
            original_columns = df.columns.tolist()
            current_columns = context.current_dataframe.columns.tolist()
            # Add any new columns at the end
            new_columns = [col for col in current_columns if col not in original_columns]
            ordered_columns = original_columns + new_columns
            # Only reorder columns that exist in the result
            ordered_columns = [col for col in ordered_columns if col in current_columns]
            context.current_dataframe = context.current_dataframe[ordered_columns]

            if failed_genes:
                logger.warning(
                    f"Failed to process {len(failed_genes)} genes: {', '.join(failed_genes[:10])}"
                )
                if len(failed_genes) > 10:
                    logger.warning(f"... and {len(failed_genes) - 10} more")
        else:
            logger.error("No genes were successfully processed")
            context.current_dataframe = df  # Keep original data

        logger.info(
            f"Parallel analysis completed: {len(successful_results)}/{num_genes} genes processed"
        )
        return context

    @staticmethod
    def _process_gene_group(
        gene_name: str,
        gene_df: pd.DataFrame,
        pedigree_data: Optional[Dict],
        scoring_config: Optional[Dict],
        analysis_cfg: Dict,
        skip_analysis: bool,
    ) -> pd.DataFrame:
        """Process a single gene's variants (runs in worker process)."""
        import pandas as pd
        from variantcentrifuge.analyze_variants import analyze_variants
        from variantcentrifuge.scoring import score_dataframe

        # Make a copy to avoid modifying the original
        result_df = gene_df.copy()

        try:
            # Apply scoring if configured
            if scoring_config:
                result_df = score_dataframe(result_df, scoring_config)

            # Apply analysis if not skipped
            if not skip_analysis:
                result_df = analyze_variants(result_df, pedigree_data or {}, analysis_cfg)

            return result_df

        except Exception as e:
            # Log error and re-raise
            import logging

            logger = logging.getLogger(__name__)
            logger.error(f"Error processing gene {gene_name}: {str(e)}")
            raise
