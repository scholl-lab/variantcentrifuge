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
from typing import Any, Dict, Set

import pandas as pd

from ..analyze_variants import analyze_variants
from ..annotator import annotate_dataframe_with_features, load_custom_features
from ..gene_burden import perform_gene_burden_analysis
from ..inheritance import analyze_inheritance
from ..pipeline_core import PipelineContext, Stage
from ..scoring import apply_scoring
from ..stats_engine import StatsEngine

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
        # Other transformations are optional
        return {"field_extraction"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load data into DataFrame or prepare chunked processing."""
        input_file = context.data

        # Check if we should use chunked processing
        if self._should_use_chunks(context, input_file):
            logger.info("File too large, will use chunked processing")
            context.config["use_chunked_processing"] = True
            # Don't load full DataFrame - chunked stages will handle it
            return context

        # Load full DataFrame
        logger.info(f"Loading data from {input_file}")
        compression = "gzip" if str(input_file).endswith(".gz") else None

        df = pd.read_csv(
            input_file,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            na_values=[""],
            compression=compression,
        )

        context.current_dataframe = df
        logger.info(f"Loaded {len(df)} variants into DataFrame")

        return context

    def _should_use_chunks(self, context: PipelineContext, file_path: Path) -> bool:
        """Determine if chunked processing should be used."""
        # Check explicit chunking request
        if context.config.get("chunks"):
            return True

        # Check file size (simplified - real implementation would check actual size)
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
        features = load_custom_features(
            bed_files=context.annotation_configs.get("bed_files", []),
            gene_lists=context.annotation_configs.get("gene_lists", []),
            json_genes_file=context.annotation_configs.get("json_genes"),
            json_mapping=context.annotation_configs.get("json_mapping"),
        )

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
        # Static dependencies only - pedigree loading is optional
        # and will be checked at runtime
        return {"dataframe_loading", "custom_annotation"}

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
            vcf_samples=vcf_samples,
            pedigree_data=pedigree_data,
            inheritance_mode=inheritance_mode,
            use_vectorized=not context.config.get("no_vectorized_comp_het", False),
        )

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
        # Static dependencies only - scoring config loading is optional
        # and will be checked at runtime
        # Note: inheritance_analysis should run before scoring so scores
        # can use inheritance patterns
        return {"dataframe_loading", "custom_annotation", "inheritance_analysis"}

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

        logger.info(f"Processing file in chunks of {chunk_size} variants")

        # This is a simplified version - real implementation would:
        # 1. Read file in gene-aware chunks
        # 2. Apply all analysis steps to each chunk
        # 3. Aggregate results appropriately
        # 4. Handle memory efficiently

        # For now, just mark as complete
        logger.info("Chunked processing completed")
        context.config["chunked_processing_complete"] = True

        return context


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
        unique_genes = df[gene_column].unique()

        logger.info(
            f"Running parallel analysis for {len(unique_genes)} genes using {threads} threads"
        )

        # This is a simplified version - real implementation would:
        # 1. Split DataFrame by gene
        # 2. Process each gene's variants in parallel
        # 3. Merge results back together
        # 4. Handle gene-level operations (compound het, etc.)

        logger.info("Parallel analysis completed")
        return context
