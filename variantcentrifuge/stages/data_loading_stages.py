"""DataFrame loading and custom annotation stages."""

from __future__ import annotations

import logging
from pathlib import Path

from ..dataframe_optimizer import load_optimized_dataframe, should_use_memory_passthrough
from ..pipeline_core import PipelineContext, Stage
from .genotype_matrix_builder import _create_standard_column_aliases

logger = logging.getLogger(__name__)


class DataFrameLoadingStage(Stage):
    """Load TSV data into DataFrame or prepare for chunked processing."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "dataframe_loading"

    @property
    def checkpoint_resume_policy(self) -> str:
        """DataFrame state can be restored from TSV artifacts."""
        return "restore"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load data into DataFrame for analysis"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return set()

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        return {
            "field_extraction",
            "parallel_complete_processing",
            "genotype_replacement",
            "phenotype_integration",
            "extra_column_removal",
        }

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system."""
        input_file = self._find_input_file(context)

        if not input_file or not Path(input_file).exists():
            raise RuntimeError(
                "Cannot restore dataframe_loading from checkpoint: no TSV input file found"
            )

        dataframe_chunks_done = context.config.get("dataframe_chunked_analysis_complete", False)
        if self._should_use_chunks(context, input_file) and not dataframe_chunks_done:
            logger.info("File too large, will use chunked processing (checkpoint skip)")
            context.config["use_chunked_processing"] = True
            return context

        logger.info(f"Restoring DataFrame from {input_file} (checkpoint skip)")
        compression = "gzip" if str(input_file).endswith(".gz") else None

        df, rename_map = load_optimized_dataframe(
            str(input_file),
            sep="\t",
            compression=compression,
            on_bad_lines="warn",
        )

        context.current_dataframe = df
        context.column_rename_map = rename_map
        logger.info(f"Restored {len(df)} variants into DataFrame after checkpoint skip")
        _create_standard_column_aliases(df)

        if should_use_memory_passthrough(df):
            context.variants_df = df
            logger.info("DataFrame stored for in-memory pass-through (checkpoint skip)")

        return context

    def _find_input_file(self, context: PipelineContext) -> Path:
        """Find the most recent TSV file in the pipeline."""
        if hasattr(context, "extra_columns_removed_tsv") and context.extra_columns_removed_tsv:
            return context.extra_columns_removed_tsv
        elif hasattr(context, "phenotypes_added_tsv") and context.phenotypes_added_tsv:
            return context.phenotypes_added_tsv
        elif hasattr(context, "genotype_replaced_tsv") and context.genotype_replaced_tsv:
            return context.genotype_replaced_tsv
        elif hasattr(context, "extracted_tsv") and context.extracted_tsv:
            return context.extracted_tsv
        else:
            input_file = context.data
            if not input_file:
                raise RuntimeError(
                    "Cannot restore dataframe_loading from checkpoint: no TSV input file found"
                )
            if input_file and (
                str(input_file).endswith(".vcf.gz") or str(input_file).endswith(".vcf")
            ):
                logger.warning(
                    f"DataFrameLoadingStage received VCF file instead of TSV: {input_file}"
                )
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
                            return Path(tsv_file)
                else:
                    raise ValueError(
                        f"DataFrameLoadingStage requires a TSV file, but got VCF: {input_file}. "
                        f"No TSV files found in context."
                    )
            return Path(input_file)

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load data into DataFrame or prepare chunked processing."""
        input_file = self._find_input_file(context)

        dataframe_chunks_done = context.config.get("dataframe_chunked_analysis_complete", False)
        if self._should_use_chunks(context, input_file) and not dataframe_chunks_done:
            logger.info("File too large, will use chunked processing")
            context.config["use_chunked_processing"] = True
            return context

        logger.info(f"Loading data from {input_file}")
        compression = "gzip" if str(input_file).endswith(".gz") else None

        df, rename_map = load_optimized_dataframe(
            str(input_file),
            sep="\t",
            compression=compression,
            on_bad_lines="warn",
        )

        context.current_dataframe = df
        context.column_rename_map = rename_map
        logger.info(f"Loaded {len(df)} variants into DataFrame")

        _create_standard_column_aliases(df)

        if should_use_memory_passthrough(df):
            context.variants_df = df
            logger.info("DataFrame stored for in-memory pass-through")

        return context

    def _should_use_chunks(self, context: PipelineContext, file_path: Path) -> bool:
        """Determine if chunked processing should be used."""
        if context.config.get("no_chunked_processing"):
            return False

        if context.config.get("force_chunked_processing"):
            return True

        try:
            sample_count = len(context.vcf_samples) if context.vcf_samples else 0
            file_size_mb = file_path.stat().st_size / (1024 * 1024)
            estimated_variants = max(1, int(file_size_mb * 1000))

            if sample_count > 0:
                estimated_memory_mb = (estimated_variants * sample_count * 8) / (1024 * 1024)

                if sample_count > 1000:
                    threshold_mb = context.config.get("chunk_threshold_mb", 50)
                    memory_threshold_mb = 500
                elif sample_count > 100:
                    threshold_mb = context.config.get("chunk_threshold_mb", 100)
                    memory_threshold_mb = 1000
                else:
                    threshold_mb = context.config.get("chunk_threshold_mb", 500)
                    memory_threshold_mb = 2000

                if file_size_mb > threshold_mb:
                    logger.info(
                        f"Using chunked processing: file size {file_size_mb:.1f}MB > "
                        f"{threshold_mb}MB threshold"
                    )
                    return True

                if estimated_memory_mb > memory_threshold_mb:
                    logger.info(
                        f"Using chunked processing: estimated memory {estimated_memory_mb:.1f}MB > "
                        f"{memory_threshold_mb}MB threshold"
                    )
                    logger.info(
                        f"Memory estimation: {estimated_variants} variants x {sample_count} samples"
                    )
                    return True
            else:
                threshold_mb = context.config.get("chunk_threshold_mb", 500)
                return bool(file_size_mb > threshold_mb)

        except Exception as e:
            logger.debug(f"Error in chunking decision: {e}")
            return False

        return False

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        for attr in [
            "extra_columns_removed_tsv",
            "phenotypes_added_tsv",
            "genotype_replaced_tsv",
            "extracted_tsv",
        ]:
            if hasattr(context, attr):
                tsv_file = getattr(context, attr)
                if tsv_file and tsv_file.exists():
                    return [tsv_file]
        if hasattr(context, "data") and context.data:
            return [context.data]
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        return []


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
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Run after chunked analysis when that stage is active."""
        return {"chunked_analysis"}

    def _has_annotation_config(self) -> bool:
        """Check if annotation config stage exists."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply custom annotations."""
        if not context.annotation_configs:
            logger.debug("No custom annotations configured")
            return context

        if context.config.get("use_chunked_processing"):
            logger.info("Annotations will be applied during chunked processing")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for annotation")
            return context

        annotation_cfg = {
            "annotate_bed_files": context.annotation_configs.get("bed_files", []),
            "annotate_gene_lists": context.annotation_configs.get("gene_lists", []),
            "annotate_json_genes": context.annotation_configs.get("json_genes", []),
            "json_gene_mapping": context.annotation_configs.get("json_mapping", ""),
            "json_genes_as_columns": context.config.get("json_genes_as_columns", False),
        }
        from . import analysis_stages

        features = analysis_stages.load_custom_features(annotation_cfg)

        logger.info(f"Applying {len(features)} custom annotation features")
        df = analysis_stages.annotate_dataframe_with_features(df, features)

        context.current_dataframe = df
        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        input_files = []
        if hasattr(context, "annotation_configs") and context.annotation_configs:
            for bed_file in context.annotation_configs.get("bed_files", []):
                if Path(bed_file).exists():
                    input_files.append(Path(bed_file))
            for gene_list in context.annotation_configs.get("gene_lists", []):
                if Path(gene_list).exists():
                    input_files.append(Path(gene_list))
            for json_file in context.annotation_configs.get("json_genes", []):
                if Path(json_file).exists():
                    input_files.append(Path(json_file))
        return input_files

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        return []
