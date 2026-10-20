"""Statistics generation and ClinVar PM5 annotation stages."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..checkpoint_models import AtomicFileOperation
from ..pipeline_core import PipelineContext, Stage
from ..sample_subset import vcf_sample_count

logger = logging.getLogger(__name__)


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
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        return {"chunked_analysis", "inheritance_analysis", "variant_scoring"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system."""
        if context.config.get("no_stats"):
            return context

        stats_output = context.config.get("stats_output_file")
        if not stats_output:
            default_stats_path = context.workspace.get_intermediate_path("statistics.tsv")
            if default_stats_path.exists():
                stats_output = str(default_stats_path)
                context.config["stats_output_file"] = stats_output
                logger.info(f"Restored statistics file path: {stats_output}")
            else:
                logger.warning("Statistics file not found during checkpoint skip")

        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Generate statistics."""
        if context.config.get("no_stats"):
            logger.debug("Statistics generation disabled")
            return context

        if context.config.get("use_chunked_processing") and not context.config.get(
            "dataframe_chunked_analysis_complete"
        ):
            logger.info("Statistics will be generated after chunked processing completes")
            return context

        df = context.current_dataframe
        if df is None:
            raise RuntimeError("No DataFrame loaded for statistics")

        stats_config_path = context.config.get("stats_config")
        config = stats_config_path or self._get_default_stats_config(vcf_sample_count(context))

        from . import analysis_stages

        engine = analysis_stages.StatsEngine(config)
        stats = engine.compute(df)

        context.statistics = stats

        stats_output = context.config.get("stats_output_file")
        if not stats_output and not context.config.get("no_stats"):
            stats_output = context.workspace.get_intermediate_path("statistics.tsv")
            logger.debug(f"Created default statistics output path: {stats_output}")

        if stats_output:
            context.config["stats_output_file"] = str(stats_output)
            context.stats_output_tsv = Path(stats_output)
            self._write_statistics(stats, stats_output)
            logger.debug(f"Statistics written to: {stats_output}")

        return context

    def _write_statistics(self, stats: dict[str, Any], output_file: str) -> None:
        """Write statistics to file in a consistent format."""
        all_data = []

        if "dataset" in stats and not stats["dataset"].empty:
            dataset_df = stats["dataset"].copy()
            dataset_df["category"] = "dataset"
            all_data.append(dataset_df)

        if "genes" in stats and not stats["genes"].empty:
            genes_df = stats["genes"].copy()
            genes_df["category"] = "gene"
            all_data.append(genes_df)

        for key, value in stats.items():
            if key not in ["dataset", "genes"] and hasattr(value, "empty") and not value.empty:
                other_df = value.copy()
                other_df["category"] = key
                all_data.append(other_df)

        with AtomicFileOperation(output_file) as temp_out:
            if all_data:
                combined_df = pd.concat(all_data, ignore_index=True, sort=False)
                combined_df.to_csv(temp_out, sep="\t", index=False)
            else:
                with open(temp_out, "w") as f:
                    f.write("statistic\tvalue\tcategory\n")
                    f.write("total_variants\t0\tdataset\n")

    def _get_default_stats_config(self, n_vcf_samples: int | None = None) -> dict[str, Any]:
        """Get default statistics configuration using basic columns.

        ``n_vcf_samples``: VCF sample count when extraction kept a subset (#212), so
        ``total_samples`` still reports the VCF, not the extracted, sample count.
        """
        n_gt_cols = "len([col for col in df.columns if '_GT' in col])"
        return {
            "stats_version": "1.0",
            "description": "Default statistics configuration using basic columns",
            "dataset_stats": [
                {
                    "name": "total_variants",
                    "expression": "len(df)",
                    "description": "Total number of variants",
                },
                {
                    "name": "total_samples",
                    "expression": n_gt_cols if n_vcf_samples is None else str(n_vcf_samples),
                    "description": "Total number of samples",
                },
                {
                    "name": "unique_genes",
                    "expression": "df['GENE'].nunique() if 'GENE' in df.columns else 0",
                    "description": "Number of unique genes",
                },
                {
                    "name": "high_impact_variants",
                    "expression": "(df['IMPACT'] == 'HIGH').sum() if 'IMPACT' in df.columns else 0",
                    "description": "Number of HIGH impact variants",
                },
                {
                    "name": "moderate_impact_variants",
                    "expression": (
                        "(df['IMPACT'] == 'MODERATE').sum() if 'IMPACT' in df.columns else 0"
                    ),
                    "description": "Number of MODERATE impact variants",
                },
            ],
            "gene_stats": [
                {
                    "name": "variant_count",
                    "expression": "len(group_df)",
                    "groupby": "GENE",
                    "description": "Number of variants per gene",
                },
                {
                    "name": "high_impact_count",
                    "expression": (
                        "(group_df['IMPACT'] == 'HIGH').sum() "
                        "if 'IMPACT' in group_df.columns else 0"
                    ),
                    "groupby": "GENE",
                    "description": "Number of HIGH impact variants per gene",
                },
                {
                    "name": "moderate_impact_count",
                    "expression": (
                        "(group_df['IMPACT'] == 'MODERATE').sum() "
                        "if 'IMPACT' in group_df.columns else 0"
                    ),
                    "groupby": "GENE",
                    "description": "Number of MODERATE impact variants per gene",
                },
            ],
        }


class ClinVarPM5Stage(Stage):
    """Annotate variants with ACMG PM5 criterion from ClinVar lookup."""

    @property
    def name(self) -> str:
        return "clinvar_pm5"

    @property
    def description(self) -> str:
        return "Annotate variants with ClinVar PM5 criterion"

    @property
    def dependencies(self) -> set[str]:
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        return {"variant_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        lookup_path = context.config.get("clinvar_pm5_lookup")
        if not lookup_path:
            logger.debug("No --clinvar-pm5-lookup provided, skipping PM5 annotation")
            return context

        if context.current_dataframe is None:
            logger.warning("No DataFrame available for PM5 annotation")
            return context

        from ..clinvar_pm5 import annotate_pm5, load_pm5_lookup

        lookup = load_pm5_lookup(lookup_path)

        df = context.current_dataframe
        gene_col = "GENE" if "GENE" in df.columns else "ANN[0].GENE"
        hgvs_col = "HGVS_P" if "HGVS_P" in df.columns else "ANN[0].HGVS_P"

        if gene_col not in df.columns:
            logger.warning("Gene column not found in DataFrame, skipping PM5")
            return context

        if hgvs_col not in df.columns:
            logger.warning("HGVS_P column not found in DataFrame, skipping PM5")
            return context

        logger.info("Running PM5 annotation on %d variants", len(df))
        context.current_dataframe = annotate_pm5(df, lookup, gene_col=gene_col, hgvs_col=hgvs_col)
        logger.info("PM5 annotation complete")
        return context
