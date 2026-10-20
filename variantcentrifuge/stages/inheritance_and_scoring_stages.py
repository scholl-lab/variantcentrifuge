"""Inheritance analysis and variant scoring stages."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ..gt_reconstruction import _find_per_sample_gt_columns
from ..pipeline_core import PipelineContext, Stage
from .sample_column_helpers import (
    cleanup_sample_columns,
    create_sample_columns_from_gt_intelligent,
)

logger = logging.getLogger(__name__)


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
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        return {"custom_annotation", "chunked_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Calculate inheritance patterns."""
        if context.config.get("dataframe_chunked_analysis_complete"):
            logger.debug("Inheritance already computed during chunked analysis")
            return context
        df = context.current_dataframe
        if df is None and not context.config.get("use_chunked_processing"):
            logger.warning("No DataFrame loaded for inheritance analysis")
            return context

        inheritance_mode = context.config.get("inheritance_mode", "simple")

        if context.config.get("use_chunked_processing"):
            logger.info(
                f"Inheritance analysis ({inheritance_mode} mode) "
                "will be applied during chunked processing"
            )
            preserve_details = self._check_if_scoring_needs_details(context)
            if preserve_details:
                logger.debug("Preserving Inheritance_Details for chunked processing scoring")

            context.config["inheritance_analysis_config"] = {
                "inheritance_mode": inheritance_mode,
                "vcf_samples": context.vcf_samples,
                "pedigree_data": context.pedigree_data,
                "use_vectorized_comp_het": not context.config.get("no_vectorized_comp_het", False),
                "preserve_details_for_scoring": preserve_details,
            }
            return context

        vcf_samples = context.vcf_samples or []
        pedigree_data: dict[str, Any] | None = context.pedigree_data

        if not vcf_samples:
            logger.warning("No VCF samples available for inheritance analysis")
            return context

        if isinstance(vcf_samples, set):
            vcf_samples = sorted(vcf_samples)
            logger.warning(
                "VCF samples was a set - converted to sorted list for deterministic order"
            )

        assert df is not None, "DataFrame must not be None for inheritance analysis"

        logger.info(
            f"Calculating inheritance patterns ({inheritance_mode} mode) "
            f"for {len(vcf_samples)} samples"
        )

        sample_columns_exist = any(s in df.columns for s in vcf_samples)

        if not sample_columns_exist and "GT" in df.columns:
            logger.debug("Preparing sample columns from GT field for inheritance analysis")
            prep_start = self._start_subtask("sample_column_preparation")

            sample_column_method = context.config.get("sample_column_creation_method", "auto")
            df = create_sample_columns_from_gt_intelligent(
                df=df,
                vcf_samples=vcf_samples,
                separator=context.config.get("separator", ";"),
                snpsift_sep=context.config.get("extract_fields_separator", ","),
                method=sample_column_method,
            )
            self._end_subtask("sample_column_preparation", prep_start)

        elif not sample_columns_exist:
            gt_cols = _find_per_sample_gt_columns(df)
            if gt_cols:
                n = min(len(gt_cols), len(vcf_samples))
                logger.info(
                    f"Creating sample columns from {n} per-sample GT columns "
                    f"for inheritance analysis"
                )
                prep_start = self._start_subtask("sample_column_preparation")
                for i in range(n):
                    df[vcf_samples[i]] = df[gt_cols[i]]
                self._end_subtask("sample_column_preparation", prep_start)

        analysis_start = self._start_subtask("inheritance_calculation")

        try:
            threads = context.config.get("threads", 1)
            min_variants_for_parallel = context.config.get(
                "min_variants_for_parallel_inheritance", 100
            )

            from . import analysis_stages

            if threads > 1 and len(df) >= min_variants_for_parallel:
                logger.info(f"Using parallel inheritance analyzer with {threads} workers")
                import time

                from ..inheritance.parallel_analyzer import analyze_inheritance_parallel

                comp_het_start = time.time()

                df = analyze_inheritance_parallel(
                    df=df,
                    sample_list=vcf_samples,
                    pedigree_data=pedigree_data or {},
                    use_vectorized_comp_het=not context.config.get("no_vectorized_comp_het", False),
                    n_workers=threads,
                    min_variants_for_parallel=min_variants_for_parallel,
                    reference=context.config.get("reference"),
                )

                comp_het_time = time.time() - comp_het_start
                if comp_het_time > 1.0:
                    self._subtask_times["compound_het_analysis"] = comp_het_time
            else:
                df = analysis_stages.analyze_inheritance(
                    df=df,
                    sample_list=vcf_samples,
                    pedigree_data=pedigree_data or {},
                    use_vectorized_comp_het=not context.config.get("no_vectorized_comp_het", False),
                    reference=context.config.get("reference"),
                )

        except Exception as e:
            logger.error(f"Error in inheritance analysis: {e}")
            raise

        self._end_subtask("inheritance_calculation", analysis_start)

        output_start = self._start_subtask("output_processing")
        from ..inheritance.analyzer import process_inheritance_output

        preserve_details = self._check_if_scoring_needs_details(context)
        if preserve_details:
            logger.debug("Preserving Inheritance_Details for scoring configuration")

        df = process_inheritance_output(
            df, inheritance_mode, preserve_details_for_scoring=preserve_details
        )
        self._end_subtask("output_processing", output_start)

        if vcf_samples:
            cleanup_start = self._start_subtask("column_cleanup")
            df = cleanup_sample_columns(df, vcf_samples, preserve_columns=[])
            self._end_subtask("column_cleanup", cleanup_start)

        context.current_dataframe = df
        return context

    def _check_if_scoring_needs_details(self, context: PipelineContext) -> bool:
        """Check if scoring configuration requires Inheritance_Details column."""
        if not context.scoring_config:
            return False

        if "formulas" in context.scoring_config:
            for formula_dict in context.scoring_config["formulas"]:
                for formula_name, formula_expr in formula_dict.items():
                    if isinstance(formula_expr, str) and any(
                        var in formula_expr for var in ["details", "segregation_p_value"]
                    ):
                        logger.info(
                            f"Scoring formula '{formula_name}' requires Inheritance_Details"
                        )
                        return True

        if "variables" in context.scoring_config:
            variables = context.scoring_config["variables"]
            if "Inheritance_Details" in variables:
                logger.info("Scoring configuration maps Inheritance_Details column")
                return True

        return False

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        input_files = []
        if hasattr(context, "ped_file") and context.ped_file:
            input_files.append(Path(context.ped_file))
        return input_files

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        return []


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
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        return {
            "custom_annotation",
            "inheritance_analysis",
            "chunked_analysis",
            "variant_analysis",
            "variant_identifier",
        }

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply variant scoring."""
        if not context.scoring_config:
            logger.debug("No scoring configuration loaded")
            return context

        if context.config.get("use_chunked_processing") and not context.config.get(
            "dataframe_chunked_analysis_complete"
        ):
            logger.info("Scoring will be applied during chunked processing")
            return context

        df = context.current_dataframe
        if df is None:
            logger.warning("No DataFrame loaded for scoring")
            return context

        from . import analysis_stages

        logger.info("Applying variant scoring formulas")
        df = analysis_stages.apply_scoring(df, context.scoring_config)

        context.current_dataframe = df
        return context
