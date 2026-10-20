"""Variant-level analysis stage."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..analyze_variants import analyze_variants_df
from ..gt_reconstruction import _find_per_sample_gt_columns, reconstruct_gt_column
from ..pipeline_core import PipelineContext, Stage

logger = logging.getLogger(__name__)


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
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        return {"custom_annotation", "inheritance_analysis", "chunked_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Run variant analysis."""
        if context.current_dataframe is None:
            logger.warning("No DataFrame available for variant analysis")
            return context

        df = context.current_dataframe
        should_have_inheritance = self._should_calculate_inheritance(context)

        inheritance_cols_present = [
            col for col in df.columns if col in ["Inheritance_Pattern", "Inheritance_Details"]
        ]
        if inheritance_cols_present:
            logger.debug(f"Input DataFrame has inheritance columns: {inheritance_cols_present}")
            for col in inheritance_cols_present:
                non_null_count = df[col].notna().sum()
                logger.debug(f"Column '{col}' has {non_null_count} non-null values")
        else:
            if should_have_inheritance:
                logger.warning(
                    "No inheritance columns found in input DataFrame - "
                    "they were lost before VariantAnalysisStage"
                )
            else:
                logger.debug(
                    "No inheritance columns found - inheritance analysis was not requested"
                )

        analysis_config: dict[str, Any] = {
            "vcf_sample_names": context.vcf_samples or [],
            "sample_list": ",".join(context.vcf_samples or []),
            "no_db_links": context.config.get("no_db_links", False),
            "base_name": context.workspace.base_name,
            "reference": context.config.get("reference", "GRCh37"),
            "case_samples": context.config.get("case_samples") or [],
            "control_samples": context.config.get("control_samples") or [],
            "case_phenotypes": context.config.get("case_phenotypes") or [],
            "control_phenotypes": context.config.get("control_phenotypes") or [],
            "no_stats": True,
        }
        _original_gt_cols = _find_per_sample_gt_columns(df)
        if "GT" not in df.columns and context.vcf_samples and _original_gt_cols:
            logger.info("Reconstructing GT column for variant analysis (local copy only)")
            reconstructed = reconstruct_gt_column(df.copy(deep=False), context.vcf_samples)
            df = df.copy(deep=False)
            df["GT"] = reconstructed["GT"].values

        df_to_analyze = df
        if context.column_rename_map:
            reverse_map = {v: k for k, v in context.column_rename_map.items()}
            df_to_analyze = df.rename(columns=reverse_map, copy=False)

        logger.info("Running variant-level analysis")
        analysis_df = analyze_variants_df(df_to_analyze, analysis_config)

        if analysis_df is not None:
            if context.column_rename_map:
                analysis_df = analysis_df.rename(columns=context.column_rename_map, copy=False)

            key_cols = ["CHROM", "POS", "REF", "ALT"] + [
                c
                for c in ("FEATUREID", "ANN[0].FEATUREID", "ANN_0__FEATUREID")
                if c in df.columns and c in analysis_df.columns
            ]

            orig_key_cols_present = [col for col in key_cols if col in df.columns]
            analysis_key_cols_present = [col for col in key_cols if col in analysis_df.columns]
            logger.debug(f"Original DataFrame key columns: {orig_key_cols_present}")
            logger.debug(f"Analysis DataFrame key columns: {analysis_key_cols_present}")

            if all(col in df.columns for col in key_cols) and all(
                col in analysis_df.columns for col in key_cols
            ):
                original_only_cols = [col for col in df.columns if col not in analysis_df.columns]

                if original_only_cols:
                    inheritance_cols = [
                        col
                        for col in original_only_cols
                        if col in ["Inheritance_Pattern", "Inheritance_Details"]
                    ]
                    if inheritance_cols:
                        logger.debug(f"Preserving inheritance columns: {inheritance_cols}")

                    try:
                        if len(analysis_df) == len(df):
                            for col in original_only_cols:
                                analysis_df[col] = df[col].values
                        else:
                            merge_df = df[key_cols + original_only_cols].copy()
                            analysis_df = pd.merge(analysis_df, merge_df, on=key_cols, how="left")

                        logger.debug(
                            f"Successfully preserved columns from original dataframe: "
                            f"{original_only_cols}"
                        )

                        if inheritance_cols:
                            preserved_inheritance = [
                                col for col in inheritance_cols if col in analysis_df.columns
                            ]
                            if preserved_inheritance:
                                logger.debug(
                                    f"Inheritance columns successfully preserved: "
                                    f"{preserved_inheritance}"
                                )
                                for col in preserved_inheritance:
                                    non_null_count = analysis_df[col].notna().sum()
                                    logger.debug(
                                        f"Column '{col}' has {non_null_count} non-null values"
                                    )
                            else:
                                msg = f"Inheritance columns: {inheritance_cols}"
                                if should_have_inheritance:
                                    logger.error(f"Failed to preserve {msg}")
                                else:
                                    logger.debug(f"No {msg} to preserve")
                    except Exception as e:
                        logger.error(f"Error preserving columns: {e}")
                        logger.warning(
                            "Inheritance columns may be lost due to preservation failure"
                        )

                    if "VAR_ID" in analysis_df.columns:
                        cols = ["VAR_ID", *[c for c in analysis_df.columns if c != "VAR_ID"]]
                        analysis_df = analysis_df[cols]

                    if "Custom_Annotation" in analysis_df.columns and "GT" in analysis_df.columns:
                        cols = analysis_df.columns.tolist()
                        cols.remove("Custom_Annotation")
                        gt_idx = cols.index("GT")
                        cols.insert(gt_idx + 1, "Custom_Annotation")
                        analysis_df = analysis_df[cols]

            else:
                logger.warning(
                    "Could not preserve columns from original dataframe - missing key columns "
                    "for merging"
                )
                inheritance_cols = [
                    col
                    for col in df.columns
                    if col in ["Inheritance_Pattern", "Inheritance_Details"]
                ]
                if inheritance_cols and should_have_inheritance:
                    logger.error(
                        f"CRITICAL: Inheritance columns {inheritance_cols} will be lost "
                        f"due to missing key columns for merge"
                    )

            context.current_dataframe = analysis_df
            logger.info(f"Variant analysis complete: {len(analysis_df)} variants analyzed")

            final_inheritance_cols = [
                c
                for c in analysis_df.columns
                if c in ("Inheritance_Pattern", "Inheritance_Details")
            ]
            if final_inheritance_cols:
                logger.debug(f"Final DataFrame has inheritance columns: {final_inheritance_cols}")
                for col in final_inheritance_cols:
                    n_vals = analysis_df[col].notna().sum()
                    logger.debug("Column '%s' has %d non-null values", col, n_vals)
            else:
                if should_have_inheritance:
                    logger.error(
                        "INHERITANCE COLUMNS LOST - Final DataFrame missing inheritance columns!"
                    )
                else:
                    logger.debug("Final DataFrame has no inheritance columns - none were expected")
        else:
            logger.warning("No analysis results generated - DataFrame unchanged")
            inheritance_cols_present = [
                col for col in df.columns if col in ["Inheritance_Pattern", "Inheritance_Details"]
            ]
            if inheritance_cols_present:
                logger.debug(
                    f"Unchanged DataFrame still has inheritance columns: {inheritance_cols_present}"
                )
            else:
                if should_have_inheritance:
                    logger.warning("Unchanged DataFrame also missing inheritance columns")
                else:
                    logger.debug(
                        "Unchanged DataFrame has no inheritance columns - none were expected"
                    )

        return context

    def _should_calculate_inheritance(self, context: PipelineContext) -> bool:
        """Check if inheritance analysis should be calculated."""
        has_ped_file = context.pedigree_data is not None and len(context.pedigree_data) > 0
        has_inheritance_mode = context.config.get("inheritance_mode")
        has_calculate_inheritance_config = context.config.get("calculate_inheritance", False)
        requires_inheritance_for_scoring = self._check_if_scoring_needs_details(context)

        should_calculate = (
            has_ped_file
            or has_inheritance_mode
            or has_calculate_inheritance_config
            or requires_inheritance_for_scoring
        )

        logger.debug(
            f"Should calculate inheritance: {should_calculate} (ped={has_ped_file}, "
            f"mode={has_inheritance_mode}, config={has_calculate_inheritance_config}, "
            f"scoring={requires_inheritance_for_scoring})"
        )
        return bool(should_calculate)

    def _check_if_scoring_needs_details(self, context: PipelineContext) -> bool:
        """Check if scoring configuration requires Inheritance_Details."""
        if not context.scoring_config:
            return False

        scoring_formulas = context.scoring_config.get("scoring_formulas", {})
        for formula_name, formula_data in scoring_formulas.items():
            formula_str = formula_data.get("formula", "")
            if "Inheritance_Details" in formula_str:
                logger.debug(f"Scoring formula '{formula_name}' requires Inheritance_Details")
                return True

        return False

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        return []
