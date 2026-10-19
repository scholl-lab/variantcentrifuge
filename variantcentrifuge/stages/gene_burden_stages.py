"""Gene burden analysis stage."""

from __future__ import annotations

import logging
import time as _time
from pathlib import Path
from typing import Any, cast

import pandas as pd

from ..gt_reconstruction import _find_per_sample_gt_columns, reconstruct_gt_column
from ..pipeline_core import PipelineContext, Stage

logger = logging.getLogger(__name__)


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
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"dataframe_loading", "sample_config_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        return {"custom_annotation", "chunked_analysis"}

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system."""
        if not context.config.get("perform_gene_burden"):
            return context

        burden_output = context.config.get("gene_burden_output")
        if not burden_output:
            output_dir = context.config.get("output_dir", "output")
            base_name = context.config.get("output_file_base", "gene_burden_results")

            default_path = Path(output_dir) / f"{base_name}.gene_burden.tsv"
            default_path_gz = Path(str(default_path) + ".gz")

            if default_path_gz.exists():
                burden_output = str(default_path_gz)
                context.config["gene_burden_output"] = burden_output
                logger.info(f"Restored gene burden file path: {burden_output}")
            elif default_path.exists():
                burden_output = str(default_path)
                context.config["gene_burden_output"] = burden_output
                logger.info(f"Restored gene burden file path: {burden_output}")
            else:
                logger.warning("Gene burden output file not found during checkpoint skip")

        return context

    def _resolve_gene_burden_output(self, context: PipelineContext) -> tuple[str, str | None]:
        """Resolve and persist the gene burden sidecar output path."""
        burden_output = context.config.get("gene_burden_output")
        if not burden_output:
            output_dir = context.config.get("output_dir", "output")
            output_file = context.config.get("output_file")
            if output_file and output_file != "-":
                output_path = Path(output_file)
                if not output_path.is_absolute():
                    output_path = Path(output_dir) / output_path
                name = output_path.name
                if name.endswith(".tsv.gz"):
                    base_name = name[: -len(".tsv.gz")]
                elif name.endswith(".tsv"):
                    base_name = name[: -len(".tsv")]
                else:
                    base_name = output_path.stem
                burden_output = str(output_path.with_name(f"{base_name}.gene_burden.tsv"))
            else:
                base_name = context.config.get("output_file_base", "gene_burden_results")
                burden_output = str(Path(output_dir) / f"{base_name}.gene_burden.tsv")
            context.config["gene_burden_output"] = burden_output

        compression = "gzip" if str(burden_output).endswith(".gz") else None
        return burden_output, compression

    def _gene_burden_output_columns(self, context: PipelineContext) -> list[str]:
        """Return the header schema for empty gene burden output."""
        mode = context.config.get("gene_burden_mode", "alleles")
        if mode == "samples":
            return [
                "GENE",
                "proband_count",
                "control_count",
                "proband_carrier_count",
                "control_carrier_count",
                "raw_p_value",
                "corrected_p_value",
                "odds_ratio",
                "or_ci_lower",
                "or_ci_upper",
                "n_qualifying_variants",
            ]
        return [
            "GENE",
            "proband_count",
            "control_count",
            "proband_allele_count",
            "control_allele_count",
            "raw_p_value",
            "corrected_p_value",
            "odds_ratio",
            "or_ci_lower",
            "or_ci_upper",
        ]

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Perform gene burden analysis."""
        if not context.config.get("perform_gene_burden"):
            logger.debug("Gene burden analysis not requested")
            return context

        if context.config.get("use_chunked_processing") and not context.config.get(
            "dataframe_chunked_analysis_complete"
        ):
            logger.info("Gene burden analysis will be performed after chunked processing completes")
            return context

        case_samples = context.config.get("case_samples", [])
        control_samples = context.config.get("control_samples", [])

        logger.debug(
            f"Gene burden analysis context check: "
            f"case_samples={len(case_samples) if case_samples else 0}, "
            f"control_samples={len(control_samples) if control_samples else 0}"
        )

        config_keys = list(context.config.keys())
        sample_related_keys = [k for k in config_keys if "sample" in k.lower()]
        logger.debug(f"All sample-related config keys: {sample_related_keys}")

        if not case_samples and not control_samples:
            raise ValueError(
                "Gene burden analysis requested but both case and control cohorts are missing."
            )
        if not case_samples:
            raise ValueError("Gene burden analysis requested but case cohort is missing.")
        if not control_samples:
            raise ValueError("Gene burden analysis requested but control cohort is missing.")

        df = context.current_dataframe
        if df is None:
            raise RuntimeError("No DataFrame loaded for gene burden analysis")
        if df.empty:
            burden_output, compression = self._resolve_gene_burden_output(context)
            empty = pd.DataFrame(columns=self._gene_burden_output_columns(context))
            empty.to_csv(burden_output, sep="\t", index=False, compression=cast(Any, compression))
            context.gene_burden_results = empty
            return context

        if "GENE" not in df.columns:
            logger.error("DataFrame missing required 'GENE' column for gene burden analysis")
            logger.debug(f"Available columns: {list(df.columns)[:20]}...")
            return context

        df_for_burden = df
        if "GT" not in df.columns and context.vcf_samples:
            gt_cols = _find_per_sample_gt_columns(df)
            if gt_cols:
                logger.info("Reconstructing GT column for gene burden analysis (local copy only)")
                df_for_burden = reconstruct_gt_column(df.copy(), context.vcf_samples)

        if "GT" not in df_for_burden.columns:
            logger.error("DataFrame missing required 'GT' column for gene burden analysis")
            return context

        logger.debug(
            f"Performing gene burden analysis: {len(case_samples)} cases, "
            f"{len(control_samples)} controls"
        )

        burden_config = {
            "gene_burden_mode": context.config.get("gene_burden_mode", "alleles"),
            "correction_method": context.config.get("correction_method", "fdr"),
            "confidence_interval_method": context.config.get(
                "confidence_interval_method", "normal_approx"
            ),
            "confidence_interval_alpha": context.config.get("confidence_interval_alpha", 0.05),
        }

        all_vcf_samples = set(context.vcf_samples) if context.vcf_samples else set()

        logger.info(f"Using {len(all_vcf_samples)} VCF samples for gene burden analysis")

        from ..helpers import assign_case_control_counts

        _t_cc = _time.monotonic()
        df_with_counts = assign_case_control_counts(
            df=df_for_burden,
            case_samples=set(case_samples),
            control_samples=set(control_samples),
            all_samples=all_vcf_samples,
        )
        _t_cc_elapsed = _time.monotonic() - _t_cc
        logger.info(
            f"Case/control count assignment completed in {_t_cc_elapsed:.2f}s "
            f"({len(df_with_counts)} variants, {len(all_vcf_samples)} samples)"
        )

        _t_burden = _time.monotonic()
        from . import analysis_stages

        burden_results = analysis_stages.perform_gene_burden_analysis(
            df=df_with_counts,
            cfg=burden_config,
            case_samples=set(case_samples),
            control_samples=set(control_samples),
            vcf_samples=list(context.vcf_samples) if context.vcf_samples else None,
        )
        _t_burden_elapsed = _time.monotonic() - _t_burden
        n_genes = len(burden_results) if burden_results is not None else 0
        logger.info(
            f"Gene burden analysis completed in {_t_burden_elapsed:.2f}s ({n_genes} genes tested)"
        )

        context.gene_burden_results = burden_results

        burden_output, compression = self._resolve_gene_burden_output(context)

        burden_results.to_csv(
            burden_output, sep="\t", index=False, compression=cast(Any, compression)
        )
        logger.info(f"Wrote gene burden results to {burden_output}")

        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        if hasattr(context, "config") and context.config.get("gene_burden_output"):
            burden_output = context.config["gene_burden_output"]
            if Path(burden_output).exists():
                return [Path(burden_output)]
        return []
