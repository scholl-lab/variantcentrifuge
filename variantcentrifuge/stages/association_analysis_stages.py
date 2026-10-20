"""Association analysis stage using AssociationEngine framework."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from ..association.engine import AssociationEngine
from ..association.weights import resolve_score_weight_column
from ..gene_burden import (
    _aggregate_gene_burden_from_columns,
    _aggregate_gene_burden_from_gt,
    _aggregate_gene_burden_legacy,
    drop_transcript_duplicates,
)
from ..gt_reconstruction import _find_per_sample_gt_columns, reconstruct_gt_column
from ..pipeline_core import PipelineContext, Stage
from .analysis_stage_helpers import (
    _build_assoc_config_from_context,
    _count_significant_association_results,
)
from .association_outcome import build_quantitative_outcome
from .genotype_matrix_builder import (
    _find_cadd_column,
    _find_effect_column,
    _find_revel_column,
    _GenotypeMatrixBuilder,
)

logger = logging.getLogger(__name__)


def _emit_score_artifact(config: dict, path: str, gene_burden_data: list[dict]) -> None:
    """Write the meta-analysis score artifact from this run's per-gene matrix builders."""
    import hashlib
    import json

    from ..association.meta.artifact import emit_stratum_artifact
    from ..links import normalize_reference

    reference = str(config.get("reference") or "")
    # Strata are pooled only when they selected variants the same way.
    selection = {k: config.get(k) for k in ("filters", "final_filter", "bcftools_prefilter")}
    mask_hash = hashlib.sha256(json.dumps(selection, sort_keys=True).encode()).hexdigest()
    emit_stratum_artifact(
        path,
        {
            d["GENE"]: d["_genotype_matrix_builder"]
            for d in gene_burden_data
            if "_genotype_matrix_builder" in d
        },
        stratum_id=Path(path).name.split(".")[0],
        genome_build=normalize_reference(reference) or reference,
        annotation_version=reference,
        mask_id=mask_hash[:12],
        mask_hash=mask_hash,
    )


def _cohort_size_ok(n_cases_total: int, n_controls_total: int) -> bool:
    """Log binary cohort size/balance warnings; False when too few cases to analyse."""
    if n_cases_total < 10:
        logger.error(
            f"Association analysis: only {n_cases_total} case(s) found — "
            "fewer than 10 cases produces invalid results. Aborting."
        )
        return False
    if n_cases_total < 50:
        logger.warning(
            f"Association analysis: only {n_cases_total} case(s) — "
            "fewer than 50 cases provides no practical power for regression tests"
        )
    elif n_cases_total < 200:
        logger.warning(
            f"Association analysis: {n_cases_total} case(s) — "
            "fewer than 200 cases is underpowered for SKAT; interpret results with caution"
        )
    if n_controls_total > 0 and n_cases_total > 0:
        ratio = n_controls_total / n_cases_total
        if ratio > 20:
            logger.warning(
                f"Association analysis: case:control ratio is 1:{ratio:.0f} — "
                "exceeds 1:20; Type I error inflation risk without SPA/Firth correction"
            )
    return True


class AssociationAnalysisStage(Stage):
    """Perform association analysis using the modular AssociationEngine framework."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "association_analysis"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Perform association analysis"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"dataframe_loading", "sample_config_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        return {
            "custom_annotation",
            "pca_computation",
            "chunked_analysis",
            "scoring_config_loading",
            "variant_scoring",
        }

    @property
    def parallel_safe(self) -> bool:
        """Return False — rpy2/R SKAT backend is not thread-safe (SKAT-08)."""
        return False

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system."""
        if not context.config.get("perform_association"):
            return context

        assoc_output = context.config.get("association_output")
        if not assoc_output:
            output_dir = context.config.get("output_dir", "output")
            base_name = context.config.get("output_file_base", "association_results")

            default_path = Path(output_dir) / f"{base_name}.association.tsv"
            default_path_gz = Path(str(default_path) + ".gz")

            if default_path_gz.exists():
                assoc_output = str(default_path_gz)
                context.config["association_output"] = assoc_output
                logger.info(f"Restored association output file path: {assoc_output}")
            elif default_path.exists():
                assoc_output = str(default_path)
                context.config["association_output"] = assoc_output
                logger.info(f"Restored association output file path: {assoc_output}")
            else:
                logger.warning("Association output file not found during checkpoint skip")

        return context

    def _association_output_columns(self, test_names: list[str]) -> list[str]:
        """Return the header schema for empty association output."""
        columns = [
            "gene",
            "n_cases",
            "n_controls",
            "n_variants",
        ]
        for test_name in test_names:
            columns.extend([f"{test_name}_pvalue", f"{test_name}_qvalue"])
        columns.extend(["primary_test", "primary_pvalue", "primary_qvalue", "warnings"])
        return list(dict.fromkeys(columns))

    def _resolve_association_output(self, context: PipelineContext) -> tuple[str, str | None]:
        """Resolve and persist the association sidecar output path."""
        assoc_output = context.config.get("association_output")
        if not assoc_output:
            output_dir = context.config.get("output_dir", "output")
            base_name = context.config.get("output_file_base", "association_results")
            assoc_output = str(Path(output_dir) / f"{base_name}.association.tsv")
            context.config["association_output"] = assoc_output

        use_compression = context.config.get("gzip_intermediates", True)
        if use_compression and not str(assoc_output).endswith(".gz"):
            assoc_output = str(assoc_output) + ".gz"
            context.config["association_output"] = assoc_output
            return assoc_output, "gzip"
        return assoc_output, None

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Perform association analysis."""
        if not context.config.get("perform_association"):
            logger.debug("Association analysis not requested")
            return context

        pca_result = context.get_result("pca_computation")
        if pca_result and isinstance(pca_result, dict):
            pca_file_from_stage = pca_result.get("pca_file")
            if pca_file_from_stage:
                context.config.setdefault("pca_file", pca_file_from_stage)
                logger.debug(f"Association: PCA file from pipeline stage: {pca_file_from_stage}")

        assoc_config = _build_assoc_config_from_context(context)
        logger.info(f"Association analysis: trait type = {assoc_config.trait_type}")
        quantitative = assoc_config.trait_type == "quantitative"

        _json_assoc = context.config.get("association", {}) or {}
        test_names: list[str] = (
            context.config.get("association_tests")
            or _json_assoc.get("association_tests")
            or ["fisher"]
        )

        case_samples = context.config.get("case_samples", [])
        control_samples = context.config.get("control_samples", [])
        if isinstance(case_samples, str):
            case_samples = [s.strip() for s in case_samples.split(",") if s.strip()]
            context.config["case_samples"] = case_samples
        if isinstance(control_samples, str):
            control_samples = [s.strip() for s in control_samples.split(",") if s.strip()]
            context.config["control_samples"] = control_samples

        if quantitative:
            # The outcome is the numeric phenotype column; case/control labels play no role.
            case_samples, control_samples = [], []
            binary_only = sorted({"fisher", "logistic_burden"} & set(test_names))
            if binary_only:
                raise ValueError(
                    f"Association tests {binary_only} need a binary trait; "
                    "use linear_burden/skat for --trait-type quantitative."
                )
            if context.config.get("score_artifact"):
                raise ValueError("--score-artifact supports binary traits only")
        elif not case_samples and not control_samples:
            raise ValueError(
                "Association analysis requested but both case and control cohorts are missing."
            )
        elif not case_samples:
            raise ValueError("Association analysis requested but case cohort is missing.")
        elif not control_samples:
            raise ValueError("Association analysis requested but control cohort is missing.")

        df = context.current_dataframe
        if df is None:
            df = context.variants_df
        if df is None:
            raise RuntimeError("No DataFrame loaded for association analysis")
        df = drop_transcript_duplicates(df)

        if df.empty:
            from ..checkpoint_models import AtomicFileOperation

            assoc_output, compression = self._resolve_association_output(context)
            empty = pd.DataFrame(columns=self._association_output_columns(test_names))
            with AtomicFileOperation(assoc_output) as temp_out:
                empty.to_csv(temp_out, sep="\t", index=False, compression=cast(Any, compression))
            context.association_results = empty
            return context

        engine = AssociationEngine.from_names(test_names, assoc_config)

        if "GENE" not in df.columns:
            logger.error("DataFrame missing required 'GENE' column for association analysis")
            return context

        score_weight_column = resolve_score_weight_column(
            assoc_config.variant_weights,
            assoc_config.variant_weight_column,
        )
        if score_weight_column is not None and score_weight_column not in df.columns:
            raise ValueError(
                f"Requested variant weight column '{score_weight_column}' is missing from "
                "the association input DataFrame."
            )

        needs_functional_annotations = assoc_config.variant_weights in {
            "cadd",
            "revel",
            "combined",
        }

        if "coast" in test_names:
            _coast_base_required = [
                ("ANN_0__EFFECT", "ANN[0].EFFECT"),
                ("ANN_0__IMPACT", "ANN[0].IMPACT"),
            ]
            _coast_model = getattr(assoc_config, "coast_classification", None) or ""
            _is_sift_model = (
                not _coast_model
                or _coast_model.endswith("sift_polyphen")
                or _coast_model == "sift_polyphen"
            )
            if _is_sift_model:
                _coast_base_required.extend(
                    [
                        ("dbNSFP_SIFT_pred", "dbNSFP_SIFT_pred"),
                        ("dbNSFP_Polyphen2_HDIV_pred", "dbNSFP_Polyphen2_HDIV_pred"),
                    ]
                )
            elif _coast_model.endswith("cadd") or _coast_model == "cadd":
                _coast_base_required.extend(
                    [
                        ("dbNSFP_CADD_phred", "dbNSFP_CADD_phred"),
                    ]
                )
            missing = [
                san
                for san, alt in _coast_base_required
                if san not in df.columns and alt not in df.columns
            ]
            if missing:
                logger.error(
                    "COAST allelic series test requires columns %s but they are "
                    "missing from the data. Use --coast-classification to select "
                    "a model, or add required fields to --fields. Skipping COAST.",
                    missing,
                )
                test_names = [t for t in test_names if t != "coast"]
                if not test_names:
                    return context
                engine = AssociationEngine.from_names(test_names, assoc_config)

        score_artifact = context.config.get("score_artifact")
        # The score artifact needs the same per-gene genotype matrices as regression tests.
        regression_tests = any(
            t in test_names
            for t in ("logistic_burden", "linear_burden", "skat", "skat_python", "coast")
        )
        needs_regression = regression_tests or bool(score_artifact)

        df_with_per_sample_gt: pd.DataFrame | None = None
        gt_cols = _find_per_sample_gt_columns(df)
        if needs_regression and gt_cols:
            df_with_per_sample_gt = df
        elif needs_regression and "GT" in df.columns and context.vcf_samples:
            from .sample_column_helpers import create_sample_columns_from_gt_intelligent

            unpacked_samples = list(context.vcf_samples)
            df_with_per_sample_gt = create_sample_columns_from_gt_intelligent(
                df.copy(), vcf_samples=unpacked_samples
            )
            gt_series_str = df_with_per_sample_gt["GT"].fillna("").astype(str)
            for s_col in unpacked_samples:
                if s_col not in df_with_per_sample_gt.columns:
                    df_with_per_sample_gt[s_col] = gt_series_str.str.extract(
                        rf"(?:^|;)\s*{re.escape(s_col)}\(([^:)]+)", expand=False
                    ).fillna("./.")
                listed = gt_series_str.str.contains(rf"(?:^|;)\s*{re.escape(s_col)}\(", regex=True)
                fill_ref = (df_with_per_sample_gt[s_col] == "./.") & ~listed
                df_with_per_sample_gt.loc[fill_ref, s_col] = "0/0"

        if "GT" not in df.columns and context.vcf_samples and gt_cols:
            logger.info("Reconstructing GT column for association aggregation (local copy only)")
            df = reconstruct_gt_column(df.copy(), context.vcf_samples)

        case_set = set(case_samples)
        control_set = set(control_samples)
        vcf_samples_list = list(context.vcf_samples) if context.vcf_samples else None
        quantitative_outcome = None
        if quantitative:
            if not vcf_samples_list:
                raise ValueError("Quantitative association analysis requires VCF sample IDs.")
            assigned_samples, quantitative_outcome = build_quantitative_outcome(
                context.phenotype_data,
                vcf_samples_list,
                context.config.get("phenotype_value_column"),
            )
        else:
            assigned_samples = (
                [s for s in vcf_samples_list if s in case_set or s in control_set]
                if vcf_samples_list
                else []
            )

        n_cases_total = len(case_samples)
        n_controls_total = len(control_samples)
        if not quantitative and not _cohort_size_ok(n_cases_total, n_controls_total):
            return context

        covariate_matrix = None
        covariate_col_names: list[str] = []
        if assoc_config.covariate_file and assigned_samples:
            from ..association.covariates import load_covariates

            covariate_matrix, covariate_col_names = load_covariates(
                assoc_config.covariate_file,
                assigned_samples,
                assoc_config.covariate_columns,
                assoc_config.categorical_covariates,
            )
            logger.info(
                f"Association analysis: loaded {covariate_matrix.shape[1]} "
                f"covariate(s): {covariate_col_names}"
            )

        pca_file = assoc_config.pca_file
        if pca_file and assigned_samples:
            from ..association.pca import load_pca_file, merge_pca_covariates

            n_pcs = assoc_config.pca_components
            pca_matrix, pca_col_names = load_pca_file(
                pca_file,
                assigned_samples,
                n_components=n_pcs,
                remove_sample_substring=context.config.get("remove_sample_substring"),
            )
            covariate_matrix, covariate_col_names = merge_pca_covariates(
                pca_matrix,
                pca_col_names,
                covariate_matrix,
                covariate_col_names,
            )
            logger.info(
                f"PCA: merged {len(pca_col_names)} PCs with existing covariates "
                f"-> {len(covariate_col_names)} total columns"
            )

        phenotype_vector = None
        if quantitative_outcome is not None:
            phenotype_vector = quantitative_outcome
            logger.info(
                f"Association analysis: quantitative outcome "
                f"'{context.config.get('phenotype_value_column')}' for "
                f"{len(assigned_samples)} samples"
            )
        elif assigned_samples:
            phenotype_vector = np.array(
                [1.0 if s in case_set else 0.0 for s in assigned_samples],
                dtype=float,
            )
            n_pv_cases = int(phenotype_vector.sum())
            n_pv_controls = int((1.0 - phenotype_vector).sum())
            logger.info(
                f"Association analysis: phenotype vector built — "
                f"{n_pv_cases} cases, {n_pv_controls} controls"
            )

        has_case_ctrl = True
        gt_columns = _find_per_sample_gt_columns(df)
        use_column_aggregation = bool(
            has_case_ctrl
            and gt_columns
            and vcf_samples_list
            and len(gt_columns) <= len(vcf_samples_list)
        )
        use_gt_aggregation = has_case_ctrl and not use_column_aggregation and "GT" in df.columns

        if use_column_aggregation:
            logger.info(
                f"Association analysis: using column-based aggregation "
                f"({len(gt_columns)} GT columns, "
                f"{len(case_set)} cases, {len(control_set)} controls)"
            )
            assert vcf_samples_list is not None
            gene_burden_data = _aggregate_gene_burden_from_columns(
                df, case_set, control_set, vcf_samples_list, gt_columns
            )
        elif use_gt_aggregation:
            logger.info(
                "Association analysis: using packed GT string collapsing "
                f"({len(case_set)} cases, {len(control_set)} controls)"
            )
            gene_burden_data = _aggregate_gene_burden_from_gt(df, case_set, control_set)
        else:
            logger.info("Association analysis: using pre-computed per-variant counts (legacy mode)")
            gene_burden_data = _aggregate_gene_burden_legacy(df)

        if not gene_burden_data:
            logger.warning("No genes found with variant data for association analysis.")
            return context

        gt_source_df = df_with_per_sample_gt if df_with_per_sample_gt is not None else df
        gt_columns_for_matrix = _find_per_sample_gt_columns(gt_source_df)
        if not gt_columns_for_matrix and assigned_samples:
            gt_columns_for_matrix = [s for s in assigned_samples if s in gt_source_df.columns]
        if needs_regression and gt_columns_for_matrix and assigned_samples:
            is_binary = assoc_config.trait_type == "binary"
            _gene_groups = gt_source_df.groupby("GENE")

            if vcf_samples_list and len(gt_columns_for_matrix) == len(vcf_samples_list):
                sample_to_gt = dict(zip(vcf_samples_list, gt_columns_for_matrix, strict=True))
                assigned_gt_cols = [sample_to_gt[s] for s in assigned_samples if s in sample_to_gt]
            else:
                assigned_gt_cols = [c for c in assigned_samples if c in gt_source_df.columns]
                if not assigned_gt_cols:
                    assigned_gt_cols = gt_columns_for_matrix

            for gene_data in gene_burden_data:
                gene_name = gene_data.get("GENE", "")
                try:
                    gene_df = _gene_groups.get_group(gene_name)
                except KeyError:
                    gene_df = gt_source_df.iloc[0:0]

                cadd_column = _find_cadd_column(gene_df) if needs_functional_annotations else None
                revel_column = _find_revel_column(gene_df) if needs_functional_annotations else None
                effect_column = (
                    _find_effect_column(gene_df)
                    if needs_functional_annotations or score_weight_column is not None
                    else None
                )

                builder = _GenotypeMatrixBuilder(
                    gene_df=gene_df,
                    vcf_samples=assigned_samples,
                    gt_columns=assigned_gt_cols,
                    is_binary=is_binary,
                    missing_site_threshold=assoc_config.missing_site_threshold,
                    missing_sample_threshold=assoc_config.missing_sample_threshold,
                    phenotype_vector=phenotype_vector,
                    covariate_matrix=covariate_matrix,
                    score_column=score_weight_column,
                    cadd_column=cadd_column,
                    revel_column=revel_column,
                    effect_column=effect_column,
                )
                gene_data["_genotype_matrix_builder"] = builder
                gene_data["vcf_samples"] = assigned_samples
                gene_data["gene_df"] = gene_df.reset_index(drop=True)

            del df_with_per_sample_gt, gt_source_df

        gene_data_by_gene = {d.get("GENE", d.get("gene", "")): d for d in gene_burden_data}

        if score_artifact:  # before run_all, which consumes the matrix builders
            _emit_score_artifact(context.config, score_artifact, gene_burden_data)
            if not regression_tests:  # don't make the engine build matrices nobody tests
                for d in gene_burden_data:
                    d.pop("_genotype_matrix_builder", None)

        results_df = engine.run_all(gene_burden_data)

        if results_df.empty:
            logger.warning("Association analysis produced no results.")
            return context

        from ..association.diagnostics import compute_per_gene_warnings

        warnings_by_gene: dict[str, str] = {}
        # Case-carrier warnings are a binary-trait concept.
        for gene in [] if quantitative else results_df["gene"].unique():
            gdata = gene_data_by_gene.get(gene, {})
            case_carriers = gdata.get("proband_carrier_count", 0)
            gene_warnings = compute_per_gene_warnings(gene, case_carriers, assoc_config)
            if gene_warnings:
                warnings_by_gene[gene] = ";".join(gene_warnings)
        results_df["warnings"] = results_df["gene"].map(warnings_by_gene).fillna("")

        assoc_output, compression = self._resolve_association_output(context)

        from ..checkpoint_models import AtomicFileOperation

        with AtomicFileOperation(assoc_output) as temp_out:
            results_df.to_csv(temp_out, sep="\t", index=False, compression=cast(Any, compression))
        logger.info(f"Wrote association results to {assoc_output}")

        if assoc_config.diagnostics_output:
            from ..association.diagnostics import emit_sample_size_warnings, write_diagnostics

            cohort_warnings = (
                []
                if quantitative
                else emit_sample_size_warnings(n_cases_total, n_controls_total, assoc_config)
            )
            write_diagnostics(
                results_df=results_df,
                diagnostics_dir=assoc_config.diagnostics_output,
                test_names=test_names,
                n_cases=n_cases_total,
                n_controls=n_controls_total,
                cohort_warnings=cohort_warnings,
            )

        context.association_results = results_df

        n_genes = len(results_df)
        n_sig = _count_significant_association_results(results_df)
        logger.info(
            f"Association analysis: {n_genes} genes tested, {n_sig} significant (FDR < 0.05)"
        )

        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        if hasattr(context, "config") and context.config.get("association_output"):
            assoc_output = context.config["association_output"]
            if Path(assoc_output).exists():
                return [Path(assoc_output)]
        return []
