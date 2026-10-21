"""Association analysis stage using AssociationEngine framework."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from ..association.aggregation import _aggregate_gene_burden_legacy
from ..association.covariate_checks import check_convergence_failures, preflight_covariates
from ..association.engine import AssociationEngine
from ..association.engine_workers import available_cpus, gc_paused
from ..association.gene_sets import add_gene_set_units, load_gene_sets, run_gene_and_set_units
from ..association.gt_codes import aggregate_gene_counts, parse_genotype_codes
from ..association.timing import log_stage_summary, substep
from ..association.weights import resolve_score_weight_column
from ..gene_burden import drop_transcript_duplicates
from ..gt_reconstruction import _find_per_sample_gt_columns
from ..pipeline_core import PipelineContext, Stage
from ..sample_subset import subset_summary
from .analysis_stage_helpers import (
    _build_assoc_config_from_context,
    _count_significant_association_results,
)
from .association_outcome import build_quantitative_outcome
from .association_stage_helpers import (
    _coast_missing_columns,
    _cohort_size_ok,
    _emit_score_artifact,
    _load_covariates_and_pcs,
    _merge_carrier_warnings,
    _report_missingness_qc,
    _score_artifact_hook,
)
from .genotype_matrix_builder import (
    _find_cadd_column,
    _find_effect_column,
    _find_revel_column,
    _GenotypeMatrixBuilder,
)

logger = logging.getLogger(__name__)
_PER_SAMPLE_COLUMN = re.compile(r"^GEN[_\[](\d+)[_\]\.]+([A-Za-z0-9_]+)$")


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

        timings: dict[str, float] = {}  # sub-step seconds, read by benchmarks and progress logs
        context.stage_results[self.name] = {"substep_seconds": timings}

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

        set_ids: list[str] = []
        if context.config.get("gene_set_file"):
            df, set_ids = add_gene_set_units(df, load_gene_sets(context.config["gene_set_file"]))

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
            missing = _coast_missing_columns(df, assoc_config)
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

        with substep(timings, "sample_setup"):
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
                    context.config.get("remove_sample_substring"),
                )
            else:
                assigned_samples = (
                    [s for s in vcf_samples_list if s in case_set or s in control_set]
                    if vcf_samples_list
                    else []
                )

        n_cases_total = len(case_samples)
        n_controls_total = len(control_samples)
        guard_ratio = (
            None if assoc_config.diagnostics_output else assoc_config.max_case_control_ratio
        )
        if not quantitative and not _cohort_size_ok(n_cases_total, n_controls_total, guard_ratio):
            return context

        with substep(timings, "covariates_phenotype"):
            raw_categorical: dict[str, Any] = {}
            covariate_matrix, covariate_col_names = _load_covariates_and_pcs(
                assoc_config,
                assigned_samples,
                context.config.get("remove_sample_substring"),
                raw_categorical,
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
            if needs_regression:  # Fisher ignores covariates; only their consumers are checked
                covariate_matrix, covariate_col_names = preflight_covariates(
                    covariate_matrix, covariate_col_names, phenotype_vector, binary=not quantitative
                )

        gt_cols = _find_per_sample_gt_columns(df)
        # Samples whose genotypes are parsed: the assigned VCF samples, or (packed GT
        # without VCF sample IDs) the case/control lists themselves.
        samples = (
            assigned_samples
            if vcf_samples_list
            else list(dict.fromkeys([*case_samples, *control_samples]))
        )
        # Set-up forks only where an engine loop (genes, or gene sets) does; never with R tests.
        n_workers, n_sets = assoc_config.association_workers, len(set_ids)
        setup_workers = (
            (available_cpus() if n_workers == -1 else n_workers)
            if engine.runs_parallel(df["GENE"].nunique() - n_sets) or engine.runs_parallel(n_sets)
            else 1
        )
        genotypes = None
        if "GT" in df.columns or (gt_cols and vcf_samples_list):
            with substep(timings, "gt_parse"):  # once per run, assigned samples only (#219)
                genotypes = parse_genotype_codes(
                    df, samples, vcf_samples_list, gt_cols, setup_workers
                )
            # Per-sample GEN[i].* columns are no longer needed: keep only annotations.
            # Packed GT too: genotypes come from the codes from here on.
            df = df.drop(
                columns=[c for c in df.columns if c == "GT" or _PER_SAMPLE_COLUMN.match(c)]
            )

        with substep(timings, "aggregation"):
            gene_rows = df.groupby("GENE", observed=True).indices
            if genotypes is not None:
                logger.info(
                    f"Association analysis: genotype-code aggregation "
                    f"({len(case_set)} cases, {len(control_set)} controls)"
                )
                gene_burden_data = aggregate_gene_counts(
                    genotypes,
                    gene_rows,
                    np.array([s in case_set for s in samples], dtype=bool),
                    np.array([s in control_set for s in samples], dtype=bool),
                    assoc_config.missing_as_ref,
                    setup_workers,
                )
            else:
                logger.info(
                    "Association analysis: using pre-computed per-variant counts (legacy mode)"
                )
                gene_burden_data = _aggregate_gene_burden_legacy(df)

        if not gene_burden_data:
            logger.warning("No genes found with variant data for association analysis.")
            return context

        if needs_regression and genotypes is not None and assigned_samples:
            is_binary = assoc_config.trait_type == "binary"
            cadd_column = _find_cadd_column(df) if needs_functional_annotations else None
            revel_column = _find_revel_column(df) if needs_functional_annotations else None
            effect_column = (
                _find_effect_column(df)
                if needs_functional_annotations or score_weight_column is not None
                else None
            )
            with substep(timings, "matrix_builders"), gc_paused():
                # One take in gene order, then a contiguous slice per gene: the same rows
                # and index labels as df.iloc[rows], without a pandas take per gene.
                ordered_rows = [gene_rows[d["GENE"]] for d in gene_burden_data]
                by_gene = df.take(np.concatenate(ordered_rows))
                ends = np.cumsum([len(r) for r in ordered_rows])
                for gene_data, rows, end in zip(gene_burden_data, ordered_rows, ends, strict=True):
                    gene_df = by_gene.iloc[end - len(rows) : end]
                    gene_data["_genotype_matrix_builder"] = _GenotypeMatrixBuilder(
                        gene_df=gene_df,
                        vcf_samples=assigned_samples,
                        gt_columns=[],
                        is_binary=is_binary,
                        missing_site_threshold=assoc_config.missing_site_threshold,
                        missing_sample_threshold=assoc_config.missing_sample_threshold,
                        phenotype_vector=phenotype_vector,
                        covariate_matrix=covariate_matrix,
                        score_column=score_weight_column,
                        cadd_column=cadd_column,
                        revel_column=revel_column,
                        effect_column=effect_column,
                        missing_as_ref=assoc_config.missing_as_ref,
                        genotypes=genotypes,
                        rows=rows,
                    )
                    gene_data["vcf_samples"] = assigned_samples
                    gene_data["gene_df"] = gene_df.reset_index(drop=True)

        gene_data_by_gene = {d.get("GENE", d.get("gene", "")): d for d in gene_burden_data}

        # The score artifact is computed per gene inside the engine's workers, from the
        # matrix built for the tests (one build per gene, #212); builders are saved first
        # because run_all consumes them.
        artifact_builders = {
            d["GENE"]: d["_genotype_matrix_builder"]
            for d in gene_burden_data
            if score_artifact and "_genotype_matrix_builder" in d
        }
        gene_hook = (
            _score_artifact_hook(artifact_builders, assoc_config.skat_robust, df)
            if artifact_builders
            else None
        )

        with substep(timings, "engine"):
            results_df = (
                run_gene_and_set_units(
                    engine, test_names, assoc_config, gene_burden_data, set_ids, gene_hook
                )
                if set_ids
                else engine.run_all(gene_burden_data, gene_hook=gene_hook)
            )
        if score_artifact:
            with substep(timings, "score_artifact"):
                _emit_score_artifact(
                    context.config,
                    score_artifact,
                    artifact_builders,
                    getattr(engine, "gene_hook_results", None),
                    assoc_config.skat_robust,
                )
            if not regression_tests:  # matrices were built only for the artifact
                for d in gene_burden_data:
                    d.pop("n_variants_qc", None)

        if results_df.empty:
            logger.warning("Association analysis produced no results.")
            return context

        _merge_carrier_warnings(results_df, gene_data_by_gene, assoc_config, quantitative)
        check_convergence_failures(
            results_df, covariate_col_names, phenotype_vector, covariate_matrix
        )
        _report_missingness_qc(results_df, gene_data_by_gene, assoc_config)

        assoc_output, compression = self._resolve_association_output(context)

        from ..checkpoint_models import AtomicFileOperation

        with substep(timings, "write_output"), AtomicFileOperation(assoc_output) as temp_out:
            results_df.to_csv(temp_out, sep="\t", index=False, compression=cast(Any, compression))
        logger.info(f"Wrote association results to {assoc_output}")

        if assoc_config.diagnostics_output:
            from ..association.diagnostics import (
                emit_sample_size_warnings,
                write_covariate_balance,
                write_diagnostics,
            )

            cohort_warnings = (
                []
                if quantitative
                else emit_sample_size_warnings(n_cases_total, n_controls_total, assoc_config)
            )
            write_diagnostics(  # QQ/lambda_GC describe the gene family only
                results_df=results_df[results_df["unit_type"] == "gene"] if set_ids else results_df,
                diagnostics_dir=assoc_config.diagnostics_output,
                test_names=test_names,
                n_cases=n_cases_total,
                n_controls=n_controls_total,
                cohort_warnings=cohort_warnings,
            )
            if raw_categorical and not quantitative:
                write_covariate_balance(
                    raw_categorical, case_set, control_set, assoc_config.diagnostics_output
                )

        context.association_results = results_df

        n_genes = len(results_df)
        log_stage_summary(
            timings,
            [
                f"aggregation path: {'genotype-code' if genotypes is not None else 'legacy counts'}"
                f"; genotypes parsed for {len(samples)} assigned samples of "
                f"{len(vcf_samples_list or [])} extracted samples",
                subset_summary(context),
                *getattr(engine, "run_info", []),
            ],
        )
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
