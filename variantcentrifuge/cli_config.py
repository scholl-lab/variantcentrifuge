"""Configuration post-processing and validation helpers for VariantCentrifuge CLI."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

from .cli_memory import calculate_sort_memory_limit as _calculate_sort_memory_limit
from .config import normalize_annotation_config
from .vcf_header_parser import resolve_tumor_normal_indices
from .version import __version__

logger = logging.getLogger("variantcentrifuge")


def _expand_tumor_normal_vars(expr: str, tn_vars: dict[str, str]) -> str:
    """Replace supported tumor/normal placeholders while leaving other braces intact."""
    for key, val in tn_vars.items():
        expr = expr.replace(f"{{{key}}}", val)
    return expr


def resolve_presets_and_tumor_normal(
    args: argparse.Namespace,
    cfg: dict[str, Any],
    filters: str | None,
) -> tuple[str | None, dict[str, str]]:
    """Resolve tumor-normal sample indices and expand preset filter templates."""
    preset_val = cfg.get("presets")
    preset_dict = preset_val if isinstance(preset_val, dict) else {}

    if not preset_dict:
        from .config import load_config as _load_default_config
        from .field_profile import resolve_profile

        try:
            default_cfg = _load_default_config()
            profile_name = getattr(args, "field_profile", None) or cfg.get("field_profile")
            if profile_name:
                default_cfg["field_profile"] = profile_name
            resolve_profile(default_cfg)
            preset_dict = default_cfg.get("presets", {})
        except Exception:
            preset_dict = {}

    raw_presets = getattr(args, "preset", None) or cfg.get("preset")
    if not raw_presets and isinstance(preset_val, (list, str)):
        raw_presets = preset_val
    raw_presets = [raw_presets] if isinstance(raw_presets, str) else (raw_presets or [])
    raw_presets = [p.strip() for item in raw_presets for p in str(item).split(",") if p.strip()]
    raw_final_filter = (
        args.final_filter if args.final_filter is not None else cfg.get("final_filter")
    )
    expr_texts = [preset_dict.get(p.strip(), "") for p in raw_presets]
    if filters:
        expr_texts.append(filters)
    if raw_final_filter:
        expr_texts.append(raw_final_filter)
    tumor_idx_arg = getattr(args, "tumor_sample_index", None)
    normal_idx_arg = getattr(args, "normal_sample_index", None)
    vcf_arg = getattr(args, "vcf_file", None)
    needs_tumor_idx = bool(
        tumor_idx_arg is not None
        or cfg.get("tumor_sample_index") is not None
        or any("{tumor_idx}" in e for e in expr_texts)
    )
    needs_normal_idx = bool(
        normal_idx_arg is not None
        or cfg.get("normal_sample_index") is not None
        or any("{normal_idx}" in e for e in expr_texts)
    )
    needs_tumor_normal = needs_tumor_idx or needs_normal_idx
    t_idx, n_idx = resolve_tumor_normal_indices(
        vcf_arg if needs_tumor_normal else None,
        cli_tumor_idx=tumor_idx_arg,
        cli_normal_idx=normal_idx_arg,
        cfg_tumor_idx=cfg.get("tumor_sample_index"),
        cfg_normal_idx=cfg.get("normal_sample_index"),
        require_tumor=needs_tumor_idx,
        require_normal=needs_normal_idx,
    )
    tumor_normal_vars = {"tumor_idx": t_idx, "normal_idx": n_idx}
    for _k, _d in (
        ("tumor_dp_min", 20),
        ("normal_dp_min", 20),
        ("tumor_af_min", 0.05),
        ("normal_af_max", 0.03),
    ):
        _v = getattr(args, _k, None)
        tumor_normal_vars[_k] = str(cfg.get(_k, _d) if _v is None else _v)

    # Handle presets (from CLI --preset or config "preset")
    if raw_presets:
        cfg["preset"] = [p.strip() for p in raw_presets if p.strip()]
        chosen_presets = []
        for p in raw_presets:
            p_str = p.strip()
            if p_str in preset_dict:
                expr = _expand_tumor_normal_vars(preset_dict[p_str], tumor_normal_vars)
                chosen_presets.append(f"({expr})")
            else:
                msg = f"Preset '{p_str}' not found in config."
                logger.error(msg)
                raise ValueError(msg)
        if chosen_presets:
            combined_preset_filter = " & ".join(chosen_presets)
            filters = (
                f"({combined_preset_filter}) & ({filters})" if filters else combined_preset_filter
            )
    if filters and ("{tumor_" in filters or "{normal_" in filters):
        filters = _expand_tumor_normal_vars(filters, tumor_normal_vars)

    return filters, tumor_normal_vars


def populate_cli_config(
    args: argparse.Namespace,
    cfg: dict[str, Any],
    parser: argparse.ArgumentParser | None,
    reference: str | None,
    filters: str | None,
    fields: str | None,
    tumor_normal_vars: dict[str, str] | None = None,
) -> None:
    """Populate and validate pipeline configuration from CLI arguments."""

    def _fail(msg: str) -> None:
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    tumor_normal_vars = tumor_normal_vars or {}
    cfg["perform_gene_burden"] = bool(
        getattr(args, "perform_gene_burden", False) or cfg.get("perform_gene_burden", False)
    )
    args.perform_gene_burden = cfg["perform_gene_burden"]
    cfg["no_stats"] = bool(getattr(args, "no_stats", False) or cfg.get("no_stats", False))
    args.no_stats = cfg["no_stats"]
    cfg["gene_burden_mode"] = (
        getattr(args, "gene_burden_mode", None) or cfg.get("gene_burden_mode") or "alleles"
    )
    cfg["correction_method"] = (
        getattr(args, "correction_method", None) or cfg.get("correction_method") or "fdr"
    )

    # Association analysis configuration
    cfg["perform_association"] = bool(
        getattr(args, "perform_association", False) or cfg.get("perform_association", False)
    )
    args.perform_association = cfg["perform_association"]
    assoc_tests = getattr(args, "association_tests", None)
    if assoc_tests:
        if isinstance(assoc_tests, str):
            cfg["association_tests"] = [t.strip() for t in assoc_tests.split(",")]
        else:
            cfg["association_tests"] = list(assoc_tests)
    elif cfg.get("association_tests"):
        if isinstance(cfg["association_tests"], str):
            cfg["association_tests"] = [t.strip() for t in cfg["association_tests"].split(",")]
    else:
        cfg["association_tests"] = ["fisher"] if cfg["perform_association"] else []
    args.association_tests = cfg["association_tests"]
    cfg["skat_backend"] = getattr(args, "skat_backend", "python")
    cfg["coast_backend"] = getattr(args, "coast_backend", "python")
    cfg["coast_classification"] = None

    # Phase 19: Covariate and regression configuration
    cfg["covariate_file"] = getattr(args, "covariate_file", None)
    covariates_arg = getattr(args, "covariates", None)
    cfg["covariate_columns"] = (
        [c.strip() for c in covariates_arg.split(",")] if covariates_arg else None
    )
    cat_cov_arg = getattr(args, "categorical_covariates", None)
    cfg["categorical_covariates"] = (
        [c.strip() for c in cat_cov_arg.split(",")] if cat_cov_arg else None
    )
    cfg["trait_type"] = getattr(args, "trait_type", "binary")
    cfg["variant_weights"] = getattr(args, "variant_weights", "beta:1,25")
    cfg["variant_weight_column"] = getattr(args, "variant_weight_column", None)
    # Phase 23: Parse --variant-weight-params JSON string
    _vwp_raw = getattr(args, "variant_weight_params", None)
    if _vwp_raw:
        import json as _json

        try:
            parsed_params = _json.loads(_vwp_raw)
            if not isinstance(parsed_params, dict):
                raise TypeError(f"expected a JSON object, got {type(parsed_params).__name__}")
            cfg["variant_weight_params"] = parsed_params
        except (ValueError, TypeError) as _e:
            _fail(f"--variant-weight-params: invalid JSON: {_e}")
    else:
        cfg["variant_weight_params"] = None
    cfg["diagnostics_output"] = getattr(args, "diagnostics_output", None)
    cfg["score_artifact"] = getattr(args, "score_artifact", None)
    # Phase 32: Region restriction
    cfg["regions_bed"] = getattr(args, "regions_bed", None)
    # Phase 32: Unified PCA configuration (--pca replaces --pca-file/--pca-tool)
    cfg["pca"] = getattr(args, "pca", None)
    cfg["pca_components"] = getattr(args, "pca_components", 10)
    cfg["pca_sites"] = getattr(args, "pca_sites", None)
    # Phase 23: COAST weights — parse comma-separated floats
    _coast_weights_raw = getattr(args, "coast_weights", None)
    if _coast_weights_raw:
        try:
            _coast_weights_parsed = [float(x) for x in _coast_weights_raw.split(",")]
            if len(_coast_weights_parsed) != 3:
                _fail("--coast-weights must have exactly 3 values (BMV, DMV, PTV)")
            cfg["coast_weights"] = _coast_weights_parsed
        except ValueError as _e:
            _fail(f"--coast-weights values must be numeric: {_e}")
    else:
        cfg["coast_weights"] = None
    # Association workers: 0/None = auto-allocate from --threads (default)
    _aw = getattr(args, "association_workers", 0)
    cfg["association_workers"] = _aw if _aw else 0
    # Phase 28: SKAT method and diagnostic thresholds
    cfg["skat_method"] = getattr(args, "skat_method", "SKAT")
    # Phase 33: Gene-level FDR weighting
    _gpw = getattr(args, "gene_prior_weights", None)
    if _gpw and not Path(_gpw).is_file():
        _fail(f"--gene-prior-weights file not found: {_gpw}")
    cfg["gene_prior_weights"] = _gpw
    cfg["gene_prior_weight_column"] = getattr(args, "gene_prior_weight_column", "weight")
    cfg["association_min_cases"] = getattr(args, "min_cases", 200)
    cfg["association_max_case_control_ratio"] = getattr(args, "max_case_control_ratio", 20.0)
    cfg["association_min_case_carriers"] = getattr(args, "min_case_carriers", 10)

    # Handle add_chr configuration
    if args.add_chr:
        cfg["add_chr"] = True
    # If flag is not provided, use the value from config file (already loaded)

    # Tumor-normal parameters (store resolved values for downstream use)
    cfg["tumor_sample_index"] = tumor_normal_vars["tumor_idx"]
    cfg["normal_sample_index"] = tumor_normal_vars["normal_idx"]

    # IGV parameters
    cfg["igv_enabled"] = bool(args.igv or cfg.get("igv_enabled", False))
    for opt in (
        "bam_mapping_file",
        "igv_reference",
        "igv_fasta",
        "igv_ideogram",
        "igv_max_allele_len_filename",
        "igv_hash_len_filename",
        "igv_max_variant_part_filename",
    ):
        v = getattr(args, opt, None) or cfg.get(opt)
        cfg[opt] = v
        if v is not None:
            setattr(args, opt, v)

    cfg["igv_flanking"] = (
        args.igv_flanking
        if getattr(args, "igv_flanking", None) is not None
        else cfg.get("igv_flanking", 50)
    )

    # Handle precedence and validation
    if args.igv:
        if not args.bam_mapping_file:
            msg = "For IGV integration, --bam-mapping-file must be provided."
            logger.error(msg)
            if parser is not None:
                parser.error(msg)
            raise ValueError(msg)

        if not args.igv_fasta and not args.igv_reference:
            msg = "For IGV integration, either --igv-reference or --igv-fasta must be provided."
            logger.error(msg)
            if parser is not None:
                parser.error(msg)
            raise ValueError(msg)

        if args.igv_fasta and args.igv_reference:
            logger.warning(
                "Both --igv-fasta and --igv-reference provided. "
                "Local FASTA file will take precedence."
            )

        cfg["igv_reference"] = args.igv_reference
    # MODIFIED: End of local IGV FASTA feature

    # COAST auto-field injection (COAST-02, COAST-05)
    # When COAST is one of the selected association tests, auto-inject annotation fields
    # required by the selected classification model into the field extraction list.
    _association_tests = cfg.get("association_tests", [])
    if "coast" in _association_tests:
        _coast_model = getattr(args, "coast_classification", "sift_polyphen")
        _scoring_base = Path(__file__).resolve().parent.parent / "scoring"
        if not _scoring_base.is_dir():
            _scoring_base = Path.cwd() / "scoring"
        _coast_model_dir = str(_scoring_base / "coast_classification" / _coast_model)
        if not Path(_coast_model_dir).is_dir():
            msg = (
                f"COAST classification model not found: {_coast_model_dir}. "
                f"Available models are in scoring/coast_classification/."
            )
            logger.error(msg)
            if parser is not None:
                parser.error(msg)
            raise FileNotFoundError(msg)

        _var_config_path = Path(_coast_model_dir) / "variable_assignment_config.json"
        try:
            import json as _json

            with open(_var_config_path, encoding="utf-8") as _f:
                _var_config = _json.load(_f)
        except (OSError, ValueError) as _e:
            msg = f"COAST: failed to read variable_assignment_config: {_e}"
            logger.error(msg)
            if parser is not None:
                parser.error(msg)
            raise ValueError(msg) from _e

        # Read vcf_fields list (actual VCF annotation field names for auto-injection).
        # The "variables" keys are internal COAST_* names used by the formula engine,
        # not VCF field names, so we use the explicit vcf_fields list instead.
        _vcf_fields = _var_config.get("vcf_fields", [])

        if _vcf_fields:
            if fields:
                # fields_to_extract uses space-separated format; --fields uses commas.
                # Support both by splitting on commas first, then whitespace.
                if "," in fields:
                    _field_list = [f.strip() for f in fields.split(",")]
                else:
                    _field_list = fields.split()
                _injected = []
                for _rf in _vcf_fields:
                    if _rf not in _field_list:
                        _field_list.append(_rf)
                        _injected.append(_rf)
                if _injected:
                    logger.info(
                        f"COAST auto-injected fields for '{_coast_model}' model: {_injected}"
                    )
                # Rejoin using the same delimiter as the input
                _sep = "," if "," in fields else " "
                fields = _sep.join(_field_list)
            else:
                fields = " ".join(_vcf_fields)
                logger.info(
                    f"COAST set extraction fields for '{_coast_model}' model: {_vcf_fields}"
                )

        # Store resolved absolute path in cfg for AssociationConfig propagation
        cfg["coast_classification"] = _coast_model_dir
        logger.debug(f"COAST classification model resolved to: {_coast_model_dir}")

    # Update reference/filters/fields
    cfg["reference"] = reference
    cfg["filters"] = filters
    cfg["fields_to_extract"] = fields

    # Remove sample substring (optional)
    cfg["remove_sample_substring"] = args.remove_sample_substring

    # Toggle link columns
    cfg["no_links"] = args.no_links

    # Resolve build-aware external link templates
    from .links import resolve_links

    resolve_links(cfg)

    from .pipeline_validation import resolve_threads

    raw_threads = getattr(args, "threads", None) or cfg.get("threads", "auto")
    resolved_threads = resolve_threads(raw_threads, cfg)
    cfg["threads"] = resolved_threads
    args.threads = resolved_threads

    # Output file - important to override default config value
    if args.output_file is not None:
        cfg["output_file"] = args.output_file

    # Transcript list/file, genotype filter, splitting mode, PM5, scoring, phenotype
    for _cfg_k, _arg_v in (
        ("transcript_list", args.transcript_list),
        ("transcript_file", args.transcript_file),
        ("fields_file", getattr(args, "fields_file", None)),
        ("filter_exception_file", args.filter_exception_file),
        ("genotype_filter", args.genotype_filter),
        ("gene_genotype_file", args.gene_genotype_file),
        ("snpeff_splitting_mode", args.split_snpeff_lines),
        ("clinvar_pm5_lookup", args.clinvar_pm5_lookup),
        ("scoring_config_path", args.scoring_config_path),
        ("phenotype_file", args.phenotype_file),
        ("phenotype_sample_column", args.phenotype_sample_column),
        ("phenotype_value_column", args.phenotype_value_column),
        ("json_gene_mapping", args.json_gene_mapping),
    ):
        resolved = _arg_v if _arg_v is not None else cfg.get(_cfg_k)
        cfg[_cfg_k] = resolved
        if resolved is not None:
            setattr(args, _cfg_k, resolved)

    if cfg.get("filter_exception_file") and cfg.get("fields_to_extract"):
        raw_fields = str(cfg["fields_to_extract"])
        sep = "," if "," in raw_fields and " " not in raw_fields else " "
        parts = [
            p.strip()
            for p in (raw_fields.split(",") if sep == "," else raw_fields.split())
            if p.strip()
        ]
        if "RAW_FILTER" not in parts and "INFO/RAW_FILTER" not in parts:
            parts.append("RAW_FILTER")
            cfg["fields_to_extract"] = sep.join(parts)

    # Unified extra sample fields logic:
    if args.append_extra_sample_fields is not None:
        # Flag used => enable genotype+field appending
        cfg["append_extra_sample_fields"] = True
        # Possibly empty list if used with no args
        cfg["extra_sample_fields"] = args.append_extra_sample_fields
    else:
        # Flag not used at all => preserve config file settings if present
        cfg["extra_sample_fields"] = cfg.get("extra_sample_fields") or []
        cfg["append_extra_sample_fields"] = bool(
            cfg.get("append_extra_sample_fields", False) or cfg["extra_sample_fields"]
        )

    if args.extra_sample_field_delimiter is not None:
        cfg["extra_sample_field_delimiter"] = args.extra_sample_field_delimiter
    else:
        cfg.setdefault("extra_sample_field_delimiter", ":")

    # Unified annotation configuration
    for opt, (k1, k2) in {
        "annotate_bed": ("annotate_bed_files", "annotate_bed"),
        "annotate_gene_list": ("annotate_gene_lists", "annotate_gene_list_files"),
    }.items():
        v = getattr(args, opt, None) or cfg.get(k1) or cfg.get(k2)
        if v:
            setattr(args, opt, v)
            cfg[k1] = v
            cfg[k2] = v

    ann_json = getattr(args, "annotate_json_genes", None) or cfg.get("annotate_json_genes")
    if ann_json:
        args.annotate_json_genes = ann_json
        cfg["annotate_json_genes"] = ann_json
    else:
        cfg.setdefault("annotate_json_genes", [])
    cfg["json_genes_as_columns"] = bool(
        getattr(args, "json_genes_as_columns", False) or cfg.get("json_genes_as_columns", False)
    )
    normalize_annotation_config(cfg)

    # Validate unified annotation arguments
    if args.annotate_json_genes and not args.json_gene_mapping:
        msg = "--annotate-json-genes requires --json-gene-mapping to be provided."
        logger.error(msg)
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)
    if args.json_gene_mapping and not args.annotate_json_genes:
        msg = "--json-gene-mapping requires --annotate-json-genes to be provided."
        logger.error(msg)
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    # Validate --regions-bed file existence
    if getattr(args, "regions_bed", None) and not Path(args.regions_bed).is_file():
        _fail(f"--regions-bed file not found: {args.regions_bed}")

    # Validate association analysis arguments
    if args.association_tests and not args.perform_association:
        _fail("--association-tests requires --perform-association to be set")

    if getattr(args, "variant_weight_column", None) and not args.perform_association:
        _fail("--variant-weight-column requires --perform-association to be set")

    if args.perform_association:
        variant_weights_arg = getattr(args, "variant_weights", "beta:1,25")
        if variant_weights_arg == "score_column" and not getattr(
            args, "variant_weight_column", None
        ):
            _fail("--variant-weight-column is required when --variant-weights score_column")
        if (
            variant_weights_arg.startswith("column:")
            and not variant_weights_arg[len("column:") :].strip()
        ):
            _fail("--variant-weights column:<name> requires a non-empty column name")

    # These outputs only make sense with association analysis
    for dest in ("covariate_file", "diagnostics_output", "score_artifact"):
        if getattr(args, dest, None) and not args.perform_association:
            _fail(f"--{dest.replace('_', '-')} requires --perform-association to be set")
    if getattr(args, "score_artifact", None) and getattr(args, "trait_type", "binary") != "binary":
        _fail("--score-artifact supports binary traits only")

    # PCA arg only makes sense with association analysis
    if getattr(args, "pca", None) and not args.perform_association:
        _fail("--pca requires --perform-association to be set")

    # Trait type / test compatibility check
    trait_type = getattr(args, "trait_type", "binary")
    raw_tests = cfg.get("association_tests") or getattr(args, "association_tests", None)
    if trait_type == "quantitative" and raw_tests:
        test_list = [
            t.strip() for t in (raw_tests.split(",") if isinstance(raw_tests, str) else raw_tests)
        ]
        if "logistic_burden" in test_list:
            _fail(
                "--trait-type quantitative is incompatible with logistic_burden "
                "(logistic regression is for binary traits only). "
                "Use linear_burden for quantitative traits."
            )

    # Inheritance analysis configuration
    cfg["ped_file"] = args.ped if args.ped is not None else cfg.get("ped_file")
    # Enable inheritance analysis if inheritance mode is specified OR if PED file is provided
    # This allows single-sample analysis without PED
    cfg["calculate_inheritance"] = bool(
        args.inheritance_mode is not None
        or args.ped is not None
        or cfg.get("calculate_inheritance", False)
    )
    # If inheritance is calculated but no mode specified, default to "simple"
    cfg["inheritance_mode"] = (
        args.inheritance_mode
        or cfg.get("inheritance_mode")
        or ("simple" if cfg["calculate_inheritance"] else None)
    )
    cfg[
        "use_vectorized_comp_het"
    ] = not args.no_vectorized_comp_het  # Default to True unless disabled

    # Late filtering configuration
    cfg["late_filtering"] = bool(
        args.late_filtering or cfg.get("late_filtering") or cfg.get("Late_filtering") or False
    )

    # Bcftools pre-filtering configuration
    cfg["bcftools_prefilter"] = (
        args.bcftools_prefilter
        if args.bcftools_prefilter is not None
        else cfg.get("bcftools_prefilter")
    )

    # Final filter configuration
    final_filter = args.final_filter if args.final_filter is not None else cfg.get("final_filter")
    cfg["final_filter"] = (
        _expand_tumor_normal_vars(final_filter, tumor_normal_vars) if final_filter else final_filter
    )
    cfg["fail_on_empty"] = bool(args.fail_on_empty or cfg.get("fail_on_empty", False))

    # Chunked processing configuration
    cfg["no_chunked_processing"] = args.no_chunked_processing
    cfg["force_chunked_processing"] = args.force_chunked_processing
    sort_limit_arg = getattr(args, "sort_memory_limit", None)
    if sort_limit_arg is not None:
        sort_limit = sort_limit_arg
    elif "sort_memory_limit" in cfg and cfg["sort_memory_limit"] is not None:
        sort_limit = cfg["sort_memory_limit"]
    else:
        sort_limit = "auto"
    cfg["sort_memory_limit"] = _calculate_sort_memory_limit(sort_limit, cfg.get("max_memory_gb"))
    cfg["sort_parallel"] = args.sort_parallel

    # Inheritance memory configuration
    cfg["force_inheritance_processing"] = args.force_inheritance_processing
    cfg["sample_column_creation_method"] = args.sample_column_creation_method
    cfg["memory_safety_factor"] = args.memory_safety_factor
    cfg["inheritance_memory_fraction"] = args.inheritance_memory_fraction

    # Statistics configuration
    cfg["stats_config"] = args.stats_config

    # Checkpoint configuration
    cfg["enable_checkpoint"] = args.enable_checkpoint
    cfg["resume"] = args.resume
    cfg["checkpoint_checksum"] = args.checkpoint_checksum
    cfg["resume_from"] = (
        args.resume_from
        if getattr(args, "resume_from", None) is not None
        else cfg.get("resume_from")
    )
    cfg["list_stages"] = args.list_stages
    cfg["list_checkpoints"] = args.list_checkpoints
    cfg["interactive_resume"] = args.interactive_resume
    cfg["pipeline_version"] = __version__  # Add pipeline version for compatibility checking


from .cli_fields import load_fields_from_file

__all__ = [
    "load_fields_from_file",
    "populate_cli_config",
    "resolve_presets_and_tumor_normal",
]
