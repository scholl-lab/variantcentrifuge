"""Configuration post-processing and validation helpers for VariantCentrifuge CLI."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

from .config import normalize_annotation_config
from .vcf_header_parser import resolve_tumor_normal_indices
from .version import __version__

logger = logging.getLogger("variantcentrifuge")


def _calculate_sort_memory_limit(sort_limit: str, max_memory_gb: float | None = None) -> str:
    """Calculate intelligent sort memory limit based on available memory."""
    if sort_limit != "auto":
        return sort_limit

    try:
        import psutil

        # Use provided max memory or detect available
        if max_memory_gb:
            available_gb = max_memory_gb
        else:
            memory_info = psutil.virtual_memory()
            available_gb = memory_info.available / (1024**3)

        # Use 10% of available memory for sorting, with reasonable bounds
        sort_memory_gb = max(1, min(available_gb * 0.10, 50))  # Between 1GB and 50GB

        if sort_memory_gb >= 1:
            return f"{int(sort_memory_gb)}G"
        else:
            return f"{int(sort_memory_gb * 1024)}M"

    except Exception as e:
        logger.debug(f"Could not calculate sort memory limit: {e}. Using default 4G")
        return "4G"  # Better default than 2G for modern systems


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
    preset_dict = cfg.get("presets", {})
    raw_presets = args.preset or cfg.get("preset") or []
    if isinstance(raw_presets, str):
        raw_presets = [raw_presets]
    raw_final_filter = (
        args.final_filter if args.final_filter is not None else cfg.get("final_filter")
    )
    expr_texts = [preset_dict.get(p.strip(), "") for p in raw_presets]
    if filters:
        expr_texts.append(filters)
    if raw_final_filter:
        expr_texts.append(raw_final_filter)
    needs_tumor_idx = bool(
        args.tumor_sample_index is not None
        or cfg.get("tumor_sample_index") is not None
        or any("{tumor_idx}" in e for e in expr_texts)
    )
    needs_normal_idx = bool(
        args.normal_sample_index is not None
        or cfg.get("normal_sample_index") is not None
        or any("{normal_idx}" in e for e in expr_texts)
    )
    needs_tumor_normal = needs_tumor_idx or needs_normal_idx
    t_idx, n_idx = resolve_tumor_normal_indices(
        args.vcf_file if needs_tumor_normal else None,
        cli_tumor_idx=args.tumor_sample_index,
        cli_normal_idx=args.normal_sample_index,
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
        chosen_presets = []
        for p in raw_presets:
            p_str = p.strip()
            if p_str in preset_dict:
                expr = _expand_tumor_normal_vars(preset_dict[p_str], tumor_normal_vars)
                chosen_presets.append(f"({expr})")
            else:
                logger.error(f"Preset '{p_str}' not found in config.")
                sys.exit(1)
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
    parser: argparse.ArgumentParser,
    reference: str | None,
    filters: str | None,
    fields: str | None,
    tumor_normal_vars: dict[str, str],
) -> None:
    """Populate and validate pipeline configuration from CLI arguments."""
    cfg["perform_gene_burden"] = args.perform_gene_burden
    cfg["no_stats"] = args.no_stats
    cfg["gene_burden_mode"] = args.gene_burden_mode
    cfg["correction_method"] = args.correction_method

    # Association analysis configuration
    cfg["perform_association"] = args.perform_association
    if args.association_tests:
        cfg["association_tests"] = [t.strip() for t in args.association_tests.split(",")]
    else:
        cfg["association_tests"] = ["fisher"] if args.perform_association else []
    cfg["skat_backend"] = getattr(args, "skat_backend", "python")
    cfg["coast_backend"] = getattr(args, "coast_backend", "python")
    # coast_classification: resolved to absolute path later during auto-field injection
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
            import sys as _sys

            print(
                f"error: --variant-weight-params: invalid JSON: {_e}",
                file=_sys.stderr,
            )
            _sys.exit(2)
    else:
        cfg["variant_weight_params"] = None
    cfg["diagnostics_output"] = getattr(args, "diagnostics_output", None)
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
                import sys as _sys

                print(
                    "error: --coast-weights must have exactly 3 values (BMV, DMV, PTV)",
                    file=_sys.stderr,
                )
                _sys.exit(2)
            cfg["coast_weights"] = _coast_weights_parsed
        except ValueError as _e:
            import sys as _sys

            print(
                f"error: --coast-weights values must be numeric: {_e}",
                file=_sys.stderr,
            )
            _sys.exit(2)
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
        parser.error(f"--gene-prior-weights file not found: {_gpw}")
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
    cfg["igv_enabled"] = args.igv
    cfg["bam_mapping_file"] = args.bam_mapping_file

    # MODIFIED: Start of local IGV FASTA feature
    # Store local FASTA-related parameters
    cfg["igv_fasta"] = args.igv_fasta
    cfg["igv_ideogram"] = args.igv_ideogram

    # MODIFIED: IGV filename shortening parameters
    cfg["igv_max_allele_len_filename"] = args.igv_max_allele_len_filename
    cfg["igv_hash_len_filename"] = args.igv_hash_len_filename
    cfg["igv_max_variant_part_filename"] = args.igv_max_variant_part_filename

    # MODIFIED: Start of IGV flanking feature
    # Use CLI argument if provided, otherwise use config value, with a default of 50
    cfg["igv_flanking"] = args.igv_flanking or cfg.get("igv_flanking", 50)
    # MODIFIED: End of IGV flanking feature

    # Handle precedence and validation
    if args.igv:
        if not args.bam_mapping_file:
            logger.error("For IGV integration, --bam-mapping-file must be provided.")
            sys.exit(1)

        if not args.igv_fasta and not args.igv_reference:
            logger.error(
                "For IGV integration, either --igv-reference or --igv-fasta must be provided."
            )
            sys.exit(1)

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
        import json as _json
        import os as _os

        _coast_model = getattr(args, "coast_classification", "sift_polyphen")
        # Resolve scoring base dir: package is at variantcentrifuge/, scoring/ is one level up
        _pkg_dir = _os.path.dirname(_os.path.abspath(__file__))
        _scoring_base = _os.path.join(_os.path.dirname(_pkg_dir), "scoring")
        if not _os.path.isdir(_scoring_base):
            _scoring_base = _os.path.join(_os.getcwd(), "scoring")
        _coast_model_dir = _os.path.join(_scoring_base, "coast_classification", _coast_model)
        if not _os.path.isdir(_coast_model_dir):
            logger.error(
                f"COAST classification model not found: {_coast_model_dir}. "
                f"Available models are in scoring/coast_classification/."
            )
            sys.exit(1)

        # Read required VCF fields from variable_assignment_config
        _var_config_path = _os.path.join(_coast_model_dir, "variable_assignment_config.json")
        try:
            with open(_var_config_path, encoding="utf-8") as _f:
                _var_config = _json.load(_f)
        except (OSError, _json.JSONDecodeError) as _e:
            logger.error(f"COAST: failed to read variable_assignment_config: {_e}")
            sys.exit(1)

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

    # Threads — resolve "auto" to detected CPU count
    if str(args.threads).lower() == "auto":
        from .memory.resource_manager import ResourceManager

        rm = ResourceManager(config=cfg)
        resolved_threads = rm.cpu_cores
        logger.info(f"Auto-detected {resolved_threads} CPU cores for --threads")
        cfg["threads"] = resolved_threads
        args.threads = resolved_threads
    else:
        cfg["threads"] = int(args.threads)
        args.threads = int(args.threads)

    # Output file - important to override default config value
    if args.output_file is not None:
        cfg["output_file"] = args.output_file

    # Transcript list/file, genotype filter, splitting mode, PM5, scoring, phenotype
    for _cfg_k, _arg_v in (
        ("transcript_list", args.transcript_list),
        ("transcript_file", args.transcript_file),
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
        cfg[_cfg_k] = _arg_v if _arg_v is not None else cfg.get(_cfg_k)

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
    if args.annotate_bed:
        cfg["annotate_bed_files"] = args.annotate_bed
    if args.annotate_json_genes:
        cfg["annotate_json_genes"] = args.annotate_json_genes
    else:
        cfg.setdefault("annotate_json_genes", [])
    cfg["json_genes_as_columns"] = bool(
        args.json_genes_as_columns or cfg.get("json_genes_as_columns", False)
    )
    if args.annotate_gene_list:
        cfg["annotate_gene_list_files"] = args.annotate_gene_list
        cfg["annotate_gene_lists"] = args.annotate_gene_list
    normalize_annotation_config(cfg)

    # Validate unified annotation arguments
    if args.annotate_json_genes and not args.json_gene_mapping:
        logger.error("--annotate-json-genes requires --json-gene-mapping to be provided.")
        sys.exit(1)
    if args.json_gene_mapping and not args.annotate_json_genes:
        logger.error("--json-gene-mapping requires --annotate-json-genes to be provided.")
        sys.exit(1)

    # Validate --regions-bed file existence
    if getattr(args, "regions_bed", None) and not Path(args.regions_bed).is_file():
        parser.error(f"--regions-bed file not found: {args.regions_bed}")

    # Validate association analysis arguments
    if args.association_tests and not args.perform_association:
        parser.error("--association-tests requires --perform-association to be set")

    if getattr(args, "variant_weight_column", None) and not args.perform_association:
        parser.error("--variant-weight-column requires --perform-association to be set")

    if args.perform_association:
        variant_weights_arg = getattr(args, "variant_weights", "beta:1,25")
        if variant_weights_arg == "score_column" and not getattr(
            args, "variant_weight_column", None
        ):
            parser.error("--variant-weight-column is required when --variant-weights score_column")
        if (
            variant_weights_arg.startswith("column:")
            and not variant_weights_arg[len("column:") :].strip()
        ):
            parser.error("--variant-weights column:<name> requires a non-empty column name")

    # Covariate file only makes sense with association analysis
    if getattr(args, "covariate_file", None) and not args.perform_association:
        parser.error("--covariate-file requires --perform-association to be set")

    # Diagnostics output only makes sense with association analysis
    if getattr(args, "diagnostics_output", None) and not args.perform_association:
        parser.error("--diagnostics-output requires --perform-association to be set")

    # PCA arg only makes sense with association analysis
    if getattr(args, "pca", None) and not args.perform_association:
        parser.error("--pca requires --perform-association to be set")

    # Trait type / test compatibility check
    trait_type = getattr(args, "trait_type", "binary")
    if trait_type == "quantitative" and args.association_tests:
        test_list = [t.strip() for t in args.association_tests.split(",")]
        if "logistic_burden" in test_list:
            parser.error(
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
    cfg["late_filtering"] = bool(args.late_filtering or cfg.get("late_filtering", False))

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
    cfg["sort_memory_limit"] = _calculate_sort_memory_limit(
        args.sort_memory_limit, cfg.get("max_memory_gb")
    )
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
    cfg["resume_from"] = args.resume_from
    cfg["list_stages"] = args.list_stages
    cfg["list_checkpoints"] = args.list_checkpoints
    cfg["interactive_resume"] = args.interactive_resume
    cfg["pipeline_version"] = __version__  # Add pipeline version for compatibility checking
