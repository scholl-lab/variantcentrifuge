"""Association/statistics part of the CLI args-to-config mapping."""

from __future__ import annotations

import argparse
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def populate_association_config(
    args: argparse.Namespace, cfg: dict[str, Any], _fail: Callable[[str], None]
) -> None:
    """Map statistics and association arguments into *cfg*."""
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
    cfg["no_sample_subset"] = bool(
        getattr(args, "no_sample_subset", False) or cfg.get("no_sample_subset", False)
    )
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
    cfg["mask_label"] = getattr(args, "mask_label", None)
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
    _gsf = getattr(args, "gene_set_file", None)
    if _gsf and not Path(_gsf).is_file():
        _fail(f"--gene-set-file file not found: {_gsf}")
    cfg["gene_set_file"] = _gsf
    cfg["gene_prior_weight_column"] = getattr(args, "gene_prior_weight_column", "weight")
    cfg["association_min_cases"] = getattr(args, "min_cases", 200)
    cfg["association_max_case_control_ratio"] = getattr(args, "max_case_control_ratio", 20.0)
    cfg["association_min_case_carriers"] = getattr(args, "min_case_carriers", 10)
    # Unset on the CLI (None): the config file / JSON "association" section or default applies.
    for key in ("missing_site_threshold", "missing_sample_threshold", "missing_as_ref"):
        if getattr(args, key, None) is not None:
            cfg[key] = getattr(args, key)


def inject_coast_fields(
    args: argparse.Namespace,
    cfg: dict[str, Any],
    parser: argparse.ArgumentParser | None,
    fields: str | None,
) -> str | None:
    """Resolve the COAST model and inject its required fields; returns the updated fields."""
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
    return fields


def validate_association_args(
    args: argparse.Namespace, cfg: dict[str, Any], _fail: Callable[[str], None]
) -> None:
    """Validate association-related argument combinations."""
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
    for dest in (
        "covariate_file",
        "diagnostics_output",
        "score_artifact",
        "mask_label",
        "gene_set_file",
        "missing_site_threshold",
        "missing_sample_threshold",
        "missing_as_ref",
    ):
        if getattr(args, dest, None) is not None and not args.perform_association:
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
        binary_only = sorted({"fisher", "logistic_burden"} & set(test_list))
        if binary_only:
            _fail(
                f"--trait-type quantitative is incompatible with {', '.join(binary_only)} "
                "(binary-trait tests; fisher is the default when --association-tests is not "
                "given). Use linear_burden and/or skat for quantitative traits."
            )
