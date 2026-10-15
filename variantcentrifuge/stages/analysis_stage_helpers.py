"""
Association configuration validation and resolution helpers for analysis stages.
"""

import logging
import os
from typing import TYPE_CHECKING, Any

import pandas as pd

from ..association.base import AssociationConfig

if TYPE_CHECKING:
    from ..pipeline_core import PipelineContext

logger = logging.getLogger(__name__)

# Valid keys for the "association" section in config.json
VALID_ASSOCIATION_KEYS: frozenset[str] = frozenset(
    {
        "correction_method",
        "gene_burden_mode",
        "trait_type",
        "variant_weights",
        "variant_weight_params",
        "variant_weight_column",
        "skat_backend",
        "skat_method",
        "coast_backend",
        "covariate_file",
        "covariate_columns",
        "categorical_covariates",
        "pca",
        "pca_file",
        "pca_tool",
        "pca_components",
        "coast_weights",
        "association_tests",
        "min_cases",
        "max_case_control_ratio",
        "min_case_carriers",
        "diagnostics_output",
        "confidence_interval_method",
        "confidence_interval_alpha",
        "continuity_correction",
        "missing_site_threshold",
        "missing_sample_threshold",
        "firth_max_iter",
        "association_workers",
        "gene_prior_weights",
        "gene_prior_weight_column",
    }
)


def _validate_association_config_dict(d: dict) -> None:
    """Validate association config dict. Raises ValueError listing ALL invalid keys/values.

    Performs two levels of validation:
    1. Unknown keys — any key not in VALID_ASSOCIATION_KEYS is rejected.
    2. Type and enum validation — known keys must have the correct Python type.

    All errors are collected before raising so the user sees everything wrong
    in a single error message (fail-fast but complete).

    Parameters
    ----------
    d : dict
        The "association" sub-dict from config.json.

    Raises
    ------
    ValueError
        If any validation errors are found; message lists all errors.
    """
    errors: list[str] = []

    unknown = sorted(set(d) - VALID_ASSOCIATION_KEYS)
    if unknown:
        errors.append(f"Unknown keys: {unknown}")

    # Type validation sets (intersection with actual keys avoids KeyError)
    str_keys = {
        "correction_method",
        "gene_burden_mode",
        "trait_type",
        "variant_weights",
        "variant_weight_column",
        "skat_backend",
        "skat_method",
        "coast_backend",
        "covariate_file",
        "pca_file",
        "pca_tool",
        "confidence_interval_method",
        "diagnostics_output",
        "gene_prior_weights",
        "gene_prior_weight_column",
    }
    int_keys = {
        "pca_components",
        "min_cases",
        "min_case_carriers",
        "firth_max_iter",
        "association_workers",
    }
    float_keys = {
        "max_case_control_ratio",
        "confidence_interval_alpha",
        "continuity_correction",
        "missing_site_threshold",
        "missing_sample_threshold",
    }
    list_str_keys = {"covariate_columns", "categorical_covariates", "association_tests"}
    list_float_keys = {"coast_weights"}
    dict_keys = {"variant_weight_params"}

    for key in str_keys & set(d):
        if d[key] is not None and not isinstance(d[key], str):
            errors.append(f"'{key}' must be a string, got {type(d[key]).__name__}")

    for key in int_keys & set(d):
        if not isinstance(d[key], int):
            errors.append(f"'{key}' must be an integer, got {type(d[key]).__name__}")

    for key in float_keys & set(d):
        if not isinstance(d[key], (int, float)):
            errors.append(f"'{key}' must be a number, got {type(d[key]).__name__}")

    for key in list_str_keys & set(d):
        if d[key] is not None and not isinstance(d[key], list):
            errors.append(f"'{key}' must be a list, got {type(d[key]).__name__}")

    for key in list_float_keys & set(d):
        if d[key] is not None and not isinstance(d[key], list):
            errors.append(f"'{key}' must be a list, got {type(d[key]).__name__}")

    for key in dict_keys & set(d):
        if d[key] is not None and not isinstance(d[key], dict):
            errors.append(f"'{key}' must be an object, got {type(d[key]).__name__}")

    # Enum-like value validation
    if "correction_method" in d and d["correction_method"] not in ("fdr", "bonferroni"):
        errors.append(
            f"'correction_method' must be 'fdr' or 'bonferroni', got '{d['correction_method']}'"
        )
    if "trait_type" in d and d["trait_type"] not in ("binary", "quantitative"):
        errors.append(f"'trait_type' must be 'binary' or 'quantitative', got '{d['trait_type']}'")
    if "skat_backend" in d and d["skat_backend"] not in ("auto", "r", "python"):
        errors.append(f"'skat_backend' must be 'auto', 'r', or 'python', got '{d['skat_backend']}'")
    if "coast_backend" in d and d["coast_backend"] not in ("auto", "python"):
        errors.append(f"'coast_backend' must be 'auto' or 'python', got '{d['coast_backend']}'")

    if errors:
        raise ValueError(
            f"Invalid 'association' config section ({len(errors)} error(s)):\n"
            + "\n".join(f"  - {e}" for e in errors)
        )


def _resolve_association_workers(cfg: dict, _get) -> int:
    """Resolve association_workers from CLI/JSON/auto.

    0 (default) means auto-allocate from --threads. -1 means all CPU cores.
    Explicit positive values are used directly.
    """
    raw = _get("association_workers", default=0, nullable=False)
    if raw == 0:
        # Auto: derive from --threads
        threads = int(cfg.get("threads", 1))
        logger.info(f"Association workers auto-allocated from --threads: {threads}")
        return threads
    if raw == -1:
        cpus = os.cpu_count() or 1
        logger.info(f"Association workers auto-detected: {cpus} CPU cores")
        return cpus
    return int(raw)


def _build_assoc_config_from_context(context: "PipelineContext") -> AssociationConfig:
    """Build AssociationConfig from PipelineContext with JSON config + CLI override precedence.

    Precedence (highest to lowest):
    1. CLI-set values — keys set directly in context.config by cli.py
    2. JSON association section — context.config["association"] from config.json
    3. AssociationConfig field defaults

    The "association" section in config.json is validated at call time with
    _validate_association_config_dict(). If the section is invalid, a ValueError
    is raised immediately (fail-fast before any processing).

    Parameters
    ----------
    context : PipelineContext
        Pipeline context carrying context.config (the merged config dict from
        config.json + CLI arguments, as built by setup_stages.ConfigurationLoadingStage
        and cli.py).

    Returns
    -------
    AssociationConfig
        Fully populated configuration with correct precedence applied.
    """
    cfg = context.config

    # Load and validate the JSON "association" sub-section (may be absent)
    json_assoc: dict = cfg.get("association", {}) or {}
    if json_assoc:
        _validate_association_config_dict(json_assoc)

    # Helper: resolve field value with CLI > JSON > default precedence.
    # CLI keys are stored directly in context.config by cli.py. If the CLI key
    # is present (not None for nullable fields, not equal to default for typed
    # fields), the CLI value wins. The JSON association section fills gaps.
    # Finally, AssociationConfig field defaults apply.
    #
    # CLI-set detection strategy: if the key exists in cfg AND the value is
    # not None (for optional fields) or is non-None, it was set by the CLI.
    # For fields with non-None defaults (e.g. variant_weights="beta:1,25"),
    # we cannot distinguish CLI-set from JSON-filled, so CLI always wins when
    # present in cfg (the CLI parser always writes the key, even for defaults).
    # This is correct: if a user explicitly runs the CLI they get CLI semantics.

    def _get(
        cli_key: str,
        json_key: str | None = None,
        default: Any = None,
        nullable: bool = True,
    ) -> Any:
        """Resolve a config field with CLI > JSON > default precedence.

        Parameters
        ----------
        cli_key : str
            Key in context.config set by cli.py.
        json_key : str | None
            Key in json_assoc (association section). If None, uses cli_key.
        default : Any
            Default value if neither CLI nor JSON provides a value.
        nullable : bool
            If True, None in context.config means "not set by CLI" and JSON
            can override. If False (non-nullable typed field), the CLI value
            always wins when the key is present.
        """
        jk = json_key if json_key is not None else cli_key
        cli_val = cfg.get(cli_key)
        json_val = json_assoc.get(jk)

        if nullable:
            # For nullable fields: None means "not set"; prefer first non-None value
            if cli_val is not None:
                return cli_val
            if json_val is not None:
                return json_val
            return default
        else:
            # For non-nullable fields: CLI key being present means CLI wins
            if cli_key in cfg:
                return cli_val
            if jk in json_assoc:
                return json_val
            return default

    return AssociationConfig(
        correction_method=_get("correction_method", default="fdr", nullable=False),
        gene_burden_mode=_get("gene_burden_mode", default="samples", nullable=False),
        confidence_interval_method=_get(
            "confidence_interval_method", default="normal_approx", nullable=False
        ),
        confidence_interval_alpha=_get("confidence_interval_alpha", default=0.05, nullable=False),
        continuity_correction=_get("continuity_correction", default=0.5, nullable=False),
        covariate_file=_get("covariate_file", default=None, nullable=True),
        covariate_columns=_get("covariate_columns", default=None, nullable=True),
        categorical_covariates=_get("categorical_covariates", default=None, nullable=True),
        trait_type=_get("trait_type", default="binary", nullable=False),
        variant_weights=_get("variant_weights", default="beta:1,25", nullable=False),
        variant_weight_params=_get("variant_weight_params", default=None, nullable=True),
        variant_weight_column=_get("variant_weight_column", default=None, nullable=True),
        missing_site_threshold=_get("missing_site_threshold", default=0.10, nullable=False),
        missing_sample_threshold=_get("missing_sample_threshold", default=0.80, nullable=False),
        firth_max_iter=_get("firth_max_iter", default=25, nullable=False),
        skat_backend=_get("skat_backend", default="python", nullable=False),
        skat_method=_get("skat_method", default="SKAT", nullable=False),
        coast_backend=_get("coast_backend", default="python", nullable=False),
        min_cases=_get("association_min_cases", json_key="min_cases", default=200, nullable=False),
        max_case_control_ratio=_get(
            "association_max_case_control_ratio",
            json_key="max_case_control_ratio",
            default=20.0,
            nullable=False,
        ),
        min_case_carriers=_get(
            "association_min_case_carriers",
            json_key="min_case_carriers",
            default=10,
            nullable=False,
        ),
        diagnostics_output=_get("diagnostics_output", default=None, nullable=True),
        pca_file=_get("pca_file", default=None, nullable=True),
        pca_tool=_get("pca_tool", default=None, nullable=True),
        pca_components=_get("pca_components", default=10, nullable=False),
        coast_weights=_get("coast_weights", default=None, nullable=True),
        association_workers=_resolve_association_workers(cfg, _get),
        coast_classification=_get("coast_classification", default=None, nullable=True),
        gene_prior_weights=_get("gene_prior_weights", default=None, nullable=True),
        gene_prior_weight_column=_get("gene_prior_weight_column", default="weight", nullable=False),
    )


def _count_significant_association_results(
    results_df: pd.DataFrame,
    threshold: float = 0.05,
) -> int:
    """Count significant association rows using the run's primary corrected p-value."""
    if "primary_qvalue" in results_df.columns:
        return int((results_df["primary_qvalue"] < threshold).sum())
    corr_cols = [c for c in results_df.columns if c.endswith("_qvalue")]
    return int((results_df[corr_cols].min(axis=1) < threshold).sum()) if corr_cols else 0
