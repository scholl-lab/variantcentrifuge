"""Reusable analysis profile resolution, precedence merging, and persistence."""

from __future__ import annotations

import argparse
import copy
import json
import logging
import os
import tempfile
from pathlib import Path
from typing import Any

logger = logging.getLogger("variantcentrifuge")

EXCLUDED_PROFILE_KEYS: frozenset[str] = frozenset(
    {
        # Sample-specific and run-specific I/O paths
        "vcf_file",
        "from_tsv",
        "output_file",
        "output_dir",
        "log_file",
        "ped",
        "ped_file",
        "phenotype_file",
        "bam_mapping_file",
        "sample_sheet",
        "case_samples",
        "control_samples",
        "case_samples_file",
        "control_samples_file",
        "case_phenotypes_file",
        "control_phenotypes_file",
        "sample_file",
        "sample_list",
        "pseudonymize_table",
        "pseudonymize_metadata",
        "pseudonymize_ped",
        "stats_output_file",
        "diagnostics_output",
        "covariate_file",
        "pca_file",
        "resume",
        "resume_from",
        "interactive_resume",
        "list_stages",
        "list_checkpoints",
        "show_checkpoint_status",
        "save_profile",
        "list_profiles",
        "profile",
        "config",
        # Structural config sections that live in base config.json
        "profiles",
        "field_profiles",
        "presets",
        "links",
        "links_by_build",
        "igv_links",
        "pipeline_version",
    }
)


def get_explicit_cli_dests(
    parser: argparse.ArgumentParser | None,
    argv: list[str] | None,
) -> set[str]:
    """Return the set of parser action destinations explicitly passed in argv."""
    if parser is None or not argv:
        return set()

    opt_to_dest: dict[str, str] = {}
    for action in getattr(parser, "_actions", []):
        for opt in getattr(action, "option_strings", []):
            opt_to_dest[opt] = action.dest

    explicit: set[str] = set()
    for token in argv:
        if not token.startswith("-"):
            continue
        opt_name = token.split("=", 1)[0]
        if opt_name in opt_to_dest:
            explicit.add(opt_to_dest[opt_name])
    return explicit


def list_analysis_profiles(config: dict[str, Any]) -> list[dict[str, str]]:
    """Return summary of available built-in or configured analysis profiles."""
    profiles = config.get("profiles", {})
    if not isinstance(profiles, dict) or not profiles:
        from .config import load_config

        default_cfg = load_config()
        profiles = default_cfg.get("profiles", {})

    result: list[dict[str, str]] = []
    if isinstance(profiles, dict):
        for name, prof in profiles.items():
            desc = prof.get("description", "") if isinstance(prof, dict) else ""
            result.append({"name": str(name), "description": str(desc)})
    return result


def resolve_analysis_profile(profile_spec: str, config: dict[str, Any]) -> dict[str, Any]:
    """Resolve a named analysis profile from config or load from a JSON file path."""
    if not profile_spec or not str(profile_spec).strip():
        raise ValueError("Analysis profile name or path must not be empty.")

    spec = str(profile_spec).strip()
    profiles = config.get("profiles", {})
    if not isinstance(profiles, dict) or spec not in profiles:
        from .config import load_config

        try:
            default_profiles = load_config().get("profiles", {})
            if isinstance(default_profiles, dict):
                merged = dict(default_profiles)
                if isinstance(profiles, dict):
                    merged.update(profiles)
                profiles = merged
        except Exception:
            if not isinstance(profiles, dict):
                profiles = {}

    if isinstance(profiles, dict) and spec in profiles:
        raw_profile = profiles[spec]
        if not isinstance(raw_profile, dict):
            raise ValueError(f"Analysis profile '{spec}' in configuration must be a JSON object.")
        resolved = copy.deepcopy(raw_profile)
        resolved.pop("description", None)
        return resolved

    cand_path = Path(spec)
    if cand_path.exists():
        if not cand_path.is_file():
            raise ValueError(f"Analysis profile path is not a file: {spec}")
        try:
            loaded = json.loads(cand_path.read_text(encoding="utf-8"))
        except Exception as exc:
            raise ValueError(f"Failed to parse analysis profile JSON '{spec}': {exc}") from exc
        if not isinstance(loaded, dict):
            raise ValueError(f"Analysis profile '{spec}' must be a JSON object.")
        if (
            "profile" in loaded
            and isinstance(loaded["profile"], dict)
            and len(loaded) <= 2
            and all(k in ("profile", "description") for k in loaded)
        ):
            loaded = loaded["profile"]
        resolved = copy.deepcopy(loaded)
        resolved.pop("description", None)
        for excluded in ("profiles", "field_profiles", "presets"):
            resolved.pop(excluded, None)
        return resolved

    available = ", ".join(sorted(profiles.keys())) if isinstance(profiles, dict) else "none"
    raise ValueError(
        f"Analysis profile '{spec}' not found. Available built-in profiles: {available}"
    )


def apply_analysis_profile(
    args: argparse.Namespace,
    cfg: dict[str, Any],
    parser: argparse.ArgumentParser | None = None,
    argv: list[str] | None = None,
) -> dict[str, Any]:
    """Apply analysis profile with precedence: explicit CLI > --profile > base config > defaults."""
    explicit_dests = get_explicit_cli_dests(parser, argv)
    explicit_values = {
        dest: copy.deepcopy(getattr(args, dest)) for dest in explicit_dests if hasattr(args, dest)
    }

    # 1. Propagate base config values into non-explicit CLI args so base config beats defaults
    for key, val in cfg.items():
        if key in EXCLUDED_PROFILE_KEYS and key not in (
            "vcf_file",
            "from_tsv",
            "output_file",
            "output_dir",
            "ped",
            "ped_file",
            "phenotype_file",
            "bam_mapping_file",
            "sample_sheet",
            "case_samples",
            "control_samples",
            "case_samples_file",
            "control_samples_file",
            "case_phenotypes_file",
            "control_phenotypes_file",
        ):
            continue
        if hasattr(args, key) and key not in explicit_dests and val is not None:
            setattr(args, key, copy.deepcopy(val))

    # 2. Resolve and apply --profile on top of base config, respecting explicit CLI flags
    profile_spec = getattr(args, "profile", None) or cfg.get("profile")
    if profile_spec:
        profile_settings = resolve_analysis_profile(str(profile_spec), cfg)
        logger.info("Applying analysis profile: %s", profile_spec)
        for key, val in profile_settings.items():
            if key in explicit_dests:
                # Explicit CLI flag overrides profile
                cfg[key] = copy.deepcopy(explicit_values[key])
            else:
                cfg[key] = copy.deepcopy(val)
                if hasattr(args, key):
                    setattr(args, key, copy.deepcopy(val))
                if key in ("ped", "ped_file") and "ped" not in explicit_dests:
                    args.ped = copy.deepcopy(val)
                    args.ped_file = copy.deepcopy(val)

    # 3. Ensure all explicit CLI values are reflected in args and cfg
    for dest, val in explicit_values.items():
        setattr(args, dest, copy.deepcopy(val))
        if dest not in ("config", "profile", "save_profile", "list_profiles") and val is not None:
            cfg[dest] = copy.deepcopy(val)

    return cfg


def save_analysis_profile(
    cfg: dict[str, Any],
    save_path: str | Path,
    base_fields_to_extract: str | None = None,
) -> Path:
    """Save reusable analysis configuration to a JSON profile file atomically."""
    target = Path(save_path)
    target.parent.mkdir(parents=True, exist_ok=True)

    export_data: dict[str, Any] = {}
    for key, val in cfg.items():
        if key.startswith("_") or key in EXCLUDED_PROFILE_KEYS or val is None:
            continue
        export_data[key] = copy.deepcopy(val)

    if base_fields_to_extract is not None:
        export_data["fields_to_extract"] = base_fields_to_extract

    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{target.name}.",
        suffix=".tmp",
        dir=str(target.parent),
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(export_data, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(tmp_name, target)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise

    logger.info("Saved reusable analysis profile to %s", target)
    return target
