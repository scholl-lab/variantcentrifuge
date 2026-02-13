"""Field profile resolution for version-specific annotation field names.

Field profiles handle differences between annotation database versions
(e.g., dbNSFP v4.x vs v5.x) by defining parameterized filter fragments
that get expanded at runtime based on the active profile.

Presets use ``{{fragment_name:param}}`` syntax. At resolution time, the
fragment is looked up in the active profile and ``{0}`` placeholders in
the fragment template are replaced with the parameter value.

Presets without ``{{...}}`` pass through unchanged (backward compatible).
"""

import logging
import re
from typing import Any

logger = logging.getLogger("variantcentrifuge")

# Matches {{fragment_name:param}} patterns in preset expressions.
# fragment_name: word chars (letters, digits, underscore)
# param: everything up to the closing }}
FRAGMENT_PATTERN = re.compile(r"\{\{(\w+):([^}]+)\}\}")


def resolve_profile(config: dict[str, Any]) -> dict[str, Any]:
    """Resolve field profile templates in config.

    Expands ``{{fragment:param}}`` patterns in presets and merges
    profile-specific fields into ``fields_to_extract`` and
    ``html_report_default_hidden_columns``.

    The active profile is determined by (in order of precedence):
    1. ``config["field_profile"]`` (set by ``--field-profile`` CLI arg)
    2. ``config["default_field_profile"]``
    3. Falls back to ``"dbnsfp4"``

    Parameters
    ----------
    config : dict
        The full configuration dictionary (modified in place).

    Returns
    -------
    dict
        The same config dict, with templates expanded and fields merged.

    Raises
    ------
    ValueError
        If the requested profile is not found or a preset references
        an unknown fragment name.
    """
    profiles = config.get("field_profiles", {})
    if not profiles:
        # No profiles defined — nothing to resolve (backward compatible)
        return config

    profile_name = config.get("field_profile") or config.get("default_field_profile", "dbnsfp4")

    if profile_name not in profiles:
        available = ", ".join(sorted(profiles.keys())) or "none defined"
        raise ValueError(f"Field profile '{profile_name}' not found. Available: {available}")

    profile = profiles[profile_name]
    fragments = profile.get("fragments", {})

    logger.debug("Resolving field profile '%s'", profile_name)

    # Expand templates in presets
    presets = config.get("presets", {})
    for name, expr in presets.items():
        presets[name] = _expand_fragments(expr, fragments, name)

    # Merge profile fields into fields_to_extract
    base_fields = config.get("fields_to_extract", "")
    profile_fields = profile.get("fields_to_extract", "")
    if profile_fields:
        config["fields_to_extract"] = f"{base_fields} {profile_fields}".strip()

    # Merge profile hidden columns
    base_hidden: list[str] = config.get("html_report_default_hidden_columns", [])
    profile_hidden: list[str] = profile.get("hidden_columns", [])
    for col in profile_hidden:
        if col not in base_hidden:
            base_hidden.append(col)

    logger.debug(
        "Profile '%s' resolved: %d fragments expanded, %d fields added, %d hidden columns added",
        profile_name,
        len(fragments),
        len(profile_fields.split()) if profile_fields else 0,
        len(profile_hidden),
    )

    return config


def _expand_fragments(expr: str, fragments: dict[str, str], preset_name: str) -> str:
    """Replace all ``{{fragment_name:param}}`` patterns in an expression.

    Parameters
    ----------
    expr : str
        The preset filter expression, possibly containing templates.
    fragments : dict
        Mapping of fragment name to template string with ``{0}`` placeholders.
    preset_name : str
        Name of the preset (used in error messages).

    Returns
    -------
    str
        The expression with all templates expanded.

    Raises
    ------
    ValueError
        If a referenced fragment name is not defined in the active profile.
    """

    def _replacer(match: re.Match[str]) -> str:
        frag_name = match.group(1)
        param = match.group(2)
        if frag_name not in fragments:
            raise ValueError(
                f"Preset '{preset_name}' references unknown fragment "
                f"'{frag_name}'. Available fragments: "
                f"{', '.join(sorted(fragments.keys())) or 'none'}"
            )
        return fragments[frag_name].replace("{0}", param)

    return FRAGMENT_PATTERN.sub(_replacer, expr)


def list_profiles(config: dict[str, Any]) -> list[dict[str, str]]:
    """Return summary of available field profiles.

    Parameters
    ----------
    config : dict
        The full configuration dictionary.

    Returns
    -------
    list of dict
        Each dict has ``"name"`` and ``"description"`` keys.
    """
    profiles = config.get("field_profiles", {})
    default = config.get("default_field_profile", "dbnsfp4")
    result = []
    for name, profile in profiles.items():
        desc = profile.get("description", "")
        if name == default:
            desc = f"{desc} [default]" if desc else "[default]"
        result.append({"name": name, "description": desc})
    return result
