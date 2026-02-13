"""Tests for the field_profile module.

Tests cover fragment expansion, profile resolution, error handling,
backward compatibility, and the list_profiles helper.
"""

import pytest

from variantcentrifuge.field_profile import (
    FRAGMENT_PATTERN,
    _expand_fragments,
    list_profiles,
    resolve_profile,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def dbnsfp4_config():
    """Minimal config with dbnsfp4 and dbnsfp5 profiles."""
    return {
        "default_field_profile": "dbnsfp4",
        "field_profiles": {
            "dbnsfp4": {
                "description": "dbNSFP v4.x (separate exome+genome)",
                "fields_to_extract": "dbNSFP_gnomAD_exomes_AF dbNSFP_gnomAD_genomes_AF",
                "hidden_columns": ["dbNSFP_gnomAD_exomes_AC"],
                "fragments": {
                    "gnomad_af_below": (
                        "((dbNSFP_gnomAD_exomes_AF[0] < {0}) "
                        "| (na dbNSFP_gnomAD_exomes_AC[0])) "
                        "& ((dbNSFP_gnomAD_genomes_AF[0] < {0}) "
                        "| (na dbNSFP_gnomAD_genomes_AC[0]))"
                    ),
                    "gnomad_af_below_or": (
                        "(dbNSFP_gnomAD_exomes_AF[0] < {0}) "
                        "| (na dbNSFP_gnomAD_exomes_AC[0]) "
                        "| (dbNSFP_gnomAD_genomes_AF[0] < {0}) "
                        "| (na dbNSFP_gnomAD_genomes_AC[0])"
                    ),
                    "gnomad_ac_atmost": (
                        "((dbNSFP_gnomAD_exomes_AC[0] <= {0}) "
                        "| (na dbNSFP_gnomAD_exomes_AC[0])) "
                        "& ((dbNSFP_gnomAD_genomes_AC[0] <= {0}) "
                        "| (na dbNSFP_gnomAD_genomes_AC[0]))"
                    ),
                },
            },
            "dbnsfp5": {
                "description": "dbNSFP v5.x (joint exome+genome)",
                "fields_to_extract": "dbNSFP_gnomAD4.1_joint_AF",
                "hidden_columns": ["dbNSFP_gnomAD4.1_joint_AC"],
                "fragments": {
                    "gnomad_af_below": (
                        "((dbNSFP_gnomAD4.1_joint_AF[0] < {0}) | (na dbNSFP_gnomAD4.1_joint_AC[0]))"
                    ),
                    "gnomad_af_below_or": (
                        "(dbNSFP_gnomAD4.1_joint_AF[0] < {0}) | (na dbNSFP_gnomAD4.1_joint_AC[0])"
                    ),
                    "gnomad_ac_atmost": (
                        "((dbNSFP_gnomAD4.1_joint_AC[0] <= {0}) "
                        "| (na dbNSFP_gnomAD4.1_joint_AC[0]))"
                    ),
                },
            },
        },
        "fields_to_extract": "CHROM POS REF ALT",
        "html_report_default_hidden_columns": ["QUAL", "AC"],
        "presets": {
            "rare": "({{gnomad_af_below:0.0001}})",
            "coding": "((ANN[ANY].IMPACT has 'HIGH'))",
            "super_rare": "({{gnomad_ac_atmost:2}})",
            "mixed": "(({{gnomad_af_below_or:0.05}}) | (pathogenic))",
        },
    }


# ---------------------------------------------------------------------------
# FRAGMENT_PATTERN regex tests
# ---------------------------------------------------------------------------


class TestFragmentPattern:
    """Tests for the FRAGMENT_PATTERN regex."""

    def test_matches_simple_fragment(self):
        match = FRAGMENT_PATTERN.search("({{gnomad_af_below:0.05}})")
        assert match is not None
        assert match.group(1) == "gnomad_af_below"
        assert match.group(2) == "0.05"

    def test_matches_integer_param(self):
        match = FRAGMENT_PATTERN.search("({{gnomad_ac_atmost:2}})")
        assert match is not None
        assert match.group(1) == "gnomad_ac_atmost"
        assert match.group(2) == "2"

    def test_matches_multiple_fragments(self):
        text = "(({{frag1:a}}) & ({{frag2:b}}))"
        matches = FRAGMENT_PATTERN.findall(text)
        assert len(matches) == 2
        assert matches[0] == ("frag1", "a")
        assert matches[1] == ("frag2", "b")

    def test_no_match_plain_text(self):
        assert FRAGMENT_PATTERN.search("no templates here") is None

    def test_no_match_single_braces(self):
        assert FRAGMENT_PATTERN.search("{not_a_template:val}") is None


# ---------------------------------------------------------------------------
# _expand_fragments tests
# ---------------------------------------------------------------------------


class TestExpandFragments:
    """Tests for the _expand_fragments function."""

    def test_single_fragment(self):
        fragments = {"af": "AF < {0}"}
        result = _expand_fragments("({{af:0.01}})", fragments, "test")
        assert result == "(AF < 0.01)"

    def test_multiple_fragments(self):
        fragments = {"af": "AF < {0}", "ac": "AC <= {0}"}
        result = _expand_fragments("({{af:0.01}}) & ({{ac:5}})", fragments, "test")
        assert result == "(AF < 0.01) & (AC <= 5)"

    def test_same_fragment_different_params(self):
        fragments = {"af": "AF < {0}"}
        result = _expand_fragments("({{af:0.05}}) | ({{af:0.01}})", fragments, "test")
        assert result == "(AF < 0.05) | (AF < 0.01)"

    def test_no_fragments(self):
        fragments = {"af": "AF < {0}"}
        expr = "((ANN[ANY].IMPACT has 'HIGH'))"
        result = _expand_fragments(expr, fragments, "test")
        assert result == expr

    def test_unknown_fragment_raises(self):
        fragments = {"af": "AF < {0}"}
        with pytest.raises(ValueError, match="unknown fragment 'missing'"):
            _expand_fragments("({{missing:0.01}})", fragments, "test")

    def test_unknown_fragment_error_includes_preset_name(self):
        fragments = {}
        with pytest.raises(ValueError, match="Preset 'rare'"):
            _expand_fragments("({{af:0.01}})", fragments, "rare")

    def test_param_with_decimal(self):
        fragments = {"af": "AF[0] < {0}"}
        result = _expand_fragments("({{af:0.0001}})", fragments, "test")
        assert result == "(AF[0] < 0.0001)"

    def test_multiple_placeholders_in_fragment(self):
        """A single {0} placeholder used twice in the template."""
        fragments = {"dual": "(A < {0}) & (B < {0})"}
        result = _expand_fragments("({{dual:0.05}})", fragments, "test")
        assert result == "((A < 0.05) & (B < 0.05))"


# ---------------------------------------------------------------------------
# resolve_profile tests
# ---------------------------------------------------------------------------


class TestResolveProfile:
    """Tests for the resolve_profile function."""

    def test_default_profile_dbnsfp4(self, dbnsfp4_config):
        result = resolve_profile(dbnsfp4_config)
        assert result is dbnsfp4_config  # mutates in place

        # Presets should be expanded with dbnsfp4 fragments
        assert "dbNSFP_gnomAD_exomes_AF" in result["presets"]["rare"]
        assert "dbNSFP_gnomAD_genomes_AF" in result["presets"]["rare"]
        assert "{{" not in result["presets"]["rare"]

    def test_explicit_profile_dbnsfp5(self, dbnsfp4_config):
        dbnsfp4_config["field_profile"] = "dbnsfp5"
        result = resolve_profile(dbnsfp4_config)

        # Presets should be expanded with dbnsfp5 fragments
        assert "dbNSFP_gnomAD4.1_joint_AF" in result["presets"]["rare"]
        assert "dbNSFP_gnomAD_exomes_AF" not in result["presets"]["rare"]
        assert "{{" not in result["presets"]["rare"]

    def test_non_template_preset_unchanged(self, dbnsfp4_config):
        resolve_profile(dbnsfp4_config)
        assert dbnsfp4_config["presets"]["coding"] == "((ANN[ANY].IMPACT has 'HIGH'))"

    def test_fields_merged(self, dbnsfp4_config):
        resolve_profile(dbnsfp4_config)
        fields = dbnsfp4_config["fields_to_extract"]
        assert "CHROM" in fields
        assert "dbNSFP_gnomAD_exomes_AF" in fields
        assert "dbNSFP_gnomAD_genomes_AF" in fields

    def test_fields_merged_dbnsfp5(self, dbnsfp4_config):
        dbnsfp4_config["field_profile"] = "dbnsfp5"
        resolve_profile(dbnsfp4_config)
        fields = dbnsfp4_config["fields_to_extract"]
        assert "CHROM" in fields
        assert "dbNSFP_gnomAD4.1_joint_AF" in fields
        assert "dbNSFP_gnomAD_exomes_AF" not in fields

    def test_hidden_columns_merged(self, dbnsfp4_config):
        resolve_profile(dbnsfp4_config)
        hidden = dbnsfp4_config["html_report_default_hidden_columns"]
        assert "QUAL" in hidden  # base
        assert "dbNSFP_gnomAD_exomes_AC" in hidden  # profile

    def test_hidden_columns_no_duplicates(self, dbnsfp4_config):
        # Pre-populate with a profile column to verify no duplicate
        dbnsfp4_config["html_report_default_hidden_columns"].append("dbNSFP_gnomAD_exomes_AC")
        resolve_profile(dbnsfp4_config)
        hidden = dbnsfp4_config["html_report_default_hidden_columns"]
        assert hidden.count("dbNSFP_gnomAD_exomes_AC") == 1

    def test_missing_profile_raises(self, dbnsfp4_config):
        dbnsfp4_config["field_profile"] = "nonexistent"
        with pytest.raises(ValueError, match="not found"):
            resolve_profile(dbnsfp4_config)

    def test_missing_profile_lists_available(self, dbnsfp4_config):
        dbnsfp4_config["field_profile"] = "nonexistent"
        with pytest.raises(ValueError, match="dbnsfp4"):
            resolve_profile(dbnsfp4_config)

    def test_no_profiles_section_passthrough(self):
        """Config without field_profiles should pass through unchanged."""
        cfg = {
            "presets": {"rare": "(AF < 0.01)"},
            "fields_to_extract": "CHROM POS",
        }
        original_presets = cfg["presets"].copy()
        resolve_profile(cfg)
        assert cfg["presets"] == original_presets

    def test_empty_profiles_section_passthrough(self):
        """Config with empty field_profiles should pass through unchanged."""
        cfg = {
            "field_profiles": {},
            "presets": {"rare": "(AF < 0.01)"},
            "fields_to_extract": "CHROM POS",
        }
        original_presets = cfg["presets"].copy()
        resolve_profile(cfg)
        assert cfg["presets"] == original_presets

    def test_profile_with_empty_fields(self):
        """Profile with empty fields_to_extract should not add trailing space."""
        cfg = {
            "default_field_profile": "minimal",
            "field_profiles": {
                "minimal": {
                    "description": "Minimal profile",
                    "fields_to_extract": "",
                    "hidden_columns": [],
                    "fragments": {},
                },
            },
            "fields_to_extract": "CHROM POS",
            "html_report_default_hidden_columns": [],
            "presets": {},
        }
        resolve_profile(cfg)
        assert cfg["fields_to_extract"] == "CHROM POS"

    def test_mixed_preset_expansion(self, dbnsfp4_config):
        """Mixed preset with template and non-template parts."""
        resolve_profile(dbnsfp4_config)
        mixed = dbnsfp4_config["presets"]["mixed"]
        assert "dbNSFP_gnomAD_exomes_AF" in mixed
        assert "(pathogenic)" in mixed
        assert "{{" not in mixed


class TestResolveProfileWithRealConfig:
    """Tests using the actual config.json to verify backward compatibility."""

    def test_load_and_resolve_default_profile(self):
        """Default profile resolution produces valid expanded presets."""
        from variantcentrifuge.config import load_config

        cfg = load_config()
        resolve_profile(cfg)

        # All template markers should be expanded
        for name, expr in cfg["presets"].items():
            assert "{{" not in expr, f"Preset '{name}' still has unexpanded templates"

    def test_default_profile_rare_preset_structure(self):
        """The 'rare' preset under dbnsfp4 has the expected gnomAD fields."""
        from variantcentrifuge.config import load_config

        cfg = load_config()
        resolve_profile(cfg)

        rare = cfg["presets"]["rare"]
        assert "dbNSFP_gnomAD_exomes_AF[0] < 0.0001" in rare
        assert "dbNSFP_gnomAD_genomes_AF[0] < 0.0001" in rare
        assert "na dbNSFP_gnomAD_exomes_AC[0]" in rare

    def test_dbnsfp5_rare_preset_structure(self):
        """The 'rare' preset under dbnsfp5 has joint gnomAD fields."""
        from variantcentrifuge.config import load_config

        cfg = load_config()
        cfg["field_profile"] = "dbnsfp5"
        resolve_profile(cfg)

        rare = cfg["presets"]["rare"]
        assert "dbNSFP_gnomAD4.1_joint_AF[0] < 0.0001" in rare
        assert "dbNSFP_gnomAD_exomes_AF" not in rare

    def test_default_profile_fields_include_gnomad(self):
        """Default profile fields include gnomAD exome/genome fields."""
        from variantcentrifuge.config import load_config

        cfg = load_config()
        resolve_profile(cfg)

        fields = cfg["fields_to_extract"]
        assert "dbNSFP_gnomAD_exomes_AF" in fields
        assert "dbNSFP_gnomAD_genomes_AF" in fields
        assert "dbNSFP_gnomAD_exomes_AC" in fields
        assert "dbNSFP_gnomAD_genomes_AC" in fields

    def test_dbnsfp5_profile_fields_include_joint(self):
        """dbnsfp5 profile fields include joint gnomAD fields."""
        from variantcentrifuge.config import load_config

        cfg = load_config()
        cfg["field_profile"] = "dbnsfp5"
        resolve_profile(cfg)

        fields = cfg["fields_to_extract"]
        assert "dbNSFP_gnomAD4.1_joint_AF" in fields
        assert "dbNSFP_gnomAD4.1_joint_AC" in fields
        assert "dbNSFP_gnomAD_exomes_AF" not in fields

    def test_default_profile_hidden_columns(self):
        """Default profile adds gnomAD AC columns to hidden list."""
        from variantcentrifuge.config import load_config

        cfg = load_config()
        resolve_profile(cfg)

        hidden = cfg["html_report_default_hidden_columns"]
        assert "dbNSFP_gnomAD_exomes_AC" in hidden
        assert "dbNSFP_gnomAD_genomes_AC" in hidden

    def test_all_presets_have_balanced_parens(self):
        """Expanded presets with templates should have balanced parentheses.

        Note: Some presets have pre-existing unbalanced parentheses that are
        not caused by the field profile system. These are excluded here.
        """
        from variantcentrifuge.config import load_config

        cfg = load_config()
        resolve_profile(cfg)

        # Pre-existing paren issues (not from field profiles, not template presets)
        known_unbalanced = {
            "moderate_and_high_prediction",
            "high_or_lof_or_nmd",
            "high_or_pathogenic",
        }

        for name, expr in cfg["presets"].items():
            if name in known_unbalanced:
                continue
            open_count = expr.count("(")
            close_count = expr.count(")")
            assert open_count == close_count, (
                f"Preset '{name}' has unbalanced parentheses: "
                f"{open_count} open, {close_count} close"
            )

    def test_mutect2_preset_expansion(self):
        """mutect2 presets expand gnomAD AC fragments correctly."""
        from variantcentrifuge.config import load_config

        cfg = load_config()
        resolve_profile(cfg)

        for preset_name in ("mutect2_TvsN_pass", "mutect2_TvsN", "mutect2_To_pass"):
            expr = cfg["presets"][preset_name]
            assert "dbNSFP_gnomAD_exomes_AC[0] <= 2" in expr
            assert "dbNSFP_gnomAD_genomes_AC[0] <= 2" in expr
            assert "{{" not in expr

    def test_aif_preset_expansion(self):
        """aif preset expands both 0.05 and 0.01 gnomAD OR fragments."""
        from variantcentrifuge.config import load_config

        cfg = load_config()
        resolve_profile(cfg)

        aif = cfg["presets"]["aif"]
        assert "dbNSFP_gnomAD_exomes_AF[0] < 0.05" in aif
        assert "dbNSFP_gnomAD_exomes_AF[0] < 0.01" in aif
        assert "{{" not in aif


# ---------------------------------------------------------------------------
# list_profiles tests
# ---------------------------------------------------------------------------


class TestListProfiles:
    """Tests for the list_profiles function."""

    def test_list_profiles(self, dbnsfp4_config):
        profiles = list_profiles(dbnsfp4_config)
        assert len(profiles) == 2
        names = {p["name"] for p in profiles}
        assert names == {"dbnsfp4", "dbnsfp5"}

    def test_default_profile_marked(self, dbnsfp4_config):
        profiles = list_profiles(dbnsfp4_config)
        default = next(p for p in profiles if p["name"] == "dbnsfp4")
        assert "[default]" in default["description"]

    def test_non_default_profile_not_marked(self, dbnsfp4_config):
        profiles = list_profiles(dbnsfp4_config)
        other = next(p for p in profiles if p["name"] == "dbnsfp5")
        assert "[default]" not in other["description"]

    def test_empty_profiles(self):
        profiles = list_profiles({"field_profiles": {}})
        assert profiles == []

    def test_no_profiles_key(self):
        profiles = list_profiles({})
        assert profiles == []

    def test_list_profiles_from_real_config(self):
        """Real config should have at least dbnsfp4 and dbnsfp5."""
        from variantcentrifuge.config import load_config

        cfg = load_config()
        profiles = list_profiles(cfg)
        names = {p["name"] for p in profiles}
        assert "dbnsfp4" in names
        assert "dbnsfp5" in names
