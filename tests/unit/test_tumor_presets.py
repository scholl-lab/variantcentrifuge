"""Unit tests for tumor-normal filtering presets and template expansion."""

from __future__ import annotations

from typing import ClassVar

import pytest

from variantcentrifuge.cli import create_parser
from variantcentrifuge.config import load_config
from variantcentrifuge.field_profile import resolve_profile


@pytest.fixture()
def config() -> dict:
    """Load the default config with profile resolution."""
    cfg = load_config(None)
    resolve_profile(cfg)
    return cfg


# ---------------------------------------------------------------------------
# Preset existence
# ---------------------------------------------------------------------------
class TestTumorPresetsDefined:
    EXPECTED: ClassVar[list[str]] = [
        "somatic",
        "somatic_pass",
        "somatic_strict",
        "loh",
        "germline_shared",
        "tumor_only",
    ]

    def test_all_somatic_presets_exist(self, config: dict) -> None:
        presets = config.get("presets", {})
        for name in self.EXPECTED:
            assert name in presets, f"Preset '{name}' missing from config"

    def test_legacy_mutect2_presets_preserved(self, config: dict) -> None:
        presets = config.get("presets", {})
        assert "mutect2_TvsN_pass" in presets
        assert "mutect2_TvsN" in presets
        assert "mutect2_To_pass" in presets


# ---------------------------------------------------------------------------
# Template variable expansion
# ---------------------------------------------------------------------------
class TestPresetExpansion:
    def _expand(self, expr: str, **overrides: str) -> str:
        """Expand tumor-normal template variables in a preset expression."""
        defaults = {
            "tumor_idx": "1",
            "normal_idx": "0",
            "tumor_dp_min": "20",
            "normal_dp_min": "20",
            "tumor_af_min": "0.05",
            "normal_af_max": "0.03",
        }
        defaults.update(overrides)
        return expr.format_map(defaults)

    def test_somatic_default_indices(self, config: dict) -> None:
        expr = config["presets"]["somatic"]
        expanded = self._expand(expr)
        assert "GEN[0].DP >= 20" in expanded  # normal
        assert "GEN[1].DP >= 20" in expanded  # tumor
        assert "GEN[0].AF < 0.03" in expanded  # normal AF max
        assert "GEN[1].AF >= 0.05" in expanded  # tumor AF min

    def test_somatic_custom_indices(self, config: dict) -> None:
        expr = config["presets"]["somatic"]
        expanded = self._expand(expr, tumor_idx="0", normal_idx="1")
        assert "GEN[1].DP >= 20" in expanded  # normal at index 1
        assert "GEN[0].DP >= 20" in expanded  # tumor at index 0
        assert "GEN[1].AF < 0.03" in expanded  # normal at index 1
        assert "GEN[0].AF >= 0.05" in expanded  # tumor at index 0

    def test_somatic_custom_thresholds(self, config: dict) -> None:
        expr = config["presets"]["somatic"]
        expanded = self._expand(
            expr, tumor_dp_min="50", normal_dp_min="30", tumor_af_min="0.01", normal_af_max="0.01"
        )
        assert "GEN[0].DP >= 30" in expanded
        assert "GEN[1].DP >= 50" in expanded
        assert "GEN[0].AF < 0.01" in expanded
        assert "GEN[1].AF >= 0.01" in expanded

    def test_somatic_pass_includes_filter(self, config: dict) -> None:
        expr = config["presets"]["somatic_pass"]
        expanded = self._expand(expr)
        assert "FILTER = 'PASS'" in expanded

    def test_somatic_strict_hardcoded_thresholds(self, config: dict) -> None:
        expr = config["presets"]["somatic_strict"]
        expanded = self._expand(expr)
        assert "GEN[0].DP >= 50" in expanded  # hardcoded in strict
        assert "GEN[1].DP >= 50" in expanded

    def test_loh_genotype_pattern(self, config: dict) -> None:
        expr = config["presets"]["loh"]
        expanded = self._expand(expr)
        assert "GEN[0].GT = '0/1'" in expanded  # normal het
        assert "GEN[1].GT = '1/1'" in expanded  # tumor hom-alt
        assert "GEN[1].AF >= 0.85" in expanded  # allelic imbalance

    def test_germline_shared_pattern(self, config: dict) -> None:
        expr = config["presets"]["germline_shared"]
        expanded = self._expand(expr)
        assert "GEN[0].AF >= 0.20" in expanded  # normal
        assert "GEN[1].AF >= 0.20" in expanded  # tumor

    def test_tumor_only_no_normal_refs(self, config: dict) -> None:
        expr = config["presets"]["tumor_only"]
        expanded = self._expand(expr)
        # tumor_only should not reference normal sample index
        assert "GEN[0]" not in expanded or "{normal_idx}" not in expr
        assert "GEN[1].DP >= 20" in expanded
        assert "GEN[1].AF >= 0.05" in expanded

    def test_legacy_presets_no_template_vars(self, config: dict) -> None:
        """Legacy mutect2 presets should have no {tumor_idx} variables."""
        for name in ("mutect2_TvsN_pass", "mutect2_TvsN", "mutect2_To_pass"):
            expr = config["presets"][name]
            assert "{tumor_idx}" not in expr
            assert "{normal_idx}" not in expr


# ---------------------------------------------------------------------------
# CLI argument parsing
# ---------------------------------------------------------------------------
class TestTumorCLIArgs:
    def test_tumor_sample_index_default(self) -> None:
        parser = create_parser()
        args = parser.parse_args(["-v", "test.vcf"])
        assert args.tumor_sample_index is None

    def test_tumor_sample_index_custom(self) -> None:
        parser = create_parser()
        args = parser.parse_args(["-v", "test.vcf", "--tumor-sample-index", "0"])
        assert args.tumor_sample_index == 0

    def test_normal_sample_index_custom(self) -> None:
        parser = create_parser()
        args = parser.parse_args(["-v", "test.vcf", "--normal-sample-index", "2"])
        assert args.normal_sample_index == 2

    def test_tumor_dp_min(self) -> None:
        parser = create_parser()
        args = parser.parse_args(["-v", "test.vcf", "--tumor-dp-min", "50"])
        assert args.tumor_dp_min == 50

    def test_normal_dp_min(self) -> None:
        parser = create_parser()
        args = parser.parse_args(["-v", "test.vcf", "--normal-dp-min", "30"])
        assert args.normal_dp_min == 30

    def test_tumor_af_min(self) -> None:
        parser = create_parser()
        args = parser.parse_args(["-v", "test.vcf", "--tumor-af-min", "0.01"])
        assert args.tumor_af_min == 0.01

    def test_normal_af_max(self) -> None:
        parser = create_parser()
        args = parser.parse_args(["-v", "test.vcf", "--normal-af-max", "0.05"])
        assert args.normal_af_max == 0.05
