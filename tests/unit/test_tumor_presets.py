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

    def test_mutect2_presets_use_template_vars(self, config: dict) -> None:
        """mutect2_* presets use {tumor_idx}/{normal_idx} and explicit per-sample DP."""
        for name in ("mutect2_TvsN_pass", "mutect2_TvsN"):
            expr = config["presets"][name]
            assert "{tumor_idx}" in expr
            assert "{normal_idx}" in expr
            expanded = self._expand(expr, tumor_idx="1", normal_idx="0")
            assert "GEN[1].DP >= 50" in expanded
            assert "GEN[0].DP >= 20" in expanded
            assert "GEN[1].AF >= 0.05" in expanded
            assert "GEN[0].AF < 0.03" in expanded

        to_expr = config["presets"]["mutect2_To_pass"]
        assert "{tumor_idx}" in to_expr
        assert "{normal_idx}" not in to_expr
        to_expanded = self._expand(to_expr, tumor_idx="0")
        assert "GEN[0].DP >= 50" in to_expanded
        assert "GEN[0].AF >= 0.05" in to_expanded

    def test_mutect2_tvsn_pass_parity_with_somatic_pass(self, config: dict) -> None:
        """mutect2_TvsN_pass and somatic_pass target the same tumor/normal sample indices."""
        for t_idx, n_idx in (("1", "0"), ("0", "1")):
            m_exp = self._expand(
                config["presets"]["mutect2_TvsN_pass"],
                tumor_idx=t_idx,
                normal_idx=n_idx,
            )
            s_exp = self._expand(
                config["presets"]["somatic_pass"],
                tumor_idx=t_idx,
                normal_idx=n_idx,
                tumor_dp_min="50",
                normal_dp_min="20",
            )
            for clause in (
                f"GEN[{t_idx}].DP >= 50",
                f"GEN[{n_idx}].DP >= 20",
                f"GEN[{t_idx}].AF >= 0.05",
                f"GEN[{n_idx}].AF < 0.03",
            ):
                assert clause in m_exp
                assert clause in s_exp


# ---------------------------------------------------------------------------
# VCF header ##tumor_sample / ##normal_sample resolution
# ---------------------------------------------------------------------------
class TestHeaderSampleResolution:
    def test_normal_first_vcf_header(self) -> None:
        from variantcentrifuge.vcf_header_parser import resolve_tumor_normal_indices

        hdr = (
            "##fileformat=VCFv4.2\n"
            "##tumor_sample=APA1_T\n"
            "##normal_sample=APA1_N\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAPA1_N\tAPA1_T\n"
        )
        assert resolve_tumor_normal_indices(header_text=hdr) == ("1", "0")

    def test_tumor_first_vcf_header(self) -> None:
        from variantcentrifuge.vcf_header_parser import resolve_tumor_normal_indices

        hdr = (
            "##fileformat=VCFv4.2\n"
            "##tumor_sample=APA1_T\n"
            "##normal_sample=APA1_N\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAPA1_T\tAPA1_N\n"
        )
        assert resolve_tumor_normal_indices(header_text=hdr) == ("0", "1")

    def test_tumor_only_single_sample_vcf_header(self) -> None:
        from variantcentrifuge.vcf_header_parser import resolve_tumor_normal_indices

        hdr = (
            "##fileformat=VCFv4.2\n"
            "##tumor_sample=APA51_T\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAPA51_T\n"
        )
        t_idx, _ = resolve_tumor_normal_indices(header_text=hdr)
        assert t_idx == "0"

    def test_mismatched_tumor_sample_raises(self) -> None:
        from variantcentrifuge.vcf_header_parser import resolve_tumor_normal_indices

        hdr = (
            "##fileformat=VCFv4.2\n"
            "##tumor_sample=WRONG_T\n"
            "##normal_sample=APA1_N\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAPA1_N\tAPA1_T\n"
        )
        with pytest.raises(ValueError, match="##tumor_sample"):
            resolve_tumor_normal_indices(header_text=hdr)

    def test_conflicting_cli_index_and_header_raises(self) -> None:
        from variantcentrifuge.vcf_header_parser import resolve_tumor_normal_indices

        hdr = (
            "##fileformat=VCFv4.2\n"
            "##tumor_sample=APA1_T\n"
            "##normal_sample=APA1_N\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAPA1_N\tAPA1_T\n"
        )
        with pytest.raises(ValueError, match="conflicts with"):
            resolve_tumor_normal_indices(header_text=hdr, cli_tumor_idx=0)

    def test_out_of_range_index_raises(self) -> None:
        from variantcentrifuge.vcf_header_parser import resolve_tumor_normal_indices

        hdr = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAPA1_N\tAPA1_T\n"
        )
        with pytest.raises(ValueError, match="out of range"):
            resolve_tumor_normal_indices(header_text=hdr, cli_tumor_idx=2)

    def test_get_vcf_names_raises_on_bcftools_failure(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from variantcentrifuge.helpers import get_vcf_names

        class DummyProc:
            returncode = 1

            def communicate(self) -> tuple[str, str]:
                return "", "Failed to open file"

        monkeypatch.setattr("subprocess.Popen", lambda *_a, **_kw: DummyProc())
        with pytest.raises(RuntimeError, match="bcftools query -l failed"):
            get_vcf_names("nonexistent.vcf")

    def test_cli_expands_tumor_idx_in_final_filter(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: pytest.TempPathFactory
    ) -> None:
        from variantcentrifuge import cli as cli_module

        vcf = tmp_path / "pair.vcf"  # type: ignore[operator]
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##tumor_sample=TUM\n"
            "##normal_sample=NRM\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUM\tNRM\n"
        )
        captured: list[dict] = []
        monkeypatch.setattr(cli_module, "run_pipeline", lambda a: captured.append(a.config.copy()))
        monkeypatch.setattr(
            "sys.argv",
            [
                "variantcentrifuge",
                "-v",
                str(vcf),
                "--preset",
                "somatic_pass",
                "--final-filter",
                "GEN[{tumor_idx}].AF >= 0.1 & GEN[{normal_idx}].AF <= 0.01",
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert cli_module.main() == 0
        assert captured[0]["tumor_sample_index"] == "0"
        assert captured[0]["normal_sample_index"] == "1"
        assert captured[0]["final_filter"] == "GEN[0].AF >= 0.1 & GEN[1].AF <= 0.01"


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
