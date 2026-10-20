"""Regression tests for Issue #168: optional subcommands and verified reference downloads."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from pathlib import Path
from typing import Any

import pytest

from variantcentrifuge import build_pm5_lookup as pm5_module
from variantcentrifuge import cli as cli_module
from variantcentrifuge.cli import main, parse_args


class TestSubcommandAndLegacyCompatibility:
    """Compare legacy and 'run' CLI invocations and verify option edge cases."""

    @pytest.mark.parametrize(
        "legacy_argv",
        [
            ["-v", "cohort.vcf.gz", "-g", "BRCA1", "--preset", "rare", "--xlsx"],
            ["--vcf-f", "abbreviated.vcf.gz", "--output-d", "out_abbr", "--html-r"],
            ["-vattached.vcf.gz", "-gTP53", "-oattached_out.tsv"],
            ["-v", "report", "-o", "run", "--output-dir", "validate-config"],
            ["-v", "build-pm5", "-g", "report"],
        ],
    )
    def test_parse_args_legacy_matches_run_subcommand(self, legacy_argv: list[str]) -> None:
        """Legacy invocations and 'run' subcommand invocations produce identical Namespace."""
        legacy_ns = parse_args(legacy_argv)
        run_ns = parse_args(["run", *legacy_argv])
        assert vars(legacy_ns) == vars(run_ns)

    def test_main_legacy_vs_run_invokes_pipeline_identically(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """main() with legacy args and 'run' subcommand pass equivalent args to run_pipeline."""
        captured: list[dict[str, Any]] = []

        monkeypatch.setattr(cli_module, "validate_required_tools", lambda *a, **kw: None)
        monkeypatch.setattr(cli_module, "validate_vcf_file", lambda *a, **kw: None)

        def fake_run_pipeline(args: Any) -> None:
            data = vars(args).copy()
            data.pop("start_time", None)
            captured.append(data)

        monkeypatch.setattr(cli_module, "run_pipeline", fake_run_pipeline)

        base_args = [
            "-v",
            "report",
            "-g",
            "BRCA1",
            "--filters",
            "PASS",
            "--output-dir",
            str(tmp_path / "out"),
        ]

        rc_legacy = main(base_args)
        rc_run = main(["run", *base_args])

        assert rc_legacy == 0
        assert rc_run == 0
        assert len(captured) == 2
        assert captured[0]["vcf_file"] == "report"
        assert captured[0] == captured[1]

    def test_utility_modes_preserve_behavior(
        self,
        monkeypatch: pytest.MonkeyPatch,
        capsys: pytest.CaptureFixture[str],
        tmp_path: Path,
    ) -> None:
        """Help, version, field-profile listing, and checkpoint-status keep exit behavior."""
        # --help
        with pytest.raises(SystemExit) as exc_help:
            main(["--help"])
        assert exc_help.value.code == 0
        out_help = capsys.readouterr().out
        assert "--vcf-file" in out_help

        # --version
        with pytest.raises(SystemExit) as exc_ver:
            main(["--version"])
        assert exc_ver.value.code == 0
        out_ver = capsys.readouterr().out
        assert "variantcentrifuge" in out_ver

        # --list-field-profiles
        with pytest.raises(SystemExit) as exc_prof:
            main(["--list-field-profiles"])
        assert exc_prof.value.code == 0
        out_prof = capsys.readouterr().out
        assert "dbnsfp4" in out_prof

        # --show-checkpoint-status
        with pytest.raises(SystemExit) as exc_ckpt:
            main(["--show-checkpoint-status", "--output-dir", str(tmp_path)])
        assert exc_ckpt.value.code == 0
        out_ckpt = capsys.readouterr().out
        assert "checkpoint" in out_ckpt.lower()

        # --list-stages
        with pytest.raises(SystemExit) as exc_stages:
            main(["--list-stages"])
        assert exc_stages.value.code == 0


class TestReportAndValidateConfigSubcommands:
    """Verify 'report' and 'validate-config' subcommands operate without the VCF pipeline."""

    def test_report_subcommand_generates_html_and_xlsx_without_pipeline(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """'variantcentrifuge report' builds HTML and Excel from TSV without run_pipeline."""

        def forbid_pipeline(*args: Any, **kwargs: Any) -> None:
            raise AssertionError("run_pipeline must not be called by report subcommand")

        monkeypatch.setattr(cli_module, "run_pipeline", forbid_pipeline)

        tsv_path = tmp_path / "variants.tsv"
        tsv_path.write_text(
            "CHROM\tPOS\tREF\tALT\tGENE\tIMPACT\tGT\nchr1\t100\tA\tG\tBRCA1\tHIGH\tSample1(0/1)\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "report_out"

        rc = main(
            [
                "report",
                "--input-tsv",
                str(tsv_path),
                "--output-dir",
                str(out_dir),
                "--html",
                "--xlsx",
            ]
        )
        assert rc == 0
        assert (out_dir / "report" / "index.html").is_file()
        assert (out_dir / "report" / "variants.json").is_file()
        assert (out_dir / "variants.xlsx").is_file()

    def test_report_subcommand_missing_tsv_returns_error(self, tmp_path: Path) -> None:
        """'variantcentrifuge report' returns non-zero when input TSV does not exist."""
        rc = main(
            [
                "report",
                "--input-tsv",
                str(tmp_path / "missing.tsv"),
                "--output-dir",
                str(tmp_path / "out"),
            ]
        )
        assert rc != 0

    def test_validate_config_subcommand_valid_and_invalid(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """'variantcentrifuge validate-config' validates configs and companion files."""

        def forbid_pipeline(*args: Any, **kwargs: Any) -> None:
            raise AssertionError("run_pipeline must not be called by validate-config")

        monkeypatch.setattr(cli_module, "run_pipeline", forbid_pipeline)

        # 1. Default package config is valid
        rc_default = main(["validate-config"])
        assert rc_default == 0
        assert "valid" in capsys.readouterr().out.lower()

        # 2. Invalid JSON config fails cleanly
        bad_cfg = tmp_path / "bad.json"
        bad_cfg.write_text("{broken json", encoding="utf-8")
        rc_bad = main(["validate-config", "-c", str(bad_cfg)])
        assert rc_bad != 0

        # 3. Unknown field profile fails cleanly
        rc_bad_prof = main(["validate-config", "--field-profile", "nonexistent_profile"])
        assert rc_bad_prof != 0

        # 4. Missing scoring config path fails cleanly
        rc_bad_score = main(
            ["validate-config", "--scoring-config-path", str(tmp_path / "no_scoring_dir")]
        )
        assert rc_bad_score != 0

        # 5. Unknown preset via CLI and via config file fails cleanly
        assert main(["validate-config", "--preset", "unknown_preset_xyz"]) == 1
        bad_preset_ref_cfg = tmp_path / "bad_preset_ref.json"
        bad_preset_ref_cfg.write_text(
            json.dumps({"preset": ["unknown_preset_in_cfg"]}), encoding="utf-8"
        )
        assert main(["validate-config", "-c", str(bad_preset_ref_cfg)]) == 1

        # 6. Config-defined sample_sheet missing or conflicting with ped_file fails cleanly
        missing_ss_cfg = tmp_path / "missing_ss.json"
        missing_ss_cfg.write_text(
            json.dumps({"sample_sheet": str(tmp_path / "no_such_sheet.tsv")}),
            encoding="utf-8",
        )
        assert main(["validate-config", "-c", str(missing_ss_cfg)]) == 1

        valid_ss = tmp_path / "valid_ss.tsv"
        valid_ss.write_text("sample_id\tsex\nS1\t1\n", encoding="utf-8")
        conflict_ss_cfg = tmp_path / "conflict_ss.json"
        conflict_ss_cfg.write_text(
            json.dumps(
                {
                    "sample_sheet": str(valid_ss),
                    "ped_file": str(tmp_path / "some.ped"),
                }
            ),
            encoding="utf-8",
        )
        assert main(["validate-config", "-c", str(conflict_ss_cfg)]) == 1

        # 7. Invalid scoring schema (non-dict variables) and invalid formula AST syntax fail cleanly
        bad_var_dir = tmp_path / "bad_scoring_vars"
        bad_var_dir.mkdir(parents=True, exist_ok=True)
        (bad_var_dir / "variable_assignment_config.json").write_text(
            json.dumps({"variables": "not_a_dict"}), encoding="utf-8"
        )
        (bad_var_dir / "formula_config.json").write_text(
            json.dumps({"formulas": [{"score": "a + 1"}]}), encoding="utf-8"
        )
        assert main(["validate-config", "--scoring-config-path", str(bad_var_dir)]) == 1

        bad_ast_dir = tmp_path / "bad_scoring_ast"
        bad_ast_dir.mkdir(parents=True, exist_ok=True)
        (bad_ast_dir / "variable_assignment_config.json").write_text(
            json.dumps({"variables": {"IMPACT": "a|default:0"}}), encoding="utf-8"
        )
        (bad_ast_dir / "formula_config.json").write_text(
            json.dumps({"formulas": [{"score": "a + * b"}]}), encoding="utf-8"
        )
        assert main(["validate-config", "--scoring-config-path", str(bad_ast_dir)]) == 1

    def test_build_pm5_subcommand_delegates_to_build_pm5_main(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """'variantcentrifuge build-pm5' forwards arguments to build_pm5_lookup.main."""
        forwarded: list[list[str] | None] = []

        def fake_pm5_main(args: list[str] | None = None) -> int:
            forwarded.append(args)
            return 0

        monkeypatch.setattr(pm5_module, "main", fake_pm5_main)
        out_json = str(tmp_path / "pm5.json")
        rc = main(["build-pm5", "--download", "--output", out_json])
        assert rc == 0
        assert forwarded == [["--download", "--output", out_json]]

        # download-references automatically prepends --download when not already present
        forwarded.clear()
        rc_dl = main(["download-references", "--output", out_json])
        assert rc_dl == 0
        assert forwarded == [["--download", "--output", out_json]]

    def test_validate_config_respects_configured_phenotype_columns(self, tmp_path: Path) -> None:
        """validate-config uses phenotype_sample_column and phenotype_value_column from config."""
        pheno_file = tmp_path / "custom_pheno.tsv"
        pheno_file.write_text("SampleID\tHPO_Terms\nS1\tHP:0001250\n", encoding="utf-8")
        cfg_file = tmp_path / "custom_pheno_cfg.json"
        cfg_file.write_text(
            json.dumps(
                {
                    "phenotype_file": str(pheno_file),
                    "phenotype_sample_column": "SampleID",
                    "phenotype_value_column": "HPO_Terms",
                }
            ),
            encoding="utf-8",
        )
        assert main(["validate-config", "-c", str(cfg_file)]) == 0

    def test_comma_separated_preset_in_run_and_resolve(self) -> None:
        """--preset rare,coding splits on commas in resolve_presets_and_tumor_normal."""
        from variantcentrifuge.cli_config import resolve_presets_and_tumor_normal
        from variantcentrifuge.config import load_config
        from variantcentrifuge.field_profile import resolve_profile

        parser = cli_module.create_parser()
        args = parser.parse_args(["-v", "in.vcf.gz", "--preset", "rare,coding"])
        cfg = load_config()
        resolve_profile(cfg)
        filt, _ = resolve_presets_and_tumor_normal(args, cfg, None)
        assert filt is not None
        assert cfg["preset"] == ["rare", "coding"]

    def test_shortened_utility_flags_and_interactive_resume_validation(
        self, tmp_path: Path
    ) -> None:
        """Shortened utility flags and unvalidated interactive-resume abort via parser.error."""
        with pytest.raises(SystemExit) as exc_short_vcf:
            main(["-v", "x.vcf.gz", "--preset", "rare", "--show-vcf"])
        assert exc_short_vcf.value.code == 2

        with pytest.raises(SystemExit) as exc_short_prof:
            main(["--list-prof"])
        assert exc_short_prof.value.code == 2

        # --interactive-resume without --enable-checkpoint fails validation with code 2
        with pytest.raises(SystemExit) as exc_ir_no_ckpt:
            main(["-v", "x.vcf.gz", "--interactive-resume", "--output-dir", str(tmp_path)])
        assert exc_ir_no_ckpt.value.code == 2


def _fake_urlopen(payload: bytes, error: BaseException | None = None) -> Any:
    """Return a urlopen stand-in streaming payload, then raising error if given."""

    class _Stream(io.BytesIO):
        def read(self, size: int | None = -1) -> bytes:
            data = super().read(size)
            if not data and error is not None:
                raise error
            return data

    def fake(url: str, timeout: float | None = None) -> io.BytesIO:
        assert timeout == 60
        return _Stream(payload)

    return fake


class TestVerifiedReferenceDownloads:
    """Verify reference download integrity, atomic publication, and provenance."""

    @staticmethod
    def _make_clinvar_gzip_bytes() -> bytes:
        header = (
            "#AlleleID\tType\tName\tGeneSymbol\tClinicalSignificance\t"
            "Origin\tAssembly\tReviewStatus\n"
        )
        row = (
            "1\tsingle nucleotide variant\tNM_000059.4(BRCA2):c.8167G>C (p.Asp2723His)\t"
            "BRCA2\tPathogenic\tgermline\tGRCh38\t"
            "criteria provided, multiple submitters, no conflicts\n"
        )
        return gzip.compress((header + row).encode("utf-8"))

    def test_matching_sha256_publishes_and_records_provenance(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Matching SHA-256 succeeds and records digest, URL, and timestamp in _meta."""
        payload = self._make_clinvar_gzip_bytes()
        expected_digest = hashlib.sha256(payload).hexdigest()

        monkeypatch.setattr(pm5_module.urllib.request, "urlopen", _fake_urlopen(payload))

        dest_gz = tmp_path / "variant_summary.txt.gz"
        out_json = tmp_path / "pm5_lookup.json"

        rc = pm5_module.main(
            [
                "--download",
                "--input",
                str(dest_gz),
                "--expected-sha256",
                expected_digest,
                "--output",
                str(out_json),
            ]
        )
        assert rc == 0
        assert dest_gz.read_bytes() == payload
        data = json.loads(out_json.read_text(encoding="utf-8"))
        meta = data["_meta"]
        assert meta["sha256"] == expected_digest
        assert meta["source_url"] == pm5_module.CLINVAR_FTP_URL
        assert "timestamp" in meta
        assert "provenance" not in meta

    def test_mismatching_sha256_preserves_existing_valid_file(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Mismatching SHA-256 raises ValueError and leaves pre-existing target file intact."""
        original_payload = self._make_clinvar_gzip_bytes()
        dest_gz = tmp_path / "variant_summary.txt.gz"
        dest_gz.write_bytes(original_payload)

        corrupt_payload = gzip.compress(b"corrupted new payload")

        monkeypatch.setattr(pm5_module.urllib.request, "urlopen", _fake_urlopen(corrupt_payload))

        with pytest.raises(ValueError, match="SHA-256"):
            pm5_module.download_variant_summary(
                str(dest_gz),
                expected_sha256="0" * 64,
            )

        # Pre-existing valid file must remain untouched and no leftover temp files remain
        assert dest_gz.read_bytes() == original_payload
        leftover_tmps = [p for p in tmp_path.iterdir() if p != dest_gz]
        assert leftover_tmps == []

    def test_interrupted_download_preserves_existing_valid_file(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Interrupted download cleans up temp file and preserves pre-existing target file."""
        original_payload = self._make_clinvar_gzip_bytes()
        dest_gz = tmp_path / "variant_summary.txt.gz"
        dest_gz.write_bytes(original_payload)

        monkeypatch.setattr(
            pm5_module.urllib.request,
            "urlopen",
            _fake_urlopen(b"\x1f\x8b\x08\x00partial_stream", KeyboardInterrupt()),
        )

        with pytest.raises(KeyboardInterrupt):
            pm5_module.download_variant_summary(str(dest_gz))

        assert dest_gz.read_bytes() == original_payload
        leftover_tmps = [p for p in tmp_path.iterdir() if p != dest_gz]
        assert leftover_tmps == []

    def test_truncated_gzip_download_rejected_and_preserves_existing_file(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Truncated gzip stream fails verification and preserves pre-existing target file."""
        original_payload = self._make_clinvar_gzip_bytes()
        dest_gz = tmp_path / "variant_summary.txt.gz"
        dest_gz.write_bytes(original_payload)

        monkeypatch.setattr(
            pm5_module.urllib.request, "urlopen", _fake_urlopen(original_payload[:-10])
        )

        with pytest.raises(ValueError, match=r"gzip|corrupt|truncated"):
            pm5_module.download_variant_summary(str(dest_gz))

        assert dest_gz.read_bytes() == original_payload
        leftover_tmps = [p for p in tmp_path.iterdir() if p != dest_gz]
        assert leftover_tmps == []

    def test_local_input_sha256_and_stats_and_invalid_gzip_headers(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Verify local --input with --expected-sha256, --stats, and invalid gzip checks."""
        payload = self._make_clinvar_gzip_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        local_gz = tmp_path / "local.txt.gz"
        local_gz.write_bytes(payload)
        out_gz = tmp_path / "lookup.json.gz"

        # Matching local sha256 + --stats + .json.gz output
        rc = pm5_module.main(
            [
                "--input",
                str(local_gz),
                "--expected-sha256",
                digest,
                "--output",
                str(out_gz),
                "--stats",
            ]
        )
        assert rc == 0
        assert "PM5 Lookup Table Statistics" in capsys.readouterr().out
        with gzip.open(out_gz, "rt", encoding="utf-8") as fh:
            assert json.load(fh)["_meta"]["source_url"] is None

        # Mismatching local sha256
        with pytest.raises(ValueError, match="SHA-256"):
            pm5_module.main(
                [
                    "--input",
                    str(local_gz),
                    "--expected-sha256",
                    "f" * 64,
                    "--output",
                    str(out_gz),
                ]
            )

        # verify_gzip_integrity edge cases: missing/empty and non-gzip magic
        with pytest.raises(ValueError, match="empty or missing"):
            pm5_module.verify_gzip_integrity(tmp_path / "does_not_exist.gz")

        not_gz = tmp_path / "plain.gz"
        not_gz.write_bytes(b"not a gzip header")
        with pytest.raises(ValueError, match="not a valid gzip"):
            pm5_module.verify_gzip_integrity(not_gz)

    def test_subcommand_edge_branches(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Cover validate-config ped/phenotype/presets and utility options."""
        from variantcentrifuge import vcf_header_parser
        from variantcentrifuge.cli_subcommands import (
            dispatch_subcommand,
            handle_early_utility_options,
            handle_listing_and_resume_options,
        )

        # Unsupported subcommand raises ValueError
        with pytest.raises(ValueError, match="Unsupported subcommand"):
            dispatch_subcommand("unknown-cmd", [])

        # validate-config with valid ped and phenotype file
        ped_file = tmp_path / "fam.ped"
        ped_file.write_text("FAM1\tS1\t0\t0\t1\t2\n", encoding="utf-8")
        pheno_file = tmp_path / "pheno.tsv"
        pheno_file.write_text("sample_id\tphenotype\nS1\tHP:0001250\n", encoding="utf-8")

        assert (
            main(
                [
                    "validate-config",
                    "--ped",
                    str(ped_file),
                    "--phenotype-file",
                    str(pheno_file),
                ]
            )
            == 0
        )

        # validate-config with missing ped file
        assert main(["validate-config", "--ped", str(tmp_path / "missing.ped")]) != 0

        # validate-config with non-dict presets or field_profiles
        bad_presets = tmp_path / "bad_presets.json"
        bad_presets.write_text(json.dumps({"presets": ["not_a_dict"]}), encoding="utf-8")
        assert main(["validate-config", "-c", str(bad_presets)]) != 0

        bad_fp = tmp_path / "bad_fp.json"
        bad_fp.write_text(json.dumps({"field_profiles": ["not_a_dict"]}), encoding="utf-8")
        assert main(["validate-config", "-c", str(bad_fp)]) != 0

        # report with explicit --output-file on .tsv.gz
        tsv_gz = tmp_path / "input.tsv.gz"
        with gzip.open(tsv_gz, "wt", encoding="utf-8") as f:
            f.write("CHROM\tPOS\tREF\tALT\tGENE\nchr1\t10\tA\tT\tBRCA1\n")
        explicit_xlsx = tmp_path / "custom.xlsx"
        assert (
            main(
                [
                    "report",
                    "-i",
                    str(tsv_gz),
                    "-o",
                    str(tmp_path / "rep_out"),
                    "--output-file",
                    str(explicit_xlsx),
                    "--xlsx",
                ]
            )
            == 0
        )
        assert explicit_xlsx.is_file()

        # --show-vcf-annotations (table and json) with mocked _read_header (no external bcftools)
        vcf_txt = tmp_path / "header.vcf"
        vcf_header_text = (
            "##fileformat=VCFv4.2\n"
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        vcf_txt.write_text(vcf_header_text, encoding="utf-8")
        monkeypatch.setattr(vcf_header_parser, "_read_header", lambda _p: vcf_header_text)

        with pytest.raises(SystemExit) as exc_ann1:
            handle_early_utility_options(["--show-vcf-annotations", "-v", str(vcf_txt)])
        assert exc_ann1.value.code == 0

        with pytest.raises(SystemExit) as exc_ann2:
            handle_early_utility_options(
                [
                    "--show-vcf-annotations",
                    "-v",
                    str(vcf_txt),
                    "--annotation-format",
                    "json",
                ]
            )
        assert exc_ann2.value.code == 0

        # --list-checkpoints when no checkpoint exists
        import argparse

        ns_ckpt = argparse.Namespace(
            list_stages=False,
            list_checkpoints=True,
            interactive_resume=False,
            output_dir=str(tmp_path / "no_ckpt"),
        )
        with pytest.raises(SystemExit) as exc_lc:
            handle_listing_and_resume_options(ns_ckpt, {})
        assert exc_lc.value.code == 0

        # --interactive-resume when no checkpoint exists
        ns_ir = argparse.Namespace(
            list_stages=False,
            list_checkpoints=False,
            interactive_resume=True,
            output_dir=str(tmp_path / "no_ckpt"),
        )
        with pytest.raises(SystemExit) as exc_ir:
            handle_listing_and_resume_options(ns_ir, {})
        assert exc_ir.value.code == 1
