"""Regression tests for Issue #167:
- Organize variantcentrifuge --help into task-oriented argument groups
- Demote deprecated/internal flags to argparse.SUPPRESS or Advanced group
- Add curated examples epilog
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from variantcentrifuge.checkpoint_models import AtomicFileOperation
from variantcentrifuge.cli import create_parser
from variantcentrifuge.converter import append_tsv_as_sheet, convert_to_excel
from variantcentrifuge.generate_html_report import generate_html_report
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.run_summary import (
    SUMMARY_SCHEMA_VERSION,
    write_run_summary,
)
from variantcentrifuge.stages.analysis_stages import (
    AssociationAnalysisStage,
    GeneBurdenAnalysisStage,
)
from variantcentrifuge.stages.tsv_output_stages import TSVOutputStage


@pytest.mark.unit
class TestIssue167HelpOrganization:
    """Verify task-oriented argument groups and examples epilog in --help."""

    def test_help_contains_task_oriented_groups(self) -> None:
        """--help output must organize arguments into task-oriented groups."""
        parser = create_parser()
        help_text = parser.format_help()

        expected_group_keywords = [
            "Core Input/Output",
            "Gene",
            "Filtering",
            "Inheritance",
            "Reporting & Visualization",
            "Checkpoint & Resume",
            "Privacy",
            "Performance",
            "Advanced & Tuning Options",
        ]

        for keyword in expected_group_keywords:
            assert any(
                keyword.lower() in group.title.lower()
                for group in parser._action_groups
                if group.title
            ), f"Expected group matching '{keyword}' in argument groups"
            assert keyword.lower() in help_text.lower(), f"Expected '{keyword}' in help text"

    def test_help_contains_curated_examples_epilog(self) -> None:
        """--help output must include a curated examples epilog."""
        parser = create_parser()
        help_text = parser.format_help()

        assert "Examples:" in help_text
        assert "variantcentrifuge -g" in help_text
        assert "--preset rare,coding" in help_text

    def test_internal_tuning_flags_demoted_to_advanced_group(self) -> None:
        """Internal tuning flags must be located in the Advanced/Internal group."""
        parser = create_parser()
        advanced_groups = [g for g in parser._action_groups if g.title and "Advanced" in g.title]
        assert len(advanced_groups) == 1, "Expected exactly one Advanced/Internal options group"

        advanced_group = advanced_groups[0]
        advanced_flags = {
            opt for action in advanced_group._group_actions for opt in action.option_strings
        }

        # Check internal knobs are in the advanced group
        internal_knobs = [
            "--force-chunked-processing",
            "--force-inheritance-processing",
            "--memory-safety-factor",
            "--inheritance-memory-fraction",
            "--sample-column-creation-method",
        ]
        for knob in internal_knobs:
            assert knob in advanced_flags, f"Expected internal flag '{knob}' in Advanced group"

    def test_deprecated_flags_are_suppressed_from_help(self) -> None:
        """Deprecated flags must be hidden from --help but remain functional."""
        parser = create_parser()
        help_text = parser.format_help()

        # Deprecated aliases should NOT appear in help
        assert "--pca-file" not in help_text
        assert "--pca-tool" not in help_text
        assert "--ped-file" not in help_text

        # But they must still parse without error (backward compatibility)
        parsed_pca = parser.parse_args(["-v", "test.vcf", "--pca-file", "my.eigenvec"])
        assert parsed_pca.pca == "my.eigenvec"

        parsed_ped = parser.parse_args(["-v", "test.vcf", "--ped-file", "family.ped"])
        assert parsed_ped.ped == "family.ped"


@pytest.mark.unit
class TestIssue167AtomicOutputsAndRunSummary:
    """Verify atomic publication of final outputs and consolidated run_summary.json."""

    def test_atomic_file_operation_failure_preserves_existing_file(self, tmp_path: Path) -> None:
        target = tmp_path / "result.tsv"
        target.write_text("original valid content")

        with (
            pytest.raises(RuntimeError, match="crash during write"),
            AtomicFileOperation(target) as tmp,
            open(tmp, "w", encoding="utf-8") as f,
        ):
            f.write("corrupt partial content")
            raise RuntimeError("crash during write")

        assert target.read_text() == "original valid content"
        tmp_files = list(tmp_path.glob("result.tsv.tmp.*"))
        assert len(tmp_files) == 0

    def test_atomic_file_operation_failure_on_new_file_leaves_no_file(self, tmp_path: Path) -> None:
        target = tmp_path / "new_result.tsv"

        with (
            pytest.raises(RuntimeError, match="disk error"),
            AtomicFileOperation(target) as tmp,
            open(tmp, "w", encoding="utf-8") as f,
        ):
            f.write("incomplete content")
            raise RuntimeError("disk error")

        assert not target.exists()
        tmp_files = [f for f in tmp_path.glob("*") if ".tmp." in f.name]
        assert len(tmp_files) == 0

    def test_atomic_file_operation_subprocess_kill_leaves_no_partial_artifact(
        self, tmp_path: Path
    ) -> None:
        target_file = tmp_path / "never_published.tsv"
        ready_file = tmp_path / "ready.flag"

        code = (
            "import time\n"
            "from pathlib import Path\n"
            "from variantcentrifuge.checkpoint_models import AtomicFileOperation\n"
            f"target = Path({str(target_file)!r})\n"
            f"ready = Path({str(ready_file)!r})\n"
            "with AtomicFileOperation(target) as tmp:\n"
            "    with open(tmp, 'w') as f:\n"
            "        f.write('partial content')\n"
            "    ready.write_text('ready')\n"
            "    time.sleep(30)\n"
        )
        proc = subprocess.Popen([sys.executable, "-c", code])
        try:
            for _ in range(50):
                if ready_file.exists():
                    break
                time.sleep(0.1)
            assert ready_file.exists(), "Subprocess did not start AtomicFileOperation in time"

            proc.kill()
            proc.wait(timeout=5)

            assert not target_file.exists()
        finally:
            if proc.poll() is None:
                proc.kill()

    def test_tsv_output_stage_atomic_write_failure_preserves_original(self, tmp_path: Path) -> None:
        output_file = tmp_path / "output.tsv"
        output_file.write_text("existing header\n1\t100\tA\tT\n")

        stage = TSVOutputStage()
        context = PipelineContext(
            args=argparse.Namespace(),
            config={"final_output": str(output_file)},
            workspace=Workspace(tmp_path, "test"),
        )
        context.current_dataframe = pd.DataFrame([{"CHROM": "chr1", "POS": 100}])

        with (
            patch.object(pd.DataFrame, "to_csv", side_effect=OSError("Write failed mid-stream")),
            pytest.raises(IOError, match="Write failed mid-stream"),
        ):
            stage._process(context)

        assert output_file.read_text() == "existing header\n1\t100\tA\tT\n"
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_excel_output_atomic_failure_preserves_original(self, tmp_path: Path) -> None:
        tsv_file = tmp_path / "input.tsv"
        tsv_file.write_text("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tT\n")
        xlsx_file = tmp_path / "results.xlsx"
        convert_to_excel(str(tsv_file), cfg={}, output_file=str(xlsx_file))
        orig_bytes = xlsx_file.read_bytes()

        with (
            patch("pandas.DataFrame.to_excel", side_effect=ValueError("Workbook failure")),
            pytest.raises(ValueError, match="Workbook failure"),
        ):
            convert_to_excel(str(tsv_file), cfg={}, output_file=str(xlsx_file))

        assert xlsx_file.read_bytes() == orig_bytes
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_excel_sheet_addition_failure_preserves_original(self, tmp_path: Path) -> None:
        tsv_file = tmp_path / "input.tsv"
        tsv_file.write_text("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tT\n")
        xlsx_file = tmp_path / "results.xlsx"
        convert_to_excel(str(tsv_file), cfg={}, output_file=str(xlsx_file))
        orig_bytes = xlsx_file.read_bytes()

        with (
            patch("pandas.DataFrame.to_excel", side_effect=ValueError("Sheet write failed")),
            pytest.raises(ValueError, match="Sheet write failed"),
        ):
            append_tsv_as_sheet(str(xlsx_file), str(tsv_file), sheet_name="Metadata")

        assert xlsx_file.read_bytes() == orig_bytes
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_html_report_atomic_failure_leaves_target_intact(self, tmp_path: Path) -> None:
        html_file = tmp_path / "index.html"
        html_file.write_text("<html>original</html>")

        variants_json = tmp_path / "variants.json"
        variants_json.write_text("[]")
        summary_json = tmp_path / "summary.json"
        summary_json.write_text("{}")

        with (
            patch("jinja2.Environment.get_template", side_effect=RuntimeError("Template error")),
            pytest.raises(RuntimeError, match="Template error"),
        ):
            generate_html_report(str(variants_json), str(summary_json), str(tmp_path), cfg={})

        assert html_file.read_text() == "<html>original</html>"
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_analysis_stages_atomic_failure_preserves_original(self, tmp_path: Path) -> None:
        burden_out = tmp_path / "burden.tsv"
        burden_out.write_text("gene\tp_val\nBRCA1\t0.01\n")

        burden_stage = GeneBurdenAnalysisStage()
        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={
                "perform_gene_burden": True,
                "gene_burden_output": str(burden_out),
                "vcf_samples": ["s1", "s2"],
                "case_samples": ["s1"],
                "control_samples": ["s2"],
            },
            workspace=Workspace(tmp_path, "test"),
        )
        ctx.current_dataframe = pd.DataFrame()

        with (
            patch.object(pd.DataFrame, "to_csv", side_effect=OSError("CSV write error")),
            pytest.raises(IOError, match="CSV write error"),
        ):
            burden_stage._process(ctx)

        assert burden_out.read_text() == "gene\tp_val\nBRCA1\t0.01\n"
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

        assoc_out = tmp_path / "assoc.tsv"
        assoc_out.write_text("gene\tmethod\tp_val\nPKD1\tfisher\t0.05\n")

        assoc_stage = AssociationAnalysisStage()
        ctx_assoc = PipelineContext(
            args=argparse.Namespace(),
            config={
                "perform_association": True,
                "gzip_intermediates": False,
                "association_output": str(assoc_out),
                "association_tests": ["fisher"],
                "case_samples": ["s1"],
                "control_samples": ["s2"],
            },
            workspace=Workspace(tmp_path, "test"),
        )
        ctx_assoc.current_dataframe = pd.DataFrame()

        with (
            patch.object(pd.DataFrame, "to_csv", side_effect=OSError("Assoc write error")),
            pytest.raises(IOError, match="Assoc write error"),
        ):
            assoc_stage._process(ctx_assoc)

        assert assoc_out.read_text() == "gene\tmethod\tp_val\nPKD1\tfisher\t0.05\n"
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_successful_outputs_plain_and_gzip_tsv(self, tmp_path: Path) -> None:
        plain_tsv = tmp_path / "plain.tsv"
        with AtomicFileOperation(plain_tsv) as tmp:
            pd.DataFrame([{"A": 1, "B": 2}]).to_csv(tmp, sep="\t", index=False)

        assert plain_tsv.exists()
        loaded = pd.read_csv(plain_tsv, sep="\t")
        assert len(loaded) == 1
        assert loaded.iloc[0]["A"] == 1

        gzip_tsv = tmp_path / "gzip.tsv.gz"
        with AtomicFileOperation(gzip_tsv) as tmp:
            pd.DataFrame([{"A": 3, "B": 4}]).to_csv(tmp, sep="\t", index=False, compression="gzip")

        assert gzip_tsv.exists()
        loaded_gz = pd.read_csv(gzip_tsv, sep="\t", compression="gzip")
        assert len(loaded_gz) == 1
        assert loaded_gz.iloc[0]["A"] == 3

    def test_run_summary_json_schema_and_stored_sha256(self, tmp_path: Path) -> None:
        vcf = tmp_path / "input.vcf.gz"
        with gzip.open(vcf, "wb") as f:
            f.write(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        out_tsv = tmp_path / "final.tsv"
        out_tsv.write_text("CHROM\tPOS\nchr1\t100\n")

        out_gz = tmp_path / "final.tsv.gz"
        with gzip.open(out_gz, "wb") as f:
            f.write(b"CHROM\tPOS\nchr1\t200\n")

        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={"vcf_file": str(vcf)},
            workspace=Workspace(tmp_path, "test"),
        )
        ctx.vcf_samples = ["sample1", "sample2"]
        ctx.final_output_path = out_tsv
        ctx.report_paths["compressed_tsv"] = out_gz
        ctx.current_dataframe = pd.DataFrame([{"CHROM": "chr1", "POS": 100, "GENE": "PKD1"}])
        ctx.mark_complete("vcf_setup")
        ctx.mark_complete("dataframe_loading")
        ctx.mark_complete("tsv_output")

        summary_file = write_run_summary(ctx)
        assert summary_file.is_file()

        with open(summary_file, encoding="utf-8") as f:
            summary = json.load(f)

        assert summary["schema_version"] == SUMMARY_SCHEMA_VERSION
        assert summary["execution_status"] == "success"
        assert "variantcentrifuge" in summary["versions"]
        assert "python" in summary["versions"]
        assert len(summary["stages"]) == 3
        assert summary["tier_counts"]["final_variants"] == 1
        assert summary["tier_counts"]["unique_genes"] == 1

        for category in ("inputs", "outputs"):
            for item in summary["files"][category]:
                filepath = Path(item["path"])
                assert filepath.is_file()
                expected_sha = hashlib.sha256(filepath.read_bytes()).hexdigest()
                assert item["sha256"] == expected_sha
                assert item["size_bytes"] == os.path.getsize(filepath)

    def test_run_summary_works_without_checkpoint_mode(self, tmp_path: Path) -> None:
        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={"vcf_file": "sample.vcf"},
            workspace=Workspace(tmp_path, "test"),
        )
        summary_file = write_run_summary(ctx)
        assert summary_file.is_file()
        with open(summary_file, encoding="utf-8") as f:
            summary = json.load(f)
        assert summary["schema_version"] == SUMMARY_SCHEMA_VERSION
        assert summary["execution_status"] == "success"
