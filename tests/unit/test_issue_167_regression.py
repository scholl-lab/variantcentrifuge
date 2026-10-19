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
        def fake_to_csv(path_or_buf, *args, **kwargs):
            if path_or_buf is not None:
                Path(path_or_buf).write_text("partial_corrupt_data")
            raise OSError("Write failed mid-stream")

        # 1. Existing destination remains intact
        output_file = tmp_path / "output.tsv"
        output_file.write_text("existing header\n1\t100\tA\tT\n")

        stage = TSVOutputStage()
        context = PipelineContext(
            args=argparse.Namespace(),
            config={"output_file": str(output_file)},
            workspace=Workspace(tmp_path, "test"),
        )
        context.current_dataframe = pd.DataFrame([{"CHROM": "chr1", "POS": 100}])

        with (
            patch.object(pd.DataFrame, "to_csv", side_effect=fake_to_csv),
            pytest.raises(IOError, match="Write failed mid-stream"),
        ):
            stage._process(context)

        assert output_file.read_text() == "existing header\n1\t100\tA\tT\n"
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

        # 2. Non-existent destination remains absent
        output_new = tmp_path / "output_new.tsv"
        context.config["output_file"] = str(output_new)
        with (
            patch.object(pd.DataFrame, "to_csv", side_effect=fake_to_csv),
            pytest.raises(IOError, match="Write failed mid-stream"),
        ):
            stage._process(context)

        assert not output_new.exists()
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

        # 3. Gzipped TSV non-existent destination remains absent
        output_gz = tmp_path / "output_new.tsv.gz"
        context.config["output_file"] = str(output_gz)
        with (
            patch.object(pd.DataFrame, "to_csv", side_effect=fake_to_csv),
            pytest.raises(IOError, match="Write failed mid-stream"),
        ):
            stage._process(context)

        assert not output_gz.exists()
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_chunk_copy_atomic_failure_preserves_destination(self, tmp_path: Path) -> None:
        """TSVOutputStage._process chunk copy failure with partial bytes leaves existing file
        intact and absent file absent, without leaving any .tmp. files.
        """
        from variantcentrifuge.stages.tsv_output_stages import TSVOutputStage

        chunk_tsv = tmp_path / "chunk_results.tsv"
        chunk_tsv.write_text("CHROM\tPOS\nchr1\t100\n")

        dst_existing = tmp_path / "final_chunk_existing.tsv"
        dst_existing.write_text("original final data")

        def fake_copy2(src, dst, *args, **kwargs):
            Path(dst).write_bytes(b"corrupt partial chunk bytes")
            raise OSError("Copy aborted mid-stream")

        # 1. Destination exists: remains intact
        ws1 = Workspace(tmp_path / "ws_chunk1", "test")
        ctx1 = PipelineContext(
            args=argparse.Namespace(),
            config={"output_file": str(dst_existing), "use_chunked_processing": True},
            workspace=ws1,
        )
        ctx1.chunked_analysis_tsv = chunk_tsv

        with (
            patch("shutil.copy2", side_effect=fake_copy2),
            pytest.raises(OSError, match="Copy aborted mid-stream"),
        ):
            TSVOutputStage()._process(ctx1)

        assert dst_existing.read_text() == "original final data"
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

        # 2. Destination did not exist: remains absent
        dst_new = tmp_path / "new_chunk_absent.tsv"
        ws2 = Workspace(tmp_path / "ws_chunk2", "test")
        ctx2 = PipelineContext(
            args=argparse.Namespace(),
            config={"output_file": str(dst_new), "use_chunked_processing": True},
            workspace=ws2,
        )
        ctx2.chunked_analysis_tsv = chunk_tsv

        with (
            patch("shutil.copy2", side_effect=fake_copy2),
            pytest.raises(OSError, match="Copy aborted mid-stream"),
        ):
            TSVOutputStage()._process(ctx2)

        assert not dst_new.exists()
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_excel_output_atomic_failure_preserves_original(self, tmp_path: Path) -> None:
        tsv_file = tmp_path / "input.tsv"
        tsv_file.write_text("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tT\n")
        xlsx_file = tmp_path / "results.xlsx"
        convert_to_excel(str(tsv_file), cfg={}, output_file=str(xlsx_file))
        orig_bytes = xlsx_file.read_bytes()

        def fake_to_excel(excel_writer, *args, **kwargs):
            if hasattr(excel_writer, "handles") and hasattr(excel_writer.handles, "handle"):
                excel_writer.handles.handle.write(b"corrupt partial workbook")
            raise ValueError("Workbook failure mid-stream")

        with (
            patch("pandas.DataFrame.to_excel", side_effect=fake_to_excel),
            pytest.raises(ValueError, match="Workbook failure mid-stream"),
        ):
            convert_to_excel(str(tsv_file), cfg={}, output_file=str(xlsx_file))

        assert xlsx_file.read_bytes() == orig_bytes
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

        # Destination does not exist initially
        xlsx_new = tmp_path / "new_results.xlsx"
        with (
            patch("pandas.DataFrame.to_excel", side_effect=fake_to_excel),
            pytest.raises(ValueError, match="Workbook failure mid-stream"),
        ):
            convert_to_excel(str(tsv_file), cfg={}, output_file=str(xlsx_new))

        assert not xlsx_new.exists()
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_excel_sheet_addition_failure_preserves_original(self, tmp_path: Path) -> None:
        tsv_file = tmp_path / "input.tsv"
        tsv_file.write_text("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tT\n")
        xlsx_file = tmp_path / "results.xlsx"
        convert_to_excel(str(tsv_file), cfg={}, output_file=str(xlsx_file))
        orig_bytes = xlsx_file.read_bytes()

        def fake_sheet_to_excel(excel_writer, *args, **kwargs):
            raise ValueError("Sheet write failed mid-stream")

        with (
            patch("pandas.DataFrame.to_excel", side_effect=fake_sheet_to_excel),
            pytest.raises(ValueError, match="Sheet write failed mid-stream"),
        ):
            append_tsv_as_sheet(str(xlsx_file), str(tsv_file), sheet_name="Metadata")

        assert xlsx_file.read_bytes() == orig_bytes
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_finalize_excel_file_failure_preserves_original(self, tmp_path: Path) -> None:
        """finalize_excel_file failure during wb.save preserves original workbook without .tmp."""
        from variantcentrifuge.converter import finalize_excel_file

        tsv_file = tmp_path / "input_finalize.tsv"
        tsv_file.write_text("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tT\n")
        xlsx_file = tmp_path / "results_finalize.xlsx"
        convert_to_excel(str(tsv_file), cfg={}, output_file=str(xlsx_file))
        orig_bytes = xlsx_file.read_bytes()

        def fake_wb_save(filename, *args, **kwargs):
            Path(filename).write_bytes(b"corrupted partial openpyxl save")
            raise OSError("Disk error during openpyxl wb.save")

        with (
            patch("openpyxl.Workbook.save", side_effect=fake_wb_save),
            pytest.raises(OSError, match=r"Disk error during openpyxl wb\.save"),
        ):
            finalize_excel_file(str(xlsx_file), cfg={})

        assert xlsx_file.read_bytes() == orig_bytes
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_html_report_atomic_failure_leaves_target_intact(self, tmp_path: Path) -> None:
        html_file = tmp_path / "index.html"
        html_file.write_text("<html>original</html>")

        variants_json = tmp_path / "variants.json"
        variants_json.write_text("[]")
        summary_json = tmp_path / "summary.json"
        summary_json.write_text(
            '{"num_variants": 0, "num_genes": 0, "num_samples": 0, "impact_distribution": {}}'
        )

        with (
            patch("os.replace", side_effect=OSError("Replace failed")),
            pytest.raises(OSError, match="Replace failed"),
        ):
            generate_html_report(str(variants_json), str(summary_json), str(tmp_path), cfg={})

        assert html_file.read_text() == "<html>original</html>"
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_html_report_stage_process_failure_when_absent(self, tmp_path: Path) -> None:
        """HTMLReportStage._process failure when index.html does not exist leaves index.html absent
        and leaves no .tmp. files.
        """
        from variantcentrifuge.stages.output_stages import HTMLReportStage

        ws = Workspace(tmp_path / "ws_html", "test")
        report_dir = ws.output_dir / "report"
        report_dir.mkdir(parents=True, exist_ok=True)
        html_file = report_dir / "index.html"
        assert not html_file.exists()

        variants_json = report_dir / "variants.json"
        variants_json.write_text("[]")
        summary_json = report_dir / "summary.json"
        summary_json.write_text("{}")

        input_tsv = ws.output_dir / "final_variants.tsv"
        input_tsv.write_text("CHROM\tPOS\nchr1\t100\n")

        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={"html_report": True},
            workspace=ws,
        )
        ctx.final_output_path = input_tsv
        ctx.current_dataframe = pd.DataFrame({"CHROM": ["chr1"], "POS": [100]})

        def fake_produce_report_json(variant_tsv: str, out_dir: str) -> None:
            variants_json.write_text("[]")
            summary_json.write_text("{}")

        def fake_generate_report(*args, **kwargs):
            # simulate writing partial tmp file
            for f in report_dir.glob("*.tmp*"):
                f.write_text("partial html")
            raise RuntimeError("Report generation failed mid-stream")

        with (
            patch(
                "variantcentrifuge.stages.output_stages.produce_report_json",
                side_effect=fake_produce_report_json,
            ),
            patch(
                "variantcentrifuge.stages.output_stages.generate_html_report",
                side_effect=fake_generate_report,
            ),
            pytest.raises(RuntimeError, match="Report generation failed mid-stream"),
        ):
            HTMLReportStage()._process(ctx)

        assert not html_file.exists()
        assert len([f for f in report_dir.glob("*") if ".tmp." in f.name]) == 0

    def test_analysis_stages_atomic_failure_preserves_original(self, tmp_path: Path) -> None:
        def fake_analysis_to_csv(path_or_buf, *args, **kwargs):
            if path_or_buf is not None:
                Path(path_or_buf).write_text("partial_corrupt_analysis")
            raise OSError("CSV write error")

        # 1. Existing burden output remains intact
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
            patch.object(pd.DataFrame, "to_csv", side_effect=fake_analysis_to_csv),
            pytest.raises(IOError, match="CSV write error"),
        ):
            burden_stage._process(ctx)

        assert burden_out.read_text() == "gene\tp_val\nBRCA1\t0.01\n"
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

        # Non-existent burden output remains absent
        burden_new = tmp_path / "burden_new.tsv"
        ctx.config["gene_burden_output"] = str(burden_new)
        with (
            patch.object(pd.DataFrame, "to_csv", side_effect=fake_analysis_to_csv),
            pytest.raises(IOError, match="CSV write error"),
        ):
            burden_stage._process(ctx)
        assert not burden_new.exists()
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

        # 2. Association output
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
            patch.object(pd.DataFrame, "to_csv", side_effect=fake_analysis_to_csv),
            pytest.raises(IOError, match="CSV write error"),
        ):
            assoc_stage._process(ctx_assoc)

        assert assoc_out.read_text() == "gene\tmethod\tp_val\nPKD1\tfisher\t0.05\n"
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

        # Non-existent association output remains absent
        assoc_new = tmp_path / "assoc_new.tsv"
        ctx_assoc.config["association_output"] = str(assoc_new)
        with (
            patch.object(pd.DataFrame, "to_csv", side_effect=fake_analysis_to_csv),
            pytest.raises(IOError, match="CSV write error"),
        ):
            assoc_stage._process(ctx_assoc)
        assert not assoc_new.exists()
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

    def test_run_summary_atomic_failure_preserves_destination(self, tmp_path: Path) -> None:
        summary_path = tmp_path / "run_summary.json"
        summary_path.write_text('{"original": true}')

        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={"vcf_file": "test.vcf"},
            workspace=Workspace(tmp_path, "test"),
        )

        def fake_dump(obj, fp, *args, **kwargs):
            fp.write("partial summary json")
            raise OSError("Summary write disk full")

        # Destination exists: remains intact
        with (
            patch("json.dump", side_effect=fake_dump),
            pytest.raises(OSError, match="Summary write disk full"),
        ):
            write_run_summary(ctx, output_path=summary_path)

        assert summary_path.read_text() == '{"original": true}'
        assert len([f for f in tmp_path.glob("*") if ".tmp." in f.name]) == 0

        # Destination did not exist: remains absent
        summary_new = tmp_path / "new_run_summary.json"
        with (
            patch("json.dump", side_effect=fake_dump),
            pytest.raises(OSError, match="Summary write disk full"),
        ):
            write_run_summary(ctx, output_path=summary_new)

        assert not summary_new.exists()
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
        ctx.current_dataframe = pd.DataFrame(
            [{"CHROM": "chr1", "POS": 100, "GENE": "PKD1", "IMPACT": "HIGH"}]
        )
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
        assert summary["tier_counts"]["samples"] == 2
        assert summary["tier_counts"]["final_variants"] == 1
        assert summary["tier_counts"]["final_genes"] == 1
        assert summary["tier_counts"]["unique_genes"] == 1
        assert summary["tier_counts"]["extracted_variants"] == 1
        assert summary["tier_counts"]["pre_final_filter_variants"] == 1
        assert summary["tier_counts"]["impact_tiers"]["HIGH"] == 1
        assert summary["tier_counts"]["impact_tiers"]["MODERATE"] == 0
        assert summary["tier_counts"]["impact_tiers"]["LOW"] == 0
        assert summary["tier_counts"]["impact_tiers"]["MODIFIER"] == 0

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

    def test_html_report_output_validation(self, tmp_path: Path) -> None:
        from variantcentrifuge.pipeline_validation import validate_requested_outputs

        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={"html_report": True},
            workspace=Workspace(tmp_path, "test"),
        )
        with pytest.raises(RuntimeError, match="Requested HTML report output was not written"):
            validate_requested_outputs(ctx)

        html_file = tmp_path / "index.html"
        html_file.write_text("<html></html>")
        ctx.report_paths["html"] = html_file
        # Now should pass
        validate_requested_outputs(ctx)

    def test_run_summary_excludes_itself_and_records_durations_and_gz_rows(
        self, tmp_path: Path
    ) -> None:
        import gzip

        out_gz = tmp_path / "final.tsv.gz"
        with gzip.open(out_gz, "wt", encoding="utf-8") as f:
            f.write("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tT\nchr1\t200\tG\tC\n")

        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={"output_file": str(out_gz)},
            workspace=Workspace(tmp_path, "test"),
        )
        ctx.final_output_path = out_gz
        ctx.stage_durations["tsv_output"] = 1.2345
        ctx.mark_complete("tsv_output")

        summary_file = write_run_summary(ctx)
        with open(summary_file, encoding="utf-8") as f:
            summary = json.load(f)

        # Gzip row count should be 2 (3 lines minus header)
        assert summary["tier_counts"]["final_variants"] == 2

        # Stage durations
        assert summary["stages"][0]["name"] == "tsv_output"
        assert summary["stages"][0]["duration_seconds"] == 1.2345

        # Manifest must NOT include run_summary.json itself
        output_paths = [item["path"] for item in summary["files"]["outputs"]]
        assert str(summary_file.resolve()) not in output_paths
        assert summary["execution_status"] == "success"

    def test_collect_metadata_rows(self, tmp_path: Path) -> None:
        from variantcentrifuge.run_summary import _collect_metadata_rows

        tx_file = tmp_path / "tx.txt"
        tx_file.write_text("ENST1\n")
        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={
                "preset": ["somatic_pass"],
                "tumor_sample": "TUMOR",
                "normal_sample": "NORMAL",
                "transcript_file": str(tx_file),
                "git_sha": "abc1234",
            },
            workspace=Workspace(tmp_path, "test"),
        )
        ctx.vcf_samples = ["TUMOR", "NORMAL"]
        rows = _collect_metadata_rows(ctx, lambda tool: "1.0")
        row_dict = dict(rows)
        assert row_dict["Tool"] == "variantcentrifuge"
        assert row_dict["Git SHA"] == "abc1234"
        assert "tumor=TUMOR, normal=NORMAL" in row_dict["Tumor/Normal Samples"]
        assert row_dict["tool.bcftools_version"] == "1.0"

        ctx2 = PipelineContext(
            args=argparse.Namespace(),
            config={"tumor_sample_index": "0", "normal_sample_index": "1"},
            workspace=Workspace(tmp_path, "test2"),
        )
        ctx2.vcf_samples = ["S_TUMOR", "S_NORMAL"]
        rows2 = dict(_collect_metadata_rows(ctx2, lambda tool: "1.0"))
        assert "tumor=S_TUMOR, normal=S_NORMAL" in rows2["Tumor/Normal Samples"]

    def test_pipeline_validation_comprehensive(self, tmp_path: Path) -> None:
        from variantcentrifuge.pipeline_validation import (
            resolve_threads,
            validate_requested_outputs,
        )

        # validate_requested_outputs: TSV output missing
        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={"output_file": "final.tsv"},
            workspace=Workspace(tmp_path, "test"),
        )
        with pytest.raises(RuntimeError, match="Requested TSV output was not written"):
            validate_requested_outputs(ctx)

        # stdout is ignored
        ctx.config["output_file"] = "stdout"
        validate_requested_outputs(ctx)

        # Association output missing
        ctx.config["output_file"] = None
        ctx.config["perform_association"] = True
        ctx.config["association_output"] = str(tmp_path / "missing_assoc.tsv")
        with pytest.raises(RuntimeError, match="Requested association output was not written"):
            validate_requested_outputs(ctx)

        # Association output exists
        assoc_file = tmp_path / "missing_assoc.tsv"
        assoc_file.touch()
        validate_requested_outputs(ctx)

        # Gene burden output missing
        ctx.config["perform_association"] = False
        ctx.config["perform_gene_burden"] = True
        ctx.config["gene_burden_output"] = str(tmp_path / "missing_burden.tsv")
        with pytest.raises(RuntimeError, match="Requested gene burden output was not written"):
            validate_requested_outputs(ctx)

        # Gene burden output exists
        burden_file = tmp_path / "missing_burden.tsv"
        burden_file.touch()
        validate_requested_outputs(ctx)

        # Excel output missing
        ctx.config["perform_gene_burden"] = False
        ctx.config["xlsx"] = True
        with pytest.raises(RuntimeError, match="Requested Excel output was not written"):
            validate_requested_outputs(ctx)

        excel_file = tmp_path / "report.xlsx"
        excel_file.touch()
        ctx.report_paths["excel"] = excel_file
        validate_requested_outputs(ctx)

        # resolve_threads
        assert resolve_threads("4") == 4
        assert resolve_threads(8) == 8
        assert resolve_threads("invalid") == 1
        assert resolve_threads("auto") >= 1

    def test_final_filtering_stage_records_pre_final_filter_variants(self, tmp_path: Path) -> None:
        """FinalFilteringStage records pre_final_filter_variants accurately when filtering."""
        from variantcentrifuge.pipeline_core.context import PipelineContext
        from variantcentrifuge.pipeline_core.workspace import Workspace
        from variantcentrifuge.run_summary import collect_tier_counts
        from variantcentrifuge.stages.tsv_output_stages import FinalFilteringStage

        ws = Workspace(tmp_path / "ws_filter", "test")
        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={"final_filter": "QUAL > 20"},
            workspace=ws,
        )
        # 3 variants, 2 satisfy QUAL > 20
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1", "chr2"],
                "POS": [100, 200, 300],
                "QUAL": [10, 30, 40],
                "GENE": ["G1", "G1", "G2"],
            }
        )
        ctx.current_dataframe = df
        ctx.vcf_samples = ["sample1"]

        stage = FinalFilteringStage()
        stage._process(ctx)

        assert len(ctx.current_dataframe) == 2
        tier_counts = collect_tier_counts(ctx)
        assert tier_counts["pre_final_filter_variants"] == 3
        assert tier_counts["final_variants"] == 2

    def test_file_manifest_includes_stats_output_and_report_summary_json(
        self, tmp_path: Path
    ) -> None:
        """collect_file_manifest includes stats_output_file and report/summary.json."""
        from variantcentrifuge.pipeline_core.context import PipelineContext
        from variantcentrifuge.pipeline_core.workspace import Workspace
        from variantcentrifuge.run_summary import collect_file_manifest, compute_stored_sha256

        ws = Workspace(tmp_path / "ws_manifest", "test")
        report_dir = ws.output_dir / "report"
        report_dir.mkdir(parents=True, exist_ok=True)

        summary_json = report_dir / "summary.json"
        summary_json.write_text('{"summary": "data"}')

        stats_file = ws.output_dir / "stats.tsv"
        stats_file.write_text("metric\tval\ncount\t10\n")

        ctx = PipelineContext(
            args=argparse.Namespace(),
            config={
                "stats_output_file": "stats.tsv",
                "html_report": True,
            },
            workspace=ws,
        )
        ctx.stats_output_tsv = stats_file
        ctx.report_paths["summary_json"] = summary_json

        manifest = collect_file_manifest(ctx)
        output_paths = {item["path"]: item for item in manifest["outputs"]}

        stats_resolved = str(stats_file.resolve())
        assert stats_resolved in output_paths
        assert output_paths[stats_resolved]["sha256"] == compute_stored_sha256(str(stats_file))

        summary_resolved = str(summary_json.resolve())
        assert summary_resolved in output_paths
        assert output_paths[summary_resolved]["sha256"] == compute_stored_sha256(str(summary_json))
