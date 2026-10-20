"""Regression tests for Issue #163:
- Module-level logging.getLogger(__name__) in runner
- Guard psutil gracefully for memory tracking in runner (import failure and runtime errors)
"""

from __future__ import annotations

import logging
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

import variantcentrifuge.pipeline_core.runner as runner_module
from variantcentrifuge.pipeline_core import PipelineContext, PipelineRunner, Stage


class SimpleTestStage(Stage):
    """Simple stage for runner execution tests."""

    @property
    def name(self) -> str:
        return "simple_stage"

    @property
    def dependencies(self) -> set[str]:
        return set()

    def _process(self, context: PipelineContext) -> PipelineContext:
        context.mark_complete(self.name)
        return context


@pytest.mark.unit
class TestRunnerLogger:
    """Verify module-level logger conventions in runner.py."""

    def test_runner_module_logger(self) -> None:
        """runner.py must have a module-level logger named after its module (__name__)."""
        assert hasattr(runner_module, "logger")
        assert isinstance(runner_module.logger, logging.Logger)
        assert runner_module.logger.name == "variantcentrifuge.pipeline_core.runner"


@pytest.mark.unit
class TestRunnerPsutilGracefulGuarding:
    """Verify psutil is guarded gracefully when unavailable or failing."""

    def test_runner_executes_when_psutil_none(self) -> None:
        """PipelineRunner should complete execution without error when psutil is None."""
        runner = PipelineRunner(enable_checkpoints=False, max_workers=1)
        stage = SimpleTestStage()
        context = PipelineContext(args=Mock(), config={}, workspace=Mock())

        with patch.object(runner_module, "psutil", None):
            result_context = runner.run([stage], context)

        assert result_context.is_complete("simple_stage")
        # Memory metrics should either be empty or contain graceful defaults
        assert stage.name not in runner._stage_metrics or runner._stage_metrics[stage.name] == {}

    def test_runner_handles_psutil_runtime_exception(self) -> None:
        """PipelineRunner should gracefully catch runtime exceptions during memory measurement."""
        runner = PipelineRunner(enable_checkpoints=False, max_workers=1)
        stage = SimpleTestStage()
        context = PipelineContext(args=Mock(), config={}, workspace=Mock())

        mock_psutil = Mock()
        mock_psutil.Process.side_effect = RuntimeError("Process memory info unavailable")

        with patch.object(runner_module, "psutil", mock_psutil):
            result_context = runner.run([stage], context)

        assert result_context.is_complete("simple_stage")

    def test_runner_tracks_memory_when_psutil_available(self) -> None:
        """PipelineRunner tracks memory when psutil works normally."""
        runner = PipelineRunner(enable_checkpoints=False, max_workers=1)
        stage = SimpleTestStage()
        context = PipelineContext(args=Mock(), config={}, workspace=Mock())

        mock_process = Mock()
        mock_process.memory_info.return_value = Mock(rss=100 * 1024 * 1024)  # 100 MB
        mock_psutil = Mock()
        mock_psutil.Process.return_value = mock_process

        with patch.object(runner_module, "psutil", mock_psutil):
            result_context = runner.run([stage], context)

        assert result_context.is_complete("simple_stage")
        assert stage.name in runner._stage_metrics
        assert runner._stage_metrics[stage.name]["mem_before_mb"] == pytest.approx(100.0)
        assert runner._stage_metrics[stage.name]["mem_after_mb"] == pytest.approx(100.0)


@pytest.mark.unit
class TestDocCliSynchronization:
    """Verify synchronization between CLI options and documentation."""

    def test_documented_flags_exist_in_usage_md(self) -> None:
        """usage.md must document --coast-classification, --filter-exception-file, and --pca."""
        usage_text = Path("docs/source/usage.md").read_text(encoding="utf-8")
        assert "--coast-classification" in usage_text
        assert "--filter-exception-file" in usage_text
        assert "`--pca`" in usage_text or "--pca PATH" in usage_text

    def test_coast_backend_documentation_omits_r(self) -> None:
        """COAST backend documentation must not recommend or mention the unsupported 'r' backend."""
        docs_dir = Path("docs/source")
        for doc_file in docs_dir.rglob("*.md"):
            content = doc_file.read_text(encoding="utf-8")
            assert "--coast-backend r" not in content, (
                f"{doc_file} contains deprecated/unsupported '--coast-backend r'"
            )

    def test_pca_sites_documentation_omits_bed(self) -> None:
        """--pca-sites documentation must specify sites VCF/BCF, not BED."""
        docs_dir = Path("docs/source")
        for doc_file in docs_dir.rglob("*.md"):
            content = doc_file.read_text(encoding="utf-8")
            assert "--pca-sites sites.bed" not in content, (
                f"{doc_file} contains incorrect '--pca-sites sites.bed'"
            )

    def test_parse_representative_cli_argument_lists(self) -> None:
        """Registered parser successfully parses representative documented argument lists."""
        from variantcentrifuge.cli import create_parser

        parser = create_parser()

        # 1. COAST classification options
        args1 = parser.parse_args(["-v", "test.vcf.gz", "--coast-classification", "sift_polyphen"])
        assert args1.coast_classification == "sift_polyphen"

        args2 = parser.parse_args(["-v", "test.vcf.gz", "--coast-classification", "cadd"])
        assert args2.coast_classification == "cadd"

        # 2. Filter exception file
        args3 = parser.parse_args(["-v", "test.vcf.gz", "--filter-exception-file", "hotspots.tsv"])
        assert args3.filter_exception_file == "hotspots.tsv"

        # 3. PCA unified flag and backward compatibility aliases
        args4 = parser.parse_args(["-v", "test.vcf.gz", "--pca", "eigenvec.tsv"])
        assert args4.pca == "eigenvec.tsv"

        args5 = parser.parse_args(["-v", "test.vcf.gz", "--pca", "akt"])
        assert args5.pca == "akt"

        args6 = parser.parse_args(["-v", "test.vcf.gz", "--pca-file", "eigenvec.tsv"])
        assert args6.pca == "eigenvec.tsv"

        args7 = parser.parse_args(["-v", "test.vcf.gz", "--pca-tool", "akt"])
        assert args7.pca == "akt"

        # 4. COAST backend valid choices
        args_py = parser.parse_args(["-v", "test.vcf.gz", "--coast-backend", "python"])
        assert args_py.coast_backend == "python"

        args_auto = parser.parse_args(["-v", "test.vcf.gz", "--coast-backend", "auto"])
        assert args_auto.coast_backend == "auto"

        # 5. Invalid COAST backend 'r' rejected
        with pytest.raises(SystemExit) as exc:
            parser.parse_args(["-v", "test.vcf.gz", "--coast-backend", "r"])
        assert exc.value.code == 2
