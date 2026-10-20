"""--perform-gene-burden is deprecated (#217)."""

import argparse
import contextlib
import logging
from unittest.mock import MagicMock

import pytest

from variantcentrifuge.cli_association_groups import add_association_argument_group
from variantcentrifuge.stages.gene_burden_stages import GeneBurdenAnalysisStage


def _run(config, caplog):
    ctx = MagicMock()
    ctx.config = config
    with (
        caplog.at_level(logging.WARNING, logger="variantcentrifuge"),
        contextlib.suppress(ValueError),  # missing cohorts: stops right after the warning
    ):
        GeneBurdenAnalysisStage()._process(ctx)
    return [r.getMessage() for r in caplog.records if "deprecated" in r.getMessage()]


def test_warns_when_requested(caplog):
    cfg = {"perform_gene_burden": True}
    with pytest.warns(DeprecationWarning, match="perform-association"):
        msgs = _run(cfg, caplog)
    assert len(msgs) == 1 and "--missing-as-ref" not in msgs[0]


def test_mentions_missing_as_ref(caplog):
    cfg = {"perform_gene_burden": True, "missing_as_ref": True}
    with pytest.warns(DeprecationWarning):
        msgs = _run(cfg, caplog)
    assert "--missing-as-ref is ignored" in msgs[0]


def test_silent_when_not_requested(caplog):
    assert _run({"perform_gene_burden": False}, caplog) == []


def test_help_says_deprecated():
    p = argparse.ArgumentParser()
    add_association_argument_group(p)
    assert "DEPRECATED" in p.format_help()
