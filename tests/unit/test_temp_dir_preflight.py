"""--temp-dir wiring and free-space preflight (#213)."""

import contextlib
import os
import tempfile
from collections import namedtuple
from pathlib import Path
from unittest.mock import patch

from variantcentrifuge import disk_preflight
from variantcentrifuge.disk_preflight import check_free_space
from variantcentrifuge.filters import _java_tmp_flags
from variantcentrifuge.pipeline_core.workspace import Workspace

Usage = namedtuple("Usage", "total used free")


def test_workspace_default_unchanged(tmp_path):
    ws = Workspace(tmp_path / "out", "b")
    assert ws.temp_dir.parent == ws.intermediate_dir
    assert _java_tmp_flags({}) == []


def test_workspace_temp_base_leaves_process_tempdir_alone(tmp_path):
    before = (tempfile.gettempdir(), os.environ.get("TMPDIR"))
    base = tmp_path / ("long_scratch_dir_" * 6) / "t"  # >58 chars would break AF_UNIX if global
    ws = Workspace(tmp_path / "out", "b", base)
    assert ws.temp_dir.parent == base
    assert (tempfile.gettempdir(), os.environ.get("TMPDIR")) == before
    ws.cleanup()
    os.close(tempfile.mkstemp()[0])


def test_java_flag_follows_cfg():
    assert _java_tmp_flags({"temp_dir": "/x"}) == ["-Djava.io.tmpdir=/x"]


def test_filter_and_extract_mkstemp_use_cfg_temp_dir(tmp_path):
    from variantcentrifuge import filters

    cfg = {"temp_dir": str(tmp_path)}
    with (
        patch.object(filters, "validate_bed_vcf_contigs"),
        patch.object(filters.tempfile, "mkstemp", side_effect=RuntimeError("stop")) as mk,
        contextlib.suppress(RuntimeError),
    ):
        filters.extract_variants("v.vcf.gz", "b.bed", cfg, "o.vcf.gz")
    assert mk.call_args.kwargs["dir"] == str(tmp_path)


def _vcf(tmp_path):
    f = tmp_path / "x.vcf.gz"
    f.write_bytes(b"x" * 1000)
    return f


def test_preflight_warns_never_fails(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(disk_preflight.shutil, "disk_usage", lambda d: Usage(0, 0, 500))
    check_free_space(_vcf(tmp_path), [tmp_path, Path(tmp_path)])
    assert len(caplog.records) == 1  # below even 1x: still only a warning; duplicates once
    caplog.clear()
    monkeypatch.setattr(disk_preflight.shutil, "disk_usage", lambda d: Usage(0, 0, 10**9))
    check_free_space(_vcf(tmp_path), [tmp_path])
    assert not caplog.records


def test_preflight_skips_missing_input(tmp_path):
    check_free_space(tmp_path / "nope.vcf.gz", [tmp_path])
