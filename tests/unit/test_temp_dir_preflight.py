"""--temp-dir wiring and free-space preflight (#213)."""

import os
import tempfile
from collections import namedtuple
from pathlib import Path

import pytest

from variantcentrifuge import disk_preflight
from variantcentrifuge.disk_preflight import check_free_space
from variantcentrifuge.filters import _java_tmp_flags
from variantcentrifuge.pipeline_core.workspace import Workspace

Usage = namedtuple("Usage", "total used free")


@pytest.fixture(autouse=True)
def _restore_tempdir(monkeypatch):
    old = tempfile.tempdir
    monkeypatch.delenv("VC_TEMP_DIR", raising=False)
    monkeypatch.setenv("TMPDIR", os.environ.get("TMPDIR", "/tmp"))  # restored on teardown
    yield
    tempfile.tempdir = old


def test_workspace_default_unchanged(tmp_path):
    ws = Workspace(tmp_path / "out", "b")
    assert ws.temp_dir.parent == ws.intermediate_dir
    assert _java_tmp_flags() == []


def test_workspace_temp_base(tmp_path):
    ws = Workspace(tmp_path / "out", "b", tmp_path / "scratch" / "t")
    assert ws.temp_dir.parent == tmp_path / "scratch" / "t"
    assert tempfile.gettempdir() == str(ws.temp_dir)
    assert _java_tmp_flags() == [f"-Djava.io.tmpdir={ws.temp_dir}"]


def _vcf(tmp_path):
    f = tmp_path / "x.vcf.gz"
    f.write_bytes(b"x" * 1000)
    return f


def test_preflight_fails_when_below_input_size(tmp_path, monkeypatch):
    monkeypatch.setattr(disk_preflight.shutil, "disk_usage", lambda d: Usage(0, 0, 500))
    with pytest.raises(OSError, match="--temp-dir"):
        check_free_space(_vcf(tmp_path), [tmp_path])


def test_preflight_warns_then_passes(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(disk_preflight.shutil, "disk_usage", lambda d: Usage(0, 0, 5000))
    check_free_space(_vcf(tmp_path), [tmp_path, Path(tmp_path)])
    assert len(caplog.records) == 1  # duplicate dirs checked once
    caplog.clear()
    monkeypatch.setattr(disk_preflight.shutil, "disk_usage", lambda d: Usage(0, 0, 10**9))
    check_free_space(_vcf(tmp_path), [tmp_path])
    assert not caplog.records


def test_preflight_skips_missing_input(tmp_path):
    check_free_space(tmp_path / "nope.vcf.gz", [tmp_path])


def test_cleanup_restores_process_tempdir(tmp_path):
    before = (tempfile.tempdir, os.environ.get("TMPDIR"), os.environ.get("VC_TEMP_DIR"))
    ws = Workspace(tmp_path / "out", "b", tmp_path / "scratch")
    assert os.environ["VC_TEMP_DIR"] == str(ws.temp_dir)
    ws.cleanup()
    assert (tempfile.tempdir, os.environ.get("TMPDIR"), os.environ.get("VC_TEMP_DIR")) == before
    os.close(tempfile.mkstemp()[0])  # would raise if tempdir pointed at the deleted dir
