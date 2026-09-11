"""Regression coverage for editable installs without a working C++ compiler."""

import runpy
from pathlib import Path

import pytest
import setuptools
from setuptools import Distribution, Extension
from setuptools.command.build_ext import build_ext


@pytest.fixture
def optional_build(tmp_path, monkeypatch):
    """Load the real build command without invoking package installation."""
    monkeypatch.setattr(setuptools, "setup", lambda **kwargs: None)
    namespace = runpy.run_path(str(Path(__file__).parents[2] / "setup.py"))
    monkeypatch.chdir(tmp_path)
    package = tmp_path / "variantcentrifuge"
    package.mkdir()
    extension = Extension("variantcentrifuge._qfc", sources=["qfc.cpp"])
    distribution = Distribution({"packages": ["variantcentrifuge"], "ext_modules": [extension]})
    command = namespace["OptionalBuildExt"](distribution)
    command.build_lib = str(tmp_path / "build")
    command.inplace = True
    command.ensure_finalized()
    return command, extension, package


def test_failed_optional_build_does_not_copy_missing_binary(optional_build, monkeypatch):
    """An unavailable compiler must not break the editable copy phase."""
    command, extension, package = optional_build

    def fail_compile(self, ext):
        raise RuntimeError("C++ compiler unavailable")

    monkeypatch.setattr(build_ext, "build_extension", fail_compile)
    with pytest.warns(UserWarning, match="Davies C extension build failed"):
        command.build_extension(extension)

    command.copy_extensions_to_source()
    assert list(package.iterdir()) == []


def test_successful_optional_build_copies_binary(optional_build, monkeypatch):
    """The fallback must preserve installation of successfully built binaries."""
    command, extension, package = optional_build
    filename = Path(command.get_ext_filename(extension.name))

    def compile_binary(self, ext):
        artifact = Path(self.build_lib) / filename
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_bytes(b"compiled extension")

    monkeypatch.setattr(build_ext, "build_extension", compile_binary)
    command.build_extension(extension)
    command.copy_extensions_to_source()
    assert (package / filename.name).read_bytes() == b"compiled extension"
