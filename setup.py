"""
Minimal setup.py for CFFI C extension support.

All package metadata is in pyproject.toml. This file exists solely to pass
``cffi_modules`` to setuptools (a setup() parameter not supported by
pyproject.toml's [tool.setuptools] table) and to install OptionalBuildExt
so that C extension build failures are non-fatal.

The Davies qfc() C++ extension is optional: when the build fails, SKAT
will use Kuonen saddlepoint or Liu moment-matching as the p-value fallback.
"""

import os

from setuptools import setup
from setuptools.command.build_ext import build_ext


class OptionalBuildExt(build_ext):
    """Allow the Davies C extension build to fail gracefully.

    If the C++ compiler is unavailable or compilation fails for any reason,
    the package installs successfully and SKAT uses the Liu moment-matching
    fallback for p-value computation.
    """

    def build_extension(self, ext):
        try:
            super().build_extension(ext)
        except Exception as e:
            import warnings

            warnings.warn(
                f"Davies C extension build failed: {e}. "
                "SKAT will use Kuonen saddlepoint / Liu moment-matching fallback "
                "for p-values (no accuracy loss for p > 1e-5).",
                stacklevel=2,
            )


# cffi_modules only processed when cffi is available and C extensions are enabled
cffi_modules = []
if not os.environ.get("VARIANTCENTRIFUGE_NO_C_EXT"):
    try:
        import cffi  # noqa: F401

        cffi_modules = ["variantcentrifuge/_davies_build.py:ffibuilder"]
    except ImportError:
        pass

setup(
    cffi_modules=cffi_modules,
    cmdclass={"build_ext": OptionalBuildExt},
)
