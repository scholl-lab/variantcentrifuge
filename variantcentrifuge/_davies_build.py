"""
CFFI out-of-line API build script for the Davies qfc C++ extension.

Usage (for direct compilation):
    python variantcentrifuge/_davies_build.py

The compiled module is installed automatically when the package is built
via setuptools + cffi_modules (see setup.py).
"""

import os

from cffi import FFI

ffibuilder = FFI()

# C declaration for the qfc() function exposed via extern "C" in qfc.cpp.
# Signature from Davies (1980) / CompQuadForm package src/qfc.cpp.
ffibuilder.cdef(
    """
    void qfc(double* lb1, double* nc1, int* n1, int *r1, double *sigma,
             double *c1, int *lim1, double *acc, double* trace,
             int* ifault, double *res);
"""
)

_src = os.path.join(
    os.path.dirname(__file__),
    "association",
    "data",
    "qfc.cpp",
)

ffibuilder.set_source(
    "variantcentrifuge._qfc",
    # Empty C header: the function is declared via cdef above.
    # The actual implementation is in qfc.cpp (C++).
    "",
    sources=[_src],
    source_extension=".cpp",
    libraries=[],
    extra_compile_args=["-std=c++11", "-O2"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
