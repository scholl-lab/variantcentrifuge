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
    # Forward declaration of qfc() with explicit C linkage so the CFFI-generated
    # C++ wrapper can call the C-exported symbol from qfc.cpp.
    # The extern "C" is required because source_extension='.cpp' causes CFFI to
    # compile the wrapper as C++, which would otherwise produce C++ name-mangled
    # symbol references incompatible with the C-linked export in qfc.cpp.
    """
    #ifdef __cplusplus
    extern "C" {
    #endif

    void qfc(double* lb1, double* nc1, int* n1, int *r1, double *sigma,
             double *c1, int *lim1, double *acc, double* trace,
             int* ifault, double *res);

    #ifdef __cplusplus
    }
    #endif
    """,
    sources=[_src],
    source_extension=".cpp",
    libraries=[],
    extra_compile_args=["-std=c++11", "-O2"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
