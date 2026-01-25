# -*- coding: utf-8 -*-
"""
CFFI library compilation for VECTOR-CASE.

This module handles the compilation of the VECTOR-CASE C library using CFFI,
with automatic platform detection and optimization.
"""

from __future__ import annotations

import importlib.util
import os
import platform
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from .platform_utils import detect_simd_capability, get_compile_flags

if TYPE_CHECKING:
    from types import ModuleType

# Module name for the compiled C library
MODULE_NAME = "vector_case"

# CFFI C definitions
CFFI_CDEF = """
typedef struct {
    double r_1;
    double R_1;
    double L_1;
    double N_1;

    double r_2;
    double R_2;
    double L_2;
    double N_2;

} CoilCalculationData;

typedef struct {
    uint32_t k_terms;
    uint32_t l_terms;
    uint32_t n_terms;

} SumPrecisionData;

double calculate_mutual_inductance_far(
    const CoilCalculationData data,
    const SumPrecisionData precision,
    const double d
);

double calculate_mutual_inductance_near(
    const CoilCalculationData data,
    const SumPrecisionData precision,
    const double d,
    const double Z
);

static double calculate_mutual_inductance_radial(
    const CoilCalculationData data,
    const double offset_1,
    const double offset_2,
    const double relative_tol,
    const bool verbose
);

double calculate_mutual_inductance(
    const CoilCalculationData data,
    const SumPrecisionData precision,
    const double d,
    const bool verbose
);

static double calculate_mutual_inductance_raw(
    const double r_1,
    const double R_1,
    const double L_1,
    const double N_1,
    const double r_2,
    const double R_2,
    const double L_2,
    const double N_2,
    const double d,
    const uint32_t k_terms,
    const uint32_t l_terms,
    const uint32_t n_terms,
    const bool verbose
);

static void benchmark_mutual_inductance_far(
    const SumPrecisionData precision,
    const uint32_t n_repeats
);

static void benchmark_mutual_inductance_near(
    const SumPrecisionData precision,
    const uint32_t n_repeats
);

static void benchmark_mutual_inductance_radial(
    const uint32_t n_repeats, 
    const bool fast, 
    const double relative_tol
);
"""


def get_library_paths() -> tuple[Path, Path]:
    """
    Get paths for the C source and output directories.

    Returns
    -------
    tuple[Path, Path]
        Tuple of (source_directory, output_directory).

    Notes
    -----
    The source directory contains the C header files.
    The output directory is where the compiled library will be placed.
    """
    # This file is in demo/vector_case_py/
    package_dir = Path(__file__).parent.resolve()
    demo_dir = package_dir.parent
    src_dir = demo_dir.parent / "src"
    return src_dir, demo_dir


def is_library_compiled() -> bool:
    """
    Check if the VECTOR-CASE library has been compiled.

    Returns
    -------
    bool
        True if the library is already compiled and importable.
    """
    return importlib.util.find_spec(MODULE_NAME) is not None


def compile_library(force: bool = False, verbose: bool = True) -> None:
    """
    Compile the VECTOR-CASE C library using CFFI.

    Parameters
    ----------
    force : bool, optional
        If True, recompile even if the library exists. Default is False.
    verbose : bool, optional
        If True, print compilation progress. Default is True.

    Raises
    ------
    ImportError
        If CFFI is not installed.
    RuntimeError
        If compilation fails.

    Notes
    -----
    The compiled library is placed in the demo directory and can be imported
    as a Python module. The function automatically detects the platform and
    applies appropriate compiler flags for optimal performance.

    Examples
    --------
    >>> compile_library(verbose=True)  # doctest: +SKIP
    Compiling 'vector_case' library...
    Detected SIMD capability: AVX512
    Compilation successful!
    """
    try:
        from cffi import FFI
    except ImportError as e:
        raise ImportError(
            "CFFI is required for compilation. Install with: pip install cffi"
        ) from e

    src_dir, output_dir = get_library_paths()

    # Check if already compiled
    if not force and is_library_compiled():
        if verbose:
            print(f"Module '{MODULE_NAME}' already compiled.")
        return

    if verbose:
        print(f"Compiling '{MODULE_NAME}' library...")
        simd = detect_simd_capability()
        print(f"Detected SIMD capability: {simd.name}")

    ffi = FFI()
    ffi.cdef(CFFI_CDEF)

    # Platform-specific library dependencies
    system = platform.system()
    libraries = ["m"] if system != "Windows" else []

    compile_flags = get_compile_flags()

    if verbose:
        print(f"Using compile flags: {compile_flags}")

    # Use forward slashes for path (works on all platforms)
    include_path = str(src_dir / "inductance.h").replace("\\", "/")

    ffi.set_source(
        MODULE_NAME,
        f'#include "{include_path}"',
        libraries=libraries,
        extra_compile_args=compile_flags,
    )

    # Change to output directory for compilation
    original_dir = os.getcwd()
    try:
        os.chdir(output_dir)
        ffi.compile(verbose=verbose)
    finally:
        os.chdir(original_dir)

    if verbose:
        print("Compilation successful!")


def load_library() -> "ModuleType":
    """
    Load the compiled VECTOR-CASE library.

    This function will compile the library if it hasn't been compiled yet.

    Returns
    -------
    ModuleType
        The loaded vector_case module with ``ffi`` and ``lib`` attributes.

    Raises
    ------
    ImportError
        If the library cannot be loaded after compilation.

    Examples
    --------
    >>> module = load_library()  # doctest: +SKIP
    >>> hasattr(module, 'lib')
    True
    >>> hasattr(module, 'ffi')
    True
    """
    compile_library(force=False, verbose=False)

    spec = importlib.util.find_spec(MODULE_NAME)
    if spec is None:
        raise ImportError(f"Failed to load '{MODULE_NAME}' module after compilation")

    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise ImportError(f"No loader found for '{MODULE_NAME}' module")

    spec.loader.exec_module(module)
    return module


# Cached module reference
_cached_module: Optional["ModuleType"] = None


def get_library() -> "ModuleType":
    """
    Get the VECTOR-CASE library, using cached reference if available.

    This is the recommended way to access the library, as it avoids
    repeated loading overhead.

    Returns
    -------
    ModuleType
        The loaded vector_case module.

    Examples
    --------
    >>> lib = get_library()  # doctest: +SKIP
    >>> lib is get_library()  # Same instance
    True
    """
    global _cached_module
    if _cached_module is None:
        _cached_module = load_library()
    return _cached_module
