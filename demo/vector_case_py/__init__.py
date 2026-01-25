# -*- coding: utf-8 -*-
"""
VECTOR-CASE Python Interface
============================

A Pythonic interface for the VECTOR-CASE mutual inductance calculation library.

This package provides an object-oriented interface to the underlying C library
for calculating mutual inductance between coaxial coils. The library is
automatically compiled using CFFI if not already available.

Quick Start
-----------
Basic usage::

    from vector_case_py import Coil, SumPrecision, MutualInductanceCalculator

    # Define two coils
    coil1 = Coil(inner_radius=0.1, outer_radius=0.2, length=0.1, num_turns=100)
    coil2 = Coil(inner_radius=0.3, outer_radius=0.4, length=0.1, num_turns=100)

    # Create calculator and compute inductance
    calculator = MutualInductanceCalculator()
    inductance = calculator.calculate(coil1, coil2, distance=0.1)

Modules
-------
structs
    Data classes for coils and precision parameters.
platform_utils
    Platform detection and compiler flag generation.
compiler
    CFFI library compilation utilities.
calculator
    Main mutual inductance calculator.

Notes
-----
Supported platforms: Linux, macOS (Intel and Apple Silicon), Windows.
The module automatically detects CPU capabilities and applies appropriate
optimization flags (SSE, AVX2, AVX512, NEON).

Author
------
Davor Dobrota
"""

__version__ = "1.0.0"
__author__ = "Davor Dobrota"

# Core classes - main public API
from .structs import Coil, SumPrecision
from .calculator import MutualInductanceCalculator

# Platform utilities
from .platform_utils import (
    SIMDCapability,
    detect_simd_capability,
    get_compile_flags,
    get_platform_info,
)

# Compilation utilities
from .compiler import (
    compile_library,
    load_library,
    get_library,
    is_library_compiled,
)

# Define public API
__all__ = [
    # Version info
    "__version__",
    "__author__",
    # Core classes
    "Coil",
    "SumPrecision",
    "MutualInductanceCalculator",
    # Platform utilities
    "SIMDCapability",
    "detect_simd_capability",
    "get_compile_flags",
    "get_platform_info",
    # Compilation utilities
    "compile_library",
    "load_library",
    "get_library",
    "is_library_compiled",
]
