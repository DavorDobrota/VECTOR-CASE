# -*- coding: utf-8 -*-
"""
Platform detection and compiler flag generation for VECTOR-CASE.

This module handles cross-platform detection of SIMD capabilities and
generates appropriate compiler flags for each platform.

Supported Platforms
-------------------
- **Linux**: GCC/Clang with detection via ``/proc/cpuinfo``
- **macOS**: Clang with detection via ``sysctl`` (Intel) or assumes NEON (Apple Silicon)
- **Windows**: MSVC with limited detection (defaults to AVX2 on 64-bit)
"""

from __future__ import annotations

import platform
import subprocess
from enum import Enum, auto
from typing import Optional


class SIMDCapability(Enum):
    """
    Enumeration of supported SIMD instruction sets.

    Attributes
    ----------
    NONE : auto
        No SIMD acceleration.
    SSE : auto
        SSE instruction set (x86/x64).
    AVX2 : auto
        AVX2 instruction set (x86/x64).
    AVX512 : auto
        AVX512 instruction set (x86/x64).
    NEON : auto
        ARM NEON instruction set (Apple Silicon, ARM64).
    """

    NONE = auto()
    SSE = auto()
    AVX2 = auto()
    AVX512 = auto()
    NEON = auto()


def detect_simd_capability() -> SIMDCapability:
    """
    Detect the highest SIMD capability available on the current CPU.

    Returns
    -------
    SIMDCapability
        The highest supported SIMD instruction set.

    Notes
    -----
    Detection methods vary by platform:

    - **Linux**: Reads ``/proc/cpuinfo``
    - **macOS (Intel)**: Uses ``sysctl`` to query CPU features
    - **macOS (Apple Silicon)**: Returns NEON (always available on ARM64)
    - **Windows**: Limited detection; assumes AVX2 on 64-bit systems

    Examples
    --------
    >>> capability = detect_simd_capability()
    >>> capability in SIMDCapability
    True
    """
    system = platform.system()

    if system == "Linux":
        return _detect_simd_linux()
    elif system == "Darwin":
        return _detect_simd_macos()
    elif system == "Windows":
        return _detect_simd_windows()
    else:
        return SIMDCapability.NONE


def _detect_simd_linux() -> SIMDCapability:
    """
    Detect SIMD capability on Linux systems.

    Reads ``/proc/cpuinfo`` to determine available instruction sets.

    Returns
    -------
    SIMDCapability
        The highest detected SIMD capability.
    """
    try:
        with open("/proc/cpuinfo", "r") as f:
            cpuinfo = f.read().lower()

        if "avx512" in cpuinfo:
            return SIMDCapability.AVX512
        elif "avx2" in cpuinfo:
            return SIMDCapability.AVX2
        elif "sse" in cpuinfo:
            return SIMDCapability.SSE
    except (IOError, OSError):
        pass
    return SIMDCapability.NONE


def _detect_simd_macos() -> SIMDCapability:
    """
    Detect SIMD capability on macOS systems.

    For Intel Macs, uses ``sysctl`` to query CPU features.
    For Apple Silicon (ARM64), returns NEON which is always available.

    Returns
    -------
    SIMDCapability
        The highest detected SIMD capability.

    Notes
    -----
    Apple Silicon Macs (M1, M2, etc.) use ARM architecture with NEON
    instructions. The ``-march=native`` flag handles optimization
    automatically on these systems.
    """
    try:
        # Check if running on Apple Silicon (ARM64)
        machine = platform.machine().lower()
        if machine in ("arm64", "aarch64"):
            return SIMDCapability.NEON

        # Intel Mac: Query specific CPU features via sysctl
        # Check AVX512
        try:
            result = subprocess.run(
                ["sysctl", "-n", "hw.optional.avx512f"],
                capture_output=True,
                text=True,
                check=False,
                timeout=5,
            )
            if result.returncode == 0 and result.stdout.strip() == "1":
                return SIMDCapability.AVX512
        except (subprocess.SubprocessError, OSError):
            pass

        # Check AVX2
        try:
            result = subprocess.run(
                ["sysctl", "-n", "hw.optional.avx2_0"],
                capture_output=True,
                text=True,
                check=False,
                timeout=5,
            )
            if result.returncode == 0 and result.stdout.strip() == "1":
                return SIMDCapability.AVX2
        except (subprocess.SubprocessError, OSError):
            pass

        # Fallback: check machdep.cpu.features for older detection method
        try:
            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.features"],
                capture_output=True,
                text=True,
                check=False,
                timeout=5,
            )
            if result.returncode == 0:
                features = result.stdout.lower()
                if "avx2" in features:
                    return SIMDCapability.AVX2
                elif "avx" in features:
                    return SIMDCapability.SSE  # AVX but not AVX2
                elif "sse" in features:
                    return SIMDCapability.SSE
        except (subprocess.SubprocessError, OSError):
            pass

        # SSE is standard on all Intel Macs
        return SIMDCapability.SSE

    except Exception:
        pass

    return SIMDCapability.NONE


def _detect_simd_windows() -> SIMDCapability:
    """
    Detect SIMD capability on Windows systems.

    Currently uses a simplified detection that assumes AVX2 on 64-bit
    systems, as proper CPUID detection requires additional libraries.

    Returns
    -------
    SIMDCapability
        The detected or assumed SIMD capability.

    Notes
    -----
    For more accurate detection on Windows, consider using a library
    like ``cpuinfo`` (``pip install py-cpuinfo``).
    """
    try:
        import struct

        # 64-bit Windows systems generally support at least AVX2
        if struct.calcsize("P") * 8 == 64:
            return SIMDCapability.AVX2
        return SIMDCapability.SSE
    except Exception:
        pass
    return SIMDCapability.SSE


def get_compile_flags(simd: Optional[SIMDCapability] = None) -> list[str]:
    """
    Generate compiler flags for the current platform and SIMD capability.

    Parameters
    ----------
    simd : SIMDCapability, optional
        Target SIMD instruction set. If None, auto-detected.

    Returns
    -------
    list[str]
        List of compiler flags appropriate for the platform.

    Notes
    -----
    Flags are optimized for each platform:

    - **Linux/macOS (GCC/Clang)**: ``-O3``, ``-march=native``, ``-ffast-math``
    - **Windows (MSVC)**: ``/O2``, ``/fp:fast``, ``/arch:AVX2``
    - **Apple Silicon**: Uses ``-march=native`` (NEON auto-enabled)

    Examples
    --------
    >>> flags = get_compile_flags()
    >>> "-O3" in flags or "/O2" in flags
    True
    """
    if simd is None:
        simd = detect_simd_capability()

    system = platform.system()

    if system == "Windows":
        return _get_msvc_flags(simd)
    else:
        return _get_gcc_clang_flags(simd)


def _get_gcc_clang_flags(simd: SIMDCapability) -> list[str]:
    """
    Generate GCC/Clang compiler flags.

    Parameters
    ----------
    simd : SIMDCapability
        Target SIMD instruction set.

    Returns
    -------
    list[str]
        Compiler flags for GCC or Clang.
    """
    flags = ["-std=c11", "-O3", "-ffast-math"]

    # Use native architecture for best performance
    # This works on both Intel and Apple Silicon
    flags.extend(["-march=native", "-mtune=native"])

    # Add explicit SIMD flags as hints (may be redundant with -march=native)
    if simd == SIMDCapability.AVX512:
        flags.append("-mavx512f")
    elif simd == SIMDCapability.AVX2:
        flags.append("-mavx2")
    elif simd == SIMDCapability.SSE:
        flags.append("-msse4.2")
    # NEON is automatically enabled with -march=native on ARM64

    return flags


def _get_msvc_flags(simd: SIMDCapability) -> list[str]:
    """
    Generate MSVC compiler flags.

    Parameters
    ----------
    simd : SIMDCapability
        Target SIMD instruction set.

    Returns
    -------
    list[str]
        Compiler flags for MSVC.
    """
    flags = ["/O2", "/fp:fast"]

    if simd == SIMDCapability.AVX512:
        flags.append("/arch:AVX512")
    elif simd == SIMDCapability.AVX2:
        flags.append("/arch:AVX2")
    elif simd == SIMDCapability.SSE:
        flags.append("/arch:SSE2")

    return flags


def get_platform_info() -> dict[str, str]:
    """
    Get information about the current platform.

    Returns
    -------
    dict[str, str]
        Dictionary containing platform information including system,
        machine architecture, processor, Python version, and detected
        SIMD capability.

    Examples
    --------
    >>> info = get_platform_info()
    >>> "system" in info
    True
    >>> "simd_capability" in info
    True
    """
    return {
        "system": platform.system(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python_version": platform.python_version(),
        "simd_capability": detect_simd_capability().name,
    }
