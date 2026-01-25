# VECTOR-CASE Python Interface

A Pythonic interface for the VECTOR-CASE mutual inductance calculation library.

## Overview

This package provides an object-oriented interface to the underlying C library for calculating mutual inductance between coaxial coils. The library is automatically compiled using CFFI if not already available.

## Requirements

- Python 3.9+
- CFFI (`pip install cffi`)
- A C compiler (GCC, Clang, or MSVC)

## Installation

1. Install the required Python packages:
Using a virtual environment is recommended.
```bash
python -m venv venv
source venv/bin/activate
pip install cffi
```
or (if the uv package manager is installed)
```bash 
uv sync
```
This handles both installing the required packages and activating the virtual environment.
Optionally, the library can also be installed
```bash
pip install vector-case-py
```

2. The C library will be automatically compiled on first import.

## Quick Start

```python
from vector_case_py import Coil, SumPrecision, MutualInductanceCalculator

# Define two coils
coil1 = Coil(inner_radius=0.1, outer_radius=0.2, length=0.1, num_turns=100)
coil2 = Coil(inner_radius=0.3, outer_radius=0.4, length=0.1, num_turns=100)

# Create calculator and compute inductance
calculator = MutualInductanceCalculator()
inductance = calculator.calculate(coil1, coil2, distance=0.1)
print(f'Mutual inductance: {inductance:.10e} H')
```

## Package Structure

```
vector_case_py/
├── __init__.py       # Main exports and demo function
├── structs.py        # Coil and SumPrecision dataclasses
├── platform_utils.py # SIMD detection and compiler flags
├── calculator.py     # MutualInductanceCalculator class
├── compiler.py       # CFFI library compilation
└── py.typed          # PEP 561 marker for type checking
```

## API Reference

### Classes

#### `Coil`

Represents a coaxial coil for inductance calculations.

**Parameters:**

| Parameter      | Type  | Description                                |
|----------------|-------|--------------------------------------------|
| `inner_radius` | float | Inner radius of the coil (meters)          |
| `outer_radius` | float | Outer radius of the coil (meters)          |
| `length`       | float | Axial length of the coil (meters)          |
| `num_turns`    | float | Number of wire turns                       |

#### `SumPrecision`

Precision parameters for the series expansion computation.

**Parameters:**

| Parameter | Type | Default | Description                       |
|-----------|------|---------|-----------------------------------|
| `k_terms` | int  | 32      | Terms for k-index series          |
| `l_terms` | int  | 32      | Terms for l-index series          |
| `n_terms` | int  | 64      | Terms for n-index series          |

#### `MutualInductanceCalculator`

Main calculator class for mutual inductance computations.

**Methods:**

- `calculate(coil1, coil2, distance, verbose=False)` - Calculate mutual inductance (auto-selects method)
- `calculate_near(coil1, coil2, distance, z_parameter)` - Near-field calculation
- `calculate_far(coil1, coil2, distance)` - Far-field calculation
- `benchmark_near(n_repeats=10000)` - Benchmark near-field method
- `benchmark_far(n_repeats=10000)` - Benchmark far-field method

### Functions

#### `compile_library(force=False, verbose=True)`

Manually compile the C library.

#### `detect_simd_capability()`

Detect the highest SIMD capability available on the current CPU.

#### `get_compile_flags(simd=None)`

Generate compiler flags for the current platform and SIMD capability.

#### `get_platform_info()`

Get information about the current platform (system, machine, SIMD capability).

## Platform Support

The module automatically detects the operating system and CPU capabilities to apply optimal compiler flags:

| Platform       | Compiler | SIMD Detection                      |
|----------------|----------|--------------------------------------|
| Linux          | GCC      | Via `/proc/cpuinfo`                  |
| macOS (Intel)  | Clang    | Via `sysctl` CPU feature queries     |
| macOS (Apple Silicon) | Clang | Auto-detects ARM64, uses NEON   |
| Windows        | MSVC     | Defaults to AVX2 on 64-bit systems   |

### Apple Silicon Support

On Apple Silicon Macs (M1, M2, M3, etc.), the library:
- Automatically detects ARM64 architecture
- Uses NEON SIMD instructions (always available on ARM64)
- Compiles with `-march=native` for optimal performance

## Performance Tips

1. **SIMD Optimization**: The library automatically detects and uses the highest available SIMD instruction set (SSE, AVX2, AVX512).

2. **Precision Selection**: For optimal accuracy:
   - Use `k_terms == l_terms` in general
   - Use `n_terms == 2 * k_terms` for near-field calculations
   - Use `n_terms == k_terms` for far-field calculations

3. **Force Recompilation**: If you've upgraded your CPU or changed compiler settings:
   ```python
   from vector_case_py import compile_library
   compile_library(force=True)
   ```

## Example: Calculating Inductance at Multiple Distances

```python
import numpy as np
from vector_case_py import Coil, MutualInductanceCalculator

# Define coils
coil1 = Coil(0.1, 0.2, 0.1, 100)
coil2 = Coil(0.3, 0.4, 0.1, 100)

calculator = MutualInductanceCalculator()

# Calculate at multiple distances
distances = np.linspace(0.05, 0.5, 20)
inductances = [calculator.calculate(coil1, coil2, d) for d in distances]

# Plot results (requires matplotlib)
import matplotlib.pyplot as plt
plt.plot(distances, inductances)
plt.xlabel('Distance (m)')
plt.ylabel('Mutual Inductance (H)')
plt.title('Mutual Inductance vs. Distance')
plt.show()
```

Note: `numpy` and `matplotlib` are not required dependencies, so they must be installed separately.

## License

MIT License - see the main project LICENSE file.
