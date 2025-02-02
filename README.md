# VECTOR-CASE
Vectorized Efficient C Tool for Analytical Series Expansion of Mutual Inductance 
for Circular Coils with Rectangular Cross-Section in Coaxial Configuration

## Introduction
This is the code that accompanies the paper "Calculation of mutual inductance of two coaxial thick coils with
rectangular cross-section by using cylindrical multipole expansion", by Filip Vučič, and Davor Dobrota.
The paper available at IEEE Xplore (https://ieeexplore.ieee.org/document/10856207), and can be cited as (preprint):
```
@ARTICLE{10856207,
  author={Vučić, Filip and Dobrota, Davor},
  journal={IEEE Transactions on Magnetics}, 
  title={Calculation of mutual inductance of two coaxial thick coils with rectangular cross section by using cylindrical multipole expansion}, 
  year={2025},
  volume={},
  number={},
  pages={1-1},
  doi={10.1109/TMAG.2025.3535634}}
```

The code is meant to be a header-only C library. This makes it easy to include in other projects, and to use in
conjunction with other libraries as well as from other languages. For an example of this, see how to use the code
from Python.

## Use

### Base use-case
The code is written in C, and is intended to be used as header-only library. The code is vectorized, and 
is intended to be used with SSE, AVX, or AVX-512 instruction sets. To enable them, uncomment the appropriate line
in src/settings.h or pass an appropriate definition to the compiler by some other means.

One can generate lookup tables for the sums using the functions found in generate_lookup_tables.py 


### Python
There is also an implementation of the method in Python, using the Decimal class to achieve adaptive precision.

There is a demonstration of how the methods from the header files can be used in Python bz means of the CFFI
library. Additionally one might consider using cppyy, but at a potential loss of performance due to more limited
control of the compiler and compilation flags.

### Mathematica and implementation of external method
Finally, an implementation of method by Župan et al (
T. Župan, Ž. Štih, and B. Trkulja, “Fast and precise method for  inductance calculation of coaxial circular coils 
with rectangular cross section using the one-dimensional integration of elementary functions  applicable 
to superconducting magnets,” IEEE Transactions on Applied Superconductivity, vol. 24, no. 2, pp. 81–89, 2014.)
is implemented in C (but does not converge with Gauss-Legendre) and Mathematica. 

