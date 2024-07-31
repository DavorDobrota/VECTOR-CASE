# VECTOR-CASE
Vectorized Efficient C Tool for Analytical Series Expansion of Mutual Inductance 
for Circular Coils with Rectangular Cross-Section in Coaxial Configuration

## Introduction
This is the code that accompanies the paper "Calculation of mutual inductance of two coaxial thick coils with
rectangular cross section by using cylindrical multipole expansion", by Filip Vučič, and Davor Dobrota.
The paper available at (to be determined), and can be cited as (to be determined).

## Use
The code is written in C, and is intended to be used as header-only library. The code is vectorized, and 
is intended to be used with SSE, AVX, or AVX-512 instruction sets. To enable them, uncomment the appropriate line
in src/settings.h or pass an appropriate definition to the compiler by some other means.

One cn generate lookup tables for the sums using the functions found in generate_lookup_tables.py 

There is also an implementation of the method in Python, using the Decimal class to achieve adaptive precision.

Finally, an implementation of method by Župan et al (
T. Župan, Ž. Štih, and B. Trkulja, “Fast and precise method for  inductance calculation of coaxial circular coils 
with rectangular cross section using the one-dimensional integration of elementary functions  applicable 
to superconducting magnets,” IEEE Transactions on Applied Superconductivity, vol. 24, no. 2, pp. 81–89, 2014.)
is implemented in C (but does not converge with Gauss-Legendre) and Mathematica. 

There is a demonstration of how the methods from the header files can be used in Python bz means of the CFFI 
library. Additionally one might consider using cppyy, but at a potential loss of performance due to more limited
control of the compiler and compilation flags. 