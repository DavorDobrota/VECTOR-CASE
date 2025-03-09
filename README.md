# VECTOR-CASE
Vectotized Efficient C Tool fOr Rectangular-coil Coaxial Analytical Series Expansion of mutual inductance

## Introduction
This is the code that accompanies the paper "Calculation of Mutual Inductance of Two Coaxial Thick Coils With
Rectangular Cross Section by Using Cylindrical Multipole Expansion",
The paper is available at [IEEE Xplore](https://ieeexplore.ieee.org/document/10856207), and can be cited as:
```
@ARTICLE{10856207,
    author={Vučić, Filip and Dobrota, Davor},
    journal={IEEE Transactions on Magnetics}, 
    title={Calculation of Mutual Inductance of Two Coaxial Thick Coils With Rectangular Cross Section by Using Cylindrical Multipole Expansion}, 
    year={2025},
    volume={61},
    number={3},
    pages={1-16},
    doi={10.1109/TMAG.2025.3535634}
}
```

The code is meant to be a header-only C library. This makes it easy to include in other projects, and to use in
conjunction with other libraries as well as from other languages. For an example of this, see how to use the code
from Python.

## Use Case
The code is written in C, and is intended to be used as header-only library. It is vectorized, and 
intended to be used with SSE, AVX, or AVX-512 instruction sets. To enable them, uncomment the appropriate line
in src/settings.h or pass an appropriate definition to the compiler by some other means.
Note that the most aggressive of the uncommented options is AVX-512 and it will be used even if other lines are
uncommented. Be mindful whether your CPU supports the instruction set you are targeting.

One can generate lookup tables for the sums using the functions found in 
[generate_lookup_tables.py](./src/generate_lookup_tables.py). We provide headers with 64 terms along each 
direction as we find this a good balance between size and precision .


## Python Codes
There is also an implementation of the method in Python, using the Decimal class to achieve adaptive precision.
This is at the cost of some performance, but eliminates all fears of overflow when using a large number of terms. 

There is a demonstration of how the methods from the header files can be used in Python by means of the CFFI
library in [C_interface_demo.py](./demo/C_interface_demo.py).
Additionally one might consider using cppyy, but at a potential loss of performance due to more limited
control of the compiler and compilation flags. As far as we are aware it is not simple to pass 
`-march=native` or equivalent flags. However, performance will still be much better than any other
approach even without the most aggresive optimization (see the paper for more details).

## Validation by External Method
Finally, we implemented the method by Župan et al for comparison:

T. Župan, Ž. Štih, and B. Trkulja, “Fast and precise method for  inductance calculation of coaxial circular coils 
with rectangular cross section using the one-dimensional integration of elementary functions  applicable 
to superconducting magnets” ([link to the paper](https://ieeexplore.ieee.org/document/6718014))

We provide implementations in C and Mathematica. Note that the C implementation is very fast but due to the 
oscillatory nature of the integrand the Gauss-Legendre quadrature does not converge. The Gauss-Kronrod 
adaptive quadrature used by Mathematica seems to be more robust, but it too fails for some configurations.

