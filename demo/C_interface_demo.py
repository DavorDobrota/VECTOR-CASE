# -*- coding: utf-8 -*-
"""Calling C code from Python using CFFI
@Author: Davor Dobrota

This script demonstrates how to use the C interface of the inductance calculation library
from Python. The library is compiled using the CFFI library, which is a foreign function
interface for Python calling C code. The code is compiled automatically if vector_case
module is not found in the current directory. CFFI can easily be installed using pip.

Note that the performance flags are commented out in the CFFI compilation, but they can be
introduced for better performance. If your processors supports SSE, AVX2 or AVX512, you can
go to settings.py and uncomment the corresponding (highest capability) flag. If you do this
after you have already compiled the library, you need delete the generated files and recompile
the library.

Feel free to use this as a template for your own projects.
"""

# Make sure the package cffi is installed using pip or other means
from cffi import FFI
import importlib.util

module_name = "vector_case"


if not (importlib.util.find_spec(module_name) is not None):

    ffi = FFI()

    ffi.cdef(
        """
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
    
    double calculate_mutual_inductance(
        const CoilCalculationData data,
        const SumPrecisionData precision,
        const double d,
        const bool verbose
    );
    
    double calculate_mutual_inductance_raw(
        const double  r_1,
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
    
    void benchmark_mutual_inductance_far(const SumPrecisionData precision, const uint32_t n_repeats);
    void benchmark_mutual_inductance_near(const SumPrecisionData precision, const uint32_t n_repeats);
    """
    )

    ffi.set_source(
        module_name,        # Mind the path here if you change it
        """
        #include "../src/inductance.h" 
        """,
        libraries=["m"],    # Link with math library, needed for linux
        extra_compile_args=[
            # "-std=c11",
            # "-O3",
            # "-march=native",
            # "-mtune=native"
        ],
    )

    ffi.compile(verbose=True)

spec = importlib.util.find_spec("vector_case")
vector_case = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vector_case)


if __name__ == "__main__":

    r_1 = 0.1
    R_1 = 0.2
    L_1 = 0.1
    N_1 = 100

    r_2 = 0.3
    R_2 = 0.4
    L_2 = 0.1
    N_2 = 100

    l_terms = 32
    k_terms = 32
    n_terms = 64

    d = 0.1

    # Using directly the C function
    inductance = vector_case.lib.calculate_mutual_inductance_raw(
        r_1, R_1, L_1, N_1, r_2, R_2, L_2, N_2, d, k_terms, l_terms, n_terms, True
    )
    print("Calculated Inductance:", inductance)

    # Using the struct for a more elaborate example
    coil_data = vector_case.ffi.new(
        "CoilCalculationData *",
        {
            'r_1': r_1, 'R_1': R_1, 'L_1': L_1, 'N_1': N_1,
            'r_2': r_2, 'R_2': R_2, 'L_2': L_2, 'N_2': N_2
        }
    )
    precision_data = vector_case.ffi.new(
        "SumPrecisionData *", {'n_terms': n_terms, 'l_terms': l_terms, 'k_terms': k_terms}
    )
    Z = 1.4
    inductance = vector_case.lib.calculate_mutual_inductance_near(coil_data[0], precision_data[0], d, Z)
    print("Calculated Inductance manually:", inductance)

    # Finally, a benchmark
    vector_case.lib.benchmark_mutual_inductance_near(precision_data[0], 10000)
