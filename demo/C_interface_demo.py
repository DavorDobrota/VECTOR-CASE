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
    
    void benchmark_mutual_inductance_far(const SumPrecisionData precision, const uint32_t n_repeats);
    void benchmark_mutual_inductance_near(const SumPrecisionData precision, const uint32_t n_repeats);
    """
    )

    ffi.set_source(
        module_name,
        """
        #include "inductance_near.h" 
        #include "inductance_far.h"
        """,
        libraries=["m"],
        extra_compile_args=["-std=c11", "-O3", "-march=native", "-mtune=native"],
    )

    ffi.compile(verbose=True)

spec = importlib.util.find_spec("vector_case")
vector_case = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vector_case)


if __name__ == "__main__":

    coil_data = vector_case.ffi.new(
        "CoilCalculationData *",
        {'r_1': 0.1, 'R_1': 0.2, 'L_1': 0.1, 'N_1': 100,
         'r_2': 0.3, 'R_2': 0.4, 'L_2': 0.1, 'N_2': 100}
    )
    precision_data = vector_case.ffi.new(
        "SumPrecisionData *", {'n_terms': 32, 'l_terms': 32, 'k_terms': 64}
    )

    d = 0.1
    Z = 1.0
    inductance = vector_case.lib.calculate_mutual_inductance_near(coil_data[0], precision_data[0], d, Z)
    print("Calculated Inductance:", inductance)

    vector_case.lib.benchmark_mutual_inductance_near(precision_data[0], 10000)
