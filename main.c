#include <stdio.h>

#include "inductance_naive.h"
#include "inductance_far.h"
#include "inductance_zupan.h"

int main() {
    CoilCalculationData data;

    // Input parameters for the problem
    data.r_1 = 0.1;
    data.R_1 = 0.2;
    data.L_1 = 0.1;
    data.N_1 = 100.0;

    data.r_2 = 0.3;
    data.R_2 = 0.4;
    data.L_2 = 0.1;
    data.N_2 = 100.0;

    const FP_TYPE d = 0.9;

    // Precision of the sum
    SumPrecisionData precision;
    precision.k_terms = 32;
    precision.l_terms = 32;
    precision.n_terms = 32;

    // Calculate the mutual inductance
    FP_TYPE M_12 = calculate_mutual_inductance_far(data, precision, d);
    FP_TYPE M_12_unoptimized = calculate_mutual_inductance_far_unoptimized(data, precision, d, true);
//    FP_TYPE M_12_zupan = calculate_mutual_inductance_zupan(data, d, 1, 15);

    benchmark_mutual_inductance_far(precision, 10000);

    // Print to 16 decimal places
    printf("M_12              = %.16g\n", M_12);
    printf("M_12 Unoptimized  = %.16g\n", M_12_unoptimized);
//    printf("M_12 Zupan        = %.16g\n", M_12_zupan);

    return 0;
}
