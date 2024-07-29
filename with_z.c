#pragma clang diagnostic push
#pragma ide diagnostic ignored "portability-simd-intrinsics"

#include <stdio.h>

#include "structs.h"

#include "inductance_near.h"
#include "inductance_fast.h"
#include "inductance_remainder_bound.h"


int main() {
    // Input parameters for the problem
    CoilCalculationData data;

    data.r_1 = 0.5;
    data.R_1 = 1.5;
    data.L_1 = 1.0;
    data.N_1 = 1.0;

    data.r_2 = 0.5;
    data.R_2 = 1.5;
    data.L_2 = 1.0;
    data.N_2 = 1.0;

    const FP_TYPE d = 1.0;

    // Data structure to hold the input parameters

    SumPrecisionData precision;
    precision.k_terms = 24;
    precision.l_terms = 24;
    precision.n_terms = 48;

    for (int i = -10; i < 200; i++) {
        FP_TYPE M_12 = calculate_mutual_inductance_near(data, precision, d, (FP_TYPE) i * 0.01);
        printf("%f\t%.16g\n", (double) i * 0.01, M_12);
    }

//    for (int i = 4; i < 32; ++i) {
//        for (int j = 4; j < 32; ++j) {
//            SumPrecisionData prec{(uint32_t) i, (uint32_t) j, (uint32_t) std::round(1.4 * std::max(i, j))};
//            FP_TYPE M_12 = calculate_mutual_inductance(data, prec, 0.4, 0.6);
//            std::cout << std::setprecision(16) << M_12 << "\t";
//        }
//        std::cout << std::endl;
//    }

//    for (int i = 1; i <= 43; ++i) {
//        SumPrecisionData loc_prec{18, 18, (uint32_t) i};
//        FP_TYPE M_12 = calculate_mutual_inductance(data, loc_prec, d, 0.6);
//        std::cout << std::setprecision(16) << M_12 << std::endl;
//    }

    benchmark_mutual_inductance_near(precision, 10000);

    // Calculate the mutual inductance
    FP_TYPE M_12 = calculate_mutual_inductance_near(data, precision, d, 6.0 * data.R_1);
    FP_TYPE M_12_1 = guess_best_inductance_near(data, precision, d, 0.0, 6 * data.R_1, true, 1e-5);

    // Print to 16 decimal places
    printf("M_12 = %.16g\n", M_12);
    printf("M_12_1 = %.16g\n", M_12_1);

    printf("Remainder: %.16g\n", calculate_inductance_remainder_unoptimized(data, precision, d, 3 * data.R_2, true));

    return 0;
}

#pragma clang diagnostic pop