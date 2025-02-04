#include <stdio.h>

#include "inductance_naive.h"
#include "inductance.h"
#include "inductance_zupan.h"
#include "inductance_remainder_bound.h"

/**
 * @brief Main function to test the inductance calculation functions. It is not a
 * main function in a classical sense, but rather meant to showcase the usage of the
 * library.
 */
int main() {
    CoilCalculationData data;

    // Input parameters for the problem
    data.r_1 = 1.0;
    data.R_1 = 2.0;
    data.L_1 = 1.0;
    data.N_1 = 100.0;

    data.r_2 = 1.5;
    data.R_2 = 4.0;
    data.L_2 = 1.0;
    data.N_2 = 100.0;

    const FP_TYPE d = 4e-6;

    // Precision of the sum
    SumPrecisionData precision;
    precision.k_terms = 12;
    precision.l_terms = 12;
    precision.n_terms = 12;

    // Calculate the mutual inductance
    FP_TYPE M_12 = calculate_mutual_inductance_far(data, precision, d);
    FP_TYPE M_12_unoptimized = calculate_mutual_inductance_far_unoptimized(data, precision, d, true);
//    FP_TYPE M_12_zupan = calculate_mutual_inductance_zupan(data, d, 1, 15);

    benchmark_mutual_inductance_far(precision, 10000);

    // Print to 16 decimal places
    printf("M_12              = %.16g\n", M_12);
    printf("M_12 Unoptimized  = %.16g\n", M_12_unoptimized);
//    printf("M_12 Zupan        = %.16g\n", M_12_zupan);

    // Input parameters for the problem

//    data.r_1 = 0.1;
//    data.R_1 = 0.2;
//    data.L_1 = 0.1;
//    data.N_1 = 100.0;
//
//    data.r_2 = 0.3;
//    data.R_2 = 0.4;
//    data.L_2 = 0.1;
//    data.N_2 = 100.0;

    const FP_TYPE d_near = d;

    // Data structure to hold the input parameters

    precision.k_terms = 32;
    precision.l_terms = 32;
    precision.n_terms = 64;

    for (int i = -10; i < 150; i++) {
        FP_TYPE M_12 = calculate_mutual_inductance_near(data, precision, d_near, (FP_TYPE) i * 0.1);
        printf("%f\t%.16g\n", (double) i * 0.1, M_12);
    }

//    for (int i = 1; i <= 43; ++i) {
//        SumPrecisionData loc_prec{18, 18, (uint32_t) i};
//        FP_TYPE M_12 = calculate_mutual_inductance(data, loc_prec, d, 0.6);
//        std::cout << std::setprecision(16) << M_12 << std::endl;
//    }

    FP_TYPE R = data.R_1 > data.R_2 ? data.R_1 : data.R_2;
    benchmark_mutual_inductance_near(precision, 100000);

    // Calculate the mutual inductance
    M_12 = calculate_mutual_inductance_near(data, precision, d, 6.0 * data.R_1);
    M_12_unoptimized = calculate_mutual_inductance_near_unoptimized(
        data, precision, d_near, 4.0 * R, true
    );
//    FP_TYPE M_12_1 = guess_best_inductance_near(
//        data, precision, d_near, 0.0, 5.0 * R, true, 1e-6
//    );
    FP_TYPE M_12_1 = calculate_mutual_inductance(data, precision, d, true);

    // Print to 16 decimal places
    printf("M_12 = %.16g\n", M_12);
    printf("M_12_1 = %.16g\n", M_12_1);
    printf("M_12 Unoptimized = %.16g\n", M_12_unoptimized);
    printf("Remainder: %.16g\n",
           calculate_inductance_remainder_unoptimized(data, precision, d, 3.0 * data.R_2, true)
    );

    return 0;
}
