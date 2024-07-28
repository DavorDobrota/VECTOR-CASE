#pragma clang diagnostic push
#pragma ide diagnostic ignored "portability-simd-intrinsics"

#include <iostream>
#include <iomanip>
#include <chrono>

#include "structs.h"
#include "factorial_lookup.h"

#include "inductance_near.hpp"
#include "inductance_fast.hpp"
#include "remainder_estimate.hpp"


//FP_TYPE calculate_mutual_inductance(
//        const CoilCalculationData& data,
//        const SumPrecisionData &precision,
//        const FP_TYPE d,
//        const FP_TYPE Z,
//        bool timing = false
//) {
//    // Useful calculations that can be performed at compile time, before the main loop
//    const FP_TYPE denom_1 = (FP_TYPE) 1.0 / (Z + data.L_1 + d);
//    const FP_TYPE denom_2 = (FP_TYPE) 1.0 / (Z + data.L_1 + data.L_2 + d);
//
//    const FP_TYPE R_1_sq = data.R_1 * data.R_1;
//    const FP_TYPE r_1_sq = data.r_1 * data.r_1;
//
//    const FP_TYPE R_2_sq = data.R_2 * data.R_2;
//    const FP_TYPE r_2_sq = data.r_2 * data.r_2;
//
//    const FP_TYPE denom_1_sq = denom_1 * denom_1;
//    const FP_TYPE denom_2_sq = denom_2 * denom_2;
//
//    const FP_TYPE L_1_plus_Z = data.L_1 + Z;
//
//    // Initialize timer
//    std::chrono::high_resolution_clock::time_point begin_time;
//
//    if (timing) {
//        begin_time = std::chrono::high_resolution_clock::now();
//    }
//
//    // Variable which will store the result of the main loop
//    FP_TYPE M_12 = 0.0;
//
//    // These three variables need to be remembered for each iteration
//    // of the inner loops to be able to restore index values
//    FP_TYPE loop_denom_1 = denom_1_sq;
//    FP_TYPE loop_denom_2 = denom_2_sq;
//    FP_TYPE power_2 = 4.0;
//
//    // Values that have to be set for the first loop
//    FP_TYPE loop_R_1 = R_1_sq * data.R_1;
//    FP_TYPE loop_r_1 = r_1_sq * data.r_1;
//
//    // Main loop
//    for (uint32_t l = 0; l < precision.l_terms; ++l) {
//
//        if (l > 0) {
//            loop_R_1 *= R_1_sq;
//            loop_r_1 *= r_1_sq;
//
//            loop_denom_1 *= denom_1_sq;
//            loop_denom_2 *= denom_2_sq;
//
//            power_2 *= 4.0;
//        }
//
//        FP_TYPE loop_R_1_sub_r_1 = loop_R_1 - loop_r_1;
//
//        // Restore point for first loop
//        FP_TYPE save_first_loop_denom_1 = loop_denom_1;
//        FP_TYPE save_first_loop_denom_2 = loop_denom_2;
//        FP_TYPE save_first_power_2 = power_2;
//
//        // Values that have to be set for the second loop
//        FP_TYPE loop_R_2 = R_2_sq * data.R_2;
//        FP_TYPE loop_r_2 = r_2_sq * data.r_2;
//
//        for (uint32_t k = 0; k < precision.k_terms; ++k) {
//            if (k > 0) {
//                loop_R_2 *= R_2_sq;
//                loop_r_2 *= r_2_sq;
//
//                loop_denom_1 *= denom_1_sq;
//                loop_denom_2 *= denom_2_sq;
//
//                power_2 *= 4.0;
//            }
//
//            FP_TYPE loop_R_2_sub_r_2 = loop_R_2 - loop_r_2;
//
//            FP_TYPE sign = 1.0;
//            if ((l + k) % 2 == 1) {
//                sign = -1.0;
//            }
//            auto factor = (FP_TYPE) ((2 * l + 3) * (2 * k + 3));
//            FP_TYPE table_value = factorial_array[l] * factorial_array[l + 1]
//                                * factorial_array[k] * factorial_array[k + 1];
//
//            // Restore point for second loop
//            FP_TYPE save_second_loop_denom_1 = loop_denom_1;
//            FP_TYPE save_second_loop_denom_2 = loop_denom_2;
//            FP_TYPE save_second_power_2 = power_2;
//
//            // Values that have to be set for the third loop
//            FP_TYPE loop_L1_plus_Z = L_1_plus_Z;
//            FP_TYPE loop_Z = Z;
//
//            for (uint32_t n = 0; n < precision.max_n_term; ++n) {
//
//                if (n > 0) {
//                    loop_denom_1 *= denom_1;
//                    loop_denom_2 *= denom_2;
//
//                    loop_L1_plus_Z *= L_1_plus_Z;
//                    loop_Z *= Z;
//                }
//
//                FP_TYPE numerator = sign
//                                  * factorial_array[2 * l + 2 * k + n + 1]
//                                  * (loop_L1_plus_Z - loop_Z)
//                                  * loop_R_1_sub_r_1
//                                  * loop_R_2_sub_r_2
//                                  * (loop_denom_1 - loop_denom_2);
//
//                FP_TYPE denominator = power_2
//                                    * factor
//                                    * table_value
//                                    * factorial_array[n + 1];
//
//                M_12 += numerator / denominator;
//            }
//
//            loop_denom_1 = save_second_loop_denom_1;
//            loop_denom_2 = save_second_loop_denom_2;
//            power_2 = save_second_power_2;
//        }
//
//        loop_denom_1 = save_first_loop_denom_1;
//        loop_denom_2 = save_first_loop_denom_2;
//        power_2 = save_first_power_2;
//    }
//
//    M_12 *= 4.0 * local_pi * local_pi * 1e-7 * data.N_1 * data.N_2
//          / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));
//
//    if (timing) {
//        double interval = duration_cast<std::chrono::duration<double>>(
//                std::chrono::high_resolution_clock::now() - begin_time).count();
//        std::cout << "Optimized Time = " << interval << " s" << std::endl;
//    }
//
//    return M_12;
//}

//FP_TYPE calculate_mutual_inductance(
//        const CoilCalculationData& data,
//        const SumPrecisionData &precision,
//        const FP_TYPE d,
//        const FP_TYPE Z,
//        bool timing = false
//) {
//    // Useful calculations that can be performed at compile time, before the main loop
//    const FP_TYPE denom_1 = (FP_TYPE) 1.0 / (Z + data.L_1 + d);
//    const FP_TYPE denom_2 = (FP_TYPE) 1.0 / (Z + data.L_1 + data.L_2 + d);
//
//    const FP_TYPE R_1_sq = data.R_1 * data.R_1;
//    const FP_TYPE r_1_sq = data.r_1 * data.r_1;
//
//    const FP_TYPE R_2_sq = data.R_2 * data.R_2;
//    const FP_TYPE r_2_sq = data.r_2 * data.r_2;
//
//    const FP_TYPE denom_1_sq = denom_1 * denom_1;
//    const FP_TYPE denom_2_sq = denom_2 * denom_2;
//
//    const FP_TYPE L_1_plus_Z = data.L_1 + Z;
//
//    // Initialize timer
//    std::chrono::high_resolution_clock::time_point begin_time;
//
//    if (timing) {
//        begin_time = std::chrono::high_resolution_clock::now();
//    }
//
//    // Variable which will store the result of the main loop
//    FP_TYPE M_12 = 0.0;
//
//    // These three variables need to be remembered for each iteration
//    // of the inner loops to be able to restore index values
//    FP_TYPE loop_denom_1 = denom_1_sq;
//    FP_TYPE loop_denom_2 = denom_2_sq;
//
//    // Values that have to be set for the first loop
//    FP_TYPE loop_R_1 = R_1_sq * data.R_1;
//    FP_TYPE loop_r_1 = r_1_sq * data.r_1;
//
//    // Main loop
//    for (uint32_t l = 0; l < precision.l_terms; ++l) {
//
//        if (l > 0) {
//            loop_R_1 *= R_1_sq;
//            loop_r_1 *= r_1_sq;
//
//            loop_denom_1 *= denom_1_sq;
//            loop_denom_2 *= denom_2_sq;
//        }
//
//        FP_TYPE loop_R_1_sub_r_1 = loop_R_1 - loop_r_1;
//
//        // Restore point for first loop
//        FP_TYPE save_first_loop_denom_1 = loop_denom_1;
//        FP_TYPE save_first_loop_denom_2 = loop_denom_2;
//
//        // Values that have to be set for the second loop
//        FP_TYPE loop_R_2 = R_2_sq * data.R_2;
//        FP_TYPE loop_r_2 = r_2_sq * data.r_2;
//
//        for (uint32_t k = 0; k < precision.k_terms; ++k) {
//            if (k > 0) {
//                loop_R_2 *= R_2_sq;
//                loop_r_2 *= r_2_sq;
//
//                loop_denom_1 *= denom_1_sq;
//                loop_denom_2 *= denom_2_sq;
//            }
//
//            FP_TYPE loop_R_2_sub_r_2 = loop_R_2 - loop_r_2;
//
//            // Restore point for second loop
//            FP_TYPE save_second_loop_denom_1 = loop_denom_1;
//            FP_TYPE save_second_loop_denom_2 = loop_denom_2;
//
//            // Values that have to be set for the third loop
//            FP_TYPE loop_L1_plus_Z = L_1_plus_Z;
//            FP_TYPE loop_Z = Z;
//
//            for (uint32_t n = 0; n < precision.max_n_term; ++n) {
//
//                if (n > 0) {
//                    loop_denom_1 *= denom_1;
//                    loop_denom_2 *= denom_2;
//
//                    loop_L1_plus_Z *= L_1_plus_Z;
//                    loop_Z *= Z;
//                }
//
//                FP_TYPE numerator = (loop_L1_plus_Z - loop_Z)
//                                  * loop_R_1_sub_r_1
//                                  * loop_R_2_sub_r_2
//                                  * (loop_denom_1 - loop_denom_2);
//
//                M_12 += numerator * lookup_table_near[l][k][n];
//            }
//
//            loop_denom_1 = save_second_loop_denom_1;
//            loop_denom_2 = save_second_loop_denom_2;
//        }
//
//        loop_denom_1 = save_first_loop_denom_1;
//        loop_denom_2 = save_first_loop_denom_2;
//    }
//
//    M_12 *= 4.0 * local_pi * local_pi * 1e-7 * data.N_1 * data.N_2
//          / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));
//
//    if (timing) {
//        double interval = std::chrono::duration_cast<std::chrono::duration<double>>(
//                std::chrono::high_resolution_clock::now() - begin_time).count();
//        std::cout << "Lookup optimized time = " << interval << " s" << std::endl;
//    }
//
//    return M_12;
//}


int main() {
    // Input parameters for the problem
    const FP_TYPE r_1 = 1.2;
    const FP_TYPE R_1 = 1.25;
    const FP_TYPE L_1 = 0.05;
    const FP_TYPE N_1 = 100.0;

    const FP_TYPE r_2 = 0.05;
    const FP_TYPE R_2 = 0.1;
    const FP_TYPE L_2 = 0.05;
    const FP_TYPE N_2 = 100.0;

    const FP_TYPE d = 1.225;

    // Data structure to hold the input parameters
    CoilCalculationData data{N_1, L_1, R_1, r_1, N_2, L_2, R_2, r_2};
    SumPrecisionData precision{8, 8, 16};

    for (int i = -10; i < 200; i++) {
        // Calculate the mutual inductance
        FP_TYPE M_12 = calculate_mutual_inductance_near(data, precision, d, (FP_TYPE) i * 0.01);

        // Print to 16 decimal places
        std::cout << (float) i * 0.01 << "\t" << std::setprecision(15) << M_12 << std::endl;
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

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    FP_TYPE volatile value;
//    float volatile value_float;
    for (int i = 0; i < 10000; ++i) {
//        value = calculate_mutual_inductance_double_fast(
//                data.N_1, data.L_1, data.R_1, data.r_1,
//                data.N_2, data.L_2, data.R_2, data.r_2,
//                (double) i * 0.0001, 1.0);
        value = calculate_mutual_inductance_near(data, precision, (FP_TYPE) i * 0.0001f, 1.0);
//        value_float = calculate_mutual_inductance_float_fast_avx(
//                data.N_1, data.L_1, data.R_1, data.r_1,
//                data.N_2, data.L_2, data.R_2, data.r_2,
//                (float) i * 0.0001f, 1.0);
    }

    double interval = std::chrono::duration_cast<std::chrono::duration<double>>(
            std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Value = " << value << std::endl;
//    std::cout << "Float Value = " << value_float << std::endl;
    std::cout << "Fast Time average = " << interval / 10000 << " s" << std::endl;

    // Calculate the mutual inductance
    FP_TYPE M_12 = calculate_mutual_inductance_near(data, precision, d, 6 * data.R_1);
    double M_12_fast = calculate_mutual_inductance_double_fast(
            data.N_1, data.L_1, data.R_1, data.r_1,
            data.N_2, data.L_2, data.R_2, data.r_2,
            d, 6 * data.R_1);
    FP_TYPE M_12_1 = guess_best_inductance_near(data, precision, d, 0.0, 6 * data.R_1, true);

    // Print to 16 decimal places
    std::cout << std::endl << "M_12 =   " << std::setprecision(16) << M_12 << std::endl;
    std::cout << "M_12_fast = " << std::setprecision(16) << M_12_fast << std::endl;
    std::cout << "M_12_1 = " << std::setprecision(16) << M_12_1 << std::endl;

    std::cout << "Remainder: " << calculate_inductance_remainder_unoptimized(data, precision, d, 3 * data.R_2, true) << std::endl;

    return 0;
}

#pragma clang diagnostic pop