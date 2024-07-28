#ifndef VECTOR_CASE_INDUCTANCE_NAIVE_HPP
#define VECTOR_CASE_INDUCTANCE_NAIVE_HPP

#include <iostream>
#include <chrono>
#include <cmath>

#include "structs.h"
#include "factorial_lookup.h"


// TODO - fix whatever the bug is here, no idea right now, optimized method works
double calculate_mutual_inductance_far_unoptimized(
        const CoilCalculationData& data,
        const SumPrecisionData &precision,
        const FP_TYPE d,
        bool timing = false
) {
    FP_TYPE M_12 = 0.0;

    std::chrono::high_resolution_clock::time_point begin_time;

    if (timing) {
        begin_time = std::chrono::high_resolution_clock::now();
    }

    for (uint32_t n = 0; n < precision.n_terms; ++n) {
        for (uint32_t l = 0; l < precision.l_terms; ++l) {
            for (uint32_t k = 0; k < precision.k_terms; ++k) {

                FP_TYPE sign = (l + k) % 2 == 0 ? 1.0 : -1.0;

                FP_TYPE numerator = sign
                                  * factorial_array[2 * n + 2 * l + 2 * k + 1]
                                  * std::pow(data.L_1, 2 * n)
                                  * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                                  * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                                  * (1.0 / std::pow(0.5 * data.L_1 + d, 2 * n + 2 * l + 2 * k + 2)
                                     - 1.0 / std::pow(0.5 * data.L_1 + d + data.L_2, 2 * n + 2 * l + 2 * k + 2)
                                    );

                FP_TYPE denominator = std::pow(2.0, 2 * n + 2 * l + 2 * k + 2)
                                    * ((FP_TYPE) ((2 * l + 3) * (2 * k + 3)))
                                    * factorial_array[l] * factorial_array[l + 1]
                                    * factorial_array[k] * factorial_array[k + 1]
                                    * factorial_array[2 * n + 1];

                FP_TYPE term = numerator / denominator;

                if (!std::isfinite(term)) break;

                M_12 += term;
            }
        }
    }

    M_12 *= (FP_TYPE) (4.0e-7) * local_pi * local_pi * data.N_1 * data.N_2
            / (data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    if (timing) {
        double interval = std::chrono::duration_cast<std::chrono::duration<double>>(
                std::chrono::high_resolution_clock::now() - begin_time).count();
        std::cout << "Unoptimized Time = " << interval << " s" << std::endl;
    }

    return M_12;
}


double calculate_mutual_inductance_near_unoptimized_normed(
        const CoilCalculationData& data,
        const SumPrecisionData &precision,
        const FP_TYPE d,
        const FP_TYPE Z,
        bool timing = false
) {
    FP_TYPE M_12 = 0.0;

    std::chrono::high_resolution_clock::time_point begin_time;

    if (timing) {
        begin_time = std::chrono::high_resolution_clock::now();
    }

    for (uint32_t n = 0; n < precision.n_terms; ++n) {
        for (uint32_t l = 0; l < precision.l_terms; ++l) {
            for (uint32_t k = 0; k < precision.k_terms; ++k) {

                FP_TYPE sign = (l + k) % 2 == 0 ? 1.0 : -1.0;

                FP_TYPE numerator = sign
                                    * factorial_array[2 * l + 2 * k + n + 1]
                                    * (std::pow(1.0 + data.L_1 / Z, n + 1) - 1.0)
                                    * (std::pow(data.R_1 / Z, 2 * l + 3) - std::pow(data.r_1 / Z, 2 * l + 3))
                                    * (std::pow(data.R_2 / Z, 2 * k + 3) - std::pow(data.r_2 / Z, 2 * k + 3))
                                    * (std::pow(1.0 + (data.L_1 + data.L_2 + d) / Z, 2 * l + 2 * k + n + 2)
                                       - std::pow(1.0 + (data.L_1 + d) / Z, 2 * l + 2 * k + n + 2)
                                    );

                FP_TYPE denominator = std::pow(2.0, 2 * l + 2 * k + 2)
                                      * ((FP_TYPE) ((2 * l + 3) * (2 * k + 3)))
                                      * factorial_array[l] * factorial_array[l + 1]
                                      * factorial_array[k] * factorial_array[k + 1]
                                      * factorial_array[n + 1]
                                      * std::pow((1.0 + (data.L_1 + d) / Z) * (1.0 + (data.L_1 + data.L_2 + d) / Z),
                                                 2 * l + 2 * k + n + 2);

                FP_TYPE term = numerator / denominator;

                if (!std::isfinite(term)) break;

                M_12 += term;
            }
        }
    }

    M_12 *= (FP_TYPE) (4e-7) * local_pi * local_pi * data.N_1 * data.N_2 * std::pow(Z, 5)
            / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    if (timing) {
        double interval = std::chrono::duration_cast<std::chrono::duration<double>>(
                std::chrono::high_resolution_clock::now() - begin_time).count();
        std::cout << "Unoptimized Time = " << interval << " s" << std::endl;
    }

    return M_12;
}

double calculate_mutual_inductance_near_unoptimized(
        const CoilCalculationData& data,
        const SumPrecisionData &precision,
        const FP_TYPE d,
        const FP_TYPE Z,
        bool timing = false
) {
    FP_TYPE M_12 = 0.0;

    std::chrono::high_resolution_clock::time_point begin_time;

    if (timing) {
        begin_time = std::chrono::high_resolution_clock::now();
    }

    for (uint32_t n = 0; n < precision.n_terms; ++n) {
        for (uint32_t l = 0; l < precision.l_terms; ++l) {
            for (uint32_t k = 0; k < precision.k_terms; ++k) {

                FP_TYPE sign = (l + k) % 2 == 0 ? 1.0 : -1.0;

                FP_TYPE numerator = sign
                                    * factorial_array[2 * l + 2 * k + n + 1]
                                    * (std::pow(Z + data.L_1, n + 1) - std::pow(Z, n + 1))
                                    * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                                    * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                                    * (1.0 / std::pow(Z + data.L_1 + d, 2 * l + 2 * k + n + 2)
                                       - 1.0 / std::pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * k + n + 2)
                                    );

                FP_TYPE denominator = std::pow(2.0, 2 * l + 2 * k + 2)
                                      * ((FP_TYPE) ((2 * l + 3) * (2 * k + 3)))
                                      * factorial_array[l] * factorial_array[l + 1]
                                      * factorial_array[k] * factorial_array[k + 1]
                                      * factorial_array[n + 1];

                FP_TYPE term = numerator / denominator;

                if (!std::isfinite(term)) continue;

                M_12 += term;
            }
        }
    }

    M_12 *= 4.0 * local_pi * local_pi * 1e-7 * data.N_1 * data.N_2
            / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    if (timing) {
        double interval = std::chrono::duration_cast<std::chrono::duration<double>>(
                std::chrono::high_resolution_clock::now() - begin_time).count();
        std::cout << "Unoptimized Time = " << interval << " s" << std::endl;
    }

    return M_12;
}

#endif //VECTOR_CASE_INDUCTANCE_NAIVE_HPP