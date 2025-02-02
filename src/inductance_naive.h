#ifndef VECTOR_CASE_INDUCTANCE_NAIVE_H
#define VECTOR_CASE_INDUCTANCE_NAIVE_H

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#include "structs.h"
#include "factorial_lookup.h"


/**
 * @brief Naive implementation of the mutual inductance calculation for coils that are far apart.
 *
 * This function is meant to show how slow C code can be when not optimized for the hardware nor
 * algorithmically (using pow). It is not recommended to use this function for any practical
 * application.
 *
 * @param data The coil calculation data containing the physical properties of the coils.
 * @param precision The precision data specifying the number of terms in the series expansion.
 * @param d The distance between the coils.
 * @param timing Whether to print the time taken for the calculation.
 * @return
 */
double calculate_mutual_inductance_far_unoptimized(
        const CoilCalculationData data,
        const SumPrecisionData precision,
        const double d,
        const bool timing
) {
    double M_12 = 0.0;

    struct timespec start_time;

    if (timing) {
        timespec_get(&start_time, TIME_UTC);
    }

    for (uint32_t n = 0; n < precision.n_terms; ++n) {
        for (uint32_t l = 0; l < precision.l_terms; ++l) {
            for (uint32_t k = 0; k < precision.k_terms; ++k) {

                double sign = (l + k) % 2 == 0 ? 1.0 : -1.0;

                double numerator = sign
                                 * factorial_array[2 * n + 2 * l + 2 * k + 1]
                                 * pow(data.L_1, 2 * n)
                                 * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                                 * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                                 * (1.0 / pow(0.5 * data.L_1 + d, 2 * n + 2 * l + 2 * k + 2)
                                    - 1.0 / pow(0.5 * data.L_1 + d + data.L_2, 2 * n + 2 * l + 2 * k + 2)
                                   );

                double denominator = pow(2.0, 2 * n + 2 * l + 2 * k + 2)
                                   * ((double) ((2 * l + 3) * (2 * k + 3)))
                                   * factorial_array[l] * factorial_array[l + 1]
                                   * factorial_array[k] * factorial_array[k + 1]
                                   * factorial_array[2 * n + 1];

                double term = numerator / denominator;

                M_12 += term;
            }
        }
    }

    M_12 *= 4.0e-7 * local_pi * local_pi * data.N_1 * data.N_2
          / (data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    if (timing) {
        struct timespec end_time;
        timespec_get(&end_time, TIME_UTC);

        double interval = (double) (end_time.tv_sec - start_time.tv_sec)
                        + (double) (end_time.tv_nsec - start_time.tv_nsec) * 1e-9;
        printf("Unoptimized Time = %g s\n", interval);
    }

    return M_12;
}


/**
 * @brief Naive implementation of the mutual inductance calculation for coils that are near each other.
 * Normed version.
 *
 * This function is meant to show how slow C code can be when not optimized for the hardware nor
 * algorithmically (using pow). It is not recommended to use this function for any practical
 * application.
 *
 * The normed version was introduced hoping to avoid numerical issues, but in the end it turned
 * out that there are really no numerical issues.
 *
 * @param data The coil calculation data containing the physical properties of the coils.
 * @param precision The precision data specifying the number of terms in the series expansion.
 * @param d The distance between the coils.
 * @param Z The distance between the coils and the point of interest.
 * @param timing Whether to print the time taken for the calculation.
 * @return
 */
double calculate_mutual_inductance_near_unoptimized_normed(
        const CoilCalculationData data,
        const SumPrecisionData precision,
        const double d,
        const double Z,
        const bool timing
) {
    double M_12 = 0.0;

    struct timespec start_time;

    if (timing) {
        timespec_get(&start_time, TIME_UTC);
    }

    for (uint32_t n = 0; n < precision.n_terms; ++n) {
        for (uint32_t l = 0; l < precision.l_terms; ++l) {
            for (uint32_t k = 0; k < precision.k_terms; ++k) {

                double sign = (l + k) % 2 == 0 ? 1.0 : -1.0;

                double numerator = sign
                                    * factorial_array[2 * l + 2 * k + n + 1]
                                    * (pow(1.0 + data.L_1 / Z, n + 1) - 1.0)
                                    * (pow(data.R_1 / Z, 2 * l + 3) - pow(data.r_1 / Z, 2 * l + 3))
                                    * (pow(data.R_2 / Z, 2 * k + 3) - pow(data.r_2 / Z, 2 * k + 3))
                                    * (pow(1.0 + (data.L_1 + data.L_2 + d) / Z, 2 * l + 2 * k + n + 2)
                                       - pow(1.0 + (data.L_1 + d) / Z, 2 * l + 2 * k + n + 2)
                                    );

                double denominator = pow(2.0, 2 * l + 2 * k + 2)
                                      * ((double) ((2 * l + 3) * (2 * k + 3)))
                                      * factorial_array[l] * factorial_array[l + 1]
                                      * factorial_array[k] * factorial_array[k + 1]
                                      * factorial_array[n + 1]
                                      * pow((1.0 + (data.L_1 + d) / Z) * (1.0 + (data.L_1 + data.L_2 + d) / Z),
                                                 2 * l + 2 * k + n + 2);

                double term = numerator / denominator;

                M_12 += term;
            }
        }
    }

    M_12 *= (double) (4e-7) * local_pi * local_pi * data.N_1 * data.N_2 * pow(Z, 5)
          / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    if (timing) {
        struct timespec end_time;
        timespec_get(&end_time, TIME_UTC);

        double interval = (double) (end_time.tv_sec - start_time.tv_sec)
                          + (double) (end_time.tv_nsec - start_time.tv_nsec) * 1e-9;
        printf("Unoptimized Time = %g s\n", interval);
    }

    return M_12;
}


/**
 * @brief Naive implementation of the mutual inductance calculation for coils that are near each other.
 *
 * This function is meant to show how slow C code can be when not optimized for the hardware nor
 * algorithmically (using pow). It is not recommended to use this function for any practical
 * application.
 *
 * @param data The coil calculation data containing the physical properties of the coils.
 * @param precision The precision data specifying the number of terms in the series expansion.
 * @param d The distance between the coils.
 * @param Z The value of the free parameter, modulates the convergence of the series.
 * @param timing Whether to print the time taken for the calculation.
 * @return
 */
double calculate_mutual_inductance_near_unoptimized(
        const CoilCalculationData data,
        const SumPrecisionData precision,
        const double d,
        const double Z,
        const bool timing
) {
    double M_12 = 0.0;

    struct timespec start_time;

    if (timing) {
        timespec_get(&start_time, TIME_UTC);
    }

    for (uint32_t n = 0; n < precision.n_terms; ++n) {
        for (uint32_t l = 0; l < precision.l_terms; ++l) {
            for (uint32_t k = 0; k < precision.k_terms; ++k) {

                double sign = (l + k) % 2 == 0 ? 1.0 : -1.0;

                double numerator = sign
                                    * factorial_array[2 * l + 2 * k + n + 1]
                                    * (pow(Z + data.L_1, n + 1) - pow(Z, n + 1))
                                    * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                                    * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                                    * (1.0 / pow(Z + data.L_1 + d, 2 * l + 2 * k + n + 2)
                                       - 1.0 / pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * k + n + 2)
                                    );

                double denominator = pow(2.0, 2 * l + 2 * k + 2)
                                      * ((double) ((2 * l + 3) * (2 * k + 3)))
                                      * factorial_array[l] * factorial_array[l + 1]
                                      * factorial_array[k] * factorial_array[k + 1]
                                      * factorial_array[n + 1];

                double term = numerator / denominator;

                M_12 += term;
            }
        }
    }

    M_12 *= 4.0 * local_pi * local_pi * 1e-7 * data.N_1 * data.N_2
          / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    if (timing) {
        struct timespec end_time;
        timespec_get(&end_time, TIME_UTC);

        double interval = (double) (end_time.tv_sec - start_time.tv_sec)
                          + (double) (end_time.tv_nsec - start_time.tv_nsec) * 1e-9;
        printf("Unoptimized Time = %g s\n", interval);
    }

    return M_12;
}

#endif //VECTOR_CASE_INDUCTANCE_NAIVE_H