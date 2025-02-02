#ifndef VECTOR_CASE_INDUCTANCE_REMAINDER_BOUND_H
#define VECTOR_CASE_INDUCTANCE_REMAINDER_BOUND_H

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "structs.h"
#include "factorial_lookup.h"


/**
 * @brief Remainder bound for the triple series expansion of coaxial mutual inductance. This version
 * is unoptimized but the more optimized one is not yet done.
 *
 * The paper outlined the estimate of the remainder term for the triple series expansion in the form
 * of a bound. This function calculates the bound for the remainder for given parameters. The bound is
 * usually quite pessimistic, but it can be a good indicator of the error if the coils are sufficiently
 * separated.
 *
* @param data The coil calculation data containing the physical properties of the coils.
 * @param precision The precision data specifying the number of terms in the series expansion.
 * @param d The distance between the coils.
 * @param Z The value of free parameter, modulates the convergence of the series.
 * @param timing Whether to print the time taken for the calculation.
 * @return The calculated remainder bound.
 */
FP_TYPE calculate_inductance_remainder_unoptimized(
        const CoilCalculationData data,
        const SumPrecisionData precision,
        const FP_TYPE d,
        const FP_TYPE Z,
        bool timing
) {
    struct timespec start_time;

    if (timing) {
        timespec_get(&start_time, TIME_UTC);
    }

    FP_TYPE abs_Z = fabs(Z);

    int L = (int) precision.l_terms - 1;
    int K = (int) precision.k_terms - 1;
    int N = (int) precision.n_terms - 1;

    FP_TYPE u_1 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int k = 0; k <= K; ++k) {
            FP_TYPE num = factorial_array[2 * l + 2 * k + N + 2]
                          * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                          * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                          * pow(abs_Z, N + 2)
                          * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * l + 2 * k + N + 3)
                             - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * k + N + 3));

            FP_TYPE den = factorial_array[N + 2]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[l] * factorial_array[l + 1]
                          * pow(2.0, 2 * l + 2 * k + 2)
                          * (FP_TYPE)((2 * k + 3) * (2 * l + 3));

            u_1 += num / den;
        }
    }

    FP_TYPE u_2 = 0.0;

    for (int k = 0; k <= K; ++k) {
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * L + 2 * k + n + 3]
                          * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                          * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                          * pow(abs_Z, n + 1)
                          * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * k + n + 4)
                             - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[2 * L + 3]
                          * pow(2.0, 2 * k + 1)
                          * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

            u_2 += num / den;
        }
    }

    FP_TYPE u_3 = 0.0;

    for (int l = 0; l <= L; ++l) {
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * l + 2 * K + n + 3]
                          * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                          * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                          * pow(abs_Z, n + 1)
                          * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * l + 2 * K + n + 4)
                             - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[2 * K + 3]
                          * factorial_array[l] * factorial_array[l + 1]
                          * pow(2.0, 2 * l + 1)
                          * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

            u_3 += num / den;
        }
    }

    FP_TYPE u_4 = 0.0;

    for (int l = 0; l <= L; ++l) {
        FP_TYPE num = factorial_array[2 * l + 2 * K + N + 4]
                      * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                      * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                      * pow(abs_Z, N + 2)
                      * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * l + 2 * K + N + 5)
                         - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + N + 5));

        FP_TYPE den = factorial_array[N + 2]
                      * factorial_array[2 * K + 3]
                      * factorial_array[l] * factorial_array[l + 1]
                      * pow(2.0, 2 * l + 1)
                      * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

        u_4 += num / den;
    }

    FP_TYPE u_5 = 0.0;

    for (int k = 0; k <= K; ++k) {
        FP_TYPE num = factorial_array[2 * L + 2 * k + N + 4]
                      * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                      * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                      * pow(abs_Z, N + 2)
                      * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * k + N + 5)
                         - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + N + 5));

        FP_TYPE den = factorial_array[N + 2]
                      * factorial_array[k] * factorial_array[k + 1]
                      * factorial_array[2 * L + 3]
                      * pow(2.0, 2 * k + 1)
                      * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

        u_5 += num / den;
    }

    FP_TYPE u_6 = 0.0;

    for (int n = 0; n <= N; ++n) {
        FP_TYPE num = factorial_array[2 * L + 2 * K + n + 5]
                      * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                      * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                      * pow(abs_Z, n + 1)
                      * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * K + n + 6)
                         - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + n + 6));

        FP_TYPE den = factorial_array[n + 1]
                      * factorial_array[2 * K + 3]
                      * factorial_array[2 * L + 3]
                      * (FP_TYPE)((2 * K + 5) * (2 * L + 5));

        u_6 += num / den;
    }

    FP_TYPE u_7 = (factorial_array[2 * L + 2 * K + N + 6]
                   * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                   * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                   * pow(abs_Z, N + 2)
                   * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * K + N + 7)
                      - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + N + 7))
                  ) / (factorial_array[N + 2]
                       * factorial_array[2 * K + 3]
                       * factorial_array[2 * L + 3]
                       * (FP_TYPE)((2 * K + 5) * (2 * L + 5)));

    FP_TYPE u = u_1 + u_2 + u_3 + u_4 + u_5 + u_6 + u_7;

    FP_TYPE q_1 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int k = 0; k <= K; ++k) {
            FP_TYPE num = exp(lgamma((FP_TYPE) (2 * l + 2 * k + N) + (FP_TYPE) 2.5))
                          * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                          * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                          * pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                          * ((FP_TYPE) (1.0) / pow(d, (FP_TYPE) (2 * l + 2 * k + N) + 2.5)
                             - (FP_TYPE) (1.0) / pow(data.L_2 + d, (FP_TYPE) (2 * l + 2 * k + N) + 2.5));

            FP_TYPE den = factorial_array[N + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[l] * factorial_array[l + 1]
                          * pow(2.0, 2 * l + 2 * k + 2)
                          * (FP_TYPE)((2 * k + 3) * (2 * l + 3))
                          * sqrt((FP_TYPE) (2 * (2 * N + 3)));

            q_1 += num / den;
        }
    }

    FP_TYPE q_2 = 0.0;

    for (int k = 0; k <= K; ++k){
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * L + 2 * k + n + 3]
                          * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                          * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                          * pow(Z + data.L_1, n + 1)
                          * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * k + n + 4)
                             - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[2 * L + 3]
                          * pow(2.0, 2 * k + 1)
                          * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

            q_2 += num / den;
        }
    }

    FP_TYPE q_3 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * l + 2 * K + n + 3]
                          * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                          * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                          * pow(Z + data.L_1, n + 1)
                          * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * l + 2 * K + n + 4)
                             - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[2 * K + 3]
                          * factorial_array[l] * factorial_array[l + 1]
                          * pow(2.0, 2 * l + 1)
                          * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

            q_3 += num / den;
        }
    }

    FP_TYPE q_4 = 0.0;


    for (int k = 0; k <= K; ++k) {
        FP_TYPE num = exp(lgamma((FP_TYPE) (2 * L + 2 * k + N) + (FP_TYPE) 4.5))
                      * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                      * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                      * pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                      * ((FP_TYPE) (1.0) / pow(d, (FP_TYPE) (2 * L + 2 * k + N) + 4.5)
                         - (FP_TYPE) (1.0) / pow(data.L_2 + d, (FP_TYPE) (2 * L + 2 * k + N) + 4.5));

        FP_TYPE den = factorial_array[N + 1]
                      * factorial_array[k] * factorial_array[k + 1]
                      * factorial_array[2 * L + 3]
                      * pow(2.0, 2 * k + 1)
                      * (FP_TYPE)((2 * k + 3) * (2 * L + 5))
                      * sqrt((FP_TYPE) (2 * (2 * N + 3)));

        q_4 += num / den;
    }

    FP_TYPE q_5 = 0.0;

    for (int l = 0; l <= L; ++l) {
        FP_TYPE num = exp(lgamma((FP_TYPE) (2 * l + 2 * K + N) + (FP_TYPE) 4.5))
                      * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                      * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                      * pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                      * ((FP_TYPE) (1.0) / pow(d, (FP_TYPE) (2 * l + 2 * K + N) + 4.5)
                         - (FP_TYPE) (1.0) / pow(data.L_2 + d, (FP_TYPE) (2 * l + 2 * K + N) + 4.5));

        FP_TYPE den = factorial_array[N + 1]
                      * factorial_array[l] * factorial_array[l + 1]
                      * factorial_array[2 * K + 3]
                      * pow(2.0, 2 * l + 1)
                      * (FP_TYPE)((2 * K + 5) * (2 * l + 3))
                      * sqrt((FP_TYPE) (2 * (2 * N + 3)));

        q_5 += num / den;
    }

    FP_TYPE q_6 = 0.0;

    for (int n = 0; n <= N; ++n) {
        FP_TYPE num = factorial_array[2 * L + 2 * K + n + 5]
                      * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                      * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                      * pow(Z + data.L_1, n + 1)
                      * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * K + n + 6)
                         - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + n + 6));

        FP_TYPE den = factorial_array[n + 1]
                      * factorial_array[2 * K + 3]
                      * factorial_array[2 * L + 3]
                      * (FP_TYPE)((2 * K + 5) * (2 * L + 5));

        q_6 += num / den;
    }

    FP_TYPE q_7 = (exp(lgamma((FP_TYPE) (2 * L + 2 * K + N) + (FP_TYPE) 6.5))
                   * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                   * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                   * pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                   * ((FP_TYPE) (1.0) / pow(d, (FP_TYPE) (2 * L + 2 * K + N) + 6.5)
                      - (FP_TYPE) (1.0) / pow(data.L_2 + d, (FP_TYPE) (2 * L + 2 * K + N) + 6.5))
                  ) / (factorial_array[N + 1]
                       * factorial_array[2 * K + 3]
                       * factorial_array[2 * L + 3]
                       * (FP_TYPE)((2 * K + 5) * (2 *L + 5))
                       * sqrt((FP_TYPE) (2 * (2 * N + 3))));

    FP_TYPE q = q_1 + q_2 + q_3 + q_4 + q_5 + q_6 + q_7;

    FP_TYPE r = q + u;

    r *= (FP_TYPE) (4.0e-7) * local_pi * local_pi * data.N_1 * data.N_2
         / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    if (timing) {
        struct timespec end_time;
        timespec_get(&end_time, TIME_UTC);

        double interval = (double) (end_time.tv_sec - start_time.tv_sec)
                        + (double) (end_time.tv_nsec - start_time.tv_nsec) * 1e-9;
        printf("Remainder Time = %f s\n", interval);
    }

    return r;
}


/**
 * Work in progress: Optimized version of the remainder bound calculation for the triple series expansion
 *
 * @param data The coil calculation data containing the physical properties of the coils
 * @param precision The precision data specifying the number of terms in the series expansion
 * @param d The distance between the coils
 * @param Z The value of the free parameter, modulates the convergence of the series
 * @param timing Whether to print the time taken for the calculation
 * @return The calculated remainder bound
 */
FP_TYPE calculate_inductance_remainder(
        const CoilCalculationData data,
        const SumPrecisionData precision,
        const FP_TYPE d,
        const FP_TYPE Z,
        bool timing
) {
    struct timespec start_time;

    if (timing) {
        timespec_get(&start_time, TIME_UTC);
    }

    int L = (int) precision.l_terms - 1;
    int K = (int) precision.k_terms - 1;
    int N = (int) precision.n_terms - 1;

    FP_TYPE R_1_sub_r_1_arr[MAX_TERMS_NEAR + 1];
    FP_TYPE R_2_sub_r_2_arr[MAX_TERMS_NEAR + 1];
    FP_TYPE abs_Z_arr[MAX_TERMS_NEAR + 1];
    FP_TYPE L_1_plus_Z_arr[MAX_TERMS_NEAR + 1];

    FP_TYPE denom_first_1_arr[MAX_TERMS_NEAR + 1];
    FP_TYPE denom_first_2_arr[MAX_TERMS_NEAR + 1];
    FP_TYPE denom_second_1_arr[MAX_TERMS_NEAR + 1];
    FP_TYPE denom_second_2_arr[MAX_TERMS_NEAR + 1];

    FP_TYPE R_1_sq = data.R_1 * data.R_1;
    FP_TYPE r_1_sq = data.r_1 * data.r_1;

    FP_TYPE R_2_sq = data.R_2 * data.R_2;
    FP_TYPE r_2_sq = data.r_2 * data.r_2;

    FP_TYPE denom_first_1 = (FP_TYPE) 1.0 / (Z + data.L_1 + d);
    FP_TYPE denom_first_2 = (FP_TYPE) 1.0 / (Z + data.L_1 + data.L_2 + d);
    FP_TYPE denom_second_1 = (FP_TYPE) 1.0 / (d);
    FP_TYPE denom_second_2 = (FP_TYPE) 1.0 / (data.L_2 + d);

    FP_TYPE denom_first_1_sq = denom_first_1 * denom_first_1;
    FP_TYPE denom_first_2_sq = denom_first_2 * denom_first_2;
    FP_TYPE denom_second_1_sq = denom_second_1 * denom_second_1;
    FP_TYPE denom_second_2_sq = denom_second_2 * denom_second_2;

    FP_TYPE abs_Z = fabs(Z);
    FP_TYPE L_1_plus_Z = data.L_1 + Z;

    FP_TYPE temp_R_1 = R_1_sq * data.R_1;
    FP_TYPE temp_r_1 = r_1_sq * data.r_1;

    for (int l = 0; l <= precision.l_terms; ++l) {
        R_1_sub_r_1_arr[l] = temp_R_1 - temp_r_1;
        temp_R_1 *= R_1_sq;
        temp_r_1 *= r_1_sq;
    }

    FP_TYPE temp_R_2 = R_2_sq * data.R_2;
    FP_TYPE temp_r_2 = r_2_sq * data.r_2;

    for (int k = 0; k <= precision.k_terms; ++k) {
        R_2_sub_r_2_arr[k] = temp_R_2 - temp_r_2;
        temp_R_2 *= R_2_sq;
        temp_r_2 *= r_2_sq;
    }

    FP_TYPE temp_abs_Z = abs_Z;
    FP_TYPE temp_L_1_plus_Z = L_1_plus_Z;

    for (int n = 0; n <= precision.n_terms; ++n) {
        abs_Z_arr[n] = temp_abs_Z;
        L_1_plus_Z_arr[n] = temp_L_1_plus_Z;

        temp_abs_Z *= abs_Z;
        temp_L_1_plus_Z *= L_1_plus_Z;
    }

    FP_TYPE u_1 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int k = 0; k <= K; ++k) {
            FP_TYPE num = factorial_array[2 * l + 2 * k + N + 2]
                          * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                          * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                          * pow(abs_Z, N + 2)
                          * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * l + 2 * k + N + 3)
                             - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * k + N + 3));

            FP_TYPE den = factorial_array[N + 2]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[l] * factorial_array[l + 1]
                          * pow(2.0, 2 * l + 2 * k + 2)
                          * (FP_TYPE)((2 * k + 3) * (2 * l + 3));

            u_1 += num / den;
        }
    }

    FP_TYPE u_2 = 0.0;

    for (int k = 0; k <= K; ++k) {
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * L + 2 * k + n + 3]
                          * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                          * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                          * pow(abs_Z, n + 1)
                          * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * k + n + 4)
                             - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[2 * L + 3]
                          * pow(2.0, 2 * k + 1)
                          * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

            u_2 += num / den;
        }
    }

    FP_TYPE u_3 = 0.0;

    for (int l = 0; l <= L; ++l) {
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * l + 2 * K + n + 3]
                          * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                          * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                          * pow(abs_Z, n + 1)
                          * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * l + 2 * K + n + 4)
                             - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[2 * K + 3]
                          * factorial_array[l] * factorial_array[l + 1]
                          * pow(2.0, 2 * l + 1)
                          * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

            u_3 += num / den;
        }
    }

    FP_TYPE u_4 = 0.0;

    for (int l = 0; l <= L; ++l) {
        FP_TYPE num = factorial_array[2 * l + 2 * K + N + 4]
                      * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                      * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                      * pow(abs_Z, N + 2)
                      * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * l + 2 * K + N + 5)
                         - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + N + 5));

        FP_TYPE den = factorial_array[N + 2]
                      * factorial_array[2 * K + 3]
                      * factorial_array[l] * factorial_array[l + 1]
                      * pow(2.0, 2 * l + 1)
                      * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

        u_4 += num / den;
    }

    FP_TYPE u_5 = 0.0;

    for (int k = 0; k <= K; ++k) {
        FP_TYPE num = factorial_array[2 * L + 2 * k + N + 4]
                      * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                      * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                      * pow(abs_Z, N + 2)
                      * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * k + N + 5)
                         - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + N + 5));

        FP_TYPE den = factorial_array[N + 2]
                      * factorial_array[k] * factorial_array[k + 1]
                      * factorial_array[2 * L + 3]
                      * pow(2.0, 2 * k + 1)
                      * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

        u_5 += num / den;
    }

    FP_TYPE u_6 = 0.0;

    for (int n = 0; n <= N; ++n) {
        FP_TYPE num = factorial_array[2 * L + 2 * K + n + 5]
                      * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                      * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                      * pow(abs_Z, n + 1)
                      * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * K + n + 6)
                         - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + n + 6));

        FP_TYPE den = factorial_array[n + 1]
                      * factorial_array[2 * K + 3]
                      * factorial_array[2 * L + 3]
                      * (FP_TYPE)((2 * K + 5) * (2 * L + 5));

        u_6 += num / den;
    }

    FP_TYPE u_7 = (factorial_array[2 * L + 2 * K + N + 6]
                   * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                   * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                   * pow(abs_Z, N + 2)
                   * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * K + N + 7)
                      - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + N + 7))
                  ) / (factorial_array[N + 2]
                       * factorial_array[2 * K + 3]
                       * factorial_array[2 * L + 3]
                       * (FP_TYPE)((2 * K + 5) * (2 * L + 5)));

    FP_TYPE u = u_1 + u_2 + u_3 + u_4 + u_5 + u_6 + u_7;

    FP_TYPE q_1 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int k = 0; k <= K; ++k) {
            FP_TYPE num = exp(lgamma((FP_TYPE) (2 * l + 2 * k + N) + (FP_TYPE) 2.5))
                          * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                          * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                          * pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                          * ((FP_TYPE) (1.0) / pow(d, (FP_TYPE) (2 * l + 2 * k + N) + 2.5)
                             - (FP_TYPE) (1.0) / pow(data.L_2 + d, (FP_TYPE) (2 * l + 2 * k + N) + 2.5));

            FP_TYPE den = factorial_array[N + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[l] * factorial_array[l + 1]
                          * pow(2.0, 2 * l + 2 * k + 2)
                          * (FP_TYPE)((2 * k + 3) * (2 * l + 3))
                          * sqrt((FP_TYPE) (2 * (2 * N + 3)));

            q_1 += num / den;
        }
    }

    FP_TYPE q_2 = 0.0;

    for (int k = 0; k <= K; ++k){
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * L + 2 * k + n + 3]
                          * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                          * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                          * pow(Z + data.L_1, n + 1)
                          * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * k + n + 4)
                             - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[2 * L + 3]
                          * pow(2.0, 2 * k + 1)
                          * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

            q_2 += num / den;
        }
    }

    FP_TYPE q_3 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * l + 2 * K + n + 3]
                          * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                          * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                          * pow(Z + data.L_1, n + 1)
                          * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * l + 2 * K + n + 4)
                             - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[2 * K + 3]
                          * factorial_array[l] * factorial_array[l + 1]
                          * pow(2.0, 2 * l + 1)
                          * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

            q_3 += num / den;
        }
    }

    FP_TYPE q_4 = 0.0;


    for (int k = 0; k <= K; ++k) {
        FP_TYPE num = exp(lgamma((FP_TYPE) (2 * L + 2 * k + N) + (FP_TYPE) 4.5))
                      * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                      * (pow(data.R_2, 2 * k + 3) - pow(data.r_2, 2 * k + 3))
                      * pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                      * ((FP_TYPE) (1.0) / pow(d, (FP_TYPE) (2 * L + 2 * k + N) + 4.5)
                         - (FP_TYPE) (1.0) / pow(data.L_2 + d, (FP_TYPE) (2 * L + 2 * k + N) + 4.5));

        FP_TYPE den = factorial_array[N + 1]
                      * factorial_array[k] * factorial_array[k + 1]
                      * factorial_array[2 * L + 3]
                      * pow(2.0, 2 * k + 1)
                      * (FP_TYPE)((2 * k + 3) * (2 * L + 5))
                      * sqrt((FP_TYPE) (2 * (2 * N + 3)));

        q_4 += num / den;
    }

    FP_TYPE q_5 = 0.0;

    for (int l = 0; l <= L; ++l) {
        FP_TYPE num = exp(lgamma((FP_TYPE) (2 * l + 2 * K + N) + (FP_TYPE) 4.5))
                      * (pow(data.R_1, 2 * l + 3) - pow(data.r_1, 2 * l + 3))
                      * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                      * pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                      * ((FP_TYPE) (1.0) / pow(d, (FP_TYPE) (2 * l + 2 * K + N) + 4.5)
                         - (FP_TYPE) (1.0) / pow(data.L_2 + d, (FP_TYPE) (2 * l + 2 * K + N) + 4.5));

        FP_TYPE den = factorial_array[N + 1]
                      * factorial_array[l] * factorial_array[l + 1]
                      * factorial_array[2 * K + 3]
                      * pow(2.0, 2 * l + 1)
                      * (FP_TYPE)((2 * K + 5) * (2 * l + 3))
                      * sqrt((FP_TYPE) (2 * (2 * N + 3)));

        q_5 += num / den;
    }

    FP_TYPE q_6 = 0.0;

    for (int n = 0; n <= N; ++n) {
        FP_TYPE num = factorial_array[2 * L + 2 * K + n + 5]
                      * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                      * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                      * pow(Z + data.L_1, n + 1)
                      * ((FP_TYPE) (1.0) / pow(Z + data.L_1 + d, 2 * L + 2 * K + n + 6)
                         - (FP_TYPE) (1.0) / pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + n + 6));

        FP_TYPE den = factorial_array[n + 1]
                      * factorial_array[2 * K + 3]
                      * factorial_array[2 * L + 3]
                      * (FP_TYPE)((2 * K + 5) * (2 * L + 5));

        q_6 += num / den;
    }

    FP_TYPE q_7 = (exp(lgamma((FP_TYPE) (2 * L + 2 * K + N) + (FP_TYPE) 6.5))
                   * (pow(data.R_1, 2 * L + 5) - pow(data.r_1, 2 * L + 5))
                   * (pow(data.R_2, 2 * K + 5) - pow(data.r_2, 2 * K + 5))
                   * pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                   * ((FP_TYPE) (1.0) / pow(d, (FP_TYPE) (2 * L + 2 * K + N) + 6.5)
                      - (FP_TYPE) (1.0) / pow(data.L_2 + d, (FP_TYPE) (2 * L + 2 * K + N) + 6.5))
                  ) / (factorial_array[N + 1]
                       * factorial_array[2 * K + 3]
                       * factorial_array[2 * L + 3]
                       * (FP_TYPE)((2 * K + 5) * (2 *L + 5))
                       * sqrt((FP_TYPE) (2 * (2 * N + 3))));

    FP_TYPE q = q_1 + q_2 + q_3 + q_4 + q_5 + q_6 + q_7;

    FP_TYPE r = q + u;

    r *= (FP_TYPE) (4.0e-7) * local_pi * local_pi * data.N_1 * data.N_2
         / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    return r;
}

#endif //VECTOR_CASE_INDUCTANCE_REMAINDER_BOUND_H
