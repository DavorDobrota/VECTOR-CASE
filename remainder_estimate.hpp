#ifndef VECTOR_CASE_REMAINDER_ESTIMATE_HPP
#define VECTOR_CASE_REMAINDER_ESTIMATE_HPP

#include <cmath>
#include <chrono>

#include "structs.h"
#include "factorial_lookup.h"


FP_TYPE calculate_inductance_remainder_unoptimized(
        const CoilCalculationData& data,
        const SumPrecisionData &precision,
        const FP_TYPE d,
        const FP_TYPE Z,
        bool timing = false
) {
    std::chrono::high_resolution_clock::time_point begin_time;

    if (timing) {
        begin_time = std::chrono::high_resolution_clock::now();
    }

    FP_TYPE abs_Z = std::abs(Z);

    auto L = (int) precision.l_terms - 1;
    auto K = (int) precision.k_terms - 1;
    auto N = (int) precision.n_terms - 1;

    FP_TYPE u_1 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int k = 0; k <= K; ++k) {
            FP_TYPE num = factorial_array[2 * l + 2 * k + N + 2]
                          * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                          * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                          * std::pow(abs_Z, N + 2)
                          * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * l + 2 * k + N + 3)
                             - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * k + N + 3));

            FP_TYPE den = factorial_array[N + 2]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[l] * factorial_array[l + 1]
                          * std::pow(2.0, 2 * l + 2 * k + 2)
                          * (FP_TYPE)((2 * k + 3) * (2 * l + 3));

            u_1 += num / den;
        }
    }

    FP_TYPE u_2 = 0.0;

    for (int k = 0; k <= K; ++k) {
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * L + 2 * k + n + 3]
                          * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                          * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                          * std::pow(abs_Z, n + 1)
                          * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * k + n + 4)
                             - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[2 * L + 3]
                          * std::pow(2.0, 2 * k + 1)
                          * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

            u_2 += num / den;
        }
    }

    FP_TYPE u_3 = 0.0;

    for (int l = 0; l <= L; ++l) {
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * l + 2 * K + n + 3]
                          * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                          * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                          * std::pow(abs_Z, n + 1)
                          * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * l + 2 * K + n + 4)
                             - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[2 * K + 3]
                          * factorial_array[l] * factorial_array[l + 1]
                          * std::pow(2.0, 2 * l + 1)
                          * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

            u_3 += num / den;
        }
    }

    FP_TYPE u_4 = 0.0;

    for (int l = 0; l <= L; ++l) {
        FP_TYPE num = factorial_array[2 * l + 2 * K + N + 4]
                      * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                      * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                      * std::pow(abs_Z, N + 2)
                      * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * l + 2 * K + N + 5)
                         - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + N + 5));

        FP_TYPE den = factorial_array[N + 2]
                      * factorial_array[2 * K + 3]
                      * factorial_array[l] * factorial_array[l + 1]
                      * std::pow(2.0, 2 * l + 1)
                      * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

        u_4 += num / den;
    }

    FP_TYPE u_5 = 0.0;

    for (int k = 0; k <= K; ++k) {
        FP_TYPE num = factorial_array[2 * L + 2 * k + N + 4]
                      * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                      * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                      * std::pow(abs_Z, N + 2)
                      * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * k + N + 5)
                         - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + N + 5));

        FP_TYPE den = factorial_array[N + 2]
                      * factorial_array[k] * factorial_array[k + 1]
                      * factorial_array[2 * L + 3]
                      * std::pow(2.0, 2 * k + 1)
                      * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

        u_5 += num / den;
    }

    FP_TYPE u_6 = 0.0;

    for (int n = 0; n <= N; ++n) {
        FP_TYPE num = factorial_array[2 * L + 2 * K + n + 5]
                      * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                      * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                      * std::pow(abs_Z, n + 1)
                      * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * K + n + 6)
                         - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + n + 6));

        FP_TYPE den = factorial_array[n + 1]
                      * factorial_array[2 * K + 3]
                      * factorial_array[2 * L + 3]
                      * (FP_TYPE)((2 * K + 5) * (2 * L + 5));

        u_6 += num / den;
    }

    FP_TYPE u_7 = (factorial_array[2 * L + 2 * K + N + 6]
                   * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                   * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                   * std::pow(abs_Z, N + 2)
                   * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * K + N + 7)
                      - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + N + 7))
                  ) / (factorial_array[N + 2]
                       * factorial_array[2 * K + 3]
                       * factorial_array[2 * L + 3]
                       * (FP_TYPE)((2 * K + 5) * (2 * L + 5)));

    FP_TYPE u = u_1 + u_2 + u_3 + u_4 + u_5 + u_6 + u_7;

    FP_TYPE q_1 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int k = 0; k <= K; ++k) {
            FP_TYPE num = std::exp(std::lgamma((FP_TYPE) (2 * l + 2 * k + N) + (FP_TYPE) 2.5))
                          * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                          * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                          * std::pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                          * ((FP_TYPE) (1.0) / std::pow(d, (FP_TYPE) (2 * l + 2 * k + N) + 2.5)
                             - (FP_TYPE) (1.0) / std::pow(data.L_2 + d, (FP_TYPE) (2 * l + 2 * k + N) + 2.5));

            FP_TYPE den = factorial_array[N + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[l] * factorial_array[l + 1]
                          * std::pow(2.0, 2 * l + 2 * k + 2)
                          * (FP_TYPE)((2 * k + 3) * (2 * l + 3))
                          * std::sqrt((FP_TYPE) (2 * (2 * N + 3)));

            q_1 += num / den;
        }
    }

    FP_TYPE q_2 = 0.0;

    for (int k = 0; k <= K; ++k){
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * L + 2 * k + n + 3]
                          * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                          * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                          * std::pow(Z + data.L_1, n + 1)
                          * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * k + n + 4)
                             - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[2 * L + 3]
                          * std::pow(2.0, 2 * k + 1)
                          * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

            q_2 += num / den;
        }
    }

    FP_TYPE q_3 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * l + 2 * K + n + 3]
                          * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                          * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                          * std::pow(Z + data.L_1, n + 1)
                          * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * l + 2 * K + n + 4)
                             - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[2 * K + 3]
                          * factorial_array[l] * factorial_array[l + 1]
                          * std::pow(2.0, 2 * l + 1)
                          * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

            q_3 += num / den;
        }
    }

    FP_TYPE q_4 = 0.0;


    for (int k = 0; k <= K; ++k) {
        FP_TYPE num = std::exp(std::lgamma((FP_TYPE) (2 * L + 2 * k + N) + (FP_TYPE) 4.5))
                      * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                      * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                      * std::pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                      * ((FP_TYPE) (1.0) / std::pow(d, (FP_TYPE) (2 * L + 2 * k + N) + 4.5)
                         - (FP_TYPE) (1.0) / std::pow(data.L_2 + d, (FP_TYPE) (2 * L + 2 * k + N) + 4.5));

        FP_TYPE den = factorial_array[N + 1]
                      * factorial_array[k] * factorial_array[k + 1]
                      * factorial_array[2 * L + 3]
                      * std::pow(2.0, 2 * k + 1)
                      * (FP_TYPE)((2 * k + 3) * (2 * L + 5))
                      * std::sqrt((FP_TYPE) (2 * (2 * N + 3)));

        q_4 += num / den;
    }

    FP_TYPE q_5 = 0.0;

    for (int l = 0; l <= L; ++l) {
        FP_TYPE num = std::exp(std::lgamma((FP_TYPE) (2 * l + 2 * K + N) + (FP_TYPE) 4.5))
                      * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                      * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                      * std::pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                      * ((FP_TYPE) (1.0) / std::pow(d, (FP_TYPE) (2 * l + 2 * K + N) + 4.5)
                         - (FP_TYPE) (1.0) / std::pow(data.L_2 + d, (FP_TYPE) (2 * l + 2 * K + N) + 4.5));

        FP_TYPE den = factorial_array[N + 1]
                      * factorial_array[l] * factorial_array[l + 1]
                      * factorial_array[2 * K + 3]
                      * std::pow(2.0, 2 * l + 1)
                      * (FP_TYPE)((2 * K + 5) * (2 * l + 3))
                      * std::sqrt((FP_TYPE) (2 * (2 * N + 3)));

        q_5 += num / den;
    }

    FP_TYPE q_6 = 0.0;

    for (int n = 0; n <= N; ++n) {
        FP_TYPE num = factorial_array[2 * L + 2 * K + n + 5]
                      * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                      * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                      * std::pow(Z + data.L_1, n + 1)
                      * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * K + n + 6)
                         - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + n + 6));

        FP_TYPE den = factorial_array[n + 1]
                      * factorial_array[2 * K + 3]
                      * factorial_array[2 * L + 3]
                      * (FP_TYPE)((2 * K + 5) * (2 * L + 5));

        q_6 += num / den;
    }

    FP_TYPE q_7 = (std::exp(std::lgamma((FP_TYPE) (2 * L + 2 * K + N) + (FP_TYPE) 6.5))
                   * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                   * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                   * std::pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                   * ((FP_TYPE) (1.0) / std::pow(d, (FP_TYPE) (2 * L + 2 * K + N) + 6.5)
                      - (FP_TYPE) (1.0) / std::pow(data.L_2 + d, (FP_TYPE) (2 * L + 2 * K + N) + 6.5))
                  ) / (factorial_array[N + 1]
                       * factorial_array[2 * K + 3]
                       * factorial_array[2 * L + 3]
                       * (FP_TYPE)((2 * K + 5) * (2 *L + 5))
                       * std::sqrt((FP_TYPE) (2 * (2 * N + 3))));

    FP_TYPE q = q_1 + q_2 + q_3 + q_4 + q_5 + q_6 + q_7;

    FP_TYPE r = q + u;

    r *= (FP_TYPE) (4.0e-7) * local_pi * local_pi * data.N_1 * data.N_2
         / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    if (timing) {
        double interval = std::chrono::duration_cast<std::chrono::duration<double>>(
                std::chrono::high_resolution_clock::now() - begin_time).count();
        std::cout << "Time = " << interval << " s" << std::endl;
    }

    return r;
}

FP_TYPE calculate_inductance_remainder(
        const CoilCalculationData& data,
        const SumPrecisionData &precision,
        const FP_TYPE d,
        const FP_TYPE Z,
        bool timing = false
) {
    std::chrono::high_resolution_clock::time_point begin_time;

    if (timing) {
        begin_time = std::chrono::high_resolution_clock::now();
    }

    auto L = (int) precision.l_terms - 1;
    auto K = (int) precision.k_terms - 1;
    auto N = (int) precision.n_terms - 1;

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

    FP_TYPE abs_Z = std::abs(Z);
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
                          * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                          * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                          * std::pow(abs_Z, N + 2)
                          * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * l + 2 * k + N + 3)
                             - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * k + N + 3));

            FP_TYPE den = factorial_array[N + 2]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[l] * factorial_array[l + 1]
                          * std::pow(2.0, 2 * l + 2 * k + 2)
                          * (FP_TYPE)((2 * k + 3) * (2 * l + 3));

            u_1 += num / den;
        }
    }

    FP_TYPE u_2 = 0.0;

    for (int k = 0; k <= K; ++k) {
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * L + 2 * k + n + 3]
                          * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                          * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                          * std::pow(abs_Z, n + 1)
                          * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * k + n + 4)
                             - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[2 * L + 3]
                          * std::pow(2.0, 2 * k + 1)
                          * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

            u_2 += num / den;
        }
    }

    FP_TYPE u_3 = 0.0;

    for (int l = 0; l <= L; ++l) {
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * l + 2 * K + n + 3]
                          * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                          * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                          * std::pow(abs_Z, n + 1)
                          * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * l + 2 * K + n + 4)
                             - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[2 * K + 3]
                          * factorial_array[l] * factorial_array[l + 1]
                          * std::pow(2.0, 2 * l + 1)
                          * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

            u_3 += num / den;
        }
    }

    FP_TYPE u_4 = 0.0;

    for (int l = 0; l <= L; ++l) {
        FP_TYPE num = factorial_array[2 * l + 2 * K + N + 4]
                      * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                      * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                      * std::pow(abs_Z, N + 2)
                      * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * l + 2 * K + N + 5)
                         - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + N + 5));

        FP_TYPE den = factorial_array[N + 2]
                      * factorial_array[2 * K + 3]
                      * factorial_array[l] * factorial_array[l + 1]
                      * std::pow(2.0, 2 * l + 1)
                      * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

        u_4 += num / den;
    }

    FP_TYPE u_5 = 0.0;

    for (int k = 0; k <= K; ++k) {
        FP_TYPE num = factorial_array[2 * L + 2 * k + N + 4]
                      * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                      * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                      * std::pow(abs_Z, N + 2)
                      * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * k + N + 5)
                         - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + N + 5));

        FP_TYPE den = factorial_array[N + 2]
                      * factorial_array[k] * factorial_array[k + 1]
                      * factorial_array[2 * L + 3]
                      * std::pow(2.0, 2 * k + 1)
                      * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

        u_5 += num / den;
    }

    FP_TYPE u_6 = 0.0;

    for (int n = 0; n <= N; ++n) {
        FP_TYPE num = factorial_array[2 * L + 2 * K + n + 5]
                      * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                      * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                      * std::pow(abs_Z, n + 1)
                      * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * K + n + 6)
                         - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + n + 6));

        FP_TYPE den = factorial_array[n + 1]
                      * factorial_array[2 * K + 3]
                      * factorial_array[2 * L + 3]
                      * (FP_TYPE)((2 * K + 5) * (2 * L + 5));

        u_6 += num / den;
    }

    FP_TYPE u_7 = (factorial_array[2 * L + 2 * K + N + 6]
                   * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                   * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                   * std::pow(abs_Z, N + 2)
                   * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * K + N + 7)
                      - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + N + 7))
                  ) / (factorial_array[N + 2]
                       * factorial_array[2 * K + 3]
                       * factorial_array[2 * L + 3]
                       * (FP_TYPE)((2 * K + 5) * (2 * L + 5)));

    FP_TYPE u = u_1 + u_2 + u_3 + u_4 + u_5 + u_6 + u_7;

    FP_TYPE q_1 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int k = 0; k <= K; ++k) {
            FP_TYPE num = std::exp(std::lgamma((FP_TYPE) (2 * l + 2 * k + N) + (FP_TYPE) 2.5))
                          * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                          * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                          * std::pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                          * ((FP_TYPE) (1.0) / std::pow(d, (FP_TYPE) (2 * l + 2 * k + N) + 2.5)
                             - (FP_TYPE) (1.0) / std::pow(data.L_2 + d, (FP_TYPE) (2 * l + 2 * k + N) + 2.5));

            FP_TYPE den = factorial_array[N + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[l] * factorial_array[l + 1]
                          * std::pow(2.0, 2 * l + 2 * k + 2)
                          * (FP_TYPE)((2 * k + 3) * (2 * l + 3))
                          * std::sqrt((FP_TYPE) (2 * (2 * N + 3)));

            q_1 += num / den;
        }
    }

    FP_TYPE q_2 = 0.0;

    for (int k = 0; k <= K; ++k){
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * L + 2 * k + n + 3]
                          * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                          * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                          * std::pow(Z + data.L_1, n + 1)
                          * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * k + n + 4)
                             - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * k + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[k] * factorial_array[k + 1]
                          * factorial_array[2 * L + 3]
                          * std::pow(2.0, 2 * k + 1)
                          * (FP_TYPE)((2 * k + 3) * (2 * L + 5));

            q_2 += num / den;
        }
    }

    FP_TYPE q_3 = 0.0;

    for (int l = 0; l <= L; ++l){
        for (int n = 0; n <= N; ++n) {
            FP_TYPE num = factorial_array[2 * l + 2 * K + n + 3]
                          * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                          * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                          * std::pow(Z + data.L_1, n + 1)
                          * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * l + 2 * K + n + 4)
                             - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * l + 2 * K + n + 4));

            FP_TYPE den = factorial_array[n + 1]
                          * factorial_array[2 * K + 3]
                          * factorial_array[l] * factorial_array[l + 1]
                          * std::pow(2.0, 2 * l + 1)
                          * (FP_TYPE)((2 * K + 5) * (2 * l + 3));

            q_3 += num / den;
        }
    }

    FP_TYPE q_4 = 0.0;


    for (int k = 0; k <= K; ++k) {
        FP_TYPE num = std::exp(std::lgamma((FP_TYPE) (2 * L + 2 * k + N) + (FP_TYPE) 4.5))
                      * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                      * (std::pow(data.R_2, 2 * k + 3) - std::pow(data.r_2, 2 * k + 3))
                      * std::pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                      * ((FP_TYPE) (1.0) / std::pow(d, (FP_TYPE) (2 * L + 2 * k + N) + 4.5)
                         - (FP_TYPE) (1.0) / std::pow(data.L_2 + d, (FP_TYPE) (2 * L + 2 * k + N) + 4.5));

        FP_TYPE den = factorial_array[N + 1]
                      * factorial_array[k] * factorial_array[k + 1]
                      * factorial_array[2 * L + 3]
                      * std::pow(2.0, 2 * k + 1)
                      * (FP_TYPE)((2 * k + 3) * (2 * L + 5))
                      * std::sqrt((FP_TYPE) (2 * (2 * N + 3)));

        q_4 += num / den;
    }

    FP_TYPE q_5 = 0.0;

    for (int l = 0; l <= L; ++l) {
        FP_TYPE num = std::exp(std::lgamma((FP_TYPE) (2 * l + 2 * K + N) + (FP_TYPE) 4.5))
                      * (std::pow(data.R_1, 2 * l + 3) - std::pow(data.r_1, 2 * l + 3))
                      * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                      * std::pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                      * ((FP_TYPE) (1.0) / std::pow(d, (FP_TYPE) (2 * l + 2 * K + N) + 4.5)
                         - (FP_TYPE) (1.0) / std::pow(data.L_2 + d, (FP_TYPE) (2 * l + 2 * K + N) + 4.5));

        FP_TYPE den = factorial_array[N + 1]
                      * factorial_array[l] * factorial_array[l + 1]
                      * factorial_array[2 * K + 3]
                      * std::pow(2.0, 2 * l + 1)
                      * (FP_TYPE)((2 * K + 5) * (2 * l + 3))
                      * std::sqrt((FP_TYPE) (2 * (2 * N + 3)));

        q_5 += num / den;
    }

    FP_TYPE q_6 = 0.0;

    for (int n = 0; n <= N; ++n) {
        FP_TYPE num = factorial_array[2 * L + 2 * K + n + 5]
                      * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                      * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                      * std::pow(Z + data.L_1, n + 1)
                      * ((FP_TYPE) (1.0) / std::pow(Z + data.L_1 + d, 2 * L + 2 * K + n + 6)
                         - (FP_TYPE) (1.0) / std::pow(Z + data.L_1 + data.L_2 + d, 2 * L + 2 * K + n + 6));

        FP_TYPE den = factorial_array[n + 1]
                      * factorial_array[2 * K + 3]
                      * factorial_array[2 * L + 3]
                      * (FP_TYPE)((2 * K + 5) * (2 * L + 5));

        q_6 += num / den;
    }

    FP_TYPE q_7 = (std::exp(std::lgamma((FP_TYPE) (2 * L + 2 * K + N) + (FP_TYPE) 6.5))
                   * (std::pow(data.R_1, 2 * L + 5) - std::pow(data.r_1, 2 * L + 5))
                   * (std::pow(data.R_2, 2 * K + 5) - std::pow(data.r_2, 2 * K + 5))
                   * std::pow(Z + data.L_1, (FP_TYPE) N + 1.5)
                   * ((FP_TYPE) (1.0) / std::pow(d, (FP_TYPE) (2 * L + 2 * K + N) + 6.5)
                      - (FP_TYPE) (1.0) / std::pow(data.L_2 + d, (FP_TYPE) (2 * L + 2 * K + N) + 6.5))
                  ) / (factorial_array[N + 1]
                       * factorial_array[2 * K + 3]
                       * factorial_array[2 * L + 3]
                       * (FP_TYPE)((2 * K + 5) * (2 *L + 5))
                       * std::sqrt((FP_TYPE) (2 * (2 * N + 3))));

    FP_TYPE q = q_1 + q_2 + q_3 + q_4 + q_5 + q_6 + q_7;

    FP_TYPE r = q + u;

    r *= (FP_TYPE) (4.0e-7) * local_pi * local_pi * data.N_1 * data.N_2
         / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    return r;
}

#endif //VECTOR_CASE_REMAINDER_ESTIMATE_HPP
