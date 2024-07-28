#pragma clang diagnostic push
#pragma ide diagnostic ignored "portability-simd-intrinsics"
#ifndef VECTOR_CASE_FAR_INDUCTANCE_HPP
#define VECTOR_CASE_FAR_INDUCTANCE_HPP

#include <iostream>
#include <chrono>
#include <cmath>

#include "structs.h"
#include "sum_lookup_table_far.h"

#if defined(USE_AVX) || defined(USE_AVX512)
#include "immintrin.h"
#endif


double calculate_mutual_inductance_far(
        const CoilCalculationData& data,
        const SumPrecisionData &precision,
        const FP_TYPE d
) {
    if (d + 0.5 * data.L_1 <= std::max(data.R_1, data.R_2))
        throw std::invalid_argument(
            "The distance d is too small, the method will not converge"
        );

    // Useful calculations that can be performed at compile time, before the main loop
    const FP_TYPE denom_1 = (FP_TYPE) (1.0) / ((FP_TYPE) 0.5 * data.L_1 + d);
    const FP_TYPE denom_2 = (FP_TYPE) (1.0) / ((FP_TYPE) 0.5 * data.L_1 + d + data.L_2);

    const FP_TYPE L_1_sq = data.L_1 * data.L_1;
    const FP_TYPE R_1_sq = data.R_1 * data.R_1;
    const FP_TYPE r_1_sq = data.r_1 * data.r_1;

    const FP_TYPE R_2_sq = data.R_2 * data.R_2;
    const FP_TYPE r_2_sq = data.r_2 * data.r_2;

    const FP_TYPE denom_1_sq = denom_1 * denom_1;
    const FP_TYPE denom_2_sq = denom_2 * denom_2;

    // Variable which will store the result of the main loop
    FP_TYPE M_12 = 0.0;

#ifdef USE_AVX
    __m256d M_12_vec = _mm256_set1_pd(0.0);
#elif defined(USE_AVX512)
    __m512d M_12_vec = _mm512_set1_pd(0.0);
#endif

    // Inner loop lookup table for efficient vectorization
    FP_TYPE inner_denom_1_arr[MAX_TERMS_FAR]{};
    FP_TYPE inner_denom_2_arr[MAX_TERMS_FAR]{};
    FP_TYPE inner_L_1_arr[MAX_TERMS_FAR]{};

    auto temp_denom_1 = (FP_TYPE) 1.0;
    auto temp_denom_2 = (FP_TYPE) 1.0;
    auto temp_L_1_plus_Z = (FP_TYPE) 1.0;

    for (uint32_t n = 0; n < precision.n_terms; ++n) {
        if (n > 0) {
            temp_denom_1 *= denom_1_sq;
            temp_denom_2 *= denom_2_sq;
            temp_L_1_plus_Z *= L_1_sq;
        }

        inner_denom_1_arr[n] = temp_denom_1;
        inner_denom_2_arr[n] = temp_denom_2;
        inner_L_1_arr[n] = temp_L_1_plus_Z;
    }

    // These three variables need to be remembered for each iteration
    // of the inner loops to be able to restore index values
    FP_TYPE loop_denom_1 = denom_1_sq;
    FP_TYPE loop_denom_2 = denom_2_sq;

    // Value that has to be set for the first loop
    FP_TYPE loop_R_1 = R_1_sq * data.R_1;
    FP_TYPE loop_r_1 = r_1_sq * data.r_1;

    // Main loop
    for (uint32_t l = 0; l < precision.l_terms; ++l) {

        if (l > 0) {
            loop_R_1 *= R_1_sq;
            loop_r_1 *= r_1_sq;

            loop_denom_1 *= denom_1_sq;
            loop_denom_2 *= denom_2_sq;
        }

        FP_TYPE loop_R_1_sub_r_1 = loop_R_1 - loop_r_1;

        // Restore point for first loop
        FP_TYPE save_first_loop_denom_1 = loop_denom_1;
        FP_TYPE save_first_loop_denom_2 = loop_denom_2;

        // Values that have to be set for the second loop
        FP_TYPE loop_R_2 = R_2_sq * data.R_2;
        FP_TYPE loop_r_2 = r_2_sq * data.r_2;

        for (uint32_t k = 0; k < precision.k_terms; ++k) {

            if (k > 0) {
                loop_R_2 *= R_2_sq;
                loop_r_2 *= r_2_sq;

                loop_denom_1 *= denom_1_sq;
                loop_denom_2 *= denom_2_sq;
            }

            FP_TYPE loop_R_2_sub_r_2 = loop_R_2 - loop_r_2;

            // Restore point for second loop
            FP_TYPE save_second_loop_denom_1 = loop_denom_1;
            FP_TYPE save_second_loop_denom_2 = loop_denom_2;

        #ifdef USE_AVX
            __m256d loop_denom_1_vec = _mm256_set1_pd(loop_denom_1);
            __m256d loop_denom_2_vec = _mm256_set1_pd(loop_denom_2);

            __m256d loop_R_1_sub_r_1_vec = _mm256_set1_pd(loop_R_1_sub_r_1);
            __m256d loop_R_2_sub_r_2_vec = _mm256_set1_pd(loop_R_2_sub_r_2);

            for (uint32_t n = 0; n < precision.n_terms; n += 4) {
                __m256d inner_L_1_vec = _mm256_loadu_pd(&inner_L_1_arr[n]);
                __m256d inner_loop_denom_1_vec = _mm256_loadu_pd(&inner_denom_1_arr[n]);
                __m256d inner_loop_denom_2_vec = _mm256_loadu_pd(&inner_denom_2_arr[n]);

                inner_loop_denom_1_vec = _mm256_mul_pd(inner_loop_denom_1_vec, loop_denom_1_vec);
                inner_loop_denom_2_vec = _mm256_mul_pd(inner_loop_denom_2_vec, loop_denom_2_vec);

                __m256d numerator = _mm256_mul_pd(
                    _mm256_mul_pd(inner_L_1_vec, loop_R_1_sub_r_1_vec),
                    _mm256_mul_pd(
                        loop_R_2_sub_r_2_vec,
                        _mm256_sub_pd(inner_loop_denom_1_vec, inner_loop_denom_2_vec)
                    )
                );

                __m256d lookup_table_vec = _mm256_loadu_pd(&lookup_table_far[l][k][n]);

                M_12_vec = _mm256_add_pd(M_12_vec, _mm256_mul_pd(lookup_table_vec, numerator));
            }
        #elif defined(USE_AVX512)
            __m512d loop_denom_1_vec = _mm512_set1_pd(loop_denom_1);
            __m512d loop_denom_2_vec = _mm512_set1_pd(loop_denom_2);

            __m512d loop_R_1_sub_r_1_vec = _mm512_set1_pd(loop_R_1_sub_r_1);
            __m512d loop_R_2_sub_r_2_vec = _mm512_set1_pd(loop_R_2_sub_r_2);

            for (uint32_t n = 0; n < precision.n_terms; n += 8) {
                __m512d inner_loop_denom_1_vec = _mm512_loadu_pd(&inner_denom_1_arr[n]);
                __m512d inner_loop_denom_2_vec = _mm512_loadu_pd(&inner_denom_2_arr[n]);
                __m512d inner_L_1_vec = _mm512_loadu_pd(&inner_L_1_arr[n]);

                inner_loop_denom_1_vec = _mm512_mul_pd(inner_loop_denom_1_vec, loop_denom_1_vec);
                inner_loop_denom_2_vec = _mm512_mul_pd(inner_loop_denom_2_vec, loop_denom_2_vec);

                __m512d numerator = _mm512_mul_pd(
                    _mm512_mul_pd(loop_R_2_sub_r_2_vec, loop_R_1_sub_r_1_vec),
                    _mm512_mul_pd(
                        inner_L_1_vec,
                        _mm512_sub_pd(inner_loop_denom_1_vec, inner_loop_denom_2_vec)
                    )
                );

                __m512d lookup_table_vec = _mm512_loadu_pd(&lookup_table_far[l][k][n]);

                M_12_vec = _mm512_add_pd(M_12_vec, _mm512_mul_pd(lookup_table_vec, numerator));
            }
        #else
            for (uint32_t n = 0; n < precision.n_terms; ++n) {

                FP_TYPE inner_loop_denom_1 = loop_denom_1 * inner_denom_1_arr[n];
                FP_TYPE inner_loop_denom_2 = loop_denom_2 * inner_denom_2_arr[n];

                FP_TYPE numerator = inner_L_1_arr[n]
                                  * loop_R_1_sub_r_1
                                  * loop_R_2_sub_r_2
                                  * (inner_loop_denom_1 - inner_loop_denom_2);

                M_12 += lookup_table_far[l][k][n] * numerator;
            }
        #endif // USE_AVX

            loop_denom_1 = save_second_loop_denom_1;
            loop_denom_2 = save_second_loop_denom_2;
        }

        loop_denom_1 = save_first_loop_denom_1;
        loop_denom_2 = save_first_loop_denom_2;
    }

#ifdef USE_AVX
    double temp[4];
    _mm256_storeu_pd(temp, M_12_vec);
    M_12 += temp[0] + temp[1] + temp[2] + temp[3];
#elif defined(USE_AVX512)
    M_12 += _mm512_reduce_add_pd(M_12_vec);
#endif

    M_12 *= (FP_TYPE) (4.0e-7) * local_pi * local_pi * data.N_1 * data.N_2
            / (data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    return M_12;
}

#endif //VECTOR_CASE_FAR_INDUCTANCE_HPP
#pragma clang diagnostic pop