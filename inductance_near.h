#ifndef VECTOR_CASE_INDUCTANCE_NEAR_H
#define VECTOR_CASE_INDUCTANCE_NEAR_H

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#include "structs.h"
#include "sum_lookup_table_near.h"

#if defined(USE_SSE) || defined(USE_AVX) || defined(USE_AVX512)
#include "immintrin.h"
#endif


FP_TYPE calculate_mutual_inductance_near(
        const CoilCalculationData data,
        const SumPrecisionData precision,
        const FP_TYPE d,
        const FP_TYPE Z
) {
    // Useful calculations that can be performed at compile time, before the main loop
    const FP_TYPE denom_1 = (FP_TYPE) 1.0 / (Z + data.L_1 + d);
    const FP_TYPE denom_2 = (FP_TYPE) 1.0 / (Z + data.L_1 + data.L_2 + d);

    const FP_TYPE R_1_sq = data.R_1 * data.R_1;
    const FP_TYPE r_1_sq = data.r_1 * data.r_1;

    const FP_TYPE R_2_sq = data.R_2 * data.R_2;
    const FP_TYPE r_2_sq = data.r_2 * data.r_2;

    const FP_TYPE denom_1_sq = denom_1 * denom_1;
    const FP_TYPE denom_2_sq = denom_2 * denom_2;

    const FP_TYPE L_1_plus_Z = data.L_1 + Z;

    // Variable which will store the result of the main loop
    FP_TYPE M_12 = 0.0;

#if defined(USE_SSE)
    __m128d M_12_vec = _mm_set1_pd(0.0);
#elif defined(USE_AVX)
    __m256d M_12_vec = _mm256_set1_pd(0.0);
#elif defined(USE_AVX512)
    __m512d M_12_vec = _mm512_set1_pd(0.0);
#endif

    // Inner loop lookup table for efficient vectorization
    FP_TYPE inner_denom_1_arr[MAX_TERMS_NEAR] = {0.0};
    FP_TYPE inner_denom_2_arr[MAX_TERMS_NEAR] = {0.0};
    FP_TYPE inner_L_1_plus_Z_sub_Z_arr[MAX_TERMS_NEAR] = {0.0};

    FP_TYPE temp_denom_1 = 1.0;
    FP_TYPE temp_denom_2 = 1.0;
    FP_TYPE temp_L_1_plus_Z = L_1_plus_Z;
    FP_TYPE temp_Z = Z;

    for (uint32_t n = 0; n < precision.n_terms; ++n) {
        if (n > 0) {
            temp_denom_1 *= denom_1;
            temp_denom_2 *= denom_2;

            temp_L_1_plus_Z *= L_1_plus_Z;
            temp_Z *= Z;
        }

        inner_denom_1_arr[n] = temp_denom_1;
        inner_denom_2_arr[n] = temp_denom_2;
        inner_L_1_plus_Z_sub_Z_arr[n] = temp_L_1_plus_Z - temp_Z;
    }

    // These three variables need to be remembered for each iteration
    // of the inner loops to be able to restore index values
    FP_TYPE loop_denom_1 = denom_1_sq;
    FP_TYPE loop_denom_2 = denom_2_sq;

    // Values that have to be set for the first loop
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

        #if defined(USE_SSE)
            __m128d loop_denom_1_vec = _mm_set1_pd(loop_denom_1);
            __m128d loop_denom_2_vec = _mm_set1_pd(loop_denom_2);

            __m128d loop_R_1_sub_r_1_vec = _mm_set1_pd(loop_R_1_sub_r_1);
            __m128d loop_R_2_sub_r_2_vec = _mm_set1_pd(loop_R_2_sub_r_2);

            for (uint32_t n = 0; n < precision.n_terms; n += 2) {

                __m128d loop_L_1_plus_Z_sub_Z_vec = _mm_loadu_pd(&inner_L_1_plus_Z_sub_Z_arr[n]);
                __m128d inner_loop_denom_1_vec = _mm_loadu_pd(&inner_denom_1_arr[n]);
                __m128d inner_loop_denom_2_vec = _mm_loadu_pd(&inner_denom_2_arr[n]);

                inner_loop_denom_1_vec = _mm_mul_pd(loop_denom_1_vec, inner_loop_denom_1_vec);
                inner_loop_denom_2_vec = _mm_mul_pd(loop_denom_2_vec, inner_loop_denom_2_vec);

                __m128d numerator = _mm_mul_pd(
                    _mm_mul_pd(loop_R_2_sub_r_2_vec, loop_R_1_sub_r_1_vec),
                    _mm_mul_pd(
                        loop_L_1_plus_Z_sub_Z_vec,
                        _mm_sub_pd(inner_loop_denom_1_vec, inner_loop_denom_2_vec)
                    )
                );

                __m128d lookup_table_vec = _mm_loadu_pd(&lookup_table_near[l][k][n]);

                M_12_vec = _mm_add_pd(M_12_vec, _mm_mul_pd(numerator, lookup_table_vec));
            }
        #elif defined(USE_AVX)
            __m256d loop_denom_1_vec = _mm256_set1_pd(loop_denom_1);
            __m256d loop_denom_2_vec = _mm256_set1_pd(loop_denom_2);

            __m256d loop_R_1_sub_r_1_vec = _mm256_set1_pd(loop_R_1_sub_r_1);
            __m256d loop_R_2_sub_r_2_vec = _mm256_set1_pd(loop_R_2_sub_r_2);

            for (uint32_t n = 0; n < precision.n_terms; n += 4) {

                __m256d loop_L_1_plus_Z_sub_Z_vec = _mm256_loadu_pd(&inner_L_1_plus_Z_sub_Z_arr[n]);
                __m256d inner_loop_denom_1_vec = _mm256_loadu_pd(&inner_denom_1_arr[n]);
                __m256d inner_loop_denom_2_vec = _mm256_loadu_pd(&inner_denom_2_arr[n]);

                inner_loop_denom_1_vec = _mm256_mul_pd(loop_denom_1_vec, inner_loop_denom_1_vec);
                inner_loop_denom_2_vec = _mm256_mul_pd(loop_denom_2_vec, inner_loop_denom_2_vec);

                __m256d numerator = _mm256_mul_pd(
                    _mm256_mul_pd(loop_R_2_sub_r_2_vec, loop_R_1_sub_r_1_vec),
                    _mm256_mul_pd(
                        loop_L_1_plus_Z_sub_Z_vec,
                        _mm256_sub_pd(inner_loop_denom_1_vec, inner_loop_denom_2_vec)
                    )
                );

                __m256d lookup_table_vec = _mm256_loadu_pd(&lookup_table_near[l][k][n]);

                M_12_vec = _mm256_add_pd(M_12_vec, _mm256_mul_pd(numerator, lookup_table_vec));
            }
        #elif defined(USE_AVX512)
            __m512d loop_denom_1_vec = _mm512_set1_pd(loop_denom_1);
            __m512d loop_denom_2_vec = _mm512_set1_pd(loop_denom_2);

            __m512d loop_R_1_sub_r_1_vec = _mm512_set1_pd(loop_R_1_sub_r_1);
            __m512d loop_R_2_sub_r_2_vec = _mm512_set1_pd(loop_R_2_sub_r_2);

            for (uint32_t n = 0; n < precision.n_terms; n += 8) {

                __m512d loop_L_1_plus_Z_sub_Z_vec = _mm512_loadu_pd(&inner_L_1_plus_Z_sub_Z_arr[n]);
                __m512d inner_loop_denom_1_vec = _mm512_loadu_pd(&inner_denom_1_arr[n]);
                __m512d inner_loop_denom_2_vec = _mm512_loadu_pd(&inner_denom_2_arr[n]);

                inner_loop_denom_1_vec = _mm512_mul_pd(loop_denom_1_vec, inner_loop_denom_1_vec);
                inner_loop_denom_2_vec = _mm512_mul_pd(loop_denom_2_vec, inner_loop_denom_2_vec);

                __m512d numerator = _mm512_mul_pd(
                    _mm512_mul_pd(loop_R_2_sub_r_2_vec, loop_R_1_sub_r_1_vec),
                    _mm512_mul_pd(
                        loop_L_1_plus_Z_sub_Z_vec,
                        _mm512_sub_pd(inner_loop_denom_1_vec,inner_loop_denom_2_vec)
                    )
                );

                __m512d lookup_table_vec = _mm512_loadu_pd(&lookup_table_near[l][k][n]);

                M_12_vec = _mm512_add_pd(M_12_vec, _mm512_mul_pd(numerator, lookup_table_vec));
            }
        #else

            for (uint32_t n = 0; n < precision.n_terms; ++n) {

                FP_TYPE inner_loop_denom_1 = loop_denom_1 * inner_denom_1_arr[n];
                FP_TYPE inner_loop_denom_2 = loop_denom_2 * inner_denom_2_arr[n];

                FP_TYPE numerator = inner_L_1_plus_Z_sub_Z_arr[n]
                                    * loop_R_1_sub_r_1
                                    * loop_R_2_sub_r_2
                                    * (inner_loop_denom_1 - inner_loop_denom_2);

                M_12 += lookup_table_near[l][k][n] * numerator;
            }
        #endif // USE_AVX

            loop_denom_1 = save_second_loop_denom_1;
            loop_denom_2 = save_second_loop_denom_2;
        }

        loop_denom_1 = save_first_loop_denom_1;
        loop_denom_2 = save_first_loop_denom_2;
    }

#if defined(USE_SSE)
    double temp[2];
    _mm_storeu_pd(temp, M_12_vec);
    M_12 += temp[0] + temp[1];
#elif defined(USE_AVX)
    double temp[4];
    _mm256_storeu_pd(temp, M_12_vec);
    M_12 += temp[0] + temp[1] + temp[2] + temp[3];
#elif defined(USE_AVX512)
    M_12 += _mm512_reduce_add_pd(M_12_vec);
#endif

    M_12 *= 4.0 * local_pi * local_pi * 1e-7 * data.N_1 * data.N_2
          / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    return M_12;
}

FP_TYPE calculate_mutual_inductance_near_dz(
        const CoilCalculationData data,
        const SumPrecisionData precision,
        const FP_TYPE d,
        const FP_TYPE Z
) {
    // Useful calculations that can be performed at compile time, before the main loop
    const FP_TYPE denom_1 = 1.0 / (Z + data.L_1 + d);
    const FP_TYPE denom_2 = 1.0 / (Z + data.L_1 + data.L_2 + d);

    const FP_TYPE R_1_sq = data.R_1 * data.R_1;
    const FP_TYPE r_1_sq = data.r_1 * data.r_1;

    const FP_TYPE R_2_sq = data.R_2 * data.R_2;
    const FP_TYPE r_2_sq = data.r_2 * data.r_2;

    const FP_TYPE denom_1_sq = denom_1 * denom_1;
    const FP_TYPE denom_2_sq = denom_2 * denom_2;

    const FP_TYPE L_1_plus_Z = data.L_1 + Z;

    // Variable which will store the result of the main loop
    FP_TYPE M_12 = 0.0;

#if defined(USE_SSE)
    __m128d M_12_vec = _mm_set1_pd(0.0);
    __m128d denom_1_vec = _mm_set1_pd(denom_1);
    __m128d denom_2_vec = _mm_set1_pd(denom_2);
#elif defined(USE_AVX)
    __m256d M_12_vec = _mm256_set1_pd(0.0);
    __m256d denom_1_vec = _mm256_set1_pd(denom_1);
    __m256d denom_2_vec = _mm256_set1_pd(denom_2);
#elif defined(USE_AVX512)
    __m512d M_12_vec = _mm512_set1_pd(0.0);
    __m512d denom_1_vec = _mm512_set1_pd(denom_1);
    __m512d denom_2_vec = _mm512_set1_pd(denom_2);
#endif

    // Inner loop lookup table for efficient vectorization
    FP_TYPE inner_denom_1_arr[MAX_TERMS_NEAR] = {0.0};
    FP_TYPE inner_denom_2_arr[MAX_TERMS_NEAR] = {0.0};
    FP_TYPE inner_L_1_plus_Z_sub_Z_arr[MAX_TERMS_NEAR] = {0.0};
    FP_TYPE inner_L_1_plus_Z_sub_Z_arr_2[MAX_TERMS_NEAR] = {0.0};

    FP_TYPE temp_denom_1 = 1.0;
    FP_TYPE temp_denom_2 = 1.0;
    FP_TYPE temp_L_1_plus_Z = 1.0;
    FP_TYPE temp_Z = 1.0;

    for (uint32_t n = 0; n < precision.n_terms; ++n) {
        if (n > 0) {
            temp_denom_1 *= denom_1;
            temp_denom_2 *= denom_2;

            temp_L_1_plus_Z *= L_1_plus_Z;
            temp_Z *= Z;
        }

        inner_denom_1_arr[n] = temp_denom_1;
        inner_denom_2_arr[n] = temp_denom_2;
        inner_L_1_plus_Z_sub_Z_arr[n] = temp_L_1_plus_Z * L_1_plus_Z - temp_Z * Z;
        inner_L_1_plus_Z_sub_Z_arr_2[n] = (FP_TYPE) (n + 1) * (temp_L_1_plus_Z - temp_Z);
    }

    // These three variables need to be remembered for each iteration
    // of the inner loops to be able to restore index values
    FP_TYPE loop_denom_1 = denom_1_sq;
    FP_TYPE loop_denom_2 = denom_2_sq;

    // Values that have to be set for the first loop
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

        #if defined(USE_SSE)
            __m128d loop_denom_1_vec = _mm_set1_pd(loop_denom_1);
            __m128d loop_denom_2_vec = _mm_set1_pd(loop_denom_2);

            __m128d loop_R_1_sub_r_1_vec = _mm_set1_pd(loop_R_1_sub_r_1);
            __m128d loop_R_2_sub_r_2_vec = _mm_set1_pd(loop_R_2_sub_r_2);

            for (uint32_t n = 0; n < precision.n_terms; n += 2) {

                __m128d loop_L_1_plus_Z_sub_Z_vec = _mm_loadu_pd(&inner_L_1_plus_Z_sub_Z_arr[n]);
                __m128d loop_L_1_plus_Z_sub_Z_vec_2 = _mm_loadu_pd(&inner_L_1_plus_Z_sub_Z_arr_2[n]);

                __m128d inner_loop_denom_1_vec = _mm_loadu_pd(&inner_denom_1_arr[n]);
                __m128d inner_loop_denom_2_vec = _mm_loadu_pd(&inner_denom_2_arr[n]);

                inner_loop_denom_1_vec = _mm_mul_pd(loop_denom_1_vec, inner_loop_denom_1_vec);
                inner_loop_denom_2_vec = _mm_mul_pd(loop_denom_2_vec, inner_loop_denom_2_vec);

                // Mind the endian order
                __m128d factor = _mm_set_pd(
                        -(double) (2 * k + 2 * l + n + 2 + 1),
                        -(double) (2 * k + 2 * l + n + 2 + 0)
                );

                __m128d first = _mm_mul_pd(loop_L_1_plus_Z_sub_Z_vec_2,
                                           _mm_sub_pd(inner_loop_denom_1_vec, inner_loop_denom_2_vec)
                );

                inner_loop_denom_1_vec = _mm_mul_pd(inner_loop_denom_1_vec, denom_1_vec);
                inner_loop_denom_2_vec = _mm_mul_pd(inner_loop_denom_2_vec, denom_2_vec);

                __m128d second = _mm_mul_pd(
                        _mm_mul_pd(loop_L_1_plus_Z_sub_Z_vec, factor),
                        _mm_sub_pd(inner_loop_denom_1_vec, inner_loop_denom_2_vec)
                );

                __m128d numerator = _mm_mul_pd(
                        _mm_mul_pd(loop_R_2_sub_r_2_vec, loop_R_1_sub_r_1_vec),
                        _mm_add_pd(first, second)
                );

                __m128d lookup_table_vec = _mm_loadu_pd(&lookup_table_near[l][k][n]);
                M_12_vec = _mm_add_pd(M_12_vec, _mm_mul_pd(numerator, lookup_table_vec));
            }
        #elif defined(USE_AVX)
            __m256d loop_denom_1_vec = _mm256_set1_pd(loop_denom_1);
            __m256d loop_denom_2_vec = _mm256_set1_pd(loop_denom_2);

            __m256d loop_R_1_sub_r_1_vec = _mm256_set1_pd(loop_R_1_sub_r_1);
            __m256d loop_R_2_sub_r_2_vec = _mm256_set1_pd(loop_R_2_sub_r_2);

            for (uint32_t n = 0; n < precision.n_terms; n += 4) {

                __m256d loop_L_1_plus_Z_sub_Z_vec = _mm256_loadu_pd(&inner_L_1_plus_Z_sub_Z_arr[n]);
                __m256d loop_L_1_plus_Z_sub_Z_vec_2 = _mm256_loadu_pd(&inner_L_1_plus_Z_sub_Z_arr_2[n]);

                __m256d inner_loop_denom_1_vec = _mm256_loadu_pd(&inner_denom_1_arr[n]);
                __m256d inner_loop_denom_2_vec = _mm256_loadu_pd(&inner_denom_2_arr[n]);

                inner_loop_denom_1_vec = _mm256_mul_pd(loop_denom_1_vec, inner_loop_denom_1_vec);
                inner_loop_denom_2_vec = _mm256_mul_pd(loop_denom_2_vec, inner_loop_denom_2_vec);

                // Mind the endian order
                __m256d factor = _mm256_set_pd(
                    -(double) (2 * k + 2 * l + n + 2 + 3),
                    -(double) (2 * k + 2 * l + n + 2 + 2),
                    -(double) (2 * k + 2 * l + n + 2 + 1),
                    -(double) (2 * k + 2 * l + n + 2 + 0)
                );

                __m256d first = _mm256_mul_pd(loop_L_1_plus_Z_sub_Z_vec_2,
                    _mm256_sub_pd(inner_loop_denom_1_vec, inner_loop_denom_2_vec)
                );

                inner_loop_denom_1_vec = _mm256_mul_pd(inner_loop_denom_1_vec, denom_1_vec);
                inner_loop_denom_2_vec = _mm256_mul_pd(inner_loop_denom_2_vec, denom_2_vec);

                __m256d second = _mm256_mul_pd(
                    _mm256_mul_pd(loop_L_1_plus_Z_sub_Z_vec, factor),
                    _mm256_sub_pd(inner_loop_denom_1_vec, inner_loop_denom_2_vec)
                );

                __m256d numerator = _mm256_mul_pd(
                    _mm256_mul_pd(loop_R_2_sub_r_2_vec, loop_R_1_sub_r_1_vec),
                    _mm256_add_pd(first, second)
                );

                __m256d lookup_table_vec = _mm256_loadu_pd(&lookup_table_near[l][k][n]);
                M_12_vec = _mm256_add_pd(M_12_vec, _mm256_mul_pd(numerator, lookup_table_vec));
            }
        #elif defined(USE_AVX512)
            __m512d loop_denom_1_vec = _mm512_set1_pd(loop_denom_1);
            __m512d loop_denom_2_vec = _mm512_set1_pd(loop_denom_2);

            __m512d loop_R_1_sub_r_1_vec = _mm512_set1_pd(loop_R_1_sub_r_1);
            __m512d loop_R_2_sub_r_2_vec = _mm512_set1_pd(loop_R_2_sub_r_2);

            for (uint32_t n = 0; n < precision.n_terms; n += 8) {

                __m512d loop_L_1_plus_Z_sub_Z_vec = _mm512_loadu_pd(&inner_L_1_plus_Z_sub_Z_arr[n]);
                __m512d loop_L_1_plus_Z_sub_Z_vec_2 = _mm512_loadu_pd(&inner_L_1_plus_Z_sub_Z_arr_2[n]);

                __m512d inner_loop_denom_1_vec = _mm512_loadu_pd(&inner_denom_1_arr[n]);
                __m512d inner_loop_denom_2_vec = _mm512_loadu_pd(&inner_denom_2_arr[n]);

                inner_loop_denom_1_vec = _mm512_mul_pd(loop_denom_1_vec, inner_loop_denom_1_vec);
                inner_loop_denom_2_vec = _mm512_mul_pd(loop_denom_2_vec, inner_loop_denom_2_vec);

                // Mind the endian order
                __m512d factor = _mm512_set_pd(
                        -(double) (2 * k + 2 * l + n + 2 + 7),
                        -(double) (2 * k + 2 * l + n + 2 + 6),
                        -(double) (2 * k + 2 * l + n + 2 + 5),
                        -(double) (2 * k + 2 * l + n + 2 + 4),
                        -(double) (2 * k + 2 * l + n + 2 + 3),
                        -(double) (2 * k + 2 * l + n + 2 + 2),
                        -(double) (2 * k + 2 * l + n + 2 + 1),
                        -(double) (2 * k + 2 * l + n + 2 + 0)
                );

                __m512d first = _mm512_mul_pd(
                    loop_L_1_plus_Z_sub_Z_vec_2,
                    _mm512_sub_pd(inner_loop_denom_1_vec, inner_loop_denom_2_vec)
                );

                inner_loop_denom_1_vec = _mm512_mul_pd(inner_loop_denom_1_vec, denom_1_vec);
                inner_loop_denom_2_vec = _mm512_mul_pd(inner_loop_denom_2_vec, denom_2_vec);

                __m512d second = _mm512_mul_pd(
                    _mm512_mul_pd(loop_L_1_plus_Z_sub_Z_vec, factor),
                    _mm512_sub_pd(inner_loop_denom_1_vec, inner_loop_denom_2_vec)
                );

                __m512d numerator = _mm512_mul_pd(
                    _mm512_mul_pd(loop_R_2_sub_r_2_vec, loop_R_1_sub_r_1_vec),
                    _mm512_add_pd(first, second)
                );

                __m512d lookup_table_vec = _mm512_loadu_pd(&lookup_table_near[l][k][n]);
                M_12_vec = _mm512_add_pd(M_12_vec, _mm512_mul_pd(numerator, lookup_table_vec));
            }
        #else

            for (uint32_t n = 0; n < precision.n_terms; ++n) {

                FP_TYPE inner_loop_denom_1_1 = loop_denom_1 * inner_denom_1_arr[n];
                FP_TYPE inner_loop_denom_2_1 = loop_denom_2 * inner_denom_2_arr[n];

                FP_TYPE inner_loop_denom_1_2 = inner_loop_denom_1_1 * denom_1;
                FP_TYPE inner_loop_denom_2_2 = inner_loop_denom_2_1 * denom_2;
                auto factor = -(FP_TYPE) (2 * k + 2 * l + n + 2);

                FP_TYPE numerator = loop_R_1_sub_r_1
                                    * loop_R_2_sub_r_2
                                    * (inner_L_1_plus_Z_sub_Z_arr_2[n]
                                       * (inner_loop_denom_1_1 - inner_loop_denom_2_1)
                                     + inner_L_1_plus_Z_sub_Z_arr[n] * factor
                                       * (inner_loop_denom_1_2 - inner_loop_denom_2_2)
                                      );

                M_12 += numerator * lookup_table_near[l][k][n];
            }
        #endif

            loop_denom_1 = save_second_loop_denom_1;
            loop_denom_2 = save_second_loop_denom_2;
        }

        loop_denom_1 = save_first_loop_denom_1;
        loop_denom_2 = save_first_loop_denom_2;
    }
#if defined(USE_SSE)
    double temp[2];
    _mm_storeu_pd(temp, M_12_vec);
    M_12 += temp[0] + temp[1];
#elif defined(USE_AVX)
    double temp[4];
    _mm256_storeu_pd(temp, M_12_vec);
    M_12 += temp[0] + temp[1] + temp[2] + temp[3];
#elif defined(USE_AVX512)
    M_12 += _mm512_reduce_add_pd(M_12_vec);
#endif

    M_12 *= 4.0 * local_pi * local_pi * 1e-7 * data.N_1 * data.N_2
            / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2));

    return M_12;
}


FP_TYPE guess_best_inductance_near(
        const CoilCalculationData data,
        const SumPrecisionData precision,
        const FP_TYPE d,
        const FP_TYPE Z_start,
        const FP_TYPE Z_end,
        bool verbose,
        const FP_TYPE r_tol
){
    struct timespec start_time;

    if (verbose) {
        timespec_get(&start_time, TIME_UTC);
    }

    FP_TYPE characteristic_length = data.R_1 > data.R_2 ? data.R_1 : data.R_2;

    const FP_TYPE inv_phi = (sqrt(5.0) - 1.0) / 2.0;  // 1/phi
    const FP_TYPE inv_phi_sq = (3.0 - sqrt(5.0)) / 2.0;  // 1/phi^2
    FP_TYPE Z_low = Z_start;
    FP_TYPE Z_high = Z_end;

    FP_TYPE h = Z_high - Z_low;
    if (h <= 0) {
        fprintf(stderr, "Z_start must be less than Z_end\n");
        return -1.0;
    }

    FP_TYPE Z1 = Z_low + inv_phi_sq * h;
    FP_TYPE Z2 = Z_low + inv_phi * h;

    FP_TYPE F1 = fabs(calculate_mutual_inductance_near_dz(data, precision, d, Z1));
    FP_TYPE F2 = fabs(calculate_mutual_inductance_near_dz(data, precision, d, Z2));

    uint32_t counter = 2;

    while (h > r_tol * characteristic_length) {
        h = inv_phi * h;

        if (F1 < F2) {
            Z_high = Z2;
            Z2 = Z1;
            Z1 = Z_low + inv_phi_sq * h;

            F2 = F1;
            F1 = fabs(calculate_mutual_inductance_near_dz(data, precision, d, Z1));
        } else {
            Z_low = Z1;
            Z1 = Z2;
            Z2 = Z_low + inv_phi * h;

            F1 = F2;
            F2 = fabs(calculate_mutual_inductance_near_dz(data, precision, d, Z2));
        }
        counter++;
    }

    FP_TYPE best_Z = (Z_low + Z_high) / 2;

    FP_TYPE best_value = calculate_mutual_inductance_near(data, precision, d, best_Z);

    if (verbose) {
        struct timespec end_time;
        timespec_get(&end_time, TIME_UTC);

        double interval = (double) (end_time.tv_sec - start_time.tv_sec)
                        + (double) (end_time.tv_nsec - start_time.tv_nsec) / 1e9;

        printf("Optimized with search time =    %g s\n", interval);
        printf("Number of iterations =          %u\n", counter);
        printf("Time per iteration =            %g s\n", interval / (double) counter);
        printf("Best mutual inductance at   Z = %.15g is %.15g\n", best_Z, best_value);
    }

    return best_value;
}

void benchmark_mutual_inductance_near(const SumPrecisionData precision, const uint32_t n_repeats) {
    struct timespec start_time;
    struct timespec end_time;

    CoilCalculationData data = {
        .R_1 = 0.1,
        .r_1 = 0.05,
        .N_1 = 10,
        .L_1 = 0.05,
        .R_2 = 0.1,
        .r_2 = 0.05,
        .N_2 = 10,
        .L_2 = 0.05
    };

    // Volatile to prevent the compiler optimizing out the for loop
    volatile FP_TYPE result;

    timespec_get(&start_time, TIME_UTC);

    for (uint32_t i = 0; i < n_repeats; ++i) {
        result = calculate_mutual_inductance_near(data, precision, (FP_TYPE) i * 0.0001, 1.0);
    }

    timespec_get(&end_time, TIME_UTC);

    double interval = (double) (end_time.tv_sec - start_time.tv_sec)
                    + (double) (end_time.tv_nsec - start_time.tv_nsec) * 1e-9;

    printf("\nBenchmarked mutual inductance near with %u repeats\n", n_repeats);
    printf("Precision: L = %u, K = %u, N = %u\n", precision.l_terms, precision.k_terms, precision.n_terms);
    printf("Elapsed time =          %g s\n", interval);
    printf("Time per iteration =    %g s\n", interval / (double) n_repeats);
    printf("Result (printed to prevent compiler optimization) = %.15g\n\n", result);
}


#endif //VECTOR_CASE_INDUCTANCE_NEAR_H