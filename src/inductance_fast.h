#ifndef VECTOR_CASE_INDUCTANCE_FAST_H
#define VECTOR_CASE_INDUCTANCE_FAST_H

#include "settings.h"

#if defined(USE_AVX) || defined(USE_AVX512)
#include "immintrin.h"
#endif

#include "structs.h"
#include "fast_lookup_tables.h"

const double local_pi_double = 3.14159265358979323846;
const float local_pi_float = 3.14159265358979323846f;


double calculate_mutual_inductance_near_float(
        const double N_1,
        const double L_1,
        const double R_1,
        const double r_1,
        const double N_2,
        const double L_2,
        const double R_2,
        const double r_2,
        const double d,
        const double Z
) {
    const double denom_1 = 1.0 / (Z + L_1 + d);
    const double denom_2 = 1.0 / (Z + L_1 + L_2 + d);

    const double R_1_sq = R_1 * R_1;
    const double r_1_sq = r_1 * r_1;

    const double R_2_sq = R_2 * R_2;
    const double r_2_sq = r_2 * r_2;

    const double denom_1_sq = denom_1 * denom_1;
    const double denom_2_sq = denom_2 * denom_2;

    const double L_1_plus_Z = L_1 + Z;

    double R_1_sub_r_1_arr[8];

    double temp_R_1 = R_1_sq * R_1;
    double temp_r_1 = r_1_sq * r_1;

    R_1_sub_r_1_arr[0] = temp_R_1 - temp_r_1;
    for (int i = 1; i < 8; i++) {
        temp_R_1 *= R_1_sq;
        temp_r_1 *= r_1_sq;

        R_1_sub_r_1_arr[i] = temp_R_1 - temp_r_1;
    }

    double R_2_sub_r_2_arr[8];

    double temp_R_2 = R_2_sq * R_2;
    double temp_r_2 = r_2_sq * r_2;

    R_2_sub_r_2_arr[0] = temp_R_2 - temp_r_2;
    for (int i = 1; i < 8; i++) {
        temp_R_2 *= R_2_sq;
        temp_r_2 *= r_2_sq;

        R_2_sub_r_2_arr[i] = temp_R_2 - temp_r_2;
    }

    double denom_1_arr[8];
    double denom_2_arr[8];

    double temp_denom_1 = 1.0;
    double temp_denom_2 = 1.0;

    denom_1_arr[0] = temp_denom_1;
    denom_2_arr[0] = temp_denom_2;

    for (int i = 1; i < 8; i++) {
        temp_denom_1 *= denom_1_sq;
        temp_denom_2 *= denom_2_sq;

        denom_1_arr[i] = temp_denom_1;
        denom_2_arr[i] = temp_denom_2;
    }

    double inner_denom_1_arr[16];
    double inner_denom_2_arr[16];
    double L_1_plus_Z_sub_Z_arr[16];

    double temp_inner_denom_1 = 1.0;
    double temp_inner_denom_2 = 1.0;
    double temp_L_1_plus_Z = L_1_plus_Z;
    double temp_Z = Z;

    inner_denom_1_arr[0] = temp_inner_denom_1;
    inner_denom_2_arr[0] = temp_inner_denom_2;
    L_1_plus_Z_sub_Z_arr[0] = temp_L_1_plus_Z - temp_Z;

    for (int i = 1; i < 16; i++) {
        temp_inner_denom_1 *= denom_1;
        temp_inner_denom_2 *= denom_2;
        temp_L_1_plus_Z *= L_1_plus_Z;
        temp_Z *= Z;

        inner_denom_1_arr[i] = temp_inner_denom_1;
        inner_denom_2_arr[i] = temp_inner_denom_2;
        L_1_plus_Z_sub_Z_arr[i] = temp_L_1_plus_Z - temp_Z;
    }

    double M_12 = 0.0;

#if defined(USE_AVX)
    __m256d M_12_vect = _mm256_set1_pd(0.0);

    for (int l = 0; l < 8; ++l) {
        __m256d R_2_sub_r_2 = _mm256_set1_pd(R_2_sub_r_2_arr[l]);

        for (int k = 0; k < 8; ++k) {
            __m256d R_1_sub_r_1 = _mm256_set1_pd(R_1_sub_r_1_arr[k]);

            __m256d loop_denom_1 = _mm256_set1_pd(denom_1_arr[k] * denom_1_arr[l] * denom_1_sq);
            __m256d loop_denom_2 = _mm256_set1_pd(denom_2_arr[k] * denom_2_arr[l] * denom_2_sq);

            for (int n = 0; n < 16; n += 4) {
                __m256d L_1_plus_Z_sub_Z = _mm256_loadu_pd(&L_1_plus_Z_sub_Z_arr[n]);

                __m256d inner_loop_denom_1 = _mm256_loadu_pd(&inner_denom_1_arr[n]);
                __m256d inner_loop_denom_2 = _mm256_loadu_pd(&inner_denom_2_arr[n]);

                inner_loop_denom_1 = _mm256_mul_pd(inner_loop_denom_1, loop_denom_1);
                inner_loop_denom_2 = _mm256_mul_pd(inner_loop_denom_2, loop_denom_2);

                __m256d element = _mm256_mul_pd(
                    _mm256_mul_pd(R_1_sub_r_1, R_2_sub_r_2),
                    _mm256_mul_pd(
                        L_1_plus_Z_sub_Z,
                        _mm256_sub_pd(inner_loop_denom_1, inner_loop_denom_2)
                    )
                );

                __m256d lookup = _mm256_loadu_pd(&near_table_double_fast[k][l][n]);

                M_12_vect = _mm256_add_pd(M_12_vect, _mm256_mul_pd(element, lookup));
            }
        }
    }

    double temp[4];
    _mm256_storeu_pd(temp, M_12_vect);
    M_12 += temp[0] + temp[1] + temp[2] + temp[3];
#elif defined(USE_AVX512)
    __m512d M_12_vec = _mm512_set1_pd(0.0);

    for (int l = 0; l < 8; ++l) {
        __m512d R_2_sub_r_2 = _mm512_set1_pd(R_2_sub_r_2_arr[l]);

        for (int k = 0; k < 8; ++k) {
            __m512d R_1_sub_r_1 = _mm512_set1_pd(R_1_sub_r_1_arr[k]);

            __m512d loop_denom_1 = _mm512_set1_pd(denom_1_arr[k] * denom_1_arr[l] * denom_1_sq);
            __m512d loop_denom_2 = _mm512_set1_pd(denom_2_arr[k] * denom_2_arr[l] * denom_2_sq);

            for (int n = 0; n < 16; n += 8) {
                __m512d L_1_plus_Z_sub_Z = _mm512_loadu_pd(&L_1_plus_Z_sub_Z_arr[n]);

                __m512d inner_loop_denom_1 = _mm512_loadu_pd(&inner_denom_1_arr[n]);
                __m512d inner_loop_denom_2 = _mm512_loadu_pd(&inner_denom_2_arr[n]);

                inner_loop_denom_1 = _mm512_mul_pd(inner_loop_denom_1, loop_denom_1);
                inner_loop_denom_2 = _mm512_mul_pd(inner_loop_denom_2, loop_denom_2);

                __m512d element = _mm512_mul_pd(
                    _mm512_mul_pd(R_1_sub_r_1, R_2_sub_r_2),
                    _mm512_mul_pd(
                        L_1_plus_Z_sub_Z,
                        _mm512_sub_pd(inner_loop_denom_1, inner_loop_denom_2)
                    )
                );

                __m512d lookup = _mm512_loadu_pd(&near_table_double_fast[k][l][n]);

                M_12_vec = _mm512_add_pd(M_12_vec, _mm512_mul_pd(element, lookup));
            }
        }
    }
    M_12 += _mm512_reduce_add_pd(M_12_vec);
#else
    for (int l = 0; l < 8; ++l) {
        for (int k = 0; k < 8; ++k) {

            double loop_denom_1 = denom_1_arr[k] * denom_1_arr[l] * denom_1_sq;
            double loop_denom_2 = denom_2_arr[k] * denom_2_arr[l] * denom_2_sq;

            for (int n = 0; n < 16; ++n) {

                double inner_loop_denom_1 = loop_denom_1 * inner_denom_1_arr[n];
                double inner_loop_denom_2 = loop_denom_2 * inner_denom_2_arr[n];

                M_12 += near_table_double_fast[k][l][n]
                        * R_1_sub_r_1_arr[k]
                        * R_2_sub_r_2_arr[l]
                        * L_1_plus_Z_sub_Z_arr[n]
                        * (inner_loop_denom_1 - inner_loop_denom_2);
            }
        }
    }
#endif

    M_12 *= 4.0e-7 * local_pi_double * local_pi_double * N_1 * N_2
          / (L_1 * L_2 * (R_1 - r_1) * (R_2 - r_2));

    return M_12;
}

float calculate_mutual_inductance_float_fast(
        const float N_1,
        const float L_1,
        const float R_1,
        const float r_1,
        const float N_2,
        const float L_2,
        const float R_2,
        const float r_2,
        const float d,
        const float Z
) {
    const float denom_1 = 1.0f / (Z + L_1 + d);
    const float denom_2 = 1.0f / (Z + L_1 + L_2 + d);

    const float R_1_sq = R_1 * R_1;
    const float r_1_sq = r_1 * r_1;

    const float R_2_sq = R_2 * R_2;
    const float r_2_sq = r_2 * r_2;

    const float denom_1_sq = denom_1 * denom_1;
    const float denom_2_sq = denom_2 * denom_2;

    const float L_1_plus_Z = L_1 + Z;

    float R_1_sub_r_1_arr[8];
    float R_2_sub_r_2_arr[8];
    float L_1_plus_Z_sub_Z_arr[8];

    float temp_R_1 = R_1_sq * R_1;
    float temp_r_1 = r_1_sq * r_1;

    R_1_sub_r_1_arr[0] = temp_R_1 - temp_r_1;
    for (int i = 1; i < 8; i++) {
        temp_R_1 *= R_1_sq;
        temp_r_1 *= r_1_sq;

        R_1_sub_r_1_arr[i] = temp_R_1 - temp_r_1;
    }

    float temp_R_2 = R_2_sq * R_2;
    float temp_r_2 = r_2_sq * r_2;

    R_2_sub_r_2_arr[0] = temp_R_2 - temp_r_2;
    for (int i = 1; i < 8; i++) {
        temp_R_2 *= R_2_sq;
        temp_r_2 *= r_2_sq;

        R_2_sub_r_2_arr[i] = temp_R_2 - temp_r_2;
    }

    float temp_L_1_plus_Z = L_1_plus_Z;
    float temp_Z = Z;

    L_1_plus_Z_sub_Z_arr[0] = temp_L_1_plus_Z - temp_Z;
    for (int i = 1; i < 8; i++) {
        temp_L_1_plus_Z *= L_1_plus_Z;
        temp_Z *= Z;

        L_1_plus_Z_sub_Z_arr[i] = temp_L_1_plus_Z - temp_Z;
    }

    float denom_1_arr[8];
    float denom_2_arr[8];
    float inner_denom_1_arr[8];
    float inner_denom_2_arr[8];

    float temp_denom_1 = 1.0f;
    float temp_denom_2 = 1.0f;
    float temp_inner_denom_1 = 1.0f;
    float temp_inner_denom_2 = 1.0f;

    denom_1_arr[0] = temp_denom_1;
    denom_2_arr[0] = temp_denom_2;
    inner_denom_1_arr[0] = temp_inner_denom_1;
    inner_denom_2_arr[0] = temp_inner_denom_2;

    for (int i = 1; i < 8; i++) {
        temp_denom_1 *= denom_1_sq;
        temp_denom_2 *= denom_2_sq;
        temp_inner_denom_1 *= denom_1;
        temp_inner_denom_2 *= denom_2;

        denom_1_arr[i] = temp_denom_1;
        denom_2_arr[i] = temp_denom_2;
        inner_denom_1_arr[i] = temp_inner_denom_1;
        inner_denom_2_arr[i] = temp_inner_denom_2;
    }

    float M_12 = 0.0f;

    for (int l = 0; l < 8; ++l) {
        for (int k = 0; k < 8; ++k) {

            float loop_denom_1 = denom_1_arr[k] * denom_1_arr[l] * denom_1_sq;
            float loop_denom_2 = denom_2_arr[k] * denom_2_arr[l] * denom_2_sq;

            for (int n = 0; n < 8; ++n) {

                loop_denom_1 = loop_denom_1 * inner_denom_1_arr[n];
                loop_denom_2 = loop_denom_2 * inner_denom_2_arr[n];

                M_12 += near_table_float_fast[k][l][n]
                        * R_1_sub_r_1_arr[k]
                        * R_2_sub_r_2_arr[l]
                        * L_1_plus_Z_sub_Z_arr[n]
                        * (loop_denom_1 - loop_denom_2);
            }
        }
    }

    M_12 *= 4.0e-7f * local_pi_float * local_pi_float * N_1 * N_2
            / (L_1 * L_2 * (R_1 - r_1) * (R_2 - r_2));

    return M_12;
}

#if defined(USE_AVX)
float calculate_mutual_inductance_float_fast_avx(
        const float N_1,
        const float L_1,
        const float R_1,
        const float r_1,
        const float N_2,
        const float L_2,
        const float R_2,
        const float r_2,
        const float d,
        const float Z
) {
    const float denom_1 = 1.0f / (Z + L_1 + d);
    const float denom_2 = 1.0f / (Z + L_1 + L_2 + d);

    const float R_1_sq = R_1 * R_1;
    const float r_1_sq = r_1 * r_1;

    const float R_2_sq = R_2 * R_2;
    const float r_2_sq = r_2 * r_2;

    const float denom_1_sq = denom_1 * denom_1;
    const float denom_2_sq = denom_2 * denom_2;

    const float L_1_plus_Z = L_1 + Z;

    float R_1_sub_r_1_arr[8];
    float R_2_sub_r_2_arr[8];
    float L_1_plus_Z_sub_Z_arr[8];

    float temp_R_1 = R_1_sq * R_1;
    float temp_r_1 = r_1_sq * r_1;

    R_1_sub_r_1_arr[0] = temp_R_1 - temp_r_1;
    for (int i = 1; i < 8; i++) {
        temp_R_1 *= R_1_sq;
        temp_r_1 *= r_1_sq;

        R_1_sub_r_1_arr[i] = temp_R_1 - temp_r_1;
    }

    float temp_R_2 = R_2_sq * R_2;
    float temp_r_2 = r_2_sq * r_2;

    R_2_sub_r_2_arr[0] = temp_R_2 - temp_r_2;
    for (int i = 1; i < 8; i++) {
        temp_R_2 *= R_2_sq;
        temp_r_2 *= r_2_sq;

        R_2_sub_r_2_arr[i] = temp_R_2 - temp_r_2;
    }

    float temp_L_1_plus_Z = L_1_plus_Z;
    float temp_Z = Z;

    L_1_plus_Z_sub_Z_arr[0] = temp_L_1_plus_Z - temp_Z;
    for (int i = 1; i < 8; i++) {
        temp_L_1_plus_Z *= L_1_plus_Z;
        temp_Z *= Z;

        L_1_plus_Z_sub_Z_arr[i] = temp_L_1_plus_Z - temp_Z;
    }

    float denom_1_arr[8];
    float denom_2_arr[8];
    float inner_denom_1_arr[8];
    float inner_denom_2_arr[8];

    float temp_denom_1 = 1.0;
    float temp_denom_2 = 1.0;
    float temp_inner_denom_1 = 1.0;
    float temp_inner_denom_2 = 1.0;

    denom_1_arr[0] = temp_denom_1;
    denom_2_arr[0] = temp_denom_2;
    inner_denom_1_arr[0] = temp_inner_denom_1;
    inner_denom_2_arr[0] = temp_inner_denom_2;

    for (int i = 1; i < 8; i++) {
        temp_denom_1 *= denom_1_sq;
        temp_denom_2 *= denom_2_sq;
        temp_inner_denom_1 *= denom_1;
        temp_inner_denom_2 *= denom_2;

        denom_1_arr[i] = temp_denom_1;
        denom_2_arr[i] = temp_denom_2;
        inner_denom_1_arr[i] = temp_inner_denom_1;
        inner_denom_2_arr[i] = temp_inner_denom_2;
    }

    __m256 M_12_vect = _mm256_set1_ps(0.0f);

    for (int l = 0; l < 8; ++l) {

        __m256 R_1_sub_r_1 = _mm256_set1_ps(R_1_sub_r_1_arr[l]);

        for (int k = 0; k < 8; ++k) {

            __m256 loop_denom_1 = _mm256_set1_ps(denom_1_arr[k] * denom_1_arr[l] * denom_1_sq);
            __m256 loop_denom_2 = _mm256_set1_ps(denom_2_arr[k] * denom_2_arr[l] * denom_2_sq);
            __m256 R_2_sub_r_2 = _mm256_set1_ps(R_2_sub_r_2_arr[k]);

            __m256 element = _mm256_mul_ps(R_1_sub_r_1, R_2_sub_r_2);

            __m256 vec_denom_1 = _mm256_loadu_ps(&inner_denom_1_arr[0]);
            __m256 vec_denom_2 = _mm256_loadu_ps(&inner_denom_2_arr[0]);

            loop_denom_1 = _mm256_mul_ps(loop_denom_1, vec_denom_1);
            loop_denom_2 = _mm256_mul_ps(loop_denom_2, vec_denom_2);

            __m256 L_1_plus_Z_sub_Z = _mm256_loadu_ps(&L_1_plus_Z_sub_Z_arr[0]);

            element = _mm256_mul_ps(element, L_1_plus_Z_sub_Z);
            element = _mm256_mul_ps(element, _mm256_sub_ps(loop_denom_1, loop_denom_2));

            __m256 lookup = _mm256_loadu_ps(&near_table_float_fast[k][l][0]);
            element = _mm256_mul_ps(element, lookup);

            M_12_vect = _mm256_add_ps(M_12_vect, element);
        }
    }

    float M_12 = M_12_vect[0] + M_12_vect[1] + M_12_vect[2] + M_12_vect[3]
               + M_12_vect[4] + M_12_vect[5] + M_12_vect[6] + M_12_vect[7];

    M_12 *= 4.0e-7f * local_pi_float * local_pi_float * N_1 * N_2
            / (L_1 * L_2 * (R_1 - r_1) * (R_2 - r_2));

    return M_12;
}
#endif

#endif //VECTOR_CASE_INDUCTANCE_FAST_H