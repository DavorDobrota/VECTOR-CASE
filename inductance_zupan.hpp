#ifndef VECTOR_CASE_INDUCTANCE_ZUPAN_HPP
#define VECTOR_CASE_INDUCTANCE_ZUPAN_HPP

#include <cmath>

#include "structs.h"
#include "gauss_legendre_lookup.h"


FP_TYPE calculate_kernel_zupan(
        const FP_TYPE r,
        const FP_TYPE R,
        const FP_TYPE z,
        const FP_TYPE Z,
        const FP_TYPE phi
) {
    FP_TYPE cos_phi = std::cos(phi);
    FP_TYPE sin_phi = std::sin(phi);

    FP_TYPE sin_2_phi = std::sin(2.0 * phi);
    FP_TYPE cos_2_phi = std::cos(2.0 * phi);
    FP_TYPE cos_3_phi = std::cos(3.0 * phi);

    FP_TYPE Y = Z - z;
    FP_TYPE X = std::sqrt(r * r + R * R + Y * Y - 2.0 * r * R * cos_phi);

    FP_TYPE r_sq = r * r;
    FP_TYPE R_sq = R * R;
    FP_TYPE Y_sq = Y * Y;

    FP_TYPE temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8;

    temp_1 = -10.0 * R * r * (r_sq * r + R_sq * R) * cos_phi * cos_phi;
    temp_2 = 5.0 * (r_sq * r + R_sq * R) * (r_sq + R_sq + 3.0 * Y_sq);
    temp_3 = -2.0 * X * r * R * (r_sq + R_sq);

    FP_TYPE first = temp_1 + cos_phi * (temp_2 + temp_3);

    temp_1 = -2.0 * r_sq * r_sq;
    temp_2 = r_sq * (9.0 * Y_sq - 16.0 * R_sq);
    temp_3 = R_sq * (9.0 * Y_sq - 2.0 * R_sq);
    temp_4 = 6.0 * (r_sq * r_sq + R_sq * R_sq) * cos_2_phi;
    temp_5 = 2.0 * Y_sq * Y_sq / (sin_phi * sin_phi);

    FP_TYPE second = -0.5 * X * (temp_1 + temp_2 + temp_3 + temp_4 + temp_5);

    temp_1 = -2.0 * Y_sq * Y_sq * Y * cos_phi / (sin_phi * sin_phi * sin_phi);
    temp_2 = std::atan2(Y_sq * cos_phi + r * R * sin_phi * sin_phi, Y * X * sin_phi);
    temp_3 = -30.0 * Y * r_sq * R_sq * std::log(X + Y);
    temp_4 = 15.0 * Y * (r_sq * r_sq + R_sq * R_sq) * cos_2_phi * std::log(X + Y);
    temp_5 = 5.0 * r_sq * r * (r_sq - 4.0 * Y_sq) * std::log(R - r * cos_phi + X);
    temp_6 = 5.0 * R_sq * R * (R_sq - 4.0 * Y_sq) * std::log(r - R * cos_phi + X);
    temp_7 = 2.0 * r_sq * r_sq * r * std::log(
            std::sqrt(Y_sq + r_sq * sin_phi * sin_phi) / (R - r * cos_phi + X));
    temp_8 = 2.0 * R_sq * R_sq * R * std::log(
            std::sqrt(Y_sq + R_sq * sin_phi * sin_phi) / (r - R * cos_phi + X));

    FP_TYPE third = 0.5 * (temp_1 * temp_2 + temp_3 + temp_4
                         + cos_phi * (temp_5 + temp_6 + temp_7 + temp_8));

    temp_1 = r_sq * r_sq * r;
    temp_2 = 5.0 * std::log(R - r * cos_phi + X);
    temp_3 = 2.0 * std::log(std::sqrt(Y_sq + r_sq * sin_phi * sin_phi) / (R - r * cos_phi + X));
    temp_4 = R_sq * R_sq * R;
    temp_5 = 5.0 * std::log(r - R * cos_phi + X);
    temp_6 = 2.0 * std::log(std::sqrt(Y_sq + R_sq * sin_phi * sin_phi) / (r - R * cos_phi + X));

    FP_TYPE fourth = -0.5 * cos_3_phi * (temp_1 * (temp_2 + temp_3) + temp_4 * (temp_5 + temp_6));

    temp_1 = r_sq * r_sq;
    temp_2 = 4.0 * std::atan2(Y, r * sin_phi);
    temp_3 = 3.0 * std::atan2(Y * (r * cos_phi - R), r * sin_phi * X);
    temp_4 = R_sq * R_sq;
    temp_5 = 4.0 * std::atan2(Y, R * sin_phi);
    temp_6 = 3.0 * std::atan2(Y * (R * cos_phi - r), R * sin_phi * X);

    FP_TYPE fifth = -2.5 * Y * sin_2_phi * (temp_1 * (temp_2 + temp_3) + temp_4 * (temp_5 + temp_6));

    FP_TYPE result = (FP_TYPE) (1.0 / 60.0) * (first + second + third + fourth + fifth);

    return result;
}

FP_TYPE calculate_mutual_inductance_zupan(
        const CoilCalculationData &data,
        const FP_TYPE d,
        const uint32_t sub_intervals,
        const uint32_t num_gl_points
) {
    FP_TYPE R_1 = data.r_1;
    FP_TYPE R_2 = data.R_1;
    FP_TYPE R_3 = data.r_2;
    FP_TYPE R_4 = data.R_2;

    FP_TYPE Z_1 = 0.0;
    FP_TYPE Z_2 = data.L_1;
    FP_TYPE Z_3 = data.L_1 + d;
    FP_TYPE Z_4 = data.L_1 + d + data.L_2;

    FP_TYPE M_12 = 0.0;

    FP_TYPE interval_len = local_pi / ((FP_TYPE) sub_intervals);

    for (uint32_t i = 0; i < sub_intervals; ++i) {

        FP_TYPE local_offset = i * interval_len;

        for (uint32_t k = 0; k < num_gl_points; ++k) {
            FP_TYPE phi = local_offset + 0.5 * interval_len * (1.0 + gl_positions[num_gl_points - 1][k]);
            FP_TYPE weight = gl_weights[num_gl_points - 1][k] * local_pi * 0.5;
            FP_TYPE cos_phi = std::cos(phi);

            FP_TYPE G = calculate_kernel_zupan(R_1, R_3, Z_1, Z_3, phi)
                      - calculate_kernel_zupan(R_1, R_3, Z_1, Z_4, phi)
                      - calculate_kernel_zupan(R_1, R_3, Z_2, Z_3, phi)
                      + calculate_kernel_zupan(R_1, R_3, Z_2, Z_4, phi)
                      - calculate_kernel_zupan(R_1, R_4, Z_1, Z_3, phi)
                      + calculate_kernel_zupan(R_1, R_4, Z_1, Z_4, phi)
                      + calculate_kernel_zupan(R_1, R_4, Z_2, Z_3, phi)
                      - calculate_kernel_zupan(R_1, R_4, Z_2, Z_4, phi)
                      - calculate_kernel_zupan(R_2, R_3, Z_1, Z_3, phi)
                      + calculate_kernel_zupan(R_2, R_3, Z_1, Z_4, phi)
                      + calculate_kernel_zupan(R_2, R_3, Z_2, Z_3, phi)
                      - calculate_kernel_zupan(R_2, R_3, Z_2, Z_4, phi)
                      + calculate_kernel_zupan(R_2, R_4, Z_1, Z_3, phi)
                      - calculate_kernel_zupan(R_2, R_4, Z_1, Z_4, phi)
                      - calculate_kernel_zupan(R_2, R_4, Z_2, Z_3, phi)
                      + calculate_kernel_zupan(R_2, R_4, Z_2, Z_4, phi);

            M_12 += weight * cos_phi * G;
        }
    }

    M_12 *= (FP_TYPE) (4.0e-7) * local_pi * data.N_1 * data.N_2
          / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2) * (FP_TYPE) sub_intervals);

    return M_12;
}

#endif //VECTOR_CASE_INDUCTANCE_ZUPAN_HPP
