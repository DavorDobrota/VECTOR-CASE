#ifndef VECTOR_CASE_INDUCTANCE_H
#define VECTOR_CASE_INDUCTANCE_H

#include "inductance_near.h"
#include "inductance_far.h"


/**
 * @brief Calculate the mutual inductance between two coils for any configuration. It resorts to
 * finding the best value of Z to calculate the mutual inductance if the convergence criterion
 * is not met.
 *
 * @details The function does not assume correct input and checks for pathological cases where the coils
 * are not properly defined. It also checks if the number of terms in the series expansion is
 * within the allowed limits. Returns -1.0 if the input is invalid.
 *
 * @param data The coil calculation data containing the physical properties of the coils.
 * @param precision The precision data specifying the number of terms in the series expansion.
 * @param d The distance between the coils, d > 0.
 * @param verbose Whether to print the results and performance of the optimization process.
 * @return The calculated mutual inductance.
 */
FP_TYPE calculate_mutual_inductance(
        const CoilCalculationData data,
        const SumPrecisionData precision,
        const FP_TYPE d,
        const bool verbose
) {
    if (data.r_1 <= 0.0 || data.R_1 <= 0.0 || data.L_1 <= 0.0 || data.N_1 <= 0.0 ||
        data.r_2 <= 0.0 || data.R_2 <= 0.0 || data.L_2 <= 0.0 || data.N_2 <= 0.0) {
        fprintf(stderr, "Error: All parameters must be positive\n");
        return -1.0;
    }
    if (data.r_1 >= data.R_1 || data.r_2 >= data.R_2) {
        fprintf(stderr, "Error: The inner radius must be less than the outer radius\n");
        return -1.0;
    }
    if (precision.k_terms == 0 || precision.l_terms == 0 || precision.n_terms == 0) {
        fprintf(stderr, "Error: Number of terms must be positive\n");
        return -1.0;
    }
    if (d <= 0.0) {
        fprintf(stderr, "Error: The distance between the coils must be positive\n");
        return -1.0;
    }

    FP_TYPE R = data.R_1 > data.R_2 ? data.R_1 : data.R_2;
    FP_TYPE L = data.L_1;

    if (verbose && (precision.k_terms != precision.l_terms)) {
        printf("Warning: precision.k_terms == precision.l_terms is recommended in general\n");
    }

    if (data.L_1 <= 2 * R) {
        if (d > 3.0 * R - L) {
            if (precision.k_terms > MAX_TERMS_FAR || precision.l_terms > MAX_TERMS_FAR || precision.n_terms > MAX_TERMS_FAR) {
                fprintf(stderr, "Error: Number of terms exceeds the maximum allowed value of %d\n", MAX_TERMS_FAR);
                return -1.0;
            }
            if (verbose && (precision.n_terms != precision.k_terms)) {
                printf("Warning: precision.n_terms == precision.k_terms is recommended for long distances\n");
            }

            return calculate_mutual_inductance_far(data, precision, d);
        }
        else {
            if (precision.k_terms > MAX_TERMS_NEAR || precision.l_terms > MAX_TERMS_NEAR || precision.n_terms > MAX_TERMS_NEAR) {
                fprintf(stderr, "Error: Number of terms exceeds the maximum allowed value of %d\n", MAX_TERMS_NEAR);
                return -1.0;
            }
            if (verbose && (precision.n_terms != 2 * precision.k_terms)) {
                printf("Warning: precision.n_terms == 2 * precision.k_terms is recommended in general\n");
            }

            return guess_best_inductance_near(data, precision, d,
                                              fmax(-0.5 * L, R - d), 6.0 * R - d,
                                              verbose, 1e-6);
        }
    }
    else {
        if (d > data.L_1) {
            if (precision.k_terms > MAX_TERMS_FAR || precision.l_terms > MAX_TERMS_FAR || precision.n_terms > MAX_TERMS_FAR) {
                fprintf(stderr, "Error: Number of terms exceeds the maximum allowed value of %d\n", MAX_TERMS_FAR);
                return -1.0;
            }
            if (verbose && precision.n_terms != precision.k_terms) {
                printf("Warning: precision.n_terms == precision.k_terms is recommended for long distances\n");
            }

            return calculate_mutual_inductance_far(data, precision, d);
        }
        else {
            if (precision.k_terms > MAX_TERMS_NEAR || precision.l_terms > MAX_TERMS_NEAR || precision.n_terms > MAX_TERMS_NEAR) {
                fprintf(stderr, "Error: Number of terms exceeds the maximum allowed value of %d\n", MAX_TERMS_NEAR);
                return -1.0;
            }
            if (verbose && (precision.n_terms != 2 * precision.k_terms)) {
                printf("Warning: precision.n_terms == 2 * precision.k_terms is recommended in general\n");
            }

            return guess_best_inductance_near(data, precision, d,
                                              fmax(-0.5 * L, -d), 2.0 * L - d,
                                              verbose, 1e-6);
        }
    }
}


/**
 * @brief A version of calculate_mutual_inductance that takes raw input values. Meant as a
 * minimal interface for the mutual inductance calculation.
 *
 * @param r_1 The inner radius of the first coil, r_1 > 0.
 * @param R_1 The outer radius of the first coil, R_1 > 0.
 * @param L_1 The length of the first coil, L_1 > 0.
 * @param N_1 The number of turns in the first coil, N_1 > 0.
 * @param r_2 The inner radius of the second coil, r_2 > 0.
 * @param R_2 The outer radius of the second coil, R_2 > 0.
 * @param L_2 The length of the second coil, L_2 > 0.
 * @param N_2 The number of turns in the second coil, N_2 > 0.
 * @param d The distance between the coils, d > 0.
 * @param k_terms The number of terms in the series expansion for the k index.
 * @param l_terms The number of terms in the series expansion for the l index.
 * @param n_terms The number of terms in the series expansion for the n index.
 * @param verbose Whether to print the results and performance of the optimization process.
 * @return The calculated mutual inductance.
 */
FP_TYPE calculate_mutual_inductance_raw(
        const FP_TYPE r_1,
        const FP_TYPE R_1,
        const FP_TYPE L_1,
        const FP_TYPE N_1,
        const FP_TYPE r_2,
        const FP_TYPE R_2,
        const FP_TYPE L_2,
        const FP_TYPE N_2,
        const FP_TYPE d,
        const uint32_t k_terms,
        const uint32_t l_terms,
        const uint32_t n_terms,
        const bool verbose
) {
    const CoilCalculationData data = {r_1, R_1, L_1, N_1, r_2, R_2, L_2, N_2};

    const SumPrecisionData precision = {k_terms, l_terms, n_terms};

    return calculate_mutual_inductance(data, precision, d, verbose);
}

#endif //VECTOR_CASE_INDUCTANCE_H