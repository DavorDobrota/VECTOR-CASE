#ifndef VECTOR_CASE_INDUCTANCE_RADIAL_H
#define VECTOR_CASE_INDUCTANCE_RADIAL_H

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "factorial_lookup.h"
#include "structs.h"
#include "timing.h"
#include "hypergeometric.h"


#define MAX_NUM_TERMS 85
#define N_INTERPOLATION_VALUES 20


/**
 * @brief The zeroth term function for the radial series expansion.
 * Includes compensation for small parameters to avoid loss of precision
 *
 * @param x Evaluation parameter for the function
 * @param r Inner radius
 * @param R Outer radius
 * @return Function value at given parameters
 */
static FP_TYPE f_0_func(const FP_TYPE x, const FP_TYPE r, const FP_TYPE R) {
    FP_TYPE result = 0.0;

    if (R * 1e-12 > r) {
        return -1.0;
    }

    if (fabs(x) < 1e-2 * r) {
        const FP_TYPE x_sq = x * x;
        const FP_TYPE R_sq = R * R;
        const FP_TYPE r_sq = r * r;

        const FP_TYPE inv_R_2_sq = 1.0 / R_sq;
        const FP_TYPE inv_r_2_sq = 1.0 / r_sq;

        result += log(r / R) * x_sq;

        FP_TYPE temp_x = x_sq * x_sq;
        FP_TYPE temp_inv_R = inv_R_2_sq * inv_R_2_sq;
        FP_TYPE temp_inv_r = inv_r_2_sq * inv_r_2_sq;
        result += 0.25 * (temp_inv_r - temp_inv_R) * temp_x;

        temp_x *= x_sq;
        temp_inv_R *= inv_R_2_sq;
        temp_inv_r *= inv_r_2_sq;
        result -= 0.09375 * (temp_inv_r - temp_inv_R) * temp_x;

        temp_inv_R *= inv_R_2_sq;
        temp_inv_r *= inv_r_2_sq;
        temp_x *= x_sq;
        result += 0.05208333333333333 * (temp_inv_r - temp_inv_R) * temp_x;

        result += r * sqrt(r_sq + x_sq) - R * sqrt(R_sq + x_sq);
    }
    else {
        const FP_TYPE sqrt_R_sq_add_x_sq = sqrt(R * R + x * x);
        const FP_TYPE sqrt_r_sq_add_x_sq = sqrt(r * r + x * x);

        const FP_TYPE term_1 = x * x * atanh(R / sqrt_R_sq_add_x_sq) + R * sqrt_R_sq_add_x_sq;
        const FP_TYPE term_2 = x * x * atanh(r / sqrt_r_sq_add_x_sq) + r * sqrt_r_sq_add_x_sq;

        result += term_2 - term_1;
    }

    return result;
}

/**
 * @brief The first term function for the radial series expansion.
 * Includes compensation for small parameters to avoid loss of precision
 *
 * @param x Evaluation parameter for the function
 * @param r Inner radius
 * @param R Outer radius
 * @return Function value at given parameters
 */
static FP_TYPE f_1_func(const FP_TYPE x, const FP_TYPE r, const FP_TYPE R) {
    FP_TYPE result = 0.0;

    if (R * 1e-12 > r) {
        return -1.0;
    }

    if (fabs(x) < 1e-2 * r) {
        const FP_TYPE x_sq = x * x;
        const FP_TYPE R_sq = R * R;
        const FP_TYPE r_sq = r * r;

        const FP_TYPE inv_R_2_sq = 1.0 / R_sq;
        const FP_TYPE inv_r_2_sq = 1.0 / r_sq;

        result += log(R / r);

        FP_TYPE temp_x = x_sq;
        FP_TYPE temp_inv_R = inv_R_2_sq;
        FP_TYPE temp_inv_r = inv_r_2_sq;
        result -= 0.25 * (temp_inv_r - temp_inv_R) * temp_x;

        temp_x *= x_sq;
        temp_inv_R *= inv_R_2_sq;
        temp_inv_r *= inv_r_2_sq;
        result += 0.09375 * (temp_inv_r - temp_inv_R) * temp_x;

        temp_x *= x_sq;
        temp_inv_R *= inv_R_2_sq;
        temp_inv_r *= inv_r_2_sq;
        result -= 0.052083333333333 * (temp_inv_r - temp_inv_R) * temp_x;

        temp_x *= x_sq;
        temp_inv_R *= inv_R_2_sq;
        temp_inv_r *= inv_r_2_sq;
        result += 0.0341796875 * (temp_inv_r - temp_inv_R) * temp_x;

        result += r / sqrt(r_sq + x_sq) - R / sqrt(R_sq + x_sq);
    }
    else {
        const FP_TYPE sqrt_R_sq_add_x_sq = sqrt(R * R + x * x);
        const FP_TYPE sqrt_r_sq_add_x_sq = sqrt(r * r + x * x);

        const FP_TYPE inv_sqrt_R_2_sq_add_x_sq = 1.0 / sqrt_R_sq_add_x_sq;
        const FP_TYPE inv_sqrt_r_2_sq_add_x_sq = 1.0 / sqrt_r_sq_add_x_sq;

        const FP_TYPE term_1 = atanh(R * inv_sqrt_R_2_sq_add_x_sq) - R * inv_sqrt_R_2_sq_add_x_sq;
        const FP_TYPE term_2 = atanh(r * inv_sqrt_r_2_sq_add_x_sq) - r * inv_sqrt_r_2_sq_add_x_sq;

        result += term_1 - term_2;
    }

    return result;
}

/**
 * @brief The n-th term function of the radial series expansion for coaxial mutual inductance.
 *
 * @param x Evaluation parameter for the function
 * @param r_2 Inner radius of the secondary coil
 * @param R_2 Outer radius of the secondary coil
 * @param n Term index, n >= 2
 * @return Function value at given parameters
 */
static FP_TYPE f_n_func(const FP_TYPE x, const FP_TYPE r_2, const FP_TYPE R_2, const uint32_t n) {
    const FP_TYPE argument_1 = -x * x / (r_2 * r_2);
    const FP_TYPE argument_2 = -x * x / (R_2 * R_2);

    assert(n >= 2);

    const FP_TYPE r_2_pow = pow(r_2, 2 * n - 2);
    const FP_TYPE R_2_pow = pow(R_2, 2 * n - 2);

    const FP_TYPE term_1 = hypergeometric3F2(n + 0.5, n - 0.5, n - 1.0, 0.5, n, argument_1) / r_2_pow;
    const FP_TYPE term_2 = hypergeometric3F2(n + 0.5, n - 0.5, n - 1.0, 0.5, n, argument_2) / R_2_pow;

    return term_1 - term_2;
}

/**
 * @brief Calculate the mutual inductance between two coils which are radially separated.
 *
 * This function calculates the mutual inductance between two coils using a single series expansion
 * based on hypergeometric functions. Depending on the geometric properties, the 3F2 hypergeometric
 * function is evaluated either in terms of its series representation or via numerical integration.
 * Therefore, the execution time heavily depends on the coil geometry.
 *
 * The number of terms in the sum is not fixed; it is determined dynamically based on the relative
 * tolerance provided. If the series does not converge, i.e., the error does not fall below the
 * provided tolerance within MAX_NUM_TERMS, a truncated sum estimation is performed to approximate
 * the result (if the last N_INTERPOLATION_VALUES do not change sign).
 *
 * @param data The coil calculation data containing the physical properties of the coils.
 * @param offset_1 The offset of the first coil along the z-axis.
 * @param offset_2 The offset of the second coil along the z-axis.
 * @param relative_tol The relative tolerance for convergence of the series expansion.
 * @param verbose Whether to print convergence information.
 * @return Mutual inductance of radially separated coils.
 */
static FP_TYPE calculate_mutual_inductance_radial(
        const CoilCalculationData data,
        const FP_TYPE offset_1,
        const FP_TYPE offset_2,
        const FP_TYPE relative_tol,
        const bool verbose
) {
    FP_TYPE M_12 = 0.0;

    const FP_TYPE Z_1 = offset_1;
    const FP_TYPE Z_2 = Z_1 + data.L_1;
    const FP_TYPE Z_3 = offset_2;
    const FP_TYPE Z_4 = Z_3 + data.L_2;

    const FP_TYPE R_1_sq = data.R_1 * data.R_1;
    const FP_TYPE r_1_sq = data.r_1 * data.r_1;

    FP_TYPE R_1_pow = R_1_sq * data.R_1;
    FP_TYPE r_1_pow = r_1_sq * data.r_1;

    FP_TYPE terms[MAX_NUM_TERMS] = {0.0};

    const FP_TYPE zeroth_term = f_0_func(Z_4 - Z_2, data.r_2, data.R_2) + f_0_func(Z_3 - Z_1, data.r_2, data.R_2)
                              - f_0_func(Z_4 - Z_1, data.r_2, data.R_2) - f_0_func(Z_3 - Z_2, data.r_2, data.R_2);

    M_12 += (R_1_pow - r_1_pow) * zeroth_term / 12.0;
    terms[0] = (R_1_pow - r_1_pow) * zeroth_term / 12.0;

    R_1_pow *= R_1_sq;
    r_1_pow *= r_1_sq;

    const FP_TYPE first_term = f_1_func(Z_4 - Z_2, data.r_2, data.R_2) + f_1_func(Z_3 - Z_1, data.r_2, data.R_2)
                             - f_1_func(Z_4 - Z_1, data.r_2, data.R_2) - f_1_func(Z_3 - Z_2, data.r_2, data.R_2);

    M_12 += (R_1_pow - r_1_pow) * first_term * 0.0125; // 1/80
    terms[1] = (R_1_pow - r_1_pow) * first_term * 0.0125; // 1/80

    bool converged = false;

    for (uint32_t n = 2; n < MAX_NUM_TERMS; ++n) {
        const FP_TYPE two_to_two_n = pow(2.0, 2 * n);

        const FP_TYPE first_frac = factorial_array[2 * n - 3] / (factorial_array[n - 1] * factorial_array[n] * two_to_two_n);
        const FP_TYPE second_frac = factorial_array[2 * n] / (factorial_array[n] * factorial_array[n + 1] * two_to_two_n);
        const FP_TYPE third_frac = 1.0 / ((double) n * 2.0 + 3.0);

        const FP_TYPE body = f_n_func(Z_4 - Z_2, data.r_2, data.R_2, n) + f_n_func(Z_3 - Z_1, data.r_2, data.R_2, n)
                           - f_n_func(Z_4 - Z_1, data.r_2, data.R_2, n) - f_n_func(Z_3 - Z_2, data.r_2, data.R_2, n);

        R_1_pow *= R_1_sq;
        r_1_pow *= r_1_sq;

        const FP_TYPE term = (R_1_pow - r_1_pow) * first_frac * second_frac * third_frac * body;
        terms[n] = term;

        M_12 += term;

        // This is "regular" convergence, because the term did not underflow, no need for truncated sum estimation.
        if (fabs(term) < relative_tol * M_12) {
            converged = true;
            if (verbose)
                printf("Converged after %u terms with residual of (relative) magnitude %.3e\n", n, fabs(term) / M_12);
            break;
        }
    }

    if (!converged) {
        if (verbose) {
            printf(
                "Series did not converge after %d terms, residual (relative) magnitude %.3e\n",
                MAX_NUM_TERMS,
                fabs(terms[MAX_NUM_TERMS - 1]) / M_12
            );
        }

        // Extract last N_INTERPOLATION_VALUES terms and their indices
        FP_TYPE x_values[N_INTERPOLATION_VALUES];
        FP_TYPE y_values[N_INTERPOLATION_VALUES];
        int count = 0;
        bool oscillating_sign = false;

        for (uint32_t i = MAX_NUM_TERMS - N_INTERPOLATION_VALUES; i < MAX_NUM_TERMS; i++) {
            if (fabs(terms[i]) > 1e-100) {
                if (terms[i] * terms[i - 1] < 0.0) {
                    oscillating_sign = true; // Sign change detected
                    break;
                }
                x_values[count] = log((FP_TYPE)(i + 1)); // log of n (1-indexed)
                y_values[count] = log(fabs(terms[i]));  // log of absolute term value
                count++;
            }
        }

        if (!oscillating_sign && count > 2) {
            // Simple linear regression on log-log scale
            FP_TYPE sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
            for (int i = 0; i < count; i++) {
                sum_x += x_values[i];
                sum_y += y_values[i];
                sum_xy += x_values[i] * y_values[i];
                sum_xx += x_values[i] * x_values[i];
            }

            // Calculate slope and intercept
            const FP_TYPE slope = (count * sum_xy - sum_x * sum_y) / (count * sum_xx - sum_x * sum_x);
            const FP_TYPE intercept = (sum_y - slope * sum_x) / count;
            // printf("Power law fit: slope = %g, intercept = %g\n", slope, intercept);

            if (slope <= -1.0 && isfinite(slope)) {
                // Start from the next term after MAX_NUM_TERMS
                const FP_TYPE a = exp(intercept);
                const FP_TYPE b = slope;
                const FP_TYPE remaining_sum = a * pow(MAX_NUM_TERMS, b + 1) / (-b - 1);

                // Add the estimated remaining sum to M_12
                if (terms[MAX_NUM_TERMS - 1] < 0.0) {
                    M_12 -= remaining_sum; // Adjust sign if necessary
                }
                else {
                    M_12 += remaining_sum; // Adjust sign if necessary
                }

                if (verbose)
                    printf("Estimated sum remainder using power-law fit;     (relative) magnitude %.3e\n",
                        remaining_sum / M_12
                    );
            }
        }
    }

    M_12 *= (4.0e-7 * LOCAL_M_PI * LOCAL_M_PI * data.N_1 * data.N_2)
          / (data.L_1 * data.L_2 * (data.R_2 - data.r_2) * (data.R_1 - data.r_1));

    return M_12;
}

/** @brief Benchmark the mutual inductance radial function.
 *
 * @param n_repeats The number of times to repeat the calculation
 * @param fast Whether the coil configuration tested admits a series expansion for the 3F2
 * hypergeometric function or numerical integration has to be used
 * @param relative_tol The relative tolerance for the series expansion convergence
 */
static void benchmark_mutual_inductance_radial(const uint32_t n_repeats, const bool fast, const FP_TYPE relative_tol) {
    const CoilCalculationData data = {
        .R_1 = 0.4,
        .r_1 = 0.3,
        .L_1 = fast ? 0.05 : 0.5,
        .N_1 = 10,

        .R_2 = 0.6,
        .r_2 = 0.5,
        .L_2 = fast ? 0.05 : 0.5,
        .N_2 = 10,
    };

    // Volatile to prevent the compiler optimizing out the for loop
    volatile FP_TYPE result = 0.0;

    Timer timer;
    timer_start(&timer);

    for (uint32_t i = 0; i < n_repeats; ++i) {
        result += calculate_mutual_inductance_radial(
            data,
            -0.5 * data.L_1 + (FP_TYPE) i * 0.00001,
            -0.5 * data.L_1 + (FP_TYPE) i * 0.00001,
            relative_tol,
            false
        );
    }

    const double interval = timer_elapsed(&timer);

    printf("\nBenchmarked mutual inductance radial with %u repeats\n", n_repeats);
    printf("Using the series representation: %s\n", fast ? "YES" : "NO");
    printf("Elapsed time =          %g s\n", interval);
    printf("Time per iteration =    %.3e s\n", interval / (double) n_repeats);
    printf("Result (printed to prevent compiler optimization) = %.15g\n\n", result);
}

#endif //VECTOR_CASE_INDUCTANCE_RADIAL_H