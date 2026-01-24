#ifndef VECTOR_CASE_HYPERGEOMETRIC_H
#define VECTOR_CASE_HYPERGEOMETRIC_H

#include <math.h>
#include <stdint.h>

#include "settings.h"
#include "gauss_kronrod_lookup.h"

// Tolerance settings
#define SERIES_TOL 1e-16
#define MAX_RECURSION_DEPTH 5
#define MAX_SERIES_ITER 300


// Series evaluation for {}_2F_1 using its Taylor expansion.
static FP_TYPE hypergeometric2F1_series(const FP_TYPE a, const FP_TYPE b, const FP_TYPE c, const FP_TYPE x) {
    FP_TYPE term = 1.0;
    FP_TYPE sum = term;

    for (uint32_t n = 0; n < MAX_SERIES_ITER; n += 4) {
        FP_TYPE n_fp = n;
        term *= (a + n_fp) * (b + n_fp) / ((c + n_fp) * (n_fp + 1.0)) * x;
        sum += term;

        n_fp += 1;
        term *= (a + n_fp) * (b + n_fp) / ((c + n_fp) * (n_fp + 1.0)) * x;
        sum += term;

        n_fp += 1;
        term *= (a + n_fp) * (b + n_fp) / ((c + n_fp) * (n_fp + 1.0)) * x;
        sum += term;

        n_fp += 1;
        term *= (a + n_fp) * (b + n_fp) / ((c + n_fp) * (n_fp + 1.0)) * x;
        sum += term;

        if (fabs(term) < SERIES_TOL * fabs(sum)) {
            // printf("2F1 series iterations: %u\n", n + 4);
            // printf("Value of 2F1(%f, %f, %f, %f) = %.16g\n", a, b, c, x, sum);
            break;
        }
    }
    return sum;
}

// {}_2F_1 function that uses a transformation if |x| is too large.
static FP_TYPE hypergeometric2F1(const FP_TYPE a, const FP_TYPE b, const FP_TYPE c, const FP_TYPE x) {
    if (x < 0) {
        const FP_TYPE prefactor = pow(1.0 - x, c - a - b);
        const FP_TYPE result = hypergeometric2F1_series(c - a, c - b, c, x);

        // Avoid returning infinity for large negative x
        if (!isfinite(result)) {
            return 0.0;
        }
        return prefactor * result;
    }

    // Use series representation for small |x|
    if (fabs(x) < 0.5 && fabs(x) * a * b < 50.0 * c) {
        return hypergeometric2F1_series(a, b, c, x);
    }

    return -FP_INFINITE;
}

static FP_TYPE hypergeometric_3F2_series(
        const FP_TYPE a_1,
        const FP_TYPE a_2,
        const FP_TYPE a_3,
        const FP_TYPE b_1,
        const FP_TYPE b_2,
        const FP_TYPE z,
        const uint32_t max_num_terms, const FP_TYPE thresh
) {
    // The first term is included
    FP_TYPE term = 1.0;
    FP_TYPE sum = term;

    for (uint32_t n = 0; n < max_num_terms; n += 4) {
        FP_TYPE n_fp = n;
        term *= ((a_1 + n_fp) * (a_2 + n_fp) * (a_3 + n_fp)) / ((b_1 + n_fp) * (b_2 + n_fp) * (n_fp + 1.0)) * z;
        sum += term;

        n_fp += 1;
        term *= ((a_1 + n_fp) * (a_2 + n_fp) * (a_3 + n_fp)) / ((b_1 + n_fp) * (b_2 + n_fp) * (n_fp + 1.0)) * z;
        sum += term;

        n_fp += 1;
        term *= ((a_1 + n_fp) * (a_2 + n_fp) * (a_3 + n_fp)) / ((b_1 + n_fp) * (b_2 + n_fp) * (n_fp + 1.0)) * z;
        sum += term;

        n_fp += 1;
        term *= ((a_1 + n_fp) * (a_2 + n_fp) * (a_3 + n_fp)) / ((b_1 + n_fp) * (b_2 + n_fp) * (n_fp + 1.0)) * z;
        sum += term;

        // Stop if the term becomes negligible
        if (fabs(term) < thresh * fabs(sum)) {
            break;
        }
    }
    return sum;
}

static FP_TYPE integrand_3F2(
        const FP_TYPE t,
        const FP_TYPE a, const FP_TYPE b, const FP_TYPE c, const FP_TYPE d, const FP_TYPE e,
        const FP_TYPE z
) {
    return pow(t, c - 1.0) * pow(1.0 - t, e - c - 1.0) * hypergeometric2F1(a, b, d, t * z);
}

static FP_TYPE hypergeometric3F2_gauss_kronrod_recursive(
    const FP_TYPE A, const FP_TYPE B,
    const FP_TYPE a, const FP_TYPE b, const FP_TYPE c, const FP_TYPE d, const FP_TYPE e,
    const FP_TYPE z,
    const FP_TYPE tol, const uint32_t depth
){
    const FP_TYPE mid = (A + B) / 2.0;
    const FP_TYPE half_length = (B - A) / 2.0;

    // Evaluate f(mid) exactly once:
    const FP_TYPE f_mid = integrand_3F2(mid, a, b, c, d, e, z);

    // Kronrod estimate: start with center contribution
    FP_TYPE kronrod_sum = kronrod_weights_21_10[0] * f_mid;
    FP_TYPE gauss_sum = 0.0;

    // Loop i = 1..10 for the positive Kronrod nodes
    for (int i = 1; i < 11; i++) {
        // scaled abscissa on [A,B]:
        const FP_TYPE x = kronrod_nodes_21_10[i];
        const FP_TYPE abscissa = half_length * x;

        // evaluate at mid ± abscissa
        const FP_TYPE f_left  = integrand_3F2(mid - abscissa, a, b, c, d, e, z);
        const FP_TYPE f_right = integrand_3F2(mid + abscissa, a, b, c, d, e, z);

        // add to Kronrod sum
        kronrod_sum += kronrod_weights_21_10[i] * (f_left + f_right);

        // check if this Kronrod index corresponds to a Gauss node:
        // the positive Gauss indices are exactly 1,3,5,7,9 → i % 2 == 1
        if (i % 2 == 1) {
            const int gi = i / 2;  // 1→0, 3→1, 5→2, 7→3, 9→4
            gauss_sum += gauss_weights_21_10[gi] * (f_left + f_right);
        }
    }

    // Scale sums to the interval length:
    const FP_TYPE integral_kronrod = half_length * kronrod_sum;
    const FP_TYPE integral_gauss   = half_length * gauss_sum;
    const FP_TYPE error_est = fabs(integral_kronrod - integral_gauss);

    // Accept or recurse:
    if (depth >= MAX_RECURSION_DEPTH || error_est < tol) {
        // printf("Kronrod 21-10 estimate accepted at depth %d: %.16g (error: %.3g)\n", depth, integral_kronrod, error_est);
        return integral_kronrod;
    }

    const FP_TYPE left  = hypergeometric3F2_gauss_kronrod_recursive(A,   mid, a, b, c, d, e, z, 0.5 * tol, depth + 1);
    const FP_TYPE right = hypergeometric3F2_gauss_kronrod_recursive(mid, B,   a, b, c, d, e, z, 0.5 * tol, depth + 1);

    return left + right;
}

// The {}_3F_2 function evaluation using the integral representation.
static FP_TYPE hypergeometric3F2(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE z) {

    if (fabs(z) < 0.5 && fabs(z) * a * b * c < 50.0 * d * e) {
        // printf("Using series representation for {}_3F_2 with z = %f\n", z);
        return hypergeometric_3F2_series(a, b, c, d, e, z, MAX_SERIES_ITER, 1e-16);
    }

    // printf("Using adaptive Gauss-Kronrod integration for {}_3F_2 with z = %f\n", z);
    const FP_TYPE I = hypergeometric3F2_gauss_kronrod_recursive(0.0, 1.0, a, b, c, d, e, z, 1e-9, 0);

    // Multiply by the gamma factor.
    const FP_TYPE prefactor = tgamma(e) / (tgamma(c) * tgamma(e - c));
    return prefactor * I;
}

#endif //VECTOR_CASE_HYPERGEOMETRIC_H