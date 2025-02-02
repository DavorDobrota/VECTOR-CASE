#include "inductance.h"

#define NUM_TEST_CASES 29


/**
 * This function was used to generate all the test cases showcase in the paper. It is not
 * meant to be used in any other context.
 *
 * The output is very verbose and is meant to provide detailed test reports.
 */
int main() {
    const CoilCalculationData coil_configurations[NUM_TEST_CASES] = {
        // Small separation
        {0.1, 0.2, 0.1, 100.0, 0.3, 0.4, 0.1, 100.0},
        {0.5, 1.5, 1.0, 1.0, 0.5, 1.5, 1.0, 1.0},
        {0.0875, 0.1125, 0.025, 200.0, 0.0875, 0.1125, 0.025, 200.0},
        {0.0287, 0.0507, 0.0302, 10.0, 0.0287, 0.0507, 0.0302, 10.0},
        {0.0287, 0.0507, 0.0302, 10.0, 0.0287, 0.0507, 0.0302, 10.0},
        {1.0, 2.0, 1.0, 100.0, 2.0, 4.0, 1.0, 100.0},
        {1.0, 2.0, 5.0, 100.0, 1.0, 2.0, 1.0, 100.0},
        {0.5, 1.0, 5.0, 100.0, 1.0, 2.0, 1.0, 100.0},
        // Medium separation
        {0.1, 0.2, 0.1, 100.0, 0.3, 0.4, 0.1, 100.0},
        {0.1, 0.2, 0.1, 100.0, 0.3, 0.4, 0.1, 100.0},
        {0.1, 0.2, 0.1, 100.0, 0.3, 0.4, 0.1, 100.0},
        {0.5, 1.5, 1.0, 100.0, 0.5, 1.5, 1.0, 100.0},
        {0.5, 1.5, 1.0, 100.0, 0.5, 1.5, 1.0, 100.0},
        {1.0, 2.0, 1.0, 100.0, 2.0, 4.0, 1.0, 100.0},
        {1.0, 2.0, 5.0, 100.0, 1.0, 2.0, 1.0, 100.0},
        {0.5, 1.0, 5.0, 100.0, 1.0, 2.0, 1.0, 100.0},
        // Large separation
        {0.1, 0.2, 0.1, 100.0, 0.3, 0.4, 0.1, 100.0},
        {0.5, 1.5, 1.0, 1.0, 0.5, 1.5, 1.0, 1.0},
        {1.0, 2.0, 1.0, 100.0, 2.0, 4.0, 1.0, 100.0},
        {1.0, 2.0, 5.0, 100.0, 1.0, 2.0, 1.0, 100.0},
        {0.5, 1.0, 5.0, 100.0, 1.0, 2.0, 1.0, 100.0},
        // Additional test cases
        {2.0, 3.0, 1.0, 100.0, 1.0, 2.0, 1.0, 100.0},
        {2.0, 3.0, 1.0, 100.0, 2.0, 3.0, 1.0, 100.0},
        {2.0, 3.0, 1.0, 100.0, 3.0, 4.0, 1.0, 100.0},
        {2.0, 3.0, 1.0, 100.0, 3.0, 4.0, 1.0, 100.0},
        {3.0, 4.0, 1.0, 100.0, 2.0, 3.0, 1.0, 100.0},
        {2.0, 3.0, 1.0, 100.0, 2.0, 3.0, 1.0, 100.0},
        {1.0, 2.0, 1.0, 100.0, 2.0, 3.0, 1.0, 100.0},
        {1.0, 2.0, 1.0, 100.0, 2.0, 3.0, 1.0, 100.0},
    };

    const FP_TYPE distances[NUM_TEST_CASES] = {
        0.1, 0.001, 0.035, 1.0e-5, 0.0098, 4.0e-6, 0.05, 5.0e-6,
        0.2, 0.4, 0.8, 1.0, 2.0, 3.0, 2.0, 2.0,
        1.1501, 4.501, 12.0, 5.001, 5.001,
        1e-6, 1e-6, 1e-6, -1.0, 1e-6, 1e-6, 1e-6, -1.0
    };

    const SumPrecisionData precision = {32, 32, 64};

    const FP_TYPE ref_values[NUM_TEST_CASES] = {
        // Small separation
        0.000845445761529684,
        0.000000539456018627887,
        0.003737280536931003,
        0.000003108645278459047,
        0.000002117124447694834,
        0.013581683586999548,
        0.0041692155709304005,
        0.0014354523410881740,
        // Medium separation
        0.0005486756620356269,
        0.00023561180514417563,
        0.00006107543644963594,
        0.0016091670033369391,
        0.0006331887757973942,
        0.0030260151749000406,
        0.0009581220397186723,
        0.00027690370848532092,
        // Large separation
        0.00002546539035257568,
        0.000000012588109863169327,
        0.00017594575209956676,
        0.00022333836706366272,
        0.00005889356521217686,
        // Additional test cases
        0.01501732168425381,
        0.03398632339317788,
        0.03478116653856178,
        0.04580625811880162,
        0.03478116653856178,
        0.03398632339317788,
        0.01501732168425381,
        0.021453866158409343
    };

    FP_TYPE values[NUM_TEST_CASES] = {};
    FP_TYPE errors[NUM_TEST_CASES] = {};

    for (int i = 0; i < NUM_TEST_CASES; i++) {
        printf("Test case %d\n", i + 1);
        printf("r_1 = %.16g, R_1 = %.16g, L_1 = %.16g, N_1 = %.16g\n",
               coil_configurations[i].r_1, coil_configurations[i].R_1, coil_configurations[i].L_1, coil_configurations[i].N_1
        );
        printf("r_2 = %.16g, R_2 = %.16g, L_2 = %.16g, N_2 = %.16g\n",
               coil_configurations[i].r_2, coil_configurations[i].R_2, coil_configurations[i].L_2, coil_configurations[i].N_2
        );
        printf("d   = %.16g\n", distances[i]);

        FP_TYPE M_12 = calculate_mutual_inductance(coil_configurations[i], precision, distances[i], true);
        values[i] = M_12;
        errors[i] = fabs(M_12 - ref_values[i]) / ref_values[i];

        printf("M_12    = %.16g\n", M_12);
        printf("M_ref   = %.16g\n", ref_values[i]);
        printf("Rel err = %.3e\n\n", errors[i]);
    }

    printf("Summary\n");
    for (int i = 0; i < NUM_TEST_CASES; i++) {
        printf("Test case %d: value = %21.16g, error = %.2e\n", i + 1, values[i], errors[i]);
    }

    // Calculate geomean of the errors
    FP_TYPE geomean = 1.0;
    for (int i = 0; i < NUM_TEST_CASES; i++) {
        geomean *= errors[i];
    }
    geomean = pow(geomean, 1.0 / NUM_TEST_CASES);
    printf("Geometric mean of the errors = %.3e\n", geomean);

    return 0;
}