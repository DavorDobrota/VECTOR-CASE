#include <stdio.h>

#include "inductance_radial.h"
#include "timing.h"

#define NUM_TEST_CASES 28
#define TIMING_RUNS 30


int main() {
    const CoilCalculationData coil_configurations[NUM_TEST_CASES] = {
        // Mixed test cases
        {0.1, 0.3, 0.2, 100.0, 0.5, 0.6, 0.5, 100.0},
        {0.1, 0.3, 0.6, 100.0, 0.5, 0.6, 0.2, 100.0},
        {0.1, 0.3, 0.2, 100.0, 0.7, 1.0, 0.5, 100.0},
        {0.1, 0.3, 0.6, 100.0, 0.7, 1.0, 0.2, 100.0},
        {0.1, 0.3, 0.2, 100.0, 1.0, 1.1, 0.5, 100.0},
        {0.1, 0.3, 0.6, 100.0, 1.0, 1.1, 0.2, 100.0},
        // More challenging tests
        // {0.1, 0.2, 0.1, 100.0, 0.3, 0.4, 0.1, 100.0},
        // {0.2, 0.4, 0.2, 500.0, 0.6, 0.8, 0.2, 500.0},
        // {1.0, 2.0, 1.0, 100.0, 2.0, 3.0, 1.0, 100.0},
        // {2.0, 3.0, 1.0, 100.0, 3.0, 4.0, 1.0, 100.0},
        // {2.0, 3.0, 1.0, 100.0, 3.0, 4.0, 1.0, 100.0},
        // {1.0, 2.0, 1.0, 100.0, 2.0, 3.0, 1.0, 100.0},
        // {0.1, 0.3, 0.1, 100.0, 0.3, 0.6, 0.1, 100.0},
        // {0.1, 0.3, 1.0, 100.0, 0.6, 1.0, 1.0, 100.0},
        // Long coils
        // {0.2, 0.3, 3.0, 100.0, 0.32, 0.4, 1.0, 100.0},
        // {0.2, 0.3, 3.0, 100.0, 0.32, 0.4, 1.0, 100.0},
        // Short coils
        // {0.1, 0.5, 0.2, 100.0, 0.6, 1.0, 0.1, 100.0},
        // {0.1, 0.5, 0.1, 100.0, 0.6, 1.0, 0.1, 100.0},

        // new test cases - short coils separated
        {0.3, 0.5, 0.1, 100.0, 1.3, 1.5, 0.1, 100.0},
        {0.3, 0.5, 0.1, 100.0, 0.8, 1.0, 0.1, 100.0},
        {0.3, 0.5, 0.1, 100.0, 0.7, 0.9, 0.1, 100.0},
        {0.3, 0.5, 0.1, 100.0, 0.6, 0.8, 0.1, 100.0},
        {0.3, 0.5, 0.1, 100.0, 0.6, 0.8, 0.1, 100.0},
        {0.1, 0.2, 0.1, 100.0, 0.3, 0.4, 0.1, 100.0},
        {0.2, 0.4, 0.2, 500.0, 0.6, 0.8, 0.2, 500.0},

        // new test cases - long coils
        {0.3, 0.5, 1.0, 100.0, 1.3, 1.5, 1.0, 100.0},
        {0.3, 0.5, 1.0, 100.0, 0.8, 1.0, 1.0, 100.0},
        {0.3, 0.5, 1.0, 100.0, 0.7, 0.9, 1.0, 100.0},
        {0.3, 0.5, 1.0, 100.0, 0.6, 0.8, 1.0, 100.0},
        {0.3, 0.5, 1.0, 100.0, 0.55, 0.75, 1.0, 100.0},
        {0.2, 0.3, 3.0, 100.0, 0.325, 0.425, 1.0, 100.0},
        {0.2, 0.3, 3.0, 100.0, 0.325, 0.425, 1.0, 100.0},

        // new test cases - touching coils
        {0.3, 0.5, 0.1, 100.0, 0.5, 0.7, 0.1, 100.0},
        {0.3, 0.5, 0.1, 100.0, 0.5, 0.7, 0.2, 100.0},
        {1.0, 2.0, 1.0, 100.0, 2.0, 3.0, 1.0, 100.0},
        {1.0, 2.0, 1.0, 100.0, 2.0, 3.0, 1.0, 100.0},
        {2.0, 3.0, 1.0, 100.0, 3.0, 4.0, 1.0, 100.0},
        {2.0, 3.0, 1.0, 100.0, 3.0, 4.0, 1.0, 100.0},
        {0.3, 0.5, 1.0, 100.0, 0.5, 0.7, 1.0, 100.0},
        {0.2, 0.3, 3.0, 100.0, 0.3, 0.4, 1.0, 100.0},
    };

    const FP_TYPE offsets[NUM_TEST_CASES][2] = {
        // Mixed test cases
        {-0.1, -0.3},
        {-0.3, -0.1},
        {-0.1, -0.3},
        {-0.3, -0.1},
        {-0.1, -0.3},
        {-0.3, -0.1},

        // More challenging cases
        // {0.1, 0.3},
        // {-0.1, -0.1},
        // {-0.5, 0.5},
        // {-0.5, 0.5},
        // {-0.5, -0.5},
        // {-0.5, -0.5},
        // Long coils
        // {-1.5, -0.5},
        // {0.0, -0.5},
        // Short coils
        // {-0.1, -0.05},
        // {-0.1, -0.05},

        // new test cases - short coils separated
        {-0.05, -0.05},
        {-0.05, -0.05},
        {-0.05, -0.05},
        {-0.05, -0.05},
        {-0.05, 0.15},
        {0.1, 0.3},
        {-0.1, -0.1},

        // new test cases - long coils
        {-0.5, -0.5},
        {-0.5, -0.5},
        {-0.5, -0.5},
        {-0.5, -0.5},
        {-0.5, -0.5},
        {-1.5, -0.5},
        {-1.5, -1.5},

        // new test cases - touching coils
        {-0.05, -0.05},
        {-0.05, -0.1},
        {-0.5, 0.5},
        {-0.5, -0.5},
        {-0.5, 0.5},
        {-0.5, -0.5},
        {-0.5, -0.5},
        {-1.5, -0.5},
    };

    const FP_TYPE expected_results[NUM_TEST_CASES] = {
        // Expected results for the test cases
        0.0014580449583581346,
        0.0014113009326591167,
        0.0009881598332602629,
        0.0009738954938843544,
        0.0008009113803134828,
        0.0007936416715207371,

        // More challenging
        // 0.0008454457615296840,
        // 0.07077152436821374,
        // 0.015017330663528129,
        // 0.03478118242027195,
        // 0.045806257720734,
        // 0.021453865861749712,
        // Long coils
        // 0.0008084046587807803,
        // 0.00041367134696835714,
        // Short coils
        // 0.0028833825704557822,
        // 0.0028833825704557822

        // new test cases - short coils separated
        0.00238701251929403961,
        0.00393474324718242046,
        0.00456050184487856620,
        0.00547322891675942978,
        0.00438524120720827724,
        0.000845445761529629619,
        0.0707715243682153494,

        // new test cases - long coils
        0.00211719165583916296,
        0.00299393081555929469,
        0.00324862005682069596,
        0.00354288692467293184,
        0.00370785555536516870,
        0.000806353259989952957,
        0.000710517132147122916,

        // new test cases - touching coils
        0.007016399870364664,
        0.006764435897446592,
        0.01501733066353079,
        0.02145386586173871,
        0.03478118242026727,
        0.04580625772073445,
        0.003887098293944512,
        0.000809661707357942008
    };

    FP_TYPE values[NUM_TEST_CASES] = {};
    FP_TYPE errors[NUM_TEST_CASES] = {};
    FP_TYPE execution_times[NUM_TEST_CASES] = {};


    for (uint32_t i = 0; i < NUM_TEST_CASES; i++) {
        printf("Test case %2d:\n", i + 1);
        printf("Primary coil:   r_1 = %7.4f, R_1 = %7.4f, Z_1 = %7.4f, Z_2 = %7.4f, N_1 = %d\n",
               coil_configurations[i].r_1,
               coil_configurations[i].R_1,
               offsets[i][0],
               offsets[i][0] + coil_configurations[i].L_1,
               (int)coil_configurations[i].N_1
        );
        printf("Secondary coil: r_2 = %7.4f, R_2 = %7.4f, Z_3 = %7.4f, Z_4 = %7.4f, N_2 = %d\n",
               coil_configurations[i].r_2,
               coil_configurations[i].R_2,
               offsets[i][1],
               offsets[i][1] + coil_configurations[i].L_2,
               (int)coil_configurations[i].N_2
        );

        // Calculate the mutual inductance
        const FP_TYPE M_12 = calculate_mutual_inductance_radial(
            coil_configurations[i],
            offsets[i][0],
            offsets[i][1],
            1e-16
        );

        printf("M_12 for test case %2d               : %.16g\n", i + 1, M_12);
        printf("Expected result for test case %2d    : %.16g\n", i + 1, expected_results[i]);
        printf("Relative error for test case %2d     : %.16g\n",
               i + 1,
               fabs((M_12 - expected_results[i]) / expected_results[i])
        );

        values[i] = M_12;
        errors[i] = fabs(M_12 - expected_results[i]) / expected_results[i];

        Timer timer;
        timer_start(&timer);

        FP_TYPE M_12_sum = 0.0;

        for (uint32_t j = 0; j < TIMING_RUNS; j++) {
            // Repeat the calculation to get a better average time
            M_12_sum += calculate_mutual_inductance_radial(
                coil_configurations[i],
                offsets[i][0],
                offsets[i][1],
                1e-16
            );
        }

        const double elapsed_seconds = timer_elapsed(&timer);
        const FP_TYPE execution_time = (elapsed_seconds * 1000.0) / TIMING_RUNS;

        printf("Elapsed time for test case %d: %.4f milliseconds\n\n", i + 1, execution_time);
        execution_times[i] = execution_time;
    }

    // Calculate geomean and max error
    FP_TYPE geomean = 1.0;
    FP_TYPE max_error = 0.0;
    for (uint32_t i = 0; i < NUM_TEST_CASES; i++) {
        geomean *= pow(errors[i], 1.0 / NUM_TEST_CASES);
        if (errors[i] > max_error) {
            max_error = errors[i];
        }
    }
    printf("Geomean error:      %.3e\n", geomean);
    printf("Maximum error:      %.3e\n", max_error);

    FP_TYPE geomean_time = 1.0;
    for (uint32_t i = 0; i < NUM_TEST_CASES; i++) {
        geomean_time *= execution_times[i];
    }
    geomean_time = pow(geomean_time, 1.0 / NUM_TEST_CASES);
    printf("Geomean execution time: %.3f milliseconds\n", geomean_time);
}