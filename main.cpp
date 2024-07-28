#include <iostream>
#include <iomanip>
#include <chrono>

#include "inductance_naive.hpp"
#include "inductance_far.hpp"
#include "inductance_zupan.hpp"

int main() {

    // Input parameters for the problem
    constexpr FP_TYPE N_1 = 100.0;
    constexpr FP_TYPE L_1 = 0.1;
    constexpr FP_TYPE R_1 = 0.2;
    constexpr FP_TYPE r_1 = 0.1;

    constexpr FP_TYPE N_2 = 100.0;
    constexpr FP_TYPE L_2 = 0.1;
    constexpr FP_TYPE R_2 = 0.4;
    constexpr FP_TYPE r_2 = 0.3;

    constexpr FP_TYPE d = 0.001;

    // Data structure to hold the input parameters
    CoilCalculationData data{N_1, L_1, R_1, r_1, N_2, L_2, R_2, r_2};

    // Calculate the mutual inductance
//    FP_TYPE M_12 = calculate_mutual_inductance_far(data, {32, 32, 32}, d);
//    FP_TYPE M_12_unoptimized = calculate_mutual_inductance_far_unoptimized(data, {8, 8, 8}, d);
    FP_TYPE M_12_zupan = calculate_mutual_inductance_zupan(data, d, 1, 15);

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    FP_TYPE volatile value;
    for (int i = 0; i < 10000; ++i) {
        value = calculate_mutual_inductance_zupan(data, d, 1, 15);
    }

    for (int i = 1; i <= 64; ++i) {
        value = calculate_mutual_inductance_zupan(data, d, 1, i);
        std::cout << std::setprecision(16) << i << "\t" << value << std::endl;
    }

    double interval = std::chrono::duration_cast<std::chrono::duration<double>>(
            std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Value = " << value << std::endl;
    std::cout << "Zupan time average = " << interval / 10000 << " s" << std::endl;

    // Print to 16 decimal places
    std::cout << std::endl;
//    std::cout << "M_12              = " << std::setprecision(16) << M_12 << std::endl;
//    std::cout << "M_12 Unoptimized  = " << std::setprecision(16) << M_12_unoptimized << std::endl;
    std::cout << "M_12 Zupan        = " << std::setprecision(16) << M_12_zupan << std::endl;

    return 0;
}
