#ifndef VECTOR_CASE_STRUCTS_H
#define VECTOR_CASE_STRUCTS_H

#include <stdint.h>

#include "settings.h"

const FP_TYPE local_pi = 3.14159265358979323846;


/**
 * @struct CoilCalculationData
 * @brief A struct that contains the data needed to calculate the mutual inductance between two coils.
 *
 * This structure holds the geometrical and winding properties of two coils,
 * which are required for computing their mutual inductance.
 *
 * @var CoilCalculationData::r_1
 * The inner radius of the first coil.
 *
 * @var CoilCalculationData::R_1
 * The outer radius of the first coil.
 *
 * @var CoilCalculationData::L_1
 * The length of the first coil.
 *
 * @var CoilCalculationData::N_1
 * The number of turns in the first coil.
 *
 * @var CoilCalculationData::r_2
 * The inner radius of the second coil.
 *
 * @var CoilCalculationData::R_2
 * The outer radius of the second coil.
 *
 * @var CoilCalculationData::L_2
 * The length of the second coil.
 *
 * @var CoilCalculationData::N_2
 * The number of turns in the second coil.
 */
typedef struct {
    FP_TYPE r_1;
    FP_TYPE R_1;
    FP_TYPE L_1;
    FP_TYPE N_1;

    FP_TYPE r_2;
    FP_TYPE R_2;
    FP_TYPE L_2;
    FP_TYPE N_2;
    
} CoilCalculationData;


/**
 * @struct SumPrecisionData
 * @brief A struct that contains the number of terms in the series expansion.
 *
 * @var SumPrecisionData::k_terms
 * The number of terms in the series expansion for the k index.
 *
 * @var SumPrecisionData::l_terms
 * The number of terms in the series expansion for the l index.
 *
 * @var SumPrecisionData::n_terms
 * The number of terms in the series expansion for the n index.
 */
typedef struct {
    uint32_t k_terms;
    uint32_t l_terms;
    uint32_t n_terms;

} SumPrecisionData;


#endif //VECTOR_CASE_STRUCTS_H