#ifndef VECTOR_CASE_STRUCTS_H
#define VECTOR_CASE_STRUCTS_H

#include <stdint.h>

#include "settings.h"

const FP_TYPE local_pi = 3.14159265358979323846;


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

typedef struct {
    uint32_t k_terms;
    uint32_t l_terms;
    uint32_t n_terms;

} SumPrecisionData;


#endif //VECTOR_CASE_STRUCTS_H
