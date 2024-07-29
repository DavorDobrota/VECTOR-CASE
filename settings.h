#ifndef VECTOR_CASE_SETTINGS_H
#define VECTOR_CASE_SETTINGS_H

// Define the floating point type to be used
#define USING_DOUBLE

//#define USE_SSE
#define USE_AVX
//#define USE_AVX512

// For now just double precision is supported because float
// precision is not enough in most application and makes
// everything harder to manipulate and maintain
#if defined(USING_DOUBLE)
#define FP_TYPE double
#else
#error "No floating point type defined"
#endif

#endif //VECTOR_CASE_SETTINGS_H
