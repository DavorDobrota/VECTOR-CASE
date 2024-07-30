cmake_minimum_required(VERSION 3.20)
project(vector_case LANGUAGES C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -march=native -mtune=native")
set(CMAKE_C_STANDARD_REQUIRED True)

add_executable(vector_case
    src/settings.h
    src/structs.h

    src/factorial_lookup.h
    src/gauss_legendre_lookup.h
    src/sum_lookup_table_near.h
    src/sum_lookup_table_far.h

    src/inductance_naive.h
    src/inductance_fast.h
    src/inductance_far.h
    src/inductance_near.h
    src/inductance_zupan.h
    src/inductance_remainder_bound.h

    src/main.c)

# Link the math library
target_link_libraries(vector_case m)