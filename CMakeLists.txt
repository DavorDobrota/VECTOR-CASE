cmake_minimum_required(VERSION 3.20)
project(vector_case LANGUAGES C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -march=native -mtune=native")
set(CMAKE_C_STANDARD_REQUIRED True)

add_executable(vector_case
        structs.h
        factorial_lookup.h
        gauss_legendre_lookup.h
        inductance_naive.h
        inductance_far.h
        inductance_zupan.h
        main.c)
add_executable(with_z
        settings.h
        structs.h
        sum_lookup_table_near.h
        factorial_lookup.h
        fast_lookup_tables.h
        inductance_near.h
        inductance_far.h
        inductance_fast.h
        inductance_remainder_bound.h
        with_z.c)

target_link_libraries(with_z m)
target_link_libraries(vector_case m)