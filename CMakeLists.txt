cmake_minimum_required(VERSION 3.20)
project(vector_case)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -mtune=native")

add_executable(vector_case
        structs.h
        factorial_lookup.h
        gauss_legendre_lookup.h
        inductance_naive.hpp
        inductance_far.hpp
        inductance_zupan.hpp
        main.cpp)
add_executable(with_z
        settings.h
        structs.h
        sum_lookup_table_near.h
        factorial_lookup.h
        fast_lookup_tables.h
        inductance_near.hpp
        inductance_far.hpp
        inductance_fast.hpp
        remainder_estimate.hpp
        with_z.cpp)