cmake_minimum_required(VERSION 3.20)
project(vector_case LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

if(WIN32)
    # Set the appropriate arc manually because MSVC does not have an equivalent to -march=native
    set(CMAKE_C_FLAGS "/std:c11 /O2 /arch:AVX2 /fp:fast")
else()
    set(CMAKE_C_FLAGS "-Ofast -march=native -mtune=native")
endif()

set(SOURCE_FILES
    src/settings.h
    src/structs.h

    src/factorial_lookup.h
    src/gauss_legendre_lookup.h
    src/sum_lookup_table_near.h
    src/sum_lookup_table_far.h

    src/inductance.h
    src/inductance_naive.h
    src/inductance_far.h
    src/inductance_near.h
    src/inductance_zupan.h
    src/inductance_remainder_bound.h
)

add_executable(vector_case_exec ${SOURCE_FILES} src/main.c)
add_executable(vector_case_test src/inductance.h src/test.c)

add_library(vector_case SHARED ${SOURCE_FILES})
set_target_properties(vector_case PROPERTIES LINKER_LANGUAGE C)

# Link the math library
if (UNIX)
    target_link_libraries(vector_case m)
    target_link_libraries(vector_case_test m)
endif()
