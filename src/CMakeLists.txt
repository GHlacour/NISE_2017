cmake_minimum_required(VERSION 3.9)

# Set paths
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Define libraries
add_library(polar_lib OBJECT polar.c polar.h)
add_library(NISE_subs_lib OBJECT NISE_subs.c NISE_subs.h types.h lapack.h)
add_library(readinput_lib OBJECT readinput.c readinput.h types.h)
add_library(1DFFT_lib OBJECT 1DFFT.c 1DFFT.h types.h)
add_library(absorption_lib OBJECT absorption.c types.h lapack.h NISE_subs.h 1DFFT.h)
add_library(population_lib OBJECT population.c types.h lapack.h NISE_subs.h)
add_library(c_absorption_lib OBJECT c_absorption.c types.h lapack.h NISE_subs.h 1DFFT.h)
add_library(calc_2DIR_lib OBJECT calc_2DIR.c types.h polar.h lapack.h NISE_subs.h)
add_library(calc_2DES_lib OBJECT calc_2DES.c types.h polar.h lapack.h NISE_subs.h)
add_library(luminescence_lib OBJECT luminescence.c absorption.h types.h lapack.h NISE_subs.h)
add_library(calc_CD_lib OBJECT calc_CD.c calc_CD.h types.h lapack.h NISE_subs.h)
add_library(analyse_lib OBJECT analyse.c types.h lapack.h NISE_subs.h)
add_library(nrutil_lib OBJECT nrutil.h nrutil.c)
add_library(random_lib OBJECT randomlib.h randomlib.c)

# Define executables
add_executable(NISE
    NISE.c NISE.h types.h lapack.h
    $<TARGET_OBJECTS:polar_lib>
    $<TARGET_OBJECTS:NISE_subs_lib>
    $<TARGET_OBJECTS:1DFFT_lib>
    $<TARGET_OBJECTS:absorption_lib>
    $<TARGET_OBJECTS:population_lib>
    $<TARGET_OBJECTS:c_absorption_lib>
    $<TARGET_OBJECTS:calc_2DIR_lib>
    $<TARGET_OBJECTS:calc_2DES_lib>
    $<TARGET_OBJECTS:luminescence_lib>
    $<TARGET_OBJECTS:calc_CD_lib>
    $<TARGET_OBJECTS:analyse_lib>
    $<TARGET_OBJECTS:readinput_lib>
    $<TARGET_OBJECTS:random_lib>
)

add_executable(2DFFT
    2DFFT.c types.h
    $<TARGET_OBJECTS:nrutil_lib>
)

add_executable(translate
    translate.c translate.h
    $<TARGET_OBJECTS:NISE_subs_lib>
    $<TARGET_OBJECTS:readinput_lib>
    $<TARGET_OBJECTS:random_lib>
)

### Link libraries
# Math lib
target_link_libraries(translate m)
target_link_libraries(NISE m)
target_link_libraries(2DFFT m)

# LAPACK
find_package(LAPACK REQUIRED)
target_link_libraries(NISE ${LAPACK_LIBRARIES})
target_link_libraries(translate ${LAPACK_LIBRARIES})

# FFTW
find_package(FFTW REQUIRED)
target_link_libraries(NISE ${FFTW_LIBRARIES})
target_link_libraries(translate ${FFTW_LIBRARIES})

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(NISE OpenMP::OpenMP_C)
target_link_libraries(calc_2DES_lib OpenMP::OpenMP_C)