cmake_minimum_required(VERSION 3.12)

# Set paths
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Define libraries
add_library(random_lib OBJECT randomlib.h randomlib.c)

# Define executables
add_executable(NISE
    NISE.c NISE.h types.h lapack.h polar.c polar.h NISE_subs.c NISE_subs.h
    1DFFT.c 1DFFT.h absorption.h absorption.c population.c c_absorption.c
    calc_2DIR.c calc_2DES.c luminescence.c calc_CD.c calc_CD.h analyse.c
    readinput.c readinput.h
    $<TARGET_OBJECTS:random_lib>
)

add_executable(2DFFT
    2DFFT.c types.h nrutil.h nrutil.c
)

add_executable(translate
    translate.c translate.h NISE_subs.c NISE_subs.h types.h lapack.h readinput.c readinput.h
    $<TARGET_OBJECTS:random_lib>
)

### Link libraries
# Math lib, not necessary on Windows
if(NOT WIN32)
    target_link_libraries(translate m)
    target_link_libraries(NISE m)
    target_link_libraries(2DFFT m)
endif()

# LAPACK
find_package(LAPACK REQUIRED)
target_link_libraries(NISE ${LAPACK_LIBRARIES})
target_link_libraries(translate ${LAPACK_LIBRARIES})

# FFTW
find_package(FFTW REQUIRED)
target_link_libraries(NISE ${FFTW_LIBRARIES})
target_include_directories(NISE PUBLIC ${FFTW_INCLUDE_DIRS})
target_link_libraries(translate ${FFTW_LIBRARIES})
target_include_directories(translate PUBLIC ${FFTW_INCLUDE_DIRS})

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(NISE OpenMP::OpenMP_C)