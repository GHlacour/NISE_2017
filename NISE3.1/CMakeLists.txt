cmake_minimum_required(VERSION 3.8)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

add_library(polar_lib OBJECT polar.c polar.h)
add_library(NISE3.1_lib OBJECT NISE3.1.c NISE3.1.h types3.1.h lapack.h)
add_library(NISE3.1subs_lib OBJECT NISE3.1subs.c NISE3.1subs.h types3.1.h lapack.h)
add_library(readinput_lib OBJECT readinput.c readinput.h types3.1.h)
add_library(1DFFT_lib OBJECT 1DFFT.c 1DFFT.h types3.1.h)
add_library(absorption_lib OBJECT absorption.c types3.1.h lapack.h NISE3.1subs.h 1DFFT.h)
add_library(population_lib OBJECT population.c types3.1.h lapack.h NISE3.1subs.h)
add_library(c_absorption_lib OBJECT c_absorption.c types3.1.h lapack.h NISE3.1subs.h 1DFFT.h)
add_library(calc_2DIR_lib OBJECT calc_2DIR.c types3.1.h polar.h lapack.h NISE3.1subs.h)
add_library(calc_2DES_lib OBJECT calc_2DES.c types3.1.h polar.h lapack.h NISE3.1subs.h)
add_library(luminescence_lib OBJECT luminescence.c absorption.h types3.1.h lapack.h NISE3.1subs.h)
add_library(calc_CD_lib OBJECT calc_CD.c calc_CD.h types3.1.h lapack.h NISE3.1subs.h)
add_library(analyse_lib OBJECT analyse.c types3.1.h lapack.h NISE3.1subs.h)
add_library(translate_lib OBJECT translate.c translate.h)
#add_library(2DFFT_lib OBJECT 2DFFT.c types.h)
add_library(nrutil_lib OBJECT nrutil.h nrutil.c)
add_library(random_lib OBJECT randomlib.h randomlib.c)

add_executable(NISE3.1
    $<TARGET_OBJECTS:polar_lib>
    $<TARGET_OBJECTS:NISE3.1_lib>
    $<TARGET_OBJECTS:NISE3.1subs_lib>
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
target_link_libraries(NISE3.1 m)

#add_executable(2DFFT
#	$<TARGET_OBJECTS:2DFFT_lib>
#	$<TARGET_OBJECTS:nrutil_lib>
#)
#target_link_libraries(2DFFT m)

add_executable(translate
    $<TARGET_OBJECTS:translate_lib>
    $<TARGET_OBJECTS:NISE3.1subs_lib>
    $<TARGET_OBJECTS:readinput_lib>
    $<TARGET_OBJECTS:random_lib>
)
target_link_libraries(translate m)

find_package(LAPACK REQUIRED)
target_compile_options(NISE3.1 PUBLIC ${LAPACK_LINKER_FLAGS})
#target_compile_options(2DFFT PUBLIC ${LAPACK_LINKER_FLAGS})
target_compile_options(translate PUBLIC ${LAPACK_LINKER_FLAGS})
target_link_libraries(NISE3.1 ${LAPACK_LIBRARIES})
#target_link_libraries(2DFFT ${LAPACK_LIBRARIES})
target_link_libraries(translate ${LAPACK_LIBRARIES})

find_package(FFTW REQUIRED)
target_link_libraries(NISE3.1 ${FFTW_LIBRARIES})
target_link_libraries(translate ${FFTW_LIBRARIES})

find_package(OpenMP REQUIRED)
target_link_libraries(NISE3.1 OpenMP::OpenMP_C)
target_link_libraries(calc_2DES_lib OpenMP::OpenMP_C)
#target_compile_options(calc_2DES_lib PUBLIC "${OpenMP_C_FLAGS}")