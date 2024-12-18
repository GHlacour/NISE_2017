# Set paths
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Define libraries
add_library(random_lib OBJECT randomlib.h randomlib.c)

# Define executables
add_executable(NISE
    NISE.c NISE.h types.h lapack.h polar.c polar.h NISE_subs.c NISE_subs.h
    MPI_subs.c MPI_subs.h 1DFFT.c 1DFFT.h absorption.h absorption.c raman.h
    raman.c population.c population.h calc_DOS.h calc_DOS.c calc_2DIR.h
    calc_2DIR.c calc_2DIRraman.h calc_2DIRraman.c calc_2DES.h calc_2DES.c
    luminescence.c calc_Diffusion.h calc_Diffusion.c anisotropy.c
    anisotropy.h calc_CD.c calc_CD.h calc_LD.c calc_LD.h analyse.c
    analyse.h readinput.c readinput.h propagate.c propagate.h
    read_trajectory.c read_trajectory.h correlate.h correlate.c
    project.c project.h types_MPI.h types_MPI.c mcfret.c mcfret.h
    calc_Redfield.h calc_Redfield.h 
    sfg.c sfg.h calc_CG_2DES.c calc_CG_2DES.h eq_den.c eq_den.h   
    $<TARGET_OBJECTS:random_lib>
)

add_executable(2DFFT
    2DFFT.c types.h nrutil.h nrutil.c
)

add_executable(translate
    translate.c translate.h NISE_subs.c NISE_subs.h types.h lapack.h
    readinput.c readinput.h read_trajectory.c read_trajectory.h
    $<TARGET_OBJECTS:random_lib>
)

### Link libraries
# LAPACK
find_package(LAPACK REQUIRED)
target_link_libraries(NISE ${LAPACK_LIBRARIES})
target_link_libraries(translate ${LAPACK_LIBRARIES})

# FFTW
# Not using the FFTW3 script, since FFTW3 installed using ./configure && make && make install
# do not properly export those scripts. Therefore we use a custom script.
find_package(FFTW REQUIRED)
target_link_libraries(NISE ${FFTW_LIBRARIES})
target_include_directories(NISE PUBLIC ${FFTW_INCLUDE_DIRS})
target_link_libraries(translate ${FFTW_LIBRARIES})
target_include_directories(translate PUBLIC ${FFTW_INCLUDE_DIRS})
target_link_libraries(2DFFT ${FFTW_LIBRARIES})
target_include_directories(2DFFT PUBLIC ${FFTW_INCLUDE_DIRS})

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(NISE OpenMP::OpenMP_C)

# MPI
find_package(MPI REQUIRED)
target_link_libraries(NISE MPI::MPI_C)

# Math lib, not necessary on Windows.
# NOTE: Should be linked AFTER linking FFTW because the order matters!
# Also sets fast math
if(NOT WIN32)
    target_link_libraries(translate m)
    target_link_libraries(NISE m)
    target_link_libraries(2DFFT m)
endif()

if(NOT WIN32 AND NOT ACCURATE_MATHS)
    target_compile_options(NISE PUBLIC -ffast-math)
elseif(NOT ACCURATE_MATHS)
    target_compile_options(NISE PUBLIC /fp:fast)
endif()

## Set features
target_compile_features(NISE PUBLIC c_std_99)
target_compile_features(2DFFT PUBLIC c_std_99)
target_compile_features(translate PUBLIC c_std_99)
