cmake_minimum_required(VERSION 3.8)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in the source directory, or any directory with a CMakeLists.txt file! Please make a subdirectory like build and run from there (cmake ..)")
endif()

# Set up findFFTW to auto-detect the FFTW library
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Project definitions
project(NISE VERSION 3.1 LANGUAGES C)
add_subdirectory("src")
add_subdirectory("doc")
add_subdirectory("example")