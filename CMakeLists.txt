cmake_minimum_required(VERSION 3.8)

configure_file(findFFTW.cmake findFFTW/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW)

if(result)
    message(FATAL_ERROR "Cmake step for findFFTW failed: ${result}")
else()
    message("CMake step for findFFTW completed (${result})")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW)

if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")

project("NISE3.1")
add_subdirectory("NISE3.1")
